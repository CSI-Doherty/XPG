#' extension of lapply function
#'
#' rewrite lapply function so that it returns a named vector
#' @param X vector names
#' @param FUN function in lapply
#' @return A named vector
#' @export
lapply2 <- function (X, FUN, ...)  {
  names <- X
  FUN <- match.fun(FUN)
  if (!is.vector(X) || is.object(X)) {
    X <- as.list(X)
  }
  names(X) <- X
  .Internal(lapply(X, FUN))
}

#' NN Connection
#'
#' Calculating Percentage Neighbourhood Connection 
#' @return The average score per cell type and average percentage per each cell  
#' @export
nnfunc = function(s, graph, label, cat){
  nn = s@graphs[[graph]]
  tmpll = data.frame(lapply2(levels(s@meta.data[,label]), function(x){
    cellid = rownames(s@meta.data[s@meta.data[,label]==x,])
    print(nn[,colnames(nn) %in% cellid])
    return(rowSums(nn[,colnames(nn) %in% cellid]))
  }), check.names = FALSE)
  s@meta.data[,label, drop = FALSE] %>% as_tibble(rownames = 'barcodes') %>% left_join(tmpll %>% as_tibble(rownames = 'barcodes')) -> tmpll

  tmpll %>% select(!barcodes) %>% group_by(!!rlang::sym(label)) %>% summarise_each(sum) %>% melt(id.vars = label) %>%
    group_by(!!rlang::sym(label)) %>% mutate(per =  100 *value/sum(value)) -> m
  l = m[m[[label]] == m$variable, ]$per

  if(cat == 'group'){
    return(c(l, mean(l)))
  }else if(group == 'cell'){
    ct_ll = as.numeric(apply(tmpll, 1, function(x) x[x[label]]))
    rs = rowSums(tmpll[,levels(s@meta.data[,label])])
    avg = sum(ct_ll/rs*100)/dim(tmpll)[1]
    return(c(l, avg))
  }
}

#' PCA and NN functions
#'
#' Constructing PCA space and drawing NN graph using genes of interest
#' @return neighbourhood percentage connection matrix
#' @export
insider = function(gi = gi, seu = seu, graph = graph, sp = sp, arr=arr, ct = ct, k = 10, cat = 'group'){
  gene = c(sp, arr, gi)
  tmpseu <- seu[Features(seu) %in% gene,]
  tmpseu <- RunPCA(tmpseu, npcs = 5, features = gene, verbose = F)
  tmpseu <- FindNeighbors(tmpseu, reduction = "pca", k.param = k, dims = 1:5, compute.SNN = TRUE, verbose = F)
  tmpseu_nn = nnfunc(tmpseu, graph, ct, cat)
  return(tmpseu_nn)
}

#' Intermediate function to prevent error from mclapply function
#'
#' tryCatch error and rerun insider functions
#' @param gi Current gene of interest
#' @return neighbourhood percentage connection matrix
#' @export
mclapplyfunc = function(gi = gi, seu = seu, graph = graph, sp = sp, arr=arr, ct = ct, k = k, cat = 'group'){
  tryCatch({
    return(insider(gi = gi, seu = seu, graph = graph, sp = sp, arr=arr, ct = ct, k = k, cat = 'group'))
  }, error = function(e) {
    return(insider(gi = gi, seu = seu, graph = graph, sp = sp, arr=arr, ct = ct, k = k, cat = 'group'))
  })
}

#' Main XPG function 
#'
#' Initialise all empty variables
#' @param g Genes of interest
#' @param sp Starting panel gene list
#' @param graph Nearest Neighbour graph in Seurat object
#' @param seu Seurat object
#' @param ct celltype column in Seurat object
#' @param core number of available cores
#' @param cat average per cell type = 'group' or per cell = 'cell'
#' @return a list of 3 objects (1. Percentage matrix 2. Ranked gene list 3. Time(s) per run)
#' @examples 
#' output = main_func(g = c('EPCAM','MS4A1','MZB1'), sp = c('KLF2','CD8A','COL1A1','EPCAM'), graph = 'RNA_nn', seu = obj, ct = 'celltype', core = 100 , k = 10, cat = 'group');
#' @export
XPG = function(g, sp, graph, seu, ct, core, k, cat){
  arr = c(); mer = data.frame(); curr = 0; outcast = matrix()
  it = c(); jt = c()
  startlen = length(g)

  tryCatch({
    for(j in 0:startlen){
      print(paste('Round', j))
      jst <- Sys.time()
      suppressMessages(suppressWarnings({
        if(j==0){
          df = insider(gi = c(), seu = seu, graph = graph, sp = sp, arr = arr, ct = ct, k = k, cat = cat)
          class(df) <- "numeric"
          df = t(data.frame(c(df,'sp',0)))
        }else{
          df = mclapply(g, FUN = mclapplyfunc, mc.cores = core, mc.preschedule = T, seu = seu, graph = graph, sp = sp, arr = arr, ct = ct, k = k, cat = cat)
          if (is.list(df)) df = do.call(rbind, df)
          class(df) <- "numeric"
          df = cbind(df,g,rep(j,dim(df)[1]))
        }
        colnames(df) = c(levels(seu@meta.data[,ct]),'average','gene','round')

      }))
      mer = rbind(mer,df)
      arr = c(arr, df[which.max(df[,'average']),'gene'])
      g = g[!g %in% arr]

      jet <- Sys.time()
      t = difftime(jet, jst, units = 'secs')
      jt = c(jt,t)

      print(sprintf(paste('Time @', j,'genes: %.2f %s'), t, units(t)))

    }
  },error=function(e) {
    message('An Error Occurred');print(e)
  },warning=function(w) {
    message('A Warning Occurred');print(w)
  })
  return(list(mer, arr, jt))
}
