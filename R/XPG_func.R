#' NN Connection
#'
#' Calculating Percentage Neighbourhood Connection
#' @return The average score per cell type and average percentage per each cell
#' @export
nnfunc = function(s, graph, label, k, cat){
  l = as.character(s@meta.data[,label])
  nn = as(s@graphs[[graph]], 'matrix')
  nn = t(rowsum(t(nn), l))
  per = mapply(function(x,y) nn[x,y], 1:dim(nn)[1], l)/k*100
  if(cat == 'group'){per = tapply(per, l, mean)}
  return(c(per, mean(per)))
}

#' PCA and NN functions
#'
#' Constructing PCA space and drawing NN graph using genes of interest
#' @return neighbourhood percentage connection matrix
#' @export
insider = function(gi = gi, seu = seu, graph = graph, sp = sp, arr=arr, ct = ct, k = 10, cat = 'group'){
  gene = c(sp, arr, gi)
  tmpseu <- seu[Features(seu) %in% gene,]
  tmpseu <- RunPCA(tmpseu, npcs = 12, features = gene, verbose = F)
  tmpseu <- FindNeighbors(tmpseu, reduction = "pca", k.param = k, dims = 1:10, compute.SNN = TRUE, verbose = F)
  tmpseu_nn = nnfunc(tmpseu, graph, ct, k, cat)
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
#' @param k number neighbour per cell
#' @param cat average per cell type = 'group' or per cell = 'cell'
#' @return a list of 3 objects (1. Percentage matrix 2. Ranked gene list 3. Time(s) per run)
#' @examples
#' output = main_func(g = c('EPCAM','MS4A1','MZB1'), sp = c('KLF2','CD8A','COL1A1','EPCAM'), graph = 'RNA_nn', seu = obj, ct = 'celltype', core = 100 , k = 10, cat = 'group');
#' @export
XPG = function(g, sp, graph, seu, ct, core, k, cat){
  arr = c(); mer = data.frame(); curr = 0; outcast = matrix()
  it = c(); jt = c()
  startlen = length(g)
  l = levels(seu@meta.data[,ct])

  tryCatch({
    for(j in 1:startlen){
      print(paste('Round', j))
      jst <- Sys.time()
      suppressMessages(suppressWarnings({
        df = parallel::mclapply(g, FUN = mclapplyfunc, mc.cores = core, mc.preschedule = FALSE, seu = seu, graph = graph, sp = sp, arr = arr, ct = ct, k = k, cat = cat)
        if (is.list(df)) df = do.call(rbind, df)
        df = cbind(df,g,rep(j,dim(df)[1]))
        colnames(df) = c(l,'average','gene','round')
        arr = c(arr, df[,'gene'][which.max(df[,'average'])])
        g = g[!g %in% arr]
      }))
      mer = rbind(mer,df)

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
