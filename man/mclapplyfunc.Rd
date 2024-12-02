% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/XPG_func.R
\name{mclapplyfunc}
\alias{mclapplyfunc}
\title{Intermediate function to prevent error from mclapply function}
\usage{
mclapplyfunc(
  gi = gi,
  seu = seu,
  graph = graph,
  sp = sp,
  arr = arr,
  ct = ct,
  k = k,
  cat = "group"
)
}
\arguments{
\item{gi}{Current gene of interest}
}
\value{
neighbourhood percentage connection matrix
}
\description{
tryCatch error and rerun insider functions
}
