#' Find the number of component with batch-effect and component with library size
#' @export
#' @param X a data matrix with n rows representing samples and m columns representing genes.
#' @param libsize numeric vector of sample library size values. Length - number of samples.
#' @param batch numeric vector of batch value for sample. Length - number of samples.
#' @param ncomp number of components
#' @return A list of the following components:
#' compBatchNum - component number with batch effect
#' compLibsizeNum - component number with library size
#' S - weight matrix
#' A - metagene matrix
findComp = function(
  X,
  libsize=0,
  batch=0,
  ncomp=5,
  fun="logcosh",
  alg.typ="parallel",
  method = c("R","C")
)
{
  set.seed(1)
  X.ica = fastICA(X, n.comp=ncomp, alg.typ=alg.typ, fun=fun, method=method)

  S = X.ica$S
  A = X.ica$A

  corLibsizeS <- cor(libsize, S)
  compLibsizeNum = which.max(abs(corLibsizeS))

  res = matrix(0, ncomp, nlevels(factor(batch)))
  for (i in 1:nlevels(factor(batch))){
    res[, i] = data.matrix(colMeans(S[batch == levels(factor(batch))[i], ]))
  }

  batchDiff <- res[,1] - res[,2]

  compBatchNum = which.max(abs(batchDiff))

  return(list(compLibsizeNum, compBatchNum, S, A))
}

