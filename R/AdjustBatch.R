#' Adjusting the technical effect of the batch effect in RNA-Seq data
#' @export
#' @param X a data matrix with n rows representing samples and m columns representing genes.
#' @param batch numeric vector of batch value for sample. Length - number of samples.
#' @return A list of the following components:
#' newX - matrix of corrected data
#' compBatchNum - component number with batch effect
adjustBatch = function(
  X,
  batch
)
{
  findcomp = findComp(X)

  S = findcomp$S
  A = findcomp$A
  compBatchNum = findcomp$compBatchNum

  res = 0
  for (i in 1:nlevels(factor(batch))){
    res[i] = colMeans(S[batch == levels(factor(batch))[i], compBatchNum])
    S[batch == levels(factor(batch))[i], compBatchNum] = S[batch == levels(factor(batch))[i], compBatchNum] / res[i]
  }

  newX = data.matrix(S) %*% data.matrix(A)

  return(list(newX, compBatchNum))
}



