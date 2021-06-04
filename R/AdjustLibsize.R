#' Adjusting the technical effect of the library size in RNA-Seq data
#' @export
#' @param X a data matrix with n rows representing samples and m columns representing genes.
#' @param libsize numeric vector of sample library size values. Length - number of samples.
#' @return A list of the following components:
#' newX - matrix of corrected data
#' compLibsizeNum - component number with library size
adjustLibsize = function(
  X,
  libsize
)
{
  findcomp = findComp(X)

  S = findcomp$S
  A = findcomp$A
  compLibsizeNum = findcomp$compLibsizeNum

  S[, compLibsizeNum] = S[, compLibsizeNum] / libsize

  newX = data.matrix(S) %*% data.matrix(A)

  return(list(newX, compLibsizeNum))
}










