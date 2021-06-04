#' Prepare data for visualize in less dimention
#' @export
#' @param data a data matrix with n rows representing samples and m columns representing genes.
#' @param alg.type if alg.type == "tSNE" to reduce the dimension of the data, the tSNE method is used (the default). if alg.type == "PCA" the PCA method is used and if alg.type == "UMAP" the UMAP method is used.
#' @return Double matrix of components (PCA) or dimensions (tSNE, UMAP)
#' @examples
#' data = rbind(matrix(rnorm(1000, mean=50), nrow = 500, ncol = 6),
#'              matrix(rnorm(100), ncol = 6, nrow = 500))
#' labels = factor(cbind(rep(1, 500), rep(2, 500)))
#
#' data.tsne <- Rtsne(data, dims = 2, perplexity=30, verbose=TRUE,
#'                    max_iter = 500)
#'
#' ggplot(data.frame(data.tsne$Y), aes(x = data.tsne$Y[,1], y = data.tsne$Y[,2])) +
#'   geom_point(aes(color = labels))+
#'   ylab("dim2")+
#'   xlab("dim1")+
#'   ggtitle("tSNE")
forvisualData = function(data, alg.type = c("tSNE", "PCA", "UMAP")){
  if (alg.type=="tSNE") {
    set.seed(1)
    tsne <- Rtsne(data, dims = 2, perplexity=30, verbose=TRUE, max_iter = 500)
    getdata = tsne$Y
  }
  else if (alg.type=="PCA"){
    pca = prcomp(data)
    getdata = pca$x
  }
  else{
    umap <- umap(data)
    getdata <- umap$layout
  }
  return(getdata)
}
