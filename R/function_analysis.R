## functional analysis functions

#' Fisher's exact test
#'
#' wrapper function for R's fisher.test function to test each cluster identified in graph for overlap with genes of interest
#'
#' @param net_clusters Named list of cluster membership for each gene
#' @param int_genes List of all genes in category of interest (any not in reference list are ignored)
#' @param min_size Min size of cluster to analyse for overlap (default = 10)
#' @param all_genes List of all genes to consider as reference list (i.e. full population under analysis, default = NULL in
#'   which case all nodes in net_clusters used)
#' @param filename Filename to write results to, default = NULL (no file saved)
#'
#' @export
#'
fisher_test_clusters <- function(net_clusters, int_genes, min_size = 10, all_genes = NULL, filename = NULL) {

  # check int_genes and discard any not in reference list
  if (is.null(all_genes)) { all_genes = names(net_clusters) }
  int_genes = int_genes[which(int_genes %in% all_genes)]

  # get list of cluster IDs above min size
  idx = which(table(net_clusters) >= min_size)

  # table for results
  mat = matrix(ncol = 3, nrow = length(idx), dimnames = list(idx, c("cluster", "p-value", "p-adjust")))
  mat[,1] = idx

  for (i in idx) {

    i_genes = names(which(net_clusters == i))

    dfx = data.frame(cluster_i = all_genes %in% i_genes, int_genes = all_genes %in% int_genes)

    x = stats::fisher.test(table(dfx), alternative = "greater")

    mat[which(idx == i),2] = x$p.value
  }

  mat[,3] = stats::p.adjust(mat[,2], method = "BH")

  if (!is.null(filename)) {
    utils::write.table(mat, file = filename, sep = "\t", quote = FALSE, row.names = FALSE)
  }

  return(mat)
}
