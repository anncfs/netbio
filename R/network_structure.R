## Functions to analyse network structure (communities and components)

#' Check network components
#'
#' Checks connected components present in network
#'
#' @param network An igraph representation of network
#' @param filename An optional filename (default NULL) to save simple stats about graph and components
#'
#' @import igraph
#'
#' @return Logical vector indicating which nodes belong to largest component (TRUE) or not (F)
#' @export
get_components <- function(network, filename = NULL) {

  netcomp = igraph::components(network)
  x = table(netcomp$csize)

  # prints results to terminal of file
  message('No. of nodes: ', gorder(network))
  message('No. of components: ', netcomp$no)
  message('Max. component size: ', max(netcomp$csize))

  # save to file
  if (is.character(filename)) {
    write(paste0('Total: ', netcomp$no, ' components, ', gorder(network), ' nodes'), filename, append = T)
    write(paste0('Size:\t',paste0(names(x), collapse = '\t')), filename, append = T)
    write(paste0('Freq:\t',paste0(as.vector(x), collapse = '\t')), filename,  append = T)
  }

  # returns logical vector indicating which nodes in largest component
  return(netcomp$membership == which.max(netcomp$csize))
}


#' Detect communities in network
#'
#' Uses selected community detection algorithm to find densely connected regions in network
#'
#' @param network An igraph representation of network
#' @param cluster_algorithm igraph community detection algorithm ('fg' = fast_greedy,
#'   'ml' = louvain (DEFAULT), 'eb' = edge betweeness, 'wt' = walktrap)
#' @param filename Filename to write summary results about community detection (default = NULL)
#' @param min_size Minimum size of cluster to record details in \code{filename} (default = 20)
#'
#' @import igraph
#'
#' @return List of community membership for nodes
#' @export
get_communities <- function(network, filename = NULL, cluster_algorithm = 'ml', min_size = 10) {

  if (cluster_algorithm == 'ml') {
    clust = igraph::cluster_louvain(network)
  } else if (cluster_algorithm == 'fg') {
    clust = igraph::cluster_fast_greedy(network)
  } else if (cluster_algorithm == 'eb') {
    clust = igraph::cluster_edge_betweenness(network)
  } else if (cluster_algorithm == 'wt') {
    clust = igraph::cluster_walktrap(network)
  }

  # save details to file
  if (is.character(filename)) {
    write('\nCommunity detection results:', filename, append = T)
    write(paste0('Method: ', igraph::algorithm(clust), '\t'), filename, append = T)
    write(paste0('Number of clusters: ', length(clust)), filename, append = T)

    clust_size = table(clust$membership)
    idx = which(clust_size >= min_size)
    write(paste0('Cluster ID:\t',paste0(names(clust_size)[idx], collapse = '\t')), filename, append = T)
    write(paste0('Cluster size:\t',paste0(clust_size[idx], collapse = '\t')), filename,  append = T)
  }

  # print results
  message('No. of clusters: ', length(clust))

  # returns named vector with community membership
  netclust = clust$membership
  names(netclust) = clust$names
  return(netclust)
}



#' Plot to compare network communities
#'
#' 2D heatmap plot to compare (using Jaccard Index) the nodes present in two sets of communities
#'
#' @param net1 Named list giving community membership for nodes (list names) in network 1
#' @param net2 Named list giving community membesrhip for nodes in network 2 (must have same nodes as network 1)
#' @param min_size Minimum size (no. of nodes) of communities to include in comparison (default = 10)
#' @param network_names List giving names for networks to use in plot labels
#'   (default = c("network 1", "network 2"))
#'
#' @importFrom graphics axis title
#' @export
#' @return NULL
plot_community_comparison <- function(net1, net2, min_size = 10, network_names = c('network 1', 'network 2')) {

  # check net1 and net2 have same nodes in same order
  if (!all(names(net1) == names(net2))) {
    message("Ensure nodes of both community lists match and are in the same order")
    return(NULL)
  }

  # get cluster IDs for clusters above min size
  id1 = as.vector(which(table(net1) >= min_size))
  id2 = as.vector(which(table(net2) >= min_size))

  # calculate jaccard index for each pair of clusters
  mat = matrix(0, length(id1), length(id2))
  rownames(mat) = id1; colnames(mat) = id2

  for (i in 1:length(id1)){
    for (j in 1:length(id2)){
      # find gene members for each cluster
      c1_i = which(net1 == id1[i])
      c2_j = which(net2 == id2[j])

      # calculate jaccard coeff
      a = intersect(c1_i, c2_j)
      b = union(c1_i, c2_j)
      mat[i,j] = length(a)/length(b)
    }
  }

  # if at least 2 clusters for each network, plot graph
  if(all(dim(mat) > 1)){
    fields::image.plot(1:ncol(mat), 1:nrow(mat), t(mat[nrow(mat):1,]),
               xaxt = 'n', yaxt = 'n', col = grDevices::gray(10:0/10), zlim = c(0,1),
               xlab = paste0('clusters (', network_names[2], ')'),
               ylab = paste0('clusters (', network_names[1], ')'))
    axis(1, at = c(1:ncol(mat)), labels = id2, las=2, cex.axis = 0.8)
    axis(2, at = c(1:nrow(mat)), labels = rev(id1), las=2, cex.axis = 0.8)
    title(paste0('Cluster comparison (min cluster size = ', min_size, ')'))
  }
}
