## functions for cluster comparisons

#' Plot to compare network communities
#'
#' Compare edge or node values between networks for communities identified
#'
#' @param network1 Data.frame with all network edges (before thresholding),
#'   columns give node1, node2, edge score
#' @param network2 As above for second network
#' @param net_clusters Named list of clusters
#' @param min_size Min size of clusters to include in analysis (default = 10)
#'
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' @importFrom stats wilcox.test
#' @importFrom rlang .data
#'
#' @export
plot_community_scores <- function(network1, network2, net_clusters, min_size = 10){

  idx = which(table(net_clusters) >= min_size)

  mat = matrix(nrow = length(idx), ncol = 4,
               dimnames = list(idx, c("mean 1", "mean 2", "p-value", "p-adjust")))
  mat[,1] = idx

  # for each cluster, test for sig difference in edge scores (via Mann Whitney)
  p = list()  # list to allow multiple subplots (TODO improve this in future - use facets?)

  for (i in idx) {

    i_genes = names(which(net_clusters == i))

    # index for edges within cluster i in each network
    i1 = which(network1[,1] %in% i_genes & network1[,2] %in% i_genes)
    i2 = which(network2[,1] %in% i_genes & network2[,2] %in% i_genes)

    # mann whitney (wilcox) test
    mw = wilcox.test(network1[i1,3], network2[i2,3])

    mat[which(idx == i), 1:3] = c(mean(network1[i1,3]),
                                  mean(network2[i2,3]),
                                  mw$p.value)

    # violin plot
    temp = data.frame(
      network = c(rep("1", length(i1)), rep("2", length(i2))),
      edge_score = c(network1[i1,3], network2[i2,3]))

    p[[which(idx == i)]] = ggplot(temp, aes(x = .data$network, y = .data$edge_score, color = .data$network)) + geom_violin() +
      geom_boxplot(width = 0.1) + labs(title = paste0('cluster ', i)) +
      theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

  }

  do.call(gridExtra::grid.arrange, p)
  mat[,4] = stats::p.adjust(mat[,3], method = "BH")

  return(mat)
}



#' Select nodes in cluster list
#'
#' Takes a network and list of cluster identities (or min size of cluster) and returns a named
#' logical list indicating whether each gene is in the selected list of clusters
#'
#' @param net_clusters Named list of clusters
#' @param min_size Min size of clusters to include in analysis (default = 10)
#' @param select_clusters List of cluster IDs to select (default NULL, in which case min_size used to select clusters)
#'
#' @return A named logical list, where names are node names, and vector indicates presence of nodes in selection
#'
#' @export
get_cluster_nodes <- function(net_clusters, select_clusters = NULL, min_size = 10) {

  if (is.null(select_clusters)) {
    select_clusters = which(table(net_clusters) >= min_size)
  }

  nodes_select = logical(length(net_clusters))
  names(nodes_select) = names(net_clusters)

  for (i in select_clusters) {
    idx_i = which(net_clusters == i)
    nodes_select[idx_i] = TRUE
  }
  return(nodes_select)
}


#' Get combined network
#'
#' Create igraph network with edge weights to indicate shared / unique edges and selected nodes
#'
#' @param net1 Network 1 (igraph format)
#' @param net2 Network 2 (igraph format)
#' @param node_list Vector of names of selected nodes to plot in graph
#'
#' @import igraph
#' @export
get_combined_network <- function(net1, net2, node_list) {

  edges1 = as_edgelist(net1)
  edges2 = as_edgelist(net2)

  # index to keep edges involving selected nodes
  idx1 = which(edges1[,1] %in% node_list & edges1[,2] %in% node_list)
  idx2 = which(edges2[,1] %in% node_list & edges2[,2] %in% node_list)

  # create graphs using edges lists and get adjacency matrices
  g1 = graph_from_data_frame(edges1[idx1,], directed = FALSE, vertices = node_list)
  g2 = graph_from_data_frame(edges2[idx2,], directed = FALSE, vertices = node_list)

  # extract adjacency matrices and add to get edge weights indicating membership
  adj1 = as_adj(g1)
  adj2 = as_adj(g2)
  adj3 = adj1 + 2*adj2
  g3 = graph_from_adjacency_matrix(adj3, mode = "undirected", weighted = TRUE)

  print_shared_properties(adj3)

  return(g3)
}



#' Print shared node / edge details
#'
#' Print message showing number of unique / shared nodes / edges
#' @param adj igraph adjacency matrix, where edge weights 1, 2, 3 indicate unique edges to networks 1 and 2
#'   and shared edges respectively
#'
#' @export
print_shared_properties <- function(adj) {
  message('Number of edges:')
  message(length(which(adj@x == 1)), ' unique to network 1')
  message(length(which(adj@x == 2)), ' unique to network 2')
  message(length(which(adj@x == 3)), ' in both networks')

  # print message about node properties
  message('Number of nodes:')
  n1 = as.integer(apply(as.matrix(adj), 1, function(x) any(x %in% c(1,3))))
  n2 = as.integer(apply(as.matrix(adj), 1, function(x) any(x %in% c(2,3))))
  n3 = n1 + 2*n2
  message(length(which(n3 == 1)), ' unique to network 1')
  message(length(which(n3 == 2)), ' unique to network 2')
  message(length(which(n3 == 3)), ' in both networks')
}
