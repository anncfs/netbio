#' Loads PIDC data
#'
#' Loads node and edge data from PIDC output file
#'
#' Reads in the PIDC output file, removes duplicate edges (network is
#' undirected so each edge is represented twice in consecutive rows).
#'
#' @param filename Path and filename for PIDC .txt output file
#'                 (each row is in format: nodeX nodeY edgescore)
#'
#' @import data.table
#'
#' @return A data.table with row per edge, columns giving nodeX nodeY score
#' @export
load_pidc_data <- function(filename){

  # import data and remove duplicate edges
  duplicate_edges = data.table::fread(filename)
  edges = duplicate_edges[seq(2, dim(duplicate_edges)[1], 2)]

  print(paste0("Imported network: ", dim(edges)[1], " nodes"))
  return(edges)
}



#' Threshold network
#'
#' Creates thresholded network by selecting top scoring edges
#'
#' @param edges A data.frame or matrix with row per edge, columns giving nodeX nodeY score
#' @param threshold Threshold, either a % of top scoring edges to include (default) or a number of edges
#' @param nodes Reference list of node names, if not provided (default = NULL), the nodes in the threshold network
#'   are used
#' @param as_percentage Logical value, TRUE (default) means \code{threshold} is given as a %, FALSE indicates
#'   \code{threshold} is the number of edges to keep
#'
#' @import igraph
#'
#' @return Thresholded network in iGraph format
#' @export
threshold_network <- function(edges, threshold, as_percentage = TRUE, nodes = NULL){

  # calculate no. of edges to keep
  if (as_percentage) {
    n_edges = ceiling(dim(edges)[1] * threshold/100)
  } else {
    n_edges = threshold
  }

  # ensure edge list sorted by decreasing edge score
  edges = edges[order(-edges[,3]),]

  # get list of nodes
  if (is.null(nodes)) {
    nodes = edges[1:n_edges,1]
  }

  # create adjacency matrix and add edges
  adj = matrix(0, length(nodes), length(nodes))
  for (e in 1:n_edges){
    node1 = which(nodes == edges[[e,1]])
    node2 = which(nodes == edges[[e,2]])
    adj[node1,node2] = 1
    adj[node2,node1] = 1
  }

  igraph::graph_from_adjacency_matrix(adj, mode = 'undirected')
}


