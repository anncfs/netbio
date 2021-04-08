## functions for plotting networks and subnetworks

#' Plot subnetwork comparison
#'
#' Plot subgraph (selected nodes) with edge colours indicating network membership
#'
#' @param g igraph network, where edge "weight" attribute is 1,2,3 indicating edges
#'   unique to network 1, unique to network 2, or shared by both networks respectively.
#'   In network plot these are coloured according to edge_col
#' @param l graph layout (default NULL, in which case layout_with_fr(g) used), matrix
#'   with two columns giving x,y coordinates for each node (in order of V(g))
#' @param edge_col vector of 3 edge colours for unique to 1, unique to 2, shared respectively
#'   default = c(rgb(0,0,1, alpha = 0.2), rgb(1,0,0,alpha = 0.2), rgb(0,0,0,alpha=0.2))
#'   i.e. blue, red, grey
#'
#' @export
#' @return returns layout
plot_subgraph_combined <- function(g, l = NULL,
                                     edge_col = c(rgb(0,0,1,alpha = 0.2),
                                                  rgb(1,0,0,alpha = 0.2),
                                                  rgb(0,0,0,alpha=0.2))) {

  if (is.null(l)) { l = layout_with_fr(g)}

  plot(g, layout = l,
       vertex.label = NA,
       vertex.size = 1,
       vertex.frame.color = NA,
       vertex.color = rgb(0,0,0,alpha = 0.4),
       edge.color = edge_col[E(g)$weight])

  return(l)
}


#' Plot subnetwork comparison by edge type
#'
#' Plot subgraph with each edge type in separate plot (cf plot_subgraph_combined)
#'
#' @param g igraph network, where edge "weight" attribute is 1,2,3 indicating edges
#'   unique to network 1, unique to network 2, or shared by both networks respectively.
#'   In network plot these are coloured according to edge_col
#' @param l graph layout (default NULL, in which case layout_with_fr(g) used), matrix
#'   with two columns giving x,y coordinates for each node (in order of V(g))
#' @param edge_col vector of 3 edge colours for unique to 1, unique to 2, shared respectively
#'   default = c(rgb(0,0,1, alpha = 0.2), rgb(1,0,0,alpha = 0.2), rgb(0,0,0,alpha=0.2))
#'   i.e. blue, red, grey
#'
#' @importFrom grDevices rgb
#' @export
#'
plot_subgraph_separated <- function(g, l = NULL,
                                    edge_col = c(rgb(0,0,1,alpha = 0.2),
                                                 rgb(1,0,0,alpha = 0.2),
                                                 rgb(0,0,0,alpha = 0.2))) {

  if (is.null(l)) { l = layout_with_fr(g)}

  edge_base = rep(rgb(1,1,1,alpha = 0), 3)

  for (i in 1:3) {
    edge_i = edge_base
    edge_i[i] = edge_col[i]
    plot(g, layout = l,
       vertex.label = NA,
       vertex.size = 1,
       vertex.frame.color = NA,
       vertex.color = rgb(0,0,0,alpha = 0.4),
       edge.color = edge_i[E(g)$weight])
  }

  return(l)
}



#' Plot degree and centrality by node
#'
#' Various plots of degree and centrality to compare (sub)networks
#'
#' @param net1 Network 1 (igraph format)
#' @param net2 Network 2 (igraph format)
#' @param node_list Vector of names of selected nodes to plot in graph
#' @param topX No. of top ranked nodes to display names for in graphs
#'
#' @export
plot_degree_centrality <- function(net1, net2, node_list, topX = 20) {

  # calculate degree / centrality for selected nodes in full graph (i.e. includes intra and inter edges)
  deg1 = degree(net1, v = node_list)
  deg2 = degree(net2, v = node_list)
  cen1 = betweenness(net1, v = node_list)
  cen2 = betweenness(net2, v = node_list)

  # plot comparisons
  plot_deg_v_cent(deg1, cen1, topX, node_list)
  title('network 1')
  plot_deg_v_cent(deg2, cen2, topX, node_list)
  title('network 2')
  plot_node_property(deg1,deg2)
  title('degree')
  plot_node_property(cen1,cen2)
  title('centrality')

}

#' Plot single node property comparison
#'
#' Plot single graph comparing node property (e.g. degree or centrality) values
#' between two networks
#'
#' @param p1 named vector of property in first network
#' @param p2 named vector of property in second network
#'
#' @importFrom graphics abline
plot_node_property <- function(p1, p2) {

  xymax = ceiling(max(c(p1,p2)))

  plot(p1, p2, pch = 16, cex = 0.5,
       xlab = 'network 1', ylab = 'network2',
       xlim = c(0,xymax), ylim = c(0,xymax))
  abline(0,1, col = 'red')

}

#' Plot single degree v centrality graph
#'
#' Plot single degree centrality graph for chosen network / node list
#'
#' @param degX named vector of node degrees
#' @param cenX named vector of node centrality
#' @param topX no. of top ranked nodes to select
#' @param node_list Vector of selected node names
#'
#' @importFrom graphics text
plot_deg_v_cent <- function(degX, cenX, topX, node_list) {

  idx = get_top_nodes(degX, cenX, topX)
  plot(degX, cenX, pch = 16, cex = 0.5,
       xlab = 'degree', ylab = 'centrality')
  text(degX[idx], cenX[idx], node_list[idx], pos = 2, cex = 0.7)
}

#' Get top nodes
#'
#' @param degX named vector of node degrees
#' @param cenX named vector of node centrality
#' @param topX no. of top ranked nodes to select
#'
#' @return index for top nodes in degX/cenX
get_top_nodes <- function(degX, cenX, topX) {

  sort_deg = sort(degX, index.return = T, decreasing = T)
  sort_cen = sort(cenX, index.return = T, decreasing = T)
  top_nodes = unique(c(sort_deg$ix[1:topX], sort_cen$ix[1:topX]))
  return(top_nodes)
}
