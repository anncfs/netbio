% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_networks.R
\name{plot_subgraph_combined}
\alias{plot_subgraph_combined}
\title{Plot subnetwork comparison}
\usage{
plot_subgraph_combined(
  g,
  l = NULL,
  edge_col = c(rgb(0, 0, 1, alpha = 0.2), rgb(1, 0, 0, alpha = 0.2), rgb(0, 0, 0, alpha
    = 0.2))
)
}
\arguments{
\item{g}{igraph network, where edge "weight" attribute is 1,2,3 indicating edges
unique to network 1, unique to network 2, or shared by both networks respectively.
In network plot these are coloured according to edge_col}

\item{l}{graph layout (default NULL, in which case layout_with_fr(g) used), matrix
with two columns giving x,y coordinates for each node (in order of V(g))}

\item{edge_col}{vector of 3 edge colours for unique to 1, unique to 2, shared respectively
default = c(rgb(0,0,1, alpha = 0.2), rgb(1,0,0,alpha = 0.2), rgb(0,0,0,alpha=0.2))
i.e. blue, red, grey}
}
\value{
returns layout
}
\description{
Plot subgraph (selected nodes) with edge colours indicating network membership
}
