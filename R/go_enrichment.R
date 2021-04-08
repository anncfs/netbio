## Functions for GO term analysis

#' GO term enrichment of network communities
#'
#' Use WebGestaltR package to do GO term over-representation analysis for each community in network
#'
#' @param net_clusters Named list of community membership
#' @param min_size minimum size of cluster to do GO enrichment for
#' @param ref_genes list of genes to use as reference for GO enrichment (default = NULL, in which case
#'   full genome used)
#' @param ontology which ontology to use (any valid WebGestalt enrichDatabase,
#'   default = 'geneontology_Biological_Process')
#' @param filename filename to record summarised results
#' @param dir_name directory for WebGestaltR to output results (default = NULL, detailed output suppressed)
#'
#' @export
#'
go_enrichment <- function(net_clusters, min_size = 10,
                          ontology = 'geneontology_Biological_Process',
                          filename = NULL,
                          ref_genes = NULL,
                          dir_name = NULL) {

  # check output
  if (is.null(dir_name)) {
    isOutput = FALSE
    dir_name = getwd()
  } else {
    isOutput = TRUE
  }

  # create directory if doesn't exist
  if (!file.exists(dir_name)) {dir.create(dir_name)}

  # get list of cluster IDs above min size
  idx = which(table(net_clusters) >= min_size)

  for (i in idx) {

    int_genes = names(which(net_clusters == i))

    enrichResult = WebGestaltR::WebGestaltR(enrichMethod = "ORA",
                                            organism = "mmusculus", enrichDatabase = ontology,
                                            interestGene = int_genes, interestGeneType = "genesymbol",
                                            referenceGene = ref_genes, referenceGeneType = "genesymbol",
                                            referenceSet = "genome",
                                            isOutput = isOutput, outputDirectory = dir_name,
                                            projectName = paste0("cluster",i))

    # (if no sig terms by fdr, then WebGestaltR returns NULL)
    if (!is.null(enrichResult)) {
      write(paste0("Cluster ", i, " (top FDR):"), file = filename, append = T)
      gdata::write.fwf(enrichResult[1:10, c("geneSet", "description", "overlap", "enrichmentRatio", "pValue", "FDR")],
                     file = filename, rownames = F, append = T)
      write(paste0("\nCluster ", i, " (top enrichmentRatio):"), file = filename, append = T)
      temp = enrichResult[order(-enrichResult$enrichmentRatio),]
      gdata::write.fwf(temp[1:10, c("geneSet", "description", "overlap", "enrichmentRatio", "pValue", "FDR")],
                     file = filename, rownames = F, append = T)
      write('\n\n', file = filename, append = T)
    }
  }
}

