
#' @import DESeq2
#' @import ggplot2
#' @import RColorBrewer
#' @import pheatmap
#' @import ggthemes
#' @import ggsci
#' @import methods
#' @import MASS
#' @import plyr
#' @import PoiClaClu
#' @import reshape2
#' @import gridExtra
#' @importFrom SummarizedExperiment colData assay
#' @importFrom grDevices chull colorRampPalette dev.off pdf png
#' @importFrom stats cmdscale cutree dist hclust sd var
#' @importFrom utils head read.csv write.table
#' @importFrom ggdendro dendro_data segment label
#' @importFrom stats wilcox.test
.DEVis_env <- new.env(parent=emptyenv())



#' Initialize cutoff values for significance and fold-change filtering.
#'
#' This function initializes the p-value cutoffs and log2foldChange values for most visualizations,
#' filtering, and aggregation steps of this package. Significant differentially expressed genes will
#' be initially determined based on the provided p-value cutoff and filtering, when appropriate, will
#' generally be based on the provided log fold-change cutoff initialized by this function.
#' @keywords cutoff filter p-value fold-change significance threshold
#' @param p_signif P-value cutoff for determining DE genes. Default: 0.05
#' @param lfc_cut Log2 fold-change cutoff. Default: 1
#' @return This function does not return a value.
#' @export
#' @examples
#' \dontrun{
#'
#' init_cutoffs(p_signif=0.01, lfc_cut=1.5)
#'
#' }
init_cutoffs <- function(p_signif=0.05, lfc_cut=1)
{
  if(!is.numeric(p_signif))
  {
    stop("p-value cutoff must be numeric.")
    return(-1)
  }
  if(!is.numeric(lfc_cut))
  {
    stop("lfc_cut cutoff must be numeric.")
    return(-1)
  }

  assign("signif_cutoff", p_signif, envir=.DEVis_env)
  assign("lfc_cutoff", lfc_cut, envir=.DEVis_env)
}


#' Determine if visualizations are written to file, printed to screen, or both.
#'
#' This function initializes the output mode for the package to determine
#' if visualizations should be written to file, printed to screen, or both.
#' @keywords output file display
#' @param toggle Desired output mode. Determines how output plot data will be displayed. Options are "screen", "file", or "both".
#' @return This function does not return a value.
#' @export
#' @examples
#' \dontrun{
#'
#' #Don't save files, just print to the screen.
#' set_output_mode("screen")
#'
#' #Print to the screen and save files to appropriate locations.
#' set_output_mode("both")
#'
#' #Save files to appropriate locations, but do not print them to screen.
#' set_output_mode("file")
#'
#' }
set_output_mode <- function(toggle="both")
{
  if(typeof(toggle) != "character")
  {
    stop("toggle must be a string.")
    return(-1)
  }
  type_options = c("screen","file", "both")
  if(!(toggle %in% type_options))
  {
    stop('Output mode must be set to "screen", "file", or "both".')
    return(-1)
  }
  assign("output_mode", toggle, envir=.DEVis_env)
  print("Output mode set.")
}

