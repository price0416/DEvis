#' An S4 class extending the DESeqResults object containing additional fields for differential expression data.
#'
#' This class is used to separate and identify relevant pieces of data by various visualizations in this package.
#'
#' @slot numDE The total number of differentially expressed genes in the result set.
#' @slot numUp The number of upregulated differentially expressed genes in the result set.
#' @slot numDown The number of downregulated differentially expressed genes in the result set.
#' @slot allNames The names of all differentially expressed genes in the result set.
#' @slot upNames The names of all upregulated differentially expressed genes in the result set.
#' @slot downNames The names of all downregulated differentially expressed genes in the result set.
#' @slot case The "case" condition that was used to make the contrast for this result set.
#' @slot control The "control" condition that was used to make the contrast for this result set.
#' @slot contrast The contrast that was made for this result set.  I.E. case_vs_control
#' @slot design_field The field in the provided target data that was used for experimental design in differential expression.
#' I.E. Condition_Time.
#' @importClassesFrom DESeq2 DESeqResults
#' @export
setClass("DESeqResMeta",
         slots = list(
           numDE = "numeric",
           numUp = "numeric",
           numDown = "numeric",
           allNames = "character",
           upNames = "character",
           downNames = "character",
           case = "character",
           control = "character",
           contrast = "character",
           design_field = "character"
         ),

         contains = "DESeqResults"
)
