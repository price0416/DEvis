
#' Read tab or comma delimited count data.
#'
#' This function takes a count data file and returns the formatted data. Requires that init_data_paths() has been run.
#' @param count_input Count data in tab-delimited or comma-delimited (CSV) format, assumes headers exist for each column.
#' Note that column names in the input count file must match row names in the input target metadata file.
#' @param delim Indicate if file is comma or tab delimited.  "t" indicated tab delimited and "c" indicates comma separated data.  Default: "t"
#' @seealso \code{\link{prep_targets}}, \code{\link{prep_dds_from_data}}, \code{\link{init_data_paths}}
#' @return This function will return a properly formatted count data table based on a provided input file.
#' @keywords counts input
#' @export
#' @examples
#' \dontrun{
#'
#' #Read a tab-delimited text file as count input.
#' myCounts <- prep_counts(count_input="master_count_data.txt", delim="t")
#'
#' }
prep_counts <- function(count_input, delim="t")
{
  if(delim != "t" && delim != "c")
  {
    print("Error: Possible delimiter values are 't' (tab-delimited), or 'c' (comma-separated).")
    return(-1)
  }
  if(!exists("counts_dir", envir=.DEVis_env))
  {
    stop("Count directory not initialized.  Run init_data_paths() first.")
    return(-1)
  }

  #Read count matrix.
  print("Reading count file...")
  setwd(.DEVis_env$counts_dir)
  if(delim == "t")
  {
    counts <- read.csv(count_input, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  }
  if(delim == "c")
  {
    counts <- read.csv(count_input, header=TRUE, sep=",", stringsAsFactors=FALSE)
  }
  print(.DEVis_env$working_dir)
  setwd(.DEVis_env$working_dir)
  rownames(counts) <- counts$X
  counts <- counts[,c(2:ncol(counts))]
  assign("countDat", counts, envir=.DEVis_env)
  return(counts)
}


#' Read tab or comma delimited target metadata file.
#'
#' This function takes a targets metadata file and returns formatted data. Requires that init_data_paths() has been run.
#' @param target_input Target data in tab separated format, assumes headers exist for each column.
#' Note that column names in the input count file must match row names in the input target metadata file.
#' @param delim Indicate if file is comma or tab delimited.  "t" indicated tab delimited and "c" indicates comma separated data.  Default: "t"
#' @seealso \code{\link{prep_counts}}, \code{\link{prep_dds_from_data}}, \code{\link{init_data_paths}}
#' @return This function will return a properly formatted targets data table based on a provided input file.
#' @keywords targets metadata input
#' @export
#' @examples
#' \dontrun{
#'
#' #Read a tab-delimited text file as target metadata input.
#' myTargets <- prep_targets(target_input="master_count_data.txt", delim="t")
#'
#' }
prep_targets <- function(target_input, delim="t")
{
  if(delim != "t" && delim != "c")
  {
    print("Error: Possible delimiter values are 't' (tab-delimited), or 'c' (comma-separated).")
    return(-1)
  }
  if(!exists("targets_dir", envir=.DEVis_env))
  {
    stop("Targets directory not initialized.  Run init_data_paths() first.")
    return(-1)
  }

  #Read target matrix.
  print("Reading targets file...")
  setwd(.DEVis_env$targets_dir)
  if(delim == "t")
  {
    targets <- read.csv(target_input, header=TRUE, sep="\t")
  }
  if(delim == "c")
  {
    targets <- read.csv(target_input, header=TRUE, sep=",")
  }
  setwd(.DEVis_env$working_dir)
  rownames(targets) <- targets$X
  targets <- targets[,c(2:ncol(targets))]
  targets$targetID <- rownames(targets)

  assign("tgt_dat", targets, envir=.DEVis_env)

  return(targets)
}


#' Prepare a DESeq2 object based on count and target data.
#'
#' This function takes count and targets data and creates a DESeq2 object based on the provided design.
#' @param count_input Count data in dataframe format as read with prep_counts().  Column names must correspond to target rownames.
#' @param target_input Target data in dataframe format as read with prep_targets().  Row names must correspond to count column names.
#' @param experiment_design Experimental design for DE comparison.  Can include multiple factors.  Design factors should correspond to column names
#' in the provided targets file.  Include the primary factor as the final factor in the design and any secondary or batch effects prior.
#' Must be prefaced with a tilde (~).   I.E. ~ Batch1 + Batch2 + SecondaryCondition + PrimaryCondition
#' @param collapseReps Boolean. Collapse technical replicates. Default: False
#' @param rep_field_vector Target metadata column identifying which samples correspond to which replicates
#' and maps samples to their replicate identifiers.  Required if collapseReps is TRUE.
#' @param stabilization Method of normalizing transformation to perform.  Possible values are "rld" and "vst".  rld will execute
#' DESeq2's rlog transformation, whereas vst will apply variance stabilizing transformation to the data.  For larger data sets
#' vst is recommended.
#' @return This function will return a properly formatted DESeq2 object based on the provided normalization method and experimental design.
#' @seealso \code{\link{prep_counts}}, \code{\link{prep_targets}}, \code{\link{init_cutoffs}}
#' @keywords dds design replicates normalization
#' @export
#' @examples
#' \dontrun{
#'
#' /*
#'  * Create a DESeq2 data object using previously read count and target data.
#'  * Differential expression groups indicated by "Infection" column of target data.
#'  * Apply a rlog transformation to the data. TMM normalization applied.
#'  */
#' dds <- prep_dds_from_data(count_input=myCounts, target_input=myTargets,
#'                            experiment_design= ~ Infection, stabilization="rld")
#'
#'
#' /*
#'  * Create a DESeq2 data object using previously read count and target data.
#'  * Differential expression groups indicated by "Infection" column of target data.
#'  * Takes into account a batch effect column of target data.
#'  * Apply a variance stabilizing transformation to the data. TMM normalization applied.
#'  */
#' dds <- prep_dds_from_data(count_input=myCounts, target_input=myTargets,
#'                            experiment_design= ~ batch + Infection, stabilization="vst")
#'
#'
#' /*
#'  * Create a DESeq2 data object using previously read count and target data.
#'  * Differential expression groups indicated by "Time" and "Infection" columns of target data.
#'  * Takes into account a batch effect column of target data.
#'  * Apply a variance stabilizing transformation to the data. TMM normalization applied.
#'  * Collapse technical replicates based on the "replicate" column of target data.
#'  */
#' dds <- prep_dds_from_data(count_input=myCounts, target_input=myTargets,
#'                            experiment_design= ~ batch + Time + Infection,
#'                            collapseReps=TRUE,
#'                            rep_field_vector="replicate",
#'                            stabilization="vst")
#'
#'  }
prep_dds_from_data <- function(count_input, target_input, experiment_design, collapseReps=FALSE, rep_field_vector="Replicate", stabilization="rld")
{
  #Confirm parameters.
  if(stabilization != "rld" && stabilization != "vst")
  {
    stop("Possibile values for stabilization are \"rld\" and \"vst\".  Use vst for larger data sets.")
    return(-1)
  }

  #Create DESeq Object and specify design parameters.
  print("Creating DESeq Object...")
  dds <- DESeqDataSetFromMatrix(countData = count_input, colData = target_input, design = experiment_design)

  if(collapseReps == TRUE)
  {
    print("Collapsing replicates...")
    sub1 <- paste("collapseReplicates(dds,target_input$", rep_field_vector, ")", sep="")
    dds <- eval(parse(text = sub1))

    #Store the column used for replicate collapse.
    collapse_column <- deparse(substitute(rep_field_vector))
    isRepCollapsed <- substring(sub("^[^$]*", "", collapse_column),2)

  }
  else
  {
    isRepCollapsed <- "NA"
  }

  #TMM normalization.
  print("Performing TMM normalization...")
  dds <- estimateSizeFactors(dds)

  #Do the necessary data transformation.
  if(stabilization == "rld")
  {
    print("Performing rLog transform...")
    assign("stabilized_data", rlog(dds), envir=.DEVis_env)
  }
  if(stabilization == "vst")
  {
    print("Performing variance stabilization.")
    assign("stabilized_data", varianceStabilizingTransformation(dds), envir=.DEVis_env)
  }

  print("DESeq Object prepared.")
  assign("dds_glbl", dds, envir=.DEVis_env)
  assign("norm_counts", counts(dds,normalized=TRUE), envir=.DEVis_env)

  return(dds)

}
