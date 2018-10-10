#' Incorporate additional data about differentially expressed genes into a DESeq2 result set.
#'
#' Given a DESeq result set, as generated with DESeq2::results(), this function adds several additional components to the object.
#' Variables created are numDE, numUp, numDown, allNames, upNames, and downNames, allDE, upDE, downDE, case, control, contrast, and designField.
#' All calculations are performed based on significance thresholds set using the init_cutoffs() function.
#' This function is typically used as a helper function for data aggregation.
#' @param res A DESeq2 result set as generated through DESeq2::results().  Must be a DESeq object of type S4.
#' @keywords DE aggregate
#' @return An expanded DESeq result set of type DESeqResMeta containing the following fields:
#' [numDE, numUp, numDown, allNames, upNames, and downNames, allDE, upDE, downDE, case, control, contrast, and designField]
#' @seealso \code{\link{init_cutoffs}}
#' @export
#' @examples
#' \dontrun{
#'
#' #Enrich a result set. This will make several new data points available.
#' res.day1 <- results(dds, contrast=c("Condition_Time", "day1_disease", "day1_control"))
#' enriched_result <- enrich_res(res.day1)
#'
#' #Print number of differentially expressed genes.
#' print(enriched_result@numDE)
#'
#' #Get the names of all up-regulated differentially expressed genes.
#' upreg_de_genes <- enriched_result@upNames
#' }
enrich_res <- function(res)
{
  #Validate parameters.
  try(if(typeof(res) != "S4") stop("Error: enrich_res requires S4 type object as result set for parameter res.", call. = FALSE))

  if(!exists("signif_cutoff", envir=.DEVis_env))
  {
    stop("Error: enrich_res requires that significance cutoffs be established first.  Please run init_cutoffs() first.")
  }

  #Create a class that inherits DESeqResults class, and contains extra slots for metadata we want for data visualization and aggregation.
  newRes = new("DESeqResMeta")

  #Populate the new object with all the values from the initial result set.
  newRes@rownames <- res@rownames
  newRes@nrows <- res@nrows
  newRes@listData <- res@listData
  newRes@elementType <- res@elementType
  newRes@elementMetadata <- res@elementMetadata
  newRes@metadata <- res@metadata
  newRes$baseMean <- res$baseMean
  newRes$log2FoldChange <- res$log2FoldChange
  newRes$lfcSE <- res$lfcSE
  newRes$stat <- res$lfcSE
  newRes$pvalue <- res$pvalue
  newRes$padj <- res$padj


  #See if the result set has any DE genes at all.
  if(!is.na(table(res$padj < .DEVis_env$signif_cutoff)[2]))
  {
    numDE = table(res$padj < .DEVis_env$signif_cutoff)[2]
  }
  else
  {
    numDE = 0
  }
  newRes@numDE = numDE

  #Extract all DE names and store.
  DE_names <- rownames(head(res[order(res$padj),], n=newRes@numDE))
  newRes@allNames <- DE_names

  #Extract number upregulated and downregulated, then store.
  if(!is.na(table(res$padj < .DEVis_env$signif_cutoff & res$log2FoldChange < 0)[2]))
  {
    numDown = table(res$padj < .DEVis_env$signif_cutoff & res$log2FoldChange < 0)[2]
  }
  else
  {
    numDown = 0
  }
  if(!is.na(table(res$padj < .DEVis_env$signif_cutoff & res$log2FoldChange > 0)[2]))
  {
    numUp = table(res$padj < .DEVis_env$signif_cutoff & res$log2FoldChange > 0)[2]
  }
  else
  {
    numUp = 0
  }
  newRes@numDown <- numDown
  newRes@numUp <- numUp

  #Extract and store the up and downregulard DE gene names.
  downNames = rownames(res)[which(res$padj < .DEVis_env$signif_cutoff & res$log2FoldChange < 0)]
  upNames   = rownames(res)[which(res$padj < .DEVis_env$signif_cutoff & res$log2FoldChange > 0)]
  newRes@downNames <- downNames
  newRes@upNames <- upNames

  #Attach metadata about comparisons.
  case = unlist(strsplit(res@elementMetadata$description[2], split=" "))[6]
  control = unlist(strsplit(res@elementMetadata$description[2], split=" "))[8]
  design_field = unlist(strsplit(res@elementMetadata$description[2], split=" "))[5]
  contrast = paste(case,"_vs_",control,sep="")

  newRes@case <- case
  newRes@control <- control
  newRes@design_field <- design_field
  newRes@contrast <- contrast

  return(newRes)

}


#' Aggregate and retrieve data from multiple differentially expressed result sets.
#'
#' This function finds the union or intersection of DE gene names for all result sets provided, extracts those rows
#' from all DE result sets, and merges the values into a single data set containing log2 fold-change or padj values.
#' @param res_list A list of DESeq result sets created with DESeq2::results(). I.E: list(res1, res2, ..., resN).
#' @param method Method for aggregating data.  "union" will gather expression data for all samples in results list
#' when a gene is differentially expressed in at least one sample.  "intersection" will gather expression data only for genes that
#' are DE in all samples of the result list.
#' @param type Type of data to retrieve.  Options are "lfc" (log2foldchange) and "padj" (adjusted p-value).
#' @param lfc_filter Filter genes based on a minimum log fold-change as initialized in init_cutoffs(). Boolean. Default=FALSE.
#' @return Returns the aggregated data frame containing fold-change or p-values for the provided result sets.
#' @seealso \code{\link{init_cutoffs}}
#' @keywords DE aggregate filter
#' @export
#' @examples
#' \dontrun{
#'
#' #Prepare a result list for aggregation.
#' res.day1 <- results(dds, contrast=c("Condition_Time", "day1_disease", "day1_control"))
#' res.day2 <- results(dds, contrast=c("Condition_Time", "day2_disease", "day2_control"))
#' res.day3 <- results(dds, contrast=c("Condition_Time", "day3_disease", "day3_control"))
#' myResList <- list(res.day1, res.day2, res.day3)
#'
#'
#' /*
#'  * Data for all result sets will be included for each gene if that gene was found to be
#'  * differentially expressed in at least one of the provided result sets.
#'  * Filter based on a minimum fold-change.
#'  */
#' aggregated_lfc_union  <- get_de_data(res_list=myResList, method="union",
#'                                       type="lfc", lfc_filter=TRUE)
#' aggregated_pval_union <- get_de_data(res_list=myResList, method="union",
#'                                       type="padj", lfc_filter=TRUE)
#'
#'
#' /*
#'  * Data for all result sets will be included for each gene only if that gene was found to
#'  * be differentially expressed in all provided result sets. Do not apply a fold-change filter.
#'  * Significance is determined only by p-value threshold.
#'  */
#' aggregated_lfc_intersect  <- get_de_data(res_list=myResList, method="intersection",
#'                                           type="lfc", lfc_filter=FALSE)
#' aggregated_pval_intersect <- get_de_data(res_list=myResList, method="intersection",
#'                                           type="padj", lfc_filter=FALSE)
#'
#' }
get_de_data <- function(res_list, method="union", type="lfc", lfc_filter=FALSE)
{
  #Parameter validation.
  if((typeof(res_list) != "list") || (typeof(res_list[[1]]) != "S4") || (length(res_list) <= 0))
  {
    stop("init_master_res(res_list) requires list type object containing DESeq result sets.")
    return(-1)
  }
  type_options = c("lfc","pval")
  if(!(type %in% type_options))
  {
    stop('type must be one of the following: "lfc","pval"')
    return(-1)
  }
  if(!(exists("lfc_cutoff", envir=.DEVis_env)))
  {
    stop("Cutoff values have not been initialized.  Run init_cutoffs() first.")
    return(-1)
  }
  if(typeof(lfc_filter) != "logical")
  {
    stop("lfc_filter must be a logical value. TRUE or FALSE. LFC cutoff can be set using the init_cutoffs() function.")
    return(-1)
  }
  method_options = c("union","intersection")
  if(!(method %in% method_options))
  {
    stop('method must be one of the following: "union","intersection"')
    return(-1)
  }

  #Enrich the data.
  enriched <- unlist(lapply(res_list, enrich_res))
  enrich <- unlist(enriched)

  #Get the union sets.
  res.names <- lapply(enrich, function(x) x@allNames)
  if(method=="union")
  {
    union.names = Reduce(union, res.names)
    master_union_names <- union.names
    master_union <- lapply(enrich, function(x) x[union.names,])
    master_lfc   <- lapply(master_union, function(x) x$log2FoldChange)
    master_colnames_lfc <- lapply(master_union, function(x) x@contrast)
    master_colnames <- master_colnames_lfc
    master_padj  <- lapply(master_union, function(x) x$padj)
    master_colnames_padj <- lapply(master_union, function(x) paste(x@contrast, ".padj", sep=""))
  }
  if(method=="intersection")
  {
    intersect.names = Reduce(intersect, res.names)
    master_intersect_names <- intersect.names
    master_intersect <- lapply(enrich, function(x) x[intersect.names,])
    master_lfc   <- lapply(master_intersect, function(x) x$log2FoldChange)
    master_colnames_lfc <- lapply(master_intersect, function(x) x@contrast)
    master_colnames <- master_colnames_lfc
    master_padj  <- lapply(master_intersect, function(x) x$padj)
    master_colnames_padj <- lapply(master_intersect, function(x) paste(x@contrast, ".padj", sep=""))
  }

  #Finalize the master data frame, returning data for p-value or log fold-change where appropriate.
  if(type == "lfc")
  {
    master_df <- data.frame(master_lfc)
    master_colnames <- master_colnames_lfc
  }
  if(type == "pval")
  {
    master_df <- data.frame(master_padj)
    master_colnames <- master_colnames_padj
  }
  colnames(master_df) <- master_colnames
  if(method=="union")
  {
    rownames(master_df) <- rownames(master_union[[1]])
  }
  if(method=="intersection")
  {
    rownames(master_df) <- rownames(master_intersect[[1]])
  }

  #Handle log fold-change filtering if necessary.
  if(lfc_filter == TRUE & type == "lfc")
  {
    lfc_df <- data.frame(master_lfc)
    colnames(lfc_df) <- master_colnames_lfc
    rownames(lfc_df) <- rownames(master_union[[1]])
    lfc_keep <- lfc_df[rowSums(abs(lfc_df) > abs(.DEVis_env$lfc_cutoff)) >= 1, ]
    lfc_keep_genes <- rownames(lfc_keep)

    master_df$gene <- rownames(master_df)

    master_df_trim <- master_df[master_df$gene %in% lfc_keep_genes,]
    master_df <- master_df_trim
    master_df$gene <- NULL
  }

  return(master_df)
}


#' Write differentially expressed gene data for multiple result sets to file.
#'
#' This function accepts a list of DE result sets and writes their outputs individually to file.
#' Relies on init_cutoffs() significance thresholds. Filenames will be generated automatically based on the contrasts performed
#' during the results() function. Output will be written to the /DE/data/ folder.
#' @param res_list A list of DESeq result sets.
#' @param lfc_filter Also impose a filter on the DE set by a minimum absolute log2foldChange as determined by init_cutoffs(). Default=FALSE.
#' @return Output for each result set will be written to file This function does not return a value.
#' @keywords DE output result contrast
#' @seealso \code{\link{init_cutoffs}}, \code{\link{create_dir_struct}}
#' @export
#' @examples
#' \dontrun{
#'
#' #Prepare a result list for aggregation.
#' res.day1 <- results(dds, contrast=c("Condition_Time", "day1_disease", "day1_control"))
#' res.day2 <- results(dds, contrast=c("Condition_Time", "day2_disease", "day2_control"))
#' res.day3 <- results(dds, contrast=c("Condition_Time", "day3_disease", "day3_control"))
#' myResList <- list(res.day1, res.day2, res.day3)
#'
#' #Write differentally expressed gene data for each contrast to file.
#' #Include a minimum fold change filter. This will output 3 files.
#' write_all_de_results(res_list=myResList, lfc_filter=TRUE)
#'
#' }
write_all_de_results <- function(res_list, lfc_filter=FALSE)
{
  #Parameter validation.
  if((typeof(res_list) != "list") || (typeof(res_list[[1]]) != "S4") || (length(res_list) <= 0))
  {
    stop("write_all_de_results(res_list) requires list type object containing DESeq result sets.")
    return(-1)
  }
  if(typeof(lfc_filter) != "logical")
  {
    stop("lfc_filter must be a logical value. TRUE or FALSE. LFC cutoff can be set using the init_cutoffs() function.")
    return(-1)
  }
  if(!exists("DE_dir", envir=.DEVis_env))
  {
    stop("Directories not initialized.  Try running create_dir_struct().")
    return(-1)
  }

  enriched <- unlist(lapply(res_list, enrich_res))
  enrich <- unlist(enriched)

  #Extract appropriate genes and write DE gene files for each contrast to file.
  setwd(.DEVis_env$DE_dir)
  for (i in 1:length(enrich))
  {
    curContrast <- enrich[[i]]@contrast
    outfile_name <- paste(curContrast, "_DEgenes.txt", sep="")

    #Remove NA values from padj (due to cooks cutoff).  Set to 1, these genes are not DE anyway.
    res_list[[i]]$padj <- ifelse(is.na(res_list[[i]]$padj), 1, res_list[[i]]$padj)

    #Remove NA values from log2Fold change column, set to 0 if NA.
    if(lfc_filter == TRUE)
    {
      res_list[[i]]$log2FoldChange <- ifelse(is.na(res_list[[i]]$log2FoldChange), 0, res_list[[i]]$log2FoldChange)
      cur_DE <- res_list[[i]][which(res_list[[i]]$padj < .DEVis_env$signif_cutoff & abs(res_list[[i]]$log2FoldChange) > .DEVis_env$lfc_cutoff),]
    }
    else
    {
      cur_DE <- res_list[[i]][which(res_list[[i]]$padj < .DEVis_env$signif_cutoff),]
    }

    #Write the current DE set to file.
    write.table(cur_DE, outfile_name, sep="\t", row.names=TRUE, col.names=NA)
  }
  print("Results written to DE/data/.")
  setwd(.DEVis_env$working_dir)
}


#' Create a data set consisting of aggregated data for multiple contrasts.
#'
#' This function creates a master result set for the provided DE result sets.  This function finds the union of
#' DE gene names and extracts those rows from all DE result sets, then merges the sets into a single
#' master DE file containing both log2foldChange and padj values.  The results are written to the DE
#' output directory and returned by the function.
#' Relies on init_cutoffs() significant thresholds.
#' @param res_list A list of DESeq result sets created with DESeq2::results(). I.E: list(res1, res2, ..., resN).
#' @param filename Destination output filename.
#' @param method Method for aggregating data.  "union" will gather expression data for all samples in results list
#' when a gene is differentially expressed in at least one sample.  "intersection" will gather expression data only for genes that
#' are DE in all samples of the result list.
#' @param lfc_filter Impose a filter on the DE set by a minimum absolute log2foldChange as determined by init_cutoffs(). Default=FALSE.
#' @return This function returns a data frame containing union-based aggregation of all provided result sets.
#' Results of this function are also written to file in the /DE/data/ directory in tab-delimited format.
#' @seealso \code{\link{init_cutoffs}}, \code{\link{create_dir_struct}}
#' @keywords DE aggregate result master filter
#' @export
#' @examples
#' \dontrun{
#'
#' #Prepare a result list for aggregation.
#' res.day1 <- results(dds, contrast=c("Condition_Time", "day1_disease", "day1_control"))
#' res.day2 <- results(dds, contrast=c("Condition_Time", "day2_disease", "day2_control"))
#' res.day3 <- results(dds, contrast=c("Condition_Time", "day3_disease", "day3_control"))
#' myResList <- list(res.day1, res.day2, res.day3)
#'
#' #Generate the master DE data frame using union-based aggregation
#' master_df <- create_master_res(res_list=myResList, filename="master_DE.txt",
#'                                method="union")
#'
#' #Generate the master DE data frame using intersection-based aggregation.
#' #Filter genes below minimum log fold-change.
#' master_df <- create_master_res(res_list=myResList, filename="master_DE.txt",
#'                                method="union", lfc_filter=TRUE)
#'
#' }
create_master_res <- function(res_list, filename="master_DE.txt", method="union", lfc_filter=FALSE)
{
  #Parameter validation.
  if((typeof(res_list) != "list") || (typeof(res_list[[1]]) != "S4") || (length(res_list) <= 0))
  {
    stop("create_master_res() requires list type object containing DESeq result sets.")
    return(-1)
  }
  if(!exists("signif_cutoff", envir=.DEVis_env))
  {
    stop("create_master_res() requires cutoffs to be initialized.  Run init_cutoffs() first.")
    return(-1)
  }
  if(typeof(filename) != "character")
  {
    stop("filename must be a string.")
    return(-1)
  }
  if(length(filename) == 0)
  {
    stop("filename cannot be empty.")
    return(-1)
  }
  method_options = c("union","intersection")
  if(!(method %in% method_options))
  {
    stop('method must be one of the following: "union","intersection"')
    return(-1)
  }
  if(typeof(lfc_filter) != "logical")
  {
    stop("lfc_filter must be a logical value. TRUE or FALSE. LFC cutoff can be set using the init_cutoffs() function.")
    return(-1)
  }

  #Enrich the data.
  enriched <- unlist(lapply(res_list, enrich_res))
  enrich <- unlist(enriched)

  #Union-based aggregation.
  if(method=="union")
  {
    res.names <- lapply(enrich, function(x) x@allNames)
    union.names = Reduce(union, res.names)
    master_union_names <- union.names
    assign("master_union_names", union.names, envir=.DEVis_env)
    assign("aggregate_names", union.names, envir=.DEVis_env)
    master_union <- lapply(enrich, function(x) x[union.names,])
    master_lfc   <- lapply(master_union, function(x) x$log2FoldChange)
    master_padj  <- lapply(master_union, function(x) x$padj)
    master_colnames_lfc <- lapply(master_union, function(x) paste(x@contrast, ".lfc", sep=""))
    master_colnames_padj <- lapply(master_union, function(x) paste(x@contrast, ".padj", sep=""))
    master_colnames <- c(master_colnames_lfc, master_colnames_padj)

    #Create the master data frame.
    master_df <- data.frame(master_lfc, master_padj)
    colnames(master_df) <- master_colnames
    rownames(master_df) <- rownames(master_union[[1]])

    #Filter by LFC if requested.
    if(lfc_filter==TRUE)
    {
      lfc_df <- data.frame(master_lfc)
      colnames(lfc_df) <- master_colnames_lfc
      rownames(lfc_df) <- rownames(master_union[[1]])
      lfc_keep <- lfc_df[rowSums(abs(lfc_df) > abs(.DEVis_env$lfc_cutoff)) >= 1, ]
      lfc_keep_genes <- rownames(lfc_keep)

      master_df$gene <- rownames(master_df)

      master_df_trim <- master_df[master_df$gene %in% lfc_keep_genes,]
      master_df <- master_df_trim
      master_df$gene <- NULL
    }
  }
  #Intersect-based aggregation.
  if(method=="intersection")
  {
    res.names <- lapply(enrich, function(x) x@allNames)
    intersect.names = Reduce(intersect, res.names)
    master_intersect_names <- intersect.names
    assign("master_intersect_names", intersect.names, envir=.DEVis_env)
    assign("aggregate_names", intersect.names, envir=.DEVis_env)
    master_intersect <- lapply(enrich, function(x) x[intersect.names,])
    master_lfc   <- lapply(master_intersect, function(x) x$log2FoldChange)
    master_colnames_lfc <- lapply(master_intersect, function(x) x@contrast)
    master_padj  <- lapply(master_intersect, function(x) x$padj)
    master_colnames_padj <- lapply(master_intersect, function(x) paste(x@contrast, ".padj", sep=""))
    master_colnames <- c(master_colnames_lfc, master_colnames_padj)

    #Create the master data frame.
    master_df <- data.frame(master_lfc, master_padj)
    colnames(master_df) <- master_colnames
    rownames(master_df) <- rownames(master_intersect[[1]])

    if(lfc_filter==TRUE)
    {
      lfc_df <- data.frame(master_lfc)
      colnames(lfc_df) <- master_colnames_lfc
      rownames(lfc_df) <- rownames(master_intersect[[1]])
      lfc_keep <- lfc_df[rowSums(abs(lfc_df) > abs(.DEVis_env$lfc_cutoff)) >= 1, ]
      lfc_keep_genes <- rownames(lfc_keep)

      master_df$gene <- rownames(master_df)

      master_df_trim <- master_df[master_df$gene %in% lfc_keep_genes,]
      master_df <- master_df_trim
      master_df$gene <- NULL
    }
  }

  #Write the data to file.
  setwd(.DEVis_env$DE_dir)
  write.table(master_df, filename, sep="\t", row.names=TRUE, col.names=NA)
  setwd(.DEVis_env$working_dir)

  return(master_df)
}


#' Rename gene IDs based on a 1-to-1 mapping file.
#'
#' This function accepts a 1-to-1 mapping file of gene ids and updates count data to reflect the new IDs.  The mapping file should contain
#' two columns "Gene stable ID" and "Gene name".  The "Gene stable id" should correspond to gene ids in your count data file, and "gene name"
#' will be used to replace the existing genes.  This function has can be used to translate difficult to understand gene identifiers, such as
#' ensembl IDs to to more human readable gene names.  A mapping file can be generated at \url{http://www.ensembl.org/biomart/martview/}.
#' This function should be run using count data read by prep_counts().  The returned data can then be used in prep_dds_from_data().
#' @param count_data Count data as read using the prep_counts() function.
#' @param meta_file File containing 1-to-1 gene id mapping. Tab-delimited format expected.
#' Column 1 should be titled, "Gene stable ID" and column 2 should be titled "Gene name".
#' @return This function returns a count data set with gene IDs replaced by those provided in the meta_file parameter.
#' @keywords id transpose rename
#' @seealso \url{http://www.ensembl.org/biomart/martview/}, \code{\link{prep_counts}}, \code{\link{prep_dds_from_data}}
#' @export
#' @examples
#' \dontrun{
#'
#' #Transpose IDs from current count data to their desired values based on the meta_file.
#' count_data  <- prep_counts(my_count_file)
#' new_count_data <- transpose_gene_ids(count_data, meta_file)
#'
#' #Updated gene ID count information can be used in DE analysis.
#' dds <- prep_dds_from_data(new_count_data, target_data, ~ Design_Type, TRUE,
#'                            "Replicate", "vst")
#'
#' }
transpose_gene_ids <- function(count_data, meta_file)
{
  #Adjust the gene names. Mappings can be obtained from http://www.ensembl.org/biomart/, duplicates removed.
  meta_data = read.csv(meta_file, header=TRUE, sep="\t")
  count_data$id <- rownames(count_data)
  count_data$new_id <- as.character(meta_data$Gene.name[match(count_data$id, meta_data$Gene.stable.ID)])
  count_data$new_id <- ifelse(is.na(count_data$new_id), count_data$id, count_data$new_id)
  rownames(count_data) <- count_data$new_id
  count_data$new_id <- NULL
  count_data$id <- NULL

  return(count_data)
}




#' Create a composite metadata field by merging existing data.
#'
#' This function, given target metadata, will create a new composite column based on two or more existing columns
#' in the target data.  The new field will be named based on the merged fields and will be delimited using the "_" character.
#' Requires that target data has been prepared with prep_targets().
#' @param merge_fields Column names in targets data to be included in the new composite column. Should be provided as a set
#' of strings.  I.E. c("field1","field2").
#' @return This function returns target data containing the newly created composite field.
#' @keywords aggregate metadata
#' @seealso \code{\link{prep_targets}}
#' @export
#' @examples
#' \dontrun{
#'
#' myCounts <- prep_counts(count_input="master_count_data.txt", delim="t")
#' myTargets <- prep_targets(target_input="master_count_data.txt", delim="t")
#'
#' #Create a composite field based on "treatment" and "time" fields.
#' myTargets <- make_composite_field(targets=myTargets,
#'                                   merge_fields=c("treatment","time"))
#'
#'
#' #Create a composite field based on "treatment", "time", and "patientID" fields.
#' myTargets <- make_composite_field(targets=myTargets,
#'                                   merge_fields=c("treatment","time","patientID"))
#'
#'
#' }
make_composite_field <- function(merge_fields)
{
  if(typeof(merge_fields) != "character")
  {
    stop("Type Mismatch.  merge_fields must be a character type.")
    return(-1)
  }
  if(length(merge_fields) <= 1)
  {
    stop("Composite fields can only be created from 2 or more fields.  Please provide at least two fields.")
    return(-1)
  }
  if(!exists("tgt_dat", envir=.DEVis_env))
  {
    stop("make_composite_field(): Target data not present.  Please run prep_targets() first.")
    return(-1)
  }
  if(!(all(merge_fields %in% colnames(.DEVis_env$tgt_dat))))
  {
    print("make_composite_field requires that annotation correspond to available metadata.")
    print("Possible values are: ")
    print(colnames(.DEVis_env$tgt_dat))
    stop("make_composite_field(): Invalid target column provided.")
    return(-1)
  }

  #Merge the target data columns.
  tgt_dat <- .DEVis_env$tgt_dat
  if(length(merge_fields) > 2)
  {
    sub1 <- paste("paste(tgt_dat$", merge_fields[1], ",", sep="")
    for (i in 2:(length(merge_fields)-1))
    {
      sub1 <- paste(sub1, "tgt_dat$", merge_fields[i], ",", sep="")
    }
    sub1 <- paste(sub1, "tgt_dat$", merge_fields[length(merge_fields)], ", sep='_')", sep="")
  }
  else
  {
    sub1 <- paste("paste(tgt_dat$", merge_fields[1], ",", "tgt_dat$", merge_fields[2], ",sep='_')", sep="")
  }

  #Get the name of the new field.
  composite_field_name <- paste(merge_fields, collapse="_")

  #Create the composite column
  composite_column <- as.factor(as.character(eval(parse(text = sub1))))
  sub1 <- paste("tgt_dat$", composite_field_name, " <- composite_column", sep="")
  as.character(eval(parse(text = sub1)))

  #Update the stored target data.
  assign("tgt_dat", tgt_dat, envir=.DEVis_env)

  return(tgt_dat)
}




#' Create a new metadata field by renaming existing levels of an existing field.
#'
#' This function, given a target metadata field, will create a new column based on that field that has different labels.
#' This is particularly useful when cleaning figures for publication quality, as often labels will contain abbreviations or
#' delimiting characters such as "_".  This function allows for a new column to be generated with more human-friendly labels
#' that can be fed into visualizations instead of the defaults.
#' Requires that target data has been prepared with prep_targets().
#' @param target_column Column name in targets data to be included in the new column. String.
#' @param new_column_name Name for the new column that will be created. String.
#' @return This function returns target data containing the newly created label field.
#' @keywords aggregate metadata
#' @seealso \code{\link{prep_targets}}
#' @export
#' @examples
#' \dontrun{
#'
#' myCounts <- prep_counts(count_input="master_count_data.txt", delim="t")
#' myTargets <- prep_targets(target_input="master_count_data.txt", delim="t")
#'
#' #Create a  field based on the "treatment_time" fields with new labels.
#' myTargets <- create_relabel_field(target_column="treatment_time",
#'                                  new_column_name="treatment_time_relabel")
#'
#' }
create_relabel_field <- function(target_column, new_column_name)
{
  if(typeof(target_column) != "character")
  {
    stop("Type Mismatch. target_column must be a character type.")
    return(-1)
  }
  if(typeof(new_column_name) != "character")
  {
    stop("Type Mismatch. new_column_name must be a character type.")
    return(-1)
  }
  if(length(target_column) < 1)
  {
    stop("Target update label fields can only be created from 1 field at a time.")
    return(-1)
  }
  if(length(new_column_name) != 1)
  {
    stop("new_column_name fields can only be created from 1 string value.  Provide a single string for this parameter.")
    return(-1)
  }
  if(!exists("tgt_dat", envir=.DEVis_env))
  {
    stop("create_relabel_field(): Target data not present.  Please run prep_targets() first.")
    return(-1)
  }
  if(!(target_column %in% colnames(.DEVis_env$tgt_dat)))
  {
    print("create_relabel_field requires that annotation correspond to available metadata.")
    print("Possible values are: ")
    print(colnames(.DEVis_env$tgt_dat))
    stop("create_relabel_field(): Invalid target column provided.")
    return(-1)
  }
  if(new_column_name %in% colnames(.DEVis_env$tgt_dat))
  {
    stop("create_relabel_field(): Column provided for new_column_name already exists.  Provide a different column name.")
    return(-1)
  }

  #Create a copy of the initial column.
  tgt_dat <- .DEVis_env$tgt_dat
  sub1 <- paste("tgt_dat$", new_column_name, " <- as.character(tgt_dat$", target_column, ")", sep="")
  as.character(eval(parse(text = sub1)))

  #Get updated labels for each level in the original target column.
  sub1 <- paste("levels(tgt_dat$", target_column, ")", sep="")
  target_col_levels <- as.character(eval(parse(text = sub1)))
  newLabs <- c()
  for(i in 1:length(target_col_levels))
  {
    curLab <- readline(prompt=paste("Enter a value corresponding to ", target_col_levels[i], ": ", sep=""))
    curLab <- as.character(curLab)
    newLabs[i] <- curLab
  }

  #Replace the new column values with the user provided verisons.
  for(i in 1:length(target_col_levels))
  {
    sub1 <- paste("tgt_dat$", new_column_name, "[grepl('", target_col_levels[i], "',tgt_dat$", new_column_name, ")] <- '", newLabs[i], "'", sep="")
    print(sub1)
    eval(parse(text = sub1))
  }

  #Update the stored target data.
  assign("tgt_dat", tgt_dat, envir=.DEVis_env)

  return(tgt_dat)

}
