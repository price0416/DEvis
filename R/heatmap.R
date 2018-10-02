#' Create heat maps of differentially expressed genes.
#'
#' This function generates heat maps for differentially expressed gene result sets based on the experimental design.
#' Sorting and filtration can be applied to visualize genes based on specific criteria, and additional annotation
#' can be included in the heat map based on target metadata columns.
#' For example, aggregated results from several contrasts can have their differentially expressed genes sorted based on
#' a specific metric, and the top N genes can be plotted based on that sorting method.
#' @param res_list A list of DESeq result sets. Results can be calculated individually using DESeq's results() function.
#' Lists of results can be created by creating a list(result1, result2 ... result_N).
#' @param filename Filename for output plot.  Valid extensions are ".pdf" and ".png".
#' File generation can be turned off using set_output_mode("screen"). Output will be written to the /DE/heatmaps/ directory.
#' @param anno_columns Annotation column names as a single string or a concatenated list of strings. Must correspond
#' to columns in the targets metadata file. Multiple columns are allowed.  I.E. c("Condition", "Timepoint")
#' @param sort_choice Sorting method for DE genes. Possible options are: "mean", "max", "min", "variance", "max_mean","min_mean", "sd".
#' "max" sorts genes based on the highest expression level of any single gene in a result set.
#' In contrast, "max_mean" first calculates mean expression across all result sets and subsequently sorts by maximum values.
#' Min and min_mean function similarly.  Variance and sd sort genes by highest gene-wise variance or standard deviation.
#' @param specific_genes A character vector of gene names can be passed to this parameter to plot the genes specified. This overrides sorting
#' and numGene parameters.
#' @param numGenes Number of genes to include in the plot. Default: 30
#' @param theme Theme for the layout and color scheme for the plot.  Valid selections are integers between 1-6.
#' @param cluster_genes Boolean.  Cluster genes (rows) of the heat map.  Default: TRUE
#' @param cluster_contrasts Boolean.  Cluster contrasts (columns) of the heat map.  Default: TRUE
#' @param customLabels If customLabels is set to TRUE, the user will be prompted to provide a custom label for each label on the x-axis.
#' @param returnData Boolean.  Determines if this visualization should return data used to generate the visualization. Default=FALSE.
#' @return If returnData is true, this function will return a data frame containing expression data for the genes visualized in
#' the heat map based on the sorting method and numGenes parameters.
#' @keywords heatmap sort filter visualization
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
#' /*
#'  * Create a heat map of the top 25 most upregulated genes based on time and condition.
#'  * Gene-wise max value based calculation.
#'  */
#' de_heat(res_list=myResList, anno_columns=c("Time", "Condition"), sort_choice="max",
#'          numGenes=25, theme=2, returnData=FALSE)
#'
#'
#' /*
#'  * Create a heat map of the top 50 most downregulated genes based on time and condition.
#'  * Mean based value calculation.
#'  */
#' de_heat(res_list=myResList, anno_columns=c("Time", "Condition"), sort_choice="min_mean",
#'          numGenes=50, theme=2, returnData=FALSE)
#'
#'
#' /*
#'  * Create a heat map of the top 100 most highly varying genes based on time and response.
#'  * Variance based value calculation.
#'  */
#' de_heat(res_list=myResList, anno_columns=c("Time", "Response"), sort_choice="variance",
#'          numGenes=100, theme=2, returnData=FALSE)
#'
#' /*
#'  * Plot 3 specific genes, dont cluster by contrast.
#'  */
#' de_heat(res_list=myResList, anno_columns=c("Time", "Response"), sort_choice="variance",
#'         specific_genes=c("GEN1", "ABC2", "FuSG2"), theme=2, cluster_contrasts=FALSE)
#'
#' }
de_heat <- function(res_list, filename="heatmap.pdf", anno_columns, sort_choice="none", specific_genes="", numGenes=30, theme=2, cluster_genes=TRUE, cluster_contrasts=TRUE, customLabels=FALSE, returnData=FALSE)
{
  #Parameter validation.
  sort_choice <- tolower(sort_choice)
  if((typeof(res_list) != "list") || (typeof(res_list[[1]]) != "S4") || (length(res_list) <= 0))
  {
    #If they arent passing a list of results, maybe its a set of genes.
    if(typeof(res_list) != "character")
    {
      stop("de_heat() requires list type object containing DESeq result sets or a character vector of gene names..")
      return(-1)
    }
  }
  if(!exists("output_mode", envir=.DEVis_env))
  {
    stop("de_heat(): Output mode not initialized.  Please run set_output_mode() first.")
    return(-1)
  }
  if(typeof(anno_columns) != "character")
  {
    stop("de_heat(): Type mismatch. Requires anno_columns parameter be a character vector.")
    return(-1)
  }
  if(typeof(filename) != "character")
  {
    stop("filename must be a string.")
    return(-1)
  }
  if(typeof(specific_genes) != "character")
  {
    stop("specific_genes must be a character vector.")
    return(-1)
  }
  specify = FALSE
  if(specific_genes != "" && length(specific_genes >= 1))
  {
    specify = TRUE
  }
  if(!is.numeric(numGenes))
  {
    stop("numGenes must be numeric.")
    return(-1)
  }
  if(typeof(returnData) != "logical")
  {
    stop("returnData must be a logical value. Please enter TRUE or FALSE.")
    return(-1)
  }
  if(typeof(cluster_genes) != "logical")
  {
    stop("cluster_genes must be a logical value. Please enter TRUE or FALSE.")
    return(-1)
  }
  if(typeof(cluster_contrasts) != "logical")
  {
    stop("cluster_contrasts must be a logical value. Please enter TRUE or FALSE.")
    return(-1)
  }
  if(!exists("DE_hm_base_dir", envir=.DEVis_env))
  {
    stop("Directories not initialized.  Try running create_dir_struct().")
    return(-1)
  }
  type_options = c(".pdf",".png")
  filetype = substr(filename, nchar(filename)-4+1, nchar(filename))
  if(!(filetype %in% type_options))
  {
    stop('filename must have extension ".pdf" or ".png".')
    return(-1)
  }
  if(!exists("tgt_dat", envir=.DEVis_env))
  {
    stop("de_heat(): Target data not present.  Please run prep_targets() first.")
    return(-1)
  }
  if(!(all(anno_columns %in% colnames(.DEVis_env$tgt_dat))))
  {
    print("de_heat() requires that annotation correspond to available metadata.")
    print("Possible values are: ")
    print(colnames(.DEVis_env$tgt_dat))
    stop("de_heat(): Invalid target column provided.")
    return(-1)
  }
  sort_options = c("mean","max","min","sd","variance","max_mean","min_mean","none")
  if(!(sort_choice %in% sort_options))
  {
    stop('sort_choice must be one of the following: "mean","max","min","sd","variance","max_mean","min_mean","none"')
    return(-1)
  }
  if(typeof(customLabels) != "logical")
  {
    stop("customLabels must be a logical value. Please enter TRUE or FALSE.")
    return(-1)
  }
  if(theme < 1 || theme > 6)
  {
    stop("theme must be a value between 1 and 6.")
    return(-1)
  }

  #Enrich the data.
  enriched <- unlist(lapply(res_list, enrich_res))
  enrich <- unlist(enriched)

  #Pick the appropriate genes.
  if(specify)
  {
    aggregate.names <- specific_genes
    numGenes = length(specific_genes)
  }
  else
  {
    aggregate.names <- .DEVis_env$aggregate_names
    if(length(aggregate.names) == 0)
    {
      stop("Error: No genes identified in aggregated data.  Consider union-based aggregation in create_master_res. ")
      return(-1)
    }
  }

  #Get the aggregated genes data.
  master_aggregate_names <- aggregate.names
  master_aggregate <- lapply(enrich, function(x) x[aggregate.names,])
  master_lfc   <- lapply(master_aggregate, function(x) x$log2FoldChange)
  master_colnames_lfc <- lapply(master_aggregate, function(x) x@contrast)
  master_colnames <- master_colnames_lfc

  #Finalize the master data frame.
  master_df <- data.frame(master_lfc)
  colnames(master_df) <- master_colnames
  rownames(master_df) <- rownames(master_aggregate[[1]])

  #Create annotation data frame.
  contrast_field <- master_aggregate[[1]]@design_field
  anno_df = data.frame(placeholder=integer(length(master_aggregate)))
  for (item in anno_columns)
  {
    cur_item_list <- list()
    for (i in 1:length(master_aggregate))
    {
      tgt_dat <- .DEVis_env$tgt_dat
      evalStr = paste("tgt_dat[which(tgt_dat$", contrast_field, '=="', master_aggregate[[i]]@case, '"),][1,]$',item, sep="")
      cur_val <- as.character(eval(parse(text = evalStr)))
      cur_item_list <- c(cur_item_list, cur_val)
    }
    anno_df[item] <- unlist(cur_item_list)
  }
  anno_df$placeholder <- NULL
  rownames(anno_df) <- colnames(master_df)

  #Handle any sorting options.
  if(sort_choice != "none")
  {
    if(sort_choice == "max")
    {
      maxCol <- apply(master_df, 1, max)
      master_df[,"max"] <- maxCol
      master_df <- master_df[with(master_df, order(-maxCol)), ]
      master_df$max <- NULL
    }
    if(sort_choice == "min")
    {
      minCol <- apply(master_df, 1, min)
      master_df[,"min"] <- minCol
      master_df <- master_df[with(master_df, order(minCol)), ]
      master_df$min <- NULL
    }
    if(sort_choice == "mean")
    {
      meanCol <- apply(master_df, 1, mean)
      master_df[,"mean"] <- minCol
      master_df <- master_df[order(master_df$mean),]
      master_df$mean <- NULL
    }
    if(sort_choice == "variance")
    {
      varCol <- apply(master_df, 1, var)
      master_df[,"variance"] <- varCol
      master_df <- master_df[with(master_df, order(-varCol)), ]
      master_df$variance <- NULL
    }
    if(sort_choice == "sd")
    {
      sdCol <- apply(master_df, 1, sd)
      master_df[,"sd"] <- sdCol
      master_df <- master_df[with(master_df, order(-sdCol)), ]
      master_df$sd <- NULL
    }
    if(sort_choice == "max_mean")
    {
      meanCol <- apply(master_df, 1, mean)
      master_df[,"mean"] <- meanCol
      master_df <- master_df[with(master_df, order(-meanCol)), ]
      master_df$mean <- NULL
    }
    if(sort_choice == "min_mean")
    {
      meanCol <- apply(master_df, 1, mean)
      master_df[,"mean"] <- meanCol
      master_df <- master_df[with(master_df, order(meanCol)), ]
      master_df$mean <- NULL
    }
  }

  #Narrow down to the number of genes desired.
  master_df <- head(master_df,numGenes)

  #Automatically scale the image so all the genes fit and labels are readable.
  imgHeight = 10
  imgWidth = 10
  if(numGenes > 50)
  {
    imgHeight = (numGenes / 10) + 8
  }
  if(dim(master_df)[2] > 15)
  {
    imgWidth = (dim(master_df)[2] / 10) + 10
  }

  paletteLength = 100

  #Select the theme.
  if(theme == 1)
  {
    colSet = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(paletteLength)
  }
  if(theme == 2)
  {
    colSet = colorRampPalette(c("navy", "white", "firebrick3"))(paletteLength)
  }
  if(theme == 3)
  {
    colSet = colorRampPalette(c("cornflowerblue", "white", "darkorange1"))(paletteLength)
  }
  if(theme == 4)
  {
    colSet = colorRampPalette(c("royalblue4", "white", "orangered4"))(paletteLength)
  }
  if(theme == 5)
  {
    colSet = colorRampPalette(c("springgreen3", "white", "orangered2"))(paletteLength)
  }
  if(theme == 6)
  {
    colSet <- colorRampPalette(c("gray99", "gray10"))(paletteLength)
  }

  #Adjust the color scale so that 0 is always the white part of the color palette.
  break_points <- c(seq(min(as.matrix(master_df)), 0, length.out=ceiling(paletteLength/2) + 1),
                seq(max(as.matrix(master_df))/paletteLength, max(as.matrix(master_df)), length.out=floor(paletteLength/2)))


  #If custom labels are requested retrieve and update them here.
  if(customLabels == TRUE)
  {
    curLabs = as.character(colnames(master_df))
    relabel <- as.character(colnames(master_df))

    #Get updated labels for each level in the current label set.
    for(i in 1:length(curLabs))
    {
      curLab <- as.character(readline(prompt=paste("Enter a custom label corresponding to ", curLabs[i], ": ", sep="")))
      sub1 <- paste("relabel[grepl('", curLabs[i], "',relabel)] <- '", curLab, "'", sep="")
      eval(parse(text = sub1))
    }
    rownames(anno_df) <- relabel
    colnames(master_df) <- relabel
  }

  #Draw the heat map.
  setwd(.DEVis_env$DE_hm_base_dir)
  if(.DEVis_env$output_mode == "both" | .DEVis_env$output_mode == "file")
  {
    if(filetype == ".pdf")
    {
      pdf(filename, width=imgWidth, height=imgHeight, onefile=FALSE)
    }
    else
    {
      png(filename, width=imgWidth, height=imgHeight, units='in', res=600)
    }
    pheatmap(master_df, annotation_col=anno_df, color=colSet, breaks=break_points, cluster_cols=cluster_contrasts, cluster_rows=cluster_genes, fontsize_row=13, fontsize_col=10)
    dev.off()
  }
  if(.DEVis_env$output_mode == "screen" | .DEVis_env$output_mode == "both")
  {
    pheatmap(master_df, annotation_col=anno_df, color = colSet, breaks=break_points, cluster_cols=cluster_contrasts, cluster_rows=cluster_genes, fontsize_row=13, fontsize_col=10)
  }
  setwd(.DEVis_env$working_dir)

  if(returnData)
  {
    return(master_df)
  }
}





