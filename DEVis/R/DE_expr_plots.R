#' Visualize fold-change divergence for differentially expressed genes.
#'
#' This function plots the log2 fold-change values for all differentially expressed genes for each contrast in
#' a result set.  This plot visualizes the distribution and strength of expression changes for all differentially expressed genes.
#' @param res_list A list of DESeq result sets. Results can be calculated individually using DESeq's results() function.
#' Lists of results can be created by creating a list(result1, result2 ... result_N).
#' @param filename Filename for output plot.  Valid extensions are ".pdf" and ".png".
#' File generation can be turned off using set_output_mode("screen"). Output will be written to the /DE/divergence/ directory.
#' @param customLabels If customLabels is set to TRUE, the user will be prompted to provide a custom label for each label.
#' @param theme Theme for the layout and color scheme for the plot.  Valid selections are integers between 1-6.
#' @param returnData Boolean.  Determines if this visualization should return data used to generate the visualization. Default=FALSE.
#' @return If returnData is true, this function will return the long-form table for differentially expressed genes
#' containing gene names, categorical variable, and expression values.
#' @keywords DE fold-change expression visualization
#' @export
#' @examples
#' \dontrun{
#'
#' #Prepare a result list.
#' res.day1 <- results(dds, contrast=c("Condition_Time", "day1_disease", "day1_control"))
#' res.day2 <- results(dds, contrast=c("Condition_Time", "day2_disease", "day2_control"))
#' res.day3 <- results(dds, contrast=c("Condition_Time", "day3_disease", "day3_control"))
#' myResList <- list(res.day1, res.day2, res.day3)
#'
#' #Create the plot.
#' de_diverge_plot(res_list=myResList, filename="DE_divergence_plot.pdf",
#'                  theme=1, returnData=FALSE)
#'
#' }
de_diverge_plot <- function(res_list, filename="de_divergence_plot.pdf", theme=1, customLabels=FALSE, returnData=FALSE)
{
  imgWidth  = 10
  imgHeight = 10

  if((typeof(res_list) != "list") || (typeof(res_list[[1]]) != "S4") || (length(res_list) <= 0))
  {
    stop("de_diverge_plot() requires list type object containing DESeq result sets.")
    return(-1)
  }
  if(!exists("output_mode", envir=.DEVis_env))
  {
    stop("de_diverge_plot(): Output mode not initialized.  Please run set_output_mode() first.")
    return(-1)
  }
  if(typeof(filename) != "character")
  {
    stop("filename must be a string.")
    return(-1)
  }
  if(!exists("DE_diverge_dir", envir=.DEVis_env))
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
  if(theme < 1 || theme > 6)
  {
    stop("theme must be a value between 1 and 6.")
    return(-1)
  }
  if(typeof(customLabels) != "logical")
  {
    stop("customLabels must be a logical value. Please enter TRUE or FALSE.")
    return(-1)
  }
  if(typeof(returnData) != "logical")
  {
    stop("returnData must be a logical value. Please enter TRUE or FALSE.")
    return(-1)
  }

  #Get DE gene set for each contrast.
  enriched <- unlist(lapply(res_list, enrich_res))
  enrich <- unlist(enriched)
  res.names <- lapply(enrich, function(x) x@allNames)
  filteredList <- list()
  for (i in 1:length(enrich))
  {
    enrich[[i]]$gene <- rownames(enrich[[i]])
    t1 <- data.frame(unlist(enrich[[i]]$log2FoldChange[enrich[[i]]$gene %in% res.names[[i]]]))
    t1$gene <- unlist(res.names[i])
    colnames(t1) <- c(enrich[[i]]@case, "gene")
    rownames(t1) <- t1$gene

    melt_t1 <- melt(t1,id.var="gene")
    filteredList <- rbind(filteredList, melt_t1)
  }

  #Satisfy check() that our ggplot aes variables are indeed not globals.
  value    = NULL
  variable = NULL

  #If custom labels are requested retrieve and update them here.
  if(customLabels == TRUE)
  {
    curLabs = as.character(levels(filteredList$variable))
    relabel <- as.character(filteredList$variable)

    #Get updated labels for each level in the current label set.
    for(i in 1:length(curLabs))
    {
      curLab <- as.character(readline(prompt=paste("Enter a custom label corresponding to ", curLabs[i], ": ", sep="")))
      sub1 <- paste("relabel[grepl('", curLabs[i], "',relabel)] <- '", curLab, "'", sep="")
      eval(parse(text = sub1))
    }
    filteredList$variable <- relabel
  }

  #Prepare plot.
  plot_out <- ggplot(filteredList, aes(x=value, y=factor(variable, levels=as.character(unique(variable)), ordered=TRUE), label=variable)) + geom_point(stat='identity', aes(colour=variable), size=3)

  #Apply the selected theme.
   if(theme == 1)
   {
     stata_long_pal = c(stata_pal("s2color")(15), stata_pal("s1rcolor")(15))
     plot_out <- plot_out + theme_stata() + scale_fill_manual(values=stata_long_pal) + theme(axis.text.y = element_text(angle = 45, hjust = 1), plot.title = element_text(size=14, face="bold", hjust=0)) + guides(col=guide_legend(ncol=length(enrich)%/%3))

   }
  if(theme == 2)
  {
    nature_pal = c(pal_npg("nrc")(10),pal_npg("nrc", alpha = 0.6)(10),pal_npg("nrc", alpha = 0.3)(10) )
    plot_out <- plot_out + theme_igray() + scale_fill_manual(values=nature_pal) + guides(col=guide_legend(ncol=length(enrich)%/%3))
  }
  if(theme == 3)
  {
    tron_pal = c(pal_tron("legacy")(7),pal_tron("legacy", alpha = 0.7)(7),pal_tron("legacy", alpha = 0.5)(7),pal_tron("legacy", alpha = 0.3)(7),pal_tron("legacy", alpha = 0.2)(2))
    plot_out <- plot_out + theme_hc(bgcolor = "darkunica") + theme(axis.text.x = element_text(colour = "white"), axis.text.y = element_text(colour = "white")) + scale_fill_manual(values=tron_pal) + guides(col=guide_legend(ncol=length(enrich)%/%3))
  }
  if(theme == 4)
  {
    gdoc_pal <- c(pal_ucscgb("default")(26), pal_ucscgb("default",alpha=.5)(4))
    plot_out <- plot_out + theme_gdocs() + scale_fill_manual(values=gdoc_pal) + guides(col=guide_legend(ncol=length(enrich)%/%3))
  }
  if(theme == 5)
  {
    d3_pal <- c(pal_d3("category20")(20), pal_d3("category10",alpha=.5)(10))
    plot_out <- plot_out + theme_solarized() + scale_fill_manual(values=d3_pal) + guides(col=guide_legend(ncol=length(enrich)%/%3))
  }
  if(theme == 6)
  {
    plot_out <- plot_out + theme_bw() + scale_colour_grey(start = 0, end = .9) + scale_fill_grey(start = 0, end = .9) + guides(col=guide_legend(ncol=length(enrich)%/%3))
  }

  plot_out <- plot_out +
              theme(legend.text = element_text(size=12)) +
              theme(legend.title=element_blank()) +
              theme(axis.text=element_text(size=12, face="bold")) +
              theme(text = element_text(size=22,margin = margin(t = 0, r = 10, b = 0, l = 0))) +
              labs(x="", y="", title="Differentially Expressed Genes", subtitle="Log2 Fold-Change")

  #Write the plot to file.
  setwd(.DEVis_env$DE_diverge_dir)
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

    print(plot_out)
    dev.off()
  }
  if(.DEVis_env$output_mode == "screen" | .DEVis_env$output_mode == "both")
  {
    print(plot_out)
  }
  setwd(.DEVis_env$working_dir)

  #Return the data if the user wants it.
  if(returnData)
  {
    return(filteredList)
  }
}


#' Identify and visualize patterns of expression between differentially expressed genes across a series of result sets.
#'
#' This function identifies patterns of expression across differentially expressed genes and clusters them into
#' groups based on similar patterns of expression across multiple series, such as time series data (I.E. day1, day2, day3).  The resulting
#' group data are plotted to show how expression of groups change over the course of the series based on a generalized
#' linear model with a displayed 95 percent confidence interval, or based on the mean expression for each point in the series.
#' Similarity is calculated by clustering and merging genes until all differentially expressed genes have been clustered into
#' the desired number of groups.
#' @param res_list A list of DESeq result sets. Results can be calculated individually using DESeq's results() function.
#' Lists of results can be created by creating a list(result1, result2 ... result_N).
#' @param filename Filename for output plot.  Valid extensions are ".pdf" and ".png".
#' File generation for plot file can be turned off using set_output_mode("screen").
#' Output will be written to the /DE/series_plots/ directory.
#' @param designVar The design field used when performing results() extractions.
#' Must correspond to a metadata column. String.
#' @param groupBy The variable by which to group data in the plot.  I.E. "Condition_Time". String.
#' @param numGroups The number of groups to split data into.  Default = 5
#' @param theme Theme for the layout and color scheme for the plot.  Valid selections are integers between 1-6.
#' @param method Method for drawing.  "glm" uses generalized linear model to show overall trends over time.  "mean" uses the mean expression
#' value for each group as intersection points.
#' @param customLabels If customLabels is set to TRUE, the user will be prompted to provide a custom label for each label.
#' @param returnData If this value is true, this function will return a list of dataframes corresponding to groupwise splits. Default: FALSE
#' @param writeData If this value is true, this function will write data to file in the DE series folder for each group of genes. Default: FALSE
#' @param fetchGroup Specifically fetch gene information for a specific group.  This will override standard returnData and instead return data for the specified
#' group number.
#' @return If returnData is true, this function will return a list of dataframes corresponding to groupwise splits.
#' @keywords DE similarity aggregate cluster visualization
#' @export
#' @examples
#' \dontrun{
#'
#'  #This example assumes an experimenal design of ~Condition_Time.
#'
#' #Prepare a result list.
#' res.day1 <- results(dds, contrast=c("Condition_Time", "day1_disease", "day1_control"))
#' res.day2 <- results(dds, contrast=c("Condition_Time", "day2_disease", "day2_control"))
#' res.day3 <- results(dds, contrast=c("Condition_Time", "day3_disease", "day3_control"))
#' myResList <- list(res.day1, res.day2, res.day3)
#'
#' /*
#'  * Cluster genes by similarity into 5 groups, then visualize their expression over the
#'  * course of the series using a generalized linear model.
#'  */
#' de_series(res_list=myResList, filename="DE_series_pattern.pdf",
#'              designVar="Condition_Time",
#'              groupBy="Time", numGroups=5, theme=1, method="glm",
#'              returnData=FALSE, writeData=FALSE)
#'
#'
#' /*
#'  * Cluster genes by similarity into 3 groups, then visualize their expression over the
#'  * course of the series using based on mean group expression values.
#'  */
#' de_series(res_list=myResList, filename="DE_series_pattern.pdf",
#'              designVar="Condition_Time",
#'              groupBy="Time", numGroups=3, theme=2, method="mean",
#'              returnData=FALSE, writeData=FALSE)
#'
#'  }
de_series <- function(res_list, filename="series_pattern.pdf", designVar, groupBy, numGroups=5, theme=1, method="mean", customLabels=FALSE, returnData=FALSE, writeData=FALSE, fetchGroup=0)
{
  imgHeight = 10
  imgWidth = 10

  #Parameter validation.
  if((typeof(res_list) != "list") || (typeof(res_list[[1]]) != "S4") || (length(res_list) <= 0))
  {
    stop("de_series() requires list type object containing DESeq result sets.")
    return(-1)
  }
  if(!exists("tgt_dat", envir=.DEVis_env))
  {
    stop("Target data not present.  Please run prep_targets() first.")
    return(-1)
  }
  if(!exists("stabilized_data", envir=.DEVis_env))
  {
    stop("de_series() requires that prep_dds_from_data() be run first.")
    return(-1)
  }
  if(!exists("output_mode", envir=.DEVis_env))
  {
    stop("de_series(): Output mode not initialized.  Please run set_output_mode() first.")
    return(-1)
  }
  if(typeof(filename) != "character")
  {
    stop("filename must be a string.")
    return(-1)
  }
  if(typeof(customLabels) != "logical")
  {
    stop("customLabels must be a logical value. Please enter TRUE or FALSE.")
    return(-1)
  }
  if(typeof(returnData) != "logical")
  {
    stop("returnData must be a logical value. Please enter TRUE or FALSE.")
    return(-1)
  }
  if(typeof(writeData) != "logical")
  {
    stop("writeData must be a logical value. Please enter TRUE or FALSE.")
    return(-1)
  }
  if(!exists("DE_series_dir", envir=.DEVis_env))
  {
    stop("Directories not initialized.  Try running create_dir_struct().")
    return(-1)
  }
  if(!is.numeric(fetchGroup))
  {
    stop("fetchGroup must be an integer.")
    return(-1)
  }
  if(fetchGroup > numGroups)
  {
    stop("fetchGroup cannot be higher than numGroups.")
    return(-1)
  }
  type_options = c(".pdf",".png")
  filetype = substr(filename, nchar(filename)-4+1, nchar(filename))
  if(!(filetype %in% type_options))
  {
    stop('filename must have extension ".pdf" or ".png".')
    return(-1)
  }
  type_options = c("glm","mean")
  if(!(method %in% type_options))
  {
    stop('method must be "glm" or "mean".')
    return(-1)
  }
  if(!(groupBy %in% colnames(colData(.DEVis_env$stabilized_data))))
  {
    print("de_series() requires grouping to correspond to colData annotation.")
    print("Possible values are: ")
    print(colnames(colData(.DEVis_env$stabilized_data)))
    stop("de_series(): invalid group variable selection.")
    return(-1)
  }
  if(!(designVar %in% colnames(colData(.DEVis_env$stabilized_data))))
  {
    print("de_series() requires grouping to correspond to colData annotation.")
    print("Possible values are: ")
    print(colnames(colData(.DEVis_env$stabilized_data)))
    stop("de_series(): invalid group variable selection.")
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

  #Get the master aggregate DE data.
  aggregate.names <- .DEVis_env$aggregate_names
  if(length(aggregate.names) == 0)
  {
    stop("Error: No genes identified in aggregated data.  Consider union-based aggregation in create_master_res. ")
    return(-1)
  }

  master_aggregate_names <- aggregate.names
  master_aggregate <- lapply(enrich, function(x) x[aggregate.names,])
  master_lfc   <- lapply(master_aggregate, function(x) x$log2FoldChange)
  master_padj  <- lapply(master_aggregate, function(x) x$padj)
  master_colnames_lfc <- lapply(master_aggregate, function(x) x@case)
  master_colnames_padj <- lapply(master_aggregate, function(x) paste(x@case, ".padj", sep=""))
  master_colnames <- c(master_colnames_lfc)

  #Create the master data frame.
  master_df <- data.frame(master_lfc)
  colnames(master_df) <- master_colnames
  rownames(master_df) <- rownames(master_aggregate[[1]])

  #Find similar groups of genes by clustering row-wise.
  lfc_data <- get_de_data(res_list,method="union")
  dist_data <- dist(lfc_data, method = "euclidean")
  fit <- hclust(dist_data, method = "ward.D")
  groups <- cutree(fit, k=numGroups)


  #For each group, extract the correct subset of data from the master_df.
  setwd(.DEVis_env$DE_series_dir)
  group_names   <- c()
  result_dfs <- list()
  final_melt <- list()
  mean_dat   <- list()
  file_prefix = substr(filename, 1, nchar(filename)-4)
  for (i in 1:numGroups)
  {
    group_names   <- names(groups[which(groups==i)])
    result_dfs[[i]] <- master_df[group_names,]
    sdCol <- apply(result_dfs[[i]], 1, sd)
    result_dfs[[i]]$sd <- sdCol

    #Trim the outlying genes or we'll have huge ranges of error.
    #Trimming is based on absolute value of the current groups standard deviation * 4.
    lfc_keep <- result_dfs[[i]][rowSums(abs(result_dfs[[i]]) > abs(mean(result_dfs[[i]]$sd)*4)) < 1, ]
    lfc_keep_genes <- rownames(lfc_keep)
    result_dfs[[i]]$gene <- rownames(result_dfs[[i]])
    df_trim <- result_dfs[[i]][result_dfs[[i]]$gene %in% lfc_keep_genes,]
    result_dfs[[i]] <- df_trim
    result_dfs[[i]]$gene <- NULL
    result_dfs[[i]]$sd <- NULL

    #Write the data for each group if the user requested it.
    if(writeData)
    {
      curGroup_filename = paste(file_prefix, "_group", i, ".txt", sep="")
      write.table(result_dfs[[i]], curGroup_filename, sep="\t", row.names=TRUE, col.names=NA)
    }

    #If they are requesting data for a specific group, store it here and return it at the end of the function.
    if(fetchGroup != 0 && i == fetchGroup)
    {
      fetchData = result_dfs[[fetchGroup]]
      print("fetch")
      print(i)

    }

    #Adjust column names to be whatever target field the user wants.  Needs to be the same for all group splits so doing it in the loop.
    tgt_dat <- .DEVis_env$tgt_dat
    sub1 <- paste("tgt_dat$", groupBy, "[match(colnames(result_dfs[[i]]),tgt_dat$", designVar, ")]", sep="")
    newColnames <- as.character(eval(parse(text = sub1)))
    colnames(result_dfs[[i]]) <- newColnames

    if(method == "mean")
    {
      nameSet  <- colnames(result_dfs[[i]])
      meanSet  <- as.numeric(lapply(result_dfs[[i]], mean))
      groupVar <- rep((paste("group", i, sep="")), length(nameSet))
      cBound   <- data.frame(nameSet,meanSet,groupVar)
      mean_dat <- rbind(mean_dat, cBound)
    }

    #Melt each data frame and assign the data to its appropriate group.
    suppressMessages({
      curMelt <- melt(result_dfs[[i]])
    })
    groupVar <- rep((paste("group", i, sep="")), nrow(curMelt))
    curMelt <- cbind(curMelt,groupVar)
    final_melt <- rbind(final_melt,curMelt)

  }

  #Satisfy check() that our ggplot aes variables are indeed not globals.
  value    = NULL
  variable = NULL

  if(method == "mean")
  {
    colnames(mean_dat) <- c("variable","value","groupVar")
    mean_dat <- data.frame(mean_dat)
  }

  #If custom labels are requested retrieve and update them here.
  if(customLabels == TRUE)
  {
    if(method == "mean")
    {
      curLabs = as.character(levels(mean_dat$variable))
      relabel <- as.character(mean_dat$variable)
    }
    if(method == "glm")
    {
      curLabs = as.character(levels(final_melt$variable))
      relabel <- as.character(final_melt$variable)
    }

    #Get updated labels for each level in the current label set.
    for(i in 1:length(curLabs))
    {
      curLab <- as.character(readline(prompt=paste("Enter a custom label corresponding to ", curLabs[i], ": ", sep="")))
      sub1 <- paste("relabel[grepl('", curLabs[i], "',relabel)] <- '", curLab, "'", sep="")
      eval(parse(text = sub1))
    }

    if(method == "mean")
    {
      mean_dat$variable <- relabel
    }
    if(method == "glm")
    {
      final_melt$labels <- relabel
    }
  }

  #If its the mean method set the plotting data to be the mean info instead of the full melt.
  if(method == "mean")
  {
    plot_out <- ggplot(mean_dat, aes(x=factor(variable, levels=as.character(unique(variable)), ordered=TRUE), y=value, colour = groupVar, group=groupVar)) + geom_line() + geom_point()
  }
  if(method == "glm")
  {
    plot_out <- ggplot(final_melt, aes(x=factor(as.numeric(variable), levels=as.character(unique(variable)), ordered=TRUE), y = value, colour = groupVar)) + geom_smooth(aes(x = as.numeric(variable), y = value), method = 'glm')
  }

  #Apply the selected theme.
  if(theme == 1)
  {
    stata_long_pal = c(stata_pal("s2color")(15), stata_pal("s1rcolor")(15))
    plot_out <- plot_out + theme_stata() + scale_fill_manual(values=stata_long_pal)
  }
  if(theme == 2)
  {
    nature_pal = c(pal_npg("nrc")(10),pal_npg("nrc", alpha = 0.6)(10),pal_npg("nrc", alpha = 0.3)(10) )
    plot_out <- plot_out + theme_igray() + scale_fill_manual(values=nature_pal)
  }
  if(theme == 3)
  {
    tron_pal = c(pal_tron("legacy")(7),pal_tron("legacy", alpha = 0.7)(7),pal_tron("legacy", alpha = 0.5)(7),pal_tron("legacy", alpha = 0.3)(7),pal_tron("legacy", alpha = 0.2)(2))
    plot_out <- plot_out + theme_hc(bgcolor = "darkunica") + theme(axis.text.x = element_text(colour = "white"), axis.text.y = element_text(colour = "white")) + scale_fill_manual(values=tron_pal)
  }
  if(theme == 4)
  {
    gdoc_pal <- c(pal_ucscgb("default")(26), pal_ucscgb("default",alpha=.5)(4))
    plot_out <- plot_out + theme_gdocs() + scale_fill_manual(values=gdoc_pal)
  }
  if(theme == 5)
  {
    d3_pal <- c(pal_d3("category20")(20), pal_d3("category10",alpha=.5)(10))
    plot_out <- plot_out + theme_solarized() + scale_fill_manual(values=d3_pal)
  }
  if(theme == 6)
  {
    plot_out <- plot_out + theme_bw() + scale_colour_grey(start = 0, end = .9) + scale_fill_grey(start = 0, end = .9) + guides(col=guide_legend(ncol=length(enrich)%/%3))
  }

  #Finish off the themes.
  if(method == "mean")
  {
    plot_out <- plot_out +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      theme(plot.title = element_text(size=14, face="bold", hjust=0)) +
      theme(legend.text = element_text(size=12)) +
      theme(legend.title=element_blank()) +
      theme(axis.text=element_text(size=12, face="bold")) +
      theme(text = element_text(size=16,margin = margin(t = 0, r = 10, b = 0, l = 0))) +
      labs(x="", y="Log Fold-Change", title="Mean Progression of Clustered Genes")
  }
  if(method == "glm")
  {
    plot_out <- plot_out +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      theme(plot.title = element_text(size=14, face="bold", hjust=0))

    if(customLabels == TRUE)
    {
      plot_out <- plot_out +
        scale_x_continuous(breaks=as.numeric(final_melt$variable),
                           labels = final_melt$labels,
                           minor_breaks = NULL)
    }
    else
    {
      plot_out <- plot_out +
        scale_x_continuous(breaks=as.numeric(final_melt$variable),
                           labels = final_melt$variable,
                           minor_breaks = NULL)
    }

    plot_out <- plot_out +
      theme(legend.text = element_text(size=12)) +
      theme(legend.title=element_blank()) +
      theme(axis.text=element_text(size=12, face="bold")) +
      theme(text = element_text(size=16,margin = margin(t = 0, r = 10, b = 0, l = 0))) +
      labs(x="", y="Log Fold-Change", title="Linear Progression of Clustered Genes")
  }

  #Write the plot.
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
    print(plot_out)
    dev.off()
  }
  if(.DEVis_env$output_mode == "screen" | .DEVis_env$output_mode == "both")
  {
    print(plot_out)
  }
  setwd(.DEVis_env$working_dir)

  #Return the data to the user if requested.
  if(fetchGroup != 0)
  {
    return(fetchData)
  }
  if(returnData)
  {
    if(method == "mean")
    {
      return(mean_dat)
    }
    else
    {
      return(final_melt)
    }
  }

}


#' Visualize gene-wise expression of differentially expressed genes.
#'
#' This function plots log2 fold-change values for differentially expressed genes for each contrast in
#' a result set.  The set of genes displayed can be selected by means of several sorting methods.  This makes
#' it possible to view expression differences from a variety of perspectives.  This function can be applied to
#' a single or multiple result sets, making it possible to compare expression changes in specific genes across
#' different experimental conditions.
#' @param res_list A list of DESeq result sets. Results can be calculated individually using DESeq's results() function.
#' Lists of results can be created by creating a list(result1, result2 ... result_N).
#' @param filename Filename for output plot.  Valid extensions are ".pdf" and ".png".
#' File generation can be turned off using set_output_mode("screen").
#' Output will be written to the /DE/profile_plots/ directory.
#' @param sort_choice Gene selection is based on sorting method in cases where not all genes are displayed.
#' sort_choice options are: "max", "min", "variance", "max_mean", "min_mean".
#' "max" sorts genes based on the highest expression level of any single gene in a result set.
#' In contrast, "max_mean" first calculates mean expression across all result sets and subsequently sorts by maximum mean values.
#' Min and min_mean function similarly.  Variance sorts genes by highest gene-wise variance in expression,
#' displaying genes that showed the highest variability across all samples.
#' @param specific_genes A character vector of gene names can be passed to this parameter to plot the genes specified. This overrides sorting
#' and numGene parameters.
#' @param numGenes The number of genes to include in this plot.
#' @param theme Theme for the layout and color scheme for the plot.  Valid selections are integers between 1-6.
#' @param customLabels If customLabels is set to TRUE, the user will be prompted to provide a custom label for each label.
#' @param returnData Boolean.  Determines if this visualization should return data used to generate the visualization. Default=FALSE.
#' @return If returnData is true, this function will return a data frame for sort-selected differentially expressed genes
#' containing gene names and log2 fold-change values relative to the experimental control condition.
#' @keywords DE expression fold-change aggregate sort visualization
#' @export
#' @examples
#' \dontrun{
#'
#'  #This example assumes an experimenal design of ~Condition_Time.
#'
#' #Prepare a result list.
#' res.day1 <- results(dds, contrast=c("Condition_Time", "day1_disease", "day1_control"))
#' res.day2 <- results(dds, contrast=c("Condition_Time", "day2_disease", "day2_control"))
#' res.day3 <- results(dds, contrast=c("Condition_Time", "day3_disease", "day3_control"))
#' myResList <- list(res.day1, res.day2, res.day3)
#'
#' /*
#'  * Sort data by the highest expression level for any individual gene in any sample.
#'  * Select the top 50 genes from this sort and visualize them in the plot.
#'  */
#' de_profile_plot(res_list=myResList, filename="DE_profile_upReg50.pdf",
#'                   sort_choice="max",
#'                   numGenes=50, theme=1, returnData=FALSE)
#'
#' /*
#'  * Calculate the mean expression for each gene across all three time points.
#'  * Sort the data by minimum mean expression, select the top 25 genes,
#'  * and visualize them in the plot.
#'  */
#' de_profile_plot(res_list=myResList, filename="DE_profile_meanDownReg25.pdf",
#'                   sort_choice="min_mean",
#'                   numGenes=25, theme=1, returnData=FALSE)
#'
#'
#' /*
#'  * Calculate the variance for each gene across all three time points.
#'  * Sort the data by the highest gene-wise variance, select 30 genes
#'  * with the highest variance, and visualize them in the plot.
#'  * Save the data used to generate the plot as highVar_df.
#'  */
#' highVar_df <- de_profile_plot(res_list=myResList, filename="DE_profile_highVar30.pdf",
#'                                 sort_choice="variance", numGenes=30, theme=1, returnData=TRUE)
#'
#' }
de_profile_plot <- function(res_list, filename="de_profile_plot.pdf", sort_choice="max", specific_genes="", numGenes=50, theme=1, customLabels=FALSE, returnData=FALSE)
{
  imgHeight = 10
  imgWidth = 10

  if((typeof(res_list) != "list") || (typeof(res_list[[1]]) != "S4") || (length(res_list) <= 0))
  {
    stop("de_profile_plot() requires list type object containing DESeq result sets.")
    return(-1)
  }
  if(!exists("output_mode", envir=.DEVis_env))
  {
    stop("de_profile_plot(): Output mode not initialized.  Please run set_output_mode() first.")
    return(-1)
  }
  if(typeof(filename) != "character")
  {
    stop("filename must be a string.")
    return(-1)
  }
  if(!is.numeric(numGenes))
  {
    stop("numGenes must be an integer.")
    return(-1)
  }
  if(!exists("DE_profile_dir", envir=.DEVis_env))
  {
    stop("Directories not initialized.  Try running create_dir_struct().")
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
  if(typeof(customLabels) != "logical")
  {
    stop("customLabels must be a logical value. Please enter TRUE or FALSE.")
    return(-1)
  }
  if(typeof(returnData) != "logical")
  {
    stop("returnData must be a logical value. Please enter TRUE or FALSE.")
    return(-1)
  }
  type_options = c(".pdf",".png")
  filetype = substr(filename, nchar(filename)-4+1, nchar(filename))
  if(!(filetype %in% type_options))
  {
    stop('filename must have extension ".pdf" or ".png".')
    return(-1)
  }
  sort_options = c("max","min","variance","max_mean","min_mean")
  if(!(sort_choice %in% sort_options))
  {
    stop('sort_choice must be one of the following: "max","min","variance","max_mean","min_mean"')
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

  #Get the aggregated DE gene names.
  master_aggregate_names <- aggregate.names
  master_aggregate <- lapply(enrich, function(x) x[aggregate.names,])
  master_lfc   <- lapply(master_aggregate, function(x) x$log2FoldChange)
  master_colnames_lfc <- lapply(master_aggregate, function(x) x@case)
  master_colnames <- master_colnames_lfc

  #Finalize the master data frame.
  master_df <- data.frame(master_lfc)
  colnames(master_df) <- master_colnames
  rownames(master_df) <- rownames(master_aggregate[[1]])


  #Handle any sorting options.
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
  if(sort_choice == "variance")
  {
    varCol <- apply(master_df, 1, var)
    master_df[,"variance"] <- varCol
    master_df <- master_df[with(master_df, order(-varCol)), ]
    master_df$variance <- NULL
  }
  if(sort_choice == "max_mean")
  {
    meanCol <- apply(master_df, 1, mean)
    master_df[,"max_mean"] <- meanCol
    master_df <- master_df[with(master_df, order(-meanCol)), ]
    master_df$max_mean <- NULL
  }
  if(sort_choice == "min_mean")
  {
    meanCol <- apply(master_df, 1, mean)
    master_df[,"min_mean"] <- meanCol
    master_df <- master_df[with(master_df, order(meanCol)), ]
    master_df$min_mean <- NULL
  }

  master_df$gene <- rownames(master_df)
  master_subset <- head(master_df, numGenes)
  melted <- melt(master_subset, id.var="gene")

  #Satisfy check() that our ggplot aes variables are indeed not globals.
  gene     = NULL
  value    = NULL
  variable = NULL

  #Automatically scale the image so all the genes fit and labels are readable.
  if(numGenes > 50)
  {
    imgHeight = (numGenes / 10)
    imgWidth  = (numGenes / 10) + 5
  }

  #If custom labels are requested retrieve and update them here.
  if(customLabels == TRUE)
  {
    curLabs = as.character(levels(melted$variable))
    relabel <- as.character(melted$variable)

    #Get updated labels for each level in the current label set.
    for(i in 1:length(curLabs))
    {
      curLab <- as.character(readline(prompt=paste("Enter a custom label corresponding to ", curLabs[i], ": ", sep="")))
      sub1 <- paste("relabel[grepl('", curLabs[i], "',relabel)] <- '", curLab, "'", sep="")
      eval(parse(text = sub1))
    }
    melted$variable <- relabel
  }

  plot_out <- ggplot(data=melted, aes(x=gene, y=value, colour=variable, group = variable, stroke=2)) + geom_line() + geom_point()

  #Apply the selected theme.
  if(theme == 1)
  {
    stata_long_pal = c(stata_pal("s2color")(15), stata_pal("s1rcolor")(15))
    plot_out <- plot_out + theme_stata() + scale_fill_manual(values=stata_long_pal) + guides(col=guide_legend(ncol=length(enrich)%/%3))
  }
  if(theme == 2)
  {
    nature_pal = c(pal_npg("nrc")(10),pal_npg("nrc", alpha = 0.6)(10),pal_npg("nrc", alpha = 0.3)(10) )
    plot_out <- plot_out + theme_igray() + scale_fill_manual(values=nature_pal) + guides(col=guide_legend(ncol=length(enrich)%/%3))
  }
  if(theme == 3)
  {
    tron_pal = c(pal_tron("legacy")(7),pal_tron("legacy", alpha = 0.7)(7),pal_tron("legacy", alpha = 0.5)(7),pal_tron("legacy", alpha = 0.3)(7),pal_tron("legacy", alpha = 0.2)(2))
    plot_out <- plot_out + theme_hc(bgcolor = "darkunica") + theme(axis.text.x = element_text(colour = "white"), axis.text.y = element_text(colour = "white")) + scale_fill_manual(values=tron_pal) + guides(col=guide_legend(ncol=length(enrich)%/%3))
  }
  if(theme == 4)
  {
    gdoc_pal <- c(pal_ucscgb("default")(26), pal_ucscgb("default",alpha=.5)(4))
    plot_out <- plot_out + theme_gdocs() + scale_fill_manual(values=gdoc_pal) + guides(col=guide_legend(ncol=length(enrich)%/%3))
  }
  if(theme == 5)
  {
    d3_pal <- c(pal_d3("category20")(20), pal_d3("category10",alpha=.5)(10))
    plot_out <- plot_out + theme_solarized() + scale_fill_manual(values=d3_pal) + guides(col=guide_legend(ncol=length(enrich)%/%3))
  }
  if(theme == 6)
  {
    plot_out <- plot_out + theme_bw() + scale_colour_grey(start = 0, end = .9) + scale_fill_grey(start = 0, end = .9) + guides(col=guide_legend(ncol=length(enrich)%/%3))
  }

  plot_out <- plot_out +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size=16, face="bold"), legend.title=element_blank(), legend.text=element_text(size=12)) +
    theme(text = element_text(size=16,margin = margin(t = 0, r = 10, b = 0, l = 0))) +
    labs(x="", y=expression(paste("Log2 Fold-Change"))) +
    theme(axis.text=element_text(size=10, face="bold")) +
    guides(shape = guide_legend(override.aes = list(size=5))) +
    guides(colour = guide_legend(override.aes = list(size=5)))

  #Write the plot.
  setwd(.DEVis_env$DE_profile_dir)
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
    print(plot_out)
    dev.off()
  }
  if(.DEVis_env$output_mode == "screen" | .DEVis_env$output_mode == "both")
  {
    print(plot_out)
  }
  setwd(.DEVis_env$working_dir)

  if(returnData)
  {
    return(master_subset)
  }
}


#' Visualize density plots of fold-change or significance values for aggregated data sets.
#'
#' This function plots log2 fold-change or adjusted p-values for all differentially expressed genes for each contrast in
#' a result set.  Data aggregation across a series of result sets entails that not every gene in an aggregated data
#' set will necessarily be differentially expressed, depending on the aggregation method.  For instance, a gene that is differentially
#' expressed in a day1 contrast that is not significant at day2 will be included in a union based aggregation of genes for
#' day1 and day2.  Visualizing the density of fold-change or p-values can reveal to what extent aggregation was consistent
#' across conditions and can inform decision making in data aggregation.
#' @param res_list A list of DESeq result sets. Results can be calculated individually using DESeq's results() function.
#' Lists of results can be created by creating a list(result1, result2 ... result_N).
#' @param filename Filename for output plot.  Valid extensions are ".pdf" and ".png".
#' File generation can be turned off using set_output_mode("screen").
#' Output will be written to the /DE/density_plots/ directory.
#' @param type The type of data to be displayed in this plot.
#' Valid selections are "lfc" (log2foldChange) and "pval".
#' @param method The method for computing overlaying results.  Valid selections are "union" or "intersection".
#' Union merges data for all result sets for genes that are differentially expressed in at least 1 result set.
#' Intersection merges data for genes that only are differentially expressed in all result sets provided.
#' @param returnData Boolean.  Determines if this visualization should return data used to generate the visualization. Default=FALSE.
#' @return If returnData is true, this function will return a data frame for aggregated differentially expressed genes
#' containing gene names and log2 fold-change or adjusted p-values relative to the experimental control condition.
#' @keywords DE density fold-change aggregate p-value visualization
#' @export
#' @examples
#' \dontrun{
#'
#'  #This example assumes an experimenal design of ~Condition_Time.
#'
#' #Prepare a result list.
#' res.day1 <- results(dds, contrast=c("Condition_Time", "day1_disease", "day1_control"))
#' res.day2 <- results(dds, contrast=c("Condition_Time", "day2_disease", "day2_control"))
#' res.day3 <- results(dds, contrast=c("Condition_Time", "day3_disease", "day3_control"))
#' myResList <- list(res.day1, res.day2, res.day3)
#'
#' /*
#'  * Aggregate data for all contrasts in the result list using union aggregation.
#'  * Display the density plot of p-values for the aggregated data.
#'  */
#' de_density_plot(res_list=myResList, filename="DE_density_union_pval.pdf",
#'                  type="pval", method="union", returnData=FALSE)
#'
#' /*
#'  * Aggregate data for all contrasts in the result list using intersection aggregation.
#'  * Display the density plot of log fold-change values for the aggregated data.
#'  * Store the aggregate data as DE_lfc_intersect_df.
#'  */
#' DE_lfc_intersect_df <- de_density_plot(res_list=myResList,
#'                                         filename="DE_density_union_pval.pdf",
#'                                         type="lfc", method="intersection",
#'                                         returnData=TRUE)
#'
#' }
de_density_plot <- function(res_list, filename="de_density_plot.pdf", type="pval", method="union", returnData=FALSE)
{
  imgHeight = 10
  imgWidth = 10

  if((typeof(res_list) != "list") || (typeof(res_list[[1]]) != "S4") || (length(res_list) <= 0))
  {
    stop("de_density_plot() requires list type object containing DESeq result sets.")
    return(-1)
  }
  if(!exists("output_mode", envir=.DEVis_env))
  {
    stop("de_density_plot(): Output mode not initialized.  Please run set_output_mode() first.")
    return(-1)
  }
  if(typeof(filename) != "character")
  {
    stop("filename must be a string.")
    return(-1)
  }
  if(typeof(returnData) != "logical")
  {
    stop("returnData must be a logical value. Please enter TRUE or FALSE.")
    return(-1)
  }
  if(!exists("DE_density_dir", envir=.DEVis_env))
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
  type_options = c("lfc","pval")
  if(!(type %in% type_options))
  {
    stop('type must be one of the following: "lfc","pval"')
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
    master_colnames_lfc <- lapply(master_union, function(x) x@case)
    master_colnames <- master_colnames_lfc
    master_padj  <- lapply(master_union, function(x) x$padj)
    master_colnames_padj <- lapply(master_union, function(x) paste(x@case, ".padj", sep=""))
  }
  if(method=="intersection")
  {
    intersect.names = Reduce(intersect, res.names)
    master_intersect_names <- intersect.names
    master_intersect <- lapply(enrich, function(x) x[intersect.names,])
    master_lfc   <- lapply(master_intersect, function(x) x$log2FoldChange)
    master_colnames_lfc <- lapply(master_intersect, function(x) x@case)
    master_colnames <- master_colnames_lfc
    master_padj  <- lapply(master_intersect, function(x) x$padj)
    master_colnames_padj <- lapply(master_intersect, function(x) paste(x@case, ".padj", sep=""))
  }

  #Finalize the master data frame.
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

  master_df$gene <- rownames(master_df)
  melted <- melt(master_df, id.var="gene")

  #Satisfy check() that our ggplot aes variables are indeed not globals.
  value    = NULL
  variable = NULL

  plot_out <- ggplot(melted, aes(value, fill = variable)) +
    geom_density(position = "stack") +
    scale_x_continuous(breaks=seq(0,1,.05)) +
    theme(legend.text = element_text(size=12)) +
    theme(axis.text = element_text(size=10, face="bold")) +
    theme(legend.title = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(text = element_text(size=12,margin = margin(t = 0, r = 10, b = 0, l = 0)))


  #Write the plot to file.
  setwd(.DEVis_env$DE_density_dir)
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
    suppressWarnings({
      print(plot_out)
    })
    dev.off()
  }
  if(.DEVis_env$output_mode == "screen" | .DEVis_env$output_mode == "both")
  {
    suppressWarnings({
      print(plot_out)
    })
  }
  setwd(.DEVis_env$working_dir)

  if(returnData)
  {
    return(master_df)
  }

}
