#' Visualize the expression of a specific gene with regard to metadata grouping.
#'
#' This function generates a box plot to display the expression of an individual gene with regard to a
#' specified grouping that can be based on any data that exists in the targets file. For example, a plot could
#' be created to view the expression of a gene as a function of different time points or experimental conditions.
#' @param filename Filename for output plot.  Valid extensions are ".pdf" and ".png".
#' File generation can be turned off using set_output_mode("screen"). Output will be written to the /geneplots/ directory.
#' @param gene_name The name of the gene to create boxplot for. Must match a rowname in count data. String.
#' @param groupBy The group from target data that should be used to split data. Must match a rowname in count data.
#' I.E. "Timepoint" or "Infection". String.  If only two groups are present in this variable a wilcox test will be performed
#' between the two groups and the p-value will be displayed on the plot as well.
#' @param theme Theme for the layout and color scheme for the plot.  Valid selections are integers between 1-6.
#' @param returnData Boolean. Determines if this visualization should return data used to generate the visualization. Default=FALSE.
#' @return If returnData is true, this function will return the long-form table of expression containing sample names,
#' categorical grouping, and sample IDs.
#' @keywords expression DE boxplot gene visualization
#' @export
#' @examples
#' \dontrun{
#'
#' #Plot the CAPS12 gene for each time point.
#' plot_gene(filename="CAPS12_time_plot.pdf", gene_name="CAPS12",
#'            groupBy="Time", theme=1, returnData=FALSE)
#'
#' #Plot the METTL25 gene for each "response" group.  Store the long-form data table.
#' mettl25_dat <- plot_gene(filename="METTL25_time_plot.pdf", gene_name="METTL25",
#'                           groupBy="response", theme=2, returnData=TRUE)
#'
#' }
plot_gene <- function(filename="gene_plot.pdf", gene_name, groupBy, theme=1, returnData=FALSE)
{
  imgWidth  = 10
  imgHeight = 10

  #Validate the data we want exists and was processed correctly.
  if(!exists("dds_glbl", envir=.DEVis_env))
  {
    stop("plot_gene() requires that prep_dds_from_data() be run first.")
    return(-1)
  }
  if(!exists("output_mode", envir=.DEVis_env))
  {
    stop("de_heat(): Output mode not initialized.  Please run set_output_mode() first.")
    return(-1)
  }
  if(typeof(filename) != "character")
  {
    stop("filename must be a string.")
    return(-1)
  }
  if(typeof(groupBy) != "character")
  {
    stop("group variable must be a string.")
    return(-1)
  }
  if(typeof(gene_name) != "character")
  {
    stop("group variable must be a string.")
    return(-1)
  }
  if(typeof(returnData) != "logical")
  {
    stop("returnData must be a logical value. Please enter TRUE or FALSE.")
    return(-1)
  }
  if(!exists("gene_plot_dir", envir=.DEVis_env))
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
  if(!(groupBy %in% colnames(colData(.DEVis_env$stabilized_data))))
  {
    print("plot_gene() requires grouping to correspond to colData annotation.")
    print("Possible values are: ")
    print(colnames(colData(.DEVis_env$stabilized_data)))
    stop("plot_gene(): invalid group variable selection.")
    return(-1)
  }
  if(!(gene_name %in% rownames(.DEVis_env$dds_glbl)))
  {
    print(gene_name)
    stop("This gene does not exist in data. Please enter a valid gene ID.")
    return(-1)
  }
  if(theme < 1 || theme > 6)
  {
    stop("theme must be a value between 1 and 6.")
    return(-1)
  }

  geneData <- plotCounts(.DEVis_env$dds_glbl, gene_name, groupBy, returnData=TRUE)
  count_data<-assay(.DEVis_env$stabilized_data)

  #Check number of groups.
  target_data <- .DEVis_env$tgt_dat
  ux <- unique(subset(target_data, select=c(groupBy)))
  sub1 <- paste("ux$", groupBy, sep="")
  groupBy_factors <- eval(parse(text = sub1))

  #If we're comparing two groups only, do a wilcox test and add to the plot.
  if (length(groupBy_factors) == 2){

    sub1 <- paste("target_data$", groupBy, sep="")
    target_data_groupBy_colname <-  eval(parse(text = sub1))
    x_samples<-target_data$targetID[which(target_data_groupBy_colname==groupBy_factors[1])]
    y_samples<-target_data$targetID[which(target_data_groupBy_colname==groupBy_factors[2])]
    x_column<-count_data[,x_samples]
    y_column<-count_data[,y_samples]
    x<-x_column[gene_name,]
    y<-y_column[gene_name,]

    pv<-wilcox.test(x,y)
  }

  sub1 <- paste("ggplot(geneData, aes(x=", groupBy, ", y=log(count), fill=", groupBy, "))", sep="")
  plot_out <- eval(parse(text = sub1))

  plot_out <- plot_out + geom_boxplot() + scale_y_continuous(name = "Expression") + ggtitle(paste("Gene Expression: ", gene_name, sep=""))

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
    plot_out <- plot_out + theme_bw() + scale_colour_grey(start = 0, end = .9) + scale_fill_grey(start = 0, end = .9)
  }

  if (length(groupBy_factors) == 2){
    plot_out <- plot_out + annotate("text", label=paste("p-value=",signif(pv$p.value, digits = 3),sep=""), x=-Inf, y=Inf, hjust=-0.2, vjust=1.5, size=4)
  }

  plot_out <- plot_out +
              theme(legend.title = element_blank()) +
              theme(legend.text = element_text(size=12)) +
              theme(text = element_text(size=16,margin = margin(t = 0, r = 2, b = 1, l = 2))) +
              theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(size=16, face="bold",hjust = 0))


  #Write the plot to file.
  setwd(.DEVis_env$gene_plot_dir)
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
    return(geneData)
  }

}


#' Visualize overall data set as a function of a metadata grouping.
#'
#' This function plots groupwise expression as with regard to a specified target data grouping. This function can be
#' used to visualize the overall data set according to metadata grouping.  For instance, overall expression at different
#' timepoints can be viewed, or infected vs control expression levels can be plotted.  This visualization is also
#' ideal for examining the impact of normalization or filtering.
#' @param filename Filename for output plot.  Valid extensions are ".pdf" and ".png".
#' File generation can be turned off using set_output_mode("screen"). Output will be written to the /group_stats/ directory.
#' @param id_field The unique ID field for this data set.  Should be the field that corresponds to count column names and target rownames.
#' @param groupBy The field by which to group samples in the boxplot. I.E. "time"
#' @param normalized Toggles between whether normalized or non-normalized data should be used. Boolean. Default=TRUE.
#' @param theme Theme for the layout and color scheme for the plot.  Valid selections are integers between 1-6.
#' @param returnData Boolean. Determines if this visualization should return data used to generate the visualization. Default=FALSE.
#' @return If returnData is true, this function will return the long-form table of expression containing
#' categorical grouping, and sample IDs, log fold-change, and gene name.
#' @keywords expression normalization boxplot group distribution visualization
#' @export
#' @examples
#' \dontrun{
#'
#' #Plot the overall data set expression levels grouped by time.
#' plot_group_stats("group_stats_plot.pdf", "SampleID", groupBy="Time",
#'                    normalized=TRUE, theme=1, returnData=FALSE)
#'
#' }
plot_group_stats <- function(filename="group_stats_plot.pdf", id_field, groupBy, normalized=TRUE, theme=1, returnData=FALSE)
{
  imgWidth  = 10
  imgHeight = 10

  #Validate the data we want exists and was processed correctly.
  if(!exists("stabilized_data", envir=.DEVis_env))
  {
    stop("plot_group_stats() requires that prep_dds_from_data() be run first.")
    return(-1)
  }
  if(!exists("tgt_dat", envir=.DEVis_env))
  {
    stop("Target data not present.  Please run prep_targets() first.")
    return(-1)
  }
  if(!exists("output_mode", envir=.DEVis_env))
  {
    stop("plot_group_stats(): Output mode not initialized.  Please run set_output_mode() first.")
    return(-1)
  }
  if(typeof(filename) != "character")
  {
    stop("filename must be a string.")
    return(-1)
  }
  if(typeof(groupBy) != "character")
  {
    stop("group variable must be a string.")
    return(-1)
  }
  if(typeof(returnData) != "logical")
  {
    stop("returnData must be a logical value. Please enter TRUE or FALSE.")
    return(-1)
  }
  if(typeof(id_field) != "character")
  {
    stop("group variable must be a string.")
    return(-1)
  }
  if(!exists("group_stat_dir", envir=.DEVis_env))
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
  if(!(groupBy %in% colnames(colData(.DEVis_env$stabilized_data))))
  {
    print("plot_group_stats() requires grouping to correspond to colData annotation.")
    print("Possible values are: ")
    print(colnames(colData(.DEVis_env$stabilized_data)))
    stop("plot_group_stats(): invalid group variable selection.")
    return(-1)
  }
  if(!(id_field %in% colnames(colData(.DEVis_env$stabilized_data))))
  {
    print("plot_group_stats() requires id_field to correspond to colData annotation.")
    print("Possible values are: ")
    print(colnames(colData(.DEVis_env$stabilized_data)))
    stop("plot_group_stats(): invalid group variable selection.")
    return(-1)
  }
  if(theme < 1 || theme > 6)
  {
    stop("theme must be a value between 1 and 6.")
    return(-1)
  }

  if(normalized)
  {
    melted_data <- melt(assay(.DEVis_env$stabilized_data))
    tgt_dat <- .DEVis_env$tgt_dat
    sub1 <- paste("tgt_dat$", groupBy, "[match(melted_data$Var2,tgt_dat$", id_field, ")]", sep="")
    melted_data$group <- group_levels <- eval(parse(text = sub1))
  }
  else
  {
    melted_data <- melt(log2(assay(.DEVis_env$dds_glbl)))
    tgt_dat <- .DEVis_env$tgt_dat
    sub1 <- paste("tgt_dat$", groupBy, "[match(melted_data$Var2,tgt_dat$", id_field, ")]", sep="")
    melted_data$group <- group_levels <- eval(parse(text = sub1))
  }

  #Satisfy check() that our ggplot aes variables are indeed not globals.
  value    = NULL
  variable = NULL
  group    = NULL

  plot_out <- ggplot(melted_data,aes(x=group,y=value,fill=group)) + geom_boxplot() + scale_y_continuous(name = "Log(Expression)") + theme(legend.title = element_blank()) + theme(legend.text = element_text(size=8))  + theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(size=16, face="bold"))

  #Apply the selected theme.
  if(theme == 1)
  {
    stata_long_pal = c(stata_pal("s2color")(15), stata_pal("s1rcolor")(15))
    plot_out <- plot_out + theme_stata() + scale_fill_manual(values=stata_long_pal)
  }
  if(theme == 2)
  {
    nature_pal = c(pal_npg("nrc")(10),pal_npg("nrc", alpha = 0.6)(10),pal_npg("nrc", alpha = 0.3)(15) )
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
    plot_out <- plot_out + theme_bw() + scale_colour_grey(start = 0, end = .9) + scale_fill_grey(start = 0, end = .9)
  }

  plot_out <- plot_out +
              theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
              theme(legend.text = element_text(size=12)) +
              theme(axis.text = element_text(size=12, face="bold")) +
              theme(text = element_text(size=16,margin = margin(t = 0, r = 2, b = 1, l = 2))) +
              labs(x="", y="log(Expression)") +
              theme(legend.title=element_blank())

  #Write the plot to file.
  setwd(.DEVis_env$group_stat_dir)
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
    return(melted_data)
  }

}


#' Visualize differentially expressed genes as a function of experimental design.
#'
#' This function plots groupwise expression for all differentially expressed genes
#' with regard to the experimental design.  For example, for a differential expression design of
#' ~Condition_Time, boxplots for differentially expressed genes for each condition_time case in
#' the targets metadata will be produced.
#' @param res_list A list of DESeq result sets. Results can be calculated individually using DESeq's results() function.
#' Lists of results can be created by creating a list(result1, result2 ... result_N).
#' @param filename Filename for output plot.  Valid extensions are ".pdf" and ".png".
#' File generation can be turned off using set_output_mode("screen"). Output will be written to the /DE/boxplot/ directory.
#' @param theme Theme for the layout and color scheme for the plot.  Valid selections are integers between 1-6.
#' @param returnData Boolean. Determines if this visualization should return data used to generate the visualization. Default=FALSE.
#' @return If returnData is true, this function will return the long-form table of expression containing
#' categorical grouping, and sample IDs, and log2 fold-change.
#' @keywords expression DE boxplot distribution visualization
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
#' #Plot differentially expressed genes for the aggregate result sets.
#' de_boxplot(res_list=myResList, filename="DE_condition_time_boxplot.pdf",
#'             theme=2, returnData=FALSE)
#'
#' }
de_boxplot <- function(res_list, filename="de_boxplot.pdf", theme=1, returnData=FALSE)
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
    stop("de_heat(): Output mode not initialized.  Please run set_output_mode() first.")
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
  if(!exists("DE_boxplot_dir", envir=.DEVis_env))
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

  plot_out <- ggplot(filteredList, aes(x=variable,y=value,fill=variable)) + geom_boxplot() + scale_y_discrete(name = "Log2 Fold-Change") + theme(legend.title = element_blank()) + theme(legend.text = element_text(size=6)) + theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(size=16, face="bold"))

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
    plot_out <- plot_out + theme_bw() + scale_colour_grey(start = 0, end = .9) + scale_fill_grey(start = 0, end = .9)
  }

  plot_out <- plot_out +
              theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
              theme(legend.text = element_text(size=14)) +
              theme(legend.title=element_blank()) +
              theme(axis.text=element_text(size=14, face="bold")) +
              theme(text = element_text(size=22,margin = margin(t = 0, r = 2, b = 1, l = 2))) +
              labs(x="", y="Log2 Fold-Change") +
              theme(legend.title=element_blank())

  #Write the plot to file.
  setwd(.DEVis_env$DE_boxplot_dir)
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
    return(filteredList)
  }
}
