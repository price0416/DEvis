#' Visualize differentially expressed gene counts as a stacked barplot.
#'
#' This function generates a stacked barplot representing differentially expressed gene counts for a result sets,
#' distinguishing up and down regulated genes across one or more result sets.
#' This function will produce a stacked barplot and an output file containing count information.
#' Count data for the provided result sets will also be returned by this function.
#' @param res_list A list of DESeq result sets created with DESeq2::results(). I.E: list(res1, res2, ..., resN).
#' @param filename Output file destination. String. Valid extensions are .pdf and .png.  The data file accompanying this
#' plot file will have the same name, but will be output as a tab-delimited text file. Output will be written to
#' the /DE/counts/ directory. Alternatively, file generation can be turned off using set_output_mode("screen").
#' @param lfc_filter Boolean.  Determines if DE counts should be displayed with or without the application of the log2FoldChange filter.
#' @param customLabels If customLabels is set to TRUE, the user will be prompted to provide a custom label for each label on the x-axis.
#' @param theme The color and design scheme for the output plot.  Valid selections are integers between 1-6.
#' @param returnData Boolean.  Determines if this visualization should return data used to generate the visualization. Default=FALSE.
#' @return A data frame containing count information for both up and down regulated gene counts for each result set
#' provided in res_list.  Output files will be written to /DE/counts/.
#' @keywords counts summary DE visualization
#' @seealso \code{\link{create_dir_struct}}, \code{\link{set_output_mode}}
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
#' #Visualize count data for the result set and save the results.
#' de_counts(res_set=myResList, filename="DE_counts.png", theme=2)
#'
#' }
de_counts <- function(res_list, filename="de_count_barplot.pdf", lfc_filter=FALSE, customLabels=FALSE, theme=1, returnData=FALSE)
{
  imgWidth  = 8
  imgHeight = 8

  #Validate expected values.
  if((typeof(res_list) != "list") || (typeof(res_list[[1]]) != "S4") || (length(res_list) <= 0))
  {
    stop("de_counts() requires list type object containing DESeq result sets.")
    return(-1)
  }
  if(!exists("output_mode", envir=.DEVis_env))
  {
    stop("de_counts(): Output mode not initialized.  Please run set_output_mode() first.")
    return(-1)
  }
  if(typeof(filename) != "character")
  {
    stop("filename must be a string.")
    return(-1)
  }
  if(typeof(lfc_filter) != "logical")
  {
    stop("lfc_filter must be a logical value. Please enter TRUE or FALSE.")
    return(-1)
  }
  if(typeof(returnData) != "logical")
  {
    stop("returnData must be a logical value. Please enter TRUE or FALSE.")
    return(-1)
  }
  if(typeof(customLabels) != "logical")
  {
    stop("customLabels must be a logical value. Please enter TRUE or FALSE.")
    return(-1)
  }
  if(!exists("DE_count_dir", envir=.DEVis_env))
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
  prefix <- substr(filename, 1, nchar(filename)-4)
  data_outfile <- paste(prefix,".txt", sep="")
  if(theme < 1 || theme > 6)
  {
    stop("theme must be a value between 1 and 6.")
    return(-1)
  }

  #Enrich the data.
  if(lfc_filter)
  {
    enriched <- unlist(lapply(res_list, enrich_res, TRUE))
  }
  else
  {
    enriched <- unlist(lapply(res_list, enrich_res))
  }
  enrich <- unlist(enriched)

  #Get the number of up / down DE genes.
  numUp_de   <- lapply(enrich, function(x) x@numUp)
  numDown_de <- lapply(enrich, function(x) x@numDown)

  #Stitch the up/down counts into a common vector.
  updown_list <- vector()
  for (i in 1:length(numUp_de))
  {
    updown_list <- c(updown_list, numUp_de[i])
    updown_list <- c(updown_list, numDown_de[i])
  }
  updown_list <- as.numeric(unlist(updown_list))
  updown_list <- as.numeric(ifelse(is.na(updown_list), 1, updown_list))

  #Prepare the name data.
  case_de      <- unlist(lapply(enrich, function(x) x@case))
  updown_names <- rep(case_de,each=2)

  #Prepare the direction groups.
  stock_direction <- c("Up", "Down")
  upDown_direction <- rep(stock_direction, length(updown_names)/2)

  #Aggregate all of the pieces and melt it down.
  upDown_df = data.frame(updown_list, updown_names, upDown_direction)
  suppressMessages({
    upDown_melt <- melt(upDown_df)
  })

  #Satisfy check() that our ggplot aes variables are indeed not globals.
  value = NULL

  #If custom labels are requested retrieve and update them here.
  if(customLabels == TRUE)
  {
    curLabs = as.character(levels(upDown_melt$updown_names))
    updown_relabel <- as.character(upDown_melt$updown_names)

    #Get updated labels for each level in the current label set.
    for(i in 1:length(curLabs))
    {
      curLab <- as.character(readline(prompt=paste("Enter a custom label corresponding to ", curLabs[i], ": ", sep="")))
      sub1 <- paste("updown_relabel[grepl('", curLabs[i], "',updown_relabel)] <- '", curLab, "'", sep="")
      eval(parse(text = sub1))
    }
    upDown_melt$updown_names <- updown_relabel
  }

  #Generate the plot.
  plot_out <- ggplot(upDown_melt, aes(x = factor(updown_names, levels=as.character(unique(updown_names)), ordered=TRUE), y = value, fill = upDown_direction)) + geom_bar(stat = "identity")

  if(theme == 1)
  {
    group.colors <- c(Up = "#e56262", Down = "#1bb2c6")
    plot_out <- plot_out + theme_stata() + scale_colour_stata()
  }
  if(theme == 2)
  {
    group.colors <- c(Up = "#FF5640", Down = "#04819E")
    plot_out <- plot_out + theme_solarized() + scale_colour_solarized("blue")
  }
  if(theme == 3)
  {
    group.colors <- c(Up = "#FFF500", Down = "#604BD8")
    plot_out <- plot_out + theme_solarized_2(light = FALSE) + scale_colour_solarized("blue")
  }
  if(theme == 4)
  {
    group.colors <- c(Up = "#0b0c0c", Down = "#7d8484")
    plot_out <- plot_out + scale_fill_excel() + theme_excel()
  }
  if(theme == 5)
  {
    group.colors <- c(Up = "#ed2121", Down = "#0451ea")
    plot_out <- plot_out + theme_wsj() + scale_colour_wsj("colors6", "")
  }
  if(theme == 6)
  {
    group.colors <- c(Up = "gray25", Down = "gray75")
    plot_out <- plot_out + theme_bw() + scale_colour_grey(start = 0, end = .9) + scale_fill_grey(start = 0, end = .9)
  }

  plot_out <- plot_out +
              theme(axis.text.x=element_text(angle=90,hjust=-0.2,vjust=0.5), legend.title=element_blank()) +
              ggtitle("") + labs(x="", y=expression(paste("DE Genes"), title="Expression")) +
              theme(text = element_text(size=16,margin = margin(t = 0, r = 10, b = 0, l = 0))) +
              scale_fill_manual(values=group.colors)

  setwd(.DEVis_env$DE_count_dir)
  if(.DEVis_env$output_mode == "both" | .DEVis_env$output_mode == "file")
  {
    #Save the plot.
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

    #Write the output data.
    write.table(upDown_df, data_outfile, sep="\t")
  }
  if(.DEVis_env$output_mode == "screen" | .DEVis_env$output_mode == "both")
  {
    print(plot_out)
  }

  setwd(.DEVis_env$working_dir)

  if(returnData)
  {
    return(upDown_df)
  }
}
