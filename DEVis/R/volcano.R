#' Visualize the expression & significance of differentially expressed genes.
#'
#' This function generates a volcano plot that displays the relation between expression and p-value for all DE
#' genes in a result set.
#' @param res_list A list of DESeq result sets. Results can be calculated individually using DESeq's results() function.
#' Lists of results can be created by creating a list(result1, result2 ... result_N).
#' @param filename Filename for output plot.  Valid extensions are ".pdf" and ".png".
#' File generation can be turned off using set_output_mode("screen"). Output will be written to the /DE/volcano/ directory.
#' @param lfc_thresh Log2 fold change threshold for highlighting genes on the volcano plot.
#' @param strict_scale Boolean. If this is true scales for all samples will be forced to be identical, otherwise scales will
#' vary depending on the dispersion of the volcano plot for each individual sample.
#' @param num_columns Number of columns to use in the grid layout. Default=3.
#' @param theme Theme for the layout and color scheme for the plot.  Valid selections are integers between 1-6.
#' @param returnData Boolean. Determines if this visualization should return data used to generate the visualization. Default=FALSE.
#' @return If returnData is true, this function will return the long-form table of expression containing sample names,
#' categorical grouping, and sample IDs.
#' @keywords expression DE volcano significance p-value visualization
#' @export
#' @examples
#' \dontrun{
#'
#' #Make a volcano plot highlighting DE genes above 1.5 LFC threshold.
#' de_volcano(result_list, filename="volcano.pdf", lfc_thresh=1.5, 
#'            strict_scale=TRUE, theme=1, returnData=FALSE)
#'
#' }
de_volcano <- function(res_list, filename="volcano_plot.pdf", lfc_thresh=-1, strict_scale=TRUE, num_columns=3, theme=1, returnData=FALSE)
{
  imgHeight = 8
  imgWidth = 8
  
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
  if(typeof(filename) != "character")
  {
    stop("filename must be a string.")
    return(-1)
  }
  if(!is.numeric(lfc_thresh))
  {
    stop("lfc_threshold must be numeric.")
    return(-1)
  }
  if(!is.numeric(num_columns))
  {
    stop("num_columns must be numeric.")
    return(-1)
  }
  if(typeof(returnData) != "logical")
  {
    stop("returnData must be a logical value. Please enter TRUE or FALSE.")
    return(-1)
  }
  if(typeof(strict_scale) != "logical")
  {
    stop("strict_scale must be a logical value. Please enter TRUE or FALSE.")
    return(-1)
  }
  if(!exists("DE_volcano_dir", envir=.DEVis_env))
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
  if(!exists("lfc_cutoff", envir=.DEVis_env) & lfc_thresh == -1)
  {
    stop("de_volcano(): No LFC threshold initialized. Please run init_cutoffs() or specify a desired LFC threshold.")
    return(-1)
  }
  #If they don't specify a threshold use the initialized one.
  if(lfc_thresh == -1)
  {
    lfc_thresh  <- .DEVis_env$lfc_cutoff
  }
  if(lfc_thresh < -1)
  {
    stop("Please provide a positive value for lfc_thresh.")
    return(-1)
  }
  
  #Get LFC and P-value info for the provided result list.
  lfc_data  <- get_de_data(res_list)
  pval_data <- get_de_data(res_list,type="pval")
  num_samples <- dim(lfc_data)[2]
  
  #Build basic versions of the plots without any themes for each sample.
  max_x = 0
  max_y = 0
  x_range_list = c()
  plot_id_list = c()
  return_data  = list()
  for (i in 1:num_samples)
  {
    curLab =  colnames(lfc_data)[i]
    
    #Prepare data for the current sample (lfc and padj data).
    sub1 <- paste("cur_df <- data.frame(lfc_data$\"",as.character(curLab),"\"[rownames(lfc_data) %in% rownames(pval_data)], -log10(pval_data$\"",as.character(curLab),".padj\"))", sep="")
    eval(parse(text = sub1))

    rownames(cur_df) <- rownames(pval_data)
    colnames(cur_df) <- c(paste(curLab,".lfc",sep=""), paste(curLab,".padj",sep=""))
    
    #Make sure first character of field is not a number, R hates that.
    suppressWarnings({
      if(!is.na(as.numeric(substring(curLab, 1, 1))))
      {
        curLab = paste("x", curLab, sep="")
        colnames(cur_df) <- c(paste(curLab,".lfc",sep=""), paste(curLab,".padj",sep=""))
      }
    })
    
    #Handle any NA values that might exist.
    sub1 <- paste("cur_df <- cur_df[!is.na(cur_df$\"",as.character(curLab),".padj\"),]", sep="")
    eval(parse(text = sub1))
    
    #Keep track of the range so we can scale all of these appropriately at the end.
    sub1 <- paste("max_x = max(abs(cur_df$",curLab,".lfc),  max_x)", sep="")
    eval(parse(text = sub1))
    sub1 <- paste("max_y = max(abs(cur_df$",curLab,".padj), max_y)", sep="")
    eval(parse(text = sub1))
    cur_x = 0
    if(strict_scale == FALSE)
    {
      sub1 <- paste("cur_x = max(abs(cur_df$",curLab,".lfc))", sep="")
      eval(parse(text = sub1))
      x_range_list = c(x_range_list, cur_x)
    }
                      
    #Identify the LFC threshold for coloring points for each gene.
    threshold = c()
    sub1 <- paste("threshold = ifelse(cur_df$\"",curLab,".lfc\" >= ", as.character(lfc_thresh), ", 'A', ifelse(cur_df$\"",curLab,".lfc\" <= -", as.character(lfc_thresh),", 'A', 'B'))", sep="")
    eval(parse(text = sub1))
    cur_df$thresh = threshold    
    
    #Set up the basic plots for each sample.
    sub1 <- paste("vp",as.character(i)," <- ggplot(cur_df, aes(x=",curLab,".lfc, y=",curLab,".padj)) + geom_point(aes(colour = thresh), size=2.0)", sep="")
    eval(parse(text = sub1))
    
    #Keep track of our IDs.
    sub1 <- paste("plot_id_list <- c(plot_id_list,\"vp",as.character(i),"\")", sep="")
    eval(parse(text = sub1))
    
    #Keep track of things if they are requesting the data be returned.
    if(returnData == TRUE)
    {
      return_data[i] = cur_df
    }
  }

  #Round up to the nearest integer for the x/y ranges.
  max_x = ceiling(max_x)
  max_y = ceiling(max_y)
  
  #Determine the number of rows/columns we will need to draw all of these plots.
  num_cols = num_columns
  num_rows = ceiling(num_samples/num_cols)
  if(num_samples < num_cols)
  {
    num_cols = num_samples
    num_rows = 1
  }
  imgHeight = (imgHeight * num_rows)/2
  
  #Now that all of our plots are generated in skeleton form, go through and add bells and whistles.
  for (i in 1:num_samples)
  {
    curLab =  colnames(lfc_data)[i]
    
    #See if they are using the strict scale and adjust accordingly.
    if(strict_scale == TRUE)
    {
      sub1 <- paste("vp",i, " <- vp",i, " + xlim(-", max_x, ",", max_x, ") + ylim(0,", max_y, ")", sep="")
      eval(parse(text = sub1))
    }
    else
    {
      #At least center the x-axis at 0.
      sub1 <- paste("vp",i, " <- vp",i, " + xlim(-", x_range_list[i], ",", x_range_list[i], ")", sep="")
      eval(parse(text = sub1))
    }
    
    #Apply the selected theme.
    suppressWarnings({
      if(theme == 1)
      {
        sub1 <- paste("vp",i, " <- vp",i, " + theme_stata() + scale_colour_manual(values=c('A'= '#DC0000FF', 'B'='black'))", sep="")
        eval(parse(text = sub1))
      }
      if(theme == 2)
      {
        sub1 <- paste("vp",i, " <- vp",i, " + theme_igray() + scale_colour_manual(values=c('A'= '#00ffff', 'B'='black'))", sep="")
        eval(parse(text = sub1))
      }
      if(theme == 3)
      {
        sub1 <- paste("vp",i, " <- vp",i, " + theme_hc(bgcolor = \"darkunica\") + theme(axis.text.x = element_text(colour = \"white\"), axis.text.y = element_text(colour = \"white\")) + scale_colour_manual(values=c('A'= '#ffbf00', 'B'='#d6d6d6'))", sep="")
        eval(parse(text = sub1))
      }
      if(theme == 4)
      {
        sub1 <- paste("vp",i, " <- vp",i, " + theme_gdocs() + scale_colour_manual(values=c('A'= 'orange', 'B'='black'))", sep="")
        eval(parse(text = sub1))
      }
      if(theme == 5)
      {
        sub1 <- paste("vp",i, " <- vp",i, " + theme_solarized() + scale_colour_manual(values=c('A'= 'green', 'B'='black'))", sep="")
        eval(parse(text = sub1))
      }
      if(theme == 6)
      {
        sub1 <- paste("vp",i, " <- vp",i, " + theme_bw() + scale_colour_manual(values=c('A'= 'grey', 'B'='black'))", sep="")
        eval(parse(text = sub1))
      }
    })
    
    #Put on the finishing touches.
    sub1 <- paste("vp",i, " <- vp",i, " + theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.margin = unit(c(1,1,1,1), \"lines\"), plot.title = element_text(size=8, face=\"bold\"), legend.position=\"none\") +
                  theme(text = element_text(size=8,margin = margin(t = 0, r = 3, b = 0, l = 0))) +
                  labs(title=\"",curLab,"\", x=\"Log2 Fold Change\", y=expression(paste(\"-log10(P-value)\"))) +
                  theme(axis.text=element_text(size=7, face=\"bold\"))",sep="")
    eval(parse(text = sub1))
  }
  
  #Prepare the grid plotting command.
  plot_str = ""
  sub1 <- paste("grid.arrange(",sep="")
  for (i in 1:num_samples)
  {
    plot_str <- paste(plot_str, ",", plot_id_list[i],sep="")
  }
  plot_str <- substr(plot_str,2,9999)   #Dirty hack to adjust commas.
  plot_str <- paste(plot_str,",",sep="")
  sub3 <- paste(sub1, plot_str, "ncol = ", num_cols, ", nrow = ", num_rows, ")")

  #Write the plot.
  setwd(.DEVis_env$DE_volcano_dir)
  if(.DEVis_env$output_mode == "screen" | .DEVis_env$output_mode == "both")
  {
    eval(parse(text = sub3))
  }
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
    eval(parse(text = sub3))
    dev.off()
  }
  
  setwd(.DEVis_env$working_dir)
  
  #Handle returning data if user requests.
  if(returnData)
  {
    return(return_data)
  }
}