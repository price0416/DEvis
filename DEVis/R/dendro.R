#' Create dendrograms based on hierarchical clustering.
#'
#' This function plots a dendrogram for the dataset as a whole using hierarchical clustering.  Metadata
#' can be applied to the dendrogram to examine how data cluster by different conditions. By examining
#' how data cluster, it is possible to identify patterns in the data and batch effects. For example,
#' creating a dendrogram that is grouped by the "case" condition in a case vs control experiment,
#' it would be possible to identify a strong effect in the overall dataset between the case and control
#' condition if most of the "case" conditions clustered strongly together. Batch effects can also be
#' identified in this way.  For instance, if samples were prepared at two different locations and this
#' information was incorporated into the metadata as a "batch" column, visualizing the data based on
#' batch data could indicate whether or not differences in the data are due to differences in the sample
#' prep location, rather than due to a biological effect, depending on how strongly the batches cluster.
#' @param filename Filename for output plot.  Valid extensions are ".pdf" and ".png".
#' File generation can be turned off using set_output_mode("screen").
#' Output will be written to /dendrograms/ directory.
#' @param id_field The unique sample ID field for your data.  Should correspond to a column in targets data.
#' @param groupBy A field in target meta data that you wish to cluster data by.
#' @return This function does not return a value.
#' @keywords dendrogram batch cluster
#' @export
#' @examples
#' \dontrun{
#'
#' #Plot a dendrogram based on an "infection" metadata field.
#' plot_dendro(filename="infection_dendro.pdf", "sample_id", "infection")
#'
#' #Plot a dendrogram based on an "batch" metadata field.
#' plot_dendro(filename="batch_dendro.pdf", "sample_id", "batch")
#'
#' }
plot_dendro <- function(filename="dendro_plot.pdf", id_field, groupBy)
{
  imgHeight = 10
  imgWidth = 10

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
  if(typeof(id_field) != "character")
  {
    stop("group variable must be a string.")
    return(-1)
  }
  if(!exists("dendro_dir", envir=.DEVis_env))
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
  if(!exists("countDat", envir=.DEVis_env))
  {
    stop("plot_dendro(): Count data not present.  Please run prep_counts() first.")
    return(-1)
  }
  if(!exists("tgt_dat", envir=.DEVis_env))
  {
    stop("plot_dendro(): Target data not present.  Please run prep_targets() first.")
    return(-1)
  }
  if(!exists("output_mode", envir=.DEVis_env))
  {
    stop("plot_dendro(): Output mode not initialized.  Please run set_output_mode() first.")
    return(-1)
  }
  if(!exists("stabilized_data", envir=.DEVis_env))
  {
    stop("plot_dendro() requires that prep_dds_from_data() be run first.")
    return(-1)
  }
  if(!(groupBy %in% colnames(colData(.DEVis_env$stabilized_data))))
  {
    print("plot_dendro() requires grouping to correspond to annotation.")
    print("Possible values are: ")
    print(colnames(colData(.DEVis_env$stabilized_data)))
    stop("plot_dendro(): invalid group variable selection.")
    return(-1)
  }
  if(!(id_field %in% colnames(colData(.DEVis_env$stabilized_data))))
  {
    print("plot_dendro() requires id_field to correspond to annotation.")
    print("Possible values are: ")
    print(colnames(colData(.DEVis_env$stabilized_data)))
    stop("plot_dendro(): invalid group variable selection.")
    return(-1)
  }

  dat.dist <- dist(t(.DEVis_env$countDat))
  hc <- hclust(dat.dist)
  dendr <- dendro_data(hc, type="rectangle")

  tgt_dat <- .DEVis_env$tgt_dat
  sub1 <- paste("tgt_dat$", groupBy, "[match(dendr$labels$label,tgt_dat$", id_field, ")]", sep="")
  dendr$labels$label <- eval(parse(text = sub1))

  #Automatically scale the image so all the genes fit and labels are readable.
  numSamples <- length(dendr$labels$label)
  if(numSamples > 50)
  {
    imgHeight = (numSamples / 10) + 5
  }

  #Satisfy check() that our ggplot aes variables are indeed not globals.
  x = NULL
  y = NULL
  xend = NULL
  yend = NULL

  #Complete the dendrogram.
  plot_out <- ggplot() +
    geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) +
    geom_text(data=label(dendr), aes(x=x, y=y, label=label, hjust=0), size=3) +
    coord_flip() + scale_y_reverse(expand=c(0.2, 0)) +
    theme_stata() + scale_colour_stata() +
    theme(axis.line.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          axis.line.x=element_blank(),
          panel.background=element_rect(fill="white"),
          panel.grid=element_blank()) +
    theme(text = element_text(size=22,margin = margin(t = 0, r = 10, b = 0, l = 0))) +
    labs(x="", y="")

  #Draw the dendrogram.
  setwd(.DEVis_env$dendro_dir)
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
}
