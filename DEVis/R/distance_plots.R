#' Visualize the Euclidian distances between samples.
#'
#' This function computes and plots sample Euclidian distances in a symmetrical heatmap. This visualization
#' can be used to examine overall similarity between samples and whether or not data are
#' behaving as expected based on the experimental design.  This view can also be used to
#' identify potential outlying samples.
#' @param row_labels Row labels for samples. This value should correspond to a column in
#' the target data and will be used to label each row in the heatmap. Viewing the plot in the
#' context of different metadata makes it possible to search for outliers and batch effects, in
#' addition to whether or not the data behave as intended. String.
#' @param filename Output filename in string format. Valid formats are .pdf and .png.
#' File generation can be turned off using set_output_mode("screen").
#' Output will be written to the /sample_distance/euclidian/ directory.
#' @param theme Determines color scheme used for this plot.  Valid options are integers, 1-6.
#' @param returnData Boolean.  Determines if this visualization should return data used to generate the visualization. Default=FALSE.
#' @return If returnData is true, this function will return the sample distance matrix of
#' euclidian distance measurements between samples.
#' @keywords distance cluster Euclidian batch outlier
#' @seealso \code{\link{plot_poisson_dist}}
#' @export
#' @examples
#' \dontrun{
#'
#' #Visualize Euclidian distances of all samples using "SampleID" target data as labels.
#' plot_euclid_dist("SampleID", filename="euclidian_distance.pdf",
#'                    theme=2, returnData=FALSE)
#'
#' #Visualize Euclidian distances of all samples using "timepoint" target data as labels.
#' #Store the resulting distance matrix data.
#' distMatrix <- plot_euclid_dist("timepoint", filename="euclidian_distance.pdf",
#'                                  theme=2, returnData=TRUE)
#'
#' }
plot_euclid_dist <- function(row_labels, filename="euclidian_distance.pdf", theme=1, returnData=FALSE)
{
  #Validate parameters and data.
  if(!exists("stabilized_data", envir=.DEVis_env))
  {
    stop("euclid_dist() requires that prep_dds_from_data() be run first.")
    return(-1)
  }
  if(!exists("output_mode", envir=.DEVis_env))
  {
    stop("euclid_dist(): Output mode not initialized.  Please run set_output_mode() first.")
    return(-1)
  }
  if(!exists("euclidian_dir", envir=.DEVis_env))
  {
    stop("Directories not initialized.  Try running create_dir_struct().")
    return(-1)
  }
  if(typeof(returnData) != "logical")
  {
    stop("returnData must be a logical value. Please enter TRUE or FALSE.")
    return(-1)
  }
  if(!(row_labels %in% colnames(colData(.DEVis_env$stabilized_data))))
  {
    print("euclid_dist() requires that labels correspond to colData annotation.")
    print("Possible values are: ")
    print(colnames(colData(.DEVis_env$stabilized_data)))
    stop("euclid_dist(): invalid variable selection.")
    return(-1)
  }
  else
  {
    stabilized_data <- .DEVis_env$stabilized_data
    sub1 <- paste("colData(stabilized_data)$", row_labels, sep="")
    row_labels <- eval(parse(text = sub1))
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

  #Auto-scale image size.
  numSamples = length(row_labels)
  imgHeight = 10
  imgWidth = 10
  if(numSamples > 50)
  {
    imgHeight = (numSamples / 10) + 5
    imgWidth = (numSamples / 10) + 5
  }

  #Calculate distance measurements.
  sampleDists <- dist( t( assay(.DEVis_env$stabilized_data) ) )
  setwd(.DEVis_env$euclidian_dir)
  sampleDistMatrix <- as.matrix( sampleDists )
  rownames(sampleDistMatrix) <- row_labels
  colnames(sampleDistMatrix) <- NULL

  #Determine color scheme.
  if(theme == 1)
  {
    colSet <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  }
  if(theme == 2)
  {
    colSet <- colorRampPalette( rev(brewer.pal(9, "YlOrRd")) )(255)
  }
  if(theme == 3)
  {
    colSet <- colorRampPalette( rev(brewer.pal(9, "YlGn")) )(255)
  }
  if(theme == 4)
  {
    colSet <- colorRampPalette( rev(brewer.pal(9, "GnBu")) )(255)
  }
  if(theme == 5)
  {
    colSet <- colorRampPalette( rev(brewer.pal(9, "PuRd")) )(255)
  }
  if(theme == 6)
  {
    colSet <- colorRampPalette( rev(brewer.pal(9, "Greys")) )(255)
  }

  if(.DEVis_env$output_mode == "both" | .DEVis_env$output_mode == "file")
  {
    #Create the output file.
    if(filetype == ".pdf")
    {
      pdf(filename, width=imgWidth, height=imgHeight, onefile=FALSE)
    }
    else
    {
      png(filename, width=imgWidth, height=imgHeight, units='in', res=600)
    }

    #Draw the plot and finish up.
    pheatmap(sampleDistMatrix,
             clustering_distance_rows=sampleDists,
             clustering_distance_cols=sampleDists,
             color=colSet,
             fontsize_row=13)

    dev.off()
  }
  if(.DEVis_env$output_mode == "screen" | .DEVis_env$output_mode == "both")
  {
    pheatmap(sampleDistMatrix,
             clustering_distance_rows=sampleDists,
             clustering_distance_cols=sampleDists,
             color=colSet,
             fontsize_row=13)
  }

  setwd(.DEVis_env$working_dir)

  if(returnData)
  {
    return(sampleDistMatrix)
  }
}


#' Visualize the poisson distances between samples.
#'
#' This function computes and plots sample poisson distances in a symmetrical heatmap. This visualization
#' can be used to examine overall similarity between samples and whether or not data are
#' behaving as expected based on the experimental design.  Unlike Euclidian distance, poisson distance
#' calculation takes into account variation between sample counts when calculating distances.
#' This view can also be used to identify potential outlying samples.
#' @param row_labels Row labels for samples. This value should correspond to a column in
#' the target data and will be used to label each row in the heatmap. Viewing the plot in the
#' context of different metadata makes it possible to search for outliers and batch effects, in
#' addition to whether or not the data behave as intended. String.
#' @param filename Output filename in string format. Valid formats are .pdf and .png.
#' File generation can be turned off using set_output_mode("screen").
#' Output will be written to the /sample_distance/poisson/ directory.
#' @param theme Determines color scheme used for this plot.  Valid options are integers, 1-6.
#' @param returnData Boolean. Determines if this visualization should return data used to generate the visualization. Default=FALSE.
#' @return If returnData is true, this function will return the sample distance matrix of
#' poisson distance measurements between samples.
#' @keywords distance cluster poisson batch outlier
#' @seealso \code{\link{plot_euclid_dist}}
#' @export
#' @examples
#' \dontrun{
#'
#' #Visualize Euclidian distances of all samples using "SampleID" target data as labels.
#' plot_poisson_dist("SampleID", filename="poisson_distance.pdf", theme=2, returnData=FALSE)
#'
#' #Visualize Euclidian distances of all samples using "timepoint" target data as labels.
#' #Store the resulting distance matrix data.
#' distMatrix <- plot_poisson_dist("timepoint", filename="poisson_distance.pdf",
#'                                  theme=2, returnData=TRUE)
#'
#' }
plot_poisson_dist <- function(row_labels, filename="poisson_distance.pdf", theme=1, returnData=FALSE)
{
  #Validate parameters and data.
  if(!exists("stabilized_data", envir=.DEVis_env))
  {
    stop("poisson_dist() requires that prep_dds_from_data() be run first.")
    return(-1)
  }
  if(!exists("output_mode", envir=.DEVis_env))
  {
    stop("poisson_dist(): Output mode not initialized.  Please run set_output_mode() first.")
    return(-1)
  }
  if(typeof(returnData) != "logical")
  {
    stop("returnData must be a logical value. Please enter TRUE or FALSE.")
    return(-1)
  }
  if(!exists("poisson_dir", envir=.DEVis_env))
  {
    stop("Directories not initialized.  Try running create_dir_struct().")
    return(-1)
  }
  if(!(row_labels %in% colnames(colData(.DEVis_env$stabilized_data))))
  {
    print("poisson_dist() requires that labels correspond to colData annotation.")
    print("Possible values are: ")
    print(colnames(colData(.DEVis_env$stabilized_data)))
    stop("poisson_dist(): invalid variable selection.")
    return(-1)
  }
  else
  {
    stabilized_data <- .DEVis_env$stabilized_data
    sub1 <- paste("colData(stabilized_data)$", row_labels, sep="")
    row_labels <- eval(parse(text = sub1))
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

  #Auto-scale image size.
  numSamples = length(row_labels)
  imgHeight = 10
  imgWidth = 10
  if(numSamples > 50)
  {
    imgHeight = (numSamples / 10) + 5
    imgWidth = (numSamples / 10) + 5
  }

  #Determine color scheme.
  if(theme == 1)
  {
    colSet <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  }
  if(theme == 2)
  {
    colSet <- colorRampPalette( rev(brewer.pal(9, "YlOrRd")) )(255)
  }
  if(theme == 3)
  {
    colSet <- colorRampPalette( rev(brewer.pal(9, "YlGn")) )(255)
  }
  if(theme == 4)
  {
    colSet <- colorRampPalette( rev(brewer.pal(9, "GnBu")) )(255)
  }
  if(theme == 5)
  {
    colSet <- colorRampPalette( rev(brewer.pal(9, "PuRd")) )(255)
  }
  if(theme == 6)
  {
    colSet <- colorRampPalette( rev(brewer.pal(9, "Greys")) )(255)
  }


  #Calculate distance measurements.
  poisd <- PoissonDistance(t(assay(.DEVis_env$stabilized_data)))
  samplePoisDistMatrix <- as.matrix( poisd$dd )
  rownames(samplePoisDistMatrix) <- row_labels
  colnames(samplePoisDistMatrix) <- NULL

  #Create the output file.
  setwd(.DEVis_env$poisson_dir)
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

    #Draw the plot and finish up.
    suppressWarnings({
      pheatmap(samplePoisDistMatrix,
               clustering_distance_rows=poisd$dd,
               clustering_distance_cols=poisd$dd,
               color=colSet,
               fontsize_row=15)
    })
    dev.off()
  }
  if(.DEVis_env$output_mode == "screen" | .DEVis_env$output_mode == "both")
  {
    suppressWarnings({
      pheatmap(samplePoisDistMatrix,
               clustering_distance_rows=poisd$dd,
               clustering_distance_cols=poisd$dd,
               color=colSet,
               fontsize_row=15)
    })
  }
  setwd(.DEVis_env$working_dir)
}
