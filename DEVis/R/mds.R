
#' Visualize multi-dimensional scaled data for all samples, with group-wise metadata incorporated as colors and shapes.
#'
#' This function performs multi-dimensional scaling (MDS) on data to visualize levels of similarity between individual samples.
#' Colors and shape variables can be specified according to available metadata, making it possible to identify patterns in the data
#' based on groups.  For instance, by coloring data points based on "Infected/Non-infected" groups, and drawing data point shapes based on time
#' point, it is possible to view the differences between infected and non-infected samples across multiple time points.
#' @param filename Filename for output plot.  Valid extensions are ".pdf" and ".png".
#' File generation for plot file can be turned off using set_output_mode("screen").
#' Output will be written to the /MDS/standard/ directory.
#' @param color_var The group from target data that should be indicated by color.
#' Any column of metadata can be used, regardless of experimental design.
#' @param shape_var The group from target data that should be indicated by shape.
#' Any column of metadata can be used, regardless of experimental design. (Optional). Default="none".
#' @param showConf Boolean.  Draw an ellipsis representing the 95 percent confidence interval around each group. Default: TRUE
#' @param theme Theme for the layout and color scheme for the MDS plot.  Valid selections are integers between 1-6.
#' @param customLabels If customLabels is set to TRUE, the user will be prompted to provide a custom label for each label.
#' @param returnData If this value is true, this function will return a sample distance matrix. Default: FALSE
#' @return If returnData is true, this function will return a matrix containing sample distances computed to create the MDS plot.
#' @seealso \code{\link{plot_mds_hulls}}
#' @keywords mds distance visualization
#' @export
#' @examples
#' \dontrun{
#'
#' #These examples assume "Time" and "Infection" columns exist in target metadata.
#'
#' #Make a MDS plot showing all samples, colored based on time point.
#' plot_mds(filename="mds_plot.pdf", color_var="Time", shape_var="none",
#'           theme=1, returnData=FALSE)
#'
#' #Make a MDS plot, colored based on time point with shapes based on infection.
#' plot_mds(filename="mds_plot.pdf", color_var="Time", shape_var="Infection",
#'           theme=2, returnData=FALSE)
#'
#' }
plot_mds <- function(filename="mds_plot.pdf", color_var, shape_var="none", showConf=TRUE, theme=1, customLabels=FALSE, returnData=FALSE)
{
  imgWidth  = 10
  imgHeight = 10
  do_shape  = FALSE

  #Validate the data we want exists and was processed correctly.
  if(!exists("stabilized_data", envir=.DEVis_env))
  {
    stop("plot_mds() requires that prep_dds_from_data() be run first.")
    return(-1)
  }
  if(!exists("output_mode", envir=.DEVis_env))
  {
    stop("plot_mds(): Output mode not initialized.  Please run set_output_mode() first.")
    return(-1)
  }
  if(typeof(filename) != "character")
  {
    stop("filename must be a string.")
    return(-1)
  }
  if(!exists("MDS_dir", envir=.DEVis_env))
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
  if(typeof(showConf) != "logical")
  {
    stop("showConf must be a logical value. Please enter TRUE or FALSE.")
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

  #Determine how to group color and shape variables.
  if(!(color_var %in% colnames(colData(.DEVis_env$stabilized_data))))
  {
    print("plot_mds() requires that color grouping to correspond to colData annotation.")
    print("Possible values are: ")
    print(colnames(colData(.DEVis_env$stabilized_data)))
    stop("plot_mds(): invalid color variable selection.")
    return(-1)
  }
  else
  {
    stabilized_data <- .DEVis_env$stabilized_data
    sub1 <- paste("colData(stabilized_data)$", color_var, sep="")
    color_var <- eval(parse(text = sub1))
  }
  if(shape_var != "none")
  {
    do_shape = TRUE
    if(!(shape_var %in% colnames(colData(.DEVis_env$stabilized_data))))
    {
      print("plot_mds() requires that shape grouping to correspond to colData annotation, or be 'none'.")
      print("Possible values are: ")
      print(colnames(colData(.DEVis_env$stabilized_data)))
      stop("plot_mds(): invalid shape variable selection.")
      return(-1)
    }
    else
    {
      stabilized_data <- .DEVis_env$stabilized_data
      sub1 <- paste("colData(stabilized_data)$", shape_var, sep="")
      shape_var <- eval(parse(text = sub1))
    }
  }

  #Calculate distances matrix.
  sampleDists <- dist( t( assay(.DEVis_env$stabilized_data) ) )
  sampleDistMatrix <- as.matrix( sampleDists )
  mds <- data.frame(cmdscale(sampleDistMatrix))
  mds <- cbind(mds, as.data.frame(colData(.DEVis_env$stabilized_data)))

  #Satisfy check() that our ggplot aes variables are indeed not globals.
  X1    = NULL
  X2    = NULL

  #If custom labels are requested retrieve and update them here.
  if(customLabels == TRUE)
  {
    curLabs = as.character(levels(color_var))
    relabel <- as.character(color_var)

    #Get updated labels for each level in the current label set.
    for(i in 1:length(curLabs))
    {
      curLab <- as.character(readline(prompt=paste("Enter a custom label corresponding to ", curLabs[i], ": ", sep="")))
      sub1 <- paste("relabel[grepl('", curLabs[i], "',relabel)] <- '", curLab, "'", sep="")
      eval(parse(text = sub1))
    }
    color_var <- relabel

    if(do_shape)
    {
      curLabs = as.character(levels(shape_var))
      relabel <- as.character(shape_var)

      #Get updated labels for each level in the current label set.
      for(i in 1:length(curLabs))
      {
        curLab <- as.character(readline(prompt=paste("Enter a custom label corresponding to ", curLabs[i], ": ", sep="")))
        sub1 <- paste("relabel[grepl('", curLabs[i], "',relabel)] <- '", curLab, "'", sep="")
        eval(parse(text = sub1))
      }
      shape_var <- relabel
    }
  }

  #Start creating the plot.
  if(do_shape)
  {
    plot_out <- ggplot(mds, aes(X1,X2,color=color_var, shape=shape_var)) + geom_point(size=3)
  }
  else
  {
    plot_out <- ggplot(mds, aes(X1,X2,color=color_var)) + geom_point(size=3)
  }

  #Draw confidence interval ellipses.
  if(showConf)
  {
    plot_out <- plot_out + stat_ellipse(aes(x=X1, y=X2, group=color_var), type = "norm")
  }

  #Apply the appropriate theme.
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

  #Finish the details.
  plot_out <- plot_out +
    labs(x="Dimension 1", y="Dimension 2") +
    theme(legend.text = element_text(size=12)) +
    theme(legend.title=element_blank()) +
    theme(axis.text=element_text(size=10, face="bold")) +
    theme(text = element_text(size=12,margin = margin(t = 0, r = 2, b = 1, l = 2)))

  #Create output file.
  setwd(.DEVis_env$MDS_dir)
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
    return(mds)
  }
}

#' Visualize multi-dimensional scaled data for all samples, with group-wise metadata incorporated as colors, shapes, and convex hulls.
#'
#' This function performs multi-dimensional scaling (MDS) on data to visualize levels of similarity between individual samples.
#' Distance measurements can be calculated based only on the differentially expressed genes for each sample.
#' Colors and shape variables can be specified according to available metadata, making it possible to identify patterns in the data
#' based on groups.  Additionally, convex hulls are computed and drawn indicating the overall clustering range of samples in a group-wise manner.
#'
#' @param filename Filename for output plot.  Valid extensions are ".pdf" and ".png".
#' File generation for plot file can be turned off using set_output_mode("screen").
#' Output will be written to the /MDS/hulls/ directory.
#' @param color_var The group from target data that should be indicated by color.
#' Any column of metadata can be used, regardless of experimental design. String.
#' @param shape_var The group from target data that should be indicated by shape.
#' Any column of metadata can be used, regardless of experimental design. String. Optional. Default="none".
#' @param deOnly Use only DE genes when computing the plot.  Requires that create_master_res() be run first. Boolean.
#' @param showLabel Show labels indicating sample identifiers on the plot. Boolean.  Default=FALSE.
#' @param hullType Determines hull type, either "solid" or "outline". Solid hulls are partially transparent overlaying polygons.
#' Outline hulls trace an outline along the outermost distant samples of each group.  Default="solid".
#' @param exclude_data Exclude some subset of data from this plot. If this is TRUE, three additional parameters must be provided:
#' idCol, excludeCol, and excludeName.  These will then be used to subset the data for the MDS as desired.  For example,
#' it is possible to remove "Day1" data from your "Timepoint" field, or "Control_Sample" from your condition field.
#' I.E. Remove control samples from the visualization. Logical. Default=FALSE.
#' @param idCol Required only if exclude_data is TRUE.  The ID column in target data that corresponds to column names in count data.
#' @param excludeCol Required only if exclude_data is TRUE. The field from which exclusion criteria will be determined.
#' Must match a column name from target data.
#' @param excludeName Required only if exclude_data is TRUE. Value to exclude from the data based on the excludeCol.
#' Must match a level in the excludeCol.
#' @param theme Theme for the layout and color scheme for the MDS plot.  Valid selections are integers between 1-6.
#' @param customLabels If customLabels is set to TRUE, the user will be prompted to provide a custom label for each label.
#' @param returnData If this value is true, this function will return a sample distance matrix. Default: FALSE
#' @return If returnData is true, this function will return a matrix containing sample distances computed to create the MDS plot.
#' @seealso \code{\link{plot_mds}}
#' @keywords mds hull distance hull visualization DE
#' @export
#' @examples
#' \dontrun{
#'
#' #These examples assume "Time", "Infection", and sampleType columns exist in target metadata.
#'
#'  /*
#'   * Create a MDS plot with convex hulls drawn based on "Time" metadata.
#'   * Shapes are based on "Infection".  Compute based on data from all genes.
#'   */
#' plot_mds_hulls("mds_hull_plot.pdf", "Time", "Infection",
#'                 deOnly=FALSE, showLabel=FALSE, hullType="solid",
#'                 theme=1)
#'
#'  /*
#'   * Create a MDS plot with convex hulls drawn based on "Time" metadata.
#'   * Shapes are based on "Infection".  Compute based on data from only
#'   * differentailly expressed genes.  Hulls drawn only as outlines.
#'   */
#' plot_mds_hulls("mds_hull_deOnly_plot.pdf", "Time", "Infection",
#'                 deOnly=TRUE, showLabel=FALSE, hullType="outline",
#'                 theme=2)
#'
#'  /*
#'   * Create a MDS plot with convex hulls drawn based on "Time" metadata.
#'   * Shapes are based on "Infection".  Compute based on data from only
#'   * differentailly expressed genes.  Do not show samples labeled
#'   * as "mock" from "sampleType" column of target metadata.
#'   */
#' plot_mds_hulls("mds_hull_deOnly_plot.pdf", "Time", "Infection",
#'                 deOnly=TRUE, showLabel=FALSE, hullType="outline",
#'                 exclude_data=TRUE, idCol="sampleID",
#'                 excludeCol="sampleType", excludeName="mock"
#'                 theme=4)
#'
#' }
plot_mds_hulls <- function(filename="mds_hulls_plot.pdf", color_var, shape_var="none", deOnly=FALSE, showLabel=FALSE, hullType="solid", exclude_data=FALSE, idCol="", excludeCol="", excludeName="", theme=1, customLabels=FALSE, returnData=FALSE)
{
  imgHeight = 10
  imgWidth  = 10

  #Validate the data we want exists and was processed correctly.
  if(!exists("stabilized_data", envir=.DEVis_env))
  {
    stop("plot_mds_hulls() requires that prep_dds_from_data() be run first.")
    return(-1)
  }
  if(!exists("output_mode", envir=.DEVis_env))
  {
    stop("plot_mds_hulls(): Output mode not initialized.  Please run set_output_mode() first.")
    return(-1)
  }
  if(typeof(filename) != "character")
  {
    stop("filename must be a string.")
    return(-1)
  }
  if(typeof(excludeCol) != "character" | typeof(excludeName) != "character" | typeof(idCol) != "character")
  {
    print("Column exclusion requires the names of the id column that corresponds to count data column names, the name of the column to select inclusion/exclusion, and the name of the items to exclude from the dataset in the provided column.")
    print('I.E. idCol="Sample_ID", excludeCol="Condition", excludeName="Mock",  would remove all samples labeled "Mock" in the "Condition" column of the target data.')
    print("Possible column values are: ")
    print(colnames(colData(.DEVis_env$stabilized_data)))
    stop("Please provide valid column names as string values.")
    return(-1)
  }
  if(exclude_data)
  {
    if(!(idCol %in% colnames(colData(.DEVis_env$stabilized_data))))
    {
      print("Invalid ID column selected.")
      print("Possible column values are: ")
      print(colnames(colData(.DEVis_env$stabilized_data)))
      stop("Please provide an ID column that corresponds to to count data column names.")
      return(-1)
    }
    if(!(excludeCol %in% colnames(colData(.DEVis_env$stabilized_data))))
    {
      print("Invalid exclusion column selected.")
      print("Possible column values are: ")
      print(colnames(colData(.DEVis_env$stabilized_data)))
      stop("Please provide a valid column name.")
      return(-1)
    }
    stabilized_data <- .DEVis_env$stabilized_data
    sub1 <- paste("colData(stabilized_data)$", excludeCol, sep="")
    exclude_col_data <- eval(parse(text = sub1))
    for(i in 1:length(excludeName))
    {
      if(nchar(excludeName[i]) == 0)
      {
        print("Exclusion name not provided.")
        stop("Please provide a valid column name.")
        return(-1)
      }
      if(!(excludeName[i] %in% exclude_col_data))
      {
        print(excludeName[i])
        print("Exclusion name does not exist in provided exclusion column.")
        print("Possible values for exclusions are: ")
        print(levels(exclude_col_data))
        stop("Please provide a valid selection.")
        return(-1)
      }
    }
  }
  if(!exists("MDS_hull_dir", envir=.DEVis_env))
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
  hullType_options = c("solid", "outline")
  if(!(hullType %in% hullType_options))
  {
    stop('hullType must be "solid" or "outline"')
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

  #Give the user the option to plot only data from DE genes.
  #Requires that create_master_res() already ran so the union_names list can exist between samples.
  if(deOnly == TRUE)
  {
    if(!exists("aggregate_names", envir=.DEVis_env))
    {
      stop("plot_mds_hulls() requires that create_master_res() be run first if only DE genes will be used.")
      return(-1)
    }
    else
    {
      data <- data.frame(assay(.DEVis_env$stabilized_data)[.DEVis_env$aggregate_names,])
    }
  }
  else
  {
    data <- data.frame(assay(.DEVis_env$stabilized_data))
  }

  #If some data is to be excluded from the plot, calculate and filter the data portion and save the exclusion list.
  includeList <- c()
  excludeList <- c()
  if(exclude_data)
  {
    stabilized_data <- .DEVis_env$stabilized_data
    sub1 <- paste("colData(stabilized_data)$", excludeCol, sep="")
    exclude_col_data <- eval(parse(text = sub1))
    sub2 <- paste("colData(stabilized_data)$", idCol, sep="")
    id_col_data <- eval(parse(text = sub2))

    includeList <- id_col_data[which(exclude_col_data != excludeName[i])]
    data <- data[,includeList]
  }

  #Determine how to group color and shape variables.
  if(!(color_var %in% colnames(colData(.DEVis_env$stabilized_data))))
  {
    print("plot_mds_hulls() requires that color grouping to correspond to colData annotation.")
    print("Possible values are: ")
    print(colnames(colData(.DEVis_env$stabilized_data)))
    stop("plot_mds_hulls(): invalid color variable selection.")
    return(-1)
  }
  else
  {
    if(exclude_data)
    {
      stabilized_data <- .DEVis_env$stabilized_data
      sub1 <- paste("colData(stabilized_data)$", color_var, "[colData(stabilized_data)$", idCol, " %in% includeList]",  sep="")
      color_var <- eval(parse(text = sub1))
    }
    else
    {
      stabilized_data <- .DEVis_env$stabilized_data
      sub1 <- paste("colData(stabilized_data)$", color_var, sep="")
      color_var <- eval(parse(text = sub1))
    }
  }
  do_shape = FALSE
  if(shape_var != "none")
  {
    do_shape = TRUE
    if(!(shape_var %in% colnames(colData(.DEVis_env$stabilized_data))))
    {
      print("plot_mds_hulls() requires that shape grouping to correspond to colData annotation.")
      print("Possible values are: ")
      print(colnames(colData(.DEVis_env$stabilized_data)))
      stop("plot_mds_hulls(): invalid shape variable selection.")
      return(-1)
    }
    else
    {
      if(exclude_data)
      {
        stabilized_data <- .DEVis_env$stabilized_data
        sub1 <- paste("colData(stabilized_data)$", shape_var, "[colData(stabilized_data)$", idCol, " %in% includeList]",  sep="")
        shape_var <- eval(parse(text = sub1))
      }
      else
      {
        stabilized_data <- .DEVis_env$stabilized_data
        sub1 <- paste("colData(stabilized_data)$", shape_var, sep="")
        shape_var <- eval(parse(text = sub1))
      }
    }
  }

  #If custom labels are requested retrieve and update them here.
  if(customLabels == TRUE)
  {
    curLabs = as.character(levels(color_var))
    relabel <- as.character(color_var)

    #Get updated labels for each level in the current label set.
    for(i in 1:length(curLabs))
    {
      curLab <- as.character(readline(prompt=paste("Enter a custom label corresponding to ", curLabs[i], ": ", sep="")))
      sub1 <- paste("relabel[grepl('", curLabs[i], "',relabel)] <- '", curLab, "'", sep="")
      eval(parse(text = sub1))
    }
    color_var <- relabel

    if(do_shape)
    {
      curLabs = as.character(levels(shape_var))
      relabel <- as.character(shape_var)

      #Get updated labels for each level in the current label set.
      for(i in 1:length(curLabs))
      {
        curLab <- as.character(readline(prompt=paste("Enter a custom label corresponding to ", curLabs[i], ": ", sep="")))
        sub1 <- paste("relabel[grepl('", curLabs[i], "',relabel)] <- '", curLab, "'", sep="")
        eval(parse(text = sub1))
      }
      shape_var <- relabel
    }
  }


  suppressMessages({
    #Transpose data, make 0 values essentially 0.
    d <- dist(t(data))
    d[d==0] = 0.0000000000000000000000000000001

    #Fit data
    fit <- isoMDS(d)
    x <- fit$points[,1]
    y <- fit$points[,2]

    min_lim = min(min(x),min(y))
    max_lim = max(max(x),max(y))

    datai = data.frame(x=x, y=y)

    datai$id = names(x)
    datai$cols = as.factor(color_var)
    if(do_shape)
    {
      datai$shape = as.factor(shape_var)
      datai$hull = as.factor(paste0(color_var))
    }
    else
    {
      datai$hull = as.factor(paste0(color_var))
    }

    #Compute the subset of points that lie on the convex hull.
    find_hull <- function(df) df[chull(df$x, df$y), ]
    hulls <- ddply(datai, "hull", find_hull)
  })

  #Satisfy check() that our ggplot aes variables are indeed not globals.
  cols    = NULL
  shape   = NULL
  hull    = NULL

  #Plot the results.
  if(do_shape)
  {
    png1 <- ggplot(datai, aes(x= x, y=y,colour=cols,shape=shape,label = id))
  }
  else
  {
    png1 <- ggplot(datai, aes(x= x, y=y,colour=cols,label = id))
  }
  #If they want labels add them.
  if(showLabel)
  {
    png1 <- png1 + geom_text(aes(label=id),hjust=0, vjust=0, size=2)
  }

  #Create the plot.
  png <- png1 +
    geom_point(alpha=1,size = 3) +
    scale_shape_manual(values=(0:(length(levels(datai$shape))-1))) +
    xlim(min_lim, max_lim)+ylim (min_lim, max_lim)

  #Set solid or outline hull type.
  if(hullType=="solid")
  {
    png <- png + geom_polygon(data = hulls,aes(group = hull, fill=cols), alpha = 0.4, size = 0.5)
  }
  else
  {
    png <- png + geom_polygon(data = hulls,aes(group = hull, fill=cols), alpha = 0.4, size = 0.5, fill=NA)
  }


  #Apply the selected theme.
  if(theme == 1)
  {
    stata_long_pal = c(stata_pal("s2color")(15), stata_pal("s1rcolor")(15))
    png <- png + theme_stata() + scale_color_manual(values=stata_long_pal) + scale_fill_manual(values=stata_long_pal)
  }
  if(theme == 2)
  {
    nature_pal = c(pal_npg("nrc")(10),pal_npg("nrc", alpha = 0.6)(10),pal_npg("nrc", alpha = 0.3)(15) )
    png <- png + theme_igray() + scale_color_manual(values=nature_pal) + scale_fill_manual(values=nature_pal)
  }
  if(theme == 3)
  {
    tron_pal = c(pal_tron("legacy")(7),pal_tron("legacy", alpha = 0.7)(7),pal_tron("legacy", alpha = 0.5)(7),pal_tron("legacy", alpha = 0.3)(7),pal_tron("legacy", alpha = 0.2)(2))
    png <- png + theme_hc(bgcolor = "darkunica") + theme(axis.text.x = element_text(colour = "white"), axis.text.y = element_text(colour = "white")) + scale_color_manual(values=tron_pal) + scale_fill_manual(values=tron_pal)
  }
  if(theme == 4)
  {
    gdoc_pal <- c(pal_ucscgb("default")(26), pal_ucscgb("default",alpha=.5)(4))
    png <- png + theme_gdocs() + scale_color_manual(values=gdoc_pal) + scale_fill_manual(values=gdoc_pal)
  }
  if(theme == 5)
  {
    d3_pal <- c(pal_d3("category20")(20), pal_d3("category10",alpha=.5)(10))
    png <- png + theme_solarized() + scale_color_manual(values=d3_pal) + scale_fill_manual(values=d3_pal)
  }
  if(theme == 6)
  {
    png <- png + theme_bw() + scale_colour_grey(start = 0, end = .9) + scale_fill_grey(start = 0, end = .9)
  }

  png <- png +
    theme(legend.text = element_text(size=14)) +
    theme(axis.text.x = element_text(angle=90, hjust=1)) +
    theme(axis.text=element_text(size=12, face="bold")) +
    theme(legend.title=element_blank()) +
    theme(text = element_text(size=12, margin = margin(t = 0, r = 2, b = 1, l = 2))) +
    labs(x="Dimension 1", y="Dimension 2") +
    guides(colour = guide_legend(override.aes = list(size=1, alpha = 1)))


  #Save the plot.
  setwd(.DEVis_env$MDS_hull_dir)
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
    print(png)
    dev.off()
  }
  if(.DEVis_env$output_mode == "screen" | .DEVis_env$output_mode == "both")
  {
    print(png)
  }
  setwd(.DEVis_env$working_dir)

  #Return the data if requested.
  if(returnData)
  {
    return(datai)
  }
}
