% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distance_plots.R
\name{plot_euclid_dist}
\alias{plot_euclid_dist}
\title{Visualize the Euclidian distances between samples.}
\usage{
plot_euclid_dist(row_labels, filename = "euclidian_distance.pdf",
  theme = 1, returnData = FALSE)
}
\arguments{
\item{row_labels}{Row labels for samples. This value should correspond to a column in
the target data and will be used to label each row in the heatmap. Viewing the plot in the
context of different metadata makes it possible to search for outliers and batch effects, in
addition to whether or not the data behave as intended. String.}

\item{filename}{Output filename in string format. Valid formats are .pdf and .png.
File generation can be turned off using set_output_mode("screen").
Output will be written to the /sample_distance/euclidian/ directory.}

\item{theme}{Determines color scheme used for this plot.  Valid options are integers, 1-6.}

\item{returnData}{Boolean.  Determines if this visualization should return data used to generate the visualization. Default=FALSE.}
}
\value{
If returnData is true, this function will return the sample distance matrix of
euclidian distance measurements between samples.
}
\description{
This function computes and plots sample Euclidian distances in a symmetrical heatmap. This visualization
can be used to examine overall similarity between samples and whether or not data are
behaving as expected based on the experimental design.  This view can also be used to
identify potential outlying samples.
}
\examples{
\dontrun{

#Visualize Euclidian distances of all samples using "SampleID" target data as labels.
plot_euclid_dist("SampleID", filename="euclidian_distance.pdf",
                   theme=2, returnData=FALSE)

#Visualize Euclidian distances of all samples using "timepoint" target data as labels.
#Store the resulting distance matrix data.
distMatrix <- plot_euclid_dist("timepoint", filename="euclidian_distance.pdf",
                                 theme=2, returnData=TRUE)

}
}
\seealso{
\code{\link{plot_poisson_dist}}
}
\keyword{Euclidian}
\keyword{batch}
\keyword{cluster}
\keyword{distance}
\keyword{outlier}
