% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/boxplot.R
\name{plot_gene}
\alias{plot_gene}
\title{Visualize the expression of a specific gene with regard to metadata grouping.}
\usage{
plot_gene(filename = "gene_plot.pdf", gene_name, groupBy, theme = 1,
  returnData = FALSE)
}
\arguments{
\item{filename}{Filename for output plot.  Valid extensions are ".pdf" and ".png".
File generation can be turned off using set_output_mode("screen"). Output will be written to the /geneplots/ directory.}

\item{gene_name}{The name of the gene to create boxplot for. Must match a rowname in count data. String.}

\item{groupBy}{The group from target data that should be used to split data. Must match a rowname in count data.
I.E. "Timepoint" or "Infection". String.  If only two groups are present in this variable a wilcox test will be performed
between the two groups and the p-value will be displayed on the plot as well.}

\item{theme}{Theme for the layout and color scheme for the plot.  Valid selections are integers between 1-6.}

\item{returnData}{Boolean. Determines if this visualization should return data used to generate the visualization. Default=FALSE.}
}
\value{
If returnData is true, this function will return the long-form table of expression containing sample names,
categorical grouping, and sample IDs.
}
\description{
This function generates a box plot to display the expression of an individual gene with regard to a
specified grouping that can be based on any data that exists in the targets file. For example, a plot could
be created to view the expression of a gene as a function of different time points or experimental conditions.
}
\examples{
\dontrun{

#Plot the CAPS12 gene for each time point.
plot_gene(filename="CAPS12_time_plot.pdf", gene_name="CAPS12",
           groupBy="Time", theme=1, returnData=FALSE)

#Plot the METTL25 gene for each "response" group.  Store the long-form data table.
mettl25_dat <- plot_gene(filename="METTL25_time_plot.pdf", gene_name="METTL25",
                          groupBy="response", theme=2, returnData=TRUE)

}
}
\keyword{DE}
\keyword{boxplot}
\keyword{expression}
\keyword{gene}
\keyword{visualization}
