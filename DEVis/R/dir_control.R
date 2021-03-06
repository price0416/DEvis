#' Initialize the directory structure for automatically storing and structuring results.
#'
#' This function initializes and, if necessary, creates the directory structure for output of
#' analysis based on a provided base directory.  A /results/ folder will be created in the base directory
#' containing a complete directory structure for all possible results that can be generated by this package.
#' All visualizations and data aggregation tools rely on this directory structure for data output.
#' Note that this function will NOT overwrite existing directories.  Initialization of
#' data directories should be performed separately using init_data_paths().
#' @param base_dir Base directory in which to build directory structure for result output.
#' @return This function does not return a value.
#' @seealso \code{\link{init_data_paths}}, \code{\link{set_output_mode}}
#' @keywords directory organization project-management output initialization
#' @export
#' @examples
#' \dontrun{
#'
#' #Initialize a DEVis directory structure at a base directory.
#' create_dir_struct(base_dir="/Users/adam/projects/ebola/")
#'
#' This command will create the following folders:
#' /Users/adam/projects/ebola/
#'   /results/
#'     /DE/
#'       /boxplot/
#'       /counts/
#'       /data/
#'       /density_plots/
#'       /divergence/
#'       /heatmaps/
#'       /profile_plots/
#'       /series_plots/
#'       /volcano/
#'    /dendrograms/
#'    /geneplots/
#'    /group_stats/
#'    /MDS/
#'      /standard/
#'      /hulls/
#'    /sample_distance/
#'      /euclidian/
#'      /poisson/
#'
#' }
create_dir_struct <- function(base_dir)
{
  if(!dir.exists(base_dir))
  {
    print("Error:")
    print(base_dir)
    stop("Provided base directory doesn't exist.")
    return(-1)
  }

  print("Initializing directory structure...")

  #Establish Paths.
  working_dir    <- base_dir
  results_dir    <- paste(working_dir,"results/", sep="")
  DE_base_dir    <- paste(results_dir,"DE/", sep="")
  DE_dir         <- paste(results_dir,"DE/data", sep="")
  DE_hm_base_dir <- paste(results_dir,"DE/heatmaps/", sep="")
  DE_count_dir   <- paste(results_dir,"DE/counts", sep="")
  MDS_base       <- paste(results_dir,"MDS/", sep="")
  MDS_dir        <- paste(results_dir,"MDS/standard/", sep="")
  MDS_hull_dir   <- paste(results_dir,"MDS/hulls/", sep="")
  distance_dir   <- paste(results_dir,"sample_distance/", sep="")
  euclidian_dir  <- paste(results_dir,"sample_distance/euclidian/", sep="")
  poisson_dir    <- paste(results_dir,"sample_distance/poisson/", sep="")
  gene_plot_dir  <- paste(results_dir,"geneplots/", sep="")
  group_stat_dir <- paste(results_dir,"group_stats/", sep="")
  DE_diverge_dir <- paste(results_dir,"DE/divergence/", sep="")
  DE_boxplot_dir <- paste(results_dir,"DE/boxplot/", sep="")
  DE_profile_dir <- paste(results_dir,"DE/profile_plots/", sep="")
  DE_density_dir <- paste(results_dir,"DE/density_plots/", sep="")
  DE_series_dir  <- paste(results_dir,"DE/series_plots/", sep="")
  DE_volcano_dir <- paste(results_dir,"DE/volcano/", sep="")
  dendro_dir     <- paste(results_dir,"dendrograms/", sep="")

  #Save paths to the environment.
  assign("working_dir", base_dir, envir=.DEVis_env)
  assign("results_dir", paste(working_dir,"results/", sep=""), envir=.DEVis_env)
  assign("DE_base_dir", paste(results_dir,"DE/", sep=""), envir=.DEVis_env)
  assign("DE_dir", paste(results_dir,"DE/data", sep=""), envir=.DEVis_env)
  assign("DE_hm_base_dir", paste(results_dir,"DE/heatmaps/", sep=""), envir=.DEVis_env)
  assign("DE_count_dir", paste(results_dir,"DE/counts", sep=""), envir=.DEVis_env)
  assign("MDS_base", paste(results_dir,"MDS/", sep=""), envir=.DEVis_env)
  assign("MDS_dir", paste(results_dir,"MDS/standard/", sep=""), envir=.DEVis_env)
  assign("MDS_hull_dir", paste(results_dir,"MDS/hulls/", sep=""), envir=.DEVis_env)
  assign("distance_dir", paste(results_dir,"sample_distance/", sep=""), envir=.DEVis_env)
  assign("euclidian_dir", paste(results_dir,"sample_distance/euclidian/", sep=""), envir=.DEVis_env)
  assign("poisson_dir", paste(results_dir,"sample_distance/poisson/", sep=""), envir=.DEVis_env)
  assign("gene_plot_dir", paste(results_dir,"geneplots/", sep=""), envir=.DEVis_env)
  assign("group_stat_dir", paste(results_dir,"group_stats/", sep=""), envir=.DEVis_env)
  assign("DE_diverge_dir", paste(results_dir,"DE/divergence/", sep=""), envir=.DEVis_env)
  assign("DE_boxplot_dir", paste(results_dir,"DE/boxplot/", sep=""), envir=.DEVis_env)
  assign("DE_profile_dir", paste(results_dir,"DE/profile_plots/", sep=""), envir=.DEVis_env)
  assign("DE_density_dir", paste(results_dir,"DE/density_plots/", sep=""), envir=.DEVis_env)
  assign("DE_series_dir", paste(results_dir,"DE/series_plots/", sep=""), envir=.DEVis_env)
  assign("DE_volcano_dir", paste(results_dir,"DE/volcano/", sep=""), envir=.DEVis_env)
  assign("dendro_dir", paste(results_dir,"dendrograms/", sep=""), envir=.DEVis_env)


  #Create directory structure if it doesn't exist.
  dir.create(working_dir, showWarnings = FALSE)
  dir.create(results_dir, showWarnings = FALSE)
  dir.create(DE_base_dir, showWarnings = FALSE)
  dir.create(DE_hm_base_dir, showWarnings = FALSE)
  dir.create(DE_dir, showWarnings = FALSE)
  dir.create(DE_count_dir, showWarnings = FALSE)
  dir.create(MDS_base, showWarnings = FALSE)
  dir.create(MDS_dir, showWarnings = FALSE)
  dir.create(MDS_hull_dir, showWarnings = FALSE)
  dir.create(distance_dir, showWarnings = FALSE)
  dir.create(euclidian_dir, showWarnings = FALSE)
  dir.create(poisson_dir, showWarnings = FALSE)
  dir.create(gene_plot_dir, showWarnings = FALSE)
  dir.create(group_stat_dir, showWarnings = FALSE)
  dir.create(DE_diverge_dir, showWarnings = FALSE)
  dir.create(DE_boxplot_dir, showWarnings = FALSE)
  dir.create(DE_profile_dir, showWarnings = FALSE)
  dir.create(DE_density_dir, showWarnings = FALSE)
  dir.create(dendro_dir, showWarnings = FALSE)
  dir.create(DE_series_dir, showWarnings = FALSE)
  dir.create(DE_volcano_dir, showWarnings = FALSE)

  setwd(.DEVis_env$working_dir)
}


#' Initialize data paths for count and target files.
#'
#' This function initializes the paths to folders containing count and target data.  This package
#' assumes correspondence between count data column names and target data row names.  In this package,
#' metadata files are referred to as target data.
#' Initialization by this function is required by most functions.
#' @param count_path Path to directory containing count data. Count data is assumed to be tab-delimited text
#' represent raw expression counts for each gene, with samples indicated by column and genes identified by row.
#' Column names in count data must correspond to row names in target files.
#' @param target_path Path to directory containing target data. Target data can be tab-delimited or comma-separated text.
#' Target rownames must correspond to count column names.  Additionally, any data can be included in target files to
#' represent experimental conditions or additional metadata for the project.
#' @return This function does not return a value.
#' @keywords initialization directory organization project-management input
#' @seealso \code{\link{create_dir_struct}}
#' @export
#' @examples
#' \dontrun{
#'
#' #Initialize paths to count and target file directories.
#' init_data_paths(count_path="/Users/adam/projects/ebola/counts/",
#'                 target_path="/Users/adam/projects/ebola/")
#'
#' }
init_data_paths <- function(count_path, target_path)
{
  if(!dir.exists(count_path))
  {
    stop("Error: Provided count directory doesn't exist.")
    return(-1)
  }
  if(!dir.exists(target_path))
  {
    stop("Error: Provided target directory doesn't exist.")
    return(-1)
  }
  print("Initializing data paths...")
  assign("counts_dir", count_path, envir=.DEVis_env)
  assign("targets_dir", target_path, envir=.DEVis_env)
}
