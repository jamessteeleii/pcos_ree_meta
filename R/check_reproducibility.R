arguments <- commandArgs(trailingOnly = TRUE)
file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_argument) != 1L) {
  stop("Run this script with Rscript.", call. = FALSE)
}

script_path <- normalizePath(sub("^--file=", "", file_argument), mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
setwd(project_root)

source(file.path("R", "functions", "reproducibility_functions.R"), encoding = "UTF-8")

check_project_reproducibility(
  root = ".",
  check_environment = "--environment" %in% arguments
)
