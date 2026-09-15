ensure_parent_directory <- function(path) {
  directory <- dirname(path)
  if (!dir.exists(directory)) {
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  }
  path
}

save_plot_tiff <- function(plot, path, width, height, dpi = 300) {
  path <- ensure_parent_directory(path)
  ggplot2::ggsave(
    filename = path,
    plot = plot,
    dpi = dpi,
    device = "tiff",
    width = width,
    height = height
  )
  path
}
