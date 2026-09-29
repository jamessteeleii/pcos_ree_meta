#!/usr/bin/env Rscript

# Recreate the PRISMA flow diagram and append the September 2026 search update.
# This script intentionally uses only base R graphics packages so the figure can
# be regenerated without fitting models or restoring additional dependencies.

output_pdf <- file.path("output", "pdf", "prisma_updated_search_2026.pdf")
output_png <- file.path("plots", "prisma_updated_search_2026.png")
output_manuscript_png <- file.path("manuscript", "prisma.png")
grid <- asNamespace("grid")

dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_png), recursive = TRUE, showWarnings = FALSE)

colours <- list(
  ink = "#202124",
  border = "#3C4043",
  original = "#FBC02D",
  original_light = "#FFF8E1",
  other = "#D9D9D9",
  other_light = "#F4F4F4",
  update = "#7CC7C4",
  update_light = "#E7F6F5",
  phase = "#B9D2EE",
  final = "#DDEEDB"
)

draw_box <- function(x, y, width, height, label, fill = "white",
                     border = colours$border, fontsize = 8.2,
                     fontface = "plain", radius = 0.012) {
  grid$grid.roundrect(
    x = x,
    y = y,
    width = width,
    height = height,
    r = grid$unit(radius, "snpc"),
    gp = grid$gpar(fill = fill, col = border, lwd = 1.05)
  )
  grid$grid.text(
    label,
    x = x,
    y = y,
    gp = grid$gpar(
      col = colours$ink,
      fontsize = fontsize,
      fontfamily = "Arial",
      fontface = fontface,
      lineheight = 0.92
    )
  )
}

draw_header <- function(x, y, width, label, fill) {
  draw_box(
    x = x,
    y = y,
    width = width,
    height = 0.056,
    label = label,
    fill = fill,
    border = fill,
    fontsize = 9.2,
    fontface = "bold",
    radius = 0.02
  )
}

draw_phase <- function(y, height, label) {
  grid$grid.roundrect(
    x = 0.021,
    y = y,
    width = 0.026,
    height = height,
    r = grid$unit(0.012, "snpc"),
    gp = grid$gpar(fill = colours$phase, col = colours$phase)
  )
  grid$grid.text(
    label,
    x = 0.021,
    y = y,
    rot = 90,
    gp = grid$gpar(
      col = colours$ink,
      fontsize = 8.5,
      fontfamily = "Arial",
      fontface = "bold"
    )
  )
}

draw_arrow <- function(x0, y0, x1, y1) {
  grid$grid.lines(
    x = c(x0, x1),
    y = c(y0, y1),
    arrow = grid$arrow(type = "closed", length = grid$unit(0.105, "inches")),
    gp = grid$gpar(col = colours$border, lwd = 1.05)
  )
}

draw_polyline_arrow <- function(x, y) {
  grid$grid.lines(
    x = x,
    y = y,
    arrow = grid$arrow(type = "closed", length = grid$unit(0.105, "inches")),
    gp = grid$gpar(col = colours$border, lwd = 1.05)
  )
}

draw_prisma <- function() {
  grid$grid.newpage()
  grid$pushViewport(grid$viewport())

  grid$grid.text(
    "PRISMA 2020 flow diagram - updated literature search",
    x = 0.5,
    y = 0.974,
    gp = grid$gpar(
      col = colours$ink,
      fontsize = 15,
      fontfamily = "Arial",
      fontface = "bold"
    )
  )
  grid$grid.text(
    "The updated search identified no additional eligible studies.",
    x = 0.5,
    y = 0.944,
    gp = grid$gpar(
      col = colours$border,
      fontsize = 9.5,
      fontfamily = "Arial"
    )
  )

  draw_header(
    x = 0.202,
    y = 0.897,
    width = 0.315,
    label = "Original identification via databases and registers",
    fill = colours$original
  )
  draw_header(
    x = 0.535,
    y = 0.897,
    width = 0.245,
    label = "Original identification via other methods",
    fill = colours$other
  )
  draw_header(
    x = 0.832,
    y = 0.897,
    width = 0.305,
    label = "Updated search: January 2025-September 2026",
    fill = colours$update
  )

  draw_phase(y = 0.782, height = 0.175, label = "Identification")
  draw_phase(y = 0.506, height = 0.325, label = "Screening")
  draw_phase(y = 0.121, height = 0.145, label = "Included")

  # Original database/register search.
  draw_box(
    x = 0.125,
    y = 0.784,
    width = 0.155,
    height = 0.145,
    label = paste(
      "Records identified from:",
      "Databases (n = 3)",
      "PubMed (n = 230)",
      "MEDLINE/EBSCO (n = 195)",
      "Web of Science (n = 71)",
      "Registers (n = 0)",
      sep = "\n"
    ),
    fill = colours$original_light,
    fontsize = 7.7
  )
  draw_box(
    x = 0.289,
    y = 0.784,
    width = 0.155,
    height = 0.122,
    label = paste(
      "Records removed before screening:",
      "Duplicate records (n = 130)",
      "Marked ineligible by automation (n = 0)",
      "Removed for other reasons (n = 0)",
      sep = "\n"
    ),
    fill = colours$original_light,
    fontsize = 7.5
  )
  draw_box(
    x = 0.125,
    y = 0.642,
    width = 0.155,
    height = 0.064,
    label = "Records screened\n(n = 366)",
    fontsize = 8.2
  )
  draw_box(
    x = 0.289,
    y = 0.642,
    width = 0.155,
    height = 0.064,
    label = "Records excluded\n(n = 337)",
    fontsize = 8.2
  )
  draw_box(
    x = 0.125,
    y = 0.521,
    width = 0.155,
    height = 0.064,
    label = "Reports sought for retrieval\n(n = 27)",
    fontsize = 8.2
  )
  draw_box(
    x = 0.289,
    y = 0.521,
    width = 0.155,
    height = 0.064,
    label = "Reports not retrieved\n(n = 0)",
    fontsize = 8.2
  )
  draw_box(
    x = 0.125,
    y = 0.395,
    width = 0.155,
    height = 0.064,
    label = "Reports assessed for eligibility\n(n = 27)",
    fontsize = 8.2
  )
  draw_box(
    x = 0.289,
    y = 0.395,
    width = 0.155,
    height = 0.082,
    label = paste(
      "Reports excluded:",
      "Incorrect outcome (n = 10)",
      "Full text unretrievable (n = 1)",
      sep = "\n"
    ),
    fontsize = 7.7
  )

  draw_arrow(0.203, 0.784, 0.211, 0.784)
  draw_arrow(0.125, 0.712, 0.125, 0.676)
  draw_arrow(0.203, 0.642, 0.211, 0.642)
  draw_arrow(0.125, 0.610, 0.125, 0.553)
  draw_arrow(0.203, 0.521, 0.211, 0.521)
  draw_arrow(0.125, 0.489, 0.125, 0.427)
  draw_arrow(0.203, 0.395, 0.211, 0.395)

  # Original records identified through other methods.
  draw_box(
    x = 0.482,
    y = 0.784,
    width = 0.126,
    height = 0.082,
    label = "Records identified from:\nWebsites (n = 2)",
    fill = colours$other_light,
    border = colours$other,
    fontsize = 8.0
  )
  draw_box(
    x = 0.482,
    y = 0.521,
    width = 0.126,
    height = 0.064,
    label = "Reports sought for retrieval\n(n = 2)",
    fill = colours$other_light,
    border = colours$other,
    fontsize = 8.0
  )
  draw_box(
    x = 0.603,
    y = 0.521,
    width = 0.103,
    height = 0.064,
    label = "Reports not retrieved\n(n = NA)",
    fill = colours$other_light,
    border = colours$other,
    fontsize = 7.7
  )
  draw_box(
    x = 0.482,
    y = 0.395,
    width = 0.126,
    height = 0.064,
    label = "Reports assessed for eligibility\n(n = 2)",
    fill = colours$other_light,
    border = colours$other,
    fontsize = 8.0
  )
  draw_box(
    x = 0.603,
    y = 0.395,
    width = 0.103,
    height = 0.064,
    label = "Reports excluded\n(n = NA)",
    fill = colours$other_light,
    border = colours$other,
    fontsize = 7.7
  )

  draw_arrow(0.482, 0.743, 0.482, 0.553)
  draw_arrow(0.545, 0.521, 0.551, 0.521)
  draw_arrow(0.482, 0.489, 0.482, 0.427)
  draw_arrow(0.545, 0.395, 0.551, 0.395)

  # Updated database search supplied by the review team.
  draw_box(
    x = 0.758,
    y = 0.784,
    width = 0.144,
    height = 0.134,
    label = paste(
      "Records identified (n = 69):",
      "PubMed (n = 38)",
      "Web of Science (n = 13)",
      "MEDLINE/EBSCO (n = 18)",
      sep = "\n"
    ),
    fill = colours$update_light,
    fontsize = 8.0
  )
  draw_box(
    x = 0.910,
    y = 0.784,
    width = 0.132,
    height = 0.104,
    label = paste(
      "Records removed before screening:",
      "Duplicate or already in the",
      "original search (n = 8)",
      sep = "\n"
    ),
    fill = colours$update_light,
    fontsize = 7.7
  )
  draw_box(
    x = 0.758,
    y = 0.642,
    width = 0.144,
    height = 0.064,
    label = "Records screened by title\n(n = 61)",
    fontsize = 8.2
  )
  draw_box(
    x = 0.910,
    y = 0.642,
    width = 0.132,
    height = 0.064,
    label = "Records excluded by title\n(n = 59)",
    fontsize = 8.0
  )
  draw_box(
    x = 0.758,
    y = 0.521,
    width = 0.144,
    height = 0.064,
    label = "Abstracts screened\n(n = 2)",
    fontsize = 8.2
  )
  draw_box(
    x = 0.910,
    y = 0.521,
    width = 0.132,
    height = 0.064,
    label = "Records excluded by abstract\n(n = 2)",
    fontsize = 8.0
  )
  draw_box(
    x = 0.758,
    y = 0.395,
    width = 0.144,
    height = 0.064,
    label = "Reports sought for retrieval\n(n = 0)",
    fontsize = 8.2
  )
  draw_box(
    x = 0.910,
    y = 0.395,
    width = 0.132,
    height = 0.064,
    label = "Reports assessed for eligibility\n(n = 0)",
    fontsize = 8.0
  )

  draw_arrow(0.830, 0.784, 0.844, 0.784)
  draw_arrow(0.758, 0.717, 0.758, 0.676)
  draw_arrow(0.830, 0.642, 0.844, 0.642)
  draw_arrow(0.758, 0.610, 0.758, 0.553)
  draw_arrow(0.830, 0.521, 0.844, 0.521)
  draw_arrow(0.758, 0.489, 0.758, 0.427)
  draw_arrow(0.830, 0.395, 0.844, 0.395)

  # Final evidence base. The original pathways and update all converge here.
  draw_box(
    x = 0.500,
    y = 0.121,
    width = 0.310,
    height = 0.116,
    label = paste(
      "Studies included in review (n = 17)",
      "Reports of included studies (n = 18)",
      "New studies included from updated search (n = 0)",
      sep = "\n"
    ),
    fill = colours$final,
    fontsize = 9.2,
    fontface = "bold"
  )

  draw_polyline_arrow(
    x = c(0.125, 0.125, 0.345),
    y = c(0.363, 0.180, 0.180)
  )
  draw_polyline_arrow(
    x = c(0.482, 0.482, 0.482),
    y = c(0.363, 0.226, 0.180)
  )
  draw_polyline_arrow(
    x = c(0.758, 0.758, 0.655),
    y = c(0.363, 0.180, 0.180)
  )

  grid$grid.text(
    "Original-search counts are reproduced from the current PRISMA diagram. Update-search counts were supplied by the review team.",
    x = 0.5,
    y = 0.030,
    gp = grid$gpar(
      col = colours$border,
      fontsize = 7.5,
      fontfamily = "Arial",
      fontface = "italic"
    )
  )

  grid$popViewport()
}

grDevices::cairo_pdf(output_pdf, width = 16.54, height = 11.69, family = "Arial")
draw_prisma()
grDevices::dev.off()

grDevices::png(
  output_png,
  width = 4961,
  height = 3508,
  units = "px",
  res = 300,
  type = "cairo"
)
draw_prisma()
grDevices::dev.off()

if (!file.copy(output_png, output_manuscript_png, overwrite = TRUE)) {
  stop("Could not copy the updated PRISMA diagram into the manuscript directory.")
}

message("Created: ", output_pdf)
message("Created: ", output_png)
message("Updated: ", output_manuscript_png)
