manifest_field <- function(lines, label) {
  pattern <- paste0("^- ", label, ":\\s*`?([^`]+)`?\\s*$")
  match <- grep(pattern, lines, value = TRUE)
  if (length(match) != 1L) {
    stop("Data manifest must contain exactly one '", label, "' field.", call. = FALSE)
  }
  sub(pattern, "\\1", match, perl = TRUE)
}

sha256_file <- function(path) {
  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE))
  }

  if (nzchar(Sys.which("certutil"))) {
    output <- system2(
      "certutil",
      c("-hashfile", shQuote(normalizePath(path, winslash = "\\", mustWork = TRUE)), "SHA256"),
      stdout = TRUE,
      stderr = TRUE
    )
    candidates <- gsub("\\s", "", output)
    candidates <- candidates[grepl("^[[:xdigit:]]{64}$", candidates)]
    if (length(candidates)) {
      return(tolower(candidates[[1L]]))
    }
  }

  for (command in c("sha256sum", "shasum")) {
    if (nzchar(Sys.which(command))) {
      arguments <- if (identical(command, "shasum")) c("-a", "256", path) else path
      output <- system2(command, arguments, stdout = TRUE, stderr = TRUE)
      candidate <- sub("\\s+.*$", "", output[[1L]])
      if (grepl("^[[:xdigit:]]{64}$", candidate)) {
        return(tolower(candidate))
      }
    }
  }

  NA_character_
}

check_data_integrity <- function(
    data_file = file.path("data", "studies_data.csv"),
    manifest_file = file.path("data", "README.md")) {
  problems <- character()
  expect_true <- function(condition, message) {
    if (!isTRUE(condition)) {
      problems <<- c(problems, message)
    }
  }

  expect_true(file.exists(data_file), paste("Missing canonical dataset:", data_file))
  expect_true(file.exists(manifest_file), paste("Missing data manifest:", manifest_file))
  if (length(problems)) {
    return(problems)
  }

  manifest <- readLines(manifest_file, encoding = "UTF-8", warn = FALSE)
  expected_rows <- suppressWarnings(as.integer(manifest_field(manifest, "Records")))
  expected_columns <- suppressWarnings(as.integer(manifest_field(manifest, "Variables")))
  expected_hash <- tolower(manifest_field(manifest, "SHA-256"))

  bytes <- readBin(data_file, what = "raw", n = file.info(data_file)$size)
  text <- rawToChar(bytes)
  expect_true(validUTF8(text), "Canonical dataset is not valid UTF-8.")
  expect_true(
    !grepl("\uFFFD|\\?2B|\\?3", text, perl = TRUE),
    "Canonical dataset contains a replacement character or a known corrupted Greek title."
  )

  data <- tryCatch(
    read.csv(
      data_file,
      colClasses = "character",
      check.names = FALSE,
      na.strings = character(),
      encoding = "UTF-8"
    ),
    error = function(error) error
  )
  if (inherits(data, "error")) {
    return(c(problems, paste("Could not parse canonical CSV:", conditionMessage(data))))
  }

  expect_true(nrow(data) == expected_rows, paste("Expected", expected_rows, "records; found", nrow(data)))
  expect_true(ncol(data) == expected_columns, paste("Expected", expected_columns, "variables; found", ncol(data)))
  expect_true(!anyDuplicated(names(data)), "Canonical dataset has duplicate column names.")

  required_columns <- c(
    "lab", "study", "arm", "authors", "title", "year", "cond", "timepoint",
    "insulin_resistant_y_n", "method", "n", "units", "mean", "sd", "comments"
  )
  expect_true(
    all(required_columns %in% names(data)),
    paste("Missing required columns:", paste(setdiff(required_columns, names(data)), collapse = ", "))
  )

  if (all(c("lab", "study", "arm", "timepoint") %in% names(data))) {
    keys <- do.call(paste, c(data[c("lab", "study", "arm", "timepoint")], sep = "\r"))
    expect_true(!anyDuplicated(keys), "Duplicate lab/study/arm/timepoint records detected.")
  }

  empty_rows <- apply(data, 1L, function(row) all(is.na(row) | !nzchar(trimws(row))))
  expect_true(!any(empty_rows), "Canonical dataset contains an empty record.")

  textual_missing <- vapply(
    data,
    function(column) any(toupper(trimws(column)) %in% c("NA", "N/A", "NULL"), na.rm = TRUE),
    logical(1)
  )
  expect_true(
    !any(textual_missing),
    paste("Textual missing-value codes found in:", paste(names(textual_missing)[textual_missing], collapse = ", "))
  )

  numeric_columns <- grep(
    "^(lab|study|arm|year|n|m|median|sd|se|lower_range|upper_range|iqr)($|_)",
    names(data),
    value = TRUE
  )
  for (column in numeric_columns) {
    present <- !is.na(data[[column]]) & nzchar(trimws(data[[column]]))
    parsed <- suppressWarnings(as.numeric(data[[column]][present]))
    expect_true(
      !anyNA(parsed),
      paste("Non-numeric value found in numeric column:", column)
    )
  }

  positive_columns <- grep("^(n|sd)($|_)", names(data), value = TRUE)
  for (column in positive_columns) {
    present <- !is.na(data[[column]]) & nzchar(trimws(data[[column]]))
    parsed <- suppressWarnings(as.numeric(data[[column]][present]))
    expect_true(
      !any(parsed <= 0, na.rm = TRUE),
      paste("Non-positive value found in:", column)
    )
  }

  allowed_values <- list(
    cond = c("PCOS", "Control"),
    timepoint = c("baseline", "mid_intervention", "post_intervention"),
    insulin_resistant_y_n = c("", "y", "n"),
    in_macrofactor_article = c("y", "n")
  )
  for (column in intersect(names(allowed_values), names(data))) {
    unexpected <- setdiff(unique(data[[column]]), allowed_values[[column]])
    expect_true(
      !length(unexpected),
      paste("Unexpected value in", column, ":", paste(unexpected, collapse = ", "))
    )
  }

  actual_hash <- sha256_file(data_file)
  if (is.na(actual_hash)) {
    problems <- c(problems, "No SHA-256 implementation is available to verify the data checksum.")
  } else {
    expect_true(
      identical(tolower(actual_hash), expected_hash),
      paste("Data checksum differs from data/README.md. Expected", expected_hash, "but found", actual_hash)
    )
  }

  problems
}

extract_bibliography_keys <- function(path) {
  lines <- readLines(path, encoding = "UTF-8", warn = FALSE)
  entries <- grep("^@[[:alpha:]]+\\s*\\{", lines, value = TRUE)
  sub("^@[[:alpha:]]+\\s*\\{\\s*([^,]+),.*$", "\\1", entries, perl = TRUE)
}

extract_citation_keys <- function(path) {
  text <- paste(readLines(path, encoding = "UTF-8", warn = FALSE), collapse = "\n")
  matches <- gregexpr("(?<![[:alnum:]_])@[[:alnum:]_.:-]+", text, perl = TRUE)
  citations <- regmatches(text, matches)[[1L]]
  citations <- sub("^@", "", citations)
  citations[!grepl("^(fig|tbl|sec|eq|lst)-", citations)]
}

check_project_structure <- function(root = ".") {
  problems <- character()
  expect_true <- function(condition, message) {
    if (!isTRUE(condition)) {
      problems <<- c(problems, message)
    }
  }

  required_files <- c(
    "AGENTS.md", "CLAUDE.md", ".gitignore", "renv.lock", "_targets.R", "_targets.yaml",
    file.path("data", "README.md"), file.path("data", "studies_data.csv"),
    file.path("docs", "REPRODUCIBILITY.md")
  )
  missing_files <- required_files[!file.exists(file.path(root, required_files))]
  expect_true(!length(missing_files), paste("Missing project files:", paste(missing_files, collapse = ", ")))

  forbidden_files <- c(
    ".Rhistory",
    file.path("data", "studies_data.xlsx"),
    file.path("data", "Extraction notes.docx"),
    file.path("manuscript", "JCEM author affirmation.docx"),
    file.path("manuscript", "JCEM author affirmation.pdf"),
    file.path("manuscript", "_targets.yaml")
  )
  present_forbidden <- forbidden_files[file.exists(file.path(root, forbidden_files))]
  expect_true(
    !length(present_forbidden),
    paste("Retired or private files are present:", paste(present_forbidden, collapse = ", "))
  )
  expect_true(!dir.exists(file.path(root, "papers")), "The non-redistributable papers directory is present.")

  ignore_path <- file.path(root, ".gitignore")
  if (file.exists(ignore_path)) {
    ignore_lines <- trimws(readLines(ignore_path, encoding = "UTF-8", warn = FALSE))
    required_ignores <- c("_targets/", "renv/library/", ".Rhistory", "papers/", "*_files/")
    missing_ignores <- setdiff(required_ignores, ignore_lines)
    expect_true(!length(missing_ignores), paste("Missing ignore rules:", paste(missing_ignores, collapse = ", ")))
  }

  lock_path <- file.path(root, "renv.lock")
  if (file.exists(lock_path)) {
    lock <- paste(readLines(lock_path, encoding = "UTF-8", warn = FALSE), collapse = "\n")
    expect_true(
      grepl('"R"\\s*:\\s*\\{\\s*"Version"\\s*:\\s*"4\\.4\\.3"', lock, perl = TRUE),
      "renv.lock must declare R 4.4.3."
    )
    locked_packages <- c("targets", "qs", "tarchetypes", "quarto")
    missing_packages <- locked_packages[!vapply(
      locked_packages,
      function(package) grepl(paste0('"', package, '"\\s*:'), lock, perl = TRUE),
      logical(1)
    )]
    expect_true(!length(missing_packages), paste("Missing lockfile packages:", paste(missing_packages, collapse = ", ")))
    expect_true(!grepl("packagemanager\\.posit\\.co/cran/latest", lock), "renv.lock uses the moving cran/latest repository.")
  }

  code_files <- c(
    file.path(root, "_targets.R"),
    file.path(root, "_targets.yaml"),
    list.files(file.path(root, "R"), pattern = "\\.[Rr]$", recursive = TRUE, full.names = TRUE),
    list.files(file.path(root, "manuscript"), pattern = "\\.(qmd|ya?ml)$", recursive = TRUE, full.names = TRUE)
  )
  code_files <- code_files[file.exists(code_files)]
  absolute_path_hits <- vapply(code_files, function(path) {
    text <- paste(readLines(path, encoding = "UTF-8", warn = FALSE), collapse = "\n")
    grepl("(^|[^[:alpha:]])[[:alpha:]]:[/\\\\]", text, perl = TRUE)
  }, logical(1))
  expect_true(
    !any(absolute_path_hits),
    paste("Windows absolute path found in:", paste(code_files[absolute_path_hits], collapse = ", "))
  )

  if (exists("lock", inherits = FALSE)) {
    code <- paste(vapply(code_files, function(path) {
      paste(readLines(path, encoding = "UTF-8", warn = FALSE), collapse = "\n")
    }, character(1)), collapse = "\n")
    namespace_calls <- regmatches(
      code,
      gregexpr("(?<![[:alnum:].])[[:alpha:]][[:alnum:].]*::", code, perl = TRUE)
    )[[1L]]
    library_calls <- regmatches(
      code,
      gregexpr("\\b(?:library|require)\\s*\\(\\s*['\"]?[[:alpha:]][[:alnum:].]*", code, perl = TRUE)
    )[[1L]]
    referenced_packages <- c(
      sub("::$", "", namespace_calls),
      sub("^.*\\(\\s*['\"]?", "", library_calls, perl = TRUE)
    )
    referenced_packages <- unique(referenced_packages[nzchar(referenced_packages)])
    base_packages <- c("base", "datasets", "graphics", "grDevices", "methods", "parallel", "stats", "tools", "utils")
    locked_entries <- regmatches(lock, gregexpr('(?m)^    "[^"]+": \\{', lock, perl = TRUE))[[1L]]
    locked_packages <- sub('^    "([^"]+)": \\{$', "\\1", locked_entries, perl = TRUE)
    unlocked_packages <- setdiff(referenced_packages, c(base_packages, locked_packages))
    expect_true(
      !length(unlocked_packages),
      paste("R packages referenced in code but missing from renv.lock:", paste(unlocked_packages, collapse = ", "))
    )
  }

  bibliography <- file.path(root, "manuscript", "mylibrary.bib")
  qmd_files <- list.files(file.path(root, "manuscript"), pattern = "\\.qmd$", full.names = TRUE)
  if (file.exists(bibliography) && length(qmd_files)) {
    bibliography_keys <- extract_bibliography_keys(bibliography)
    citation_keys <- unique(unlist(lapply(qmd_files, extract_citation_keys), use.names = FALSE))
    missing_citations <- setdiff(citation_keys, bibliography_keys)
    expect_true(
      !length(missing_citations),
      paste("Citation keys missing from manuscript/mylibrary.bib:", paste(missing_citations, collapse = ", "))
    )
  }

  problems
}

check_local_environment <- function(root = ".") {
  problems <- character()
  expected_r <- "4.4.3"
  if (!identical(as.character(getRversion()), expected_r)) {
    problems <- c(problems, paste("Expected R", expected_r, "but found", getRversion()))
  }

  required_packages <- c(
    "targets", "qs", "tarchetypes", "quarto", "tidyverse", "here", "metafor",
    "brms", "marginaleffects", "tidybayes", "patchwork", "flextable",
    "officer", "webshot2"
  )
  missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages)) {
    problems <- c(problems, paste("Packages not installed:", paste(missing_packages, collapse = ", ")))
  }

  commands <- c("git", "quarto", "xelatex")
  missing_commands <- commands[!nzchar(Sys.which(commands))]
  if (.Platform$OS.type == "windows") {
    build_commands <- c("make", "g++")
    missing_commands <- c(missing_commands, build_commands[!nzchar(Sys.which(build_commands))])
  }
  if (length(missing_commands)) {
    problems <- c(problems, paste("System commands not found:", paste(unique(missing_commands), collapse = ", ")))
  }

  if (requireNamespace("webshot2", quietly = TRUE)) {
    chrome <- tryCatch(webshot2::find_chrome(), error = function(error) "")
    if (!nzchar(chrome)) {
      problems <- c(problems, "Chrome or Chromium was not found for webshot2.")
    }
  }

  problems
}

check_project_reproducibility <- function(
    data_file = file.path("data", "studies_data.csv"),
    manifest_file = file.path("data", "README.md"),
    root = ".",
    check_environment = FALSE) {
  problems <- c(
    check_data_integrity(data_file, manifest_file),
    check_project_structure(root)
  )
  if (isTRUE(check_environment)) {
    problems <- c(problems, check_local_environment(root))
  }
  problems <- unique(problems[nzchar(problems)])

  if (length(problems)) {
    stop(
      paste(c("Reproducibility checks failed:", paste0("- ", problems)), collapse = "\n"),
      call. = FALSE
    )
  }

  message("Reproducibility checks passed.")
  invisible(TRUE)
}
