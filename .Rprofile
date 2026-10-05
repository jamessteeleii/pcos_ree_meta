if (.Platform$OS.type == "windows") {
  locale_variables <- c("LANG", "LC_ALL", "LC_CTYPE")
  locale_values <- Sys.getenv(locale_variables, unset = NA_character_)
  unsupported_values <- c("C.UTF-8", "C.utf8")
  Sys.unsetenv(locale_variables[locale_values %in% unsupported_values])
}

source("renv/activate.R")
