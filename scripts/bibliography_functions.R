# =============================================================================
# Bibliography and software-version functions
# =============================================================================

write_package_bib <- function(
    packages,
    bib_file = "references/packages.bib",
    version_file = "tables/software_versions.csv",
    include_R = TRUE
) {
  
  # ---------------------------------------------------------------------------
  # Create output directories
  # ---------------------------------------------------------------------------
  
  dir.create(
    dirname(bib_file),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  dir.create(
    dirname(version_file),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  # ---------------------------------------------------------------------------
  # Collect citations and software versions
  # ---------------------------------------------------------------------------
  
  citations <- list()
  versions <- list()
  
  
  # R itself ------------------------------------------------------------------
  
  if (include_R) {
    
    citations[["R"]] <- citation()
    
    versions[["R"]] <- tibble::tibble(
      software = "R",
      version = paste0(
        R.version$major,
        ".",
        R.version$minor
      )
    )
  }
  
  
  # R packages ----------------------------------------------------------------
  
  for (pkg in packages) {
    
    if (!requireNamespace(pkg, quietly = TRUE)) {
      warning("Package not installed: ", pkg)
      next
    }
    
    citations[[pkg]] <- citation(pkg)
    
    versions[[pkg]] <- tibble::tibble(
      software = pkg,
      version = as.character(
        utils::packageVersion(pkg)
      )
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # Write BibTeX file
  # ---------------------------------------------------------------------------
  
  bibtex <- unlist(
    lapply(citations, toBibtex),
    use.names = FALSE
  )
  
  writeLines(
    unique(as.character(bibtex)),
    con = bib_file
  )
  
  
  # ---------------------------------------------------------------------------
  # Write software-version table
  # ---------------------------------------------------------------------------
  
  version_table <- dplyr::bind_rows(versions)
  
  readr::write_csv(
    version_table,
    version_file
  )
  
  
  # ---------------------------------------------------------------------------
  # Report
  # ---------------------------------------------------------------------------
  
  message(
    "Wrote ", length(citations),
    " software citations to: ",
    bib_file
  )
  
  message(
    "Wrote software version table to: ",
    version_file
  )
  
  
  # Return version table invisibly
  invisible(version_table)
}