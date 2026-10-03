# =============================================================================
# Generate software bibliography with stable BibTeX keys
#
# Creates BibTeX citations for R and the R packages used in the manuscript.
#
# Output:
#   manuscript/references/software.bib
#
# Citation keys are generated automatically:
#   @R
#   @brms
#   @cmdstanr
#   @posterior
#   @tidybayes
#   ...
#
# If a package has multiple recommended citations, additional entries receive
# numeric suffixes, e.g. @brms2, @brms3.
# =============================================================================

packages <- c(
  "brms",
  "cmdstanr",
  "posterior",
  "tidybayes",
  "mgcv",
  "tidyverse",
  "ggplot2",
  "patchwork",
  "flextable"
)

# -----------------------------------------------------------------------------
# Output path
# -----------------------------------------------------------------------------

output_dir <- file.path(
  "manuscript",
  "references"
)

output_file <- file.path(
  output_dir,
  "software.bib"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# -----------------------------------------------------------------------------
# Check packages
# -----------------------------------------------------------------------------

missing_packages <- packages[
  !vapply(
    packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0) {
  stop(
    "These packages are not installed: ",
    paste(missing_packages, collapse = ", ")
  )
}

# -----------------------------------------------------------------------------
# Helper: generate BibTeX with stable citation keys
# -----------------------------------------------------------------------------

make_bib_entries <- function(citation_object, base_key) {
  
  # citation() returns an object that may contain one or more references
  entries <- lapply(
    seq_along(citation_object),
    function(i) {
      
      # Convert this individual citation to BibTeX
      bib <- as.character(
        toBibtex(citation_object[i])
      )
      
      # Use package name for first citation;
      # add numeric suffix for additional citations
      key <- if (i == 1L) {
        base_key
      } else {
        paste0(base_key, i)
      }
      
      # Replace whatever appears between the opening { and first comma
      # with our stable key
      bib[1] <- sub(
        "\\{[^,]*,",
        paste0("{", key, ","),
        bib[1]
      )
      
      bib
    }
  )
  
  unlist(entries, use.names = FALSE)
}

# -----------------------------------------------------------------------------
# Generate bibliography
# -----------------------------------------------------------------------------

bib <- c(
  "% ============================================================================",
  "% Software citations",
  "% Generated automatically from R and installed package citation metadata",
  "% ============================================================================",
  "",
  "% R",
  make_bib_entries(
    citation(),
    "R"
  ),
  ""
)

for (pkg in packages) {
  
  bib <- c(
    bib,
    paste0("% ", pkg),
    make_bib_entries(
      citation(pkg),
      pkg
    ),
    ""
  )
}

# -----------------------------------------------------------------------------
# Save bibliography
# -----------------------------------------------------------------------------

writeLines(
  bib,
  con = output_file,
  useBytes = TRUE
)

message(
  "Software bibliography written to: ",
  normalizePath(
    output_file,
    winslash = "/",
    mustWork = TRUE
  )
)

# -----------------------------------------------------------------------------
# Report citation keys
# -----------------------------------------------------------------------------

cat("\nSoftware bibliography generated.\n\n")
cat("Primary citation keys:\n")
cat("  @R\n")
cat(paste0("  @", packages, "\n"))

cat(
  "\nPackages with multiple recommended citations receive ",
  "numbered keys (for example, @brms2).\n"
)