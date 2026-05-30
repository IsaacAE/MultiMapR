################################################################################
# load_data.R
#
# Utilities to load and prepare phylogenetic data in MultiMapR.
#
# MAIN FUNCTION:
#   load_data(tree_path, csv_path, ...)  →  list(tree, characters)
#
# SUPPORTED TREE FORMATS:
#   .tre / .tree / .nwk  →  read.tree()   (Newick)
#   .nex / .nexus        →  read.nexus()  (NEXUS)
#   (automatic detection by extension; can be forced with tree_format=)
#
# SUPPORTED CSV FORMATS:
#   - With or without a header for characters (automatic detection)
#   - First column: species names (any column name)
#   - Numeric and text (string) characters
#   - Inapplicable "-" and unknown "?" are preserved as-is
#   - Species names with spaces → optionally normalized to "_"
#   - UTF-8 BOM (files exported from Excel) removed automatically
#
# NOTE: This file is part of the MultiMapR package. Dependencies
# (ape, tools) are declared in DESCRIPTION; library() is not used here.
################################################################################


# ============================================================================ #
# BOM HANDLING
# ============================================================================ #

#' Removes the UTF-8 BOM (EF BB BF) if present and returns the path to use
#'
#' If the file has a BOM, copies the content without the first 3 bytes to a
#' temporary file and returns its path. If there is no BOM, returns the original path.
#' The temporary file is automatically deleted when the R session ends.
#'
#' @param csv_path Path to the original file.
#' @return Path to the clean file (temporary or original).
.remove_bom <- function(csv_path) {
  con <- file(csv_path, open = "rb")
  bom <- readBin(con, raw(), n = 3)
  has_bom <- (length(bom) == 3 &&
                bom[1] == as.raw(0xEF) &&
                bom[2] == as.raw(0xBB) &&
                bom[3] == as.raw(0xBF))

  if (!has_bom) {
    close(con)
    return(csv_path)
  }

  # The 3 BOM bytes have already been consumed; read the rest
  content <- readBin(con, raw(), n = file.info(csv_path)$size)
  close(con)

  tmp <- tempfile(fileext = ".csv")
  con_out <- file(tmp, open = "wb")
  writeBin(content, con_out)
  close(con_out)

  message("UTF-8 BOM detected and removed for reading.")
  tmp
}


# ============================================================================ #
# AUTOMATIC HEADER DETECTION
# ============================================================================ #

#' Decides whether a CSV has a header row for the character columns
#'
#' Heuristic based on the species column (column 1):
#' taxon names always contain a space or underscore
#' ("Homo_sapiens", "Homo sapiens"). If the first cell has
#' neither, it is treated as a column label -> there is a header.
#'
#'   "Species"               -> header = TRUE   (no space or "_")
#'   "Taxon"                 -> header = TRUE
#'   "Saccopteryx bilineata" -> header = FALSE  (has space)
#'   "Passer_domesticus"     -> header = FALSE  (has "_")
#'
#' @param csv_path  Path to the CSV file (already BOM-free).
#' @param sep       Field separator (default ",").
#' @return TRUE if a header is detected, FALSE otherwise.
.detect_header <- function(csv_path, sep = ",") {
  first_row <- read.csv(csv_path, header = FALSE, sep = sep,
                        stringsAsFactors = FALSE, nrows = 1,
                        colClasses = "character")

  if (nrow(first_row) < 1 || ncol(first_row) < 1) return(FALSE)

  species_cell <- trimws(as.character(first_row[1, 1]))

  # A taxon name has a space or underscore ("Homo_sapiens",
  # "Homo sapiens"). A column label ("Species", "Taxon") does not.
  looks_like_taxon <- grepl("[ _]", species_cell)
  !looks_like_taxon
}


# ============================================================================ #
# COLUMN NAME NORMALIZATION
# ============================================================================ #

#' Generates R-safe character names
#'
#' If the CSV has a header, those names are used (sanitized).
#' If there is no header, generates "char1", "char2", ...
#'
#' @param raw_names  Vector with the raw names of the character columns.
#'                   NULL or empty → automatic names are generated.
#' @param n          Number of character columns.
#' @return Character vector of length n.
.normalize_char_names <- function(raw_names, n) {
  if (is.null(raw_names) || length(raw_names) == 0) {
    return(paste0("char", seq_len(n)))
  }
  names   <- make.names(trimws(raw_names), unique = TRUE)
  empty   <- nchar(names) == 0 | grepl("^V[0-9]+$", names)
  names[empty] <- paste0("char", which(empty))
  names
}


# ============================================================================ #
# TREE READING
# ============================================================================ #

#' Sanitizes a tree file for compatibility with ape
#'
#' Removes common problems produced by different phylogenetic programs:
#'   - Windows line endings (\r\n → \n)
#'   - Numbered comments in Mesquite TRANSLATE blocks (`[0]`, `[1]`, ...)
#'     that ape cannot parse
#'
#' @param tree_path Path to the original file.
#' @return Path to the clean file (temporary if changes were made, original otherwise).
.sanitize_tree <- function(tree_path) {
  lines <- readLines(tree_path, warn = FALSE)

  # 1) readLines already handles \r\n on most platforms,
  #    but we force explicit cleanup just in case
  lines <- gsub("\r", "", lines)

  # 2) Remove numeric comments at the start of lines: "[0]", "[12]", etc.
  #    Pattern: line starting with optional spaces/tabs followed by [N]
  lines <- gsub("^(\\s*)\\[[0-9]+\\]\\s*", "\\1", lines, perl = TRUE)

  tmp <- tempfile(fileext = tools::file_ext(tree_path))
  writeLines(lines, tmp)
  tmp
}


#' Reads a phylogenetic tree file (Newick or NEXUS)
#'
#' Compatible with files from Mesquite, TNT, WinClada and other programs
#' that produce non-standard variants of the NEXUS/Newick format.
#'
#' @param tree_path    Path to the tree file.
#' @param tree_format  "auto" (default), "newick" or "nexus".
#' @return phylo object.
read_tree <- function(tree_path, tree_format = "auto") {
  if (!file.exists(tree_path))
    stop(paste0("Tree file not found: '", tree_path, "'"))

  fmt <- tolower(trimws(tree_format))

  if (fmt == "auto") {
    ext <- tolower(tools::file_ext(tree_path))
    fmt <- if (ext %in% c("nex", "nexus")) "nexus" else "newick"
  }

  # Sanitize the file before parsing
  clean_path <- .sanitize_tree(tree_path)

  tree <- tryCatch({
    if (fmt == "nexus") read.nexus(clean_path) else read.tree(clean_path)
  }, error = function(e) {
    # If the detected format fails, try the other one
    alt_fmt <- if (fmt == "nexus") "newick" else "nexus"
    message(sprintf(
      "Format '%s' failed (%s). Trying '%s'...", fmt, e$message, alt_fmt
    ))
    tryCatch(
      if (alt_fmt == "nexus") read.nexus(clean_path) else read.tree(clean_path),
      error = function(e2) stop(paste0(
        "Could not read tree '", tree_path, "'.\n",
        "  As ", fmt,     ": ", e$message, "\n",
        "  As ", alt_fmt, ": ", e2$message
      ))
    )
  })

  if (is.null(tree))
    stop(paste0("Could not read tree from '", tree_path, "'."))

  tree
}


# ============================================================================ #
# CSV READING
# ============================================================================ #

#' Reads a character CSV file and prepares it for MultiMapR
#'
#' @param csv_path           Path to the CSV file.
#' @param sep                Field separator (default ",").
#' @param species_col        Index or name of the species column (default 1).
#' @param normalize_spaces   If TRUE, replaces spaces with "_" in species names
#'                           (default FALSE).
#' @return Data.frame with a "Species" column and "charN" or custom-named columns.
read_csv_characters <- function(csv_path,
                                sep               = ",",
                                species_col       = 1,
                                normalize_spaces  = FALSE) {

  if (!file.exists(csv_path))
    stop(paste0("CSV file not found: '", csv_path, "'"))

  # Remove BOM if present (Excel UTF-8 adds it); work on clean file
  clean_path  <- .remove_bom(csv_path)
  has_header  <- .detect_header(clean_path, sep = sep)

  raw_data <- read.csv(clean_path,
                       header           = has_header,
                       sep              = sep,
                       stringsAsFactors = FALSE,
                       colClasses       = "character",
                       check.names      = FALSE)

  if (ncol(raw_data) < 2)
    stop("The CSV must have at least two columns: species + one character.")

  # --- Extract species column ---
  if (is.character(species_col)) {
    if (!species_col %in% colnames(raw_data))
      stop(paste0("Species column '", species_col, "' not found."))
    sp_idx <- which(colnames(raw_data) == species_col)
  } else {
    sp_idx <- as.integer(species_col)
    if (sp_idx < 1 || sp_idx > ncol(raw_data))
      stop(paste0("Species column index out of range: ", sp_idx))
  }

  species <- trimws(raw_data[[sp_idx]])
  if (normalize_spaces)
    species <- gsub(" ", "_", species)

  # --- Character columns (everything except the species column) ---
  char_idx   <- setdiff(seq_len(ncol(raw_data)), sp_idx)
  char_data  <- raw_data[, char_idx, drop = FALSE]

  # Determine character column names
  char_names <- if (has_header) colnames(raw_data)[char_idx] else NULL
  char_names <- .normalize_char_names(char_names, ncol(char_data))

  # Build final data.frame
  result           <- as.data.frame(char_data, stringsAsFactors = FALSE)
  colnames(result) <- char_names
  result           <- data.frame(Species = species, result,
                                 stringsAsFactors = FALSE,
                                 check.names      = FALSE)
  rownames(result) <- NULL

  message(sprintf(
    "CSV loaded: %d species, %d character(s). Header %s.",
    nrow(result),
    length(char_names),
    if (has_header) "detected (custom names)" else "not detected (automatic names)"
  ))

  result
}


# ============================================================================ #
# TREE ↔ CSV COMPATIBILITY VALIDATION
# ============================================================================ #

#' Checks that the CSV names match the tree tip.labels
#'
#' @param tree     phylo object.
#' @param data     Data.frame with a "Species" column.
#' @param strict   If TRUE raises an error on any discrepancy;
#'                 if FALSE (default) only warns and returns unmatched names.
#' @return Invisible: vector of CSV names not found in the tree.
validate_compatibility <- function(tree, data, strict = FALSE) {
  tips            <- tree$tip.label
  sp              <- data$Species
  unmatched_csv   <- setdiff(sp,   tips)
  unmatched_tree  <- setdiff(tips, sp)
  ok              <- length(unmatched_csv) == 0 && length(unmatched_tree) == 0

  if (!ok) {
    msg_parts <- character(0)
    if (length(unmatched_csv) > 0)
      msg_parts <- c(msg_parts,
                     paste0("In CSV but NOT in tree (", length(unmatched_csv), "): ",
                            paste(head(unmatched_csv, 10), collapse = ", "),
                            if (length(unmatched_csv) > 10) " ..." else ""))
    if (length(unmatched_tree) > 0)
      msg_parts <- c(msg_parts,
                     paste0("In tree but NOT in CSV (", length(unmatched_tree), "): ",
                            paste(head(unmatched_tree, 10), collapse = ", "),
                            if (length(unmatched_tree) > 10) " ..." else ""))

    msg <- paste(c("Discrepancies between tree and CSV:", msg_parts), collapse = "\n  ")
    if (strict) stop(msg) else warning(msg)
  } else {
    message("Validation OK: all ", length(tips), " tips match the CSV.")
  }

  invisible(unmatched_csv)
}


# ============================================================================ #
# MAIN FUNCTION
# ============================================================================ #

#' Loads a tree and a character CSV ready to use with MultiMapR
#'
#' Automatically detects:
#'   - Tree format (Newick / NEXUS) by extension.
#'   - Presence or absence of a header in the CSV.
#'   - UTF-8 BOM in the CSV (files exported from Excel).
#'   - Character types (numeric or strings): preserved as text.
#'
#' @param tree_path          Path to the tree file (.tre, .nwk, .nex, .nexus, ...).
#' @param csv_path           Path to the character CSV file.
#' @param sep                CSV field separator (default ",").
#' @param species_col        Species column: integer index or column name (default 1).
#' @param normalize_spaces   If TRUE converts spaces to "_" in species names
#'                           (useful when the tree uses "_" and the CSV uses " ").
#' @param tree_format        "auto" (default), "newick" or "nexus".
#' @param strict             If TRUE, raises an error on tree/CSV discrepancies
#'                           (default FALSE).
#'
#' @return List with:
#'   \item{tree}{phylo object ready for ape / MultiMapR.}
#'   \item{characters}{Data.frame with a "Species" column and character columns.}
#'
#' @examples
#' \dontrun{
#'   # CSV without header (numeric characters)
#'   d1 <- load_data("tree2_jadc.tre", "Matriz_JAIR.csv")
#'   execute_phylogeny(d1$tree, d1$characters)
#'
#'   # CSV with BOM and species names with spaces
#'   d2 <- load_data("bats_tre.tre", "bats_matrix.csv",
#'                   normalize_spaces = TRUE)
#'   execute_phylogeny(d2$tree, d2$characters)
#'
#'   # Semicolon separator, species column by name
#'   d3 <- load_data("tree.nex", "data.csv",
#'                   sep = ";", species_col = "Taxon")
#'   execute_phylogeny(d3$tree, d3$characters)
#'
#'   # Polymorphic use: pass paths directly to the orchestrator
#'   execute_phylogeny("my_tree.tre", "my_characters.csv")
#' }
load_data <- function(tree_path,
                      csv_path,
                      sep               = ",",
                      species_col       = 1,
                      normalize_spaces  = FALSE,
                      tree_format       = "auto",
                      strict            = FALSE) {

  tree       <- read_tree(tree_path, tree_format = tree_format)
  characters <- read_csv_characters(csv_path,
                                    sep              = sep,
                                    species_col      = species_col,
                                    normalize_spaces = normalize_spaces)

  validate_compatibility(tree, characters, strict = strict)

  list(tree = tree, characters = characters)
}