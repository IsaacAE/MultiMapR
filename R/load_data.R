################################################################################
# load_data.R
#
# Utilities to load and prepare phylogenetic data in MultiMapR.
#
# MAIN FUNCTION:
#   load_data(tree_path, csv_path, ...)  ->  list(tree, characters)
#
# SUPPORTED TREE FORMATS:
#   .tre / .tree / .nwk  ->  read.tree()   (Newick)
#   .nex / .nexus        ->  read.nexus()  (NEXUS)
#   TNT `tread`          ->  converted to Newick (with or without commas)
#   (automatic detection by extension; the other format and a tolerant
#    NEXUS TREES-block parser are tried if the first attempt fails)
#
# SUPPORTED CHARACTER MATRIX FORMATS:
#   CSV / TSV / ";"-separated text:
#     - With or without a header for characters (automatic detection)
#     - Species column: any column (index or name)
#     - Numeric and text (string) characters
#     - Inapplicable "-" and unknown "?" are preserved as-is
#     - Trailing separators / empty columns are discarded
#     - UTF-8 BOM (files exported from Excel) removed automatically
#   TNT (`xread`, e.g. .tnt / .ss), including `cnames` character names
#   NEXUS CHARACTERS / DATA block (Mesquite, MorphoBank, PAUP...)
#
# Polymorphisms written as [01], {01} or (0 1) are converted to "0/1",
# the notation understood by the Fitch module.
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


#' Reads a text file as lines, without BOM and without carriage returns
#'
#' @param path Path to the file.
#' @return Character vector with one element per line.
.read_lines_clean <- function(path) {
  lines <- readLines(.remove_bom(path), warn = FALSE)
  gsub("\r", "", lines, fixed = TRUE)
}


#' Removes surrounding single or double quotes (and whitespace) from names
#'
#' @param x Character vector.
#' @return Character vector.
.strip_quotes <- function(x) {
  x <- trimws(x)
  x <- sub("^(['\"])(.*)\\1$", "\\2", x, perl = TRUE)
  trimws(gsub("''", "'", x, fixed = TRUE))
}


# ============================================================================ #
# AUTOMATIC HEADER DETECTION
# ============================================================================ #

# Labels commonly used for the species column
.SPECIES_LABELS <- c("species", "specie", "especie", "especies", "taxon",
                     "taxa", "taxones", "taxonname", "terminal", "terminals",
                     "terminales", "otu", "otus", "name", "names", "nombre",
                     "nombres", "sp", "spp")

#' Decides whether the first row of a raw table is a header
#'
#' Rules, in order:
#'   1. The species cell is a typical label ("Species", "Taxon", "OTU",
#'      "Especie"...)                                       -> header.
#'   2. The character cells are consecutive integers 1..k   -> header.
#'   3. Most character cells of the first row are values that also appear
#'      in the same column further down (e.g. "0", "1", "?"), or all of them
#'      look like state codes                               -> data (no header).
#'   4. Otherwise                                           -> header.
#' With a single row the old rule is used: taxon names contain a space or "_".
#'
#' @param raw     Data.frame read with \code{header = FALSE} (all character).
#' @param sp_idx  Index of the species column.
#' @return TRUE if a header is detected, FALSE otherwise.
.detect_header <- function(raw, sp_idx = 1L) {
  first        <- trimws(as.character(unlist(raw[1, ], use.names = FALSE)))
  species_cell <- .strip_quotes(first[sp_idx])

  if (tolower(gsub("[^A-Za-z]", "", species_cell)) %in% .SPECIES_LABELS)
    return(TRUE)

  char_cols <- setdiff(seq_along(first), sp_idx)
  filled    <- char_cols[nzchar(first[char_cols])]
  if (nrow(raw) < 2 || length(filled) == 0)
    return(!grepl("[ _]", species_cell))

  if (all(grepl("^[0-9]+$", first[filled])) &&
      length(filled) > 2 &&
      identical(as.integer(first[filled]), seq_along(filled)))
    return(TRUE)

  body   <- raw[-1, , drop = FALSE]
  in_col <- vapply(filled, function(j) first[j] %in% trimws(body[[j]]), logical(1))
  codes  <- grepl("^([0-9?-]|[\\[{(][^\\]})]*[\\]})])$", first[filled], perl = TRUE)

  !(mean(in_col) >= 0.5 || all(codes))
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
#'                   NULL or empty -> automatic names are generated.
#' @param n          Number of character columns.
#' @return Character vector of length n.
.normalize_char_names <- function(raw_names, n) {
  if (is.null(raw_names) || length(raw_names) == 0) {
    return(paste0("char", seq_len(n)))
  }
  names   <- make.names(trimws(raw_names), unique = TRUE)
  empty   <- nchar(trimws(raw_names)) == 0 | grepl("^V[0-9]+$", names)
  names[empty] <- paste0("char", which(empty))
  make.unique(names)
}


# ============================================================================ #
# STATE NORMALIZATION
# ============================================================================ #

#' Converts polymorphism notations to the "0/1" form used by Fitch
#'
#'   `[01]`, `{01}`, `(01)`, `[0 1]`, `(0,1)`, `{0/1}`  ->  `0/1`
#'   `[0]`                                              ->  `0`
#' Any other value is returned unchanged (trimmed).
#'
#' @param x Character vector of states.
#' @return Character vector of the same length.
.normalize_polymorphism <- function(x) {
  x    <- trimws(x)
  poly <- grepl("^[\\[{(].*[\\]})]$", x, perl = TRUE)
  if (!any(poly)) return(x)
  x[poly] <- vapply(x[poly], function(s) {
    inner <- substr(s, 2, nchar(s) - 1)
    parts <- if (grepl("[ ,/&|]", inner)) strsplit(inner, "[ ,/&|]+")[[1]]
             else strsplit(inner, "")[[1]]
    parts <- unique(parts[nzchar(parts)])
    if (length(parts) == 0) "?" else paste(parts, collapse = "/")
  }, character(1), USE.NAMES = FALSE)
  x
}


#' Splits a TNT/NEXUS state string into one token per character
#'
#' `01?[01]-{12}(34)` -> `0` `1` `?` `[01]` `-` `{12}` `(34)`
#'
#' @param s Single string without whitespace.
#' @return Character vector of tokens.
.tokenize_states <- function(s) {
  if (!nzchar(s)) return(character(0))
  regmatches(s, gregexpr("\\[[^\\]]*\\]|\\{[^}]*\\}|\\([^)]*\\)|.", s, perl = TRUE))[[1]]
}


#' Builds the final character data.frame shared by all matrix readers
#'
#' @param species           Character vector of species names.
#' @param states            Character matrix/data.frame (rows = species).
#' @param char_names        Raw character names or NULL.
#' @param normalize_spaces  Replace spaces with "_" in species names.
#' @param source_label      Text used in the summary message.
#' @return Data.frame with a "Species" column and one column per character.
.build_character_table <- function(species, states, char_names,
                                   normalize_spaces, source_label) {
  species <- .strip_quotes(species)
  if (normalize_spaces) species <- gsub("\\s+", "_", species)

  states <- as.data.frame(states, stringsAsFactors = FALSE)
  states[] <- lapply(states, function(col) .normalize_polymorphism(as.character(col)))

  result           <- states
  colnames(result) <- .normalize_char_names(char_names, ncol(states))
  result           <- data.frame(Species = species, result,
                                 stringsAsFactors = FALSE,
                                 check.names      = FALSE)
  rownames(result) <- NULL

  dups <- unique(species[duplicated(species)])
  if (length(dups) > 0)
    warning("Duplicated species in the matrix: ", paste(head(dups, 10), collapse = ", "))

  message(sprintf("Matrix loaded (%s): %d species, %d character(s).",
                  source_label, nrow(result), ncol(result) - 1))
  result
}


# ============================================================================ #
# TREE READING
# ============================================================================ #

#' Converts the text of a TNT `tread` file to a Newick string
#'
#' Handles the quoted title after `tread`, the `proc` footer, several trees
#' separated by "*" (only the first is kept) and the space-separated,
#' comma-less TNT syntax, e.g. \code{(A (B C)(D E))}.
#'
#' @param lines Lines of the file.
#' @return Newick string ending with ";".
.tnt_tread_to_newick <- function(lines) {
  txt <- paste(lines, collapse = " ")
  txt <- sub("^\\s*tread\\b", "", txt, ignore.case = TRUE, perl = TRUE)
  txt <- gsub("'[^']*'", " ", txt)                          # quoted titles
  txt <- sub("\\bproc\\b.*$", "", txt, ignore.case = TRUE, perl = TRUE)
  txt <- gsub("\\[[^\\]]*\\]", " ", txt, perl = TRUE)       # comments

  trees <- trimws(strsplit(txt, "[*;]")[[1]])
  trees <- trees[grepl("\\(", trees)]
  if (length(trees) == 0) stop("No tree found in TNT 'tread' file.")
  if (length(trees) > 1)
    message(sprintf("TNT file contains %d trees; using the first one.", length(trees)))

  nwk <- trees[1]
  tok <- "[^\\s(),;:]"
  nwk <- gsub("\\(\\s+", "(", nwk, perl = TRUE)
  nwk <- gsub("\\s+\\)", ")", nwk, perl = TRUE)
  nwk <- gsub("\\)\\s*\\(", "),(", nwk, perl = TRUE)
  nwk <- gsub(sprintf("\\)\\s*(?=%s)", tok), "),", nwk, perl = TRUE)
  nwk <- gsub(sprintf("(%s)\\s+\\(", tok), "\\1,(", nwk, perl = TRUE)
  nwk <- gsub(sprintf("(%s)\\s+(?=%s)", tok, tok), "\\1,", nwk, perl = TRUE)
  nwk <- gsub(",+", ",", nwk)
  paste0(nwk, ";")
}


#' Sanitizes a tree file for compatibility with ape
#'
#' Removes common problems produced by different phylogenetic programs:
#'   - Windows line endings (CR+LF converted to LF) and UTF-8 BOM
#'   - Numbered comments in Mesquite TRANSLATE blocks (`[0]`, `[1]`, ...)
#'     that ape cannot parse
#'   - TNT `tread` format (with or without a quoted title): rewritten as plain
#'     Newick so that read.tree() can parse it. TNT separates tips with spaces
#'     instead of commas; the missing commas are inserted.
#'
#' @param tree_path Path to the original file.
#' @return Path to the clean file (temporary).
.sanitize_tree <- function(tree_path) {
  lines <- .read_lines_clean(tree_path)

  # Remove numeric comments at the start of lines: "[0]", "[12]", etc.
  lines <- gsub("^(\\s*)\\[[0-9]+\\]\\s*", "\\1", lines, perl = TRUE)

  # TNT tread format: first non-empty line starts with "tread"
  first_content <- which(nchar(trimws(lines)) > 0)[1]
  is_tnt <- !is.na(first_content) &&
    grepl("^\\s*tread\\b", lines[first_content], ignore.case = TRUE, perl = TRUE)

  if (is_tnt) {
    message("TNT 'tread' format detected. Converting to Newick...")
    lines <- .tnt_tread_to_newick(lines)
  }

  tmp <- tempfile(fileext = paste0(".", tools::file_ext(tree_path)))
  writeLines(lines, tmp)
  tmp
}


#' Tolerant parser for the TREES block of a NEXUS file
#'
#' Used when \code{ape::read.nexus()} fails, e.g. with TNT exports where the
#' Newick string is on a different line from \code{tree name =}, with
#' \code{[&U]} comments or with \code{end ;}. Supports TRANSLATE tables.
#'
#' @param path Path to the (sanitized) NEXUS file.
#' @return phylo object (the first tree of the block).
.read_nexus_trees_fallback <- function(path) {
  txt <- paste(.read_lines_clean(path), collapse = "\n")
  m   <- regexpr("(?is)begin\\s+trees\\s*;.*?(\\bend(block)?\\s*;|$)", txt, perl = TRUE)
  if (m == -1) stop("No TREES block found.")

  block <- regmatches(txt, m)
  block <- sub("(?is)^begin\\s+trees\\s*;", "", block, perl = TRUE)
  block <- sub("(?is)\\bend(block)?\\s*;$", "", block, perl = TRUE)
  block <- gsub("\\[[^\\]]*\\]", "", block, perl = TRUE)

  stmts <- trimws(strsplit(block, ";", fixed = TRUE)[[1]])
  trees <- grep("^u?tree\\b", stmts, ignore.case = TRUE, perl = TRUE, value = TRUE)
  if (length(trees) == 0) stop("TREES block has no 'tree' statement.")
  if (length(trees) > 1)
    message(sprintf("NEXUS file contains %d trees; using the first one.", length(trees)))

  newick <- sub("^[^=]*=\\s*", "", trees[1])
  tree   <- read.tree(text = paste0(gsub("\\s*\n\\s*", "", newick), ";"))

  translate <- grep("^translate\\b", stmts, ignore.case = TRUE, perl = TRUE, value = TRUE)
  if (length(translate) > 0) {
    pairs <- trimws(strsplit(sub("(?i)^translate\\s*", "", translate[1], perl = TRUE),
                             ",")[[1]])
    pairs <- pairs[nzchar(pairs)]
    keys  <- sub("\\s.*$", "", pairs)
    vals  <- .strip_quotes(sub("^\\S+\\s+", "", pairs))
    idx   <- match(tree$tip.label, keys)
    tree$tip.label[!is.na(idx)] <- vals[idx[!is.na(idx)]]
  }
  tree
}


#' Reads a phylogenetic tree file (Newick, NEXUS or TNT)
#'
#' Compatible with files from Mesquite, TNT, WinClada and other programs
#' that produce non-standard variants of the NEXUS/Newick format.
#' TNT files using the \code{tread} / \code{proc-;} syntax are automatically
#' converted to plain Newick before parsing. Tip labels are trimmed and
#' surrounding quotes are removed. If the file holds several trees, the first
#' one is used.
#'
#' @param tree_path    Path to the tree file.
#' @param tree_format  "auto" (default), "newick" or "nexus".
#' @return phylo object.
read_tree <- function(tree_path, tree_format = "auto") {
  if (!file.exists(tree_path))
    stop(paste0("Tree file not found: '", tree_path, "'"))

  fmt <- tolower(trimws(tree_format))

  # Sanitize the file before parsing
  clean_path <- .sanitize_tree(tree_path)

  if (fmt == "auto") {
    ext       <- tolower(tools::file_ext(tree_path))
    is_nexus  <- grepl("^\\s*#nexus", readLines(clean_path, n = 1, warn = FALSE),
                       ignore.case = TRUE)
    fmt <- if (ext %in% c("nex", "nexus", "nxs") || isTRUE(is_nexus)) "nexus" else "newick"
  }

  readers <- list(
    newick   = function(p) read.tree(p),
    nexus    = function(p) read.nexus(p),
    fallback = function(p) .read_nexus_trees_fallback(p)
  )
  order  <- c(fmt, setdiff(c("newick", "nexus"), fmt), "fallback")
  errors <- character(0)
  tree   <- NULL

  for (r in order) {
    # Warnings are only relevant if this reader succeeds
    warns <- character(0)
    tree <- tryCatch(
      withCallingHandlers(readers[[r]](clean_path), warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      }),
      error = function(e) {
        errors[[r]] <<- conditionMessage(e)
        NULL
      })
    if (inherits(tree, c("phylo", "multiPhylo")))
      for (w in unique(warns)) warning(w, call. = FALSE)
    if (inherits(tree, "multiPhylo")) {
      message(sprintf("File contains %d trees; using the first one.", length(tree)))
      tree <- tree[[1]]
    }
    if (inherits(tree, "phylo")) break
    tree <- NULL
    if (r != "fallback")
      message(sprintf("Format '%s' failed (%s). Trying next reader...",
                      r, if (r %in% names(errors)) trimws(errors[[r]]) else "no tree"))
  }

  if (is.null(tree))
    stop(paste0("Could not read tree '", tree_path, "'.\n",
                paste0("  As ", names(errors), ": ", trimws(errors), collapse = "\n")))

  tree$tip.label <- .strip_quotes(tree$tip.label)
  tree
}


# ============================================================================ #
# CHARACTER MATRIX READING
# ============================================================================ #

#' Guesses the field separator of a delimited text file
#'
#' Chooses among ",", ";" and tab the one that splits the first lines into
#' the largest consistent number of fields.
#'
#' @param path Path to the file (BOM-free).
#' @return Separator character.
.detect_sep <- function(path) {
  lines <- head(readLines(path, warn = FALSE), 20)
  lines <- lines[nzchar(trimws(lines))]
  cands <- c(",", ";", "\t")
  score <- vapply(cands, function(s) {
    n <- lengths(regmatches(lines, gregexpr(s, lines, fixed = TRUE)))
    if (min(n) == 0) 0 else stats::median(n)
  }, numeric(1))
  if (all(score == 0)) "," else cands[which.max(score)]
}


#' Reads a character CSV file and prepares it for MultiMapR
#'
#' Rows may have different numbers of fields (trailing separators are common
#' in files exported from Excel); columns that are completely empty are
#' discarded. The header is detected automatically (see \code{.detect_header})
#' unless \code{header} is given.
#'
#' @param csv_path           Path to the CSV file.
#' @param sep                Field separator (default ","). Use \code{"auto"}
#'                           to detect ",", ";" or tab. If the given separator
#'                           yields a single column, the others are tried.
#' @param species_col        Index or name of the species column (default 1).
#' @param normalize_spaces   If TRUE, replaces spaces with "_" in species names
#'                           (default FALSE).
#' @param header             \code{NA} (default, automatic), \code{TRUE} or
#'                           \code{FALSE}.
#' @return Data.frame with a "Species" column and "charN" or custom-named columns.
read_csv_characters <- function(csv_path,
                                sep               = ",",
                                species_col       = 1,
                                normalize_spaces  = FALSE,
                                header            = NA) {

  if (!file.exists(csv_path))
    stop(paste0("CSV file not found: '", csv_path, "'"))

  # Remove BOM if present (Excel UTF-8 adds it); work on clean file
  clean_path <- .remove_bom(csv_path)

  read_raw <- function(s) {
    n_fields <- utils::count.fields(clean_path, sep = s, quote = "\"",
                                    blank.lines.skip = TRUE, comment.char = "")
    n_cols   <- max(n_fields, 1L, na.rm = TRUE)
    utils::read.table(clean_path, header = FALSE, sep = s, quote = "\"",
                      fill = TRUE, colClasses = "character",
                      col.names = paste0("V", seq_len(n_cols)),
                      na.strings = character(0), comment.char = "",
                      strip.white = TRUE, blank.lines.skip = TRUE,
                      stringsAsFactors = FALSE, check.names = FALSE)
  }

  if (identical(sep, "auto")) sep <- .detect_sep(clean_path)
  raw_data <- read_raw(sep)
  if (ncol(raw_data) < 2) {
    alt <- .detect_sep(clean_path)
    if (alt != sep) {
      message(sprintf("Separator '%s' gives a single column; using '%s' instead.",
                      sep, if (alt == "\t") "\\t" else alt))
      raw_data <- read_raw(alt)
    }
  }

  # Drop fully empty rows
  raw_data <- raw_data[rowSums(raw_data != "") > 0, , drop = FALSE]
  if (ncol(raw_data) < 2 || nrow(raw_data) < 1)
    stop("The CSV must have at least two columns: species + one character.")

  # --- Species column index (names only make sense with a header) ---
  first_row <- trimws(as.character(unlist(raw_data[1, ], use.names = FALSE)))
  if (is.character(species_col)) {
    sp_idx <- match(species_col, first_row)
    if (is.na(sp_idx))
      stop(paste0("Species column '", species_col, "' not found."))
    has_header <- TRUE
  } else {
    sp_idx <- as.integer(species_col)
    if (sp_idx < 1 || sp_idx > ncol(raw_data))
      stop(paste0("Species column index out of range: ", sp_idx))
    has_header <- if (is.na(header)) .detect_header(raw_data, sp_idx) else isTRUE(header)
  }

  header_row <- if (has_header) first_row else NULL
  body       <- if (has_header) raw_data[-1, , drop = FALSE] else raw_data

  # --- Character columns: all except species, discarding empty ones ---
  char_idx <- setdiff(seq_len(ncol(raw_data)), sp_idx)
  empty    <- vapply(char_idx, function(j)
    all(body[[j]] == "") && (is.null(header_row) || header_row[j] == ""), logical(1))
  if (any(empty))
    message(sprintf("%d empty column(s) discarded (trailing separators).", sum(empty)))
  char_idx <- char_idx[!empty]

  body[char_idx] <- lapply(body[char_idx], function(col) ifelse(col == "", "?", col))

  .build_character_table(
    species          = body[[sp_idx]],
    states           = body[, char_idx, drop = FALSE],
    char_names       = if (has_header) header_row[char_idx] else NULL,
    normalize_spaces = normalize_spaces,
    source_label     = sprintf("CSV, header %s",
                               if (has_header) "detected" else "not detected")
  )
}


#' Reads a TNT character matrix (`xread` command)
#'
#' Supports an optional quoted title after \code{xread}, polymorphisms in
#' square brackets, interleaved blocks and the \code{cnames} block for
#' character names. State names from \code{cnames} are stored in the
#' \code{"state_labels"} attribute of the result.
#'
#' @param path              Path to the TNT file.
#' @param normalize_spaces  Replace spaces with "_" in species names.
#' @return Data.frame with a "Species" column and one column per character.
read_tnt_characters <- function(path, normalize_spaces = FALSE) {
  lines <- .read_lines_clean(path)
  start <- grep("^\\s*xread\\b", lines, ignore.case = TRUE, perl = TRUE)[1]
  if (is.na(start)) stop("No 'xread' command found in '", path, "'.")

  body <- paste(lines[start:length(lines)], collapse = "\n")
  body <- sub("^\\s*xread", "", body, ignore.case = TRUE, perl = TRUE)
  body <- sub("^\\s*'[^']*'", "", body, perl = TRUE)          # optional title

  dims <- regmatches(body, regexec("^\\s*([0-9]+)\\s+([0-9]+)", body))[[1]]
  if (length(dims) < 3) stop("Could not read 'nchar ntax' after xread.")
  n_char <- as.integer(dims[2])
  n_tax  <- as.integer(dims[3])
  body   <- substring(body, nchar(dims[1]) + 1)

  end      <- regexpr(";", body, fixed = TRUE)
  mat_txt  <- if (end > 0) substr(body, 1, end - 1) else body
  rest_txt <- if (end > 0) substring(body, end + 1) else ""

  mlines <- trimws(strsplit(mat_txt, "\n", fixed = TRUE)[[1]])
  mlines <- mlines[nzchar(mlines) & !grepl("^&", mlines)]

  seqs <- list()
  for (ln in mlines) {
    name <- sub("\\s.*$", "", ln)
    seqs[[name]] <- paste0(seqs[[name]], gsub("\\s+", "", substring(ln, nchar(name) + 1)))
  }

  # Character names (and state names) from cnames: "{0 name state0 state1 ;"
  char_names   <- NULL
  state_labels <- NULL
  cn <- regexpr("(?is)\\bcnames\\b(.*?)\\n\\s*;", rest_txt, perl = TRUE)
  if (cn > 0) {
    entries <- regmatches(regmatches(rest_txt, cn),
                          gregexpr("\\{[^;]*;", regmatches(rest_txt, cn)))[[1]]
    idx   <- as.integer(sub("^\\{\\s*([0-9]+).*$", "\\1", entries))
    parts <- strsplit(trimws(sub(";$", "", sub("^\\{\\s*[0-9]+", "", entries))), "\\s+")
    char_names   <- rep("", n_char)
    state_labels <- vector("list", n_char)
    ok <- !is.na(idx) & idx < n_char
    char_names[idx[ok] + 1]   <- vapply(parts[ok], `[`, character(1), 1)
    state_labels[idx[ok] + 1] <- lapply(parts[ok], `[`, -1)
  }

  .matrix_from_sequences(seqs, n_char, n_tax, char_names, state_labels,
                         normalize_spaces, "TNT")
}


#' Reads the CHARACTERS / DATA block of a NEXUS file
#'
#' Handles quoted taxon names, \code{MISSING}, \code{GAP} and
#' \code{MATCHCHAR} symbols, interleaved and wrapped matrices, polymorphisms
#' \code{(01)} / \code{{01}}, and character names from \code{CHARLABELS} or
#' \code{CHARSTATELABELS} (state names are stored in the
#' \code{"state_labels"} attribute of the result).
#'
#' @param path              Path to the NEXUS file.
#' @param normalize_spaces  Replace spaces with "_" in species names.
#' @return Data.frame with a "Species" column and one column per character.
read_nexus_characters <- function(path, normalize_spaces = FALSE) {
  txt <- paste(.read_lines_clean(path), collapse = "\n")
  m   <- regexpr("(?is)begin\\s+(characters|data)\\s*;.*?\\bend(block)?\\s*;", txt, perl = TRUE)
  if (m == -1) stop("No CHARACTERS/DATA block found in '", path, "'.")
  block <- gsub("\\[[^\\]]*\\]", "", regmatches(txt, m), perl = TRUE)

  opt <- function(key, default = NULL) {
    r <- regmatches(block, regexec(sprintf("(?i)\\b%s\\s*=\\s*\"?([^\\s;\"]+)", key),
                                   block, perl = TRUE))[[1]]
    if (length(r) < 2) default else r[2]
  }
  n_char     <- as.integer(opt("nchar", NA))
  n_tax      <- as.integer(opt("ntax", NA))
  missing    <- substr(opt("missing", "?"), 1, 1)
  gap        <- substr(opt("gap", "-"), 1, 1)
  matchchar  <- opt("matchchar")
  interleave <- grepl("(?i)\\binterleave\\b(\\s*=\\s*yes)?", block, perl = TRUE) &&
                !grepl("(?i)\\binterleave\\s*=\\s*no", block, perl = TRUE)
  if (is.na(n_char)) stop("NCHAR not declared in the CHARACTERS block.")

  quoted_tokens <- function(s) {
    tk <- regmatches(s, gregexpr("'(?:[^']|'')*'|[^\\s,]+", s, perl = TRUE))[[1]]
    .strip_quotes(tk)
  }

  # Text of a command ("CHARLABELS ... ;") and its comma-separated entries
  command_text <- function(cmd) {
    r <- regmatches(block, regexec(sprintf("(?is)\\b%s\\b((?:'(?:[^']|'')*'|[^;'])*);", cmd),
                                   block, perl = TRUE))[[1]]
    if (length(r) == 2) r[2] else NULL
  }
  command_entries <- function(s) {
    e <- trimws(regmatches(s, gregexpr("(?:'(?:[^']|'')*'|[^,'])+", s, perl = TRUE))[[1]])
    e[nzchar(e)]
  }

  # --- Character and state names ---
  char_names   <- NULL
  state_labels <- vector("list", n_char)
  cl <- command_text("charlabels")
  if (!is.null(cl)) char_names <- quoted_tokens(cl)[seq_len(n_char)]

  csl <- command_text("charstatelabels")
  if (!is.null(csl)) {
    char_names <- rep("", n_char)
    for (e in command_entries(csl)) {
      halves <- strsplit(e, "/", fixed = TRUE)[[1]]
      tk     <- quoted_tokens(halves[1])
      i      <- suppressWarnings(as.integer(tk[1]))
      if (is.na(i) || i < 1 || i > n_char) next
      char_names[i] <- if (length(tk) > 1) tk[2] else ""
      if (length(halves) > 1) state_labels[[i]] <- quoted_tokens(paste(halves[-1], collapse = "/"))
    }
  }

  sl <- command_text("statelabels")
  if (!is.null(sl)) {
    for (e in command_entries(sl)) {
      tk <- quoted_tokens(e)
      i  <- suppressWarnings(as.integer(tk[1]))
      if (!is.na(i) && i >= 1 && i <= n_char) state_labels[[i]] <- tk[-1]
    }
  }

  # --- Matrix ---
  mx <- regmatches(block, regexec("(?is)\\bmatrix\\b(.*?);", block, perl = TRUE))[[1]]
  if (length(mx) < 2) stop("No MATRIX command found in the CHARACTERS block.")
  mlines <- trimws(strsplit(mx[2], "\n", fixed = TRUE)[[1]])
  mlines <- mlines[nzchar(mlines)]

  seqs    <- list()
  current <- NULL
  for (ln in mlines) {
    # Wrapped (non-interleaved) rows: keep filling the current taxon
    if (!interleave && !is.null(current) &&
        length(.tokenize_states(seqs[[current]])) < n_char) {
      seqs[[current]] <- paste0(seqs[[current]], gsub("\\s+", "", ln))
      next
    }
    name_tok <- regmatches(ln, regexpr("^('(?:[^']|'')*'|\\S+)", ln, perl = TRUE))
    name     <- .strip_quotes(name_tok)
    states   <- gsub("\\s+", "", substring(ln, nchar(name_tok) + 1))
    seqs[[name]] <- paste0(seqs[[name]], states)
    current <- name
  }

  first_seq <- if (!is.null(matchchar)) .tokenize_states(seqs[[1]]) else NULL
  seqs <- lapply(seqs, function(s) {
    tk <- .tokenize_states(s)
    tk[tk == missing] <- "?"
    tk[tk == gap]     <- "-"
    if (!is.null(matchchar)) {
      same <- which(tk == matchchar & seq_along(tk) <= length(first_seq))
      tk[same] <- first_seq[same]
    }
    tk
  })

  .matrix_from_sequences(seqs, n_char, n_tax, char_names, state_labels,
                         normalize_spaces, "NEXUS")
}


#' Turns per-taxon state sequences into the MultiMapR character table
#'
#' @param seqs         Named list: one state string (or token vector) per taxon.
#' @param n_char       Declared number of characters.
#' @param n_tax        Declared number of taxa (NA if unknown).
#' @param char_names   Character names or NULL.
#' @param state_labels List of state names per character or NULL.
#' @param normalize_spaces,source_label See \code{.build_character_table}.
#' @return Data.frame (see \code{.build_character_table}).
.matrix_from_sequences <- function(seqs, n_char, n_tax, char_names, state_labels,
                                   normalize_spaces, source_label) {
  tokens <- lapply(seqs, function(s) if (length(s) == 1) .tokenize_states(s) else s)
  lens   <- lengths(tokens)
  bad    <- lens != n_char
  if (any(bad))
    warning(sprintf("%d taxa do not have %d states (%s); padded with '?' or truncated.",
                    sum(bad), n_char,
                    paste(head(sprintf("%s = %d", names(tokens)[bad], lens[bad]), 5),
                          collapse = ", ")))
  tokens <- lapply(tokens, function(t) c(t, rep("?", n_char))[seq_len(n_char)])

  if (!is.na(n_tax) && n_tax != length(tokens))
    warning(sprintf("Declared %d taxa but found %d in the matrix.", n_tax, length(tokens)))

  states <- do.call(rbind, tokens)
  if (!is.null(char_names) && all(!nzchar(char_names))) char_names <- NULL

  result <- .build_character_table(names(tokens), states, char_names,
                                   normalize_spaces, source_label)
  if (!is.null(state_labels) && any(lengths(state_labels) > 0)) {
    names(state_labels) <- colnames(result)[-1]
    attr(result, "state_labels") <- state_labels
  }
  result
}


#' Detects the format of a character matrix file
#'
#' @param path Path to the file.
#' @return "tnt", "nexus" or "csv".
.detect_matrix_format <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("tnt", "ss", "xread")) return("tnt")
  if (ext %in% c("nex", "nexus", "nxs")) return("nexus")
  head_lines <- head(.read_lines_clean(path), 50)
  if (length(head_lines) && grepl("^\\s*#nexus", head_lines[1], ignore.case = TRUE)) return("nexus")
  if (any(grepl("^\\s*xread\\b", head_lines, ignore.case = TRUE, perl = TRUE))) return("tnt")
  "csv"
}


#' Reads a character matrix in CSV, TNT or NEXUS format
#'
#' @param path              Path to the matrix file.
#' @param format            "auto" (default), "csv", "tnt" or "nexus".
#' @param sep,species_col,header  Only used for CSV files
#'                          (see \code{read_csv_characters}).
#' @param normalize_spaces  Replace spaces with "_" in species names.
#' @return Data.frame with a "Species" column and one column per character.
read_characters <- function(path,
                            format           = "auto",
                            sep              = ",",
                            species_col      = 1,
                            normalize_spaces = FALSE,
                            header           = NA) {
  if (!file.exists(path))
    stop(paste0("Character matrix file not found: '", path, "'"))
  fmt <- tolower(trimws(format))
  if (fmt == "auto") fmt <- .detect_matrix_format(path)
  switch(fmt,
         tnt   = read_tnt_characters(path, normalize_spaces = normalize_spaces),
         nexus = read_nexus_characters(path, normalize_spaces = normalize_spaces),
         csv   = read_csv_characters(path, sep = sep, species_col = species_col,
                                     normalize_spaces = normalize_spaces,
                                     header = header),
         stop("Unknown matrix format: '", format, "'. Use 'csv', 'tnt' or 'nexus'."))
}


# ============================================================================ #
# TREE <-> CSV COMPATIBILITY VALIDATION
# ============================================================================ #

#' Reconciles species names that differ only in spaces, underscores, quotes
#' or letter case
#'
#' "Saccopteryx bilineata" (matrix) is renamed to "Saccopteryx_bilineata"
#' (tree). Only unambiguous one-to-one matches are applied.
#'
#' @param tree  phylo object.
#' @param data  Data.frame with a "Species" column.
#' @return \code{data} with harmonized species names.
.harmonize_species <- function(tree, data) {
  tips      <- tree$tip.label
  sp        <- data$Species
  unmatched <- which(!sp %in% tips)
  if (length(unmatched) == 0) return(data)

  key <- function(x) tolower(gsub("[\\s_]+", "_", .strip_quotes(x), perl = TRUE))
  tip_keys <- key(tips)
  tip_keys[tip_keys %in% tip_keys[duplicated(tip_keys)]] <- NA

  idx <- match(key(sp[unmatched]), tip_keys)
  ok  <- !is.na(idx) & !tips[idx] %in% sp & !duplicated(idx)
  if (any(ok)) {
    sp[unmatched[ok]] <- tips[idx[ok]]
    data$Species <- sp
    message(sprintf("%d species name(s) matched to the tree after normalizing spaces/underscores/case.",
                    sum(ok)))
  }
  data
}


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

#' Loads a tree and a character matrix ready to use with MultiMapR
#'
#' Automatically detects:
#'   - Tree format (Newick / NEXUS / TNT \code{tread}).
#'   - Matrix format (CSV / TNT \code{xread} / NEXUS CHARACTERS block).
#'   - Presence or absence of a header in the CSV, and its separator
#'     problems (trailing separators, ";" or tab instead of ",").
#'   - UTF-8 BOM (files exported from Excel).
#'   - Character types (numeric or strings): preserved as text.
#'   - Species names that differ from the tips only in spaces vs "_",
#'     quotes or letter case (renamed to the tree labels).
#'
#' @param tree_path          Path to the tree file (.tre, .nwk, .nex, .nexus, ...).
#' @param csv_path           Path to the character matrix (.csv, .tnt, .nex, ...).
#' @param sep                CSV field separator (default ","; \code{"auto"}
#'                           detects it).
#' @param species_col        Species column: integer index or column name (default 1).
#' @param normalize_spaces   If TRUE converts spaces to "_" in species names
#'                           (useful when the tree uses "_" and the CSV uses " ").
#' @param tree_format        "auto" (default), "newick" or "nexus".
#' @param strict             If TRUE, raises an error on tree/CSV discrepancies
#'                           (default FALSE).
#' @param matrix_format      "auto" (default), "csv", "tnt" or "nexus".
#' @param header             CSV header: \code{NA} (default, automatic),
#'                           \code{TRUE} or \code{FALSE}.
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
#'   # TNT matrix and TNT tree
#'   d4 <- load_data("tree_tnt.tre", "matrix.tnt")
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
                      strict            = FALSE,
                      matrix_format     = "auto",
                      header            = NA) {

  tree       <- read_tree(tree_path, tree_format = tree_format)
  characters <- read_characters(csv_path,
                                format           = matrix_format,
                                sep              = sep,
                                species_col      = species_col,
                                normalize_spaces = normalize_spaces,
                                header           = header)

  characters <- .harmonize_species(tree, characters)
  validate_compatibility(tree, characters, strict = strict)

  list(tree = tree, characters = characters)
}
