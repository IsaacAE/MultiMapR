################################################################################
# character_stats.R
#
# Per-character summary statistics: missing ("?") and inapplicable ("-") data
# counts. Independent of tree topology or the chosen ancestral algorithm --
# purely a property of the character matrix. Includes plotting companions
# (plot_character_stats / plot_character_completeness) for when the printed
# table is too dense to read at a glance.
################################################################################

#' Default colorblind-safe palette (Okabe-Ito) for the three data categories
#' used throughout this file: scored, missing ("?") and inapplicable ("-").
#' @keywords internal
CHAR_STATS_COLORS <- c(scored = "#009E73", missing = "#E69F00", inapplicable = "#56B4E9")


#' Classifies each cell of a character column as scored/missing/inapplicable
#'
#' @param vals  Character (or coercible) vector: one character column.
#' @return Character vector, same length as \code{vals}, with values
#'   \code{"scored"}, \code{"missing"} or \code{"inapplicable"}.
#' @keywords internal
.classify_character_values <- function(vals) {
  vals      <- as.character(vals)
  vals_trim <- trimws(vals)
  is_missing <- is.na(vals) | vals_trim %in% c("", "?")
  is_inapp   <- !is_missing & vals_trim == "-"
  ifelse(is_missing, "missing", ifelse(is_inapp, "inapplicable", "scored"))
}


#' Missing- and inapplicable-data statistics per character
#'
#' For each character column, counts how many taxa are scored \code{"?"}
#' (missing data, including blank/\code{NA} cells) and how many are scored
#' \code{"-"} (inapplicable), alongside how many taxa carry an actual
#' observed state and how many distinct states are observed. Polymorphic
#' codings (e.g. \code{"0/1"}, \code{"01"}) count as one scored taxon and
#' contribute all of their member states to the observed-state count.
#'
#' @param character_data  Data.frame with a \code{"Species"} column and one
#'                         column per character, as produced by
#'                         \code{load_data()} or passed directly to
#'                         \code{execute_phylogeny()}.
#' @param characters       Optional character vector of column names to
#'                          summarize. Defaults to every column except
#'                          \code{"Species"}.
#' @return A data.frame with one row per character and columns:
#'   \code{character}, \code{n_taxa}, \code{n_missing}, \code{pct_missing},
#'   \code{n_inapplicable}, \code{pct_inapplicable}, \code{n_scored},
#'   \code{n_states}. An extra \code{"TOTAL"} row aggregates across all
#'   summarized characters (percentages over the pooled cell count).
#' @export
character_stats <- function(character_data, characters = NULL) {
  if (!is.data.frame(character_data))
    stop("'character_data' must be a data.frame.")

  if (is.null(characters)) characters <- setdiff(colnames(character_data), "Species")
  missing_ch <- setdiff(characters, colnames(character_data))
  if (length(missing_ch) > 0)
    stop("Character(s) not found in 'character_data': ", paste(missing_ch, collapse = ", "))
  if (length(characters) == 0)
    stop("No characters to summarize.")

  rows <- lapply(characters, function(ch) {
    vals      <- as.character(character_data[[ch]])
    vals_trim <- trimws(vals)
    n_taxa    <- length(vals)

    cls        <- .classify_character_values(vals)
    n_missing  <- sum(cls == "missing")
    n_inapp    <- sum(cls == "inapplicable")
    n_scored   <- sum(cls == "scored")

    scored_vals <- vals_trim[cls == "scored"]
    poly_split <- unlist(lapply(scored_vals, function(v) {
      if (grepl("[,/&|]", v)) trimws(strsplit(v, "[,/&|]")[[1]])
      else if (grepl("^[0-9]+$", v) && nchar(v) > 1) unique(strsplit(v, "")[[1]])
      else v
    }))
    n_states <- length(unique(poly_split))

    data.frame(
      character        = ch,
      n_taxa           = n_taxa,
      n_missing        = n_missing,
      pct_missing      = round(100 * n_missing / n_taxa, 1),
      n_inapplicable   = n_inapp,
      pct_inapplicable = round(100 * n_inapp / n_taxa, 1),
      n_scored         = n_scored,
      n_states         = n_states,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)

  total_taxa  <- sum(out$n_taxa)
  total_row <- data.frame(
    character        = "TOTAL",
    n_taxa           = total_taxa,
    n_missing        = sum(out$n_missing),
    pct_missing      = round(100 * sum(out$n_missing) / total_taxa, 1),
    n_inapplicable   = sum(out$n_inapplicable),
    pct_inapplicable = round(100 * sum(out$n_inapplicable) / total_taxa, 1),
    n_scored         = sum(out$n_scored),
    n_states         = NA_integer_,
    stringsAsFactors = FALSE
  )

  out <- rbind(out, total_row)
  rownames(out) <- NULL
  out
}


#' Prints a character-statistics table to the console
#'
#' Console-formatting wrapper around \code{\link{character_stats}}, used by
#' the interactive CLI (\code{setup_mapping_config()}) to surface missing /
#' inapplicable data counts before the user picks characters to map.
#'
#' @param character_data  See \code{\link{character_stats}}.
#' @param characters      See \code{\link{character_stats}}.
#' @return Invisibly, the data.frame from \code{\link{character_stats}}.
#' @keywords internal
print_character_stats <- function(character_data, characters = NULL) {
  st <- character_stats(character_data, characters)
  body <- st[st$character != "TOTAL", , drop = FALSE]
  total <- st[st$character == "TOTAL", , drop = FALSE]

  cat("\n=== Character statistics (missing / inapplicable data) ===\n")
  name_w <- max(8L, nchar(body$character))
  for (i in seq_len(nrow(body))) {
    cat(sprintf(
      "  %-*s  states=%-2d  scored=%-3d  missing(?)=%-3d (%5.1f%%)  inapplicable(-)=%-3d (%5.1f%%)\n",
      name_w, body$character[i], body$n_states[i], body$n_scored[i],
      body$n_missing[i], body$pct_missing[i],
      body$n_inapplicable[i], body$pct_inapplicable[i]))
  }
  cat(sprintf(
    "  %-*s              scored=%-3d  missing(?)=%-3d (%5.1f%%)  inapplicable(-)=%-3d (%5.1f%%)\n",
    name_w, "TOTAL", total$n_scored, total$n_missing, total$pct_missing,
    total$n_inapplicable, total$pct_inapplicable))

  invisible(st)
}


# ==============================================================================
# PLOTS
# ==============================================================================

#' Bar chart of missing / inapplicable data per character
#'
#' Horizontal stacked bar chart, one bar per character, with segments for
#' scored / missing ("?") / inapplicable ("-") taxon counts. Visual
#' companion to \code{\link{character_stats}} for matrices where the printed
#' table is too dense to scan quickly.
#'
#' @param character_data   See \code{\link{character_stats}}.
#' @param characters       See \code{\link{character_stats}}.
#' @param sort_by          One of \code{"none"} (matrix column order),
#'                         \code{"pct_missing"}, \code{"pct_inapplicable"} or
#'                         \code{"character"}. The two percentage options
#'                         sort descending, so the most incomplete
#'                         characters appear at the top of the plot.
#' @param colors           Named vector with entries \code{scored},
#'                         \code{missing}, \code{inapplicable}. Defaults to
#'                         an Okabe-Ito colorblind-safe triplet.
#' @param export_filename  Optional path (without extension) to also save
#'                         the plot as PNG/PDF. \code{NULL} (default) draws
#'                         on the active graphics device only.
#' @param export_format    \code{"png"} (default) or \code{"pdf"}.
#' @param width,height     Device size in inches when exporting. Defaults
#'                         scale with the number of characters.
#' @return Invisibly, the \code{character_stats()} data.frame used to draw
#'   the plot (without the \code{"TOTAL"} row), in the order plotted.
#' @export
plot_character_stats <- function(character_data, characters = NULL,
                                  sort_by = c("none", "pct_missing", "pct_inapplicable", "character"),
                                  colors = CHAR_STATS_COLORS,
                                  export_filename = NULL,
                                  export_format = c("png", "pdf"),
                                  width = NULL, height = NULL) {
  sort_by       <- match.arg(sort_by)
  export_format <- match.arg(export_format)

  st <- character_stats(character_data, characters)
  st <- st[st$character != "TOTAL", , drop = FALSE]
  st <- switch(sort_by,
    pct_missing      = st[order(-st$pct_missing), ],
    pct_inapplicable = st[order(-st$pct_inapplicable), ],
    character        = st[order(st$character), ],
    st
  )
  rownames(st) <- NULL

  n_char <- nrow(st)
  mat <- t(as.matrix(st[, c("n_scored", "n_missing", "n_inapplicable")]))
  rownames(mat) <- c("scored", "missing", "inapplicable")
  colnames(mat) <- st$character
  storage.mode(mat) <- "numeric"
  # barplot(horiz=TRUE) draws its first column at the bottom; reverse the
  # columns so the first row of `st` (as sorted above) ends up on top.
  mat <- mat[, rev(seq_len(n_char)), drop = FALSE]

  draw <- function() {
    left_mar <- max(6, max(nchar(colnames(mat))) * 0.55 + 2)
    old_par  <- par(mar = c(4, left_mar, 3, 9), xpd = NA)
    on.exit(par(old_par))

    barplot(mat, horiz = TRUE, col = colors[rownames(mat)], border = NA,
            las = 1, cex.names = min(0.9, max(0.35, 30 / n_char)),
            xlab = "Number of taxa",
            main = "Character data completeness")

    usr <- par("usr")
    legend(x = usr[2], y = usr[4], xjust = 0, yjust = 1,
           legend = c("Scored", "Missing (?)", "Inapplicable (-)"),
           fill = colors[c("scored", "missing", "inapplicable")],
           bty = "n", cex = 0.85, xpd = NA)
  }

  if (!is.null(export_filename)) {
    w_in <- width  %||% 8
    h_in <- height %||% max(3, 0.28 * n_char + 1.5)
    prev_dev <- dev.cur()
    .emtree_open_device(export_filename, export_format, w_in, h_in)
    draw()
    .close_device_restore(prev_dev)
    cat("Saved:", paste0(export_filename, ".", export_format), "\n")
  }

  draw()
  invisible(st)
}


#' Taxon x character completeness heatmap
#'
#' Tile plot with one row per taxon and one column per character, colored by
#' whether that cell is scored, missing ("?") or inapplicable ("-"). Gives an
#' at-a-glance view of where missing/inapplicable data cluster across the
#' whole matrix -- complementary to the per-character summary in
#' \code{\link{plot_character_stats}}.
#'
#' @param character_data    Data.frame with a \code{"Species"} column and one
#'                          column per character.
#' @param characters        Optional subset of character columns (default:
#'                          all).
#' @param colors            Named vector with entries \code{scored},
#'                          \code{missing}, \code{inapplicable}. Defaults to
#'                          an Okabe-Ito colorblind-safe triplet.
#' @param show_char_labels  Show character names on the x-axis. Default
#'                          \code{NULL} auto-decides: \code{TRUE} when there
#'                          are \code{<= 60} characters, \code{FALSE}
#'                          otherwise (labels would overlap into
#'                          illegibility).
#' @param export_filename   Optional path (without extension) to also save
#'                          the plot as PNG/PDF. \code{NULL} (default) draws
#'                          on the active graphics device only.
#' @param export_format     \code{"png"} (default) or \code{"pdf"}.
#' @param width,height      Device size in inches when exporting. Defaults
#'                          scale with matrix dimensions.
#' @return Invisibly, \code{NULL}.
#' @export
plot_character_completeness <- function(character_data, characters = NULL,
                                         colors = CHAR_STATS_COLORS,
                                         show_char_labels = NULL,
                                         export_filename = NULL,
                                         export_format = c("png", "pdf"),
                                         width = NULL, height = NULL) {
  if (!is.data.frame(character_data))
    stop("'character_data' must be a data.frame.")
  if (is.null(characters)) characters <- setdiff(colnames(character_data), "Species")
  missing_ch <- setdiff(characters, colnames(character_data))
  if (length(missing_ch) > 0)
    stop("Character(s) not found in 'character_data': ", paste(missing_ch, collapse = ", "))
  if (length(characters) == 0)
    stop("No characters to summarize.")

  export_format <- match.arg(export_format)

  species <- if ("Species" %in% colnames(character_data)) {
    as.character(character_data$Species)
  } else {
    as.character(seq_len(nrow(character_data)))
  }
  n_taxa <- length(species)
  n_char <- length(characters)
  if (is.null(show_char_labels)) show_char_labels <- n_char <= 60

  cls_mat <- vapply(characters, function(ch) .classify_character_values(character_data[[ch]]),
                    character(n_taxa))
  # cls_mat is n_taxa x n_char (one .classify_character_values() result per
  # column); flatten column-major (matches `matrix()`'s default fill order)
  # to recover a same-shaped color matrix.
  code <- matrix(match(as.character(cls_mat), c("scored", "missing", "inapplicable")),
                 nrow = n_taxa, ncol = n_char)
  cell_colors <- matrix(colors[c("scored", "missing", "inapplicable")][code],
                        nrow = n_taxa, ncol = n_char)

  draw <- function() {
    left_mar   <- max(6, max(nchar(species)) * 0.55 + 2)
    bottom_mar <- if (show_char_labels) max(4, max(nchar(characters)) * 0.45 + 2) else 3
    old_par <- par(mar = c(bottom_mar, left_mar, 3, 9), xpd = NA)
    on.exit(par(old_par))

    plot(NA, xlim = c(0, n_char), ylim = c(0, n_taxa), xaxs = "i", yaxs = "i",
         axes = FALSE, xlab = "", ylab = "", main = "Character matrix completeness")

    row_idx  <- seq_len(n_taxa)
    y_top    <- n_taxa - row_idx + 1
    y_bottom <- n_taxa - row_idx
    for (j in seq_len(n_char)) {
      rect(xleft = j - 1, xright = j, ybottom = y_bottom, ytop = y_top,
           col = cell_colors[, j], border = NA)
    }

    axis(2, at = y_bottom + 0.5, labels = species, las = 1,
         cex.axis = min(0.9, max(0.25, 25 / n_taxa)), tick = FALSE, line = -0.5)
    if (show_char_labels) {
      axis(1, at = seq_len(n_char) - 0.5, labels = characters, las = 2,
           cex.axis = min(0.8, max(0.25, 40 / n_char)), tick = FALSE, line = -0.5)
    }

    usr <- par("usr")
    legend(x = usr[2], y = usr[4], xjust = 0, yjust = 1,
           legend = c("Scored", "Missing (?)", "Inapplicable (-)"),
           fill = colors[c("scored", "missing", "inapplicable")],
           bty = "n", cex = 0.85, xpd = NA)
  }

  if (!is.null(export_filename)) {
    w_in <- width  %||% max(6, 0.15 * n_char + 3)
    h_in <- height %||% max(4, 0.22 * n_taxa + 1.5)
    prev_dev <- dev.cur()
    .emtree_open_device(export_filename, export_format, w_in, h_in)
    draw()
    .close_device_restore(prev_dev)
    cat("Saved:", paste0(export_filename, ".", export_format), "\n")
  }

  draw()
  invisible(NULL)
}
