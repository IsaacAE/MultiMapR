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


#' Completes a (possibly partial) scored/missing/inapplicable color vector
#'
#' Entries missing from \code{colors} or not valid R colors fall back to
#' \code{\link{CHAR_STATS_COLORS}}, so callers can override just one of the
#' three categories (e.g. \code{c(missing = "grey80")}).
#' @param colors Named character vector (any subset of \code{scored},
#'   \code{missing}, \code{inapplicable}).
#' @return Named character vector with all three entries.
#' @keywords internal
.resolve_stats_colors <- function(colors) {
  out <- CHAR_STATS_COLORS
  if (is.null(colors)) return(out)
  for (k in intersect(names(colors), names(out))) {
    v <- as.character(colors[[k]])
    if (length(v) == 1 && !is.na(v) && is_valid_color(v)) out[[k]] <- v
  }
  out
}


#' Black or white, whichever reads better on top of each fill color
#' @param fill Vector of R colors.
#' @return Character vector of \code{"#000000"} / \code{"#FFFFFF"}.
#' @keywords internal
.contrast_ink <- function(fill) {
  m <- grDevices::col2rgb(fill) / 255
  lin <- ifelse(m <= 0.03928, m / 12.92, ((m + 0.055) / 1.055)^2.4)
  lum <- 0.2126 * lin[1, ] + 0.7152 * lin[2, ] + 0.0722 * lin[3, ]
  ifelse(lum > 0.179, "#000000", "#FFFFFF")
}


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
#'                         an Okabe-Ito colorblind-safe triplet; any entry
#'                         left out keeps its default.
#' @param border           Color of the outline drawn around every bar
#'                         segment, or \code{NA} for none.
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
                                  border = "grey25",
                                  export_filename = NULL,
                                  export_format = c("png", "pdf"),
                                  width = NULL, height = NULL) {
  sort_by       <- match.arg(sort_by)
  export_format <- match.arg(export_format)
  colors        <- .resolve_stats_colors(colors)

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

    barplot(mat, horiz = TRUE, col = colors[rownames(mat)], border = border,
            lwd = 0.6, las = 1, cex.names = min(0.9, max(0.35, 30 / n_char)),
            xlab = "Number of taxa",
            main = "Character data completeness")

    usr <- par("usr")
    legend(x = usr[2], y = usr[4], xjust = 0, yjust = 1,
           legend = c("Scored", "Missing (?)", "Inapplicable (-)"),
           fill = colors[c("scored", "missing", "inapplicable")],
           border = if (is.na(border)) "grey25" else border,
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


#' Page layout (in inches) of the taxon x character heatmap
#'
#' Shared by \code{\link{plot_character_completeness}} (margins and default
#' export size) and the graphical interface (preview size), so the preview
#' and the exported file are laid out the same way.
#'
#' @param species     Taxon labels (rows).
#' @param characters  Character labels (columns).
#' @param cell_in     Target cell side, in inches.
#' @return List with \code{mai} (bottom, left, top, right margins),
#'   \code{width}, \code{height} (suggested device size) and the logical
#'   flags \code{right_labels} / \code{top_labels}.
#' @keywords internal
.completeness_layout <- function(species, characters, cell_in = 0.2) {
  n_taxa <- length(species); n_char <- length(characters)
  ch_w   <- 0.075                                # ~ width of one glyph at cex 0.8
  sp_w   <- max(nchar(species))    * ch_w + 0.2
  cl_h   <- max(nchar(characters)) * ch_w + 0.2
  right_labels <- n_char > 15                    # repeat taxa on the right
  top_labels   <- n_taxa > 20                    # repeat characters on top
  head_h <- 0.7                                  # title + legend
  mai <- c(cl_h, sp_w, head_h + if (top_labels) cl_h else 0.1,
           if (right_labels) sp_w else 0.25)
  list(mai = mai,
       width  = n_char * cell_in + mai[2] + mai[4],
       height = n_taxa * cell_in + mai[1] + mai[3],
       right_labels = right_labels, top_labels = top_labels)
}


#' Taxon x character completeness heatmap
#'
#' Tile plot with one row per taxon and one column per character, colored by
#' whether that cell is scored, missing ("?") or inapplicable ("-"). Each
#' cell is outlined and, when there is room, labelled with the state coded
#' for that taxon, so the matrix can be read cell by cell. Taxon names are
#' repeated on the right and character names on top for wide / tall
#' matrices, and heavier guide lines every \code{guide_every} rows/columns
#' help follow a row or column across the plot.
#'
#' @param character_data    Data.frame with a \code{"Species"} column and one
#'                          column per character.
#' @param characters        Optional subset of character columns (default:
#'                          all).
#' @param colors            Named vector with entries \code{scored},
#'                          \code{missing}, \code{inapplicable}. Defaults to
#'                          an Okabe-Ito colorblind-safe triplet; any entry
#'                          left out keeps its default.
#' @param border            Color of the cell outlines, or \code{NA} for none.
#' @param show_values       Write the coded state inside each cell.
#'                          \code{NULL} (default) does so whenever the text
#'                          fits the cell.
#' @param guide_every       Draw a heavier guide line every this many rows
#'                          and columns (\code{0} disables them).
#' @param show_char_labels  Show character names. \code{NULL} (default)
#'                          shows them whenever they fit the column width.
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
                                         border = "grey30",
                                         show_values = NULL,
                                         guide_every = 5,
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
  colors        <- .resolve_stats_colors(colors)

  species <- if ("Species" %in% colnames(character_data)) {
    as.character(character_data$Species)
  } else {
    as.character(seq_len(nrow(character_data)))
  }
  n_taxa <- length(species)
  n_char <- length(characters)

  vals <- vapply(characters, function(ch) {
    v <- trimws(as.character(character_data[[ch]]))
    v[is.na(v) | v == ""] <- "?"
    v
  }, character(n_taxa))
  vals <- matrix(vals, nrow = n_taxa, ncol = n_char)
  cls  <- .classify_character_values(vals)          # column-major, like `vals`
  fill <- unname(colors[cls])
  ink  <- .contrast_ink(fill)

  layout <- .completeness_layout(species, characters)

  draw <- function() {
    din <- par("din")
    mai <- layout$mai
    # Shrink the label margins if the device is too small to hold them.
    for (idx in list(c(2, 4), c(1, 3))) {
      avail <- din[if (idx[1] == 2) 1 else 2] - 0.8
      if (sum(mai[idx]) > avail) mai[idx] <- mai[idx] * max(avail, 0.1) / sum(mai[idx])
    }
    old_par <- par(mai = mai, xpd = NA)
    on.exit(par(old_par))

    plot.new()
    plot.window(xlim = c(0, n_char), ylim = c(0, n_taxa), xaxs = "i", yaxs = "i")
    pin    <- par("pin")
    cell_w <- pin[1] / n_char
    cell_h <- pin[2] / n_taxa
    line_h <- par("cin")[2]                             # text line height at cex 1

    xl <- rep(seq_len(n_char) - 1, each = n_taxa)
    yb <- rep(n_taxa - seq_len(n_taxa), times = n_char)
    # Outlines on cells narrower than ~1.5 mm would just grey the whole plot.
    rect(xl, yb, xl + 1, yb + 1, col = fill,
         border = if (min(cell_w, cell_h) < 0.06) NA else border, lwd = 0.6)

    if (guide_every > 0) {
      guide_col <- "grey10"
      if (n_char > guide_every) {
        gx <- seq(guide_every, n_char - 1, by = guide_every)
        segments(gx, 0, gx, n_taxa, col = guide_col, lwd = 1.8)
      }
      if (n_taxa > guide_every) {
        gy <- n_taxa - seq(guide_every, n_taxa - 1, by = guide_every)
        segments(0, gy, n_char, gy, col = guide_col, lwd = 1.8)
      }
    }
    rect(0, 0, n_char, n_taxa, border = "grey15", lwd = 1.2)

    # State written inside each cell
    val_cex <- min(0.8, 0.75 * cell_h / line_h,
                   0.8 * cell_w / (max(nchar(vals)) * 0.6 * par("cin")[1]))
    fits <- val_cex >= 0.35 && min(cell_w, cell_h) >= 0.11
    if (isTRUE(show_values) || (is.null(show_values) && fits)) {
      text(xl + 0.5, yb + 0.5, labels = as.vector(vals), cex = max(val_cex, 0.2),
           col = ink, family = "mono")
    }

    # Taxon labels (left, and right on wide matrices)
    sp_cex <- min(0.8, 0.9 * cell_h / line_h)
    sp_at  <- n_taxa - seq_len(n_taxa) + 0.5
    mtext(species, side = 2, at = sp_at, las = 1, line = 0.3, cex = sp_cex, adj = 1)
    if (layout$right_labels)
      mtext(species, side = 4, at = sp_at, las = 1, line = 0.3, cex = sp_cex, adj = 0)

    # Character labels (bottom, and top on tall matrices)
    ch_cex <- min(0.8, 0.9 * cell_w / line_h)
    if (isTRUE(show_char_labels) || (is.null(show_char_labels) && ch_cex >= 0.3)) {
      ch_at <- seq_len(n_char) - 0.5
      mtext(characters, side = 1, at = ch_at, las = 2, line = 0.3, cex = ch_cex, adj = 1)
      if (layout$top_labels)
        mtext(characters, side = 3, at = ch_at, las = 2, line = 0.3, cex = ch_cex, adj = 0)
    }

    # Title and a horizontal legend in the top margin
    top_y  <- function(inches) grconvertY(din[2] - inches, from = "inches", to = "user")
    left_x <- grconvertX(0.15, from = "inches", to = "user")
    text(left_x, top_y(0.2), "Character matrix completeness", adj = c(0, 0.5),
         font = 2, cex = 1)
    legend(x = left_x, y = top_y(0.47), xjust = 0, yjust = 0.5, horiz = TRUE,
           legend = c("Scored", "Missing (?)", "Inapplicable (-)"),
           fill = colors[c("scored", "missing", "inapplicable")],
           border = if (is.na(border)) "grey25" else border,
           bty = "n", cex = 0.8, xpd = NA)
  }

  if (!is.null(export_filename)) {
    w_in <- width  %||% max(6, layout$width)
    h_in <- height %||% max(4, layout$height)
    prev_dev <- dev.cur()
    .emtree_open_device(export_filename, export_format, w_in, h_in)
    draw()
    .close_device_restore(prev_dev)
    cat("Saved:", paste0(export_filename, ".", export_format), "\n")
  }

  draw()
  invisible(NULL)
}
