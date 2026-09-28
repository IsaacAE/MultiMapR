################################################################################
# Universal export function for MultiMapR.
# Supports three topologies  : phylogram | cladogram | fan
# Supports two formats       : png | pdf
#
# Geometric strategy:
#   - Blank Canvas (edge.color = "transparent") to populate .PlotPhyloEnv
#     without rendering anything, then manually draw the histories with
#     segments() / lines() + offsets calculated on the fly.
#   - Fan uses true polar geometry: radius + angle -> radial segments and arcs.
#
# Dependencies: ape (>= 5.0) -- declared in DESCRIPTION, do not use library() here
# Author:       MultiMapR -- universal export module
################################################################################

# ==============================================================================
# INTERNAL VALIDATORS
# ==============================================================================

#' Checks that `tree` is a valid phylo object with edges and tip.labels
#'
#' @param tree  Object to validate; must be of class \code{phylo}.
#' @keywords internal
.emtree_validate_tree <- function(tree) {
  if (!inherits(tree, "phylo"))
    stop("`tree` must be an object of class 'phylo'.")
  if (is.null(tree$edge) || nrow(tree$edge) == 0L)
    stop("`tree$edge` is empty: the tree has no branches.")
  if (is.null(tree$tip.label))
    stop("`tree$tip.label` is NULL: the tree has no tip labels.")
}

#' Checks that `color_list` is consistent with the number of tree edges
#'
#' @param color_list  List of character vectors (one per history).
#' @param n_edges     Expected length of each color vector (\code{nrow(tree$edge)}).
#' @keywords internal
.emtree_validate_color_list <- function(color_list, n_edges) {
  if (!is.list(color_list) || length(color_list) == 0L)
    stop("`color_list` must be a non-empty list of color vectors.")
  for (i in seq_along(color_list)) {
    vec <- color_list[[i]]
    if (!is.character(vec))
      stop(sprintf("`color_list[[%d]]` must be a character vector.", i))
    if (length(vec) != n_edges)
      stop(sprintf(
        "`color_list[[%d]]` has length %d, but the tree has %d branches. ",
        i, length(vec), n_edges),
        "Each vector must have exactly one entry per edge.")
  }
}

#' Checks scalar type and format arguments
#'
#' @param type    Tree topology string: \code{"phylogram"}, \code{"cladogram"}, or \code{"fan"}.
#' @param format  Output format: \code{"png"} or \code{"pdf"}.
#' @keywords internal
.emtree_validate_type_format <- function(type, format) {
  valid_types   <- c("phylogram", "cladogram", "fan")
  valid_formats <- c("png", "pdf")
  if (!is.character(type)   || length(type) != 1L || !type   %in% valid_types)
    stop(sprintf("`type` must be one of: %s.", paste(valid_types, collapse = ", ")))
  if (!is.character(format) || length(format) != 1L || !format %in% valid_formats)
    stop(sprintf("`format` must be one of: %s.", paste(valid_formats, collapse = ", ")))
}

#' Checks the tip-label options shared by the export and screen engines
#'
#' @param show_labels  Logical scalar.
#' @param label_size   Positive multiplier of the automatic label size.
#' @param label_color  Single valid R color.
#' @keywords internal
.emtree_validate_labels <- function(show_labels, label_size, label_color) {
  if (!is.logical(show_labels) || length(show_labels) != 1L || is.na(show_labels))
    stop("`show_labels` must be TRUE or FALSE.")
  if (!is.numeric(label_size) || length(label_size) != 1L || is.na(label_size) || label_size <= 0)
    stop("`label_size` must be a positive number.")
  if (!is.character(label_color) || length(label_color) != 1L || !is_valid_color(label_color))
    stop("`label_color` must be a single valid R color name or hex code.")
}

#' Checks that `filename` is a non-empty string (no extension required)
#'
#' @param filename  Output filename without extension; must be a non-empty string.
#' @keywords internal
.emtree_validate_filename <- function(filename) {
  if (!is.character(filename) || length(filename) != 1L || !nzchar(trimws(filename)))
    stop("`filename` must be a non-empty string (without extension).")
}


# ==============================================================================
# INTERNAL GEOMETRIC HELPERS
# ==============================================================================

#' Safely retrieves .PlotPhyloEnv from the ape namespace
.emtree_get_PlotPhyloEnv <- function() {
  get(".PlotPhyloEnv", envir = asNamespace("ape"))
}

#' Opens the appropriate graphics device according to `format`
#'
#' @param filename   Output path without extension.
#' @param format     Device format: \code{"png"} or \code{"pdf"}.
#' @param width_in   Device width in inches.
#' @param height_in  Device height in inches.
#' @return Invisible NULL. The device remains open after the call.
#' @keywords internal
.emtree_open_device <- function(filename, format, width_in, height_in) {

  path <- paste0(filename, ".", format)
  if (format == "png") {
    png(filename = path,
        width    = width_in,
        height   = height_in,
        units    = "in",
        res      = 300,
        bg       = "white")
  } else {                       # pdf
    pdf(file   = path,
        width  = width_in,
        height = height_in)
  }
  invisible(NULL)
}

#' Computes the offset vector from the current physical device parameters
#' (par("usr") / par("pin")).
#'
#' @param N             Number of histories.
#' @param offset_range  Amplitude of the Y range (default 0.1).
#' @return List with `dx` and `dy` (vectors of length N).
.emtree_calc_offsets <- function(N, offset_range = 0.1) {
  if (N == 1L)
    return(list(dx = 0, dy = 0))

  # Fixed step = offset_range between consecutive layers, set centered at 0.
  # This ensures that for any N the visual separation between branches
  # is always the same (offset_range), without gaps or asymmetries:
  #   N=2 -> c(-range/2, +range/2)  separation = range [ok]
  #   N=3 -> c(-range,   0, +range) separation = range [ok]  (same as seq() gave)
  #   N=4 -> c(-3r/2, -r/2, +r/2, +3r/2) separation = range [ok]
  #
  # Original bug: seq(-r, r, N=2) -> c(-r, +r), separation = 2*range and
  # no value at 0, leaving a gap twice the step size of N=3.
  indices    <- seq_len(N) - (N + 1L) / 2   # e.g. N=2->c(-0.5,0.5); N=3->c(-1,0,1)
  offsets_y  <- indices * offset_range

  usr <- par("usr")   # c(x1, x2, y1, y2) -- coordinate system limits
  pin <- par("pin")   # c(width_in, height_in) -- physical plot area dimensions

  scale_x <- (usr[2L] - usr[1L]) / pin[1L]   # data_units / inch (X)
  scale_y <- (usr[4L] - usr[3L]) / pin[2L]   # data_units / inch (Y)

  offsets_x <- offsets_y * (scale_x / scale_y)

  list(dx = offsets_x, dy = offsets_y)
}


# ==============================================================================
# LEGEND HELPER -- EXCLUSIVE SPACE
# ==============================================================================

#' Converts the user-chosen corner into margins and legend position
#'
#' Strategy: the margin of the chosen corner is enlarged so the legend
#' never overlaps with labels or the plot. The position returned
#' by `legend()` uses `par(xpd = TRUE)` over the outer margin panel.
#'
#' @param corner        String: "topleft" | "topright" | "bottomleft" | "bottomright"
#' @param mar_base      Numeric vector c(b, l, t, r) of base margins in lines.
#' @param legend_h      Estimated legend height in lines (default 4).
#' @param legend_w      Estimated legend width in lines (default 6).
#' @param n_chars       Number of characters in the legend (proportionally expands
#'                      `legend_h` when there are multiple blocks). Default 1.
#' @return List with `mar` (new margin vector) and `corner` (string
#'         as expected by `legend()`).
.emtree_config_legend <- function(corner, mar_base, legend_h = 4,
                                  legend_w = 6, n_chars = 1L) {
  # Kept for compatibility with calls in main_MultiMapR.R (simple_mapping).
  # In export_multimapr_tree() space reservation is done with
  # .emtree_measure_legend() after the first render, expanding xlim/ylim.
  valid_corners <- c("topleft", "topright", "bottomleft", "bottomright")
  if (!corner %in% valid_corners)
    stop(sprintf("`legend_corner` must be one of: %s.",
                 paste(valid_corners, collapse = ", ")))

  legend_h_total <- legend_h * max(1L, n_chars) + max(0L, n_chars - 1L)

  mar <- mar_base
  if (grepl("top",    corner)) mar[3L] <- max(mar[3L], legend_h_total)
  if (grepl("bottom", corner)) mar[1L] <- max(mar[1L], legend_h_total)
  if (grepl("left",   corner)) mar[2L] <- max(mar[2L], legend_w)
  if (grepl("right",  corner)) mar[4L] <- max(mar[4L], legend_w)

  list(mar = mar, corner = corner)
}

#' Horizontal gap between a legend symbol and its label, in character widths
#' (\code{legend(x.intersp = )}; R's default of 1 leaves the text too close).
#' Used for drawing and measuring every legend so reserved space matches.
#' @keywords internal
.LEGEND_X_INTERSP <- 1.6

#' Legend spacing that grows with the symbol size
#'
#' \code{legend()} spaces rows and places the label from the TEXT size, so a
#' symbol much larger than the text (\code{pt.cex >> cex}, e.g. big terminal
#' figures) overlaps the next row and its own label. Both spacings are scaled
#' by the symbol/text ratio and never drop below the defaults.
#'
#' @param cex    Legend text size.
#' @param pt_cex Legend symbol size.
#' @return List with \code{x} (\code{x.intersp}) and \code{y} (\code{y.intersp}).
#' @keywords internal
.legend_intersp <- function(cex, pt_cex) {
  r <- max(1, pt_cex / cex)
  list(x = max(.LEGEND_X_INTERSP, 0.3 * r + 1.1),
       y = max(1, 0.5 * r))
}

#' Positions of the per-character legend blocks
#'
#' Measures every block with \code{legend(..., plot = FALSE)} and lays them out
#' together, so symbols, state names and titles line up across characters:
#' \itemize{
#'   \item \code{"vertical"}: blocks stacked in character order (first on top),
#'     all starting at the same x.
#'   \item \code{"horizontal"}: blocks side by side in character order, all
#'     starting at the same top y.
#' }
#' The whole group is anchored to \code{corner} of the plot region. Blocks are
#' then drawn with \code{xjust = 0, yjust = 1} at the returned coordinates.
#'
#' @param blocks  List of \code{list(title, labels)}.
#' @param corner  "topleft" | "topright" | "bottomleft" | "bottomright".
#' @param layout  "vertical" or "horizontal".
#' @param cex,pt_cex,pch Text size, symbol size and symbol used by the legend.
#' @return List: \code{x}, \code{y} (top-left of every block, data units) and
#'   \code{w}, \code{h} (size of the whole group).
#' @keywords internal
.legend_layout <- function(blocks, corner, layout = "vertical",
                           cex = 0.8, pt_cex = cex * 1.2, pch = 15) {
  layout <- match.arg(layout, c("vertical", "horizontal"))
  usr <- par("usr")
  n   <- length(blocks)
  if (n == 0L) return(list(x = numeric(0), y = numeric(0), w = 0, h = 0))
  sp <- .legend_intersp(cex, pt_cex)
  # Titles are drawn apart from legend() (see .legend_draw_blocks()): with a
  # built-in title wider than its entries, legend() centres the entries under
  # it and the symbols of different blocks would no longer line up.
  rects <- lapply(blocks, function(b)
    legend(x = usr[1L], y = usr[4L], legend = b$labels, pch = pch,
           x.intersp = sp$x, y.intersp = sp$y,
           bty = "n", cex = cex, pt.cex = pt_cex,
           xjust = 0, yjust = 1, plot = FALSE)$rect)
  title_h <- .legend_title_height(cex)
  title_w <- vapply(blocks, function(b)
    if (is.null(b$title)) 0 else strwidth(b$title, cex = cex) + .legend_title_pad(cex), numeric(1))
  w <- pmax(vapply(rects, function(r) r$w, numeric(1)), title_w)
  h <- vapply(rects, function(r) r$h, numeric(1)) + title_h
  gap_y <- strheight("M", cex = cex) * 1.4
  gap_x <- strwidth("MM", cex = cex)
  if (layout == "horizontal") {
    W <- sum(w) + gap_x * (n - 1L); H <- max(h)
    x <- c(0, cumsum(w + gap_x))[seq_len(n)]
    y <- rep(0, n)
  } else {
    W <- max(w); H <- sum(h) + gap_y * (n - 1L)
    x <- rep(0, n)
    y <- -c(0, cumsum(h + gap_y))[seq_len(n)]
  }
  x_left <- if (grepl("right", corner)) usr[2L] - W else usr[1L]
  y_top  <- if (grepl("top",   corner)) usr[4L]     else usr[3L] + H
  list(x = x_left + x, y = y_top + y, w = W, h = H)
}

#' Height reserved for a legend block title (data units, current device)
#' @param cex Legend text size.
#' @keywords internal
.legend_title_height <- function(cex) strheight("M", cex = cex) * 1.6

#' Left offset of a legend block title, so it starts where the symbols do
#' @param cex Legend text size.
#' @keywords internal
.legend_title_pad <- function(cex) strwidth("m", cex = cex) * 0.5

#' Draws legend blocks at the positions computed by \code{.legend_layout()}
#'
#' Each block is its title (left-aligned, drawn with \code{text()}) above a
#' title-less \code{legend()}, so symbols, state names and titles line up
#' across blocks.
#'
#' @param lay     Output of \code{.legend_layout()}.
#' @param blocks  List of \code{list(title, labels, colors)}.
#' @param cex,pt_cex,pch Text size, symbol size and symbol.
#' @param filled  \code{TRUE} for fillable symbols (21-25): colors go to
#'   \code{pt.bg} with a black border; otherwise to \code{col}.
#' @keywords internal
.legend_draw_blocks <- function(lay, blocks, cex, pt_cex, pch, filled = FALSE) {
  title_h <- .legend_title_height(cex)
  pad     <- .legend_title_pad(cex)
  sp      <- .legend_intersp(cex, pt_cex)
  for (k in seq_along(blocks)) {
    b <- blocks[[k]]
    if (!is.null(b$title))
      text(lay$x[k] + pad, lay$y[k], labels = b$title, adj = c(0, 1), cex = cex)
    args <- list(x = lay$x[k], y = lay$y[k] - title_h, legend = b$labels, pch = pch,
                 x.intersp = sp$x, y.intersp = sp$y, bty = "n", cex = cex, pt.cex = pt_cex,
                 horiz = FALSE, xjust = 0, yjust = 1)
    if (filled) { args$pt.bg <- b$colors; args$col <- "black" } else args$col <- b$colors
    do.call(legend, args)
  }
  invisible(NULL)
}

#' Blocks of a grouped legend in the form expected by \code{.legend_layout()}
#' @param legend_by_char Named list character -> list(labels, colors).
#' @keywords internal
.legend_blocks <- function(legend_by_char) {
  keep <- vapply(legend_by_char, function(b) length(b$labels) > 0L, logical(1))
  Map(function(nm, b) list(title = nm, labels = b$labels, colors = b$colors),
      names(legend_by_char)[keep], legend_by_char[keep])
}

#' Measures the space occupied by the legend in data coordinates using
#' `legend(..., plot = FALSE)`, which calculates the rect without rendering.
#'
#' Must be called **after** the tree has already been drawn and par("usr") /
#' par("pin") reflect the actual canvas. Returns how many additional data units
#' need to be reserved on the X axis (left or right) and on the Y axis
#' (up or down) so the legend does not overlap the plot.
#'
#' @param corner         Corner string.
#' @param legend_by_char Grouped list character -> list(labels, colors), or NULL.
#' @param legend_labels  Flat label vector (legacy mode).
#' @param legend_title   Block title (legacy mode).
#' @param cex_ley        Legend font size.
#' @param legend_layout  "vertical" or "horizontal" (grouped legend only).
#' @return List: `dx` (extra X units), `dy` (extra Y units),
#'         `going_down` (bool), `on_right` (bool).
.emtree_measure_legend <- function(corner,
                                   legend_by_char = NULL,
                                   legend_labels  = NULL,
                                   legend_title   = NULL,
                                   cex_ley        = 0.8,
                                   legend_layout  = "vertical") {

  going_down <- grepl("top",   corner)
  on_right   <- grepl("right", corner)
  usr        <- par("usr")   # c(x1, x2, y1, y2) in data coordinates

  # -- Internal function: measures a legend block with plot=FALSE -------------
  measure_block <- function(labels, title_txt, x_ref, y_ref) {
    lg <- legend(x.intersp = .LEGEND_X_INTERSP,
                 x      = x_ref,
                 y      = y_ref,
                 legend = labels,
                 pch    = 15,
                 col    = rep("black", length(labels)),
                 title  = title_txt,
                 bty    = "n",
                 cex    = cex_ley,
                 horiz  = FALSE,
                 pt.cex = cex_ley * 1.2,
                 xjust  = if (on_right) 1 else 0,
                 yjust  = if (going_down) 1 else 0,
                 plot   = FALSE)   # <<< measures without drawing
    lg$rect
  }

  # Reference position: the corner of the data area
  x_ref <- if (on_right)   usr[2L] else usr[1L]
  y_ref <- if (going_down) usr[4L] else usr[3L]

  total_w <- 0
  total_h <- 0
  line_h  <- strheight("M", cex = cex_ley) * 1.4

  if (!is.null(legend_by_char) && length(legend_by_char) > 0L) {
    lay     <- .legend_layout(.legend_blocks(legend_by_char), corner, legend_layout,
                              cex = cex_ley, pt_cex = cex_ley * 1.2, pch = 15)
    total_w <- lay$w
    total_h <- lay$h
  } else if (!is.null(legend_labels) && length(legend_labels) > 0L) {
    rect    <- measure_block(legend_labels, legend_title, x_ref, y_ref)
    total_w <- rect$w
    total_h <- rect$h
  } else {
    return(list(dx = 0, dy = 0, going_down = going_down, on_right = on_right))
  }

  # Safety margin: 8% extra so the legend does not touch the border
  list(dx         = total_w * 1.08,
       dy         = total_h * 1.08,
       going_down = going_down,
       on_right   = on_right)
}

#' Draws the legend in the reserved corner of the margin
#'
#' Supports two modes:
#'   - **Grouped** (`legend_by_char` not NULL): receives a named list
#'     character -> list(labels, colors) and draws each character as an
#'     independent block with its header in bold, without mixing states.
#'   - **Flat** (legacy): uses `legend_labels` + `legend_colors` + `legend_title`
#'     as before, to maintain compatibility with older calls.
#'
#' Called **after** rendering the tree. Uses `par(xpd = NA)` to
#' paint outside the plot area but within the device.
#'
#' @param corner         Corner string ("topleft", etc.)
#' @param legend_by_char Named list character -> list(labels, colors). If not
#'                       NULL takes precedence over `legend_labels`.
#' @param legend_labels  Flat label vector (legacy mode).
#' @param legend_colors  Flat color vector (legacy mode).
#' @param legend_title   Single block title (legacy mode, string or NULL).
#' @param pch            Legend symbol (default 15 = filled square).
#' @param cex_ley        Font size for the legend.
#' @param legend_layout  "vertical" (blocks stacked) or "horizontal" (side by side).
.emtree_draw_legend <- function(corner,
                                legend_by_char = NULL,
                                legend_labels  = NULL,
                                legend_colors  = NULL,
                                legend_title   = NULL,
                                pch            = 15,
                                cex_ley        = 0.8,
                                legend_layout  = "vertical") {

  old_xpd <- par("xpd")
  par(xpd = NA)
  on.exit(par(xpd = old_xpd))

  # -- Grouped mode: one block per character with header ---------------------
  # All blocks are laid out together (.legend_layout) and drawn left-justified
  # with left-aligned titles, so symbols and state names line up.
  if (!is.null(legend_by_char) && length(legend_by_char) > 0L) {
    blocks <- .legend_blocks(legend_by_char)
    lay    <- .legend_layout(blocks, corner, legend_layout,
                             cex = cex_ley, pt_cex = cex_ley * 1.2, pch = pch)
    .legend_draw_blocks(lay, blocks, cex = cex_ley, pt_cex = cex_ley * 1.2, pch = pch)
    return(invisible(NULL))
  }

  # -- Legacy (flat) mode: a single block ------------------------------------
  if (is.null(legend_labels) || length(legend_labels) == 0L)
    return(invisible(NULL))

  legend(corner, x.intersp = .LEGEND_X_INTERSP,
         legend = legend_labels,
         col    = legend_colors,
         pch    = pch,
         title  = legend_title,
         bty    = "n",
         cex    = cex_ley,
         horiz  = FALSE,
         pt.cex = cex_ley * 1.2)

  invisible(NULL)
}


# ==============================================================================
# TOPOLOGY RENDERERS
# ==============================================================================

# ------------------------------------------------------------------------------
# PHYLOGRAM
# Vectorized horizontals; verticals in a loop protected against degeneration.
# ------------------------------------------------------------------------------

.emtree_render_phylogram <- function(pp, tree, color_list, lwd, offsets) {
  xx <- pp$xx
  yy <- pp$yy

  edges      <- tree$edge
  parent_idx <- edges[, 1L]
  child_idx  <- edges[, 2L]
  n_tips     <- Ntip(tree)

  # Internal nodes and their metadata
  internal_nodes <- unique(edges[edges[, 1L] > n_tips, 1L])
  node_x         <- xx[internal_nodes]
  entry_idx      <- match(internal_nodes, edges[, 2L])   # NA = root
  children_of    <- lapply(internal_nodes,
                           function(nd) edges[edges[, 1L] == nd, 2L])

  N <- length(color_list)

  for (i in seq_len(N)) {
    ec    <- color_list[[i]]
    lwd_i <- lwd[min(i, length(lwd))]
    dy    <- offsets$dy[min(i, length(offsets$dy))]
    dx    <- offsets$dx[min(i, length(offsets$dx))]

    # -- Horizontals (vectorized, no loop) ------------------------------------
    # Each edge goes from x_parent to x_child, at the height y_child (ape convention).
    segments(x0   = xx[parent_idx] + dx,
             y0   = yy[child_idx]  + dy,
             x1   = xx[child_idx]  + dx,
             y1   = yy[child_idx]  + dy,
             col  = ec,
             lwd  = lwd_i,
             lend = 1L)

    # -- Verticals (loop over internal nodes) ---------------------------------
    # Connects children of a node at the node's X coordinate.
    for (m in seq_along(internal_nodes)) {
      ie         <- entry_idx[m]
      node_color <- if (!is.na(ie)) {
        ec[ie]
      } else {
        # Root has no incoming edge. Prefer the algorithm-supplied root
        # color (its own resolved state) over borrowing the first child's
        # edge color, which was the previous (incorrect) fallback.
        rc <- attr(ec, "root_color")
        if (!is.null(rc)) rc else ec[which(edges[, 1L] == internal_nodes[m])[1L]]
      }

      children_m <- children_of[[m]]
      y_children <- yy[children_m]
      y_min      <- min(y_children)
      y_max      <- max(y_children)

      if (y_min != y_max) {   # Skip degenerate segments (polytomies / length 0)
        segments(x0   = node_x[m] + dx,
                 y0   = y_min     + dy,
                 x1   = node_x[m] + dx,
                 y1   = y_max     + dy,
                 col  = node_color,
                 lwd  = lwd_i,
                 lend = 1L)
      }
    }
  }
}


# ------------------------------------------------------------------------------
# CLADOGRAM
# Diagonals parent -> child, one parallel "lane" per history.
#
# Lanes: each history is a filled band parallel to the branch (not a thick
# line). Bands are offset PERPENDICULAR to every branch, measured in inches, so
# the ribbon looks the same on rising and falling branches and for any canvas
# proportion. A global dx/dy shift would run almost parallel to rising branches
# and collapse their lanes.
#   · Band edges are lines offset along the branch's left normal (the "upper"
#     side for both slopes, so lane order is stable). Neighbouring bands share
#     an edge, so there is never a gap between lanes.
#   · Miter joins: each edge of a child band starts where it crosses the same
#     edge of the parent's band (at the root: of a sibling's band), so corners
#     close without notches.
#   · Guard: the parent band ends on the join cut (the bisector between the
#     incoming branch and each child); beyond the node it only runs under the
#     band of the collinear child, which is drawn later and hides it. It never
#     spills onto a sibling clade. Edges are drawn root-to-tips so children
#     always sit on top of their parent at the joins, and among siblings the
#     straight continuation is drawn last, covering the outer miter corner of
#     the side branch.
#   · At the tips every lane centre ends on the tip's vertical (lanes stay
#     aligned) and the band is cut perpendicular to the branch; only there
#     does `tip_end` ("round" | "butt" | "square") apply.
#
# Branch angle (`branch_angle`, degrees from the horizontal, default 45): the
# tree is compressed horizontally toward the tips so every branch runs at
# +/-branch_angle on the device; sibling branches meet at twice that angle
# (90 degrees at the default). Tips (and therefore labels) do not move. If the
# canvas is too narrow for the requested angle, the closest angle that fits is
# used. `branch_angle = NULL` (or NA) stretches the tree to the full width.
# Lane corners sit about d / tan(angle) from the node, so angles below ~30
# degrees push them away from it; 30-60 degrees is the recommended range.
#
# ------------------------------------------------------------------------------

#' Extra room (data units) above the plot for rotated column headers
#'
#' The terminal-figure overlay draws one header per character above the top
#' tip, rotated 45 degrees (see \code{plot_ancestral_with_terminals()}). The
#' canvas must reserve that height or the headers are clipped. Measured on the
#' active device right after the canvas is drawn; accounts for the fact that
#' enlarging the y range also enlarges each inch in data units.
#' @param labels  Header strings.
#' @param cex     Header text size.
#' @param yy_tips Tip y coordinates.
#' @return Non-negative amount to add to the upper y limit.
#' @keywords internal
.emtree_header_space <- function(labels, cex, yy_tips) {
  if (length(labels) == 0L) return(0)
  usr   <- par("usr")
  pin_h <- par("pin")[2L]
  rng   <- usr[4L] - usr[3L]
  w_in  <- max(strwidth(labels, units = "inches", cex = cex, font = 2L))
  h_in  <- max(strheight(labels, units = "inches", cex = cex, font = 2L))
  a_in  <- w_in * sin(pi / 4) + h_in               # rotated text height + padding
  base  <- max(yy_tips) + diff(range(yy_tips)) * 0.04
  if (a_in >= pin_h) return(0)                     # cannot fit anyway
  max(0, (base + a_in * rng / pin_h - usr[4L]) / (1 - a_in / pin_h))
}

#' Normalizes and validates a cladogram branch angle
#'
#' @param branch_angle Degrees from the horizontal (10-80), or \code{NULL} /
#'   \code{NA} to stretch the tree to the full width.
#' @return The angle as a number, or \code{NULL} for the stretched layout.
#' @keywords internal
.emtree_branch_angle <- function(branch_angle) {
  if (is.null(branch_angle) || (length(branch_angle) == 1L && is.na(branch_angle)))
    return(NULL)
  if (!is.numeric(branch_angle) || length(branch_angle) != 1L ||
      branch_angle < 10 || branch_angle > 80)
    stop("`branch_angle` must be a number of degrees between 10 and 80, or NULL.")
  branch_angle
}

#' Automatic figure width in inches
#'
#' A cladogram drawn at a fixed branch angle is narrower than it is tall (about
#' half as wide at 45 degrees), so a fixed 12 in canvas would leave most of the
#' figure empty; its width follows the height and the angle instead (tree +
#' labels). Other topologies, and the stretched cladogram, keep 12 in.
#' @param type         Tree topology.
#' @param height       Figure height in inches.
#' @param branch_angle Cladogram branch angle in degrees, or \code{NULL}.
#' @keywords internal
.emtree_auto_width <- function(type, height, branch_angle = 45) {
  branch_angle <- .emtree_branch_angle(branch_angle)
  if (!identical(type, "cladogram") || is.null(branch_angle)) return(12)
  max(6, 0.55 * height / tan(branch_angle * pi / 180) + 4)
}

.emtree_lend_code <- function(tip_end) {
  switch(tip_end %||% "round", "round" = 0L, "butt" = 1L, "square" = 2L,
         stop("`tip_end` must be one of \"round\", \"butt\" or \"square\"."))
}

#' Geometry of a cladogram line offset by `d` inches (a band edge or lane center)
#'
#' @param X,Y   Node coordinates in inches.
#' @param edges Edge matrix of the tree.
#' @param n_tips Number of tips.
#' @param d     Signed perpendicular distance in inches.
#' @return list(x0, y0, x1, y1) in inches, one element per edge.
#' @keywords internal
.emtree_clado_lanes <- function(X, Y, edges, n_tips, d) {
  par_i <- edges[, 1L]
  chi_i <- edges[, 2L]
  E     <- nrow(edges)

  ux  <- X[chi_i] - X[par_i]
  uy  <- Y[chi_i] - Y[par_i]
  len <- sqrt(ux^2 + uy^2)
  len[len == 0] <- 1
  ux <- ux / len; uy <- uy / len
  # Lane anchor: parent shifted along the left normal (-uy, ux)
  ax <- X[par_i] - d * uy
  ay <- Y[par_i] + d * ux

  # Parameter t along lane `e` where it crosses lane `f` (NA if parallel)
  cross_t <- function(e, f) {
    cr <- ux[e] * uy[f] - uy[e] * ux[f]
    if (abs(cr) < 1e-9) return(NA_real_)
    ((ax[f] - ax[e]) * uy[f] - (ay[f] - ay[e]) * ux[f]) / cr
  }

  entry  <- match(par_i, chi_i)            # incoming edge of each edge's parent
  t0 <- numeric(E)
  for (e in seq_len(E)) {
    f <- entry[e]
    if (is.na(f)) {                        # root: pair with a sibling lane
      sib <- which(par_i == par_i[e] & seq_len(E) != e)
      f   <- if (length(sib)) sib[1L] else NA_integer_
    }
    t0[e] <- if (is.na(f)) 0 else cross_t(e, f)
    # Miter limit: nearly collinear lanes cross far away -> bevel instead
    if (is.na(t0[e]) || abs(t0[e]) > 4 * abs(d)) t0[e] <- 0
  }

  t1 <- len                                # default: shifted child node
  for (e in seq_len(E)) {
    ch <- chi_i[e]
    if (ch <= n_tips) {
      # Cut on the tip vertical so every lane ends at the same x
      if (abs(ux[e]) > 1e-9) t1[e] <- (X[ch] - ax[e]) / ux[e]
    } else {
      kids <- which(par_i == ch)
      # Children start on this lane (t0 of kid measured on the kid's lane);
      # map those points back onto this lane and extend to the farthest one.
      tk <- vapply(kids, function(k) {
        (ax[k] + t0[k] * ux[k] - ax[e]) * ux[e] + (ay[k] + t0[k] * uy[k] - ay[e]) * uy[e]
      }, numeric(1))
      t1[e] <- max(c(len[e], tk))
    }
  }

  list(x0 = ax + t0 * ux, y0 = ay + t0 * uy,
       x1 = ax + t1 * ux, y1 = ay + t1 * uy)
}

.emtree_render_cladogram <- function(pp, tree, color_list, lwd, offsets,
                                     tip_end = "round", branch_angle = 45) {
  edges      <- tree$edge
  parent_idx <- edges[, 1L]
  child_idx  <- edges[, 2L]
  n_tips     <- Ntip(tree)
  lend_code  <- .emtree_lend_code(tip_end)

  # Work in inches so perpendicular distances are isotropic on the device
  usr <- par("usr"); pin <- par("pin")
  sx  <- (usr[2L] - usr[1L]) / pin[1L]
  sy  <- (usr[4L] - usr[3L]) / pin[2L]
  xx  <- pp$xx
  branch_angle <- .emtree_branch_angle(branch_angle)
  if (!is.null(branch_angle)) {
    dxe   <- xx[child_idx] - xx[parent_idx]
    slope <- stats::median(abs((pp$yy[child_idx] - pp$yy[parent_idx])[dxe != 0] / dxe[dxe != 0]))
    x_tip <- max(xx)
    span  <- x_tip - min(xx)
    if (is.finite(slope) && slope > 0 && span > 0) {
      # device slope = tan(branch_angle)
      k     <- slope * sx / sy / tan(branch_angle * pi / 180)
      k_max <- (x_tip - usr[1L]) / span            # keep the root inside the plot
      xx    <- x_tip - (x_tip - xx) * min(k, k_max)
    }
  }
  X   <- xx / sx
  Y   <- pp$yy / sy

  N <- length(color_list)

  # Band edges: lane widths are the line widths (lwd 1 = 1/96 in), stacked
  # without gaps and centered on the branch. A hair of overlap between
  # neighbouring bands hides anti-aliasing seams.
  w_in   <- lwd[pmin(seq_len(N), length(lwd))] / 96
  bounds <- c(0, cumsum(w_in)) - sum(w_in) / 2
  seam   <- 0.004                                  # inches (~1 px at 300 dpi)
  geo    <- lapply(bounds, function(o) .emtree_clado_lanes(X, Y, edges, n_tips, o))
  centre <- lapply(seq_len(N), function(i) .emtree_clado_lanes(
    X, Y, edges, n_tips, (bounds[i] + bounds[i + 1L]) / 2))

  # Unit direction and left normal of every edge (inches)
  ux <- X[child_idx] - X[parent_idx]; uy <- Y[child_idx] - Y[parent_idx]
  ul <- sqrt(ux^2 + uy^2); ul[ul == 0] <- 1
  ux <- ux / ul; uy <- uy / ul
  nx <- -uy; ny <- ux

  # Drawing order: root to tips, so a child band always covers its parent's
  # join. Among siblings, the child that continues the incoming branch in a
  # straight line is drawn last (guard): the outer miter corner of a side
  # branch lies past the node, inside the continuation's band, and must not
  # show on top of that sibling clade.
  depth     <- ape::node.depth.edgelength(ape::compute.brlen(tree, 1))
  entry     <- match(parent_idx, child_idx)
  collinear <- !is.na(entry) &
    (ux * ux[entry] + uy * uy[entry]) > cos(2 * pi / 180)
  collinear[is.na(collinear)] <- FALSE
  draw  <- order(depth[parent_idx], collinear)
  is_tip_edge <- child_idx <= n_tips
  continues   <- child_idx %in% parent_idx[collinear]   # edge has a straight continuation

  for (e in draw) {
    for (i in seq_len(N)) {
      lo <- geo[[i]]; hi <- geo[[i + 1L]]
      # widen by `seam` along the band normal so adjacent bands overlap slightly
      sx_e <- nx[e] * seam; sy_e <- ny[e] * seam
      if (is_tip_edge[e]) {
        # Tip: cut perpendicular to the branch through the lane centre, which
        # sits on the tip vertical, so caps are clean and lanes stay aligned.
        cx <- centre[[i]]$x1[e]; cy <- centre[[i]]$y1[e]
        hw <- w_in[i] / 2
        end_lo <- c(cx - nx[e] * hw, cy - ny[e] * hw)
        end_hi <- c(cx + nx[e] * hw, cy + ny[e] * hw)
      } else if (continues[e]) {
        # Guard: a straight continuation covers everything past the node, so
        # the band stops on the node's perpendicular cut and never shows next
        # to it (the side branch's outer corner lies under the continuation).
        cn <- child_idx[e]
        end_lo <- c(X[cn] + nx[e] * bounds[i],      Y[cn] + ny[e] * bounds[i])
        end_hi <- c(X[cn] + nx[e] * bounds[i + 1L], Y[cn] + ny[e] * bounds[i + 1L])
      } else {
        end_lo <- c(lo$x1[e], lo$y1[e])
        end_hi <- c(hi$x1[e], hi$y1[e])
      }
      px <- c(lo$x0[e] - sx_e, end_lo[1L] - sx_e, end_hi[1L] + sx_e, hi$x0[e] + sx_e)
      py <- c(lo$y0[e] - sy_e, end_lo[2L] - sy_e, end_hi[2L] + sy_e, hi$y0[e] + sy_e)
      polygon(px * sx, py * sy, col = color_list[[i]][e], border = NA)
    }
  }

  # Tip ends: "butt" keeps the flat perpendicular cut; "round" / "square" add a cap.
  if (lend_code != 1L) {
    eps <- 1e-4 * max(diff(range(X)), diff(range(Y)), 1e-8)
    for (i in seq_len(N)) {
      ct <- centre[[i]]
      vx <- ct$x1 - ct$x0; vy <- ct$y1 - ct$y0
      vl <- sqrt(vx^2 + vy^2); vl[vl == 0] <- 1
      e  <- which(is_tip_edge)
      segments(x0   = (ct$x1[e] - vx[e] / vl[e] * eps) * sx,
               y0   = (ct$y1[e] - vy[e] / vl[e] * eps) * sy,
               x1   = ct$x1[e] * sx, y1 = ct$y1[e] * sy,
               col  = color_list[[i]][e],
               lwd  = lwd[min(i, length(lwd))],
               lend = lend_code)
    }
  }
}


# ------------------------------------------------------------------------------
# FAN
# Concentric polar geometry with compact separation calibrated to the tree.
#
# Spread scale:
#   Layers should appear nearly parallel, like a superimposed adjusted
#   multimapping, not as separate rings. To achieve this, the total spread
#   is expressed as a small percentage of R_max (tree maximum radius):
#
#     spread  = R_max * 0.03   -> 3% of total radius, distributed among N layers
#     vect_D  = seq(-spread/2, spread/2, length.out = N)
#
#   For each history i:
#     dr       = vect_D[i]            radial displacement
#     dtheta_n = dtheta_node[n] * dr/R_max  angular rotation per node, so
#                                     that each trace endpoint uses its own
#                                     angle -> traces parallel to the original
#
# East axis crossing correction (0 deg / 360 deg):
#   Normalization to [0, 2pi) with ifelse(); if ang_max - ang_min > pi the clade
#   crosses 0 deg and angles < pi are raised to the extended range [2pi, 4pi).
# ------------------------------------------------------------------------------

.emtree_render_fan <- function(pp, tree, color_list, lwd, offsets) {
  xx     <- pp$xx
  yy     <- pp$yy
  edges  <- tree$edge
  n_tips <- Ntip(tree)
  N      <- length(color_list)

  # Extract original angles for all nodes
  node_angles <- atan2(yy, xx)

  # -- Internal node metadata ------------------------------------------------
  internal_nodes <- unique(edges[edges[, 1L] > n_tips, 1L])
  entry_idx      <- match(internal_nodes, edges[, 2L])
  children_of    <- lapply(internal_nodes, function(nd) edges[edges[, 1L] == nd, 2L])

  # -- Spread ----------------------------------------------------------------
  R_max  <- max(sqrt(xx^2 + yy^2))
  # The spread between layers in fan mode scales with branch width (lwd)
  # so that wider layers do not overlap. Empirical factor: 0.0025 * lwd,
  # normalized by R_max to be independent of tree size.
  # The lwd vector here has one element per history (may vary); we take
  # the maximum as a conservative reference.
  lwd_ref <- max(lwd)
  spread  <- R_max * max(0.01, lwd_ref * 0.0025)
  vect_D  <- if (N == 1L) 0 else {
    indices <- seq_len(N) - (N + 1L) / 2
    indices * (spread / max(1L, N - 1L))
  }

  # Auxiliary function: exact angular compensation for parallelism
  # Deltatheta = asin(D / (R + D))  -- clamp avoids NaN when R+D ~ 0
  calc_dtheta <- function(R, D) {
    R_new <- R + D
    if (abs(R_new) < 1e-8) return(pi / 2)
    ratio <- max(-1, min(1, D / R_new))
    asin(ratio)
  }

  # -- History loop ----------------------------------------------------------
  for (i in seq_len(N)) {
    ec    <- color_list[[i]]
    lwd_i <- lwd[min(i, length(lwd))]
    dr    <- vect_D[i]

    # 1. Radial segments
    for (j in seq_len(nrow(edges))) {
      p <- edges[j, 1L];  h <- edges[j, 2L]

      R_p <- sqrt(xx[p]^2 + yy[p]^2)
      R_h <- sqrt(xx[h]^2 + yy[h]^2)

      R_p_new <- R_p + dr
      R_h_new <- R_h + dr

      # Exact angular compensation at each endpoint according to its local radius
      dtheta_p <- calc_dtheta(R_p, dr)
      dtheta_h <- calc_dtheta(R_h, dr)

      # Both endpoints use the child's angle; R handles negative radii natively
      theta_p_new <- node_angles[h] + dtheta_p
      theta_h_new <- node_angles[h] + dtheta_h

      segments(x0  = R_p_new * cos(theta_p_new),
               y0  = R_p_new * sin(theta_p_new),
               x1  = R_h_new * cos(theta_h_new),
               y1  = R_h_new * sin(theta_h_new),
               col = ec[j], lwd = lwd_i, lend = 1L)
    }

    # 2. Arcs
    for (m in seq_along(internal_nodes)) {
      nd         <- internal_nodes[m]
      ie         <- entry_idx[m]
      node_color <- if (!is.na(ie)) {
        ec[ie]
      } else {
        rc <- attr(ec, "root_color")
        if (!is.null(rc)) rc else ec[which(edges[, 1L] == nd)[1L]]
      }

      children_m <- children_of[[m]]
      R_nd       <- sqrt(xx[nd]^2 + yy[nd]^2)
      R_new      <- R_nd + dr

      dtheta_nd <- calc_dtheta(R_nd, dr)

      theta_adj_list <- numeric(length(children_m))
      for (k in seq_along(children_m)) {
        h <- children_m[k]
        theta_adj_list[k] <- node_angles[h] + dtheta_nd
      }

      ang_norm <- ifelse(theta_adj_list < 0, theta_adj_list + 2 * pi, theta_adj_list)
      ang_min  <- min(ang_norm)
      ang_max  <- max(ang_norm)

      if ((ang_max - ang_min) > pi) {
        ang_norm <- ifelse(ang_norm < pi, ang_norm + 2 * pi, ang_norm)
        ang_min  <- min(ang_norm)
        ang_max  <- max(ang_norm)
      }

      ang_seq <- seq(ang_min, ang_max, length.out = 100L)

      lines(x   = R_new * cos(ang_seq),
            y   = R_new * sin(ang_seq),
            col = node_color, lwd = lwd_i, lend = 1L)
    }
  }
}


# ==============================================================================
# MAIN FUNCTION
# ==============================================================================

#' Exports a multi-history evolutionary phylogenetic tree
#'
#' Overlays N evolutionary histories (branch color vectors) onto a
#' phylogeny, supporting three topology types and two output formats.
#' Applies the "Blank Canvas" strategy with isotropic diagonal offset
#' calculated on the fly from the actual physical parameters of the device.
#'
#' @section Internal workflow:
#' \enumerate{
#'   \item Exhaustive argument validation.
#'   \item Dynamic dimension calculation (height proportional to Ntip; fan = square).
#'   \item Graphics device opening (PNG @ 300 dpi, or PDF).
#'   \item Blank Canvas: \code{plot(tree, edge.color="transparent")}
#'         populates \code{.PlotPhyloEnv} without rendering edges.
#'   \item Extraction of nodal coordinates and physical canvas parameters.
#'   \item On-the-fly calculation of isotropic offsets (scale_x/scale_y).
#'   \item Topology-based rendering (phylogram / cladogram / fan).
#'   \item Guaranteed device closure via \code{finally}.
#' }
#'
#' @param tree          \code{\link[ape]{read.tree}} object of class \code{phylo}.
#' @param color_list    List of N character vectors (valid R/CSS colors).
#'                      Each vector must have exactly \code{nrow(tree$edge)}
#'                      elements (one per edge).
#' @param filename      Output filename \strong{without extension}.
#'                      The correct extension is added automatically.
#' @param type          Tree topology: \code{"phylogram"} (default),
#'                      \code{"cladogram"} or \code{"fan"}.
#' @param format        Output format: \code{"png"} (default) or \code{"pdf"}.
#' @param lwd           Branch width (default: \code{3}).
#' @param label_offset  Distance between the tip endpoints and their
#'                      labels, in data coordinate units
#'                      (default: \code{0.3}).
#' @param offset_range  Total amplitude of the Y offset range. The vector
#'                      \code{offsets_y} is built as
#'                      \code{seq(-range, range, length.out = N)}.
#'                      With \code{N = 1} the offset is always zero (default: \code{0.1}).
#' @param use_edge_length Logical. When \code{TRUE} (default) the tree is plotted
#'                      with its original branch lengths (if present). When
#'                      \code{FALSE} all edges are drawn with uniform length
#'                      regardless of \code{tree$edge.length}.
#' @param ladderize     Controls tip ordering before rendering. \code{FALSE}
#'                      (default) preserves the original order. \code{TRUE}
#'                      calls \code{ape::ladderize(right = FALSE)} (larger
#'                      subclade at bottom). \code{"right"} calls
#'                      \code{ape::ladderize(right = TRUE)} (larger subclade
#'                      at top). Must match the value used in
#'                      \code{execute_phylogeny()} to ensure export and screen
#'                      render are identical.
#' @param mar            Plot margins as vector \code{c(b, l, t, r)}
#'                       in lines (default: \code{c(1, 1, 1, 4)}).
#' @param width          Device width in inches. \code{NULL} = automatic (proportional to tip count).
#' @param height         Device height in inches. \code{NULL} = automatic.
#' @param legend_by_char Named list \code{character -> list(labels, colors)} for grouped legend.
#'                       When not \code{NULL}, takes precedence over \code{legend_labels} /
#'                       \code{legend_colors} / \code{legend_title}.
#' @param legend_labels  Flat character vector of legend labels (legacy mode).
#' @param legend_colors  Flat character vector of legend colors (legacy mode).
#' @param legend_corner  Corner for the legend: \code{"topleft"}, \code{"topright"},
#'                       \code{"bottomleft"} (default), or \code{"bottomright"}.
#' @param terminal_stretch Multiplier applied only to the terminal (tip) branches
#'                       when the tree has no real edge lengths, so tip labels and
#'                       colors get more visual space. Internal branches stay at
#'                       length 1 and the topology is preserved. Default \code{1}.
#' @param tip_end        Cladogram only: style of the free branch ends,
#'                       \code{"round"} (default), \code{"butt"} or \code{"square"}.
#'                       Only the tips are affected: at internal nodes the
#'                       colour bands always meet in closed miter joins.
#' @param branch_angle   Cladogram only. Angle of every branch from the horizontal,
#'                       in degrees (10-80; default \code{45}, so siblings meet at
#'                       right angles). 30-60 is recommended: smaller angles move
#'                       the lane corners away from the nodes. \code{NULL}
#'                       stretches the tree to the full canvas width.
#' @param show_labels    Logical. Draw the species names at the tips (default \code{TRUE}).
#' @param label_size     Multiplier of the automatic tip-label size (default \code{1}).
#' @param label_color    Color of the tip labels (default \code{"black"}).
#' @param legend_title   Single legend block title string (legacy mode). \code{NULL} = no title.
#' @param legend_layout  Grouped legend: \code{"vertical"} (default) stacks one block
#'                       per character; \code{"horizontal"} places them side by side.
#'                       Symbols, state names and titles are aligned in both.
#' @param overlay_fn     Optional function called after the legend, inside the open device,
#'                       with \code{par("usr")} already set. Receives \code{pp}, \code{cex_aj},
#'                       \code{label_offset_aj}, \code{R_tips}, and \code{gap_u}.
#' @param hide_fan_labels When \code{TRUE}, omits radial tip labels in fan topology,
#'                        delegating their drawing to \code{overlay_fn}.
#' @param header_labels  Optional character vector of column headers that
#'                        \code{overlay_fn} draws above the top tip (rotated 45
#'                        degrees, phylogram / cladogram). Space is reserved at
#'                        the top so they are not clipped.
#'
#' @return Invisible: full path of the generated file (string).
#'
#' @examples
#' \dontrun{
#' library(ape)
#' set.seed(42)
#' tree <- rtree(30)
#' n_e  <- nrow(tree$edge)
#'
#' ec1 <- sample(c("steelblue",  "tomato",        "gray70"), n_e, replace = TRUE)
#' ec2 <- sample(c("darkorange", "mediumseagreen", "gray70"), n_e, replace = TRUE)
#' ec3 <- sample(c("mediumpurple","gold3",         "gray70"), n_e, replace = TRUE)
#'
#' # Phylogram PNG
#' export_multimapr_tree(tree, list(ec1, ec2, ec3),
#'                       "my_tree_phylogram", type = "phylogram", format = "png")
#'
#' # Cladogram PDF
#' export_multimapr_tree(tree, list(ec1, ec2),
#'                       "my_tree_cladogram", type = "cladogram", format = "pdf")
#'
#' # Fan PNG
#' export_multimapr_tree(tree, list(ec1, ec2, ec3),
#'                       "my_tree_fan",       type = "fan",       format = "png")
#' }
#'
#' @export
export_multimapr_tree <- function(tree,
                                  color_list,
                                  filename,
                                  type          = "phylogram",
                                  format        = "png",
                                  lwd           = 3,
                                  label_offset  = 0.3,
                                  offset_range  = 0.1,
                                  use_edge_length = TRUE,
                                  ladderize       = FALSE,
                                  # terminal_stretch: multiplier applied ONLY to terminal
                                  # (tip) branches when the tree has no real edge lengths
                                  # (the usual case for character-mapping cladograms, where
                                  # every edge defaults to length 1). Internal branches are
                                  # left untouched so the topology's proportions don't change
                                  # -- this only gives tip labels/colors more visual room.
                                  # 1 = no change (default, backward compatible).
                                  terminal_stretch = 1,
                                  # tip_end: cladogram branch-end style,
                                  # "round" (default) | "butt" | "square".
                                  tip_end       = "round",
                                  # branch_angle: cladogram branch angle in degrees
                                  # (45 = right-angle joints); NULL stretches to the canvas.
                                  branch_angle  = 45,
                                  # -- TIP LABELS ----------------------------------------
                                  # show_labels: draw species names at the tips.
                                  # label_size : multiplier of the automatic label size.
                                  # label_color: color of the species names.
                                  show_labels   = TRUE,
                                  label_size    = 1,
                                  label_color   = "black",
                                  mar           = c(1, 1, 1, 4),
                                  # -- CUSTOM DIMENSION PARAMETERS -----------------------
                                  # NULL -> automatic dimensions are used, calculated
                                  # from the number of tips (original behavior).
                                  # A numeric value in inches activates dynamic scaling.
                                  width         = NULL,
                                  height        = NULL,
                                  # -- LEGEND PARAMETERS ---------------------------------
                                  # legend_by_char: named list character -> list(labels, colors).
                                  #   When provided, each character is drawn as an independent
                                  #   block with its name as a header. Takes precedence over
                                  #   legend_labels / legend_colors / legend_title.
                                  # legend_labels / legend_colors: parallel vectors with
                                  #   legend labels and colors (legacy mode). NULL = no legend.
                                  # legend_corner : corner where the legend is placed;
                                  #   one of "topleft" | "topright" | "bottomleft" | "bottomright".
                                  # legend_title  : legend block title (string or NULL, legacy only).
                                  legend_by_char = NULL,
                                  legend_labels  = NULL,
                                  legend_colors  = NULL,
                                  legend_corner  = "bottomleft",
                                  legend_title   = NULL,
                                  # legend_layout: per-character blocks stacked
                                  # ("vertical") or side by side ("horizontal").
                                  legend_layout  = "vertical",
                                  # -- OVERLAY CALLBACK ----------------------------------
                                  # Optional function executed AFTER the legend,
                                  # inside the same open device and with par("usr")
                                  # already established. Receives pp, cex_aj, label_offset_aj,
                                  # R_tips and gap_u to position elements with the
                                  # exact coordinates of the export render.
                                  overlay_fn      = NULL,
                                  # hide_fan_labels: when TRUE omits step 5b
                                  # (radial fan labels), delegating their drawing
                                  # to overlay_fn which repositions them further out.
                                  hide_fan_labels = FALSE,
                                  # header_labels: column headers the overlay draws above
                                  # the top tip (rotated 45 deg); room is reserved for them.
                                  header_labels   = NULL) {

  # -- 0. Argument validation ------------------------------------------------
  .emtree_validate_tree(tree)
  .emtree_validate_color_list(color_list, nrow(tree$edge))
  .emtree_lend_code(tip_end)   # validates tip_end
  branch_angle <- .emtree_branch_angle(branch_angle)
  .emtree_validate_labels(show_labels, label_size, label_color)
  legend_layout <- match.arg(legend_layout, c("vertical", "horizontal"))
  .emtree_validate_filename(filename)
  .emtree_validate_type_format(type, format)

  # -- 0b. Apply ladderize / edge length BEFORE any rendering ----------------
  # These transform the tree object used for ALL subsequent render steps so
  # that the exported file matches what execute_phylogeny() drew on screen.
  if (identical(ladderize, TRUE)) {
    tree <- ladderize(tree, right = FALSE)
  } else if (identical(ladderize, "right")) {
    tree <- ladderize(tree, right = TRUE)
  }
  if (!isTRUE(use_edge_length)) {
    tree$edge.length <- NULL
  }
  if (is.null(tree$edge.length)) {
    if (!is.numeric(terminal_stretch) || length(terminal_stretch) != 1L || terminal_stretch <= 0)
      stop("`terminal_stretch` must be a positive number.")
    if (terminal_stretch != 1) {
      el <- rep(1, nrow(tree$edge))
      el[tree$edge[, 2] <= Ntip(tree)] <- terminal_stretch
      tree$edge.length <- el
    }
  }

  if (!is.numeric(lwd)          || length(lwd) != 1L || lwd <= 0)
    stop("`lwd` must be a positive number.")
  if (!is.numeric(label_offset) || length(label_offset) != 1L)
    stop("`label_offset` must be a numeric scalar.")
  if (!is.numeric(offset_range) || length(offset_range) != 1L || offset_range < 0)
    stop("`offset_range` must be a non-negative number.")
  if (!is.numeric(mar) || length(mar) != 4L)
    stop("`mar` must be a numeric vector of length 4.")
  if (!is.null(width)  && (!is.numeric(width)  || length(width)  != 1L || width  <= 0))
    stop("`width` must be a positive number in inches, or NULL.")
  if (!is.null(height) && (!is.numeric(height) || length(height) != 1L || height <= 0))
    stop("`height` must be a positive number in inches, or NULL.")

  # Validate legend parameters
  # Priority: legend_by_char > legend_labels/legend_colors (legacy mode)
  has_grouped_legend <- !is.null(legend_by_char) && length(legend_by_char) > 0L
  has_flat_legend    <- !has_grouped_legend &&
    !is.null(legend_labels) && length(legend_labels) > 0L
  has_legend         <- has_grouped_legend || has_flat_legend

  if (has_grouped_legend) {
    # Only validate structure; space is reserved with expanded xlim/ylim
    # after measuring the legend with .emtree_measure_legend() (Phase 1 -> Phase 2).
    for (nm in names(legend_by_char)) {
      blk <- legend_by_char[[nm]]
      if (!is.list(blk) || is.null(blk$labels) || is.null(blk$colors))
        stop(sprintf(
          "`legend_by_char[[\"%s\"]]` must be a list with `labels` and `colors` elements.", nm))
      if (length(blk$labels) != length(blk$colors))
        stop(sprintf(
          "`legend_by_char[[\"%s\"]]`: `labels` and `colors` must have the same length.", nm))
    }

  } else if (has_flat_legend) {
    if (is.null(legend_colors) || length(legend_colors) != length(legend_labels))
      stop("`legend_colors` must have the same length as `legend_labels`.")
  }

  # -- 1. Dimension calculation and device opening ---------------------------
  n_tips <- Ntip(tree)
  N      <- length(color_list)

  # Create Exports/ folder if it does not exist and redirect output there
  dir.create("Exports", showWarnings = FALSE, recursive = TRUE)
  filename <- file.path("Exports", basename(filename))

  # Default dimensions (original behavior)
  height_default <- n_tips * 0.25 + 2   # 0.25 inches per tip + base margin
  width_default  <- .emtree_auto_width(type, height_default, branch_angle)

  # Fan -> square canvas to avoid polar distortions + 1.5 in extra
  # for more breathing room for labels and legend in radial mode.
  if (type == "fan") {
    height_default <- max(height_default, width_default)
    width_default  <- height_default
  }

  # Final dimensions: user > automatic
  height_in <- if (!is.null(height)) height else height_default
  width_in  <- if (!is.null(width))  width  else width_default

  # -- SCALE FACTORS ---------------------------------------------------------
  # Measure how much each axis grows (or shrinks) relative to the reference size.
  # scale_h > 1 -> taller canvas than default; scale_h < 1 -> more compact.
  # Used to: (a) adjust cex and label_offset proportionally,
  #          (b) compensate Y offset dilation (inverse scale).
  scale_h <- height_in / height_default
  scale_w <- width_in  / width_default

  # -- 2. Render parameters (cex, offsets, displacement) ---------------------
  # Calculated BEFORE opening the final device so they can also be used
  # in the off-screen probing device (Phase 1).
  cex_base <- max(1 / (1 + n_tips / 50), 0.2)
  cex_aj   <- cex_base * scale_h
  cex_lbl  <- cex_aj * label_size   # tip labels only

  # Offset proportional to the actual geometric scale of the data (2.5%)
  phy_tmp <- tree
  if (is.null(phy_tmp$edge.length)) phy_tmp$edge.length <- rep(1, nrow(phy_tmp$edge))
  max_depth <- max(node.depth.edgelength(phy_tmp))

  label_offset_aj <- max_depth * 0.025 * scale_w

  if (offset_range == 0.1 && N > 1L) {
    lwd_factor <- switch(type,
                         "phylogram" = 0.012, "cladogram" = 0.012, "fan" = 0.003, 0.012)
    offset_range <- max(0.05, lwd * lwd_factor)
  }
  adjusted_offset_range <- offset_range / scale_h
  lwd_vec  <- lwd * (0.95 ^ (seq_len(N) - 1L))
  cex_ley  <- max(0.5, min(1.0, 10 / (n_tips + 10)))

  mar_base <- c(1, 1, 1, 4)

  # -- Internal canvas render function --------------------------------------
  # Draws the "invisible" tree (edge.color = "transparent") to populate
  # .PlotPhyloEnv and establish par("usr"). Reused in the off-screen probe
  # and in the definitive render inside the final device.
  .render_canvas <- function(xlim_extra = NULL, ylim_extra = NULL) {
    if (type == "fan") par(mar = rep(2L, 4L)) else par(mar = mar_base)
    plot(tree,
         type           = type,
         edge.color     = "transparent",
         tip.color      = if (type == "fan") "transparent" else label_color,
         edge.width     = lwd,
         cex            = cex_lbl,
         label.offset   = label_offset_aj,
         no.margin      = if (type == "fan") FALSE else TRUE,
         show.tip.label = show_labels,
         font           = 3L,    # <--- Italic enabled
         x.lim          = xlim_extra,
         y.lim          = ylim_extra)
    .PlotPhyloEnv <- .emtree_get_PlotPhyloEnv()
    get("last_plot.phylo", envir = .PlotPhyloEnv)
  }

  # -- 3. PHASE 1 -- off-screen probing device --------------------------------
  # A temporary pdf(NULL) is opened (no file on disk) to:
  #   a) establish par("usr") with the same dimensions as the final device
  #   b) measure the legend with legend(..., plot = FALSE) in data coordinates
  # When closed with dev.off() the temporary device disappears without leaving a file.
  xlim_final <- NULL
  ylim_final <- NULL

  has_headers <- type != "fan" && length(header_labels) > 0L

  if (has_legend || has_headers ||
      (type == "fan" && hide_fan_labels && !is.null(overlay_fn))) {
    prev_dev <- dev.cur()
    pdf(NULL, width = width_in, height = height_in)   # off-screen, no file
    tryCatch({
      .render_canvas()   # only to establish par("usr")

      if (has_legend) {
        med <- .emtree_measure_legend(
          corner         = legend_corner,
          legend_by_char = if (has_grouped_legend) legend_by_char else NULL,
          legend_labels  = if (!has_grouped_legend) legend_labels else NULL,
          legend_title   = if (!has_grouped_legend) legend_title  else NULL,
          cex_ley        = cex_ley,
          legend_layout  = legend_layout
        )
      } else {
        med <- list(dx = 0, dy = 0, on_right = FALSE, going_down = FALSE)
      }

      usr <- par("usr")

      if (type == "fan") {
        expansion <- max(med$dx, med$dy)

        # If the overlay will draw rings (hide_fan_labels = TRUE), the viewport
        # must also accommodate those rings + species labels.
        # Geometry identical to .overlay_tip_figures / simple mode:
        #   R_max            = max(sqrt(pp$xx^2 + pp$yy^2)) of tips
        #   base_radius      = R_max * 1.06
        #   radius_increment = R_max * 0.10
        #   name_radius      = base_radius + (n_char-1)*increment + 0.5*increment
        #                    = R_max * (1.06 + 0.10*n_char)
        #   + estimated width of the longest label
        if (hide_fan_labels && !is.null(overlay_fn)) {
          pp_tmp        <- get("last_plot.phylo", envir = .emtree_get_PlotPhyloEnv())
          xx_t          <- pp_tmp$xx[seq_len(n_tips)]
          yy_t          <- pp_tmp$yy[seq_len(n_tips)]
          R_max_ov      <- max(sqrt(xx_t^2 + yy_t^2))
          n_char_ov     <- length(color_list)   # number of characters = number of histories
          name_radius   <- R_max_ov * (1.06 + 0.10 * n_char_ov)
          max_nc        <- max(nchar(tree$tip.label))
          lbl_est       <- if (show_labels) max_nc * R_max_ov * 0.018 * cex_lbl else 0
          total_radius  <- name_radius + lbl_est + R_max_ov * 0.05
          # The fan canvas is symmetric; we want total_radius to fit
          expansion_ov  <- max(0, total_radius - usr[2L])
          expansion     <- max(expansion, expansion_ov)
        }

        # Symmetric canvas for fan: expand equally in all directions
        # so overlay rings are not clipped in any quadrant.
        # If there is an additional legend in a corner, add only on that axis.
        xlim_final <- c(usr[1L] - expansion, usr[2L] + expansion)
        ylim_final <- c(usr[3L] - expansion, usr[4L] + expansion)
        # Add extra space on the legend side if applicable
        if (has_legend && med$dx > 0) {
          if (med$on_right) xlim_final[2L] <- xlim_final[2L] + med$dx
          else              xlim_final[1L] <- xlim_final[1L] - med$dx
        }
      } else {
        if (has_legend) {
          xlim_final <- if (med$on_right) c(usr[1L], usr[2L] + med$dx)
          else              c(usr[1L] - med$dx, usr[2L])
          ylim_final <- if (med$going_down) c(usr[3L], usr[4L] + med$dy)
          else               c(usr[3L] - med$dy, usr[4L])
        }
        if (has_headers) {
          pp_h <- get("last_plot.phylo", envir = .emtree_get_PlotPhyloEnv())
          extra <- .emtree_header_space(header_labels, cex_aj * 0.9,
                                        pp_h$yy[seq_len(n_tips)])
          if (is.null(ylim_final)) ylim_final <- usr[3:4]
          ylim_final[2L] <- ylim_final[2L] + extra
        }
      }
    }, finally = .close_device_restore(prev_dev))   # always close; leaves no file on disk
  }

  # -- 4. Open final device and render --------------------------------------
  prev_dev <- dev.cur()
  .emtree_open_device(filename, format, width_in, height_in)
  output_path <- paste0(filename, ".", format)

  tryCatch({

    # Definitive canvas: with expanded xlim/ylim if there is a legend, or normal otherwise.
    pp      <- .render_canvas(xlim_extra = xlim_final, ylim_extra = ylim_final)
    offsets <- .emtree_calc_offsets(N, adjusted_offset_range)

    # -- 5. Branch rendering by topology --------------------------------------
    if (type == "phylogram") {
      .emtree_render_phylogram(pp, tree, color_list, lwd_vec, offsets)

    } else if (type == "cladogram") {
      .emtree_render_cladogram(pp, tree, color_list, lwd_vec, offsets,
                               tip_end = tip_end, branch_angle = branch_angle)

    } else {
      .emtree_render_fan(pp, tree, color_list, lwd_vec, offsets)

      # -- 5b. Fan tip labels -- exact radius --------------------------------
      xx_tips <- pp$xx[seq_len(n_tips)]
      yy_tips <- pp$yy[seq_len(n_tips)]
      angles  <- atan2(yy_tips, xx_tips)
      R_tips  <- sqrt(xx_tips^2 + yy_tips^2)
      gap_u   <- strwidth("m", cex = cex_aj) * 1.5

      if (!hide_fan_labels && show_labels) {
        old_xpd <- par("xpd")
        par(xpd = NA)
        for (j in seq_len(n_tips)) {
          ang_j  <- angles[j]
          R_j    <- R_tips[j] + gap_u
          right_side <- cos(ang_j) >= 0
          srt_j  <- ang_j * 180 / pi
          adj_j  <- if (right_side) c(0, 0.5) else { srt_j <- srt_j + 180; c(1, 0.5) }
          text(R_j * cos(ang_j),
               R_j * sin(ang_j),
               labels = tree$tip.label[j],
               adj    = adj_j,
               cex    = cex_lbl,
               col    = label_color,
               font   = 3L,        # <--- Italic enabled
               srt    = srt_j)
        }
        par(xpd = old_xpd)
      }
    }

    # -- 6. Legend -- drawn inside the expanded viewport --------------------
    if (has_legend) {
      .emtree_draw_legend(
        corner         = legend_corner,
        legend_by_char = if (has_grouped_legend) legend_by_char else NULL,
        legend_labels  = legend_labels,
        legend_colors  = legend_colors,
        legend_title   = legend_title,
        cex_ley        = cex_ley,
        legend_layout  = legend_layout
      )
    }

    # -- 7. Optional overlay -- additional figures on top of the render -----
    if (!is.null(overlay_fn)) {
      # Calculate R_tips and gap_u for fan (NULL for other topologies)
      if (type == "fan") {
        .ov_xx    <- pp$xx[seq_len(n_tips)]
        .ov_yy    <- pp$yy[seq_len(n_tips)]
        .ov_Rtips <- sqrt(.ov_xx^2 + .ov_yy^2)
        .ov_gapu  <- strwidth("m", cex = cex_aj) * 1.5
      } else {
        .ov_Rtips <- NULL
        .ov_gapu  <- NULL
      }
      overlay_fn(pp              = pp,
                 cex_aj          = cex_aj,
                 label_offset_aj = label_offset_aj,
                 R_tips          = .ov_Rtips,
                 gap_u           = .ov_gapu)
    }

    # Indicate whether dimensions are automatic or user-defined
    dims_origin <- if (!is.null(width) || !is.null(height)) " [custom]" else " [auto]"
    message(sprintf(
      "[export_multimapr_tree] File generated: %s  (%d tips, %d histories, %.1f \u00d7 %.1f in%s%s)",
      output_path, n_tips, N, width_in, height_in,
      if (format == "png") " @ 300 dpi" else "",
      dims_origin))

  }, error = function(e) {
    stop(sprintf("[export_multimapr_tree] Error during rendering (%s, %s): %s",
                 type, format, conditionMessage(e)))
  }, finally = {
    .close_device_restore(prev_dev)   # Close the device in all cases
  })

  invisible(output_path)
}


# ==============================================================================
# USAGE EXAMPLES  (change FALSE -> TRUE to run)
# ==============================================================================
if (FALSE) {

  library(ape)
  set.seed(42)

  tree <- rtree(40)
  n_e  <- nrow(tree$edge)

  palette1 <- c("steelblue",   "tomato",         "gray70")
  palette2 <- c("darkorange",  "mediumseagreen",  "gray70")
  palette3 <- c("mediumpurple","gold3",           "gray70")

  ec1 <- sample(palette1, n_e, replace = TRUE)
  ec2 <- sample(palette2, n_e, replace = TRUE)
  ec3 <- sample(palette3, n_e, replace = TRUE)

  # -- Phylogram with legend in bottom-left corner ----------------------------
  export_multimapr_tree(
    tree          = tree,
    color_list    = list(ec1, ec2, ec3),
    filename      = "test_phylogram",
    type          = "phylogram",
    format        = "png",
    legend_labels = c("State A", "State B", "Ambiguous"),
    legend_colors = palette1,
    legend_corner = "bottomleft",
    legend_title  = "Character 1"
  )

  # -- Cladogram with legend in top-right corner ------------------------------
  export_multimapr_tree(
    tree          = tree,
    color_list    = list(ec1, ec2),
    filename      = "test_cladogram",
    type          = "cladogram",
    format        = "pdf",
    legend_labels = c("State A", "State B", "Ambiguous"),
    legend_colors = palette1,
    legend_corner = "topright",
    legend_title  = "Character 1"
  )

  # -- Fan with legend in top-left corner ------------------------------------
  export_multimapr_tree(
    tree          = tree,
    color_list    = list(ec1, ec2, ec3),
    filename      = "test_fan",
    type          = "fan",
    format        = "png",
    legend_labels = c("State A", "State B", "Ambiguous"),
    legend_colors = palette1,
    legend_corner = "topleft",
    legend_title  = "Character 1"
  )

  # -- Single history without legend (original behavior) ---------------------
  export_multimapr_tree(
    tree       = tree,
    color_list = list(ec1),
    filename   = "test_single",
    type       = "phylogram",
    format     = "pdf"
  )
}

# ==============================================================================
# SCREEN RENDERING -- ADVANCED ENGINE
# ==============================================================================

#' Renders the multi-history tree directly in R's interactive window
#'
#' Runs exactly the same Blank Canvas logic used by
#' \code{export_multimapr_tree}, but operates on the active R graphics window
#' instead of opening storage devices on disk.
#' Supports all three topologies: \code{"phylogram"}, \code{"cladogram"} and
#' \code{"fan"}.
#'
#' @param tree          \code{\link[ape]{read.tree}} object of class \code{phylo}.
#' @param color_list    List of N character vectors (valid R/CSS colors).
#'                      Each vector must have exactly \code{nrow(tree$edge)}
#'                      elements (one per edge).
#' @param type          Tree topology: \code{"phylogram"} (default),
#'                      \code{"cladogram"} or \code{"fan"}.
#' @param lwd           Branch width (default: \code{3}).
#' @param label_offset  Distance between tip endpoints and their
#'                      labels (default: \code{0.3}).
#' @param offset_range  Total amplitude of the Y offset range (default: \code{0.1}).
#' @param use_edge_length Logical. When \code{TRUE} (default) branch lengths are
#'                      used if present. When \code{FALSE} all branches are drawn
#'                      with uniform length.
#' @param ladderize     Controls tip ordering. \code{FALSE} (default) preserves
#'                      original order. \code{TRUE} calls
#'                      \code{ape::ladderize(right = FALSE)}. \code{"right"}
#'                      calls \code{ape::ladderize(right = TRUE)}.
#' @param legend_by_char Named list character -> list(labels, colors) for
#'                       grouped legend. If not NULL takes precedence over
#'                       \code{legend_labels} / \code{legend_colors}.
#' @param legend_labels  Flat label vector (legacy mode).
#' @param legend_colors  Flat color vector (legacy mode).
#' @param legend_corner  Legend corner: \code{"topleft"}, \code{"topright"},
#'                       \code{"bottomleft"} (default) or \code{"bottomright"}.
#' @param terminal_stretch Multiplier applied only to the terminal (tip) branches
#'                        when the tree has no real edge lengths. Internal
#'                        branches stay at length 1. Default \code{1}.
#' @param tip_end         Cladogram only: \code{"round"} (default), \code{"butt"}
#'                        or \code{"square"} branch ends.
#' @param branch_angle    Cladogram only: branch angle in degrees (10-80, default
#'                        \code{45}); \code{NULL} stretches to the full width.
#' @param show_labels     Logical. Draw the species names at the tips (default \code{TRUE}).
#' @param label_size      Multiplier of the automatic tip-label size (default \code{1}).
#' @param label_color     Color of the tip labels (default \code{"black"}).
#' @param legend_title    Single legend block title (legacy mode).
#' @param legend_layout   Grouped legend: \code{"vertical"} (default) or
#'                        \code{"horizontal"} blocks.
#' @param overlay_fn      Optional function called after the legend, inside the
#'                        active device with \code{par("usr")} already set.
#' @param hide_fan_labels When \code{TRUE}, omits radial tip labels in fan topology,
#'                        delegating their drawing to \code{overlay_fn}.
#' @param header_labels   Optional column headers drawn by \code{overlay_fn} above
#'                        the top tip; space is reserved for them.
#' @return Invisible NULL. Draws on the active graphics device.
plot_multimapr_screen <- function(tree,
                                  color_list,
                                  type          = "phylogram",
                                  lwd           = 3,
                                  label_offset  = 0.3,
                                  offset_range  = 0.1,
                                  use_edge_length = TRUE,
                                  ladderize       = FALSE,
                                  terminal_stretch = 1,
                                  tip_end        = "round",
                                  branch_angle   = 45,
                                  show_labels    = TRUE,
                                  label_size     = 1,
                                  label_color    = "black",
                                  legend_by_char = NULL,
                                  legend_labels  = NULL,
                                  legend_colors  = NULL,
                                  legend_corner  = "bottomleft",
                                  legend_title   = NULL,
                                  legend_layout  = "vertical",
                                  # -- OVERLAY CALLBACK ----------------------------------
                                  # Optional function executed AFTER the legend,
                                  # inside the same open device and with par("usr")
                                  # already established. Receives pp, cex_aj, label_offset_aj,
                                  # R_tips and gap_u to position elements with the
                                  # exact coordinates of the screen render.
                                  overlay_fn      = NULL,
                                  # hide_fan_labels: when TRUE omits the radial fan
                                  # labels, delegating their drawing to overlay_fn.
                                  hide_fan_labels = FALSE,
                                  header_labels   = NULL) {

  .emtree_validate_tree(tree)
  .emtree_validate_color_list(color_list, nrow(tree$edge))
  .emtree_lend_code(tip_end)   # validates tip_end
  branch_angle <- .emtree_branch_angle(branch_angle)
  .emtree_validate_labels(show_labels, label_size, label_color)
  legend_layout <- match.arg(legend_layout, c("vertical", "horizontal"))

  # Apply ladderize / edge length before rendering (mirrors export behavior)
  if (identical(ladderize, TRUE)) {
    tree <- ladderize(tree, right = FALSE)
  } else if (identical(ladderize, "right")) {
    tree <- ladderize(tree, right = TRUE)
  }
  if (!isTRUE(use_edge_length)) {
    tree$edge.length <- NULL
  }
  if (is.null(tree$edge.length)) {
    if (!is.numeric(terminal_stretch) || length(terminal_stretch) != 1L || terminal_stretch <= 0)
      stop("`terminal_stretch` must be a positive number.")
    if (terminal_stretch != 1) {
      el <- rep(1, nrow(tree$edge))
      el[tree$edge[, 2] <= Ntip(tree)] <- terminal_stretch
      tree$edge.length <- el
    }
  }

  n_tips <- Ntip(tree)
  N      <- length(color_list)

  # Ensure a graphics window is open
  if (is.null(dev.list())) dev.new()

  # Neutral scale factors for interactive screen
  scale_h <- 1
  scale_w <- 1

  cex_base <- max(1 / (1 + n_tips / 50), 0.2)
  cex_aj   <- cex_base * scale_h
  cex_lbl  <- cex_aj * label_size   # tip labels only

  # Offset proportional to actual tree scale (2.5%)
  phy_tmp <- tree
  if (is.null(phy_tmp$edge.length)) phy_tmp$edge.length <- rep(1, nrow(phy_tmp$edge))
  max_depth <- max(node.depth.edgelength(phy_tmp))
  label_offset_aj <- max_depth * 0.025 * scale_w

  adjusted_offset_range <- offset_range / scale_h
  lwd_vec <- lwd * (0.95 ^ (seq_len(N) - 1L))
  cex_ley <- max(0.5, min(1.0, 10 / (n_tips + 10)))

  mar_base <- c(1, 1, 1, 4)

  # Internal function: Blank Canvas to populate .PlotPhyloEnv
  .render_canvas_screen <- function(ylim_extra = NULL) {
    if (type == "fan") par(mar = rep(2L, 4L)) else par(mar = mar_base)
    plot(tree,
         type           = type,
         edge.color     = "transparent",
         tip.color      = if (type == "fan") "transparent" else label_color,
         edge.width     = lwd,
         cex            = cex_lbl,
         label.offset   = label_offset_aj,
         no.margin      = if (type == "fan") FALSE else TRUE,
         show.tip.label = show_labels,
         font           = 3L,
         y.lim          = ylim_extra)
    pp <- get("last_plot.phylo", envir = get(".PlotPhyloEnv", envir = asNamespace("ape")))
    return(pp)
  }

  # Clear window and render invisible canvas
  plot.new()
  pp      <- .render_canvas_screen()
  if (type != "fan" && length(header_labels) > 0L) {
    extra <- .emtree_header_space(header_labels, cex_aj * 0.9, pp$yy[seq_len(n_tips)])
    if (extra > 0) {
      usr <- par("usr")
      plot.new()
      pp <- .render_canvas_screen(ylim_extra = c(usr[3L], usr[4L] + extra))
    }
  }
  offsets <- .emtree_calc_offsets(N, adjusted_offset_range)

  # Dispatch to topology-specific geometric renderer
  if (type == "phylogram") {
    .emtree_render_phylogram(pp, tree, color_list, lwd_vec, offsets)

  } else if (type == "cladogram") {
    .emtree_render_cladogram(pp, tree, color_list, lwd_vec, offsets,
                               tip_end = tip_end, branch_angle = branch_angle)

  } else if (type == "fan") {
    .emtree_render_fan(pp, tree, color_list, lwd_vec, offsets)

    # Manual italic radial labels in fan mode
    xx_tips <- pp$xx[seq_len(n_tips)]
    yy_tips <- pp$yy[seq_len(n_tips)]
    angles  <- atan2(yy_tips, xx_tips)
    R_tips  <- sqrt(xx_tips^2 + yy_tips^2)
    gap_u   <- strwidth("m", cex = cex_aj) * 1.5

    if (!hide_fan_labels && show_labels) {
      old_xpd <- par("xpd"); par(xpd = NA)
      for (j in seq_len(n_tips)) {
        ang_j      <- angles[j]
        R_j        <- R_tips[j] + gap_u
        right_side <- cos(ang_j) >= 0
        srt_j      <- ang_j * 180 / pi
        adj_j      <- if (right_side) c(0, 0.5) else { srt_j <- srt_j + 180; c(1, 0.5) }
        text(R_j * cos(ang_j), R_j * sin(ang_j),
             labels = tree$tip.label[j], adj = adj_j, cex = cex_lbl,
             col = label_color, font = 3L, srt = srt_j)
      }
      par(xpd = old_xpd)
    }
  }

  # Draw legend if applicable
  has_grouped_legend <- !is.null(legend_by_char) && length(legend_by_char) > 0L
  has_flat_legend    <- !has_grouped_legend && !is.null(legend_labels) && length(legend_labels) > 0L

  if (has_grouped_legend || has_flat_legend) {
    .emtree_draw_legend(
      corner         = legend_corner,
      legend_by_char = legend_by_char,
      legend_labels  = legend_labels,
      legend_colors  = legend_colors,
      legend_title   = legend_title,
      cex_ley        = cex_ley,
      legend_layout  = legend_layout
    )
  }

  # -- Optional overlay -- additional figures over the screen render --------
  if (!is.null(overlay_fn)) {
    if (type == "fan") {
      .ov_xx    <- pp$xx[seq_len(n_tips)]
      .ov_yy    <- pp$yy[seq_len(n_tips)]
      .ov_Rtips <- sqrt(.ov_xx^2 + .ov_yy^2)
      .ov_gapu  <- strwidth("m", cex = cex_aj) * 1.5
    } else {
      .ov_Rtips <- NULL
      .ov_gapu  <- NULL
    }
    overlay_fn(pp              = pp,
               cex_aj          = cex_aj,
               label_offset_aj = label_offset_aj,
               R_tips          = .ov_Rtips,
               gap_u           = .ov_gapu)
  }

  invisible(NULL)
}
