################################################################################
# MultiMapR — Core Rendering Controllers
# Rendering and evolutionary state assignment module
################################################################################


# ==============================================================================
# BLOCK 2 — SIMPLE MAPPING
# ==============================================================================

#' Horizontal distance between the columns of terminal figures
#'
#' At least the historical 2.5 "M" widths, and never less than the figure's
#' own diameter plus a margin, so large figures (\code{figure_size}) do not
#' touch the next column. Figures with \code{pch} 21-25 have a diameter of
#' about \code{0.375 * pt.cex} character heights.
#'
#' @param cex_txt  Text size used for the table.
#' @param tam_fig  Figure size (\code{pt.cex}).
#' @return Separation in x data units of the current plot.
#' @keywords internal
.figure_col_sep <- function(cex_txt, tam_fig) {
  x_per_in <- diff(par("usr")[1:2]) / par("pin")[1L]
  diam_in  <- 0.375 * tam_fig * par("cin")[2L]
  max(strwidth("M", cex = cex_txt) * 2.5, diam_in * 1.35 * x_per_in)
}

#' Resolves the tip figure color given its state value
#'
#' @param value          Character state value.
#' @param color_map      Named vector: state -> color.
#' @param color_na       Fallback color for NA or unrecognized states.
#' @return A single color string.
resolve_tip_color <- function(value, color_map, color_na = "gray70") {
  value <- as.character(value)
  if (is.na(value) || value == "") return(color_na)
  if (value %in% names(color_map)) return(color_map[[value]])
  return(color_na)
}


#' SIMPLE MAPPING — neutral tree + colored figure table/rings at terminals
#'
#' For phylogram and cladogram:
#'   - Draws the tree with x.lim expanded to accommodate the table and legend.
#'   - The colored figure table appears on the right with the full character
#'     name as a column header (rotated 45°).
#'   - The legend is positioned in the chosen corner, one block per character.
#' For fan:
#'   - Concentric rings of symbols around the tree.
#'   - C1, C2... abbreviations above each ring (once, anchored at the
#'     terminal with the largest angle).
#'   - Tip labels attached to the last ring, parallel to the radius.
#'   - Legend in blocks, same as multi-mapping mode.
#'
#' Spacing strategy (phylogram / cladogram):
#'   Phase 1 — off-screen probe: the tree is drawn in pdf(NULL) with the same
#'   dimensions as the final file to obtain real par("usr") and par("pin").
#'   With those values the following are calculated:
#'     a) x_max_tabla  = x of the last column + horizontal projection of the
#'                       longest header rotated 45°
#'     b) x_max_leyenda = x_max_tabla + legend width measured with
#'                        legend(..., plot=FALSE)
#'   If x_max_leyenda > `usr[2]`, x.lim is expanded by the necessary amount
#'   (in data units) and passed as x.lim to the definitive plot().
#'   Phase 2 — definitive render: the final device is opened and everything
#'   is drawn inside the already-expanded viewport.
#'
#' @param filogenia  phylo object.
#' @param config     Configuration list from setup_mapping_config().
#' @return Invisible NULL. Draws on the active graphics device.
plot_simple_mapping <- function(filogenia, config) {

  datos_ord       <- config$aligned_data
  caracteres      <- config$characters
  colores_por_car <- config$colors_by_character
  tipo_arbol      <- config$tree_type
  pch_fig         <- config$figure_pch
  tam_fig         <- config$figure_size
  simple_mode     <- config$simple_mode %||% "figures"   # "figures" | "tip_color"
  cex_aj          <- adjust_cex(filogenia, config = config)
  n_tips          <- Ntip(filogenia)
  n_car           <- length(caracteres)
  show_lbl        <- config$show_labels %||% TRUE
  cex_lbl         <- cex_aj * (config$label_size %||% 1)
  col_lbl         <- config$label_color %||% "black"
  show_leg        <- config$show_legend %||% TRUE

  resolve_color <- function(valor, colores_estado, car) {
    .ancestral_tip_color(valor, colores_estado, config, car)
  }

  exportar  <- !is.null(config$export_filename)
  fn_export <- config$export_filename

  cat("\n=== Generating plot (Simple Mapping) ===\n")
  if (exportar) cat("    \u2192 Exporting to:",
                    paste0(fn_export, ".", config$export_format %||% "png"), "\n")

  # ── Helper: draws the legend in blocks (same as plot_superimposed_characters) ──
  .draw_simple_legend <- function(pos_ley) {
    if (!show_leg) return(invisible(NULL))
    blocks <- .simple_legend_blocks()
    lay    <- .legend_layout(blocks, pos_ley, config$legend_layout %||% "vertical",
                             cex = cex_aj, pt_cex = tam_fig, pch = pch_fig)
    old_xpd <- par("xpd"); par(xpd = TRUE)
    .legend_draw_blocks(lay, blocks, cex = cex_aj, pt_cex = tam_fig, pch = pch_fig,
                        filled = pch_fig %in% 21:25)
    par(xpd = old_xpd)
    invisible(NULL)
  }

  # ── Helper: one legend block per character (fan titles carry the ring id) ────
  .simple_legend_blocks <- function() {
    lapply(seq_along(caracteres), function(i) {
      car   <- caracteres[i]
      col_e <- c(colores_por_car[[car]], .missing_legend_entries(config, car))
      list(title  = if (tipo_arbol == "fan") paste0("C", i, " (", car, ")") else car,
           labels = names(col_e), colors = unname(col_e))
    })
  }

  # ── Main render function ─────────────────────────────────────────────────────
  # xlim_extra: if not NULL, passed as x.lim to plot() to expand the viewport
  # to the right to accommodate table + legend.
  .render_simple_mapping <- function(xlim_extra = NULL, ylim_extra = NULL) {

    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    par(lend = config$tip_end %||% "round")

    pos_leyenda <- if (!is.null(config$legend_corner)) config$legend_corner
    else if (tipo_arbol == "fan") "topleft" else "topright"

    # ==========================================================================
    # FAN
    # Strategy: no.margin = TRUE + symmetric x.lim calculated post-plot.
    # 1. First plot with tip.color = "transparent" and without x.lim to populate
    #    last_plot.phylo with the real tree coordinates.
    # 2. Calculate radio_contenido = tree + rings + label gap.
    # 3. Estimate width of longest labels in data coordinates.
    # 4. Replot with x.lim = c(-radio_total, radio_total) so that all
    #    content fits inside the viewport without relying on margins.
    # ==========================================================================
    if (tipo_arbol == "fan") {

      # ── Step 1: tree coordinates — calculated off-screen before opening ──────
      # the final device (see PHASE 1 fan block below).
      # Here we only read the already-calculated results.
      obj <- .fan_coords_probe

      xx_tip <- obj$xx[seq_len(n_tips)]
      yy_tip <- obj$yy[seq_len(n_tips)]

      angulos        <- atan2(yy_tip, xx_tip)
      tip_radius     <- sqrt(xx_tip^2 + yy_tip^2)
      max_tip_radius <- max(tip_radius)

      if (simple_mode == "tip_color") {
        # ── MODE tip_color (fan): colored terminal branches, no rings ──────────
        # Edges whose child is a tip are painted with that tip's state color;
        # all other edges stay gray.
        char1     <- caracteres[1L]
        col_tips1 <- sapply(datos_ord[[char1]], resolve_color,
                            colores_estado = colores_por_car[[char1]], car = char1)

        edge_col <- rep("gray70", nrow(filogenia$edge))
        for (e in seq_len(nrow(filogenia$edge))) {
          child <- filogenia$edge[e, 2L]
          if (child <= n_tips)
            edge_col[e] <- col_tips1[child]
        }

        # Label geometry: identical to the multimapping fan (export_multimapr_tree)
        # so that the gap between tips and species names is consistent.
        # gap_u = strwidth("m") * 1.5  (evaluated after the plot establishes usr)
        # radio_nombres = max_tip_radius + gap_u
        # radio_total   = radio_nombres + estimated label width + small padding

        # First plot with show.tip.label = FALSE to get definitive tip coordinates
        # and then measure strwidth in the correct coordinate system.
        max_nchar_lbl <- max(nchar(filogenia$tip.label))

        # Rough xlim to open the canvas (will be tightened after measuring)
        rough_lbl_w <- max_tip_radius * 0.018 * cex_aj * max_nchar_lbl
        rough_total <- max_tip_radius * 1.06 + rough_lbl_w + max_tip_radius * 0.05
        xlim_fan    <- c(-rough_total, rough_total)

        plot(filogenia,
             type           = "fan",
             cex            = cex_aj,
             label.offset   = 0,           # labels drawn manually below
             edge.width     = config$branch_width,
             edge.color     = edge_col,
             tip.color      = "transparent",
             show.tip.label = FALSE,
             no.margin      = TRUE,
             font           = 3L,
             x.lim          = xlim_fan)

        # Refresh coordinates after definitive plot
        obj2       <- get("last_plot.phylo", envir = .PlotPhyloEnv)
        xx_tip2    <- obj2$xx[seq_len(n_tips)]
        yy_tip2    <- obj2$yy[seq_len(n_tips)]
        angulos2   <- atan2(yy_tip2, xx_tip2)
        tip_radius2 <- sqrt(xx_tip2^2 + yy_tip2^2)
        R_tips2     <- max(tip_radius2)

        gap_u2        <- strwidth("m", cex = cex_aj) * 1.5
        radio_nombres2 <- R_tips2 + gap_u2

        # Draw italic tip labels radially
        old_xpd <- par("xpd"); par(xpd = TRUE)
        for (j in seq_len(n_tips)[show_lbl]) {
          ang_j  <- angulos2[j]; deg_j <- ang_j * 180 / pi
          lado_d <- cos(ang_j) >= 0
          srt_j  <- if (lado_d) deg_j else deg_j + 180
          adj_j  <- if (lado_d) c(0, 0.5) else c(1, 0.5)
          text(radio_nombres2 * cos(ang_j), radio_nombres2 * sin(ang_j),
               labels = filogenia$tip.label[j], adj = adj_j, cex = cex_lbl,
               col = col_lbl, srt = srt_j, font = 3L)
        }
        par(xpd = old_xpd)

        .draw_simple_legend(pos_leyenda)

      } else {
        # ── MODE figures (fan): concentric symbol rings ─────────────────────────

        # ── Step 2: ring geometry ───────────────────────────────────────────────
        radio_base       <- max_tip_radius * 1.06
        incremento_radio <- max_tip_radius * 0.10
        radio_ultimo     <- radio_base + (n_car - 1L) * incremento_radio

        # Space between last ring and label
        gap_etiq         <- incremento_radio * 0.5
        radio_nombres    <- radio_ultimo + gap_etiq

        # ── Step 3: estimate total radius including tip labels ──────────────────
        max_nchar_lbl <- max(nchar(filogenia$tip.label))
        char_w_u      <- max_tip_radius * 0.018 * cex_aj
        lbl_w_u       <- max_nchar_lbl * char_w_u
        radio_total   <- radio_nombres + lbl_w_u + max_tip_radius * 0.05

        # ── Step 4: replot with symmetric x.lim that contains everything ────────
        xlim_fan <- c(-radio_total, radio_total)

        plot(filogenia,
             type           = "fan",
             cex            = cex_aj,
             label.offset   = config$label_offset,
             edge.width     = config$branch_width,
             tip.color      = "transparent",
             show.tip.label = FALSE,
             no.margin      = TRUE,
             x.lim          = xlim_fan)

        # Refresh coordinates from the definitive plot
        obj    <- get("last_plot.phylo", envir = .PlotPhyloEnv)
        xx_tip <- obj$xx[seq_len(n_tips)]
        yy_tip <- obj$yy[seq_len(n_tips)]
        angulos <- atan2(yy_tip, xx_tip)

        # ── Symbol rings ──────────────────────────────────────────────────────────
        old_xpd <- par("xpd"); par(xpd = TRUE)
        for (i in seq_len(n_car)) {
          col_tips     <- sapply(datos_ord[[caracteres[i]]], resolve_color,
                                 colores_estado = colores_por_car[[caracteres[i]]],
                                 car = caracteres[i])
          radio_actual <- radio_base + (i - 1L) * incremento_radio

          points(radio_actual * cos(angulos),
                 radio_actual * sin(angulos),
                 pch = pch_fig, bg = col_tips, col = "black", cex = tam_fig)

          gap_tan <- incremento_radio * 0.55
          for (j in seq_len(n_tips)) {
            ang_j   <- angulos[j]
            deg_j   <- ang_j * 180 / pi
            lado_d  <- cos(ang_j) >= 0
            srt_j   <- if (lado_d) deg_j else deg_j + 180
            tan_x   <- -sin(ang_j)
            tan_y   <-  cos(ang_j)
            cx      <- radio_actual * cos(ang_j) + tan_x * gap_tan
            cy      <- radio_actual * sin(ang_j) + tan_y * gap_tan
            text(cx, cy,
                 labels = paste0("C", i),
                 adj    = c(0.5, 0.5),
                 cex    = cex_aj * 0.65,
                 font   = 2L,
                 srt    = srt_j)
          }
        }
        par(xpd = old_xpd)

        # ── Tip labels — parallel to radius ────────────────────────────────────
        old_xpd <- par("xpd"); par(xpd = TRUE)
        for (j in seq_len(n_tips)[show_lbl]) {
          ang_j  <- angulos[j]; deg_j <- ang_j * 180 / pi
          lado_d <- cos(ang_j) >= 0
          srt_j  <- if (lado_d) deg_j else deg_j + 180
          adj_j  <- if (lado_d) c(0, 0.5) else c(1, 0.5)
          text(radio_nombres * cos(ang_j), radio_nombres * sin(ang_j),
               labels = filogenia$tip.label[j], adj = adj_j, cex = cex_lbl, srt = srt_j,
               col    = col_lbl, font   = 3L)
        }
        par(xpd = old_xpd)

        .draw_simple_legend(pos_leyenda)

      }   # end if simple_mode

      # ==========================================================================
      # PHYLOGRAM / CLADOGRAM
      # ==========================================================================
    } else {

      par(mar = c(1, 1, 2, 1), xpd = TRUE)

      if (simple_mode == "tip_color") {
        # ── MODE: colored terminal branches — one character only ───────────────
        # Edges whose child is a tip are painted with that tip's state color;
        # all other edges (internal) stay gray.
        char1     <- caracteres[1L]
        col_tips1 <- sapply(datos_ord[[char1]], resolve_color,
                            colores_estado = colores_por_car[[char1]], car = char1)

        edge_col <- rep("gray70", nrow(filogenia$edge))
        for (e in seq_len(nrow(filogenia$edge))) {
          child <- filogenia$edge[e, 2L]
          if (child <= n_tips)
            edge_col[e] <- col_tips1[child]
        }

        plot(filogenia,
             type           = tipo_arbol,
             cex            = cex_lbl,
             label.offset   = config$label_offset,
             edge.width     = config$branch_width,
             edge.color     = edge_col,
             tip.color      = col_lbl,
             show.tip.label = show_lbl,
             font           = 3L,
             x.lim          = xlim_extra,
             y.lim          = ylim_extra)

        .draw_simple_legend(pos_leyenda)

      } else {
        # ── MODE: figures (default) ────────────────────────────────────────────
        plot(filogenia,
             type           = tipo_arbol,
             cex            = cex_lbl,
             label.offset   = config$label_offset,
             edge.width     = config$branch_width,
             tip.color      = col_lbl,
             show.tip.label = show_lbl,
             font           = 3L,
             x.lim          = xlim_extra,   # expanded for the table (+ legend)
             y.lim          = ylim_extra)   # expanded for a horizontal legend

        obj    <- get("last_plot.phylo", envir = .PlotPhyloEnv)
        xx_tip <- obj$xx[seq_len(n_tips)]
        yy_tip <- obj$yy[seq_len(n_tips)]
        usr    <- par("usr")
        pin    <- par("pin")

        sep_col    <- .figure_col_sep(cex_aj, tam_fig)
        max_tip_x  <- max(xx_tip)
        max_lbl_w  <- if (show_lbl) max(strwidth(filogenia$tip.label, cex = cex_lbl)) else 0
        start_x    <- max_tip_x + config$label_offset + max_lbl_w + max_tip_x * 0.02
        x_columnas <- start_x + seq(0, n_car - 1L) * sep_col

        y_range  <- diff(range(yy_tip))
        y_header <- max(yy_tip) + y_range * 0.04
        old_xpd <- par("xpd"); par(xpd = NA)
        text(x_columnas, y_header,
             labels = caracteres,
             cex    = cex_aj * 0.9,
             srt    = 45,
             adj    = c(0, 0.5),
             font   = 2L)
        par(xpd = old_xpd)

        fig_x <- fig_y <- fig_col <- NULL
        for (i in seq_along(caracteres)) {
          col_tips <- sapply(datos_ord[[caracteres[i]]], resolve_color,
                             colores_estado = colores_por_car[[caracteres[i]]],
                             car = caracteres[i])
          fig_x  <- c(fig_x,   rep(x_columnas[i], n_tips))
          fig_y   <- c(fig_y,   yy_tip)
          fig_col <- c(fig_col, col_tips)
        }
        points(fig_x, fig_y, pch = pch_fig, bg = fig_col, col = "black", cex = tam_fig)

        .draw_simple_legend(pos_leyenda)
      }
    }

    invisible(NULL)
  }  # end .render_simple_mapping


  # ── PHASE 1: off-screen probe ───────────────────────────────────────────────
  .fan_coords_probe <- NULL
  computed_lim <- list(x = NULL, y = NULL)

  if (tipo_arbol == "fan") {
    prev_dev <- dev.cur()
    pdf(NULL, width = 7, height = 7)
    .fan_coords_probe <- tryCatch({
      plot(filogenia,
           type           = "fan",
           cex            = cex_aj,
           label.offset   = config$label_offset,
           edge.width     = config$branch_width,
           tip.color      = "transparent",
           show.tip.label = FALSE,
           no.margin      = TRUE)
      get("last_plot.phylo", envir = .PlotPhyloEnv)
    }, error = function(e) NULL,
    finally = .close_device_restore(prev_dev))
  }

  if (tipo_arbol != "fan" && simple_mode == "figures") {
    n_tips_tmp    <- n_tips
    default_height  <- n_tips_tmp * 0.25 + 2
    default_width <- 12
    height_in_tmp   <- if (!is.null(config$height)) config$height else default_height
    width_in_tmp  <- if (!is.null(config$width))  config$width  else default_width
    # On screen, probe with the size of the device that will be drawn on
    if (!exportar && dev.cur() > 1L) {
      din <- dev.size("in")
      width_in_tmp <- din[1L]; height_in_tmp <- din[2L]
    }

    prev_dev <- dev.cur()
    pdf(NULL, width = width_in_tmp, height = height_in_tmp)
    computed_lim <- tryCatch({

      par(mar = c(1, 1, 2, 1), xpd = FALSE)
      plot(filogenia,
           type           = tipo_arbol,
           cex            = cex_lbl,
           label.offset   = config$label_offset,
           edge.width     = config$branch_width,
           tip.color      = col_lbl,
           show.tip.label = show_lbl,
           font           = 3L)

      obj_s  <- get("last_plot.phylo", envir = .PlotPhyloEnv)
      xx_s   <- obj_s$xx[seq_len(n_tips)]
      yy_s   <- obj_s$yy[seq_len(n_tips)]
      usr_s  <- par("usr")
      pin_s  <- par("pin")

      sep_col_s  <- .figure_col_sep(cex_aj, tam_fig)
      max_tip_s  <- max(xx_s)
      max_lbl_s  <- if (show_lbl) max(strwidth(filogenia$tip.label, cex = cex_lbl)) else 0
      start_x_s  <- max_tip_s + config$label_offset + max_lbl_s + max_tip_s * 0.02
      x_cols_s   <- start_x_s + seq(0, n_car - 1L) * sep_col_s

      car_largo     <- caracteres[which.max(nchar(caracteres))]
      hdr_w_horiz   <- strwidth(car_largo, cex = cex_aj * 0.9) * cos(pi / 4)
      x_max_tabla   <- max(x_cols_s) + hdr_w_horiz

      pos_ley_s  <- if (!is.null(config$legend_corner)) config$legend_corner else "topright"
      on_right_s <- grepl("right", pos_ley_s)

      # Vertical legend: reserve its width beside the table. Horizontal legend:
      # reserve its height above / below the tree (it spans the bottom or top).
      horiz <- identical(config$legend_layout, "horizontal")
      lay   <- if (show_leg)
        .legend_layout(.simple_legend_blocks(), pos_ley_s,
                       config$legend_layout %||% "vertical",
                       cex = cex_aj, pt_cex = tam_fig, pch = pch_fig) else list(w = 0, h = 0)

      x_max_total <- x_max_tabla + if (horiz) 0 else lay$w * 1.1
      ylim <- if (horiz && lay$h > 0) {
        if (grepl("top", pos_ley_s)) c(usr_s[3L], usr_s[4L] + lay$h * 1.1)
        else                         c(usr_s[3L] - lay$h * 1.1, usr_s[4L])
      }

      list(x = if (x_max_total > usr_s[2L]) c(usr_s[1L], x_max_total),
           y = ylim)

    }, error = function(e) list(x = NULL, y = NULL),
    finally = .close_device_restore(prev_dev))
  }

  # ── PHASE 2: export or draw on screen ───────────────────────────────────────
  if (exportar) {
    .export_device(fn_export, config$export_format %||% "png",
                   filogenia, tipo_arbol,
                   {
                     .render_simple_mapping(xlim_extra = computed_lim$x,
                                            ylim_extra = computed_lim$y)
                   },
                   width  = config$width,
                   height = config$height)
  } else {
    .render_simple_mapping(xlim_extra = computed_lim$x, ylim_extra = computed_lim$y)
  }
}


# ==============================================================================
# BLOCK 3 — ANCESTRAL CHARACTER RECONSTRUCTION
# ==============================================================================

# ------------------------------------------------------------------------------
# 3.1  Algorithms
# ------------------------------------------------------------------------------

#' Initializes the edge color vector in gray
#'
#' @param filogenia  phylo object.
#' @return Named character vector of length nrow(filogenia$edge).
init_edge_colors <- function(filogenia) {
  ec <- rep("gray70", nrow(filogenia$edge))
  names(ec) <- paste0(filogenia$edge[, 1], "-", filogenia$edge[, 2])
  ec
}


#' Resolves the single color used for missing ("?"), inapplicable ("-") and
#' ambiguous states under the Fitch algorithm
#'
#' Centralizes the fallback chain so every call site (legend entries, branch
#' coloring in fitch.R, and tip-figure coloring) agrees on the same color:
#' \code{config$ambiguity_color} when the user supplied one, otherwise
#' \code{FITCH_AMBIG_COLOR} (magenta) from fitch.R.
#'
#' @param config Configuration list (uses config$ambiguity_color).
#' @return A single color string.
#' @keywords internal
.get_fitch_ambig_color <- function(config) {
  if (!is.null(config$ambiguity_color)) {
    config$ambiguity_color
  } else if (exists("FITCH_AMBIG_COLOR", mode = "character")) {
    get("FITCH_AMBIG_COLOR")
  } else {
    "#FF00FF"
  }
}


#' Injects the Fitch ambiguity entry into a state->color vector (Fitch only)
#'
#' When the selected algorithm is Fitch (config$algorithm == "fitch"), appends
#' an "Ambiguous / Missing" entry with the color from
#' \code{.get_fitch_ambig_color()} to \code{colores_estado} so the legend
#' always shows what the ambiguous branches/figures mean.
#' For any other algorithm the vector is returned unchanged.
#'
#' @param colores_estado Named character vector state -> color.
#' @param config         Configuration list (uses config$algorithm).
#' @return Possibly extended named character vector.
#' @keywords internal
.add_fitch_ambig_to_legend <- function(colores_estado, config) {
  if (!identical(config$algorithm, "fitch")) return(colores_estado)
  c(colores_estado, c("Ambiguous / Missing" = .get_fitch_ambig_color(config)))
}


#' Normalizes the colors chosen for missing ("?") and inapplicable ("-") tips
#'
#' Used by simple mapping and the depth-weighted majority reconstruction,
#' where each character may give "?" and "-" their own color (Fitch always uses the single
#' ambiguity color instead, see \code{.get_fitch_ambig_color()}).
#'
#' @param missing_colors \code{NULL}; a vector named \code{"?"} and/or
#'   \code{"-"} applied to every character; or a named list
#'   \emph{character} -> such a vector.
#' @param characters Names of the mapped characters.
#' @return Named list \emph{character} -> \code{c("?" = color, "-" = color)}
#'   (only the symbols supplied), or \code{NULL}.
#' @keywords internal
.normalize_missing_colors <- function(missing_colors, characters) {
  if (is.null(missing_colors) || length(missing_colors) == 0L) return(NULL)
  check <- function(v) {
    v <- unlist(v)
    if (is.null(names(v)) || !all(names(v) %in% c("?", "-")))
      stop("`missing_colors` must be named \"?\" and/or \"-\".")
    bad <- !vapply(v, is_valid_color, logical(1))
    if (any(bad))
      stop("`missing_colors`: invalid color(s): ", paste(v[bad], collapse = ", "), ".")
    v
  }
  if (is.list(missing_colors)) {
    if (is.null(names(missing_colors)))
      stop("A `missing_colors` list must be named by character.")
    out <- lapply(missing_colors, check)
    out[intersect(names(out), characters)]
  } else {
    v <- check(missing_colors)
    stats::setNames(rep(list(v), length(characters)), characters)
  }
}

#' Color of a tip state (simple mapping and ancestral reconstruction)
#'
#' Observed states take their mapped color. Missing (\code{"?"}, empty or
#' \code{NA}) and inapplicable (\code{"-"}) tips take the Fitch ambiguity color
#' under Fitch, or the per-character color from \code{config$missing_colors}
#' in simple mapping and under the depth-weighted majority algorithm (gray70
#' when none was chosen).
#' Unmapped (unticked) states stay gray70.
#'
#' @param valor          Tip state.
#' @param colores_estado Named vector state -> color of the character.
#' @param config         Configuration list.
#' @param car            Character name.
#' @return A single color string.
#' @keywords internal
.ancestral_tip_color <- function(valor, colores_estado, config, car) {
  valor <- as.character(valor)
  if (!is.na(valor) && valor %in% names(colores_estado)) return(colores_estado[[valor]])
  if (is.na(valor) || valor %in% c("", "?", "-")) {
    if (identical(config$algorithm, "fitch")) return(.get_fitch_ambig_color(config))
    sym <- if (is.na(valor) || valor == "") "?" else valor
    mc  <- config$missing_colors[[car]]
    if (!is.null(mc) && sym %in% names(mc)) return(unname(mc[[sym]]))
  }
  "gray70"
}

#' Legend entries for the missing / inapplicable colors of one character
#'
#' Simple mapping and depth-weighted majority (Fitch shows one
#' "Ambiguous / Missing" entry).
#' An entry is added when the user chose a color for that symbol and the
#' symbol occurs in the character's data.
#'
#' @param config Configuration list.
#' @param car    Character name.
#' @return Named vector label -> color, or \code{NULL}.
#' @keywords internal
.missing_legend_entries <- function(config, car) {
  if (identical(config$algorithm, "fitch")) return(NULL)
  mc <- config$missing_colors[[car]]
  if (is.null(mc)) return(NULL)
  vals    <- as.character(config$aligned_data[[car]])
  present <- c("?" = any(is.na(vals) | vals %in% c("", "?")), "-" = any(vals %in% "-"))
  syms    <- intersect(names(mc), names(present)[present])
  syms    <- setdiff(syms, names(config$colors_by_character[[car]]))
  if (length(syms) == 0L) return(NULL)
  labels  <- c("?" = "Missing (?)", "-" = "Inapplicable (-)")
  stats::setNames(unname(mc[syms]), labels[syms])
}

#' Appends the missing / inapplicable entries to grouped legend data
#' @param ley_data Output of \code{.build_legend_data()}.
#' @param config   Configuration list.
#' @return \code{ley_data} with the extra entries.
#' @keywords internal
.add_missing_to_legend_data <- function(ley_data, config) {
  for (nm in names(ley_data$by_char)) {
    extra <- .missing_legend_entries(config, nm)
    if (length(extra) == 0L) next
    ley_data$by_char[[nm]]$labels <- c(ley_data$by_char[[nm]]$labels, names(extra))
    ley_data$by_char[[nm]]$colors <- c(ley_data$by_char[[nm]]$colors, unname(extra))
  }
  ley_data
}


#' Dispatches to the ancestral reconstruction algorithm selected in config
#'
#' @param tree         phylo object.
#' @param tip_colors   Named character vector: tip label -> color.
#' @param edge_colors  Named character vector from init_edge_colors().
#' @param config       Configuration list.
#' @return Updated edge_colors vector.
apply_ancestral_algorithm <- function(tree, tip_colors, edge_colors, config) {

  if (identical(config$algorithm, "majority")) {
    if (!exists("default_algorithm", mode = "function"))
      stop("Function 'default_algorithm' is not available. ",
           "Make sure default_alg.R is part of the package.")
    default_algorithm(tree, tip_colors, edge_colors)

  } else if (identical(config$algorithm, "fitch")) {
    if (!exists("external_algorithm", mode = "function"))
      stop("Function 'external_algorithm' is not available. ",
           "Make sure fitch.R is part of the package.")
    external_algorithm(tree, tip_colors, edge_colors, config)

  } else {
    stop("Unrecognized algorithm.")
  }
}


# ------------------------------------------------------------------------------
# 3.3  Branch painting function
# ------------------------------------------------------------------------------

#' Paints the phylogeny (phylogram) with the colors of a single ancestral mapping.
#'
#' Draws one ancestral mapping on the phylogeny.
#' For fan and cladogram it preserves the original plot.phylo() behavior.
#'
#' @param filogenia      phylo object.
#' @param edge_colors    Color vector for each branch.
#' @param config         Configuration list.
#' @param titulo_leyenda String for the legend title.
#' @param colores_estado Named vector state -> color for the legend.
#' @return Invisible NULL. Draws on the active graphics device.
plot_ancestral_branches <- function(filogenia, edge_colors, config,
                                    titulo_leyenda = "", colores_estado = NULL) {
  exportar   <- !is.null(config$export_filename)
  fn_export  <- config$export_filename
  tipo_arbol <- config$tree_type

  cat("\n=== Generating plot (Ancestral Reconstruction) ===\n")

  if (exportar) {
    cat("    \u2192 Exporting to:", paste0(fn_export, ".", config$export_format %||% "png"), "\n")

    ley_data <- if (!is.null(colores_estado) && length(colores_estado) > 0)
      list(labels = names(colores_estado), colors = unname(colores_estado))
    else
      list(labels = NULL, colors = NULL)

    export_multimapr_tree(
      tree          = filogenia,
      color_list    = list(edge_colors),
      filename      = fn_export,
      type          = tipo_arbol,
      format        = config$export_format %||% "png",
      lwd           = config$branch_width,
      width         = config$width,
      height        = config$height,
       ladderize     = config$ladderize %||% FALSE,
       terminal_stretch = config$terminal_stretch %||% 1,
       tip_end          = config$tip_end %||% "round",
       branch_angle     = config$branch_angle %||% 45,
       show_labels      = config$show_labels %||% TRUE,
       label_size       = config$label_size %||% 1,
       label_color      = config$label_color %||% "black",
       legend_labels = if (config$show_legend %||% TRUE) ley_data$labels,
       legend_colors = ley_data$colors,
       legend_corner = config$legend_corner %||% "bottomleft",
       legend_layout = config$legend_layout %||% "vertical",
       node_states   = isTRUE(config$show_nodes),
       node_shape    = config$node_shape %||% "circle",
       node_size     = config$node_size %||% 1,
       legend_title  = titulo_leyenda
    )
  }

  # SCREEN RENDER — advanced visual engine for ANY topology
  plot_multimapr_screen(
    tree           = filogenia,
    color_list     = list(edge_colors),
    type           = tipo_arbol,
    lwd            = config$branch_width,
    ladderize      = config$ladderize %||% FALSE,
      terminal_stretch = config$terminal_stretch %||% 1,
      tip_end          = config$tip_end %||% "round",
      branch_angle     = config$branch_angle %||% 45,
      show_labels      = config$show_labels %||% TRUE,
      label_size       = config$label_size %||% 1,
      label_color      = config$label_color %||% "black",
    legend_labels  = if (!is.null(colores_estado) && config$show_legend %||% TRUE) names(colores_estado),
    legend_colors  = if (!is.null(colores_estado)) unname(colores_estado) else NULL,
    legend_corner  = config$legend_corner %||% "bottomleft",
    legend_layout  = config$legend_layout %||% "vertical",
    node_states    = isTRUE(config$show_nodes),
    node_shape     = config$node_shape %||% "circle",
    node_size      = config$node_size %||% 1,
    legend_title   = titulo_leyenda
  )

  invisible(NULL)
}


#' FULL ANCESTRAL RECONSTRUCTION — orchestrates 3.1 and 3.3
#'
#' Dispatches to the appropriate render function based on config$multi_function:
#'   1 -> branch superimposition only (plot_superimposed_characters /
#'         plot_ancestral_branches for single character)
#'   2 -> ancestral reconstruction on branches PLUS simple terminal figures
#'         (plot_ancestral_with_terminals)
#'
#' @param filogenia  phylo object.
#' @param config     Configuration list from setup_mapping_config().
#' @return Invisible NULL.
plot_ancestral_reconstruction <- function(filogenia, config) {
  datos_ord       <- config$aligned_data
  caracteres      <- config$characters
  colores_por_car <- config$colors_by_character

  # multi_function == 2: ancestral reconstruction + terminal figures
  if (!is.null(config$multi_function) && config$multi_function == 2) {
    plot_ancestral_with_terminals(filogenia, config)
    return(invisible(NULL))
  }

  # multi_function == 1 (or undefined): branch superimposition only
  if (length(caracteres) == 1) {
    car            <- caracteres[1]
    colores_estado <- colores_por_car[[car]]
    tip_colors     <- sapply(datos_ord[[car]], .ancestral_tip_color,
                             colores_estado = colores_estado, config = config, car = car)

    edge_colors <- init_edge_colors(filogenia)
    edge_colors <- apply_ancestral_algorithm(filogenia, tip_colors, edge_colors, config)

    # For Fitch, append the ambiguity entry to the legend (Fitch-only)
    colores_estado_ley <- c(.add_fitch_ambig_to_legend(colores_estado, config),
                            .missing_legend_entries(config, car))

    plot_ancestral_branches(filogenia, edge_colors, config,
                            titulo_leyenda = car,
                            colores_estado = colores_estado_ley)

  } else {
    plot_superimposed_characters(filogenia, config)
  }
}


# ------------------------------------------------------------------------------
# 3.4  Ancestral reconstruction + terminal figures
# ------------------------------------------------------------------------------

#' Ancestral reconstruction on branches PLUS simple terminal figures
#'
#' Combines two visualization layers in a single plot:
#'   Layer 1 — Branches painted with the ancestral reconstruction result
#'             (using the algorithm selected in config$algorithm).
#'   Layer 2 — Colored figures (circles or squares) at the tips showing the
#'             observed tip states (equivalent to the simple mapping squares).
#'
#' The function uses plot_multimapr_screen (and export_multimapr_tree for
#' file output) to draw the ancestral layer, then overlays tip figures
#' with par(new = TRUE) / points() on the same coordinate system.
#'
#' @param filogenia  phylo object.
#' @param config     Configuration list from setup_mapping_config().
#'                   Must contain: aligned_data, characters, colors_by_character,
#'                   algorithm, tree_type, branch_width, figure_pch, figure_size,
#'                   and optionally export_filename / export_format.
#' @return Invisible NULL. Draws on the active graphics device.
plot_ancestral_with_terminals <- function(filogenia, config) {

  datos_ord       <- config$aligned_data
  caracteres      <- config$characters
  colores_por_car <- config$colors_by_character
  tipo_arbol      <- config$tree_type
  pch_fig         <- config$figure_pch %||% 22L   # 22 = square
  tam_fig         <- config$figure_size %||% 1
  exportar        <- !is.null(config$export_filename)
  fn_export       <- config$export_filename
  n_tips          <- Ntip(filogenia)

  # Missing ("?") and inapplicable ("-") tips use the same color on figures and
  # branches: the Fitch ambiguity color, or the per-character colors chosen for
  # the depth-weighted majority (see .ancestral_tip_color()).
  resolve_color <- function(valor, colores_estado, car)
    .ancestral_tip_color(valor, colores_estado, config, car)

  cat("\n=== Generating plot (Ancestral Reconstruction + Terminal Figures) ===\n")

  # ── Step 1: compute ancestral edge colors for each character ────────────────
  lista_ec <- lapply(caracteres, function(car) {
    col_e      <- colores_por_car[[car]]
    tip_colors <- sapply(datos_ord[[car]], resolve_color, colores_estado = col_e, car = car)
    ec         <- init_edge_colors(filogenia)
    ec         <- apply_ancestral_algorithm(filogenia, tip_colors, ec, config)
    ec
  })

  # ── Step 2: build legend data (branches block, grouped by character) ─────────
  ley_data <- .build_legend_data(colores_por_car, caracteres)

  # For Fitch: append the ambiguity entry to every character block (Fitch-only)
  if (identical(config$algorithm, "fitch")) {
    ambig_col <- .get_fitch_ambig_color(config)
    ambig_label <- "Ambiguous / Missing"
    for (nm in names(ley_data$by_char)) {
      if (!ambig_label %in% ley_data$by_char[[nm]]$labels) {
        ley_data$by_char[[nm]]$labels <- c(ley_data$by_char[[nm]]$labels, ambig_label)
        ley_data$by_char[[nm]]$colors <- c(ley_data$by_char[[nm]]$colors, ambig_col)
      }
    }
  }
  # Majority: per-character missing / inapplicable colors (when chosen)
  ley_data <- .add_missing_to_legend_data(ley_data, config)

  # ── Helper: overlay tip figures on an already-rendered phylogeny ────────────
  # Reads tip coordinates from .PlotPhyloEnv after the base plot is drawn.
  # ── Helper: overlay tip figures ──────────────────────────────────────────────
  # pp, cex_aj, lbl_off: when called from overlay_fn (export) they are
  # received from the export_multimapr_tree engine, which already computed them
  # with the correct scale factors. When called from screen (NULL) they are
  # read from .PlotPhyloEnv and adjust_cex() respectively.
  .overlay_tip_figures <- function(pp       = NULL,
                                             cex_aj   = NULL,
                                             lbl_off  = NULL,
                                             R_tips   = NULL,   # radius of each terminal (fan)
                                             gap_u    = NULL) { # gap laboral base (fan)
    # gap_u_fan: local name to avoid conflict with R_tips/gap_u parameters
    gap_u_fan <- gap_u

    if (is.null(pp))
      pp <- get("last_plot.phylo",
                envir = get(".PlotPhyloEnv", envir = asNamespace("ape")))

    xx_tip <- pp$xx[seq_len(n_tips)]
    yy_tip <- pp$yy[seq_len(n_tips)]

    # cex_aj: use the value from the export engine if passed, otherwise recalculate
    if (is.null(cex_aj))
      cex_aj <- max(1 / (1 + n_tips / 50), 0.2)

    old_xpd <- par("xpd"); par(xpd = TRUE)

    if (tipo_arbol == "fan") {
      # Fan: exact replica of the simple mapping geometry (plot_simple_mapping).
      #
      # Geometry (identical to simple mode):
      #   max_tip_radius   = max(R_tips)
      #   radio_base       = max_tip_radius * 1.06   <- first ring
      #   incremento_radio = max_tip_radius * 0.10   <- spacing between rings
      #   radio_actual[i]  = radio_base + (i-1) * incremento_radio
      #   gap_tan          = incremento_radio * 0.55  <- tangential offset for "C_i"
      #   radio_nombres    = radio_ultimo + incremento_radio * 0.5
      #
      # R_tips may come from the export engine (overlay_fn) or be computed here.
      angulos <- atan2(yy_tip, xx_tip)
      if (is.null(R_tips)) R_tips <- sqrt(xx_tip^2 + yy_tip^2)
      max_tip_radius   <- max(R_tips)

      radio_base       <- max_tip_radius * 1.06
      incremento_radio <- max_tip_radius * 0.10
      n_car_ov         <- length(caracteres)
      radio_ultimo     <- radio_base + (n_car_ov - 1L) * incremento_radio
      gap_etiq         <- incremento_radio * 0.5
      radio_nombres    <- radio_ultimo + gap_etiq
      gap_tan          <- incremento_radio * 0.55

      for (i in seq_len(n_car_ov)) {
        col_tips     <- sapply(datos_ord[[caracteres[i]]], resolve_color,
                               colores_estado = colores_por_car[[caracteres[i]]],
                               car = caracteres[i])
        radio_actual <- radio_base + (i - 1L) * incremento_radio

        # ── Ring i figures ───────────────────────────────────────────────────
        points(radio_actual * cos(angulos),
               radio_actual * sin(angulos),
               pch = pch_fig, bg = col_tips, col = "black", cex = tam_fig)

        # ── Tangential "C_i" label (same as simple mode) ────────────────────
        for (j in seq_len(n_tips)) {
          ang_j  <- angulos[j]
          deg_j  <- ang_j * 180 / pi
          lado_d <- cos(ang_j) >= 0
          srt_j  <- if (lado_d) deg_j else deg_j + 180
          tan_x  <- -sin(ang_j)
          tan_y  <-  cos(ang_j)
          cx     <- radio_actual * cos(ang_j) + tan_x * gap_tan
          cy     <- radio_actual * sin(ang_j) + tan_y * gap_tan
          text(cx, cy,
               labels = paste0("C", i),
               adj    = c(0.5, 0.5),
               cex    = cex_aj * 0.65,
               font   = 2L,
               srt    = srt_j)
        }
      }

      # ── Species labels (parallel to radius, beyond the last ring)
      for (j in seq_len(n_tips)[config$show_labels %||% TRUE]) {
        ang_j  <- angulos[j]
        deg_j  <- ang_j * 180 / pi
        lado_d <- cos(ang_j) >= 0
        srt_j  <- if (lado_d) deg_j else deg_j + 180
        adj_j  <- if (lado_d) c(0, 0.5) else c(1, 0.5)
        text(radio_nombres * cos(ang_j),
             radio_nombres * sin(ang_j),
             labels = filogenia$tip.label[j],
             adj    = adj_j,
             cex    = cex_aj * (config$label_size %||% 1),
             col    = config$label_color %||% "black",
             srt    = srt_j,
             font   = 3L)
      }

    } else {
      # Phylogram / cladogram: figures placed to the right of each tip label.
      # lbl_off: from the export engine if available, or proportional to x_max
      max_tip_x <- max(xx_tip)
      if (is.null(lbl_off))
        lbl_off <- config$label_offset %||% (max_tip_x * 0.025)

      max_lbl_w  <- if (config$show_labels %||% TRUE)
                      max(strwidth(filogenia$tip.label,
                                   cex = cex_aj * (config$label_size %||% 1))) else 0
      sep_col    <- .figure_col_sep(cex_aj, tam_fig)
      start_x    <- max_tip_x + lbl_off + max_lbl_w + max_tip_x * 0.02
      x_columnas <- start_x + seq(0, length(caracteres) - 1L) * sep_col

      # ── Column headers rotated 45 degrees (same as simple mode) ─────────────
      y_range  <- diff(range(yy_tip))
      y_header <- max(yy_tip) + y_range * 0.04
      old_xpd2 <- par("xpd"); par(xpd = NA)
      text(x_columnas, y_header,
           labels = caracteres,
           cex    = cex_aj * 0.9,
           srt    = 45,
           adj    = c(0, 0.5),
           font   = 2L)
      par(xpd = old_xpd2)

      # ── Figures per column ────────────────────────────────────────────────────
      for (i in seq_along(caracteres)) {
        col_tips <- sapply(datos_ord[[caracteres[i]]], resolve_color,
                           colores_estado = colores_por_car[[caracteres[i]]],
                           car = caracteres[i])
        points(rep(x_columnas[i], n_tips), yy_tip,
               pch = pch_fig, bg = col_tips, col = "black", cex = tam_fig)
      }
    }

    par(xpd = old_xpd)
    invisible(NULL)
  }

  # ── Step 3a: export to file ─────────────────────────────────────────────────
  # Fully delegated to export_multimapr_tree(), which manages the device,
  # dimensions, legend, and coordinates. Terminal figures are injected via
  # overlay_fn, which receives pp, cex_aj and label_offset_aj already computed
  # with the correct export scale factors, ensuring figures are positioned
  # correctly next to labels rather than on top of them.
  if (exportar) {
    cat("    \u2192 Exporting to:", paste0(fn_export, ".", config$export_format %||% "png"), "\n")

    export_multimapr_tree(
      tree            = filogenia,
      color_list      = lista_ec,
      filename        = fn_export,
      type            = tipo_arbol,
      format          = config$export_format %||% "png",
      lwd             = config$branch_width,
      width           = config$width,
      height          = config$height,
      offset_range    = config$offset_range %||% 0.1,
      tip_end         = config$tip_end %||% "round",
      branch_angle    = config$branch_angle %||% 45,
      show_labels     = config$show_labels %||% TRUE,
      label_size      = config$label_size %||% 1,
      label_color     = config$label_color %||% "black",
      legend_by_char  = if (config$show_legend %||% TRUE) ley_data$by_char,
      legend_corner   = config$legend_corner %||% "bottomleft",
      legend_layout   = config$legend_layout %||% "vertical",
      node_states     = isTRUE(config$show_nodes),
      node_shape      = config$node_shape %||% "circle",
      node_size       = config$node_size %||% 1,
      # For fan, hide engine labels and redraw them in the overlay,
      # further out, after all figures.
      hide_fan_labels = (tipo_arbol == "fan"),
      header_labels   = caracteres,
      overlay_fn      = function(pp, cex_aj, label_offset_aj,
                                 R_tips = NULL, gap_u = NULL) {
        .overlay_tip_figures(pp      = pp,
                                       cex_aj  = cex_aj,
                                       lbl_off = label_offset_aj,
                                       R_tips  = R_tips,
                                       gap_u   = gap_u)
      }
    )
  }

  # ── Step 3b: screen render ──────────────────────────────────────────────────
  # Always render to screen regardless of export flag.
  # For fan, hide_fan_labels = TRUE delegates labels to the overlay, which
  # redraws them further out, beyond the figure rings.
  plot_multimapr_screen(
    tree            = filogenia,
    color_list      = lista_ec,
    type            = tipo_arbol,
    lwd             = config$branch_width,
    offset_range    = config$offset_range %||% 0.1,
    tip_end         = config$tip_end %||% "round",
    branch_angle    = config$branch_angle %||% 45,
    show_labels     = config$show_labels %||% TRUE,
    label_size      = config$label_size %||% 1,
    label_color     = config$label_color %||% "black",
    legend_by_char  = if (config$show_legend %||% TRUE) ley_data$by_char,
    legend_corner   = config$legend_corner %||% "bottomleft",
    legend_layout   = config$legend_layout %||% "vertical",
    node_states     = isTRUE(config$show_nodes),
    node_shape      = config$node_shape %||% "circle",
    node_size       = config$node_size %||% 1,
    hide_fan_labels = (tipo_arbol == "fan"),
    header_labels   = caracteres,
    overlay_fn      = function(pp, cex_aj, label_offset_aj,
                               R_tips = NULL, gap_u = NULL) {
      .overlay_tip_figures(pp      = pp,
                                     cex_aj  = cex_aj,
                                     lbl_off = label_offset_aj,
                                     R_tips  = R_tips,
                                     gap_u   = gap_u)
    }
  )

  invisible(NULL)
}


# ==============================================================================
# MULTI-CHARACTER
# ==============================================================================

#' Superimposition of multiple characters on branches (advanced engine)
#'
#' Uses the Blank Canvas visual engine (plot_multimapr_screen) to render
#' on screen all topologies with the same quality as file export.
#' Isotropic offsets and the character-grouped legend are calculated and
#' drawn identically to export_multimapr_tree.
#'
#' @param filogenia        \code{phylo} object.
#' @param config           Configuration list from \code{setup_mapping_config()}.
#' @param alpha1           Alpha transparency for the first character layer (default \code{1}).
#' @param alpha2           Alpha transparency for additional character layers (default \code{1}).
#' @param offset_h_map     Horizontal offset between layers in data units (default \code{0.02}).
#' @param offset_v_map     Vertical offset between layers in data units (default \code{0.02}).
#' @param extra_offset_map Extra separation added to fan layers (default \code{0.07}).
#' @param escala_fan       Scale factor for fan topology layers (default \code{0.75}).
#' @param rotacion_fan     Rotation in radians between fan layers (default \code{2}).
#' @param base_escala      Base scale factor applied before \code{escala_fan} (default \code{1}).
#' @param legend_spacing   Vertical spacing between legend blocks in data units (default \code{0.06}).
#' @param x_limit          Optional numeric; overrides the x-axis upper limit. \code{NULL} = automatic.
#' @return Invisible NULL. Draws on the active graphics device.
plot_superimposed_characters <- function(filogenia, config,
                                         alpha1           = 1,
                                         alpha2           = 1,
                                         offset_h_map     = 0.02,
                                         offset_v_map     = 0.02,
                                         extra_offset_map = 0.07,
                                         escala_fan       = 0.75,
                                         rotacion_fan     = 2,
                                         base_escala      = 1,
                                         legend_spacing   = 0.06,
                                         x_limit          = NULL) {

  datos_ord       <- config$aligned_data
  caracteres      <- config$characters
  colores_por_car <- config$colors_by_character
  tipo_arbol      <- config$tree_type
  grosor1         <- config$branch_width %||% 2

  resolve_color <- function(valor, colores_estado, car)
    .ancestral_tip_color(valor, colores_estado, config, car)

  lista_ec <- lapply(caracteres, function(car) {
    col_e      <- colores_por_car[[car]]
    tip_colors <- sapply(datos_ord[[car]], resolve_color, colores_estado = col_e, car = car)
    ec         <- init_edge_colors(filogenia)
    ec         <- apply_ancestral_algorithm(filogenia, tip_colors, ec, config)
    ec
  })

  exportar  <- !is.null(config$export_filename)
  fn_export <- config$export_filename
  ley_data  <- .build_legend_data(colores_por_car, caracteres)

  # For Fitch: append the ambiguity entry to every character block (Fitch-only)
  if (identical(config$algorithm, "fitch")) {
    ambig_col <- .get_fitch_ambig_color(config)
    ambig_label <- "Ambiguous / Missing"
    for (nm in names(ley_data$by_char)) {
      if (!ambig_label %in% ley_data$by_char[[nm]]$labels) {
        ley_data$by_char[[nm]]$labels <- c(ley_data$by_char[[nm]]$labels, ambig_label)
        ley_data$by_char[[nm]]$colors <- c(ley_data$by_char[[nm]]$colors, ambig_col)
      }
    }
  }
  # Majority: per-character missing / inapplicable colors (when chosen)
  ley_data <- .add_missing_to_legend_data(ley_data, config)

  if (exportar) {
    cat("    \u2192 Exporting to:", paste0(fn_export, ".", config$export_format %||% "png"), "\n")
    export_multimapr_tree(
      tree            = filogenia,
      color_list      = lista_ec,
      filename        = fn_export,
      type            = tipo_arbol,
      format          = config$export_format %||% "png",
      lwd             = grosor1,
      width           = config$width,
      height          = config$height,
      offset_range    = config$offset_range %||% 0.1,
      use_edge_length = isTRUE(config$use_edge_length),
       ladderize       = config$ladderize %||% FALSE,
       terminal_stretch = config$terminal_stretch %||% 1,
       tip_end          = config$tip_end %||% "round",
       branch_angle     = config$branch_angle %||% 45,
       show_labels      = config$show_labels %||% TRUE,
       label_size       = config$label_size %||% 1,
       label_color      = config$label_color %||% "black",
       legend_by_char  = if (config$show_legend %||% TRUE) ley_data$by_char,
       legend_corner   = config$legend_corner %||% "bottomleft",
       legend_layout   = config$legend_layout %||% "vertical",
       node_states     = isTRUE(config$show_nodes),
       node_shape      = config$node_shape %||% "circle",
       node_size       = config$node_size %||% 1
     )
   }

  # SCREEN RENDER — integrated advanced visual engine
  plot_multimapr_screen(
    tree            = filogenia,
    color_list      = lista_ec,
    type            = tipo_arbol,
    lwd             = grosor1,
    offset_range    = config$offset_range %||% 0.1,
    use_edge_length = isTRUE(config$use_edge_length),
    ladderize       = config$ladderize %||% FALSE,
      terminal_stretch = config$terminal_stretch %||% 1,
      tip_end          = config$tip_end %||% "round",
      branch_angle     = config$branch_angle %||% 45,
      show_labels      = config$show_labels %||% TRUE,
      label_size       = config$label_size %||% 1,
      label_color      = config$label_color %||% "black",
    legend_by_char  = if (config$show_legend %||% TRUE) ley_data$by_char,
    legend_corner   = config$legend_corner %||% "bottomleft",
    legend_layout   = config$legend_layout %||% "vertical",
    node_states     = isTRUE(config$show_nodes),
    node_shape      = config$node_shape %||% "circle",
    node_size       = config$node_size %||% 1
  )

  invisible(NULL)
}
