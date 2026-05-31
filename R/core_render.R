################################################################################
# MultiMapR — Core Rendering Controllers
# Rendering and evolutionary state assignment module
################################################################################


# ==============================================================================
# BLOCK 2 — SIMPLE MAPPING
# ==============================================================================

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

  datos_ord       <- config$datos_ord
  caracteres      <- config$caracteres
  colores_por_car <- config$colores_por_caracter
  tipo_arbol      <- config$tipo_arbol
  pch_fig         <- config$pch_figura
  tam_fig         <- config$tam_figura
  cex_aj          <- adjust_cex(filogenia, config = config)
  n_tips          <- Ntip(filogenia)
  n_car           <- length(caracteres)

  asignar_color <- function(valor, colores_estado) {
    valor <- as.character(valor)
    if (is.na(valor) || valor == "") return("gray70")
    if (valor %in% names(colores_estado)) return(colores_estado[[valor]])
    "gray70"
  }

  exportar  <- !is.null(config$export_filename)
  fn_export <- config$export_filename

  cat("\n=== Generating plot (Simple Mapping) ===\n")
  if (exportar) cat("    \u2192 Exporting to:",
                    paste0(fn_export, ".", config$export_format %||% "png"), "\n")

  # ── Helper: draws the legend in blocks (same as plot_superimposed_characters) ──
  .draw_simple_legend <- function(pos_ley) {
    usr        <- par("usr")
    going_down <- grepl("top",  pos_ley)
    x_start    <- if (grepl("left", pos_ley)) usr[1L] else usr[2L]
    y_start    <- if (going_down)              usr[4L] else usr[3L]
    line_h     <- strheight("M", cex = cex_aj) * 1.4
    y_cursor   <- y_start
    old_xpd <- par("xpd"); par(xpd = TRUE)
    for (i in seq_along(caracteres)) {
      car   <- caracteres[i]
      col_e <- colores_por_car[[car]]
      titulo_ley <- if (tipo_arbol == "fan") paste0("C", i, " (", car, ")") else car
      lg <- legend(
        x      = x_start,
        y      = y_cursor,
        legend = names(col_e),
        pt.bg  = unname(col_e),
        col    = "black",
        pch    = pch_fig,
        title  = titulo_ley,
        bty    = "n",
        cex    = cex_aj,
        horiz  = FALSE,
        pt.cex = tam_fig,
        xjust  = if (grepl("right", pos_ley)) 1 else 0,
        yjust  = if (going_down) 1 else 0
      )
      bloque_h <- lg$rect$h
      y_cursor <- if (going_down) y_cursor - bloque_h - line_h
      else            y_cursor + bloque_h + line_h
    }
    par(xpd = old_xpd)
    invisible(NULL)
  }

  # ── Helper: measures the total legend width (plot=FALSE) ────────────────────
  .measure_legend_width <- function(pos_ley) {
    usr      <- par("usr")
    on_right <- grepl("right", pos_ley)
    x_ref    <- if (on_right) usr[2L] else usr[1L]
    y_ref    <- if (grepl("top", pos_ley)) usr[4L] else usr[3L]
    max_width <- 0
    for (i in seq_along(caracteres)) {
      car    <- caracteres[i]
      col_e  <- colores_por_car[[car]]
      titulo <- if (tipo_arbol == "fan") paste0("C", i, " (", car, ")") else car
      lg <- legend(x = x_ref, y = y_ref,
                   legend = names(col_e),
                   pch    = pch_fig,
                   col    = rep("black", length(col_e)),
                   title  = titulo,
                   bty    = "n", cex = cex_aj, horiz = FALSE,
                   pt.cex = tam_fig,
                   xjust  = if (on_right) 1 else 0,
                   yjust  = 1,
                   plot   = FALSE)
      max_width <- max(max_width, lg$rect$w)
    }
    max_width
  }

  # ── Main render function ─────────────────────────────────────────────────────
  # xlim_extra: if not NULL, passed as x.lim to plot() to expand the viewport
  # to the right to accommodate table + legend.
  .render_simple_mapping <- function(xlim_extra = NULL) {

    op <- par(no.readonly = TRUE)
    on.exit(par(op))

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
      obj <- .fan_coords_explorador

      xx_tip <- obj$xx[seq_len(n_tips)]
      yy_tip <- obj$yy[seq_len(n_tips)]

      angulos        <- atan2(yy_tip, xx_tip)
      tip_radius     <- sqrt(xx_tip^2 + yy_tip^2)
      max_tip_radius <- max(tip_radius)

      # ── Step 2: ring geometry ─────────────────────────────────────────────────
      radio_base       <- max_tip_radius * 1.06
      incremento_radio <- max_tip_radius * 0.10
      radio_ultimo     <- radio_base + (n_car - 1L) * incremento_radio

      # Space between last ring and label
      gap_etiq         <- incremento_radio * 0.5
      radio_nombres    <- radio_ultimo + gap_etiq

      # ── Step 3: estimate total radius including tip labels ────────────────────
      max_nchar_lbl <- max(nchar(filogenia$tip.label))
      char_w_u      <- max_tip_radius * 0.018 * cex_aj
      lbl_w_u       <- max_nchar_lbl * char_w_u
      radio_total   <- radio_nombres + lbl_w_u + max_tip_radius * 0.05

      # ── Step 4: replot with symmetric x.lim that contains everything ─────────
      xlim_fan <- c(-radio_total, radio_total)

      plot(filogenia,
           type           = "fan",
           cex            = cex_aj,
           label.offset   = config$label_offset,
           edge.width     = config$grosor,
           tip.color      = "transparent",
           show.tip.label = FALSE,
           no.margin      = TRUE,
           x.lim          = xlim_fan)

      # Refresh coordinates from the definitive plot
      obj    <- get("last_plot.phylo", envir = .PlotPhyloEnv)
      xx_tip <- obj$xx[seq_len(n_tips)]
      yy_tip <- obj$yy[seq_len(n_tips)]
      angulos <- atan2(yy_tip, xx_tip)

      # ── Symbol rings ───────────────────────────────────────────────────────────
      old_xpd <- par("xpd"); par(xpd = TRUE)
      for (i in seq_len(n_car)) {
        col_tips     <- sapply(datos_ord[[caracteres[i]]], asignar_color,
                               colores_estado = colores_por_car[[caracteres[i]]])
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

      # ── Tip labels — parallel to radius ──────────────────────────────────────
      old_xpd <- par("xpd"); par(xpd = TRUE)
      for (j in seq_len(n_tips)) {
        ang_j  <- angulos[j]; deg_j <- ang_j * 180 / pi
        lado_d <- cos(ang_j) >= 0
        srt_j  <- if (lado_d) deg_j else deg_j + 180
        adj_j  <- if (lado_d) c(0, 0.5) else c(1, 0.5)
        text(radio_nombres * cos(ang_j), radio_nombres * sin(ang_j),
             labels = filogenia$tip.label[j], adj = adj_j, cex = cex_aj, srt = srt_j,
             font   = 3L)
      }
      par(xpd = old_xpd)

      .draw_simple_legend(pos_leyenda)

      # ==========================================================================
      # PHYLOGRAM / CLADOGRAM
      # ==========================================================================
    } else {

      par(mar = c(1, 1, 2, 1), xpd = TRUE)

      plot(filogenia,
           type           = tipo_arbol,
           cex            = cex_aj,
           label.offset   = config$label_offset,
           edge.width     = config$grosor,
           tip.color      = "black",
           show.tip.label = TRUE,
           font           = 3L,
           x.lim          = xlim_extra)   # NULL on screen; expanded on export

      obj    <- get("last_plot.phylo", envir = .PlotPhyloEnv)
      xx_tip <- obj$xx[seq_len(n_tips)]
      yy_tip <- obj$yy[seq_len(n_tips)]
      usr    <- par("usr")
      pin    <- par("pin")

      sep_col    <- strwidth("M", cex = cex_aj) * 2.5
      max_tip_x  <- max(xx_tip)
      max_lbl_w  <- max(strwidth(filogenia$tip.label, cex = cex_aj))
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
        col_tips <- sapply(datos_ord[[caracteres[i]]], asignar_color,
                           colores_estado = colores_por_car[[caracteres[i]]])
        fig_x   <- c(fig_x,   rep(x_columnas[i], n_tips))
        fig_y   <- c(fig_y,   yy_tip)
        fig_col <- c(fig_col, col_tips)
      }
      points(fig_x, fig_y, pch = pch_fig, bg = fig_col, col = "black", cex = tam_fig)

      .draw_simple_legend(pos_leyenda)
    }

    invisible(NULL)
  }  # end .render_simple_mapping


  # ── PHASE 1: off-screen probe ───────────────────────────────────────────────
  .fan_coords_explorador <- NULL
  computed_xlim <- NULL

  if (tipo_arbol == "fan") {
    pdf(NULL, width = 7, height = 7)
    .fan_coords_explorador <- tryCatch({
      plot(filogenia,
           type           = "fan",
           cex            = cex_aj,
           label.offset   = config$label_offset,
           edge.width     = config$grosor,
           tip.color      = "transparent",
           show.tip.label = FALSE,
           no.margin      = TRUE)
      get("last_plot.phylo", envir = .PlotPhyloEnv)
    }, error = function(e) NULL,
    finally = dev.off())
  }

  if (tipo_arbol != "fan") {
    n_tips_tmp    <- n_tips
    default_height  <- n_tips_tmp * 0.25 + 2
    default_width <- 12
    height_in_tmp   <- if (!is.null(config$height)) config$height else default_height
    width_in_tmp  <- if (!is.null(config$width))  config$width  else default_width

    pdf(NULL, width = width_in_tmp, height = height_in_tmp)
    computed_xlim <- tryCatch({

      par(mar = c(1, 1, 2, 1), xpd = FALSE)
      plot(filogenia,
           type           = tipo_arbol,
           cex            = cex_aj,
           label.offset   = config$label_offset,
           edge.width     = config$grosor,
           tip.color      = "black",
           show.tip.label = TRUE,
           font           = 3L)

      obj_s  <- get("last_plot.phylo", envir = .PlotPhyloEnv)
      xx_s   <- obj_s$xx[seq_len(n_tips)]
      yy_s   <- obj_s$yy[seq_len(n_tips)]
      usr_s  <- par("usr")
      pin_s  <- par("pin")

      sep_col_s  <- strwidth("M", cex = cex_aj) * 2.5
      max_tip_s  <- max(xx_s)
      max_lbl_s  <- max(strwidth(filogenia$tip.label, cex = cex_aj))
      start_x_s  <- max_tip_s + config$label_offset + max_lbl_s + max_tip_s * 0.02
      x_cols_s   <- start_x_s + seq(0, n_car - 1L) * sep_col_s

      car_largo     <- caracteres[which.max(nchar(caracteres))]
      hdr_w_horiz   <- strwidth(car_largo, cex = cex_aj * 0.9) * cos(pi / 4)
      x_max_tabla   <- max(x_cols_s) + hdr_w_horiz

      pos_ley_s  <- if (!is.null(config$legend_corner)) config$legend_corner else "topright"
      on_right_s <- grepl("right", pos_ley_s)

      legend_width <- 0
      y_ref_s   <- if (grepl("top", pos_ley_s)) usr_s[4L] else usr_s[3L]
      for (i in seq_along(caracteres)) {
        car_i   <- caracteres[i]
        col_e_i <- colores_por_car[[car_i]]
        lg_i <- legend(x      = x_max_tabla,
                       y      = y_ref_s,
                       legend = names(col_e_i),
                       pch    = pch_fig,
                       col    = rep("black", length(col_e_i)),
                       title  = car_i,
                       bty    = "n", cex = cex_aj, horiz = FALSE,
                       pt.cex = tam_fig, xjust = 0, yjust = 1,
                       plot   = FALSE)
        legend_width <- max(legend_width, lg_i$rect$w)
      }

      x_max_total <- x_max_tabla + legend_width * 1.1

      if (x_max_total > usr_s[2L]) {
        c(usr_s[1L], x_max_total)
      } else {
        NULL
      }

    }, error = function(e) NULL,
    finally = dev.off())
  }

  # ── PHASE 2: export or draw on screen ───────────────────────────────────────
  if (exportar) {
    .export_device(fn_export, config$export_format %||% "png",
                          filogenia, tipo_arbol,
                          {
                            .render_simple_mapping(xlim_extra = computed_xlim)
                          },
                          width  = config$width,
                          height = config$height)
  } else {
    .render_simple_mapping()
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


#' Dispatches to the ancestral reconstruction algorithm selected in config
#'
#' @param tree         phylo object.
#' @param tip_colors   Named character vector: tip label -> color.
#' @param edge_colors  Named character vector from init_edge_colors().
#' @param config       Configuration list.
#' @return Updated edge_colors vector.
apply_ancestral_algorithm <- function(tree, tip_colors, edge_colors, config) {

  if (config$algoritmo == 1) {
    if (!exists("default_algorithm", mode = "function"))
      stop("Function 'default_algorithm' is not available. ",
           "Make sure default_alg.R is part of the package.")
    default_algorithm(tree, tip_colors, edge_colors)

  } else if (config$algoritmo == 2) {
    if (!exists("external_algorithm", mode = "function"))
      stop("Function 'external_algorithm' is not available. ",
           "Make sure fitch.R is part of the package.")
    external_algorithm(tree, tip_colors, edge_colors, config)

  } else {
    stop("Unrecognized algorithm.")
  }
}


# ------------------------------------------------------------------------------
# 3.2  Drawing engine — Blank Canvas
# (dibujar_filograma_canvas retained as legacy name)
# ------------------------------------------------------------------------------

#' Draws a phylogram using the "Blank Canvas" strategy
#' with geometrically isotropic offset (X + Y) to visually separate
#' multi-mappings.
#'
#' Accepts a list of color vectors (one per evolutionary history).
#' If the list has 1 element -> simple mapping (no offset).
#' If it has N > 1 -> multi-mapping with zero-centered offset:
#'   · The user only provides `desfases_y` (vertical lanes).
#'   · `desfases_x` is calculated automatically just after
#'     `plot(..., type = "n")`, when the graphics engine already knows
#'     the canvas dimensions. The calculation guarantees that the physical
#'     separation (in inches) between verticals is identical to the
#'     physical separation between horizontals, regardless of window size
#'     or tree length:
#'
#'       usr      <- par("usr")   # data coord limits: c(x1,x2,y1,y2)
#'       pin      <- par("pin")   # physical dimensions in inches: c(w, h)
#'       escala_x <- (usr\[2\] - usr\[1\]) / pin\[1\]   # data units / inch on X
#'       escala_y <- (usr\[4\] - usr\[3\]) / pin\[2\]   # data units / inch on Y
#'       desfases_x <- desfases_y * (escala_x / escala_y)
#'
#' Horizontal strokes are vectorized with a single `segments()`.
#' Vertical strokes (internal node segments) are drawn by iterating over
#' each internal node, using the *incoming* edge color to guarantee a solid,
#' uniform color. `lend = 1` ensures perfectly square joins between
#' horizontals and verticals.
#' The vertical stroke is protected with `if (y_min != y_max)` to avoid
#' collisions in polytomies and zero-length branches.
#'
#' @param filogenia      phylo object.
#' @param lista_ec       List of color vectors. Each vector has length
#'                       `nrow(filogenia$edge)`.
#' @param grosor         Line width (lwd).
#' @param label_offset   Distance between tips and labels.
#' @param x_limit        Right limit of X axis (NULL = automatic).
#' @param desfases_y     Numeric vector of Y offsets (one per history).
#'                       If NULL, `seq(-0.08, 0.08, length.out = N)` is used.
#' @return Invisible NULL. Draws on the active graphics device.
dibujar_filograma_canvas <- function(filogenia, lista_ec,
                                     grosor       = 2,
                                     label_offset = 0.3,
                                     x_limit      = NULL,
                                     desfases_y   = NULL) {

  N      <- length(lista_ec)
  n_tips <- Ntip(filogenia)
  edges  <- filogenia$edge
  E      <- nrow(edges)

  cex_aj   <- adjust_cex(filogenia)
  xlim_arg <- if (is.null(x_limit)) NULL else c(0, x_limit)

  plot(filogenia,
       type           = "phylogram",
       edge.color     = "transparent",
       tip.color      = "black",
       edge.width     = grosor,
       cex            = cex_aj,
       label.offset   = label_offset,
       no.margin      = TRUE,
       show.tip.label = TRUE,
       font           = 3L,
       x.lim          = xlim_arg)

  pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  xx <- pp$xx
  yy <- pp$yy

  if (N == 1L) {
    desfases_y <- 0
    desfases_x <- 0
  } else {
    if (is.null(desfases_y)) {
      desfases_y <- seq(-0.08, 0.08, length.out = N)
    }

    usr <- par("usr")
    pin <- par("pin")

    escala_x <- (usr[2] - usr[1]) / pin[1]
    escala_y <- (usr[4] - usr[3]) / pin[2]

    desfases_x <- desfases_y * (escala_x / escala_y)
  }

  nodos_internos <- unique(edges[edges[, 1] > n_tips, 1L])

  nodo_x       <- xx[nodos_internos]
  idx_entrada  <- match(nodos_internos, edges[, 2L])
  hijos_de     <- lapply(nodos_internos, function(nd) edges[edges[, 1L] == nd, 2L])

  padre_idx <- edges[, 1L]
  hijo_idx  <- edges[, 2L]

  for (k in seq_len(N)) {

    ec <- lista_ec[[k]]
    dy <- desfases_y[k]
    dx <- desfases_x[k]

    # Horizontal strokes (vectorized)
    segments(x0   = xx[padre_idx] + dx,
             y0   = yy[hijo_idx]  + dy,
             x1   = xx[hijo_idx]  + dx,
             y1   = yy[hijo_idx]  + dy,
             col  = ec,
             lwd  = grosor,
             lend = 1)

    # Vertical strokes (one segment per internal node)
    for (m in seq_along(nodos_internos)) {

      ie         <- idx_entrada[m]
      color_nodo <- if (!is.na(ie)) ec[ie] else ec[which(edges[, 1L] == nodos_internos[m])[1L]]

      hijos  <- hijos_de[[m]]
      y_hijos <- yy[hijos]
      y_min  <- min(y_hijos)
      y_max  <- max(y_hijos)

      if (y_min != y_max) {
        segments(x0   = nodo_x[m] + dx,
                 y0   = y_min     + dy,
                 x1   = nodo_x[m] + dx,
                 y1   = y_max     + dy,
                 col  = color_nodo,
                 lwd  = grosor,
                 lend = 1)
      }
    }
  }

  invisible(NULL)
}


# ------------------------------------------------------------------------------
# 3.3  Branch painting function
# ------------------------------------------------------------------------------

#' Paints the phylogeny (phylogram) with the colors of a single ancestral mapping.
#'
#' Delegates to `dibujar_filograma_canvas()` with a list of a single vector.
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
  tipo_arbol <- config$tipo_arbol

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
      lwd           = config$grosor,
      legend_labels = ley_data$labels,
      legend_colors = ley_data$colors,
      legend_corner = config$legend_corner %||% "bottomleft",
      legend_title  = titulo_leyenda
    )
  }

  # SCREEN RENDER — advanced visual engine for ANY topology
  plot_multimapr_screen(
    tree           = filogenia,
    color_list     = list(edge_colors),
    type           = tipo_arbol,
    lwd            = config$grosor,
    legend_labels  = if (!is.null(colores_estado)) names(colores_estado) else NULL,
    legend_colors  = if (!is.null(colores_estado)) unname(colores_estado) else NULL,
    legend_corner  = config$legend_corner %||% "bottomleft",
    legend_title   = titulo_leyenda
  )

  invisible(NULL)
}


#' FULL ANCESTRAL RECONSTRUCTION — orchestrates 3.1 and 3.3
#'
#' Dispatches to the appropriate render function based on config$funcion_multi:
#'   1 -> branch superimposition only (plot_superimposed_characters /
#'         plot_ancestral_branches for single character)
#'   2 -> ancestral reconstruction on branches PLUS simple terminal figures
#'         (plot_ancestral_with_terminals)
#'
#' @param filogenia  phylo object.
#' @param config     Configuration list from setup_mapping_config().
#' @return Invisible NULL.
plot_ancestral_reconstruction <- function(filogenia, config) {
  datos_ord       <- config$datos_ord
  caracteres      <- config$caracteres
  colores_por_car <- config$colores_por_caracter

  resolver_color <- function(valor, colores_estado) {
    valor <- as.character(valor)
    if (is.na(valor) || valor == "") return("gray70")
    if (valor %in% names(colores_estado)) return(colores_estado[[valor]])
    return("gray70")
  }

  # funcion_multi == 2: ancestral reconstruction + terminal figures
  if (!is.null(config$funcion_multi) && config$funcion_multi == 2) {
    plot_ancestral_with_terminals(filogenia, config)
    return(invisible(NULL))
  }

  # funcion_multi == 1 (or undefined): branch superimposition only
  if (length(caracteres) == 1) {
    car            <- caracteres[1]
    colores_estado <- colores_por_car[[car]]
    tip_colors     <- sapply(datos_ord[[car]], resolver_color,
                             colores_estado = colores_estado)

    edge_colors <- init_edge_colors(filogenia)
    edge_colors <- apply_ancestral_algorithm(filogenia, tip_colors, edge_colors, config)

    plot_ancestral_branches(filogenia, edge_colors, config,
                            titulo_leyenda = car,
                            colores_estado = colores_estado)

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
#'             (using the algorithm selected in config$algoritmo).
#'   Layer 2 — Colored figures (circles or squares) at the tips showing the
#'             observed tip states (equivalent to the simple mapping squares).
#'
#' The function uses plot_multimapr_screen (and export_multimapr_tree for
#' file output) to draw the ancestral layer, then overlays tip figures
#' with par(new = TRUE) / points() on the same coordinate system.
#'
#' @param filogenia  phylo object.
#' @param config     Configuration list from setup_mapping_config().
#'                   Must contain: datos_ord, caracteres, colores_por_caracter,
#'                   algoritmo, tipo_arbol, grosor, pch_figura, tam_figura,
#'                   and optionally export_filename / export_format.
#' @return Invisible NULL. Draws on the active graphics device.
plot_ancestral_with_terminals <- function(filogenia, config) {

  datos_ord       <- config$datos_ord
  caracteres      <- config$caracteres
  colores_por_car <- config$colores_por_caracter
  tipo_arbol      <- config$tipo_arbol
  pch_fig         <- config$pch_figura %||% 22L   # 22 = square
  tam_fig         <- config$tam_figura %||% 1
  exportar        <- !is.null(config$export_filename)
  fn_export       <- config$export_filename
  n_tips          <- Ntip(filogenia)

  resolver_color <- function(valor, colores_estado) {
    valor <- as.character(valor)
    if (is.na(valor) || valor == "") return("gray70")
    if (valor %in% names(colores_estado)) return(colores_estado[[valor]])
    return("gray70")
  }

  cat("\n=== Generating plot (Ancestral Reconstruction + Terminal Figures) ===\n")

  # ── Step 1: compute ancestral edge colors for each character ────────────────
  lista_ec <- lapply(caracteres, function(car) {
    col_e      <- colores_por_car[[car]]
    tip_colors <- sapply(datos_ord[[car]], resolver_color, colores_estado = col_e)
    ec         <- init_edge_colors(filogenia)
    ec         <- apply_ancestral_algorithm(filogenia, tip_colors, ec, config)
    ec
  })

  # ── Step 2: build legend data (branches block, grouped by character) ─────────
  ley_data <- .build_legend_data(colores_por_car, caracteres)

  # ── Helper: overlay tip figures on an already-rendered phylogeny ────────────
  # Reads tip coordinates from .PlotPhyloEnv after the base plot is drawn.
  # ── Helper: overlay tip figures ──────────────────────────────────────────────
  # pp, cex_aj, lbl_off: when called from overlay_fn (export) they are
  # received from the export_multimapr_tree engine, which already computed them
  # with the correct scale factors. When called from screen (NULL) they are
  # read from .PlotPhyloEnv and adjust_cex() respectively.
  .superponer_figuras_terminales <- function(pp       = NULL,
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
        col_tips     <- sapply(datos_ord[[caracteres[i]]], resolver_color,
                               colores_estado = colores_por_car[[caracteres[i]]])
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
      for (j in seq_len(n_tips)) {
        ang_j  <- angulos[j]
        deg_j  <- ang_j * 180 / pi
        lado_d <- cos(ang_j) >= 0
        srt_j  <- if (lado_d) deg_j else deg_j + 180
        adj_j  <- if (lado_d) c(0, 0.5) else c(1, 0.5)
        text(radio_nombres * cos(ang_j),
             radio_nombres * sin(ang_j),
             labels = filogenia$tip.label[j],
             adj    = adj_j,
             cex    = cex_aj,
             srt    = srt_j,
             font   = 3L)
      }

    } else {
      # Phylogram / cladogram: figures placed to the right of each tip label.
      # lbl_off: from the export engine if available, or proportional to x_max
      max_tip_x <- max(xx_tip)
      if (is.null(lbl_off))
        lbl_off <- config$label_offset %||% (max_tip_x * 0.025)

      max_lbl_w  <- max(strwidth(filogenia$tip.label, cex = cex_aj))
      sep_col    <- strwidth("M", cex = cex_aj) * 2.5
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
        col_tips <- sapply(datos_ord[[caracteres[i]]], resolver_color,
                           colores_estado = colores_por_car[[caracteres[i]]])
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
      lwd             = config$grosor,
      width           = config$width,
      height          = config$height,
      offset_range    = config$rango_desfase %||% 0.1,
      legend_by_char  = ley_data$by_char,
      legend_corner   = config$legend_corner %||% "bottomleft",
      # For fan, hide engine labels and redraw them in the overlay,
      # further out, after all figures.
      hide_fan_labels = (tipo_arbol == "fan"),
      overlay_fn      = function(pp, cex_aj, label_offset_aj,
                                 R_tips = NULL, gap_u = NULL) {
        .superponer_figuras_terminales(pp      = pp,
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
    lwd             = config$grosor,
    offset_range    = config$rango_desfase %||% 0.1,
    legend_by_char  = ley_data$by_char,
    legend_corner   = config$legend_corner %||% "bottomleft",
    hide_fan_labels = (tipo_arbol == "fan"),
    overlay_fn      = function(pp, cex_aj, label_offset_aj,
                               R_tips = NULL, gap_u = NULL) {
      .superponer_figuras_terminales(pp      = pp,
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

#' Draws an offset cladogram for the second character onward
#'
#' @param tree          \code{phylo} object.
#' @param base_coords   List with \code{xx} and \code{yy} nodal coordinate vectors
#'                      from \code{ape::plotPhyloCoor()} or equivalent.
#' @param edge_colors   Character vector of colors, one per edge.
#' @param grosor2       Branch width (line width) for the offset layer.
#' @param alpha2        Alpha transparency applied to all edge colors (0-1).
#' @param offset_clado  Displacement in data units applied to each segment endpoint
#'                      to separate the cladogram layer from the base layer.
#' @return Invisible NULL. Draws on the active graphics device.
#' @keywords internal
plot_cladogram_split <- function(tree, base_coords, edge_colors,
                                 grosor2, alpha2, offset_clado) {
  temp_tree <- tree
  temp_tree$edge.color <- sapply(edge_colors,
                                 function(c) adjustcolor(c, alpha.f = alpha2))
  xx <- base_coords$xx
  yy <- base_coords$yy

  for (j in seq_len(nrow(temp_tree$edge))) {
    n1 <- temp_tree$edge[j, 1];  n2 <- temp_tree$edge[j, 2]
    x0 <- xx[n1]; y0 <- yy[n1]
    x1 <- xx[n2]; y1 <- yy[n2]
    dx <- x1 - x0;  dy <- y1 - y0
    slope <- if (dx != 0) dy / dx else Inf

    if (slope < 0) {
      x0o <- x0 + offset_clado; y0o <- y0 + offset_clado
      x1o <- x1 + offset_clado; y1o <- y1 + offset_clado
    } else if (slope > 0) {
      x0o <- x0 - offset_clado; y0o <- y0 + offset_clado
      x1o <- x1 - offset_clado; y1o <- y1 + offset_clado
    } else {
      if (dx == 0) {
        x0o <- x0 - offset_clado; y0o <- y0
        x1o <- x1 - offset_clado; y1o <- y1
      } else {
        x0o <- x0; y0o <- y0 + offset_clado
        x1o <- x1; y1o <- y1 + offset_clado
      }
    }
    segments(x0o, y0o, x1o, y1o,
             col = temp_tree$edge.color[j],
             lwd = grosor2, lend = 1, ljoin = 1)
  }
}


#' Replots the fan tree with slight rotation and displacement for additional layers
#'
#' @param filogenia       \code{phylo} object.
#' @param edge_colors     Character vector of colors, one per edge.
#' @param grosor          Branch width (line width).
#' @param alpha           Alpha transparency applied to all edge colors (0-1).
#' @param escala          Global scale factor for the fan radius.
#' @param rotacion        Rotation offset in radians applied to the whole fan.
#' @param offset_x        Horizontal displacement in data units (default \code{0}).
#' @param offset_y        Vertical displacement in data units (default \code{0}).
#' @param cex_aj          Font size for tip labels (default \code{0.5}).
#' @param escala_relativa Additional scale factor multiplied with \code{escala} (default \code{1}).
#' @return Invisible NULL. Draws on the active graphics device.
#' @keywords internal
plot_fan_capa <- function(filogenia, edge_colors, grosor, alpha,
                          escala, rotacion, offset_x = 0, offset_y = 0,
                          cex_aj = 0.5, escala_relativa = 1.0) {

  cols <- adjustcolor(edge_colors, alpha.f = alpha)
  esc_efectiva <- escala * escala_relativa

  par(new = TRUE)
  plot.phylo(filogenia,
             type           = "fan",
             edge.color     = cols,
             edge.width     = grosor,
             show.tip.label = FALSE,
             no.margin      = TRUE,
             rotate.tree    = rotacion,
             x.lim = c(-esc_efectiva + offset_x,  esc_efectiva + offset_x),
             y.lim = c(-esc_efectiva + offset_y,  esc_efectiva + offset_y))
}


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

  datos_ord       <- config$datos_ord
  caracteres      <- config$caracteres
  colores_por_car <- config$colores_por_caracter
  tipo_arbol      <- config$tipo_arbol
  grosor1         <- config$grosor %||% 2

  resolver_color <- function(valor, colores_estado) {
    valor <- as.character(valor)
    if (is.na(valor) || valor == "") return("gray70")
    if (valor %in% names(colores_estado)) return(colores_estado[[valor]])
    return("gray70")
  }

  lista_ec <- lapply(caracteres, function(car) {
    col_e      <- colores_por_car[[car]]
    tip_colors <- sapply(datos_ord[[car]], resolver_color, colores_estado = col_e)
    ec         <- init_edge_colors(filogenia)
    ec         <- apply_ancestral_algorithm(filogenia, tip_colors, ec, config)
    ec
  })

  exportar  <- !is.null(config$export_filename)
  fn_export <- config$export_filename
  ley_data  <- .build_legend_data(colores_por_car, caracteres)

  if (exportar) {
    cat("    \u2192 Exporting to:", paste0(fn_export, ".", config$export_format %||% "png"), "\n")
    export_multimapr_tree(
      tree           = filogenia,
      color_list     = lista_ec,
      filename       = fn_export,
      type           = tipo_arbol,
      format         = config$export_format %||% "png",
      lwd            = grosor1,
      offset_range   = config$rango_desfase %||% 0.1,
      legend_by_char = ley_data$by_char,
      legend_corner  = config$legend_corner %||% "bottomleft"
    )
  }

  # SCREEN RENDER — integrated advanced visual engine
  plot_multimapr_screen(
    tree           = filogenia,
    color_list     = lista_ec,
    type           = tipo_arbol,
    lwd            = grosor1,
    offset_range   = config$rango_desfase %||% 0.1,
    legend_by_char = ley_data$by_char,
    legend_corner  = config$legend_corner %||% "bottomleft"
  )

  invisible(NULL)
}


#' MULTI-CHARACTER FIGURE — dispatched when funcion_multi == 2
#'
#' Wrapper retained for backward compatibility with callers that still
#' reference the old name. Internally calls plot_superimposed_characters.
#'
#' @param filogenia  phylo object.
#' @param config     Configuration list from setup_mapping_config().
#' @return Invisible NULL.
plot_multicharacter_figure <- function(filogenia, config) {
  plot_superimposed_characters(filogenia, config)
}