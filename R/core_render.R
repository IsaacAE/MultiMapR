# Rendering controllers: simple mapping (terminal figures) and ancestral
# reconstruction (branch coloring), for the screen and for file export.


# Simple mapping

#' Horizontal distance between the columns of terminal figures
#' @noRd
.figure_col_sep <- function(cex_txt, fig_size) {
  x_per_in <- diff(par("usr")[1:2]) / par("pin")[1L]
  diam_in  <- 0.375 * fig_size * par("cin")[2L]
  max(strwidth("M", cex = cex_txt) * 2.5, diam_in * 1.35 * x_per_in)
}

#' Resolves the tip figure color given its state value
#' @noRd
resolve_tip_color <- function(value, color_map, color_na = "gray70") {
  value <- as.character(value)
  if (is.na(value) || value == "") return(color_na)
  if (value %in% names(color_map)) return(color_map[[value]])
  return(color_na)
}


#' Simple mapping: neutral tree with colored figures (table or rings) at the tips
#' @noRd
plot_simple_mapping <- function(tree, config) {

  aligned_data       <- config$aligned_data
  characters      <- config$characters
  colors_by_char <- config$colors_by_character
  tree_type      <- config$tree_type
  pch_fig         <- config$figure_pch
  fig_size         <- config$figure_size
  simple_mode     <- config$simple_mode %||% "figures"   # "figures" | "tip_color"
  cex_adj          <- adjust_cex(tree, config = config)
  n_tips          <- Ntip(tree)
  n_chars           <- length(characters)
  show_lbl        <- config$show_labels %||% TRUE
  cex_lbl         <- cex_adj * (config$label_size %||% 1)
  col_lbl         <- config$label_color %||% "black"
  show_leg        <- config$show_legend %||% TRUE

  resolve_color <- function(value, state_colors, char) {
    .ancestral_tip_color(value, state_colors, config, char)
  }

  is_exporting  <- !is.null(config$export_filename)
  fn_export <- config$export_filename

  cat("\n=== Generating plot (Simple Mapping) ===\n")
  if (is_exporting) cat("    \u2192 Exporting to:",
                    paste0(fn_export, ".", config$export_format %||% "png"), "\n")

  .draw_simple_legend <- function(legend_pos) {
    if (!show_leg) return(invisible(NULL))
    blocks <- .simple_legend_blocks()
    lay    <- .legend_layout(blocks, legend_pos, config$legend_layout %||% "vertical",
                             cex = cex_adj, pt_cex = fig_size, pch = pch_fig)
    old_xpd <- par("xpd"); par(xpd = TRUE)
    .legend_draw_blocks(lay, blocks, cex = cex_adj, pt_cex = fig_size, pch = pch_fig,
                        filled = pch_fig %in% 21:25)
    par(xpd = old_xpd)
    invisible(NULL)
  }

  .simple_legend_blocks <- function() {
    lapply(seq_along(characters), function(i) {
      char   <- characters[i]
      state_cols <- c(colors_by_char[[char]], .missing_legend_entries(config, char))
      list(title  = if (tree_type == "fan") paste0("C", i, " (", char, ")") else char,
           labels = names(state_cols), colors = unname(state_cols))
    })
  }

  # xlim_extra / ylim_extra expand the viewport for the table and legend.
  .render_simple_mapping <- function(xlim_extra = NULL, ylim_extra = NULL) {

    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    par(lend = config$tip_end %||% "round")

    legend_position <- if (!is.null(config$legend_corner)) config$legend_corner
    else if (tree_type == "fan") "topleft" else "topright"

    # Fan: no.margin plot with a symmetric x.lim computed after a first plot, so
    # rings, labels and legend fit without relying on margins.
    if (tree_type == "fan") {

      # Tree coordinates come from the off-screen probe (PHASE 1 below).
      obj <- .fan_coords_probe

      xx_tip <- obj$xx[seq_len(n_tips)]
      yy_tip <- obj$yy[seq_len(n_tips)]

      angles        <- atan2(yy_tip, xx_tip)
      tip_radius     <- sqrt(xx_tip^2 + yy_tip^2)
      max_tip_radius <- max(tip_radius)

      if (simple_mode == "tip_color") {
        # tip_color mode: edges leading to a tip take its state color; others stay gray.
        char1     <- characters[1L]
        col_tips1 <- sapply(aligned_data[[char1]], resolve_color,
                            state_colors = colors_by_char[[char1]], char = char1)

        edge_col <- rep("gray70", nrow(tree$edge))
        for (e in seq_len(nrow(tree$edge))) {
          child <- tree$edge[e, 2L]
          if (child <= n_tips)
            edge_col[e] <- col_tips1[child]
        }

        # Same label geometry as the export fan: gap_u = strwidth("m") * 1.5 and
        # radius_labels = max_tip_radius + gap_u; total radius adds label width + padding.

        # First plot with show.tip.label = FALSE to get definitive tip coordinates
        # and then measure strwidth in the correct coordinate system.
        max_nchar_lbl <- max(nchar(tree$tip.label))

        # Rough xlim to open the canvas (will be tightened after measuring)
        rough_lbl_w <- max_tip_radius * 0.018 * cex_adj * max_nchar_lbl
        rough_total <- max_tip_radius * 1.06 + rough_lbl_w + max_tip_radius * 0.05
        xlim_fan    <- c(-rough_total, rough_total)

        plot(tree,
             type           = "fan",
             cex            = cex_adj,
             label.offset   = 0,           # labels drawn manually below
             edge.width     = config$branch_width,
             edge.color     = edge_col,
             tip.color      = "transparent",
             show.tip.label = FALSE,
             no.margin      = TRUE,
             font           = 3L,
             x.lim          = xlim_fan)

        obj2       <- get("last_plot.phylo", envir = .PlotPhyloEnv)
        xx_tip2    <- obj2$xx[seq_len(n_tips)]
        yy_tip2    <- obj2$yy[seq_len(n_tips)]
        angles2   <- atan2(yy_tip2, xx_tip2)
        tip_radius2 <- sqrt(xx_tip2^2 + yy_tip2^2)
        R_tips2     <- max(tip_radius2)

        gap_u2        <- strwidth("m", cex = cex_adj) * 1.5
        radius_labels2 <- R_tips2 + gap_u2

        old_xpd <- par("xpd"); par(xpd = TRUE)
        for (j in seq_len(n_tips)[show_lbl]) {
          ang_j  <- angles2[j]; deg_j <- ang_j * 180 / pi
          is_right <- cos(ang_j) >= 0
          srt_j  <- if (is_right) deg_j else deg_j + 180
          adj_j  <- if (is_right) c(0, 0.5) else c(1, 0.5)
          text(radius_labels2 * cos(ang_j), radius_labels2 * sin(ang_j),
               labels = tree$tip.label[j], adj = adj_j, cex = cex_lbl,
               col = col_lbl, srt = srt_j, font = 3L)
        }
        par(xpd = old_xpd)

        .draw_simple_legend(legend_position)

      } else {
        # figures mode: concentric symbol rings
        radius_base       <- max_tip_radius * 1.06
        radius_step <- max_tip_radius * 0.10
        radius_last     <- radius_base + (n_chars - 1L) * radius_step

        # Space between last ring and label
        gap_label         <- radius_step * 0.5
        radius_labels    <- radius_last + gap_label

        # Step 3: estimate total radius including tip labels
        max_nchar_lbl <- max(nchar(tree$tip.label))
        char_w_u      <- max_tip_radius * 0.018 * cex_adj
        lbl_w_u       <- max_nchar_lbl * char_w_u
        radius_total   <- radius_labels + lbl_w_u + max_tip_radius * 0.05

        # Step 4: replot with symmetric x.lim that contains everything
        xlim_fan <- c(-radius_total, radius_total)

        plot(tree,
             type           = "fan",
             cex            = cex_adj,
             label.offset   = config$label_offset,
             edge.width     = config$branch_width,
             tip.color      = "transparent",
             show.tip.label = FALSE,
             no.margin      = TRUE,
             x.lim          = xlim_fan)

        obj    <- get("last_plot.phylo", envir = .PlotPhyloEnv)
        xx_tip <- obj$xx[seq_len(n_tips)]
        yy_tip <- obj$yy[seq_len(n_tips)]
        angles <- atan2(yy_tip, xx_tip)

        old_xpd <- par("xpd"); par(xpd = TRUE)
        for (i in seq_len(n_chars)) {
          col_tips     <- sapply(aligned_data[[characters[i]]], resolve_color,
                                 state_colors = colors_by_char[[characters[i]]],
                                 char = characters[i])
          radius_current <- radius_base + (i - 1L) * radius_step

          points(radius_current * cos(angles),
                 radius_current * sin(angles),
                 pch = pch_fig, bg = col_tips, col = "black", cex = fig_size)

          gap_tan <- radius_step * 0.55
          for (j in seq_len(n_tips)) {
            ang_j   <- angles[j]
            deg_j   <- ang_j * 180 / pi
            is_right  <- cos(ang_j) >= 0
            srt_j   <- if (is_right) deg_j else deg_j + 180
            tan_x   <- -sin(ang_j)
            tan_y   <-  cos(ang_j)
            cx      <- radius_current * cos(ang_j) + tan_x * gap_tan
            cy      <- radius_current * sin(ang_j) + tan_y * gap_tan
            text(cx, cy,
                 labels = paste0("C", i),
                 adj    = c(0.5, 0.5),
                 cex    = cex_adj * 0.65,
                 font   = 2L,
                 srt    = srt_j)
          }
        }
        par(xpd = old_xpd)

        # Tip labels — parallel to radius
        old_xpd <- par("xpd"); par(xpd = TRUE)
        for (j in seq_len(n_tips)[show_lbl]) {
          ang_j  <- angles[j]; deg_j <- ang_j * 180 / pi
          is_right <- cos(ang_j) >= 0
          srt_j  <- if (is_right) deg_j else deg_j + 180
          adj_j  <- if (is_right) c(0, 0.5) else c(1, 0.5)
          text(radius_labels * cos(ang_j), radius_labels * sin(ang_j),
               labels = tree$tip.label[j], adj = adj_j, cex = cex_lbl, srt = srt_j,
               col    = col_lbl, font   = 3L)
        }
        par(xpd = old_xpd)

        .draw_simple_legend(legend_position)

      }

      # PHYLOGRAM / CLADOGRAM
    } else {

      par(mar = c(1, 1, 2, 1), xpd = TRUE)

      if (simple_mode == "tip_color") {
        # tip_color mode: edges leading to a tip take its state color; others stay gray.
        char1     <- characters[1L]
        col_tips1 <- sapply(aligned_data[[char1]], resolve_color,
                            state_colors = colors_by_char[[char1]], char = char1)

        edge_col <- rep("gray70", nrow(tree$edge))
        for (e in seq_len(nrow(tree$edge))) {
          child <- tree$edge[e, 2L]
          if (child <= n_tips)
            edge_col[e] <- col_tips1[child]
        }

        plot(tree,
             type           = tree_type,
             cex            = cex_lbl,
             label.offset   = config$label_offset,
             edge.width     = config$branch_width,
             edge.color     = edge_col,
             tip.color      = col_lbl,
             show.tip.label = show_lbl,
             font           = 3L,
             x.lim          = xlim_extra,
             y.lim          = ylim_extra)

        .draw_simple_legend(legend_position)

      } else {
        # MODE: figures (default)
        plot(tree,
             type           = tree_type,
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

        sep_col    <- .figure_col_sep(cex_adj, fig_size)
        max_tip_x  <- max(xx_tip)
        max_lbl_w  <- if (show_lbl) max(strwidth(tree$tip.label, cex = cex_lbl)) else 0
        start_x    <- max_tip_x + config$label_offset + max_lbl_w + max_tip_x * 0.02
        x_columns <- start_x + seq(0, n_chars - 1L) * sep_col

        y_range  <- diff(range(yy_tip))
        y_header <- max(yy_tip) + y_range * 0.04
        old_xpd <- par("xpd"); par(xpd = NA)
        text(x_columns, y_header,
             labels = characters,
             cex    = cex_adj * 0.9,
             srt    = 45,
             adj    = c(0, 0.5),
             font   = 2L)
        par(xpd = old_xpd)

        fig_x <- fig_y <- fig_col <- NULL
        for (i in seq_along(characters)) {
          col_tips <- sapply(aligned_data[[characters[i]]], resolve_color,
                             state_colors = colors_by_char[[characters[i]]],
                             char = characters[i])
          fig_x  <- c(fig_x,   rep(x_columns[i], n_tips))
          fig_y   <- c(fig_y,   yy_tip)
          fig_col <- c(fig_col, col_tips)
        }
        points(fig_x, fig_y, pch = pch_fig, bg = fig_col, col = "black", cex = fig_size)

        .draw_simple_legend(legend_position)
      }
    }

    invisible(NULL)
  }


  # PHASE 1: off-screen probe
  .fan_coords_probe <- NULL
  computed_lim <- list(x = NULL, y = NULL)

  if (tree_type == "fan") {
    prev_dev <- dev.cur()
    pdf(NULL, width = 7, height = 7)
    .fan_coords_probe <- tryCatch({
      plot(tree,
           type           = "fan",
           cex            = cex_adj,
           label.offset   = config$label_offset,
           edge.width     = config$branch_width,
           tip.color      = "transparent",
           show.tip.label = FALSE,
           no.margin      = TRUE)
      get("last_plot.phylo", envir = .PlotPhyloEnv)
    }, error = function(e) NULL,
    finally = .close_device_restore(prev_dev))
  }

  if (tree_type != "fan" && simple_mode == "figures") {
    n_tips_tmp    <- n_tips
    default_height  <- n_tips_tmp * 0.25 + 2
    default_width <- 12
    height_in_tmp   <- if (!is.null(config$height)) config$height else default_height
    width_in_tmp  <- if (!is.null(config$width))  config$width  else default_width
    # On screen, probe with the size of the device that will be drawn on
    if (!is_exporting && dev.cur() > 1L) {
      din <- dev.size("in")
      width_in_tmp <- din[1L]; height_in_tmp <- din[2L]
    }

    prev_dev <- dev.cur()
    pdf(NULL, width = width_in_tmp, height = height_in_tmp)
    computed_lim <- tryCatch({

      par(mar = c(1, 1, 2, 1), xpd = FALSE)
      plot(tree,
           type           = tree_type,
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

      sep_col_s  <- .figure_col_sep(cex_adj, fig_size)
      max_tip_s  <- max(xx_s)
      max_lbl_s  <- if (show_lbl) max(strwidth(tree$tip.label, cex = cex_lbl)) else 0
      start_x_s  <- max_tip_s + config$label_offset + max_lbl_s + max_tip_s * 0.02
      x_cols_s   <- start_x_s + seq(0, n_chars - 1L) * sep_col_s

      longest_char     <- characters[which.max(nchar(characters))]
      hdr_w_horiz   <- strwidth(longest_char, cex = cex_adj * 0.9) * cos(pi / 4)
      x_max_table   <- max(x_cols_s) + hdr_w_horiz

      legend_pos_s  <- if (!is.null(config$legend_corner)) config$legend_corner else "topright"
      on_right_s <- grepl("right", legend_pos_s)

      # Vertical legend: reserve its width beside the table. Horizontal legend:
      # reserve its height above / below the tree (it spans the bottom or top).
      horiz <- identical(config$legend_layout, "horizontal")
      lay   <- if (show_leg)
        .legend_layout(.simple_legend_blocks(), legend_pos_s,
                       config$legend_layout %||% "vertical",
                       cex = cex_adj, pt_cex = fig_size, pch = pch_fig) else list(w = 0, h = 0)

      x_max_total <- x_max_table + if (horiz) 0 else lay$w * 1.1
      ylim <- if (horiz && lay$h > 0) {
        if (grepl("top", legend_pos_s)) c(usr_s[3L], usr_s[4L] + lay$h * 1.1)
        else                         c(usr_s[3L] - lay$h * 1.1, usr_s[4L])
      }

      list(x = if (x_max_total > usr_s[2L]) c(usr_s[1L], x_max_total),
           y = ylim)

    }, error = function(e) list(x = NULL, y = NULL),
    finally = .close_device_restore(prev_dev))
  }

  # PHASE 2: export or draw on screen
  if (is_exporting) {
    .export_device(fn_export, config$export_format %||% "png",
                   tree, tree_type,
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


# Ancestral reconstruction

# Algorithms

#' Initializes the edge color vector in gray
#' @noRd
init_edge_colors <- function(tree) {
  ec <- rep("gray70", nrow(tree$edge))
  names(ec) <- paste0(tree$edge[, 1], "-", tree$edge[, 2])
  ec
}


#' Resolves the single color used for missing ("?"), inapplicable ("-") and ambiguous states under the Fitch algorithm
#' @noRd
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
#' @noRd
.add_fitch_ambig_to_legend <- function(state_colors, config) {
  if (!identical(config$algorithm, "fitch")) return(state_colors)
  c(state_colors, c("Ambiguous / Missing" = .get_fitch_ambig_color(config)))
}


#' Normalizes the colors chosen for missing ("?") and inapplicable ("-") tips
#' @noRd
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
#' @noRd
.ancestral_tip_color <- function(value, state_colors, config, char) {
  value <- as.character(value)
  if (!is.na(value) && value %in% names(state_colors)) return(state_colors[[value]])
  if (is.na(value) || value %in% c("", "?", "-")) {
    if (identical(config$algorithm, "fitch")) return(.get_fitch_ambig_color(config))
    sym <- if (is.na(value) || value == "") "?" else value
    mc  <- config$missing_colors[[char]]
    if (!is.null(mc) && sym %in% names(mc)) return(unname(mc[[sym]]))
  }
  "gray70"
}

#' Legend entries for the missing / inapplicable colors of one character
#' @noRd
.missing_legend_entries <- function(config, char) {
  if (identical(config$algorithm, "fitch")) return(NULL)
  mc <- config$missing_colors[[char]]
  if (is.null(mc)) return(NULL)
  vals    <- as.character(config$aligned_data[[char]])
  present <- c("?" = any(is.na(vals) | vals %in% c("", "?")), "-" = any(vals %in% "-"))
  syms    <- intersect(names(mc), names(present)[present])
  syms    <- setdiff(syms, names(config$colors_by_character[[char]]))
  if (length(syms) == 0L) return(NULL)
  labels  <- c("?" = "Missing (?)", "-" = "Inapplicable (-)")
  stats::setNames(unname(mc[syms]), labels[syms])
}

#' Appends the missing / inapplicable entries to grouped legend data
#' @noRd
.add_missing_to_legend_data <- function(legend_data, config) {
  for (nm in names(legend_data$by_char)) {
    extra <- .missing_legend_entries(config, nm)
    if (length(extra) == 0L) next
    legend_data$by_char[[nm]]$labels <- c(legend_data$by_char[[nm]]$labels, names(extra))
    legend_data$by_char[[nm]]$colors <- c(legend_data$by_char[[nm]]$colors, unname(extra))
  }
  legend_data
}


#' Dispatches to the ancestral reconstruction algorithm selected in config
#' @noRd
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


# Branch painting

#' Paints the phylogeny (phylogram) with the colors of a single ancestral mapping
#' @noRd
plot_ancestral_branches <- function(tree, edge_colors, config,
                                    legend_title = "", state_colors = NULL) {
  is_exporting   <- !is.null(config$export_filename)
  fn_export  <- config$export_filename
  tree_type <- config$tree_type

  cat("\n=== Generating plot (Ancestral Reconstruction) ===\n")

  if (is_exporting) {
    cat("    \u2192 Exporting to:", paste0(fn_export, ".", config$export_format %||% "png"), "\n")

    legend_data <- if (!is.null(state_colors) && length(state_colors) > 0)
      list(labels = names(state_colors), colors = unname(state_colors))
    else
      list(labels = NULL, colors = NULL)

    export_multimapr_tree(
      tree          = tree,
      color_list    = list(edge_colors),
      filename      = fn_export,
      type          = tree_type,
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
       legend_labels = if (config$show_legend %||% TRUE) legend_data$labels,
       legend_colors = legend_data$colors,
       legend_corner = config$legend_corner %||% "bottomleft",
       legend_layout = config$legend_layout %||% "vertical",
       node_states   = isTRUE(config$show_nodes),
       node_shape    = config$node_shape %||% "circle",
       node_size     = config$node_size %||% 1,
       legend_title  = legend_title
    )
  }

  # Screen render (any topology)
  plot_multimapr_screen(
    tree           = tree,
    color_list     = list(edge_colors),
    type           = tree_type,
    lwd            = config$branch_width,
    ladderize      = config$ladderize %||% FALSE,
      terminal_stretch = config$terminal_stretch %||% 1,
      tip_end          = config$tip_end %||% "round",
      branch_angle     = config$branch_angle %||% 45,
      show_labels      = config$show_labels %||% TRUE,
      label_size       = config$label_size %||% 1,
      label_color      = config$label_color %||% "black",
    legend_labels  = if (!is.null(state_colors) && config$show_legend %||% TRUE) names(state_colors),
    legend_colors  = if (!is.null(state_colors)) unname(state_colors) else NULL,
    legend_corner  = config$legend_corner %||% "bottomleft",
    legend_layout  = config$legend_layout %||% "vertical",
    node_states    = isTRUE(config$show_nodes),
    node_shape     = config$node_shape %||% "circle",
    node_size      = config$node_size %||% 1,
    legend_title   = legend_title
  )

  invisible(NULL)
}


#' Ancestral reconstruction drawn on the branches (screen or file)
#' @noRd
plot_ancestral_reconstruction <- function(tree, config) {
  aligned_data       <- config$aligned_data
  characters      <- config$characters
  colors_by_char <- config$colors_by_character

  # multi_function == 2: ancestral reconstruction + terminal figures
  if (!is.null(config$multi_function) && config$multi_function == 2) {
    plot_ancestral_with_terminals(tree, config)
    return(invisible(NULL))
  }

  # multi_function == 1 (or undefined): branch superimposition only
  if (length(characters) == 1) {
    char            <- characters[1]
    state_colors <- colors_by_char[[char]]
    tip_colors     <- sapply(aligned_data[[char]], .ancestral_tip_color,
                             state_colors = state_colors, config = config, char = char)

    edge_colors <- init_edge_colors(tree)
    edge_colors <- apply_ancestral_algorithm(tree, tip_colors, edge_colors, config)

    # For Fitch, append the ambiguity entry to the legend (Fitch-only)
    legend_state_colors <- c(.add_fitch_ambig_to_legend(state_colors, config),
                            .missing_legend_entries(config, char))

    plot_ancestral_branches(tree, edge_colors, config,
                            legend_title = char,
                            state_colors = legend_state_colors)

  } else {
    plot_superimposed_characters(tree, config)
  }
}


# Ancestral reconstruction + terminal figures

#' Ancestral reconstruction on branches PLUS simple terminal figures
#' @noRd
plot_ancestral_with_terminals <- function(tree, config) {

  aligned_data       <- config$aligned_data
  characters      <- config$characters
  colors_by_char <- config$colors_by_character
  tree_type      <- config$tree_type
  pch_fig         <- config$figure_pch %||% 22L   # 22 = square
  fig_size         <- config$figure_size %||% 1
  is_exporting        <- !is.null(config$export_filename)
  fn_export       <- config$export_filename
  n_tips          <- Ntip(tree)

  # Missing ("?") and inapplicable ("-") tips use the same color on figures and
  # branches: the Fitch ambiguity color, or the per-character colors chosen for
  # the depth-weighted majority (see .ancestral_tip_color()).
  resolve_color <- function(value, state_colors, char)
    .ancestral_tip_color(value, state_colors, config, char)

  cat("\n=== Generating plot (Ancestral Reconstruction + Terminal Figures) ===\n")

  # Step 1: compute ancestral edge colors for each character
  edge_color_list <- lapply(characters, function(char) {
    state_cols      <- colors_by_char[[char]]
    tip_colors <- sapply(aligned_data[[char]], resolve_color, state_colors = state_cols, char = char)
    ec         <- init_edge_colors(tree)
    ec         <- apply_ancestral_algorithm(tree, tip_colors, ec, config)
    ec
  })

  # Step 2: build legend data (branches block, grouped by character)
  legend_data <- .build_legend_data(colors_by_char, characters)

  # For Fitch: append the ambiguity entry to every character block (Fitch-only)
  if (identical(config$algorithm, "fitch")) {
    ambig_col <- .get_fitch_ambig_color(config)
    ambig_label <- "Ambiguous / Missing"
    for (nm in names(legend_data$by_char)) {
      if (!ambig_label %in% legend_data$by_char[[nm]]$labels) {
        legend_data$by_char[[nm]]$labels <- c(legend_data$by_char[[nm]]$labels, ambig_label)
        legend_data$by_char[[nm]]$colors <- c(legend_data$by_char[[nm]]$colors, ambig_col)
      }
    }
  }
  # Majority: per-character missing / inapplicable colors (when chosen)
  legend_data <- .add_missing_to_legend_data(legend_data, config)

  # Overlays tip figures on an already-rendered tree, reading tip coordinates from
  # .PlotPhyloEnv. pp, cex_adj and lbl_off come from the export engine when called via
  # overlay_fn; otherwise they are recomputed (.PlotPhyloEnv, adjust_cex()).
  .overlay_tip_figures <- function(pp       = NULL,
                                             cex_adj   = NULL,
                                             lbl_off  = NULL,
                                             R_tips   = NULL,
                                             gap_u    = NULL) {
    # gap_u_fan: local name to avoid conflict with R_tips/gap_u parameters
    gap_u_fan <- gap_u

    if (is.null(pp))
      pp <- get("last_plot.phylo",
                envir = get(".PlotPhyloEnv", envir = asNamespace("ape")))

    xx_tip <- pp$xx[seq_len(n_tips)]
    yy_tip <- pp$yy[seq_len(n_tips)]

    # cex_adj: use the export engine's value if passed, otherwise recalculate
    if (is.null(cex_adj))
      cex_adj <- max(1 / (1 + n_tips / 50), 0.2)

    old_xpd <- par("xpd"); par(xpd = TRUE)

    if (tree_type == "fan") {
      # Fan: same geometry as plot_simple_mapping(): first ring at 1.06 * max tip radius,
      # ring spacing radius_step = 0.10 * max tip radius, tangential "C_i" offset
      # 0.55 * radius_step, labels beyond the last ring. R_tips may come from the
      # export engine (overlay_fn) or be computed here.
      angles <- atan2(yy_tip, xx_tip)
      if (is.null(R_tips)) R_tips <- sqrt(xx_tip^2 + yy_tip^2)
      max_tip_radius   <- max(R_tips)

      radius_base       <- max_tip_radius * 1.06
      radius_step <- max_tip_radius * 0.10
      n_chars_ov         <- length(characters)
      radius_last     <- radius_base + (n_chars_ov - 1L) * radius_step
      gap_label         <- radius_step * 0.5
      radius_labels    <- radius_last + gap_label
      gap_tan          <- radius_step * 0.55

      for (i in seq_len(n_chars_ov)) {
        col_tips     <- sapply(aligned_data[[characters[i]]], resolve_color,
                               state_colors = colors_by_char[[characters[i]]],
                               char = characters[i])
        radius_current <- radius_base + (i - 1L) * radius_step

        # Ring i figures
        points(radius_current * cos(angles),
               radius_current * sin(angles),
               pch = pch_fig, bg = col_tips, col = "black", cex = fig_size)

        # Tangential "C_i" label (same as simple mode)
        for (j in seq_len(n_tips)) {
          ang_j  <- angles[j]
          deg_j  <- ang_j * 180 / pi
          is_right <- cos(ang_j) >= 0
          srt_j  <- if (is_right) deg_j else deg_j + 180
          tan_x  <- -sin(ang_j)
          tan_y  <-  cos(ang_j)
          cx     <- radius_current * cos(ang_j) + tan_x * gap_tan
          cy     <- radius_current * sin(ang_j) + tan_y * gap_tan
          text(cx, cy,
               labels = paste0("C", i),
               adj    = c(0.5, 0.5),
               cex    = cex_adj * 0.65,
               font   = 2L,
               srt    = srt_j)
        }
      }

      # Species labels (parallel to radius, beyond the last ring)
      for (j in seq_len(n_tips)[config$show_labels %||% TRUE]) {
        ang_j  <- angles[j]
        deg_j  <- ang_j * 180 / pi
        is_right <- cos(ang_j) >= 0
        srt_j  <- if (is_right) deg_j else deg_j + 180
        adj_j  <- if (is_right) c(0, 0.5) else c(1, 0.5)
        text(radius_labels * cos(ang_j),
             radius_labels * sin(ang_j),
             labels = tree$tip.label[j],
             adj    = adj_j,
             cex    = cex_adj * (config$label_size %||% 1),
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
                      max(strwidth(tree$tip.label,
                                   cex = cex_adj * (config$label_size %||% 1))) else 0
      sep_col    <- .figure_col_sep(cex_adj, fig_size)
      start_x    <- max_tip_x + lbl_off + max_lbl_w + max_tip_x * 0.02
      x_columns <- start_x + seq(0, length(characters) - 1L) * sep_col

      # Column headers rotated 45 degrees (same as simple mode)
      y_range  <- diff(range(yy_tip))
      y_header <- max(yy_tip) + y_range * 0.04
      old_xpd2 <- par("xpd"); par(xpd = NA)
      text(x_columns, y_header,
           labels = characters,
           cex    = cex_adj * 0.9,
           srt    = 45,
           adj    = c(0, 0.5),
           font   = 2L)
      par(xpd = old_xpd2)

      # Figures per column
      for (i in seq_along(characters)) {
        col_tips <- sapply(aligned_data[[characters[i]]], resolve_color,
                           state_colors = colors_by_char[[characters[i]]],
                           char = characters[i])
        points(rep(x_columns[i], n_tips), yy_tip,
               pch = pch_fig, bg = col_tips, col = "black", cex = fig_size)
      }
    }

    par(xpd = old_xpd)
    invisible(NULL)
  }

  # Step 3a: export to file through export_multimapr_tree(); tip figures are injected
  # via overlay_fn, which receives pp, cex_adj and label_offset_adj already scaled for
  # the export device.
  if (is_exporting) {
    cat("    \u2192 Exporting to:", paste0(fn_export, ".", config$export_format %||% "png"), "\n")

    export_multimapr_tree(
      tree            = tree,
      color_list      = edge_color_list,
      filename        = fn_export,
      type            = tree_type,
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
      legend_by_char  = if (config$show_legend %||% TRUE) legend_data$by_char,
      legend_corner   = config$legend_corner %||% "bottomleft",
      legend_layout   = config$legend_layout %||% "vertical",
      node_states     = isTRUE(config$show_nodes),
      node_shape      = config$node_shape %||% "circle",
      node_size       = config$node_size %||% 1,
      # For fan, hide engine labels and redraw them in the overlay,
      # further out, after all figures.
      hide_fan_labels = (tree_type == "fan"),
      header_labels   = characters,
      overlay_fn      = function(pp, cex_adj, label_offset_adj,
                                 R_tips = NULL, gap_u = NULL) {
        .overlay_tip_figures(pp      = pp,
                                       cex_adj  = cex_adj,
                                       lbl_off = label_offset_adj,
                                       R_tips  = R_tips,
                                       gap_u   = gap_u)
      }
    )
  }

  # Step 3b: screen render (always, regardless of export). For fan, hide_fan_labels
  # delegates the labels to the overlay so they sit beyond the figure rings.
  plot_multimapr_screen(
    tree            = tree,
    color_list      = edge_color_list,
    type            = tree_type,
    lwd             = config$branch_width,
    offset_range    = config$offset_range %||% 0.1,
    tip_end         = config$tip_end %||% "round",
    branch_angle    = config$branch_angle %||% 45,
    show_labels     = config$show_labels %||% TRUE,
    label_size      = config$label_size %||% 1,
    label_color     = config$label_color %||% "black",
    legend_by_char  = if (config$show_legend %||% TRUE) legend_data$by_char,
    legend_corner   = config$legend_corner %||% "bottomleft",
    legend_layout   = config$legend_layout %||% "vertical",
    node_states     = isTRUE(config$show_nodes),
    node_shape      = config$node_shape %||% "circle",
    node_size       = config$node_size %||% 1,
    hide_fan_labels = (tree_type == "fan"),
    header_labels   = characters,
    overlay_fn      = function(pp, cex_adj, label_offset_adj,
                               R_tips = NULL, gap_u = NULL) {
      .overlay_tip_figures(pp      = pp,
                                     cex_adj  = cex_adj,
                                     lbl_off = label_offset_adj,
                                     R_tips  = R_tips,
                                     gap_u   = gap_u)
    }
  )

  invisible(NULL)
}


# MULTI-CHARACTER

#' Superimposition of multiple characters on branches (advanced engine)
#' @noRd
plot_superimposed_characters <- function(tree, config,
                                         alpha1           = 1,
                                         alpha2           = 1,
                                         offset_h_map     = 0.02,
                                         offset_v_map     = 0.02,
                                         extra_offset_map = 0.07,
                                         fan_scale       = 0.75,
                                         fan_rotation     = 2,
                                         base_scale      = 1,
                                         legend_spacing   = 0.06,
                                         x_limit          = NULL) {

  aligned_data       <- config$aligned_data
  characters      <- config$characters
  colors_by_char <- config$colors_by_character
  tree_type      <- config$tree_type
  width1         <- config$branch_width %||% 2

  resolve_color <- function(value, state_colors, char)
    .ancestral_tip_color(value, state_colors, config, char)

  edge_color_list <- lapply(characters, function(char) {
    state_cols      <- colors_by_char[[char]]
    tip_colors <- sapply(aligned_data[[char]], resolve_color, state_colors = state_cols, char = char)
    ec         <- init_edge_colors(tree)
    ec         <- apply_ancestral_algorithm(tree, tip_colors, ec, config)
    ec
  })

  is_exporting  <- !is.null(config$export_filename)
  fn_export <- config$export_filename
  legend_data  <- .build_legend_data(colors_by_char, characters)

  # For Fitch: append the ambiguity entry to every character block (Fitch-only)
  if (identical(config$algorithm, "fitch")) {
    ambig_col <- .get_fitch_ambig_color(config)
    ambig_label <- "Ambiguous / Missing"
    for (nm in names(legend_data$by_char)) {
      if (!ambig_label %in% legend_data$by_char[[nm]]$labels) {
        legend_data$by_char[[nm]]$labels <- c(legend_data$by_char[[nm]]$labels, ambig_label)
        legend_data$by_char[[nm]]$colors <- c(legend_data$by_char[[nm]]$colors, ambig_col)
      }
    }
  }
  # Majority: per-character missing / inapplicable colors (when chosen)
  legend_data <- .add_missing_to_legend_data(legend_data, config)

  if (is_exporting) {
    cat("    \u2192 Exporting to:", paste0(fn_export, ".", config$export_format %||% "png"), "\n")
    export_multimapr_tree(
      tree            = tree,
      color_list      = edge_color_list,
      filename        = fn_export,
      type            = tree_type,
      format          = config$export_format %||% "png",
      lwd             = width1,
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
       legend_by_char  = if (config$show_legend %||% TRUE) legend_data$by_char,
       legend_corner   = config$legend_corner %||% "bottomleft",
       legend_layout   = config$legend_layout %||% "vertical",
       node_states     = isTRUE(config$show_nodes),
       node_shape      = config$node_shape %||% "circle",
       node_size       = config$node_size %||% 1
     )
   }

  # Screen render
  plot_multimapr_screen(
    tree            = tree,
    color_list      = edge_color_list,
    type            = tree_type,
    lwd             = width1,
    offset_range    = config$offset_range %||% 0.1,
    use_edge_length = isTRUE(config$use_edge_length),
    ladderize       = config$ladderize %||% FALSE,
      terminal_stretch = config$terminal_stretch %||% 1,
      tip_end          = config$tip_end %||% "round",
      branch_angle     = config$branch_angle %||% 45,
      show_labels      = config$show_labels %||% TRUE,
      label_size       = config$label_size %||% 1,
      label_color      = config$label_color %||% "black",
    legend_by_char  = if (config$show_legend %||% TRUE) legend_data$by_char,
    legend_corner   = config$legend_corner %||% "bottomleft",
    legend_layout   = config$legend_layout %||% "vertical",
    node_states     = isTRUE(config$show_nodes),
    node_shape      = config$node_shape %||% "circle",
    node_size       = config$node_size %||% 1
  )

  invisible(NULL)
}
