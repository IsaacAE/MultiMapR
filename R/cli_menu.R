# Interactive console menus: color prompts, character selection, state configuration
# and the full mapping setup (setup_mapping_config()).


# Accessible color palettes

#' Predefined color-blind-friendly palettes used in automatic color assignment
#' @noRd
ACCESSIBLE_PALETTES <- list(
  okabe  = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999", "#44AA99"),
  tol    = c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499", "#DDDDDD"),
  plasma = c("#0D0887", "#46039F", "#7201A8", "#9C179E", "#BD3786", "#D8576B", "#ED7953", "#FA9E3B", "#FDC926", "#F0F921")
)

#' Sequential color-blind-safe gradients for multi-character Ancestral Reconstruction
#' @noRd
MULTIMAP_RAMPS <- list(
  # Warm range (yellows, oranges, reds, pinks): character 1
  reds = c(
    "#FF0000",
    "#660000",
    "#FF1493",
    "#FF8C00",
    "#DC143C",
    "#FFB6C1",
    "#8B0000",
    "#C71585",
    "#FF4500",
    "#330000"
  ),

  # Greens: character 2
  greens = c(
    "#00FF00",
    "#004000",
    "#00FFFF",
    "#2E8B57",
    "#ADFF2F",
    "#00FF7F",
    "#808000",
    "#228B22",
    "#9ACD32",
    "#002000"
  ),

  # Blues / purples: character 3
  blues = c(
    "#0000FF",
    "#000066",
    "#00BFFF",
    "#8A2BE2",
    "#1E90FF",
    "#4B0082",
    "#87CEFA",
    "#4169E1",
    "#9370DB",
    "#000033"
  ),

  # Golds / browns: character 4
  golds = c(
    "#FFD700",
    "#4A2C00",
    "#FFFF66",
    "#B8860B",
    "#F5DEB3",
    "#8B4513",
    "#FFC000",
    "#CD853F",
    "#FFF8B0",
    "#2B1A00"
  ),

  # Purples / magentas: character 5
  purples = c(
    "#9400D3",
    "#2A0033",
    "#FF00FF",
    "#6A0DAD",
    "#DDA0DD",
    "#800080",
    "#E6C8FF",
    "#BA55D3",
    "#4B0055",
    "#C71FCB"
  )
)

#' Maximum number of characters superimposed in one multimapping (ancestral reconstruction)
#' @noRd
MAX_MULTIMAP <- 5L

# Named predefined palettes

PRESET_PALETTES <- list(
  wc  = c("green", "blue", "black", "red", "yellow", "cyan1", "darkgreen"),
  MM1 = c("darkviolet", "orange", "brown", "gray50", "firebrick1", "forestgreen")
)


# Exit utility

#' Checks whether the user wants to exit
#' @noRd
check_exit <- function(input) {
  if (tolower(trimws(input)) == "exit") {
    cat("Exiting...\n")
    return(TRUE)
  }
  return(FALSE)
}


# Color and character selectors

#' Prompts the user for a color for a given state
#' @noRd
prompt_state_color <- function(state) {
  prompt_msg <- paste0("Color for state '", state,
                       "' (name or hex, 'exit' to quit): ")
  color <- readline(prompt = prompt_msg)
  if (check_exit(color)) return(invisible(NULL))

  while (!is_valid_color(color)) {
    cat("Invalid color. Please try again.\n")
    color <- readline(prompt = prompt_msg)
    if (check_exit(color)) return(invisible(NULL))
  }
  color
}

#' Prompts the user for colors for all states of a character
#' @noRd
prompt_character_colors <- function(sorted_states, character_name) {
  cat("\nAssign colors for character:", character_name, "\n")
  colors <- vapply(sorted_states, prompt_state_color, character(1))
  names(colors) <- sorted_states
  colors
}

#' Prompts the user to select one or several characters from a list
#' @noRd
prompt_characters <- function(available_chars, prompt_n, min_n = 1,
                              max_n = length(available_chars)) {
  cat("\nAvailable characters:\n")
  for (i in seq_along(available_chars))
    cat(sprintf("  %d: %s\n", i, available_chars[i]))

  if (min_n == max_n && max_n == 1) {
    if (length(available_chars) == 1) {
      cat("Only one character available:", available_chars[1], "\n")
      return(available_chars[1])
    }
    repeat {
      sel <- readline(prompt = prompt_n)
      if (check_exit(sel)) return(invisible(NULL))
      sel <- suppressWarnings(as.integer(sel))
      if (!is.na(sel) && sel >= 1 && sel <= length(available_chars)) break
      cat(sprintf("Invalid selection. Please enter a number between 1 and %d.\n",
                  length(available_chars)))
    }
    return(available_chars[sel])
  }

  range_txt <- if (min_n == max_n) as.character(min_n)
  else paste0(min_n, "\u2013", max_n)
  repeat {
    n_str <- readline(prompt = paste0(prompt_n, " (", range_txt, " or 'exit' to quit): "))
    if (check_exit(n_str)) return(invisible(NULL))
    n <- suppressWarnings(as.integer(n_str))
    if (!is.na(n) && n >= min_n && n <= max_n) break
    cat(sprintf("Invalid number. Must be between %d and %d.\n", min_n, max_n))
  }

  selected <- character(n)
  for (i in seq_len(n)) {
    repeat {
      s <- readline(prompt = paste0("  Character ", i, " (number, or 'exit' to quit): "))
      if (check_exit(s)) return(invisible(NULL))
      s <- suppressWarnings(as.integer(s))
      if (!is.na(s) && s >= 1 && s <= length(available_chars) &&
          !available_chars[s] %in% selected) break
      if (!is.na(s) && available_chars[s] %in% selected)
        cat("Character already selected. Please choose a different one.\n")
      else
        cat(sprintf("Invalid selection. Please enter a number between 1 and %d.\n",
                    length(available_chars)))
    }
    selected[i] <- available_chars[s]
  }
  selected
}

#' Prompts the user which states of a character to include and their colors
#' @noRd
prompt_states_and_colors <- function(aligned_data, character, auto_palette = NULL, strict = FALSE) {
  all_states <- sort_states(as.character(unique(aligned_data[[character]])))

  # State selection (always interactive)
  cat("\nAvailable states for '", character, "':\n", sep = "")
  for (i in seq_along(all_states)) cat(sprintf("  %d: %s\n", i, all_states[i]))

  cat("  0: All states\n")
  repeat {
    sel_str <- readline(prompt = "Which states to include? (0 = all, or numbers separated by commas): ")
    if (check_exit(sel_str)) return(invisible(NULL))
    if (trimws(sel_str) == "0") {
      selected_states <- all_states
      break
    }
    indices <- suppressWarnings(as.integer(unlist(strsplit(sel_str, ","))))
    if (!any(is.na(indices)) && all(indices >= 1) && all(indices <= length(all_states))) {
      selected_states <- all_states[indices]
      break
    }
    cat(sprintf("Invalid selection. Enter 0 or numbers between 1 and %d separated by commas.\n",
                length(all_states)))
  }

  # Color assignment: automatic (palette) or manual
  if (!is.null(auto_palette)) {
    if (length(selected_states) > length(auto_palette)) {
      if (strict) {
        stop(sprintf("Insufficient colors in palette: character '%s' has %d selected state(s) but only %d color(s) provided. Add more colors to the palette.",
                     character, length(selected_states), length(auto_palette)))
      }
      cat(sprintf("  [!] Warning: '%s' has more selected states than colors in the palette. Colors will be recycled.\n",
                  character))
      assigned_colors <- rep(auto_palette, length.out = length(selected_states))
    } else {
      assigned_colors <- auto_palette[seq_along(selected_states)]
    }
    names(assigned_colors) <- selected_states
    cat(sprintf("  \u2192 Colors assigned automatically to '%s': %s\n",
                character, paste(selected_states, collapse = ", ")))
    return(assigned_colors)
  }

  prompt_character_colors(selected_states, character)
}


#' Prompts for the ancestral reconstruction algorithm (and Fitch mode)
#' @noRd
prompt_ancestral_algorithm <- function() {
  cat("\nAncestral reconstruction algorithm:\n")
  cat("  1: Default (depth-weighted majority)\n")
  cat("  2: Fitch (ACCTRAN/DELTRAN/Unambiguous \u2014 loaded automatically from fitch.R)\n")
  repeat {
    algo_sel <- readline(prompt = "Select (1/2, or 'exit' to quit): ")
    if (check_exit(algo_sel)) return(invisible(NULL))
    algo_sel <- suppressWarnings(as.integer(algo_sel))
    if (!is.na(algo_sel) && algo_sel %in% 1:2) break
    cat("Invalid option. Please enter 1 or 2.\n")
  }

  fitch_mode <- NULL
  if (algo_sel == 2) {
    cat("\nFitch optimization mode:\n")
    cat("  1: ACCTRAN     (accelerated transformation \u2014 changes placed near the root)\n")
    cat("  2: DELTRAN     (delayed transformation \u2014 changes delayed toward the tips)\n")
    cat("  3: Unambiguous (only nodes with a single state in their MPR set)\n")
    mode_str <- readline(prompt = "Select (1/2/3, Enter = 1): ")
    if (check_exit(mode_str)) return(invisible(NULL))
    fitch_sel <- if (nchar(trimws(mode_str)) == 0) 1L else as.integer(mode_str)
    if (is.na(fitch_sel) || !fitch_sel %in% 1:3) {
      cat("Invalid option, using ACCTRAN.\n")
      fitch_sel <- 1L
    }
    fitch_mode <- c("acctran", "deltran", "unambiguous")[fitch_sel]
  }

  list(algorithm = c("majority", "fitch")[algo_sel], fitch_mode = fitch_mode)
}


# Post-run loop menu

#' Prompts the user for what to do after a plot has been rendered
#' @noRd
prompt_post_run_menu <- function(allow_algorithm_change = FALSE) {
  cat("\n=== What would you like to do next? ===\n")
  cat("  1: Run again with a completely new configuration\n")
  if (allow_algorithm_change) {
    cat("  2: Keep the current configuration and choose a different algorithm\n")
    cat("  3: Exit\n")
    valid <- 1:3
  } else {
    cat("  3: Exit\n")
    valid <- c(1L, 3L)
  }

  repeat {
    opt <- readline(prompt = "Select an option: ")
    if (check_exit(opt)) return("exit")
    opt <- suppressWarnings(as.integer(opt))
    if (!is.na(opt) && opt %in% valid) break
    cat("Invalid option. Please try again.\n")
  }

  switch(as.character(opt),
         "1" = "new",
         "2" = "algorithm",
         "3" = "exit")
}


# Full mapping configuration

#' Orchestrates all interactive menus and returns the configuration list
#' @noRd
setup_mapping_config <- function(phylogeny, character_data, use_palettes = FALSE,
                                 user_palette = NULL, ambiguity_color = NULL) {

  if (!inherits(phylogeny, "phylo"))
    stop("'phylogeny' must be a 'phylo' object.")
  if (!is.data.frame(character_data))
    stop("'character_data' must be a data.frame.")

  aligned_data   <- align_tree_data(phylogeny, character_data)
  available_chars <- colnames(aligned_data)[-1]

  cat("\n=== Welcome to MultiMapR ===\n")

  print_character_stats(aligned_data, available_chars)

  # Mapping type
  cat("\nMapping type:\n")
  cat("  1: Simple mapping \u2014 colored figures at terminals\n")
  cat("  2: Ancestral reconstruction \u2014 branch coloring\n")
  repeat {
    mapping_type <- readline(prompt = "Select (1/2, or 'exit' to quit): ")
    if (check_exit(mapping_type)) return(invisible(NULL))
    mapping_type <- suppressWarnings(as.integer(mapping_type))
    if (!is.na(mapping_type) && mapping_type %in% 1:2) break
    cat("Invalid option. Please enter 1 or 2.\n")
  }
  config <- list(mapping_type = mapping_type, aligned_data = aligned_data,
                 ambiguity_color = ambiguity_color)

  # SIMPLE MAPPING
  if (mapping_type == 1) {

    selected_chars <- prompt_characters(
      available_chars,
      prompt_n = "How many characters to map?",
      min_n = 1, max_n = length(available_chars))
    if (is.null(selected_chars)) return(invisible(NULL))
    config$characters <- selected_chars

    colors_by_char <- list()
    for (i in seq_along(selected_chars)) {
      char         <- selected_chars[i]
      if (!is.null(user_palette)) {
        current_palette <- user_palette
      } else {
        current_palette <- if (use_palettes && i <= length(ACCESSIBLE_PALETTES))
          ACCESSIBLE_PALETTES[[i]] else NULL
      }
      cols <- prompt_states_and_colors(aligned_data, char,
                                       auto_palette = current_palette,
                                       strict = !is.null(user_palette))
      if (is.null(cols)) return(invisible(NULL))
      colors_by_char[[char]] <- cols
    }
    config$colors_by_character <- colors_by_char

    cat("\nSimple mapping display mode:\n")
    cat("  1: Colored figures at terminals (circles, squares, etc.)\n")
    cat("  2: Colored tip labels only (no figures)\n")
    mode_sm_str <- readline(prompt = "Select (1/2, Enter = 1): ")
    if (check_exit(mode_sm_str)) return(invisible(NULL))
    config$simple_mode <- if (trimws(mode_sm_str) == "2") "tip_color" else "figures"

    if (config$simple_mode == "figures") {
      cat("\nFigure type at terminals:\n")
      cat("  1: Circle   (pch = 21)\n")
      cat("  2: Square   (pch = 22)\n")
      cat("  3: Triangle (pch = 24)\n")
      cat("  4: Diamond  (pch = 23)\n")
      fig_str <- readline(prompt = "Select (1\u20134, Enter = 1): ")
      pch_map <- c(`1` = 21L, `2` = 22L, `3` = 24L, `4` = 23L)
      config$figure_pch <- if (nchar(trimws(fig_str)) == 0 || is.na(pch_map[fig_str]))
        21L else pch_map[fig_str]

      size_str <- readline(prompt = "Figure size (positive number, Enter = 1): ")
      config$figure_size <- if (nchar(trimws(size_str)) == 0) 1 else {
        t <- suppressWarnings(as.numeric(size_str))
        if (is.na(t) || t <= 0) { cat("Invalid value, using 1.\n"); 1 } else t
      }
    } else {
      config$figure_pch <- 21L   # unused in this mode
      config$figure_size <- 1
    }

    # ANCESTRAL RECONSTRUCTION
  } else {

    cat("\nVisualization mode:\n")
    cat("  1: Branch superimposition (one or several characters)\n")
    cat("  2: Terminal figures + colored tree\n")
    repeat {
      fun_sel <- readline(prompt = "Select (1/2, or 'exit' to quit): ")
      if (check_exit(fun_sel)) return(invisible(NULL))
      fun_sel <- suppressWarnings(as.integer(fun_sel))
      if (!is.na(fun_sel) && fun_sel %in% 1:2) break
      cat("Invalid option. Please enter 1 or 2.\n")
    }
    config$multi_function <- fun_sel

    selected_chars <- prompt_characters(
      available_chars,
      prompt_n = "How many characters to map?",
      min_n = 1L, max_n = min(MAX_MULTIMAP, length(available_chars)))
    if (is.null(selected_chars)) return(invisible(NULL))
    config$characters <- selected_chars

    colors_by_char <- list()
    manual_chars   <- character(0)  # chars whose colors were TYPED, not auto-assigned

    if (length(selected_chars) == 1) {
      # One character: Okabe-Ito palette (ACCESSIBLE_PALETTES)
      char          <- selected_chars[1]
      if (!is.null(user_palette)) {
        current_palette <- user_palette
      } else {
        current_palette <- if (use_palettes) ACCESSIBLE_PALETTES[[1]] else NULL
      }
      if (is.null(current_palette)) manual_chars <- c(manual_chars, char)
      cols <- prompt_states_and_colors(aligned_data, char,
                                       auto_palette = current_palette,
                                       strict = !is.null(user_palette))
      if (is.null(cols)) return(invisible(NULL))
      colors_by_char[[char]] <- cols

      all_states <- sort_states(as.character(unique(aligned_data[[char]])))
      config$map_all <- length(cols) == length(all_states)

    } else {
      # 2-5 characters: sequential ramps (MULTIMAP_RAMPS): reds, greens, blues, golds, purples
      for (i in seq_along(selected_chars)) {
        char          <- selected_chars[i]
        if (!is.null(user_palette)) {
          current_palette <- user_palette
        } else {
          current_palette <- if (use_palettes && i <= length(MULTIMAP_RAMPS))
            MULTIMAP_RAMPS[[i]] else NULL
        }
        if (is.null(current_palette)) manual_chars <- c(manual_chars, char)
        cols <- prompt_states_and_colors(aligned_data, char,
                                         auto_palette = current_palette,
                                         strict = !is.null(user_palette))
        if (is.null(cols)) return(invisible(NULL))
        colors_by_char[[char]] <- cols
      }
      # map_all is TRUE only if ALL characters had ALL states selected;
      # otherwise some states will render in gray.
      config$map_all <- all(vapply(selected_chars, function(char) {
        all_states <- sort_states(as.character(unique(aligned_data[[char]])))
        length(colors_by_char[[char]]) == length(all_states)
      }, logical(1)))
    }
    config$colors_by_character <- colors_by_char

    # If the user typed a color for "?" or "-" while picking state colors and no global
    # ambiguity_color was given, promote it so Fitch renders it instead of magenta.
    # Auto-assigned palette colors are not eligible: landing on "?" there is incidental.
    if (is.null(config$ambiguity_color)) {
      for (char in manual_chars) {
        cols <- colors_by_char[[char]]
        if ("?" %in% names(cols)) {
          config$ambiguity_color <- unname(cols[["?"]])
          break
        } else if ("-" %in% names(cols)) {
          config$ambiguity_color <- unname(cols[["-"]])
          break
        }
      }
    }

    if (fun_sel == 2) {
      cat("\nFigure type at terminals:\n")
      cat("  1: Circle   (pch = 21)\n")
      cat("  2: Square   (pch = 22)\n")
      cat("  3: Triangle (pch = 24)\n")
      cat("  4: Diamond  (pch = 23)\n")
      fig_str <- readline(prompt = "Select (1\u20134, Enter = 1): ")
      if (check_exit(fig_str)) return(invisible(NULL))
      pch_map <- c(`1` = 21L, `2` = 22L, `3` = 24L, `4` = 23L)
      config$figure_pch <- if (nchar(trimws(fig_str)) == 0 || is.na(pch_map[fig_str]))
        21L else pch_map[fig_str]

      size_str <- readline(prompt = "Figure size (positive number, Enter = 1): ")
      if (check_exit(size_str)) return(invisible(NULL))
      config$figure_size <- if (nchar(trimws(size_str)) == 0) 1 else {
        t <- suppressWarnings(as.numeric(size_str))
        if (is.na(t) || t <= 0) { cat("Invalid value, using 1.\n"); 1 } else t
      }
    }

    algo_result <- prompt_ancestral_algorithm()
    if (is.null(algo_result)) return(invisible(NULL))
    config$algorithm  <- algo_result$algorithm
    config$fitch_mode <- algo_result$fitch_mode
  }

  # Tree type
  letter_to_type <- c("p" = "phylogram", "c" = "cladogram", "f" = "fan")
  cat("\nTree type:\n")
  cat("  p: Phylogram\n")
  cat("  c: Cladogram\n")
  cat("  f: Fan\n")
  repeat {
    tree_input <- tolower(trimws(readline(prompt = "Select (p/c/f, or 'exit' to quit): ")))
    if (check_exit(tree_input)) return(invisible(NULL))
    if (tree_input %in% names(letter_to_type)) break
    cat("Invalid option. Please enter p, c, or f.\n")
  }
  tree_type <- letter_to_type[[tree_input]]
  cat(sprintf("  \u2192 Tree type: %s\n", tree_type))
  config$tree_type <- tree_type

  # Offset between branches grows with branch width so they never overlap: ~0.012
  # units/pt for phylogram/cladogram (1 unit ~ 1 tip); for fan it is relative to R_max
  # (~0.0025 * R_max per pt), adjusted in .emtree_render_fan(). branch_width is a
  # placeholder here (default 2); execute_phylogeny() overrides it.
  offset_factor <- switch(tree_type,
                          "phylogram"  = 0.012,
                          "cladogram"  = 0.012,
                          "fan"        = 0.003,
                          0.012
  )
  config$branch_width        <- 2          # default; overridden by execute_phylogeny()
  config$offset_range <- max(0.05, config$branch_width * offset_factor)

  # Branch lengths
  if (!is.null(phylogeny$edge.length) && tree_type %in% c("phylogram", "fan")) {
    cat("\nThe tree has branch lengths.\n")
    cat("  1: Use proportional branch lengths in the plot\n")
    cat("  2: Ignore lengths (uniform branches, cladogram style)\n")
    len_str <- readline(prompt = "Select (1/2, Enter = 1): ")
    if (check_exit(len_str)) return(invisible(NULL))
    config$use_edge_length <- if (trimws(len_str) == "2") FALSE else TRUE
  } else {
    config$use_edge_length <- FALSE
  }

  # Export
  cat("\nExport the plot?\n")
  cat("  1: Yes \u2014 save to file\n")
  cat("  2: No \u2014 display in R window\n")
  exp_str <- readline(prompt = "Select (1/2, Enter = 2): ")
  if (check_exit(exp_str)) return(invisible(NULL))
  do_export <- trimws(exp_str) == "1"

  if (do_export) {

    # Output format
    cat("\nOutput format:\n")
    cat("  1: PNG (recommended for screen)\n")
    cat("  2: PDF (vector, ideal for publications)\n")
    fmt_str <- readline(prompt = "Select (1/2, Enter = 1): ")
    if (check_exit(fmt_str)) return(invisible(NULL))
    config$export_format <- if (trimws(fmt_str) == "2") "pdf" else "png"

    # Filename without extension
    default_name <- paste0("MultiMapR_",
                           paste(config$characters, collapse = "-"), "_",
                           config$tree_type)
    filename_input <- readline(prompt = paste0(
      "Filename without extension (Enter = '", default_name, "'): "))
    if (check_exit(filename_input)) return(invisible(NULL))
    filename <- trimws(filename_input)
    filename <- sub("\\.(png|pdf)$", "", filename, ignore.case = TRUE)
    config$export_filename <- if (nzchar(filename)) filename else default_name

    # Dimensions
    cat("\nExported file dimensions:\n")
    cat("  1: Automatic (recommended \u2014 height proportional to number of terminals)\n")
    cat("  2: Custom (width and height in inches)\n")
    dim_str <- readline(prompt = "Select (1/2, Enter = 1): ")
    if (check_exit(dim_str)) return(invisible(NULL))

    if (trimws(dim_str) == "2") {
      w_str <- readline(prompt = "Width in inches (positive number, e.g. 14): ")
      if (check_exit(w_str)) return(invisible(NULL))
      w_val <- suppressWarnings(as.numeric(trimws(w_str)))
      if (is.na(w_val) || w_val <= 0) {
        cat("Invalid width. Automatic dimensions will be used.\n")
        config$width  <- NULL
        config$height <- NULL
      } else {
        h_str <- readline(prompt = "Height in inches (positive number, e.g. 20): ")
        if (check_exit(h_str)) return(invisible(NULL))
        h_val <- suppressWarnings(as.numeric(trimws(h_str)))
        if (is.na(h_val) || h_val <= 0) {
          cat("Invalid height. Automatic dimensions will be used.\n")
          config$width  <- NULL
          config$height <- NULL
        } else {
          config$width  <- w_val
          config$height <- h_val
          cat(sprintf("  \u2192 Exporting at %.1f \u00d7 %.1f inches.\n", w_val, h_val))
        }
      }
    } else {
      config$width  <- NULL
      config$height <- NULL
    }

    # Legend corner (that margin is expanded so the legend never overlaps the plot)
    cat("\nIn which corner do you want to place the legend?\n")
    cat("  1: Bottom-left   (bottomleft)  \u2014 recommended for phylogram\n")
    cat("  2: Bottom-right  (bottomright)\n")
    cat("  3: Top-left      (topleft)\n")
    cat("  4: Top-right     (topright)    \u2014 useful for fan with dense branches\n")
    corner_map <- c("1" = "bottomleft", "2" = "bottomright",
                    "3" = "topleft",    "4" = "topright")
    corner_str <- readline(prompt = "Select (1\u20134, Enter = 1): ")
    if (check_exit(corner_str)) return(invisible(NULL))
    config$legend_corner <- if (nzchar(trimws(corner_str)) && trimws(corner_str) %in% names(corner_map))
      corner_map[[trimws(corner_str)]] else "bottomleft"
    cat(sprintf("  \u2192 Legend at: %s\n", config$legend_corner))

  } else {
    config$export_filename <- NULL
    config$export_format   <- NULL
    config$width           <- NULL
    config$height          <- NULL
    config$legend_corner   <- NULL
  }

  config
}

