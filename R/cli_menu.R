################################################################################
# cli_menu.R
#
# Interactive interface module for MultiMapR.
# Contains all console menu logic: color prompts, character selection,
# state configuration, and full mapping setup.
#
# Dependencies (must exist in the global environment when this module is loaded,
# provided by data_utils.R):
#   align_tree_data() — aligns tips with character data
#   sort_states()     — sorts unique states of a character
#   is_valid_color()  — validates R color names / hex codes
#
# Author: MultiMapR — user interface module
################################################################################


# ==============================================================================
# SECTION 0 — ACCESSIBLE COLOR PALETTES
# ==============================================================================

#' Predefined color-blind-friendly palettes used in automatic color assignment.
#'
#' Used for Simple Mapping and single-character Ancestral Reconstruction.
#' Each palette contains up to 10 colors.
#'   okabe  — Okabe-Ito (most widely recommended for color-blind accessibility)
#'   tol    — Paul Tol Muted (high-contrast, print-friendly)
#'   plasma — Plasma perceptual gradient (viridis family, ordered data)
PALETAS_ACCESIBLES <- list(
  okabe  = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999", "#44AA99"),
  tol    = c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499", "#DDDDDD"),
  plasma = c("#0D0887", "#46039F", "#7201A8", "#9C179E", "#BD3786", "#D8576B", "#ED7953", "#FA9E3B", "#FDC926", "#F0F921")
)

#' Sequential color-blind-safe gradients for multi-character Ancestral Reconstruction.
#'
#' One gradient per character slot (up to 3). Each gradient has 10 steps
#' sampled from continuous perceptual scales:
#'   magma  — character 1 (dark purple → cream)
#'   mako   — character 2 (near-black → pale green)
#'   plasma — character 3 (dark blue → yellow)
GAMAS_MULTIMAPEO <- list(
  # 1. WARM RANGE (Yellows, Oranges, Reds, Pinks)
  # Ideal for Character 1.
  # Initial jump: Light Yellow -> Dark Red -> Medium Orange
  rojos = c(
    "#FF0000", # 0: Rojo puro brillante
    "#660000", # 1: Granate muy oscuro
    "#FF1493", # 2: Rosa neón (DeepPink)
    "#FF8C00", # 3: Naranja oscuro
    "#DC143C", # 4: Carmesí
    "#FFB6C1", # 5: Rosa claro
    "#8B0000", # 6: Rojo sangre
    "#C71585", # 7: Violeta medio
    "#FF4500", # 8: Rojo anaranjado
    "#330000"  # 9: Rojo casi negro
  ),

  # 2. GAMA VERDES (Para Carácter 2)
  # Salto: Verde lima -> Verde noche -> Cian puro
  verdes = c(
    "#00FF00", # 0: Verde lima puro
    "#004000", # 1: Verde bosque muy oscuro
    "#00FFFF", # 2: Cian puro (contrasta brutal con el lima)
    "#2E8B57", # 3: Verde mar
    "#ADFF2F", # 4: Amarillo verdoso
    "#00FF7F", # 5: Verde primavera
    "#808000", # 6: Oliva
    "#228B22", # 7: Verde hoja
    "#9ACD32", # 8: Verde amarillento
    "#002000"  # 9: Verde casi negro
  ),

  # 3. GAMA AZULES / MORADOS (Para Carácter 3)
  # Salto: Azul puro -> Índigo oscuro -> Azul cielo brillante
  azules = c(
    "#0000FF", # 0: Azul puro
    "#000066", # 1: Azul marino muy oscuro
    "#00BFFF", # 2: Azul cielo profundo
    "#8A2BE2", # 3: Violeta azulado brillante
    "#1E90FF", # 4: Azul Dodger
    "#4B0082", # 5: Índigo
    "#87CEFA", # 6: Azul claro
    "#4169E1", # 7: Azul real
    "#9370DB", # 8: Morado medio
    "#000033"  # 9: Azul casi negro
  )
)


# ==============================================================================
# SECTION 1 — EXIT UTILITY
# ==============================================================================

#' Checks whether the user wants to exit
#'
#' @param input String entered by the user.
#' @return TRUE if the user typed 'exit', FALSE otherwise.
check_exit <- function(input) {
  if (tolower(trimws(input)) == "exit") {
    cat("Exiting...\n")
    return(TRUE)
  }
  return(FALSE)
}


# ==============================================================================
# SECTION 2 — COLOR AND CHARACTER SELECTORS
# ==============================================================================

#' Prompts the user for a color for a given state
#'
#' Repeats the prompt until the user enters a valid color recognized
#' by R (name or hex) or types 'exit' to cancel.
#'
#' @param state String with the name of the state whose color is requested.
#' @return String with the chosen color, or invisible(NULL) if the user cancels.
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
#'
#' Iterates over each state in \code{sorted_states} calling
#' \code{prompt_state_color()} and builds a named vector state -> color.
#'
#' @param sorted_states  Vector of already-sorted states.
#' @param character_name String with the character name (for the header).
#' @return Named character vector state -> color.
prompt_character_colors <- function(sorted_states, character_name) {
  cat("\nAssign colors for character:", character_name, "\n")
  colors <- vapply(sorted_states, prompt_state_color, character(1))
  names(colors) <- sorted_states
  colors
}

#' Prompts the user to select one or several characters from a list
#'
#' Displays the available characters numbered and allows choosing between
#' \code{min_n} and \code{max_n} of them. If only one character is available
#' or \code{min_n == max_n == 1}, returns it directly without asking.
#'
#' @param available_chars Vector of available character names.
#' @param prompt_n        Prompt text asking for the number of characters.
#' @param min_n           Minimum number of characters to select (default 1).
#' @param max_n           Maximum number of characters to select (default all).
#' @return Character vector with selected character names,
#'         or invisible(NULL) if the user cancels.
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
#'
#' Displays the unique states of the character, allows choosing a subset
#' (or all with option 0), then requests a color for each selected state.
#' State selection is always interactive. If \code{auto_palette} is supplied,
#' colors are assigned automatically from that palette (no color prompts);
#' otherwise the user is asked to type a color for each state.
#'
#' @param aligned_data  Data.frame aligned with the phylogeny.
#' @param character     Name of the character column to configure.
#' @param auto_palette  Optional character vector of colors (a palette).
#'                      When not NULL, colors are auto-assigned but state
#'                      selection is still shown interactively.
#' @return Named character vector state -> color, or invisible(NULL) if cancelled.
prompt_states_and_colors <- function(aligned_data, character, auto_palette = NULL) {
  all_states <- sort_states(as.character(unique(aligned_data[[character]])))

  # --- State selection (always interactive) -----------------------------------
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
  # ----------------------------------------------------------------------------

  # --- Color assignment: automatic (palette) or manual -----------------------
  if (!is.null(auto_palette)) {
    if (length(selected_states) > length(auto_palette)) {
      cat(sprintf("  [!] Warning: '%s' has more selected states than colors in the palette. Colors will be recycled.\n",
                  character))
      colores_asignados <- rep(auto_palette, length.out = length(selected_states))
    } else {
      colores_asignados <- auto_palette[seq_along(selected_states)]
    }
    names(colores_asignados) <- selected_states
    cat(sprintf("  \u2192 Colors assigned automatically to '%s': %s\n",
                character, paste(selected_states, collapse = ", ")))
    return(colores_asignados)
  }
  # ----------------------------------------------------------------------------

  prompt_character_colors(selected_states, character)
}


# ==============================================================================
# SECTION 3 — FULL MAPPING CONFIGURATION
# ==============================================================================

#' Orchestrates all interactive menus and returns the configuration list
#'
#' Walks through menus for mapping type, characters, colors, algorithm,
#' line width, tree topology, branch lengths, and export options, then
#' packages all user decisions into a homogeneous \code{config} list consumed
#' by the rendering functions in \code{main_MultiMapR.R}.
#'
#' @param phylogeny        \code{phylo} object (ape).
#' @param character_data   Data.frame with a \code{"Species"} column and characters.
#' @param use_palettes     Logical. If \code{TRUE}, colors are assigned automatically
#'                         from \code{PALETAS_ACCESIBLES}; color prompts are skipped but
#'                         state-selection prompts are still shown interactively.
#' @return Configuration list, or invisible(NULL) if the user cancels.
#' @seealso \code{\link{execute_phylogeny}}
setup_mapping_config <- function(phylogeny, character_data, use_palettes = FALSE) {

  if (!inherits(phylogeny, "phylo"))
    stop("'phylogeny' must be a 'phylo' object.")
  if (!is.data.frame(character_data))
    stop("'character_data' must be a data.frame.")

  aligned_data   <- align_tree_data(phylogeny, character_data)
  available_chars <- colnames(aligned_data)[-1]

  cat("\n=== Welcome to MultiMapR ===\n")

  # --- Mapping type ------------------------------------------------------------
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
  config <- list(mapping_type = mapping_type, aligned_data = aligned_data)

  # Alias kept for backward compatibility with rendering functions
  config$datos_ord <- aligned_data

  # ==========================================================================
  # SIMPLE MAPPING
  # ==========================================================================
  if (mapping_type == 1) {

    selected_chars <- prompt_characters(
      available_chars,
      prompt_n = "How many characters to map?",
      min_n = 1, max_n = length(available_chars))
    if (is.null(selected_chars)) return(invisible(NULL))
    config$caracteres <- selected_chars

    colors_by_char <- list()
    for (i in seq_along(selected_chars)) {
      char         <- selected_chars[i]
      paleta_actual <- if (use_palettes && i <= length(PALETAS_ACCESIBLES))
        PALETAS_ACCESIBLES[[i]] else NULL
      cols <- prompt_states_and_colors(aligned_data, char, auto_palette = paleta_actual)
      if (is.null(cols)) return(invisible(NULL))
      colors_by_char[[char]] <- cols
    }
    config$colores_por_caracter <- colors_by_char

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
      config$pch_figura <- if (nchar(trimws(fig_str)) == 0 || is.na(pch_map[fig_str]))
        21L else pch_map[fig_str]

      size_str <- readline(prompt = "Figure size (positive number, Enter = 1): ")
      config$tam_figura <- if (nchar(trimws(size_str)) == 0) 1 else {
        t <- suppressWarnings(as.numeric(size_str))
        if (is.na(t) || t <= 0) { cat("Invalid value, using 1.\n"); 1 } else t
      }
    } else {
      config$pch_figura <- 21L   # valor por defecto, no se usará en el render
      config$tam_figura <- 1
    }

    # ==========================================================================
    # ANCESTRAL RECONSTRUCTION
    # ==========================================================================
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
    config$funcion_multi <- fun_sel

    selected_chars <- prompt_characters(
      available_chars,
      prompt_n = "How many characters to map?",
      min_n = 1L, max_n = min(3L, length(available_chars)))
    if (is.null(selected_chars)) return(invisible(NULL))
    config$caracteres <- selected_chars

    colors_by_char <- list()

    if (length(selected_chars) == 1) {
      # === 1 SOLO CARÁCTER: Usar Paleta Okabe-Ito (PALETAS_ACCESIBLES) ===
      char          <- selected_chars[1]
      paleta_actual <- if (use_palettes) PALETAS_ACCESIBLES[[1]] else NULL
      cols <- prompt_states_and_colors(aligned_data, char, auto_palette = paleta_actual)
      if (is.null(cols)) return(invisible(NULL))
      colors_by_char[[char]] <- cols

      all_states <- sort_states(as.character(unique(aligned_data[[char]])))
      config$mapear_todos <- length(cols) == length(all_states)

    } else {
      # === MULTIMAPEO (2 o 3 CARACTERES): Usar Gamas Secuenciales (GAMAS_MULTIMAPEO) ===
      # Carácter 1: Magma, Carácter 2: Mako, Carácter 3: Plasma
      for (i in seq_along(selected_chars)) {
        char          <- selected_chars[i]
        paleta_actual <- if (use_palettes && i <= length(GAMAS_MULTIMAPEO))
          GAMAS_MULTIMAPEO[[i]] else NULL
        cols <- prompt_states_and_colors(aligned_data, char, auto_palette = paleta_actual)
        if (is.null(cols)) return(invisible(NULL))
        colors_by_char[[char]] <- cols
      }
      # mapear_todos is TRUE only if ALL characters had ALL states selected;
      # otherwise some states will render in gray.
      config$mapear_todos <- all(vapply(selected_chars, function(char) {
        all_states <- sort_states(as.character(unique(aligned_data[[char]])))
        length(colors_by_char[[char]]) == length(all_states)
      }, logical(1)))
    }
    config$colores_por_caracter <- colors_by_char

    if (fun_sel == 2) {
      cat("\nFigure type at terminals:\n")
      cat("  1: Circle   (pch = 21)\n")
      cat("  2: Square   (pch = 22)\n")
      cat("  3: Triangle (pch = 24)\n")
      cat("  4: Diamond  (pch = 23)\n")
      fig_str <- readline(prompt = "Select (1\u20134, Enter = 1): ")
      if (check_exit(fig_str)) return(invisible(NULL))
      pch_map <- c(`1` = 21L, `2` = 22L, `3` = 24L, `4` = 23L)
      config$pch_figura <- if (nchar(trimws(fig_str)) == 0 || is.na(pch_map[fig_str]))
        21L else pch_map[fig_str]

      size_str <- readline(prompt = "Figure size (positive number, Enter = 1): ")
      if (check_exit(size_str)) return(invisible(NULL))
      config$tam_figura <- if (nchar(trimws(size_str)) == 0) 1 else {
        t <- suppressWarnings(as.numeric(size_str))
        if (is.na(t) || t <= 0) { cat("Invalid value, using 1.\n"); 1 } else t
      }
    }

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
    config$algoritmo <- algo_sel

    if (algo_sel == 2) {
      cat("\nFitch optimization mode:\n")
      cat("  1: ACCTRAN     (accelerated transformation \u2014 toward tips)\n")
      cat("  2: DELTRAN     (delayed transformation \u2014 toward root)\n")
      cat("  3: Unambiguous (only unambiguous states after both passes)\n")
      mode_str <- readline(prompt = "Select (1/2/3, Enter = 1): ")
      if (check_exit(mode_str)) return(invisible(NULL))
      fitch_mode <- if (nchar(trimws(mode_str)) == 0) 1L else as.integer(mode_str)
      if (is.na(fitch_mode) || !fitch_mode %in% 1:3) {
        cat("Invalid option, using ACCTRAN.\n")
        fitch_mode <- 1L
      }
      config$fitch_mode <- c("acctran", "deltran", "unambiguous")[fitch_mode]
    }
  }

  # --- Tree type ---------------------------------------------------------------
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
  config$tipo_arbol <- tree_type

  # --- Dynamic offset range based on branch width ------------------------------
  # branch_width is supplied as a parameter to execute_phylogeny() and
  # recalculated there once config$grosor is known.  Here we set a placeholder
  # using the default width (2) so config is always complete before rendering.
  # The offset between branches must grow proportionally to the width so
  # branches never overlap. The conversion factor varies by tree type because
  # data coordinates differ in scale:
  #   . phylogram / cladogram: Y axis is linear (1 unit ~= 1 tip),
  #     empirical width->offset relation is ~0.012 units/pt.
  #   . fan: coordinates are polar; spread is a percentage of the radius,
  #     so a larger factor is needed (~0.0025 x R_max per pt).
  #     Since R_max is unknown here, a relative scale factor is used,
  #     adjusted automatically in .emtree_render_fan().
  offset_factor <- switch(tree_type,
                          "phylogram"  = 0.012,
                          "cladogram"  = 0.012,
                          "fan"        = 0.003,
                          0.012
  )
  config$grosor        <- 2          # default; overridden by execute_phylogeny()
  config$rango_desfase <- max(0.05, config$grosor * offset_factor)

  # --- Branch lengths ----------------------------------------------------------
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

  # --- Export ------------------------------------------------------------------
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
                           paste(config$caracteres, collapse = "-"), "_",
                           config$tipo_arbol)
    fn_str <- readline(prompt = paste0(
      "Filename without extension (Enter = '", default_name, "'): "))
    if (check_exit(fn_str)) return(invisible(NULL))
    fn <- trimws(fn_str)
    fn <- sub("\\.(png|pdf)$", "", fn, ignore.case = TRUE)
    config$export_filename <- if (nzchar(fn)) fn else default_name

    # ── DIMENSIONS MENU ───────────────────────────────────────────────────────
    # The "Automatic" option preserves the original behavior (height ∝ Ntip).
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
    # ── END DIMENSIONS MENU ───────────────────────────────────────────────────

    # ── LEGEND CORNER MENU ────────────────────────────────────────────────────
    # Space for the legend is reserved by expanding the margin of the chosen
    # corner so it never overlaps labels or the plot itself.
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
    # ── END LEGEND CORNER MENU ────────────────────────────────────────────────

  } else {
    config$export_filename <- NULL
    config$export_format   <- NULL
    config$width           <- NULL
    config$height          <- NULL
    config$legend_corner   <- NULL
  }

  config
}
