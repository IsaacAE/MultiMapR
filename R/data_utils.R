# Side-effect-free helpers: color validation, state sorting, tree/data alignment,
# font size adjustment and legend data construction.


# Utility operator

#' Null-coalescing operator: returns x if not NULL, otherwise y
#' @noRd
`%||%` <- function(x, y) if (!is.null(x)) x else y


# Validation and transformation

#' Validates whether a string is a color recognised by R
#' @noRd
is_valid_color <- function(color) {
  color %in% colors() || grepl("^#(?:[0-9a-fA-F]{3}){1,2}$", color, perl = TRUE)
}

#' Sorts the unique states of a character (letters -> numbers -> other)
#' @noRd
sort_states <- function(states) {
  states_letters <- states[grep("^[A-Za-z]+$", states)]
  states_numbers <- states[grep("^[0-9]+$",    states)]
  states_other   <- setdiff(states, c(states_letters, states_numbers))
  sorted <- c(sort(states_letters),
              sort(as.numeric(states_numbers)),
              sort(states_other))
  sorted[!is.na(sorted)]
}

#' Aligns character data with the tip order of the phylogeny
#' @noRd
align_tree_data <- function(phylogeny, char_data) {
  order_tips <- phylogeny$tip.label
  aligned    <- char_data[match(order_tips, char_data$Species), ]
  if (any(is.na(aligned$Species))) {
    stop("Species names do not match between the tree and the character data.")
  }
  aligned
}


# Font size adjustment

#' Adjusts font size (cex) based on the number of tips of the phylogeny
#' @noRd
adjust_cex <- function(phylogeny, min_cex = 0.2, config = NULL) {
  n_tips <- Ntip(phylogeny)
  cex    <- max(1 / (1 + n_tips / 50), min_cex)

  # Scale text with a custom canvas height so it is not tiny in the exported file.
  if (!is.null(config) && !is.null(config$height)) {
    height_default <- n_tips * 0.25 + 2

    if (!is.null(config$tree_type) && config$tree_type == "fan") {
      width_default  <- 12
      height_default <- max(height_default, width_default)
    }
    scale_h <- config$height / height_default
    cex     <- cex * scale_h
  }

  return(cex)
}


# Legend data construction

#' Builds structured legend data grouped by character for export
#' @noRd
.build_legend_data <- function(colors_by_char, characters) {
  by_char    <- list()
  all_labels <- character(0)
  all_colors <- character(0)

  for (char in characters) {
    char_state_colors  <- colors_by_char[[char]]
    states <- names(char_state_colors)
    by_char[[char]] <- list(labels = states, colors = unname(char_state_colors))
    all_labels <- c(all_labels, states)
    all_colors <- c(all_colors, unname(char_state_colors))
  }

  list(by_char = by_char, labels = all_labels, colors = all_colors)
}
