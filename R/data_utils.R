################################################################################
# data_utils.R
#
# Utilidades puras de datos para MultiMapR.
# Contiene funciones sin efectos secundarios para validación de colores,
# ordenación de estados, alineación de datos filogenéticos, ajuste de
# tamaño de fuente y construcción de estructuras de leyenda.
#
# Sin dependencias de otros módulos MultiMapR.
# Requiere: ape (para Ntip en adjust_cex).
#
# Autor: MultiMapR — módulo de utilidades de datos
################################################################################


# ==============================================================================
# SECCIÓN 1 — OPERADOR UTILITARIO
# ==============================================================================

#' Operador null-coalescing: devuelve x si no es NULL, y en caso contrario
`%||%` <- function(x, y) if (!is.null(x)) x else y


# ==============================================================================
# SECCIÓN 2 — VALIDACIÓN Y TRANSFORMACIÓN
# ==============================================================================

#' Validates whether a string is a color recognised by R
#'
#' Accepts both R colour names (from \code{colors()}) and CSS-style hex strings
#' of 3 or 6 hex digits (e.g. \code{"#A3F"}, \code{"#1a2b3c"}).
#'
#' @param color String to validate.
#' @return \code{TRUE} if valid, \code{FALSE} otherwise.
is_valid_color <- function(color) {
  color %in% colors() || grepl("^#(?:[0-9a-fA-F]{3}){1,2}$", color, perl = TRUE)
}

#' Sorts the unique states of a character (letters → numbers → other)
#'
#' Partitions \code{states} into three buckets — pure-letter tokens, pure-digit
#' tokens, and everything else — sorts each bucket independently, and
#' concatenates them.  NAs in the result are dropped.
#'
#' @param states Character vector of unique states.
#' @return Sorted character vector without NAs.
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
#'
#' Reorders \code{char_data} so that its rows match \code{phylogeny$tip.label}
#' exactly.  Stops with an informative message if any tip name is missing from
#' the data frame.
#'
#' @param phylogeny  A \code{phylo} object (ape).
#' @param char_data  Data frame with a \code{"Species"} column.
#' @return Reordered data frame aligned to the tip order of \code{phylogeny}.
align_tree_data <- function(phylogeny, char_data) {
  order_tips <- phylogeny$tip.label
  aligned    <- char_data[match(order_tips, char_data$Species), ]
  if (any(is.na(aligned$Species))) {
    stop("Species names do not match between the tree and the character data.")
  }
  aligned
}


# ==============================================================================
# SECCIÓN 3 — AJUSTE DE TAMAÑO DE FUENTE
# ==============================================================================

#' Adjusts font size (cex) based on the number of tips of the phylogeny
#'
#' The base formula scales down cex as tip count grows, with a configurable
#' minimum floor.  When \code{config$height} is present, an additional vertical
#' scaling factor is applied so that text remains proportional to the exported
#' canvas size.
#'
#' @param phylogeny A \code{phylo} object (ape).
#' @param min_cex   Minimum allowed cex value (default \code{0.2}).
#' @param config    Optional configuration list.  If it contains \code{$height},
#'                  the result is multiplied by the ratio
#'                  \code{height / height_default} so that text scales with a
#'                  custom canvas height.
#' @return A single numeric cex value.
adjust_cex <- function(phylogeny, min_cex = 0.2, config = NULL) {
  n_tips <- Ntip(phylogeny)
  cex    <- max(1 / (1 + n_tips / 50), min_cex)
  
  # ── Scale by custom canvas height ──────────────────────────────────────────
  # If the user set a canvas height different from the automatic default, text
  # must grow proportionally so it is not tiny in the exported file.
  if (!is.null(config) && !is.null(config$height)) {
    height_default <- n_tips * 0.25 + 2
    
    if (!is.null(config$tipo_arbol) && config$tipo_arbol == "fan") {
      width_default  <- 12
      height_default <- max(height_default, width_default)
    }
    scale_h <- config$height / height_default
    cex     <- cex * scale_h
  }
  
  return(cex)
}


# ==============================================================================
# SECCIÓN 4 — CONSTRUCCIÓN DE DATOS DE LEYENDA
# ==============================================================================

#' Builds structured legend data grouped by character for export
#'
#' Returns two parallel representations of the same colour mapping:
#' \describe{
#'   \item{\code{by_char}}{Named list \emph{character} → \code{list(labels, colors)},
#'     consumed by \code{export_multimapr_tree(legend_by_char = ...)} to draw
#'     one labelled block per character.}
#'   \item{\code{labels} / \code{colors}}{Flat parallel vectors for legacy
#'     \code{legend()} calls.}
#' }
#'
#' @param colores_por_car Named list \emph{character} → named vector
#'   \emph{state → color}.
#' @param caracteres Character vector of character names to include.
#' @return List with elements \code{by_char}, \code{labels}, and \code{colors}.
.build_legend_data <- function(colores_por_car, caracteres) {
  by_char    <- list()
  all_labels <- character(0)
  all_colors <- character(0)
  
  for (car in caracteres) {
    col_e   <- colores_por_car[[car]]
    estados <- names(col_e)
    by_char[[car]] <- list(labels = estados, colors = unname(col_e))
    all_labels <- c(all_labels, estados)
    all_colors <- c(all_colors, unname(col_e))
  }
  
  list(by_char = by_char, labels = all_labels, colors = all_colors)
}