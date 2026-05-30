#' Main orchestrator of MultiMapR
#'
#' @param phylogeny         Path to the tree file or a phylo object.
#' @param character_data    Path to the CSV file or a data.frame.
#' @param branch_width      Optional branch width.
#' @param label_offset      Dynamic label offset.
#' @import ape
#' @export
.export_device <- function(filename, format, phylogeny, tree_type, expr,
                           width = NULL, height = NULL) {
  n_tips         <- Ntip(phylogeny)
  height_default <- n_tips * 0.25 + 2
  width_default  <- 12
  if (tree_type == "fan") {
    height_default <- max(height_default, width_default)
    width_default  <- height_default
  }
  height_in <- if (!is.null(height)) height else height_default
  width_in  <- if (!is.null(width))  width  else width_default

  dir.create("Exports", showWarnings = FALSE, recursive = TRUE)
  path <- file.path("Exports", paste0(basename(filename), ".", format))
  if (format == "pdf") {
    pdf(file = path, width = width_in, height = height_in)
  } else {
    png(filename = path,
        width    = width_in,
        height   = height_in,
        units    = "in",
        res      = 300,
        bg       = "white")
  }
  tryCatch(
    force(expr),
    error   = function(e) stop(sprintf("[MultiMapR] Error exporting %s: %s",
                                       toupper(format), conditionMessage(e))),
    finally = dev.off()
  )
  dims_origin <- if (!is.null(width) || !is.null(height)) " [custom]" else " [auto]"
  message(sprintf("[MultiMapR] File exported: %s  (%.1f \u00d7 %.1f in%s%s)",
                  path, width_in, height_in,
                  if (format == "png") " @ 300 dpi" else "",
                  dims_origin))
  invisible(path)
}


# ==============================================================================
# MAIN FUNCTION — ORCHESTRATOR
# ==============================================================================

#' Main orchestrator of MultiMapR
#'
#' Supports argument polymorphism: if `phylogeny` and `character_data`
#' are strings (file paths on disk), the data are loaded automatically
#' through `load_data()` before configuring the mapping.
#'
#' @param phylogeny          phylo object — OR — string path to the tree file.
#' @param character_data     Data.frame with a "Species" column — OR — string path
#'                           to the character CSV file.
#' @param branch_width       Branch width (default NULL = chosen by the user in
#'                           the menu). A positive number overrides the menu.
#' @param label_offset       Distance between tips and labels (default 0.3).
#' @param sep                CSV field separator (default ",").
#'                           Only used when file paths are passed.
#' @param species_col        Species column in the CSV: integer index or column
#'                           name (default 1).
#'                           Only used when file paths are passed.
#' @param normalize_spaces   If TRUE converts spaces to "_" in species names in
#'                           the CSV (useful when the tree uses "_" and the CSV
#'                           uses spaces). Default FALSE.
#'                           Only used when file paths are passed.
#' @param tree_format        Tree file format: "auto" (default), "newick" or "nexus".
#'                           Only used when file paths are passed.
#' @param strict             If TRUE raises an error on tree/CSV discrepancies;
#'                           if FALSE (default) only warns.
#'                           Only used when file paths are passed.
#' @return Invisible NULL.
#' @export
execute_phylogeny <- function(phylogeny, character_data,
                              branch_width         = NULL,
                              label_offset         = 0.3,
                              sep                  = ",",
                              species_col          = 1,
                              normalize_spaces     = FALSE,
                              tree_format          = "auto",
                              strict               = FALSE) {
  tryCatch({

    # POLYMORPHISM: If file paths (character) are passed, load data automatically
    if (is.character(phylogeny) && is.character(character_data)) {
      cat("\n[MultiMapR] File paths detected. Loading data automatically...\n")
      loaded_data      <- load_data(tree_path        = phylogeny,
                                    csv_path         = character_data,
                                    sep              = sep,
                                    species_col      = species_col,
                                    normalize_spaces = normalize_spaces,
                                    tree_format      = tree_format,
                                    strict           = strict)
      phylogeny        <- loaded_data$tree
      character_data   <- loaded_data$characters
    }

    # Continue with the normal system flow using the modular interface
    config <- setup_mapping_config(phylogeny, character_data)
    if (is.null(config)) return(invisible(NULL))

    # Force uniform branches if the user requested it in the menu
    if (isFALSE(config$use_edge_length)) {
      phylogeny$edge.length <- NULL
    }

    # External branch width configuration if passed as argument
    if (!is.null(branch_width) && is.numeric(branch_width) && branch_width > 0) {
      config$grosor <- branch_width
      offset_factor <- switch(config$tipo_arbol %||% "phylogram",
                              "phylogram" = 0.012, "cladogram" = 0.012, "fan" = 0.003, 0.012)
      config$rango_desfase <- max(0.05, branch_width * offset_factor)
    }

    # Calculate tree geometric depth to assign a proportional 2.5% offset
    phy_tmp <- phylogeny
    if (is.null(phy_tmp$edge.length)) {
      phy_tmp$edge.length <- rep(1, nrow(phy_tmp$edge))
    }
    max_depth <- max(node.depth.edgelength(phy_tmp))
    config$label_offset <- max_depth * 0.025

    # Dispatch to the corresponding graphics controller (core_render.R)
    # plot_ancestral_reconstruction() internally handles multi_function:
    #   multi_function == 1 → branch overlay
    #   multi_function == 2 → ancestral reconstruction + tip figures
    if (config$mapping_type == 1) {
      plot_simple_mapping(phylogeny, config)
    } else {
      plot_ancestral_reconstruction(phylogeny, config)
    }

  }, error = function(e) {
    cat("Critical error during execution:", e$message, "\n")
    return(invisible(NULL))
  })
}