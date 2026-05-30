#' Orquestador principal de MultiMapR
#'
#' @param filogenia Ruta al archivo de árbol o un objeto phylo.
#' @param datos_caracteres Ruta al archivo CSV o un data.frame.
#' @param grosor Grosor opcional de las ramas.
#' @param label_offset Desfase dinámico de las etiquetas.
#' @import ape
#' @export
.exportar_dispositivo <- function(filename, format, filogenia, tipo_arbol, expr,
                                  width = NULL, height = NULL) {
  n_tips        <- Ntip(filogenia)
  alto_default  <- n_tips * 0.25 + 2
  ancho_default <- 12
  if (tipo_arbol == "fan") {
    alto_default  <- max(alto_default, ancho_default)
    ancho_default <- alto_default
  }
  alto_in  <- if (!is.null(height)) height else alto_default
  ancho_in <- if (!is.null(width))  width  else ancho_default

  dir.create("Exports", showWarnings = FALSE, recursive = TRUE)
  ruta <- file.path("Exports", paste0(basename(filename), ".", format))
  if (format == "pdf") {
    pdf(file = ruta, width = ancho_in, height = alto_in)
  } else {
    png(filename = ruta,
        width    = ancho_in,
        height   = alto_in,
        units    = "in",
        res      = 300,
        bg       = "white")
  }
  tryCatch(
    force(expr),
    error   = function(e) stop(sprintf("[MultiMapR] Error exportando %s: %s",
                                       toupper(format), conditionMessage(e))),
    finally = dev.off()
  )
  dims_origen <- if (!is.null(width) || !is.null(height)) " [custom]" else " [auto]"
  message(sprintf("[MultiMapR] Archivo exportado: %s  (%.1f \u00d7 %.1f in%s%s)",
                  ruta, ancho_in, alto_in,
                  if (format == "png") " @ 300 dpi" else "",
                  dims_origen))
  invisible(ruta)
}


# ==============================================================================
# MAIN FUNCTION — ORCHESTRATOR
# ==============================================================================

#' Orquestador principal de MultiMapR
#'
#' Soporta polimorfismo de argumentos: si `filogenia` y `datos_caracteres`
#' son cadenas de texto (rutas a archivos en disco), los datos se cargan
#' automáticamente a través de `load_data()` antes de configurar el mapeo.
#'
#' @param filogenia        Objeto phylo — O — ruta (string) al archivo de árbol.
#' @param datos_caracteres Data.frame con columna "Species" — O — ruta (string)
#'                         al archivo CSV de caracteres.
#' @param grosor           Grosor de ramas (default NULL = el usuario lo elige
#'                         en el menú). Un número positivo sobreescribe el menú.
#' @param label_offset     Distancia entre terminales y etiquetas (default 0.3).
#' @return Invisible NULL.
#' @export
execute_phylogeny <- function(filogenia, datos_caracteres, grosor = NULL, label_offset = 0.3) {
  tryCatch({

    # POLIMORFISMO: Si se pasan rutas de archivos (character), cargar los datos automáticamente
    if (is.character(filogenia) && is.character(datos_caracteres)) {
      cat("\n[MultiMapR] Detectadas rutas de archivos. Cargando datos automáticamente...\n")
      datos_cargados   <- load_data(ruta_arbol = filogenia, ruta_csv = datos_caracteres)
      filogenia        <- datos_cargados$arbol
      datos_caracteres <- datos_cargados$caracteres
    }

    # Continuar con el flujo normal del sistema usando la interfaz modular
    config <- setup_mapping_config(filogenia, datos_caracteres)
    if (is.null(config)) return(invisible(NULL))

    # Forzar ramas uniformes si el usuario lo solicitó en el menú
    if (isFALSE(config$use_edge_length)) {
      filogenia$edge.length <- NULL
    }

    # Configuración externa de grosor si se pasa por argumento
    if (!is.null(grosor) && is.numeric(grosor) && grosor > 0) {
      config$grosor <- grosor
      factor_offset <- switch(config$tipo_arbol %||% "phylogram",
                              "phylogram" = 0.012, "cladogram" = 0.012, "fan" = 0.003, 0.012)
      config$rango_desfase <- max(0.05, grosor * factor_offset)
    }

    # Calcular profundidad geométrica del árbol para asignar un offset proporcional del 2.5%
    phy_tmp <- filogenia
    if (is.null(phy_tmp$edge.length)) {
      phy_tmp$edge.length <- rep(1, nrow(phy_tmp$edge))
    }
    profundidad_max <- max(node.depth.edgelength(phy_tmp))
    config$label_offset <- profundidad_max * 0.025

    # Despachar al controlador gráfico correspondiente (core_render.R)
    if (config$tipo_mapeo == 1) {
      plot_simple_mapping(filogenia, config)
    } else {
      if (config$funcion_multi == 2) {
        plot_multicharacter_figure(filogenia, config)
      } else {
        plot_ancestral_reconstruction(filogenia, config)
      }
    }

  }, error = function(e) {
    cat("Error crítico en la ejecución:", e$message, "\n")
    return(invisible(NULL))
  })
}

