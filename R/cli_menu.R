################################################################################
# cli_menus.R
#
# Módulo de interfaz interactiva de MultiMapR.
# Contiene toda la lógica de menús de consola: solicitud de colores, caracteres,
# estados y configuración completa de mapeo.
#
# Dependencias (deben existir en el entorno global al cargarse este módulo,
# provistas por data_utils.R):
#   align_tree_data() — alinea tips con datos de caracteres
#   sort_states()     — ordena estados únicos de un carácter
#   is_valid_color()  — valida nombres de color R / hex
#
# Autor: MultiMapR — módulo de interfaz de usuario
################################################################################


# ==============================================================================
# SECCIÓN 1 — UTILIDAD DE SALIDA
# ==============================================================================

#' Verifica si el usuario quiere salir
#'
#' @param input String ingresado por el usuario.
#' @return TRUE si el usuario escribió 'exit', FALSE en caso contrario.
check_exit <- function(input) {
  if (tolower(trimws(input)) == "exit") {
    cat("Saliendo...\n")
    return(TRUE)
  }
  return(FALSE)
}


# ==============================================================================
# SECCIÓN 2 — SELECTORES DE COLORES Y CARACTERES
# ==============================================================================

#' Solicita al usuario un color para un estado dado
#'
#' Repite el prompt hasta que el usuario ingrese un color válido reconocido
#' por R (nombre o hexadecimal) o escriba 'exit' para cancelar.
#'
#' @param estado String con el nombre del estado cuyo color se solicita.
#' @return String con el color elegido, o invisible(NULL) si el usuario cancela.
prompt_state_color <- function(estado) {
  prompt_msg <- paste0("Color para el estado '", estado,
                       "' (nombre o hex, 'exit' para salir): ")
  color <- readline(prompt = prompt_msg)
  if (check_exit(color)) return(invisible(NULL))
  
  while (!is_valid_color(color)) {
    cat("Color no válido. Inténtelo de nuevo.\n")
    color <- readline(prompt = prompt_msg)
    if (check_exit(color)) return(invisible(NULL))
  }
  color
}

#' Solicita al usuario los colores para todos los estados de un carácter
#'
#' Itera sobre cada estado en \code{estados_ordenados} llamando a
#' \code{prompt_state_color()} y construye un named vector estado → color.
#'
#' @param estados_ordenados Vector de estados ya ordenados.
#' @param nombre_caracter   String con el nombre del carácter (para el encabezado).
#' @return Named character vector estado → color.
prompt_character_colors <- function(estados_ordenados, nombre_caracter) {
  cat("\nAsigne colores para el carácter:", nombre_caracter, "\n")
  colores <- vapply(estados_ordenados, prompt_state_color, character(1))
  names(colores) <- estados_ordenados
  colores
}

#' Solicita al usuario que seleccione uno o varios caracteres de una lista
#'
#' Muestra los caracteres disponibles numerados y permite elegir entre
#' \code{min_n} y \code{max_n} de ellos.  Si solo hay un carácter disponible
#' o \code{min_n == max_n == 1}, lo devuelve directamente sin preguntar.
#'
#' @param caracteres_disp Vector de nombres de caracteres disponibles.
#' @param prompt_n        Texto del prompt para pedir el número de caracteres.
#' @param min_n           Mínimo de caracteres a seleccionar (default 1).
#' @param max_n           Máximo de caracteres a seleccionar (default todos).
#' @return Character vector con los nombres de los caracteres seleccionados,
#'         o invisible(NULL) si el usuario cancela.
prompt_characters <- function(caracteres_disp, prompt_n, min_n = 1,
                              max_n = length(caracteres_disp)) {
  cat("\nCaracteres disponibles:\n")
  for (i in seq_along(caracteres_disp))
    cat(sprintf("  %d: %s\n", i, caracteres_disp[i]))
  
  if (min_n == max_n && max_n == 1) {
    if (length(caracteres_disp) == 1) {
      cat("Solo hay un carácter disponible:", caracteres_disp[1], "\n")
      return(caracteres_disp[1])
    }
    sel <- readline(prompt = prompt_n)
    if (check_exit(sel)) return(invisible(NULL))
    sel <- as.integer(sel)
    if (is.na(sel) || sel < 1 || sel > length(caracteres_disp))
      stop("Selección inválida.")
    return(caracteres_disp[sel])
  }
  
  rango_txt <- if (min_n == max_n) as.character(min_n)
  else paste0(min_n, "\u2013", max_n)
  n_str <- readline(prompt = paste0(prompt_n, " (", rango_txt, " o 'exit' para salir): "))
  if (check_exit(n_str)) return(invisible(NULL))
  n <- as.integer(n_str)
  if (is.na(n) || n < min_n || n > max_n)
    stop(paste0("Número inválido. Debe ser entre ", min_n, " y ", max_n, "."))
  
  seleccionados <- character(n)
  for (i in seq_len(n)) {
    s <- readline(prompt = paste0("  Carácter ", i, " (número, o 'exit' para salir): "))
    if (check_exit(s)) return(invisible(NULL))
    s <- as.integer(s)
    if (is.na(s) || s < 1 || s > length(caracteres_disp) ||
        caracteres_disp[s] %in% seleccionados)
      stop("Selección inválida o carácter repetido.")
    seleccionados[i] <- caracteres_disp[s]
  }
  seleccionados
}

#' Solicita al usuario qué estados de un carácter quiere incluir y sus colores
#'
#' Muestra los estados únicos del carácter, permite elegir un subconjunto
#' (o todos con la opción 0) y luego solicita un color para cada uno.
#'
#' @param datos_ord Data.frame de datos ya alineado con la filogenia.
#' @param caracter  Nombre de la columna del carácter a configurar.
#' @return Named character vector estado → color, o invisible(NULL) si cancela.
prompt_states_and_colors <- function(datos_ord, caracter) {
  todos <- sort_states(as.character(unique(datos_ord[[caracter]])))
  
  cat("\nEstados disponibles para '", caracter, "':\n", sep = "")
  for (i in seq_along(todos)) cat(sprintf("  %d: %s\n", i, todos[i]))
  
  cat("  0: Todos los estados\n")
  sel_str <- readline(prompt = paste0(
    "\u00bfQu\u00e9 estados incluir? (0 = todos, o n\u00fameros separados por comas): "))
  if (check_exit(sel_str)) return(invisible(NULL))
  
  if (trimws(sel_str) == "0") {
    estados_sel <- todos
  } else {
    indices <- suppressWarnings(as.integer(unlist(strsplit(sel_str, ","))))
    if (any(is.na(indices)) || any(indices < 1) || any(indices > length(todos)))
      stop("Selección de estados inválida.")
    estados_sel <- todos[indices]
  }
  
  prompt_character_colors(estados_sel, caracter)
}


# ==============================================================================
# SECCIÓN 3 — CONFIGURACIÓN COMPLETA DEL MAPEO
# ==============================================================================

#' Orquesta todos los menús interactivos y devuelve la lista de configuración
#'
#' Recorre los menús de tipo de mapeo, caracteres, colores, algoritmo, grosor,
#' topología, longitudes de rama y opciones de exportación, y empaqueta todas
#' las decisiones del usuario en una lista \code{config} homogénea que consumen
#' las funciones de renderizado de \code{main_MultiMapR.R}.
#'
#' @param filogenia        Objeto \code{phylo} (ape).
#' @param datos_caracteres Data.frame con columna \code{"Species"} y caracteres.
#' @return Lista de configuración, o invisible(NULL) si el usuario cancela.
#' @seealso \code{\link{ejecutar_filogenia}}
setup_mapping_config <- function(filogenia, datos_caracteres) {
  
  if (!inherits(filogenia, "phylo"))
    stop("'filogenia' debe ser un objeto de tipo 'phylo'.")
  if (!is.data.frame(datos_caracteres))
    stop("'datos_caracteres' debe ser un data.frame.")
  
  datos_ord       <- align_tree_data(filogenia, datos_caracteres)
  caracteres_disp <- colnames(datos_ord)[-1]
  
  cat("\n=== Bienvenido a MultiMapR ===\n")
  
  # --- Tipo de mapeo ---------------------------------------------------------
  cat("\nTipo de mapeo:\n")
  cat("  1: Mapeo simple \u2014 figuras de color en los terminales\n")
  cat("  2: Reconstrucci\u00f3n ancestral \u2014 coloreado de ramas\n")
  tipo_mapeo <- readline(prompt = "Seleccione (1/2, o 'exit' para salir): ")
  if (check_exit(tipo_mapeo)) return(invisible(NULL))
  tipo_mapeo <- as.integer(tipo_mapeo)
  if (!tipo_mapeo %in% 1:2) stop("Opción de mapeo inválida.")
  
  config <- list(tipo_mapeo = tipo_mapeo, datos_ord = datos_ord)
  
  # ==========================================================================
  # MAPEO SIMPLE
  # ==========================================================================
  if (tipo_mapeo == 1) {
    
    caracteres_sel <- prompt_characters(
      caracteres_disp,
      prompt_n = "\u00bfCu\u00e1ntos caracteres mapear?",
      min_n = 1, max_n = length(caracteres_disp))
    if (is.null(caracteres_sel)) return(invisible(NULL))
    config$caracteres <- caracteres_sel
    
    colores_por_car <- list()
    for (car in caracteres_sel) {
      cols <- prompt_states_and_colors(datos_ord, car)
      if (is.null(cols)) return(invisible(NULL))
      colores_por_car[[car]] <- cols
    }
    config$colores_por_caracter <- colores_por_car
    
    cat("\nTipo de figura en los terminales:\n")
    cat("  1: C\u00edrculo   (pch = 21)\n")
    cat("  2: Cuadrado  (pch = 22)\n")
    cat("  3: Tri\u00e1ngulo (pch = 24)\n")
    cat("  4: Rombo     (pch = 23)\n")
    fig_str <- readline(prompt = "Seleccione (1\u20134, Enter = 1): ")
    pch_map <- c(`1` = 21L, `2` = 22L, `3` = 24L, `4` = 23L)
    config$pch_figura <- if (nchar(trimws(fig_str)) == 0 || is.na(pch_map[fig_str]))
      21L else pch_map[fig_str]
    
    tam_str <- readline(prompt = "Tama\u00f1o de las figuras (n\u00famero positivo, Enter = 1): ")
    config$tam_figura <- if (nchar(trimws(tam_str)) == 0) 1 else {
      t <- suppressWarnings(as.numeric(tam_str))
      if (is.na(t) || t <= 0) { cat("Valor inv\u00e1lido, usando 1.\n"); 1 } else t
    }
    
    # ==========================================================================
    # RECONSTRUCCIÓN ANCESTRAL
    # ==========================================================================
  } else {
    
    cat("\nModo de visualizaci\u00f3n:\n")
    cat("  1: Superposici\u00f3n en ramas (uno o varios caracteres)\n")
    cat("  2: Figuras en terminales + \u00e1rbol coloreado\n")
    fun_sel <- readline(prompt = "Seleccione (1/2, o 'exit' para salir): ")
    if (check_exit(fun_sel)) return(invisible(NULL))
    fun_sel <- as.integer(fun_sel)
    if (!fun_sel %in% 1:2) stop("Opción inválida.")
    config$funcion_multi <- fun_sel
    
    caracteres_sel <- prompt_characters(
      caracteres_disp,
      prompt_n = "\u00bfCu\u00e1ntos caracteres mapear?",
      min_n = 1L, max_n = min(3L, length(caracteres_disp)))
    if (is.null(caracteres_sel)) return(invisible(NULL))
    config$caracteres <- caracteres_sel
    
    colores_por_car <- list()
    
    if (length(caracteres_sel) == 1) {
      car <- caracteres_sel[1]
      cols <- prompt_states_and_colors(datos_ord, car)
      if (is.null(cols)) return(invisible(NULL))
      colores_por_car[[car]] <- cols
      
      todos_estados <- sort_states(as.character(unique(datos_ord[[car]])))
      config$mapear_todos <- length(cols) == length(todos_estados)
      
    } else {
      for (car in caracteres_sel) {
        # Los estados no seleccionados quedan fuera del mapa de colores y
        # resolver_color() los pintará en gris (gray70) automáticamente.
        cols <- prompt_states_and_colors(datos_ord, car)
        if (is.null(cols)) return(invisible(NULL))
        colores_por_car[[car]] <- cols
      }
      # mapear_todos es TRUE solo si en TODOS los caracteres se eligieron
      # todos los estados; en caso contrario hay estados que irán en gris.
      config$mapear_todos <- all(vapply(caracteres_sel, function(car) {
        todos <- sort_states(as.character(unique(datos_ord[[car]])))
        length(colores_por_car[[car]]) == length(todos)
      }, logical(1)))
    }
    config$colores_por_caracter <- colores_por_car
    
    if (fun_sel == 2) {
      cat("\nTipo de figura en los terminales:\n")
      cat("  1: C\u00edrculo   (pch = 21)\n")
      cat("  2: Cuadrado  (pch = 22)\n")
      cat("  3: Tri\u00e1ngulo (pch = 24)\n")
      cat("  4: Rombo     (pch = 23)\n")
      fig_str <- readline(prompt = "Seleccione (1\u20134, Enter = 1): ")
      if (check_exit(fig_str)) return(invisible(NULL))
      pch_map <- c(`1` = 21L, `2` = 22L, `3` = 24L, `4` = 23L)
      config$pch_figura <- if (nchar(trimws(fig_str)) == 0 || is.na(pch_map[fig_str]))
        21L else pch_map[fig_str]
      
      tam_str <- readline(prompt = "Tama\u00f1o de las figuras (n\u00famero positivo, Enter = 1): ")
      if (check_exit(tam_str)) return(invisible(NULL))
      config$tam_figura <- if (nchar(trimws(tam_str)) == 0) 1 else {
        t <- suppressWarnings(as.numeric(tam_str))
        if (is.na(t) || t <= 0) { cat("Valor inv\u00e1lido, usando 1.\n"); 1 } else t
      }
    }
    
    cat("\nAlgoritmo de reconstrucci\u00f3n ancestral:\n")
    cat("  1: Por defecto (mayor\u00eda ponderada por profundidad)\n")
    cat("  2: Fitch (ACCTRAN/DELTRAN/Unambiguous \u2014 se carga autom\u00e1ticamente desde fitch.R)\n")
    algo_sel <- readline(prompt = "Seleccione (1/2, o 'exit' para salir): ")
    if (check_exit(algo_sel)) return(invisible(NULL))
    algo_sel <- as.integer(algo_sel)
    if (!algo_sel %in% 1:2) stop("Opción de algoritmo inválida.")
    config$algoritmo <- algo_sel
    
    if (algo_sel == 2) {
      cat("\nModo de optimizaci\u00f3n Fitch:\n")
      cat("  1: ACCTRAN     (cambios acelerados \u2014 hacia las hojas)\n")
      cat("  2: DELTRAN     (cambios retrasados \u2014 hacia la ra\u00edz)\n")
      cat("  3: Unambiguous (solo estados sin ambig\u00fcedad tras los dos pasos)\n")
      modo_str <- readline(prompt = "Seleccione (1/2/3, Enter = 1): ")
      if (check_exit(modo_str)) return(invisible(NULL))
      modo_fitch <- if (nchar(trimws(modo_str)) == 0) 1L else as.integer(modo_str)
      if (is.na(modo_fitch) || !modo_fitch %in% 1:3) {
        cat("Opci\u00f3n inv\u00e1lida, usando ACCTRAN.\n")
        modo_fitch <- 1L
      }
      config$modo_fitch <- modo_fitch
    }
  }
  
  # --- Grosor de ramas -------------------------------------------------------
  cat("\nGrosor de las ramas (n\u00famero positivo; Enter = 2):\n")
  cat("  Referencia: fino \u2248 1, normal \u2248 2, grueso \u2248 4, muy grueso \u2248 8\n")
  gros_str <- readline(prompt = "Grosor: ")
  if (check_exit(gros_str)) return(invisible(NULL))
  gros_val <- suppressWarnings(as.numeric(trimws(gros_str)))
  if (is.na(gros_val) || gros_val <= 0) {
    cat("Valor inv\u00e1lido, usando 2.\n")
    gros_val <- 2
  }
  config$grosor <- gros_val
  
  # --- Tipo de árbol ---------------------------------------------------------
  tipos_validos <- c("phylogram", "cladogram", "fan")
  cat("\nTipo de \u00e1rbol: phylogram / cladogram / fan\n")
  tipo_arbol <- tolower(trimws(readline(prompt = "Seleccione (o 'exit' para salir): ")))
  if (check_exit(tipo_arbol)) return(invisible(NULL))
  while (!tipo_arbol %in% tipos_validos) {
    cat("Tipo inv\u00e1lido.\n")
    tipo_arbol <- tolower(trimws(readline(prompt = "Seleccione: ")))
    if (check_exit(tipo_arbol)) return(invisible(NULL))
  }
  config$tipo_arbol <- tipo_arbol
  
  # --- Calcular rango_desfase dinámico en función del grosor -----------------
  # El offset entre ramas debe crecer proporcionalmente al grosor para que
  # nunca se tapen.  El factor de conversión varía por tipo de árbol porque
  # las coordenadas de datos difieren en escala:
  #   · phylogram / cladogram: el eje Y tiene escala lineal (1 unidad ≈ 1 tip),
  #     la relación empírica grosor→desfase es ~0.012 unidades/pt.
  #   · fan: las coordenadas son polares; el spread es un porcentaje del radio,
  #     por lo que se necesita un factor mayor (~0.0025 × R_max por pt).
  #     Sin embargo, como R_max se desconoce aquí, se usa un factor de escala
  #     relativo que se ajusta automáticamente en .emtree_render_fan().
  factor_offset <- switch(tipo_arbol,
                          "phylogram"  = 0.012,
                          "cladogram"  = 0.012,
                          "fan"        = 0.003,
                          0.012
  )
  config$rango_desfase <- max(0.05, gros_val * factor_offset)
  
  # --- Longitud de Ramas -----------------------------------------------------
  if (!is.null(filogenia$edge.length) && tipo_arbol %in% c("phylogram", "fan")) {
    cat("\nEl \u00e1rbol tiene longitudes de rama.\n")
    cat("  1: Usar longitudes proporcionales al ploteo\n")
    cat("  2: Ignorar longitudes (ramas uniformes, estilo cladograma)\n")
    len_str <- readline(prompt = "Seleccione (1/2, Enter = 1): ")
    if (check_exit(len_str)) return(invisible(NULL))
    config$use_edge_length <- if (trimws(len_str) == "2") FALSE else TRUE
  } else {
    config$use_edge_length <- FALSE
  }
  
  # --- Exportar --------------------------------------------------------------
  cat("\n\u00bfExportar el gr\u00e1fico?\n")
  cat("  1: S\u00ed \u2014 guardar como archivo\n")
  cat("  2: No \u2014 mostrar en la ventana de R\n")
  exp_str <- readline(prompt = "Seleccione (1/2, Enter = 2): ")
  if (check_exit(exp_str)) return(invisible(NULL))
  exportar <- trimws(exp_str) == "1"
  
  if (exportar) {
    
    # Formato de salida
    cat("\nFormato de salida:\n")
    cat("  1: PNG (recomendado para pantalla)\n")
    cat("  2: PDF (vectorial, ideal para publicaciones)\n")
    fmt_str <- readline(prompt = "Seleccione (1/2, Enter = 1): ")
    if (check_exit(fmt_str)) return(invisible(NULL))
    config$export_format <- if (trimws(fmt_str) == "2") "pdf" else "png"
    
    # Nombre de archivo sin extensión
    nombre_default <- paste0("MultiMapR_",
                             paste(config$caracteres, collapse = "-"), "_",
                             config$tipo_arbol)
    fn_str <- readline(prompt = paste0(
      "Nombre del archivo sin extensi\u00f3n (Enter = '", nombre_default, "'): "))
    if (check_exit(fn_str)) return(invisible(NULL))
    fn <- trimws(fn_str)
    fn <- sub("\\.(png|pdf)$", "", fn, ignore.case = TRUE)
    config$export_filename <- if (nzchar(fn)) fn else nombre_default
    
    # ── MENÚ DE DIMENSIONES PERSONALIZADAS ───────────────────────────────────
    # La opción "Automáticas" conserva el comportamiento original (alto ∝ Ntip).
    cat("\nDimensiones del archivo exportado:\n")
    cat("  1: Autom\u00e1ticas (recomendado \u2014 alto proporcional al n\u00famero de terminales)\n")
    cat("  2: Personalizadas (ancho y alto en pulgadas)\n")
    dim_str <- readline(prompt = "Seleccione (1/2, Enter = 1): ")
    if (check_exit(dim_str)) return(invisible(NULL))
    
    if (trimws(dim_str) == "2") {
      w_str <- readline(prompt = "Ancho en pulgadas (n\u00famero positivo, ej. 14): ")
      if (check_exit(w_str)) return(invisible(NULL))
      w_val <- suppressWarnings(as.numeric(trimws(w_str)))
      if (is.na(w_val) || w_val <= 0) {
        cat("Valor inv\u00e1lido para el ancho. Se usar\u00e1n dimensiones autom\u00e1ticas.\n")
        config$width  <- NULL
        config$height <- NULL
      } else {
        h_str <- readline(prompt = "Alto en pulgadas (n\u00famero positivo, ej. 20): ")
        if (check_exit(h_str)) return(invisible(NULL))
        h_val <- suppressWarnings(as.numeric(trimws(h_str)))
        if (is.na(h_val) || h_val <= 0) {
          cat("Valor inv\u00e1lido para el alto. Se usar\u00e1n dimensiones autom\u00e1ticas.\n")
          config$width  <- NULL
          config$height <- NULL
        } else {
          config$width  <- w_val
          config$height <- h_val
          cat(sprintf("  \u2192 Exportando a %.1f \u00d7 %.1f pulgadas.\n", w_val, h_val))
        }
      }
    } else {
      config$width  <- NULL
      config$height <- NULL
    }
    # ── FIN MENÚ DE DIMENSIONES ───────────────────────────────────────────────
    
    # ── MENÚ DE ESQUINA DE LEYENDA ────────────────────────────────────────────
    # El espacio para la leyenda se reserva ampliando el margen de la esquina
    # elegida, de modo que nunca choque con etiquetas ni con el gráfico.
    cat("\n\u00bfEn qu\u00e9 esquina desea colocar la leyenda?\n")
    cat("  1: Inferior izquierda  (bottomleft)   \u2014 recomendado para phylogram\n")
    cat("  2: Inferior derecha    (bottomright)\n")
    cat("  3: Superior izquierda  (topleft)\n")
    cat("  4: Superior derecha    (topright)      \u2014 \u00fatil en fan con ramas densas\n")
    esquina_map <- c("1" = "bottomleft", "2" = "bottomright",
                     "3" = "topleft",    "4" = "topright")
    esq_str <- readline(prompt = "Seleccione (1\u20134, Enter = 1): ")
    if (check_exit(esq_str)) return(invisible(NULL))
    config$legend_corner <- if (nzchar(trimws(esq_str)) && trimws(esq_str) %in% names(esquina_map))
      esquina_map[[trimws(esq_str)]] else "bottomleft"
    cat(sprintf("  \u2192 Leyenda en: %s\n", config$legend_corner))
    # ── FIN MENÚ DE ESQUINA DE LEYENDA ────────────────────────────────────────
    
  } else {
    config$export_filename <- NULL
    config$export_format   <- NULL
    config$width           <- NULL
    config$height          <- NULL
    config$legend_corner   <- NULL
  }
  
  config
}