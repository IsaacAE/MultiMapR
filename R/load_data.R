################################################################################
# load_data.R
#
# Utilidades para cargar y preparar datos filogenéticos en MultiMapR.
#
# FUNCIÓN PRINCIPAL:
#   load_data(ruta_arbol, ruta_csv, ...)  →  list(arbol, caracteres)
#
# FORMATOS DE ÁRBOL SOPORTADOS:
#   .tre / .tree / .nwk  →  read.tree()   (Newick)
#   .nex / .nexus        →  read.nexus()  (NEXUS)
#   (detección automática por extensión; se puede forzar con formato_arbol=)
#
# FORMATOS DE CSV SOPORTADOS:
#   - Con o sin encabezado en los caracteres (detección automática)
#   - Primera columna: nombres de especies (cualquier nombre de columna)
#   - Caracteres numéricos y de texto (strings)
#   - Inaplicables "-" y desconocidos "?" se preservan tal cual
#   - Nombres de especies con espacios → se normalizan a "_" opcionalmente
#   - BOM UTF-8 (archivos exportados desde Excel) eliminado automáticamente
#
# NOTA: Este archivo forma parte del paquete MultiMapR. Las dependencias
# (ape, tools) se declaran en DESCRIPTION; no se usa library() aquí.
################################################################################


# ============================================================================ #
# MANEJO DE BOM
# ============================================================================ #

#' Elimina el BOM UTF-8 (EF BB BF) si existe y devuelve la ruta a usar
#'
#' Si el archivo tiene BOM, copia el contenido sin los 3 bytes iniciales a un
#' archivo temporal y devuelve su ruta. Si no tiene BOM devuelve la ruta original.
#' El archivo temporal se elimina automáticamente al terminar la sesión de R.
#'
#' @param ruta_csv Ruta al archivo original.
#' @return Ruta al archivo limpio (temporal o la original).
.quitar_bom <- function(ruta_csv) {
  con <- file(ruta_csv, open = "rb")
  bom <- readBin(con, raw(), n = 3)
  tiene_bom <- (length(bom) == 3 &&
                  bom[1] == as.raw(0xEF) &&
                  bom[2] == as.raw(0xBB) &&
                  bom[3] == as.raw(0xBF))

  if (!tiene_bom) {
    close(con)
    return(ruta_csv)
  }

  # Los 3 bytes del BOM ya fueron consumidos; leer el resto
  contenido <- readBin(con, raw(), n = file.info(ruta_csv)$size)
  close(con)

  tmp <- tempfile(fileext = ".csv")
  con_out <- file(tmp, open = "wb")
  writeBin(contenido, con_out)
  close(con_out)

  message("BOM UTF-8 detectado y eliminado para la lectura.")
  tmp
}


# ============================================================================ #
# DETECCIÓN AUTOMÁTICA DE ENCABEZADO
# ============================================================================ #

#' Decide si un CSV tiene encabezado en las columnas de caracteres
#'
#' Heurística basada en la columna de especies (columna 1):
#' los nombres de taxones siempre contienen espacio o guion bajo
#' ("Homo_sapiens", "Homo sapiens"). Si la primera celda no tiene
#' ninguno de los dos se trata como etiqueta de columna -> hay encabezado.
#'
#'   "Species"               -> header = TRUE   (sin espacio ni "_")
#'   "Taxon"                 -> header = TRUE
#'   "Saccopteryx bilineata" -> header = FALSE  (tiene espacio)
#'   "Passer_domesticus"     -> header = FALSE  (tiene "_")
#'
#' @param ruta_csv  Ruta al archivo CSV (ya sin BOM).
#' @param sep       Separador de campos (default ",").
#' @return TRUE si se detecta encabezado, FALSE si no.
.detectar_header <- function(ruta_csv, sep = ",") {
  primera <- read.csv(ruta_csv, header = FALSE, sep = sep,
                      stringsAsFactors = FALSE, nrows = 1,
                      colClasses = "character")

  if (nrow(primera) < 1 || ncol(primera) < 1) return(FALSE)

  celda_sp <- trimws(as.character(primera[1, 1]))

  # Un nombre de taxon tiene espacio o guion bajo ("Homo_sapiens",
  # "Homo sapiens"). Una etiqueta de columna ("Species", "Taxon") no.
  parece_taxon <- grepl("[ _]", celda_sp)
  !parece_taxon
}


# ============================================================================ #
# NORMALIZACIÓN DE NOMBRES DE COLUMNAS
# ============================================================================ #

#' Genera nombres de caracteres seguros para R
#'
#' Si el CSV tiene encabezado usa esos nombres (sanitizados).
#' Si no tiene encabezado genera "char1", "char2", …
#'
#' @param nombres_raw  Vector con los nombres crudos de las columnas de caracteres.
#'                     NULL o vacío → se generan nombres automáticos.
#' @param n            Número de columnas de caracteres.
#' @return Character vector de longitud n.
.normalizar_nombres_char <- function(nombres_raw, n) {
  if (is.null(nombres_raw) || length(nombres_raw) == 0) {
    return(paste0("char", seq_len(n)))
  }
  nombres <- make.names(trimws(nombres_raw), unique = TRUE)
  vacios  <- nchar(nombres) == 0 | grepl("^V[0-9]+$", nombres)
  nombres[vacios] <- paste0("char", which(vacios))
  nombres
}


# ============================================================================ #
# LECTURA DEL ÁRBOL
# ============================================================================ #

#' Sanitiza un archivo de árbol para compatibilidad con ape
#'
#' Elimina problemas comunes producidos por distintos programas filogenéticos:
#'   - Saltos de línea Windows (\r\n → \n)
#'   - Comentarios numerados en TRANSLATE de Mesquite (`[0]`, `[1]`, ...)
#'     que ape no puede parsear
#'
#' @param ruta_arbol Ruta al archivo original.
#' @return Ruta al archivo limpio (temporal si hubo cambios, original si no).
.sanitizar_arbol <- function(ruta_arbol) {
  lineas <- readLines(ruta_arbol, warn = FALSE)

  # 1) readLines ya maneja \r\n en la mayoría de plataformas,
  #    pero forzamos limpieza explícita por si acaso
  lineas <- gsub("\r", "", lineas)

  # 2) Eliminar comentarios numéricos al inicio de línea: "[0]", "[12]", etc.
  #    Patrón: línea que empieza con espacios/tabs opcionales seguidos de [N]
  lineas <- gsub("^(\\s*)\\[[0-9]+\\]\\s*", "\\1", lineas, perl = TRUE)

  tmp <- tempfile(fileext = tools::file_ext(ruta_arbol))
  writeLines(lineas, tmp)
  tmp
}


#' Reads a phylogenetic tree file (Newick or NEXUS)
#'
#' Compatible con archivos de Mesquite, TNT, WinClada y otros programas
#' que producen variantes no estándar del formato NEXUS/Newick.
#'
#' @param ruta_arbol    Ruta al archivo de árbol.
#' @param formato_arbol "auto" (default), "newick" o "nexus".
#' @return Objeto phylo.
read_tree <- function(ruta_arbol, formato_arbol = "auto") {
  if (!file.exists(ruta_arbol))
    stop(paste0("No se encontró el archivo de árbol: '", ruta_arbol, "'"))

  fmt <- tolower(trimws(formato_arbol))

  if (fmt == "auto") {
    ext <- tolower(tools::file_ext(ruta_arbol))
    fmt <- if (ext %in% c("nex", "nexus")) "nexus" else "newick"
  }

  # Limpiar el archivo antes de parsear
  ruta_limpia <- .sanitizar_arbol(ruta_arbol)

  arbol <- tryCatch({
    if (fmt == "nexus") read.nexus(ruta_limpia) else read.tree(ruta_limpia)
  }, error = function(e) {
    # Si falla con el formato detectado, intentar el otro
    fmt_alt <- if (fmt == "nexus") "newick" else "nexus"
    message(sprintf(
      "Formato '%s' falló (%s). Intentando '%s'...", fmt, e$message, fmt_alt
    ))
    tryCatch(
      if (fmt_alt == "nexus") read.nexus(ruta_limpia) else read.tree(ruta_limpia),
      error = function(e2) stop(paste0(
        "No se pudo leer el árbol '", ruta_arbol, "'.\n",
        "  Como ", fmt,     ": ", e$message, "\n",
        "  Como ", fmt_alt, ": ", e2$message
      ))
    )
  })

  if (is.null(arbol))
    stop(paste0("No se pudo leer el árbol desde '", ruta_arbol, "'."))

  arbol
}


# ============================================================================ #
# LECTURA DEL CSV
# ============================================================================ #

#' Reads a character CSV file and prepares it for MultiMapR
#'
#' @param ruta_csv            Ruta al archivo CSV.
#' @param sep                 Separador de campos (default ",").
#' @param col_especies        Índice o nombre de la columna de especies (default 1).
#' @param normalizar_espacios Si TRUE, reemplaza espacios por "_" en los nombres
#'                            de especies (default FALSE).
#' @return Data.frame con columna "Species" y columnas "charN" o nombres propios.
read_csv_characters <- function(ruta_csv,
                                sep                  = ",",
                                col_especies         = 1,
                                normalizar_espacios  = FALSE) {

  if (!file.exists(ruta_csv))
    stop(paste0("No se encontró el archivo CSV: '", ruta_csv, "'"))

  # Eliminar BOM si existe (Excel UTF-8 lo agrega); trabajar sobre archivo limpio
  ruta_limpia  <- .quitar_bom(ruta_csv)
  tiene_header <- .detectar_header(ruta_limpia, sep = sep)

  datos_raw <- read.csv(ruta_limpia,
                        header           = tiene_header,
                        sep              = sep,
                        stringsAsFactors = FALSE,
                        colClasses       = "character",
                        check.names      = FALSE)

  if (ncol(datos_raw) < 2)
    stop("El CSV debe tener al menos dos columnas: especies + un carácter.")

  # --- Extraer columna de especies ---
  if (is.character(col_especies)) {
    if (!col_especies %in% colnames(datos_raw))
      stop(paste0("No se encontró la columna de especies '", col_especies, "'."))
    idx_sp <- which(colnames(datos_raw) == col_especies)
  } else {
    idx_sp <- as.integer(col_especies)
    if (idx_sp < 1 || idx_sp > ncol(datos_raw))
      stop(paste0("Índice de columna de especies fuera de rango: ", idx_sp))
  }

  especies <- trimws(datos_raw[[idx_sp]])
  if (normalizar_espacios)
    especies <- gsub(" ", "_", especies)

  # --- Columnas de caracteres (todo excepto la columna de especies) ---
  idx_char   <- setdiff(seq_len(ncol(datos_raw)), idx_sp)
  datos_char <- datos_raw[, idx_char, drop = FALSE]

  # Determinar nombres de columnas de caracteres
  nombres_char <- if (tiene_header) colnames(datos_raw)[idx_char] else NULL
  nombres_char <- .normalizar_nombres_char(nombres_char, ncol(datos_char))

  # Construir data.frame final
  resultado           <- as.data.frame(datos_char, stringsAsFactors = FALSE)
  colnames(resultado) <- nombres_char
  resultado           <- data.frame(Species = especies, resultado,
                                    stringsAsFactors = FALSE,
                                    check.names      = FALSE)
  rownames(resultado) <- NULL

  message(sprintf(
    "CSV cargado: %d especies, %d carácter(es). Encabezado %s.",
    nrow(resultado),
    length(nombres_char),
    if (tiene_header) "detectado (nombres propios)" else "no detectado (nombres automáticos)"
  ))

  resultado
}


# ============================================================================ #
# VALIDACIÓN DE COMPATIBILIDAD ÁRBOL ↔ CSV
# ============================================================================ #

#' Verifica que los nombres del CSV coincidan con los tip.labels del árbol
#'
#' @param arbol    Objeto phylo.
#' @param datos    Data.frame con columna "Species".
#' @param estricto Si TRUE lanza error ante cualquier discrepancia;
#'                 si FALSE (default) solo advierte y devuelve los no emparejados.
#' @return Invisible: vector de nombres del CSV que no están en el árbol.
validar_compatibilidad <- function(arbol, datos, estricto = FALSE) {
  tips          <- arbol$tip.label
  sp            <- datos$Species
  sin_par_csv   <- setdiff(sp,   tips)
  sin_par_arbol <- setdiff(tips, sp)
  ok            <- length(sin_par_csv) == 0 && length(sin_par_arbol) == 0

  if (!ok) {
    msg_partes <- character(0)
    if (length(sin_par_csv) > 0)
      msg_partes <- c(msg_partes,
                      paste0("En CSV pero NO en árbol (", length(sin_par_csv), "): ",
                             paste(head(sin_par_csv, 10), collapse = ", "),
                             if (length(sin_par_csv) > 10) " ..." else ""))
    if (length(sin_par_arbol) > 0)
      msg_partes <- c(msg_partes,
                      paste0("En árbol pero NO en CSV (", length(sin_par_arbol), "): ",
                             paste(head(sin_par_arbol, 10), collapse = ", "),
                             if (length(sin_par_arbol) > 10) " ..." else ""))

    msg <- paste(c("Discrepancias entre árbol y CSV:", msg_partes), collapse = "\n  ")
    if (estricto) stop(msg) else warning(msg)
  } else {
    message("Validación OK: los ", length(tips), " tips coinciden con el CSV.")
  }

  invisible(sin_par_csv)
}


# ============================================================================ #
# FUNCIÓN PRINCIPAL
# ============================================================================ #

#' Loads a tree and a character CSV ready to use with MultiMapR
#'
#' Detecta automáticamente:
#'   - Formato del árbol (Newick / NEXUS) por extensión.
#'   - Presencia o ausencia de encabezado en el CSV.
#'   - BOM UTF-8 en el CSV (archivos exportados desde Excel).
#'   - Tipos de caracteres (numéricos o strings): se preservan como texto.
#'
#' @param ruta_arbol          Ruta al archivo de árbol (.tre, .nwk, .nex, .nexus, …).
#' @param ruta_csv            Ruta al archivo CSV de caracteres.
#' @param sep                 Separador del CSV (default ",").
#' @param col_especies        Columna de especies: índice entero o nombre (default 1).
#' @param normalizar_espacios Si TRUE convierte espacios a "_" en nombres de
#'                            especies (útil cuando árbol usa "_" y CSV usa " ").
#' @param formato_arbol       "auto" (default), "newick" o "nexus".
#' @param estricto            Si TRUE, error ante discrepancias árbol/CSV (default FALSE).
#'
#' @return Lista con:
#'   \item{arbol}{Objeto phylo listo para ape / MultiMapR.}
#'   \item{caracteres}{Data.frame con columna "Species" y columnas de caracteres.}
#'
#' @examples
#' \dontrun{
#'   # CSV sin encabezado (caracteres numéricos)
#'   d1 <- load_data("tree2_jadc.tre", "Matriz_JAIR.csv")
#'   execute_phylogeny(d1$arbol, d1$caracteres)
#'
#'   # CSV con BOM y nombres de especie con espacios
#'   d2 <- load_data("bats_tre.tre", "bats_matrix.csv",
#'                   normalizar_espacios = TRUE)
#'   execute_phylogeny(d2$arbol, d2$caracteres)
#'
#'   # Separador punto y coma, columna de especies por nombre
#'   d3 <- load_data("arbol.nex", "datos.csv",
#'                   sep = ";", col_especies = "Taxon")
#'   execute_phylogeny(d3$arbol, d3$caracteres)
#'
#'   # Uso polimórfico: pasar rutas directamente al orquestador
#'   execute_phylogeny("mi_arbol.tre", "mis_caracteres.csv")
#' }
load_data <- function(ruta_arbol,
                      ruta_csv,
                      sep                  = ",",
                      col_especies         = 1,
                      normalizar_espacios  = FALSE,
                      formato_arbol        = "auto",
                      estricto             = FALSE) {

  arbol      <- read_tree(ruta_arbol, formato_arbol = formato_arbol)
  caracteres <- read_csv_characters(ruta_csv,
                                    sep                 = sep,
                                    col_especies        = col_especies,
                                    normalizar_espacios = normalizar_espacios)

  validar_compatibilidad(arbol, caracteres, estricto = estricto)

  list(arbol = arbol, caracteres = caracteres)
}
