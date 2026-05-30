################################################################################
# Función universal de exportación para MultiMapR.
# Soporta tres topologías  : phylogram | cladogram | fan
# Soporta dos formatos      : png | pdf
#
# Estrategia geométrica:
#   - Canvas en Blanco (edge.color = "transparent") para poblar .PlotPhyloEnv
#     sin renderizar nada, y luego dibujar las historias manualmente con
#     segments() / lines() + desfases calculados al vuelo.
#   - Fan usa geometría polar real: radio + ángulo → segmentos radiales y arcos.
#
# Dependencias: ape (>= 5.0) — declarada en DESCRIPTION, no usar library() aquí
# Autor:        MultiMapR — módulo de exportación universal
################################################################################

# ==============================================================================
# VALIDADORES INTERNOS
# ==============================================================================

#' Verifica que `tree` sea un objeto phylo válido con aristas y tip.labels
.emtree_validar_arbol <- function(tree) {
  if (!inherits(tree, "phylo"))
    stop("`tree` debe ser un objeto de clase 'phylo'.")
  if (is.null(tree$edge) || nrow(tree$edge) == 0L)
    stop("`tree$edge` está vacío: el árbol no tiene ramas.")
  if (is.null(tree$tip.label))
    stop("`tree$tip.label` es NULL: el árbol no tiene etiquetas de terminales.")
}

#' Verifica que `color_list` sea coherente con el número de aristas del árbol
.emtree_validar_color_list <- function(color_list, n_edges) {
  if (!is.list(color_list) || length(color_list) == 0L)
    stop("`color_list` debe ser una lista no vacía de vectores de colores.")
  for (i in seq_along(color_list)) {
    vec <- color_list[[i]]
    if (!is.character(vec))
      stop(sprintf("`color_list[[%d]]` debe ser un vector de caracteres.", i))
    if (length(vec) != n_edges)
      stop(sprintf(
        "`color_list[[%d]]` tiene longitud %d, pero el árbol tiene %d ramas. ",
        i, length(vec), n_edges),
        "Cada vector debe tener exactamente una entrada por arista.")
  }
}

#' Verifica argumentos escalares de tipo y formato
.emtree_validar_tipo_formato <- function(type, format) {
  tipos_ok   <- c("phylogram", "cladogram", "fan")
  formatos_ok <- c("png", "pdf")
  if (!is.character(type)   || length(type) != 1L || !type   %in% tipos_ok)
    stop(sprintf("`type` debe ser uno de: %s.", paste(tipos_ok, collapse = ", ")))
  if (!is.character(format) || length(format) != 1L || !format %in% formatos_ok)
    stop(sprintf("`format` debe ser uno de: %s.", paste(formatos_ok, collapse = ", ")))
}

#' Verifica que `filename` sea un string no vacío (sin extensión requerida)
.emtree_validar_filename <- function(filename) {
  if (!is.character(filename) || length(filename) != 1L || !nzchar(trimws(filename)))
    stop("`filename` debe ser un string no vacío (sin extensión).")
}


# ==============================================================================
# HELPERS GEOMÉTRICOS INTERNOS
# ==============================================================================

#' Recupera .PlotPhyloEnv desde el namespace de ape de forma segura
.emtree_get_PlotPhyloEnv <- function() {
  get(".PlotPhyloEnv", envir = asNamespace("ape"))
}

#' Abre el dispositivo gráfico adecuado según `format`
#' @return Invisible NULL.  El dispositivo queda abierto tras la llamada.
.emtree_abrir_dispositivo <- function(filename, format, ancho_in, alto_in) {

  ruta <- paste0(filename, ".", format)
  if (format == "png") {
    png(filename = ruta,
        width    = ancho_in,
        height   = alto_in,
        units    = "in",
        res      = 300,
        bg       = "white")
  } else {                       # pdf
    pdf(file   = ruta,
        width  = ancho_in,
        height = alto_in)
  }
  invisible(NULL)
}

#' Calcula el vector de desfases a partir de los parámetros físicos actuales
#' del dispositivo (par("usr") / par("pin")).
#'
#' @param N             Número de historias.
#' @param rango_desfase Amplitud del rango en Y (default 0.1).
#' @return Lista con `dx` y `dy` (vectores de longitud N).
.emtree_calcular_desfases <- function(N, rango_desfase = 0.1) {
  if (N == 1L)
    return(list(dx = 0, dy = 0))

  # Paso fijo = rango_desfase entre capas consecutivas, conjunto centrado en 0.
  # Esto garantiza que para cualquier N la separación visual entre ramas
  # sea siempre la misma (rango_desfase), sin huecos ni asimetrías:
  #   N=2 → c(-rango/2, +rango/2)  separación = rango ✓
  #   N=3 → c(-rango,   0, +rango) separación = rango ✓  (igual que seq() daba)
  #   N=4 → c(-3r/2, -r/2, +r/2, +3r/2) separación = rango ✓
  #
  # El bug original: seq(-r, r, N=2) → c(-r, +r), separación = 2*rango y
  # ningún valor en 0, dejando un hueco del doble del paso de N=3.
  indices    <- seq_len(N) - (N + 1L) / 2   # e.g. N=2→c(-0.5,0.5); N=3→c(-1,0,1)
  desfases_y <- indices * rango_desfase

  usr <- par("usr")   # c(x1, x2, y1, y2) — límites del sistema de coordenadas
  pin <- par("pin")   # c(ancho_in, alto_in) — dimensiones físicas del área de gráfico

  escala_x <- (usr[2L] - usr[1L]) / pin[1L]   # unidades_datos / pulgada (X)
  escala_y <- (usr[4L] - usr[3L]) / pin[2L]   # unidades_datos / pulgada (Y)

  desfases_x <- desfases_y * (escala_x / escala_y)

  list(dx = desfases_x, dy = desfases_y)
}


# ==============================================================================
# HELPER DE LEYENDA — ESPACIO EXCLUSIVO
# ==============================================================================

#' Convierte la esquina elegida por el usuario en márgenes y posición de leyenda
#'
#' Estrategia: se amplía el margen de la esquina elegida para que la leyenda
#' nunca se superponga con etiquetas ni con el gráfico. La posición devuelta
#' por `legend()` usa `par(xpd = TRUE)` sobre el panel externo del margen.
#'
#' @param corner        String: "topleft" | "topright" | "bottomleft" | "bottomright"
#' @param mar_base      Vector numérico c(b, l, t, r) de márgenes base en líneas.
#' @param leyenda_h     Alto estimado de la leyenda en líneas (default 4).
#' @param leyenda_w     Ancho estimado de la leyenda en líneas (default 6).
#' @param n_chars       Número de caracteres en la leyenda (amplía `leyenda_h`
#'                      proporcionalmente cuando hay varios bloques). Default 1.
#' @return Lista con `mar` (nuevo vector de márgenes) y `corner` (string
#'         tal como espera `legend()`).
.emtree_config_leyenda <- function(corner, mar_base, leyenda_h = 4,
                                   leyenda_w = 6, n_chars = 1L) {
  # Mantenido por compatibilidad con llamadas en main_MultiMapR.R (mapeo_simple).
  # En export_multimapr_tree() la reserva de espacio se hace con
  # .emtree_medir_leyenda() tras el primer render, expandiendo xlim/ylim.
  esquinas_ok <- c("topleft", "topright", "bottomleft", "bottomright")
  if (!corner %in% esquinas_ok)
    stop(sprintf("`legend_corner` debe ser una de: %s.",
                 paste(esquinas_ok, collapse = ", ")))

  leyenda_h_total <- leyenda_h * max(1L, n_chars) + max(0L, n_chars - 1L)

  mar <- mar_base
  if (grepl("top",    corner)) mar[3L] <- max(mar[3L], leyenda_h_total)
  if (grepl("bottom", corner)) mar[1L] <- max(mar[1L], leyenda_h_total)
  if (grepl("left",   corner)) mar[2L] <- max(mar[2L], leyenda_w)
  if (grepl("right",  corner)) mar[4L] <- max(mar[4L], leyenda_w)

  list(mar = mar, corner = corner)
}

#' Mide el espacio que ocupa la leyenda en coordenadas de datos usando
#' `legend(..., plot = FALSE)`, que calcula el rect sin renderizar nada.
#'
#' Debe llamarse **después** de que el árbol ya esté dibujado y par("usr") /
#' par("pin") reflejen el lienzo real.  Devuelve cuántas unidades de datos
#' adicionales hay que reservar en el eje X (izquierda o derecha) y en el eje
#' Y (arriba o abajo) para que la leyenda no se superponga al gráfico.
#'
#' @param corner         String de esquina.
#' @param legend_by_char Lista agrupada carácter → list(labels, colors), o NULL.
#' @param legend_labels  Vector de etiquetas (modo legacy).
#' @param legend_title   Título del bloque (modo legacy).
#' @param cex_ley        Tamaño de fuente de la leyenda.
#' @return Lista: `dx` (unidades extra en X), `dy` (unidades extra en Y),
#'         `going_down` (bool), `on_right` (bool).
.emtree_medir_leyenda <- function(corner,
                                  legend_by_char = NULL,
                                  legend_labels  = NULL,
                                  legend_title   = NULL,
                                  cex_ley        = 0.8) {

  going_down <- grepl("top",   corner)
  on_right   <- grepl("right", corner)
  usr        <- par("usr")   # c(x1, x2, y1, y2) en coordenadas de datos

  # ── Función interna: mide un bloque de leyenda con plot=FALSE ─────────────
  medir_bloque <- function(labels, title_txt, x_ref, y_ref) {
    lg <- legend(x      = x_ref,
                 y      = y_ref,
                 legend = labels,
                 pch    = 15,
                 col    = rep("black", length(labels)),
                 title  = title_txt,
                 bty    = "n",
                 cex    = cex_ley,
                 horiz  = FALSE,
                 pt.cex = cex_ley * 1.2,
                 xjust  = if (on_right) 1 else 0,
                 yjust  = if (going_down) 1 else 0,
                 plot   = FALSE)   # <<< mide sin dibujar
    lg$rect
  }

  # Posición de referencia: la esquina del área de datos
  x_ref <- if (on_right)   usr[2L] else usr[1L]
  y_ref <- if (going_down) usr[4L] else usr[3L]

  total_w <- 0
  total_h <- 0
  line_h  <- strheight("M", cex = cex_ley) * 1.4

  if (!is.null(legend_by_char) && length(legend_by_char) > 0L) {
    y_cursor <- y_ref
    for (nm in names(legend_by_char)) {
      blk  <- legend_by_char[[nm]]
      rect <- medir_bloque(blk$labels, nm, x_ref, y_cursor)
      total_w  <- max(total_w, rect$w)
      total_h  <- total_h + rect$h + line_h
      y_cursor <- if (going_down) y_cursor - rect$h - line_h
      else            y_cursor + rect$h + line_h
    }
  } else if (!is.null(legend_labels) && length(legend_labels) > 0L) {
    rect    <- medir_bloque(legend_labels, legend_title, x_ref, y_ref)
    total_w <- rect$w
    total_h <- rect$h
  } else {
    return(list(dx = 0, dy = 0, going_down = going_down, on_right = on_right))
  }

  # Margen de seguridad: 8 % extra para que la leyenda no toque el borde
  list(dx         = total_w * 1.08,
       dy         = total_h * 1.08,
       going_down = going_down,
       on_right   = on_right)
}

#' Dibuja la leyenda en la esquina reservada del margen
#'
#' Soporta dos modos:
#'   - **Agrupado** (`legend_by_char` no NULL): recibe una lista nombrada
#'     carácter → list(labels, colors) y dibuja cada carácter como un bloque
#'     independiente con su encabezado en negrita, sin mezclar estados.
#'   - **Plano** (legacy): usa `legend_labels` + `legend_colors` + `legend_title`
#'     como antes, para mantener compatibilidad con llamadas antiguas.
#'
#' Se llama **después** de renderizar el árbol. Usa `par(xpd = NA)` para
#' poder pintar fuera del área de graficado pero dentro del dispositivo.
#'
#' @param corner         String de esquina ("topleft", etc.)
#' @param legend_by_char Lista nombrada carácter → list(labels, colors). Si no
#'                       es NULL tiene precedencia sobre `legend_labels`.
#' @param legend_labels  Vector plano de etiquetas (modo legacy).
#' @param legend_colors  Vector plano de colores (modo legacy).
#' @param legend_title   Título del bloque único (modo legacy, string o NULL).
#' @param pch            Símbolo de leyenda (default 15 = cuadrado relleno).
#' @param cex_ley        Tamaño de fuente para la leyenda.
.emtree_dibujar_leyenda <- function(corner,
                                    legend_by_char = NULL,
                                    legend_labels  = NULL,
                                    legend_colors  = NULL,
                                    legend_title   = NULL,
                                    pch            = 15,
                                    cex_ley        = 0.8) {

  old_xpd <- par("xpd")
  par(xpd = NA)
  on.exit(par(xpd = old_xpd))

  # ── Modo agrupado: un bloque por carácter con encabezado ──────────────────
  if (!is.null(legend_by_char) && length(legend_by_char) > 0L) {

    # Determinar la posición inicial según la esquina elegida.
    # Usamos `legend()` de forma acumulativa: cada bloque devuelve su `rect`
    # y el siguiente se posiciona justo debajo (o arriba, según la esquina).
    usr <- par("usr")   # c(x1, x2, y1, y2)

    # Separación vertical entre bloques (en unidades de datos)
    line_h <- strheight("M", cex = cex_ley) * 1.4

    # Posición de arranque según la esquina
    going_down <- grepl("top", corner)
    x_start <- if (grepl("left",  corner)) usr[1L] else usr[2L]
    y_start <- if (going_down)             usr[4L] else usr[3L]

    y_cursor <- y_start

    for (nm in names(legend_by_char)) {
      bloque <- legend_by_char[[nm]]
      lbl    <- bloque$labels
      col    <- bloque$colors

      if (length(lbl) == 0L) next

      # Encabezado del carácter (dibujado como título del bloque)
      lg <- legend(
        x      = x_start,
        y      = y_cursor,
        legend = lbl,
        col    = col,
        pch    = pch,
        title  = nm,          # nombre del carácter como encabezado
        bty    = "n",
        cex    = cex_ley,
        horiz  = FALSE,
        pt.cex = cex_ley * 1.2,
        xjust  = if (grepl("right", corner)) 1 else 0,
        yjust  = if (going_down) 1 else 0
      )

      # Avanzar el cursor: alto del rect + separación extra entre bloques
      bloque_h <- lg$rect$h
      if (going_down) {
        y_cursor <- y_cursor - bloque_h - line_h
      } else {
        y_cursor <- y_cursor + bloque_h + line_h
      }
    }

    return(invisible(NULL))
  }

  # ── Modo legacy (plano): un único bloque ──────────────────────────────────
  if (is.null(legend_labels) || length(legend_labels) == 0L)
    return(invisible(NULL))

  legend(corner,
         legend = legend_labels,
         col    = legend_colors,
         pch    = pch,
         title  = legend_title,
         bty    = "n",
         cex    = cex_ley,
         horiz  = FALSE,
         pt.cex = cex_ley * 1.2)

  invisible(NULL)
}


# ==============================================================================
# RENDERIZADORES POR TOPOLOGÍA
# ==============================================================================

# ------------------------------------------------------------------------------
# PHYLOGRAM
# Horizontales vectorizadas; verticales en loop protegido contra degeneración.
# ------------------------------------------------------------------------------

.emtree_render_phylogram <- function(pp, tree, color_list, lwd, desfases) {
  xx <- pp$xx
  yy <- pp$yy

  edges     <- tree$edge
  padre_idx <- edges[, 1L]
  hijo_idx  <- edges[, 2L]
  n_tips    <- Ntip(tree)

  # Nodos internos y sus metadatos
  nodos_internos <- unique(edges[edges[, 1L] > n_tips, 1L])
  nodo_x         <- xx[nodos_internos]
  idx_entrada    <- match(nodos_internos, edges[, 2L])   # NA = raíz
  hijos_de       <- lapply(nodos_internos,
                           function(nd) edges[edges[, 1L] == nd, 2L])

  N <- length(color_list)

  for (i in seq_len(N)) {
    ec   <- color_list[[i]]
    lwd_i <- lwd[min(i, length(lwd))]
    dy   <- desfases$dy[min(i, length(desfases$dy))]
    dx   <- desfases$dx[min(i, length(desfases$dx))]

    # ── Horizontales (vectorizado, sin loop) ──────────────────────────────────
    # Cada arista va desde x_padre hasta x_hijo, a la altura y_hijo (convención ape).
    segments(x0   = xx[padre_idx] + dx,
             y0   = yy[hijo_idx]  + dy,
             x1   = xx[hijo_idx]  + dx,
             y1   = yy[hijo_idx]  + dy,
             col  = ec,
             lwd  = lwd_i,
             lend = 1L)

    # ── Verticales (loop sobre nodos internos) ────────────────────────────────
    # Une los hijos de un nodo en la X del nodo.
    # Color: arista entrante al nodo; si es raíz, primera arista saliente.
    for (m in seq_along(nodos_internos)) {
      ie         <- idx_entrada[m]
      color_nodo <- if (!is.na(ie)) {
        ec[ie]
      } else {
        ec[which(edges[, 1L] == nodos_internos[m])[1L]]
      }

      hijos_m <- hijos_de[[m]]
      y_hijos <- yy[hijos_m]
      y_min   <- min(y_hijos)
      y_max   <- max(y_hijos)

      if (y_min != y_max) {   # Omitir segmentos degenerados (politomías / long. 0)
        segments(x0   = nodo_x[m] + dx,
                 y0   = y_min     + dy,
                 x1   = nodo_x[m] + dx,
                 y1   = y_max     + dy,
                 col  = color_nodo,
                 lwd  = lwd_i,
                 lend = 1L)
      }
    }
  }
}


# ------------------------------------------------------------------------------
# CLADOGRAM
# Geometría directa: diagonales padre → hijo.  Desfase en ambos extremos.
# ------------------------------------------------------------------------------

.emtree_render_cladogram <- function(pp, tree, color_list, lwd, desfases) {
  xx <- pp$xx
  yy <- pp$yy

  edges     <- tree$edge
  padre_idx <- edges[, 1L]
  hijo_idx  <- edges[, 2L]

  N <- length(color_list)

  for (i in seq_len(N)) {
    ec    <- color_list[[i]]
    lwd_i <- lwd[min(i, length(lwd))]
    dy    <- desfases$dy[min(i, length(desfases$dy))]
    dx    <- desfases$dx[min(i, length(desfases$dx))]

    # Segmentos diagonales padre→hijo (completamente vectorizado).
    # dx / dy se suma a AMBOS extremos para desplazar el segmento entero.
    segments(x0   = xx[padre_idx] + dx,
             y0   = yy[padre_idx] + dy,
             x1   = xx[hijo_idx]  + dx,
             y1   = yy[hijo_idx]  + dy,
             col  = ec,
             lwd  = lwd_i,
             lend = 1L)
  }
}


# ------------------------------------------------------------------------------
# FAN
# Geometría polar concéntrica con separación compacta calibrada al árbol.
#
# Escala del spread:
#   Las capas deben verse casi paralelas, como en un multimapeo superpuesto
#   ajustado, no como anillos separados.  Para lograrlo el spread total se
#   expresa como un porcentaje pequeño de R_max (radio máximo del árbol):
#
#     spread  = R_max * 0.03   → 3% del radio total, repartido entre N capas
#     vect_D  = seq(-spread/2, spread/2, length.out = N)
#
#   Para cada historia i:
#     dr       = vect_D[i]            desplazamiento radial
#     dtheta_n = dtheta_nodo[n] * dr/R_max  rotación angular por nodo, de modo
#                                     que cada extremo del trazo use su propio
#                                     ángulo → trazos paralelos al original
#
# Corrección de cruce del eje Este (0° / 360°):
#   Normalización a [0, 2π) con ifelse(); si ang_max − ang_min > π el clado
#   cruza 0° y se elevan los ángulos < π al rango extendido [2π, 4π).
# ------------------------------------------------------------------------------

.emtree_render_fan <- function(pp, tree, color_list, lwd, desfases) {
  xx     <- pp$xx
  yy     <- pp$yy
  edges  <- tree$edge
  n_tips <- Ntip(tree)
  N      <- length(color_list)

  # Extraer ángulos originales de todos los nodos
  angulos_nodos <- atan2(yy, xx)

  # ── Metadatos de nodos internos ───────────────────────────────────────────
  nodos_internos <- unique(edges[edges[, 1L] > n_tips, 1L])
  idx_entrada    <- match(nodos_internos, edges[, 2L])
  hijos_de       <- lapply(nodos_internos, function(nd) edges[edges[, 1L] == nd, 2L])

  # ── Spread ────────────────────────────────────────────────────────────────
  R_max  <- max(sqrt(xx^2 + yy^2))
  # El spread entre capas en modo fan escala con el grosor (lwd) de las ramas
  # para que capas más anchas no se superpongan. Factor empírico: 0.0025 * lwd,
  # normalizado por R_max para ser independiente del tamaño del árbol.
  # El vector lwd aquí tiene un elemento por historia (puede variar); tomamos
  # el máximo como referencia conservadora.
  lwd_ref <- max(lwd)
  spread  <- R_max * max(0.01, lwd_ref * 0.0025)
  vect_D  <- if (N == 1L) 0 else {
    indices <- seq_len(N) - (N + 1L) / 2
    indices * (spread / max(1L, N - 1L))
  }

  # Función auxiliar: compensación angular exacta para paralelismo
  # Δθ = asin(D / (R + D))  — el clamp evita NaN cuando R+D ≈ 0
  calc_dtheta <- function(R, D) {
    R_new <- R + D
    if (abs(R_new) < 1e-8) return(pi / 2)
    ratio <- max(-1, min(1, D / R_new))
    asin(ratio)
  }

  # ── Bucle de historias ────────────────────────────────────────────────────
  for (i in seq_len(N)) {
    ec    <- color_list[[i]]
    lwd_i <- lwd[min(i, length(lwd))]
    dr    <- vect_D[i]

    # 1. Trazos radiales
    for (j in seq_len(nrow(edges))) {
      p <- edges[j, 1L];  h <- edges[j, 2L]

      R_p <- sqrt(xx[p]^2 + yy[p]^2)
      R_h <- sqrt(xx[h]^2 + yy[h]^2)

      R_p_nuevo <- R_p + dr
      R_h_nuevo <- R_h + dr

      # Compensación angular exacta en cada extremo según su radio local
      dtheta_p <- calc_dtheta(R_p, dr)
      dtheta_h <- calc_dtheta(R_h, dr)

      # Ambas puntas usan el ángulo del hijo; R maneja radios negativos nativamente
      theta_p_nuevo <- angulos_nodos[h] + dtheta_p
      theta_h_nuevo <- angulos_nodos[h] + dtheta_h

      segments(x0  = R_p_nuevo * cos(theta_p_nuevo),
               y0  = R_p_nuevo * sin(theta_p_nuevo),
               x1  = R_h_nuevo * cos(theta_h_nuevo),
               y1  = R_h_nuevo * sin(theta_h_nuevo),
               col = ec[j], lwd = lwd_i, lend = 1L)
    }

    # 2. Arcos
    for (m in seq_along(nodos_internos)) {
      nd <- nodos_internos[m]
      ie <- idx_entrada[m]
      color_nodo <- if (!is.na(ie)) ec[ie] else ec[which(edges[, 1L] == nd)[1L]]

      hijos_m <- hijos_de[[m]]
      R_nd    <- sqrt(xx[nd]^2 + yy[nd]^2)
      R_nuevo <- R_nd + dr

      dtheta_nd <- calc_dtheta(R_nd, dr)

      theta_adj_list <- numeric(length(hijos_m))
      for (k in seq_along(hijos_m)) {
        h <- hijos_m[k]
        theta_adj_list[k] <- angulos_nodos[h] + dtheta_nd
      }

      ang_norm <- ifelse(theta_adj_list < 0, theta_adj_list + 2 * pi, theta_adj_list)
      ang_min  <- min(ang_norm)
      ang_max  <- max(ang_norm)

      if ((ang_max - ang_min) > pi) {
        ang_norm <- ifelse(ang_norm < pi, ang_norm + 2 * pi, ang_norm)
        ang_min  <- min(ang_norm)
        ang_max  <- max(ang_norm)
      }

      ang_seq <- seq(ang_min, ang_max, length.out = 100L)

      lines(x   = R_nuevo * cos(ang_seq),
            y   = R_nuevo * sin(ang_seq),
            col = color_nodo, lwd = lwd_i, lend = 1L)
    }
  }
}


# ==============================================================================
# FUNCIÓN PRINCIPAL
# ==============================================================================

#' Exporta un árbol filogenético multi-historia evolutiva
#'
#' Superpone N historias evolutivas (vectores de colores de ramas) sobre una
#' filogenia soportando tres tipos de topología y dos formatos de salida.
#' Aplica la estrategia de "Canvas en Blanco" con desfase diagonal isótropo
#' calculado al vuelo a partir de los parámetros físicos reales del dispositivo.
#'
#' @section Flujo interno:
#' \enumerate{
#'   \item Validación exhaustiva de argumentos.
#'   \item Cálculo dinámico de dimensiones (alto ∝ Ntip; fan = cuadrado).
#'   \item Apertura del dispositivo gráfico (PNG @ 300 dpi, o PDF).
#'   \item Canvas en Blanco: \code{plot(tree, edge.color="transparent")}
#'         puebla \code{.PlotPhyloEnv} sin renderizar aristas.
#'   \item Extracción de coordenadas nodales y parámetros físicos del lienzo.
#'   \item Cálculo al vuelo de desfases isótropos (escala_x/escala_y).
#'   \item Renderizado por topología (phylogram / cladogram / fan).
#'   \item Cierre garantizado del dispositivo vía \code{finally}.
#' }
#'
#' @param tree          Objeto \code{\link[ape]{read.tree}} de clase \code{phylo}.
#' @param color_list    Lista de N vectores de caracteres (colores R/CSS válidos).
#'                      Cada vector debe tener exactamente \code{nrow(tree$edge)}
#'                      elementos (uno por arista).
#' @param filename      Nombre del archivo de salida \strong{sin extensión}.
#'                      La extensión correcta se añade automáticamente.
#' @param type          Topología del árbol: \code{"phylogram"} (default),
#'                      \code{"cladogram"} o \code{"fan"}.
#' @param format        Formato de salida: \code{"png"} (default) o \code{"pdf"}.
#' @param lwd           Grosor de las ramas (default: \code{3}).
#' @param label_offset  Distancia entre el extremo de los terminales y sus
#'                      etiquetas, en unidades de coordenadas de datos
#'                      (default: \code{0.3}).
#' @param rango_desfase Amplitud total del rango de desfase en Y.  El vector
#'                      \code{desfases_y} se construye como
#'                      \code{seq(-rango, rango, length.out = N)}.
#'                      Con \code{N = 1} el desfase es siempre cero (default: \code{0.1}).
#' @param mar           Márgenes del gráfico como vector \code{c(b, l, t, r)}
#'                      en líneas (default: \code{c(1, 1, 1, 4)}).
#'
#' @return Invisible: ruta completa del archivo generado (string).
#'
#' @examples
#' \dontrun{
#' library(ape)
#' set.seed(42)
#' tree <- rtree(30)
#' n_e  <- nrow(tree$edge)
#'
#' ec1 <- sample(c("steelblue",  "tomato",        "gray70"), n_e, replace = TRUE)
#' ec2 <- sample(c("darkorange", "mediumseagreen", "gray70"), n_e, replace = TRUE)
#' ec3 <- sample(c("mediumpurple","gold3",         "gray70"), n_e, replace = TRUE)
#'
#' # Phylogram PNG
#' export_multimapr_tree(tree, list(ec1, ec2, ec3),
#'                       "mi_arbol_phylogram", type = "phylogram", format = "png")
#'
#' # Cladogram PDF
#' export_multimapr_tree(tree, list(ec1, ec2),
#'                       "mi_arbol_cladogram", type = "cladogram", format = "pdf")
#'
#' # Fan PNG
#' export_multimapr_tree(tree, list(ec1, ec2, ec3),
#'                       "mi_arbol_fan",       type = "fan",       format = "png")
#' }
#'
#' @export
export_multimapr_tree <- function(tree,
                                  color_list,
                                  filename,
                                  type          = "phylogram",
                                  format        = "png",
                                  lwd           = 3,
                                  label_offset  = 0.3,
                                  rango_desfase = 0.1,
                                  mar           = c(1, 1, 1, 4),
                                  # ── PARÁMETROS DE DIMENSIONES PERSONALIZADAS ──────────
                                  # NULL → se usan las dimensiones automáticas calculadas
                                  # a partir del número de terminales (comportamiento original).
                                  # Un valor numérico en pulgadas activa el escalado dinámico.
                                  width         = NULL,
                                  height        = NULL,
                                  # ── PARÁMETROS DE LEYENDA ─────────────────────────────
                                  # legend_by_char: lista nombrada carácter → list(labels, colors).
                                  #   Cuando se proporciona, cada carácter se dibuja como un bloque
                                  #   independiente con su nombre como encabezado. Tiene precedencia
                                  #   sobre legend_labels / legend_colors / legend_title.
                                  # legend_labels / legend_colors: vectores paralelos con las
                                  #   etiquetas y colores de la leyenda (modo legacy).  NULL = sin leyenda.
                                  # legend_corner : esquina donde se ubica la leyenda;
                                  #   una de "topleft" | "topright" | "bottomleft" | "bottomright".
                                  # legend_title  : título del bloque de leyenda (string o NULL, solo legacy).
                                  legend_by_char = NULL,
                                  legend_labels  = NULL,
                                  legend_colors  = NULL,
                                  legend_corner  = "bottomleft",
                                  legend_title   = NULL) {

  # ── 0. Validación de argumentos ─────────────────────────────────────────────
  .emtree_validar_arbol(tree)
  .emtree_validar_color_list(color_list, nrow(tree$edge))
  .emtree_validar_filename(filename)
  .emtree_validar_tipo_formato(type, format)

  if (!is.numeric(lwd)          || length(lwd) != 1L || lwd <= 0)
    stop("`lwd` debe ser un número positivo.")
  if (!is.numeric(label_offset) || length(label_offset) != 1L)
    stop("`label_offset` debe ser un escalar numérico.")
  if (!is.numeric(rango_desfase) || length(rango_desfase) != 1L || rango_desfase < 0)
    stop("`rango_desfase` debe ser un número no negativo.")
  if (!is.numeric(mar) || length(mar) != 4L)
    stop("`mar` debe ser un vector numérico de longitud 4.")
  if (!is.null(width)  && (!is.numeric(width)  || length(width)  != 1L || width  <= 0))
    stop("`width` debe ser un número positivo en pulgadas, o NULL.")
  if (!is.null(height) && (!is.numeric(height) || length(height) != 1L || height <= 0))
    stop("`height` debe ser un número positivo en pulgadas, o NULL.")

  # Validar parámetros de leyenda
  # Prioridad: legend_by_char > legend_labels/legend_colors (modo legacy)
  tiene_leyenda_agrupada <- !is.null(legend_by_char) && length(legend_by_char) > 0L
  tiene_leyenda_plana    <- !tiene_leyenda_agrupada &&
    !is.null(legend_labels) && length(legend_labels) > 0L
  tiene_leyenda          <- tiene_leyenda_agrupada || tiene_leyenda_plana

  if (tiene_leyenda_agrupada) {
    # Solo validar estructura; el espacio se reserva con xlim/ylim expandidos
    # tras medir la leyenda con .emtree_medir_leyenda() (Fase 1 → Fase 2).
    for (nm in names(legend_by_char)) {
      blk <- legend_by_char[[nm]]
      if (!is.list(blk) || is.null(blk$labels) || is.null(blk$colors))
        stop(sprintf(
          "`legend_by_char[[\"%s\"]]` debe ser una lista con elementos `labels` y `colors`.", nm))
      if (length(blk$labels) != length(blk$colors))
        stop(sprintf(
          "`legend_by_char[[\"%s\"]]`: `labels` y `colors` deben tener la misma longitud.", nm))
    }

  } else if (tiene_leyenda_plana) {
    if (is.null(legend_colors) || length(legend_colors) != length(legend_labels))
      stop("`legend_colors` debe tener la misma longitud que `legend_labels`.")
  }

  # ── 1. Cálculo de dimensiones y apertura del dispositivo ────────────────────
  n_tips <- Ntip(tree)
  N      <- length(color_list)

  # Crear carpeta Exports/ si no existe y redirigir la salida allí
  dir.create("Exports", showWarnings = FALSE, recursive = TRUE)
  filename <- file.path("Exports", basename(filename))

  # Dimensiones por defecto (comportamiento original)
  alto_default  <- n_tips * 0.25 + 2   # 0.25 pulgadas por terminal + margen base
  ancho_default <- 12

  # Fan → lienzo cuadrado para evitar deformaciones polares + 1.5 in extra
  # para dar mayor respiración a las etiquetas y la leyenda en modo radial.
  if (type == "fan") {
    alto_default  <- max(alto_default, ancho_default)
    ancho_default <- alto_default
  }

  # Dimensiones finales: usuario > automático
  alto_in  <- if (!is.null(height)) height else alto_default
  ancho_in <- if (!is.null(width))  width  else ancho_default

  # ── FACTORES DE ESCALA ───────────────────────────────────────────────────────
  # Miden cuánto crece (o encoge) cada eje respecto al tamaño de referencia.
  # escala_h > 1 → lienzo más alto que el default; escala_h < 1 → más compacto.
  # Se usan para: (a) ajustar cex y label_offset proporcionalmente,
  #               (b) compensar la dilatación del desfase en Y (escala inversa).
  escala_h <- alto_in  / alto_default
  escala_w <- ancho_in / ancho_default

  # ── 2. Parámetros de render (cex, offsets, desfase) ─────────────────────────
  # Se calculan ANTES de abrir el dispositivo final para poder usarlos también
  # en el dispositivo off-screen de sondeo (Fase 1).
  cex_base <- max(1 / (1 + n_tips / 50), 0.2)
  cex_aj   <- cex_base * escala_h

  # Offset proporcional a la escala geométrica real de los datos (2.5%)
  phy_tmp <- tree
  if (is.null(phy_tmp$edge.length)) phy_tmp$edge.length <- rep(1, nrow(phy_tmp$edge))
  max_depth <- max(node.depth.edgelength(phy_tmp))

  label_offset_aj <- max_depth * 0.025 * escala_w

  if (rango_desfase == 0.1 && N > 1L) {
    factor_lwd <- switch(type,
                         "phylogram" = 0.012, "cladogram" = 0.012, "fan" = 0.003, 0.012)
    rango_desfase <- max(0.05, lwd * factor_lwd)
  }
  rango_desfase_ajustado <- rango_desfase / escala_h
  lwd_vec <- lwd * (0.95 ^ (seq_len(N) - 1L))
  cex_ley <- max(0.5, min(1.0, 10 / (n_tips + 10)))

  mar_base <- c(1, 1, 1, 4)

  # ── Función interna de render del canvas ─────────────────────────────────────
  # Dibuja el árbol "invisible" (edge.color = "transparent") para poblar
  # .PlotPhyloEnv y establecer par("usr"). Se reutiliza en el sondeo off-screen
  # y en el render definitivo dentro del dispositivo final.
  .render_canvas <- function(xlim_extra = NULL, ylim_extra = NULL) {
    if (type == "fan") par(mar = rep(2L, 4L)) else par(mar = mar_base)
    plot(tree,
         type           = type,
         edge.color     = "transparent",
         tip.color      = if (type == "fan") "transparent" else "black",
         edge.width     = lwd,
         cex            = cex_aj,
         label.offset   = label_offset_aj,
         no.margin      = if (type == "fan") FALSE else TRUE,
         show.tip.label = TRUE,
         font           = 3L,    # <--- Itálica activada
         x.lim          = xlim_extra,
         y.lim          = ylim_extra)
    .PlotPhyloEnv <- .emtree_get_PlotPhyloEnv()
    get("last_plot.phylo", envir = .PlotPhyloEnv)
  }

  # ── 3. FASE 1 — sondeo en dispositivo off-screen ─────────────────────────────
  # Se abre un pdf(NULL) temporal (sin archivo en disco) para:
  #   a) establecer par("usr") con las mismas dimensiones que el dispositivo final
  #   b) medir la leyenda con legend(..., plot = FALSE) en coordenadas de datos
  # Al cerrar con dev.off() el dispositivo temporal desaparece sin dejar archivo.
  xlim_final <- NULL
  ylim_final <- NULL

  if (tiene_leyenda) {
    pdf(NULL, width = ancho_in, height = alto_in)   # off-screen, sin archivo
    tryCatch({
      .render_canvas()   # solo para establecer par("usr")

      med <- .emtree_medir_leyenda(
        corner         = legend_corner,
        legend_by_char = if (tiene_leyenda_agrupada) legend_by_char else NULL,
        legend_labels  = if (!tiene_leyenda_agrupada) legend_labels else NULL,
        legend_title   = if (!tiene_leyenda_agrupada) legend_title  else NULL,
        cex_ley        = cex_ley
      )

      usr <- par("usr")

      if (type == "fan") {
        expansion <- max(med$dx, med$dy)
        xlim_final <- if (med$on_right) c(usr[1L], usr[2L] + expansion)
        else              c(usr[1L] - expansion, usr[2L])
        ylim_final <- c(usr[3L], usr[4L])
        if (med$going_down) ylim_final[2L] <- ylim_final[2L] + med$dy * 0.5
        else                ylim_final[1L] <- ylim_final[1L] - med$dy * 0.5
      } else {
        xlim_final <- if (med$on_right) c(usr[1L], usr[2L] + med$dx)
        else              c(usr[1L] - med$dx, usr[2L])
        ylim_final <- if (med$going_down) c(usr[3L], usr[4L] + med$dy)
        else               c(usr[3L] - med$dy, usr[4L])
      }
    }, finally = dev.off())   # cerrar siempre; no deja archivo en disco
  }

  # ── 4. Abrir dispositivo final y renderizar ───────────────────────────────────
  .emtree_abrir_dispositivo(filename, format, ancho_in, alto_in)
  ruta_salida <- paste0(filename, ".", format)

  tryCatch({

    # Canvas definitivo: con xlim/ylim expandidos si hay leyenda, o normal si no.
    pp      <- .render_canvas(xlim_extra = xlim_final, ylim_extra = ylim_final)
    desfases <- .emtree_calcular_desfases(N, rango_desfase_ajustado)

    # ── 5. Renderizado de ramas según topología ───────────────────────────────
    if (type == "phylogram") {
      .emtree_render_phylogram(pp, tree, color_list, lwd_vec, desfases)

    } else if (type == "cladogram") {
      .emtree_render_cladogram(pp, tree, color_list, lwd_vec, desfases)

    } else {
      .emtree_render_fan(pp, tree, color_list, lwd_vec, desfases)

      # ── 5b. Etiquetas de terminales en fan — radio exacto ─────────────────────
      xx_tips <- pp$xx[seq_len(n_tips)]
      yy_tips <- pp$yy[seq_len(n_tips)]
      angulos <- atan2(yy_tips, xx_tips)
      R_tips  <- sqrt(xx_tips^2 + yy_tips^2)
      gap_u   <- strwidth("m", cex = cex_aj) * 1.5

      old_xpd <- par("xpd")
      par(xpd = NA)
      for (j in seq_len(n_tips)) {
        ang_j  <- angulos[j]
        R_j    <- R_tips[j] + gap_u
        lado_d <- cos(ang_j) >= 0
        srt_j  <- ang_j * 180 / pi
        adj_j  <- if (lado_d) c(0, 0.5) else { srt_j <- srt_j + 180; c(1, 0.5) }
        text(R_j * cos(ang_j),
             R_j * sin(ang_j),
             labels = tree$tip.label[j],
             adj    = adj_j,
             cex    = cex_aj,
             font   = 3L,        # <--- Itálica activada
             srt    = srt_j)
      }
      par(xpd = old_xpd)
    }

    # ── 6. Leyenda — dibujada dentro del viewport expandido ───────────────────
    if (tiene_leyenda) {
      .emtree_dibujar_leyenda(
        corner         = legend_corner,
        legend_by_char = if (tiene_leyenda_agrupada) legend_by_char else NULL,
        legend_labels  = legend_labels,
        legend_colors  = legend_colors,
        legend_title   = legend_title,
        cex_ley        = cex_ley
      )
    }

    # Indicar si las dimensiones son automáticas o definidas por el usuario
    dims_origen <- if (!is.null(width) || !is.null(height)) " [custom]" else " [auto]"
    message(sprintf(
      "[export_multimapr_tree] Archivo generado: %s  (%d tips, %d historias, %.1f × %.1f in%s%s)",
      ruta_salida, n_tips, N, ancho_in, alto_in,
      if (format == "png") " @ 300 dpi" else "",
      dims_origen))

  }, error = function(e) {
    stop(sprintf("[export_multimapr_tree] Error durante el renderizado (%s, %s): %s",
                 type, format, conditionMessage(e)))
  }, finally = {
    dev.off()   # Cerrar el dispositivo en todos los casos
  })

  invisible(ruta_salida)
}


# ==============================================================================
# EJEMPLOS DE USO  (cambiar FALSE → TRUE para ejecutar)
# ==============================================================================
if (FALSE) {

  library(ape)
  set.seed(42)

  tree <- rtree(40)
  n_e  <- nrow(tree$edge)

  paleta1 <- c("steelblue",   "tomato",         "gray70")
  paleta2 <- c("darkorange",  "mediumseagreen",  "gray70")
  paleta3 <- c("mediumpurple","gold3",           "gray70")

  ec1 <- sample(paleta1, n_e, replace = TRUE)
  ec2 <- sample(paleta2, n_e, replace = TRUE)
  ec3 <- sample(paleta3, n_e, replace = TRUE)

  # ── Phylogram con leyenda en esquina inferior-izquierda ────────────────────
  export_multimapr_tree(
    tree          = tree,
    color_list    = list(ec1, ec2, ec3),
    filename      = "test_phylogram",
    type          = "phylogram",
    format        = "png",
    legend_labels = c("Estado A", "Estado B", "Ambiguo"),
    legend_colors = paleta1,
    legend_corner = "bottomleft",
    legend_title  = "Carácter 1"
  )

  # ── Cladogram con leyenda en esquina superior-derecha ──────────────────────
  export_multimapr_tree(
    tree          = tree,
    color_list    = list(ec1, ec2),
    filename      = "test_cladogram",
    type          = "cladogram",
    format        = "pdf",
    legend_labels = c("Estado A", "Estado B", "Ambiguo"),
    legend_colors = paleta1,
    legend_corner = "topright",
    legend_title  = "Carácter 1"
  )

  # ── Fan con leyenda en esquina superior-izquierda ──────────────────────────
  export_multimapr_tree(
    tree          = tree,
    color_list    = list(ec1, ec2, ec3),
    filename      = "test_fan",
    type          = "fan",
    format        = "png",
    legend_labels = c("Estado A", "Estado B", "Ambiguo"),
    legend_colors = paleta1,
    legend_corner = "topleft",
    legend_title  = "Carácter 1"
  )

  # ── Una sola historia sin leyenda (comportamiento original) ────────────────
  export_multimapr_tree(
    tree       = tree,
    color_list = list(ec1),
    filename   = "test_single",
    type       = "phylogram",
    format     = "pdf"
  )
}

# ==============================================================================
# RENDERIZADO EN PANTALLA — MOTOR AVANZADO
# ==============================================================================

#' Renderiza el árbol multi-historia directamente en la ventana interactiva de R
#'
#' Ejecuta exactamente la misma lógica del Canvas en Blanco que usa
#' \code{export_multimapr_tree}, pero opera sobre la ventana gráfica activa
#' de R en lugar de abrir celdas de almacenamiento en disco.
#' Soporta las tres topologías: \code{"phylogram"}, \code{"cladogram"} y
#' \code{"fan"}.
#'
#' @param tree          Objeto \code{\link[ape]{read.tree}} de clase \code{phylo}.
#' @param color_list    Lista de N vectores de caracteres (colores R/CSS válidos).
#'                      Cada vector debe tener exactamente \code{nrow(tree$edge)}
#'                      elementos (uno por arista).
#' @param type          Topología del árbol: \code{"phylogram"} (default),
#'                      \code{"cladogram"} o \code{"fan"}.
#' @param lwd           Grosor de las ramas (default: \code{3}).
#' @param label_offset  Distancia entre el extremo de los terminales y sus
#'                      etiquetas (default: \code{0.3}).
#' @param rango_desfase Amplitud total del rango de desfase en Y (default: \code{0.1}).
#' @param legend_by_char Lista nombrada carácter → list(labels, colors) para
#'                       leyenda agrupada. Si no es NULL tiene precedencia sobre
#'                       \code{legend_labels} / \code{legend_colors}.
#' @param legend_labels  Vector plano de etiquetas (modo legacy).
#' @param legend_colors  Vector plano de colores (modo legacy).
#' @param legend_corner  Esquina de la leyenda: \code{"topleft"}, \code{"topright"},
#'                       \code{"bottomleft"} (default) o \code{"bottomright"}.
#' @param legend_title   Título del bloque único de leyenda (modo legacy).
#' @return Invisible NULL. Dibuja en el dispositivo gráfico activo.
plot_multimapr_screen <- function(tree,
                                  color_list,
                                  type          = "phylogram",
                                  lwd           = 3,
                                  label_offset  = 0.3,
                                  rango_desfase = 0.1,
                                  legend_by_char = NULL,
                                  legend_labels  = NULL,
                                  legend_colors  = NULL,
                                  legend_corner  = "bottomleft",
                                  legend_title   = NULL) {

  .emtree_validar_arbol(tree)
  .emtree_validar_color_list(color_list, nrow(tree$edge))

  n_tips <- Ntip(tree)
  N      <- length(color_list)

  # Asegurar que haya una ventana gráfica abierta
  if (is.null(dev.list())) dev.new()

  # Factores de escala neutros para pantalla interactiva
  escala_h <- 1
  escala_w <- 1

  cex_base <- max(1 / (1 + n_tips / 50), 0.2)
  cex_aj   <- cex_base * escala_h

  # Offset proporcional a la escala real del árbol (2.5%)
  phy_tmp <- tree
  if (is.null(phy_tmp$edge.length)) phy_tmp$edge.length <- rep(1, nrow(phy_tmp$edge))
  max_depth <- max(node.depth.edgelength(phy_tmp))
  label_offset_aj <- max_depth * 0.025 * escala_w

  rango_desfase_ajustado <- rango_desfase / escala_h
  lwd_vec <- lwd * (0.95 ^ (seq_len(N) - 1L))
  cex_ley <- max(0.5, min(1.0, 10 / (n_tips + 10)))

  mar_base <- c(1, 1, 1, 4)

  # Función interna: Canvas en Blanco para poblar .PlotPhyloEnv
  .render_canvas_screen <- function() {
    if (type == "fan") par(mar = rep(2L, 4L)) else par(mar = mar_base)
    plot(tree,
         type           = type,
         edge.color     = "transparent",
         tip.color      = if (type == "fan") "transparent" else "black",
         edge.width     = lwd,
         cex            = cex_aj,
         label.offset   = label_offset_aj,
         no.margin      = if (type == "fan") FALSE else TRUE,
         show.tip.label = TRUE,
         font           = 3L)
    pp <- get("last_plot.phylo", envir = get(".PlotPhyloEnv", envir = asNamespace("ape")))
    return(pp)
  }

  # Limpiar ventana y renderizar canvas invisible
  plot.new()
  pp       <- .render_canvas_screen()
  desfases <- .emtree_calcular_desfases(N, rango_desfase_ajustado)

  # Despachar al renderizador geométrico según topología
  if (type == "phylogram") {
    .emtree_render_phylogram(pp, tree, color_list, lwd_vec, desfases)

  } else if (type == "cladogram") {
    .emtree_render_cladogram(pp, tree, color_list, lwd_vec, desfases)

  } else if (type == "fan") {
    .emtree_render_fan(pp, tree, color_list, lwd_vec, desfases)

    # Etiquetas radiales manuales en itálicas para modo fan
    xx_tips <- pp$xx[seq_len(n_tips)]
    yy_tips <- pp$yy[seq_len(n_tips)]
    angulos <- atan2(yy_tips, xx_tips)
    R_tips  <- sqrt(xx_tips^2 + yy_tips^2)
    gap_u   <- strwidth("m", cex = cex_aj) * 1.5

    old_xpd <- par("xpd"); par(xpd = NA)
    for (j in seq_len(n_tips)) {
      ang_j  <- angulos[j]
      R_j    <- R_tips[j] + gap_u
      lado_d <- cos(ang_j) >= 0
      srt_j  <- ang_j * 180 / pi
      adj_j  <- if (lado_d) c(0, 0.5) else { srt_j <- srt_j + 180; c(1, 0.5) }
      text(R_j * cos(ang_j), R_j * sin(ang_j),
           labels = tree$tip.label[j], adj = adj_j, cex = cex_aj, font = 3L, srt = srt_j)
    }
    par(xpd = old_xpd)
  }

  # Dibujar leyenda si corresponde
  tiene_leyenda_agrupada <- !is.null(legend_by_char) && length(legend_by_char) > 0L
  tiene_leyenda_plana    <- !tiene_leyenda_agrupada && !is.null(legend_labels) && length(legend_labels) > 0L

  if (tiene_leyenda_agrupada || tiene_leyenda_plana) {
    .emtree_dibujar_leyenda(
      corner         = legend_corner,
      legend_by_char = legend_by_char,
      legend_labels  = legend_labels,
      legend_colors  = legend_colors,
      legend_title   = legend_title,
      cex_ley        = cex_ley
    )
  }

  invisible(NULL)
}
