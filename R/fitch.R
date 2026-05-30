################################################################################
# fitch.R
#
# Adaptador del algoritmo de Fitch para su uso como 'algoritmo_externo'
# dentro del paquete MultiMapR.
#
# Cuando MultiMapR pregunte por el algoritmo, selecciona opción 2 (Fitch).
# El modo de optimización se elige interactivamente desde la terminal:
#   1 = ACCTRAN     (cambios acelerados — hacia las hojas)
#   2 = DELTRAN     (cambios retrasados — hacia la raíz)
#   3 = Unambiguous (solo nodos con estado único tras los dos pasos de Fitch)
#
# FIRMA que MultiMapR espera:
#   algoritmo_externo(tree, tip_colors, edge_colors, config) → edge_colors
#
#   - tree        : objeto phylo
#   - tip_colors  : vector nombrado (color por cada tip, en el orden de tip.label)
#   - edge_colors : vector inicializado con "gray70" (longitud == nrow(tree$edge))
#   - config      : lista de configuración de MultiMapR; se usa config$modo_fitch
#                   (string: "acctran", "deltran" o "unambiguous") para
#                   determinar el método de optimización. Valor por defecto: "deltran".
#
# NOTA: tip_colors ya contiene el color asignado por el usuario a cada estado;
#       esta función invierte ese mapeo color→estado para correr Fitch y luego
#       devuelve los colores de ramas resultantes en el mismo formato.
################################################################################

# ============================================================================ #
# PASO DESCENDENTE (pass de hojas → raíz)
# ============================================================================ #
.fitch_paso_descendente <- function(phylo, estados_tips) {
  n_tips       <- Ntip(phylo)
  n_nodes      <- phylo$Nnode
  total_nodes  <- n_tips + n_nodes

  # Estados únicos reales (excluye "?" y "-")
  all_states <- sort(unique(estados_tips[!estados_tips %in% c("?", "-")]))

  conjuntos   <- vector("list", total_nodes)
  operaciones <- character(total_nodes)

  # Inicializar terminales
  for (i in seq_len(n_tips)) {
    st <- estados_tips[i]
    if (is.na(st) || st == "?" || st == "-" || trimws(st) == "") {
      conjuntos[[i]] <- all_states          # desconocido = todos los estados
    } else {
      if (grepl("^[0-9]+$", st) && nchar(st) > 1) {
        # Polimorfismos numéricos sin separador (ej. "01" -> "0", "1")
        conjuntos[[i]] <- unique(strsplit(st, "")[[1]])
      } else if (grepl("[,/&|]", st)) {
        # Polimorfismos con separador explícito (ej. "0/1" o "urbano,bosque")
        conjuntos[[i]] <- unique(trimws(strsplit(st, "[,/&|]")[[1]]))
      } else {
        # Palabras enteras o estados de un solo caracter ("bosque", "granivoro", "0", "A")
        conjuntos[[i]] <- st
      }
    }
    operaciones[i] <- "terminal"
  }

  # Ordenar nodos internos de hojas a raíz (post-order por n° de descendientes)
  nodos_int <- (n_tips + 1):total_nodes
  contar_desc <- function(nd) {
    if (nd <= n_tips) return(1L)
    hijos <- phylo$edge[phylo$edge[, 1] == nd, 2]
    sum(sapply(hijos, contar_desc))
  }
  desc_counts     <- sapply(nodos_int, contar_desc)
  orden_postorder <- nodos_int[order(desc_counts)]

  for (nd in orden_postorder) {
    hijos <- phylo$edge[phylo$edge[, 1] == nd, 2]
    # Fitch generalizado: funciona con politomías (2 o más hijos)
    inter <- Reduce(intersect, lapply(hijos, function(h) conjuntos[[h]]))

    if (length(inter) > 0) {
      conjuntos[[nd]]  <- inter
      operaciones[nd]  <- "inter"
    } else {
      conjuntos[[nd]]  <- Reduce(union, lapply(hijos, function(h) conjuntos[[h]]))
      operaciones[nd]  <- "union"
    }
  }

  list(conjuntos   = conjuntos,
       operaciones = operaciones,
       n_tips      = n_tips,
       n_nodes     = n_nodes)
}


# ============================================================================ #
# PASO ASCENDENTE (pass de raíz → hojas)
#
# Implementa la regla de Swofford & Maddison (1987) para calcular los
# conjuntos MPR de cada nodo:
#   - Si el nodo fue UNIÓN en el paso descendente:
#       MPR(nd) = union(s_down[nd], MPR(padre))
#   - Si el nodo fue INTERSECCIÓN en el paso descendente:
#       MPR(nd) = union(s_down[nd], intersect(MPR(padre), estados_hijos))
#       donde estados_hijos = estados de todos los hijos en s_down
# ============================================================================ #
.fitch_paso_ascendente <- function(phylo, res_desc) {
  s_down      <- res_desc$conjuntos
  operaciones <- res_desc$operaciones
  n_tips      <- res_desc$n_tips

  raiz <- n_tips + 1L
  s_up <- s_down  # Inicializamos MPR como copia del paso descendente

  preorder_int <- function(nd) {
    res   <- nd
    hijos <- phylo$edge[phylo$edge[, 1] == nd, 2]
    for (h in hijos) if (h > n_tips) res <- c(res, preorder_int(h))
    res
  }
  orden_pre <- preorder_int(raiz)

  for (nd in orden_pre[-1]) {
    padre <- phylo$edge[phylo$edge[, 2] == nd, 1][1]

    c_nodo  <- s_down[[nd]]
    c_padre <- s_up[[padre]] # El MPR definitivo del padre

    if (operaciones[nd] == "union") {
      s_up[[nd]] <- union(c_nodo, c_padre)
    } else {
      hijos <- phylo$edge[phylo$edge[, 1] == nd, 2]
      estados_hijos <- unlist(lapply(hijos, function(h) s_down[[h]]))
      validos <- intersect(c_padre, estados_hijos)
      s_up[[nd]] <- union(c_nodo, validos)
    }
  }

  res_desc$conjuntos <- s_up
  res_desc
}


# ============================================================================ #
# OPTIMIZACIÓN (ACCTRAN / DELTRAN / Unambiguous) → un estado único por nodo
#
# modo_actual (derivado de config$modo_fitch dentro de la función):
#   "acctran"     — en ambigüedades, usa conjuntos del paso DESCENDENTE;
#                   si el estado del padre está en el conjunto del hijo,
#                   lo hereda (acelera cambios hacia las hojas).
#   "deltran"     — en ambigüedades, usa conjuntos MPR del paso ASCENDENTE;
#                   si el estado del padre está en el conjunto del hijo,
#                   lo hereda (retrasa cambios hacia la raíz).
#   "unambiguous" — un nodo es inambiguo EXCLUSIVAMENTE cuando ACCTRAN y
#                   DELTRAN le asignan exactamente el mismo estado; el resto
#                   queda como NA. La raíz es inambigua solo si su conjunto
#                   MPR tiene exactamente un elemento.
#
# REGLA DE NEGOCIO ESTRICTA (inambiguo):
#   Un nodo es inambiguo si y solo si acc[nd] == del[nd] (misma cadena, sin NA).
#   No se relaja esta condición: si comparten estado es inambiguo; si difieren
#   — aunque ambos sean válidos — el nodo queda como NA.
#
# Los terminales son siempre intocables: se copian directamente desde
# estados_tips y el recorrido pre-order opera solo sobre nodos internos.
# ============================================================================ #
.fitch_optimizar <- function(phylo, res_desc, res_asc, estados_tips,
                             config = NULL) {
  s_down  <- res_desc$conjuntos
  s_up    <- res_asc$conjuntos
  n_tips  <- res_desc$n_tips
  n_nodes <- res_desc$n_nodes
  total   <- n_tips + n_nodes
  raiz    <- n_tips + 1L

  # ------------------------------------------------------------------
  # OPCIONES LOCALES: extraer modo desde config (nunca desde un global)
  # ------------------------------------------------------------------
  modo_actual <- tolower(config$modo_fitch %||% "deltran")

  # ------------------------------------------------------------------
  # Helper: calcula la optimización para ACCTRAN o DELTRAN
  # Argumento `m`: "acctran" o "deltran"
  # ------------------------------------------------------------------
  calc_modo <- function(m) {
    opt <- character(total)
    opt[seq_len(n_tips)] <- estados_tips  # Terminales: siempre intocables

    # La raíz toma el primer elemento del conjunto MPR (único punto de arranque)
    opt[raiz] <- s_up[[raiz]][1]

    # Recorrido pre-order exclusivo sobre nodos internos
    preorder_int <- function(nd) {
      res   <- nd
      hijos <- phylo$edge[phylo$edge[, 1] == nd, 2]
      for (h in hijos) if (h > n_tips) res <- c(res, preorder_int(h))
      res
    }
    orden_pre <- preorder_int(raiz)

    for (nd in orden_pre[-1]) {
      padre    <- phylo$edge[phylo$edge[, 2] == nd, 1][1]
      st_padre <- opt[padre]

      c_nd_up   <- s_up[[nd]]
      c_nd_down <- s_down[[nd]]

      if (length(c_nd_up) == 1L) {
        # Nodo inambiguo en el MPR: asignación directa
        opt[nd] <- c_nd_up[1]
      } else if (m == "acctran") {
        # ACCTRAN: hereda padre si está en el conjunto DESCENDENTE
        if (!is.na(st_padre) && st_padre %in% c_nd_down) opt[nd] <- st_padre
        else opt[nd] <- c_nd_down[1]
      } else {
        # DELTRAN: hereda padre si está en el conjunto MPR ASCENDENTE
        if (!is.na(st_padre) && st_padre %in% c_nd_up) opt[nd] <- st_padre
        else opt[nd] <- c_nd_up[1]
      }
    }
    opt
  }

  # ------------------------------------------------------------------
  # Despacho según modo_actual
  # ------------------------------------------------------------------
  if (modo_actual %in% c("acctran", "deltran")) {
    return(calc_modo(modo_actual))
  }

  # Unambiguous: inambiguo EXCLUSIVAMENTE si ACCTRAN y DELTRAN coinciden
  acc <- calc_modo("acctran")
  del <- calc_modo("deltran")
  opt_unamb <- character(total)

  # Terminales: siempre intocables
  opt_unamb[seq_len(n_tips)] <- estados_tips

  for (nd in (n_tips + 1L):total) {
    if (nd == raiz) {
      # La raíz es inambigua solo si su conjunto MPR tiene exactamente un elemento
      if (length(s_up[[raiz]]) == 1L) opt_unamb[raiz] <- s_up[[raiz]][1]
      else                            opt_unamb[raiz] <- NA_character_
    } else {
      # Nodo interno: inambiguo si y solo si ACCTRAN == DELTRAN (sin NA)
      if (!is.na(acc[nd]) && !is.na(del[nd]) && acc[nd] == del[nd]) {
        opt_unamb[nd] <- acc[nd]
      } else {
        opt_unamb[nd] <- NA_character_
      }
    }
  }
  opt_unamb
}


# ============================================================================ #
# FUNCIÓN PÚBLICA — firma compatible con MultiMapR
#
#   algoritmo_externo(tree, tip_colors, edge_colors, config) → edge_colors
#
# tip_colors : vector nombrado (nombre = estado del tip, valor = color).
#              Si names(tip_colors) está disponible, se usa directamente como
#              paleta estado→color, incluido el color asignado a "?" o "-".
# edge_colors: vector de colores inicializado (se sobreescribe completamente).
# ============================================================================ #
algoritmo_externo <- function(tree, tip_colors, edge_colors, config = NULL) {
  n_tips <- Ntip(tree)

  estados_tips <- names(tip_colors)

  if (!is.null(estados_tips)) {
    # Extraer paleta directa desde los nombres del vector
    paleta <- setNames(as.character(tip_colors), estados_tips)
    paleta <- paleta[!duplicated(names(paleta))]

    # Conservar el color asignado a la ambigüedad / datos faltantes
    color_ambiguo <- "gray70"
    for (m in c("?", "-")) {
      if (m %in% names(paleta) && paleta[[m]] != "gray70") {
        color_ambiguo <- paleta[[m]]
        break
      }
    }
  } else {
    # Fallback por seguridad cuando tip_colors no tiene nombres
    colores_unicos <- unique(tip_colors[tip_colors != "gray70"])
    estados_etiq   <- as.character(seq_len(length(colores_unicos)) - 1L)
    paleta         <- setNames(colores_unicos, estados_etiq)
    color_a_estado <- setNames(estados_etiq, colores_unicos)

    estados_tips <- vapply(tip_colors, function(col) {
      if (col == "gray70" || is.na(col)) "?"
      else color_a_estado[col]
    }, character(1))
    names(estados_tips) <- NULL
    color_ambiguo <- "gray70"
  }

  res_desc <- .fitch_paso_descendente(tree, estados_tips)
  res_asc  <- .fitch_paso_ascendente(tree, res_desc)
  estados  <- .fitch_optimizar(tree, res_desc, res_asc, estados_tips, config = config)

  for (i in seq_len(nrow(tree$edge))) {
    hijo        <- tree$edge[i, 2]
    estado_hijo <- estados[hijo]

    if (!is.na(estado_hijo) && estado_hijo %in% names(paleta)) {
      edge_colors[i] <- paleta[estado_hijo]
    } else {
      # Usa el color preservado de "?" en lugar de forzar siempre gris
      edge_colors[i] <- color_ambiguo
    }
  }

  edge_colors
}
