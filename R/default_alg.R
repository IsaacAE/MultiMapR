################################################################################
# algoritmo_default.R
#
# Algoritmo de reconstrucción ancestral por defecto para MultiMapR.
# Asigna colores a ramas internas por mayoría ponderada por profundidad.
#
# Función exportada:
#   algoritmo_default(tree, tip_colors, edge_colors) → edge_colors actualizado
#
# Se carga automáticamente desde main_MultiMapR.R cuando el usuario
# selecciona el algoritmo 1 (por defecto).
################################################################################


#' Obtiene todos los descendientes de un nodo (iterativo)
#'
#' Versión iterativa de getDescendants() que evita problemas de recursión
#' en árboles grandes y no usa <<-.
#'
#' @param tree  Objeto phylo.
#' @param node  Número del nodo de partida.
#' @return Vector de índices de nodos descendientes.
.get_descendants <- function(tree, node) {
  edges       <- tree$edge
  pendientes  <- node
  result      <- integer(0)

  while (length(pendientes) > 0) {
    nd         <- pendientes[1]
    pendientes <- pendientes[-1]
    hijos      <- edges[edges[, 1] == nd, 2]
    result     <- c(result, hijos)
    pendientes <- c(pendientes, hijos)
  }
  result
}


#' Algoritmo por defecto: colorea ramas por mayoría ponderada por profundidad
#'
#' Pasos:
#'   1. Colorea las ramas que llegan a cada tip con el color del tip.
#'   2. Recorre los nodos internos en postorder (hojas primero).
#'   3. Para cada nodo interno, determina el color predominante entre
#'      sus descendientes directos ponderado por profundidad; ese color
#'      se asigna a la arista que ENTRA al nodo.
#'   4. Tras asignar todos los nodos, propaga el color de cada nodo interno
#'      a las aristas verticales que lo conectan con su padre, evitando
#'      que queden segmentos en gris cuando el padre aún no tenía color.
#'
#' @param tree        Objeto phylo.
#' @param tip_colors  Named character vector: `tip_colors[i]` = color del tip i.
#' @param edge_colors Named character vector inicializado con gris
#'                    (nombres: "padre-hijo").
#' @return Named character vector de longitud nrow(tree$edge) con colores
#'         actualizados para cada arista.
algoritmo_default <- function(tree, tip_colors, edge_colors) {
  n_tips   <- Ntip(tree)
  n_nodes  <- Nnode(tree)
  edges    <- tree$edge

  # ------------------------------------------------------------------
  # Paso 1: ramas que llegan a tips
  # ------------------------------------------------------------------
  for (i in seq_len(n_tips)) {
    idx              <- which(edges[, 2] == i)
    edge_colors[idx] <- tip_colors[i]
  }

  # ------------------------------------------------------------------
  # Paso 2: postorder de nodos internos (más profundos primero)
  # ------------------------------------------------------------------
  all_nodes    <- (n_tips + 1):(n_tips + n_nodes)
  node_depths  <- node.depth.edgelength(tree)
  sorted_nodes <- all_nodes[order(node_depths[all_nodes], decreasing = TRUE)]

  # Guardamos el color asignado a cada nodo para usarlo en el paso 4
  node_color <- character(n_tips + n_nodes)
  node_color[seq_len(n_tips)] <- tip_colors

  for (nd in sorted_nodes) {
    desc <- .get_descendants(tree, nd)
    if (length(desc) == 0) next

    desc_colors    <- character(0)
    internal_nodes <- integer(0)

    for (d in desc) {
      if (d <= n_tips) {
        desc_colors <- c(desc_colors, tip_colors[d])
      } else {
        idx            <- which(edges[, 2] == d)
        desc_colors    <- c(desc_colors, edge_colors[idx])
        internal_nodes <- c(internal_nodes, d)
      }
    }

    # Color prioritario: el del nodo interno descendiente con menos
    # descendientes propios (el más cercano en la jerarquía)
    if (length(internal_nodes) > 0) {
      n_div      <- sapply(internal_nodes,
                           function(x) length(.get_descendants(tree, x)))
      prio_node  <- internal_nodes[which.min(n_div)]
      prio_idx   <- which(edges[, 2] == prio_node)
      prio_color <- edge_colors[prio_idx]
    } else {
      prio_color <- NA_character_
    }

    predominant <- if (is.na(prio_color) || length(internal_nodes) == 0) {
      names(sort(table(desc_colors), decreasing = TRUE))[1]
    } else {
      prio_color
    }

    idx              <- which(edges[, 2] == nd)
    edge_colors[idx] <- predominant
    node_color[nd]   <- predominant
  }

  # ------------------------------------------------------------------
  # Paso 3: propagar colores a aristas entre nodos internos que hayan
  # quedado en gris.
  #
  # Problema: en un phylogram, la arista padre→hijo tiene un segmento
  # vertical (el "codo") que ape dibuja con el color de esa arista.
  # Si la arista quedó en gris porque el nodo hijo aún era gris cuando
  # se procesó el padre, la vertical queda sin color.
  #
  # Solución: para cada arista padre→hijo donde AMBOS son nodos internos,
  # si el color es gris pero el nodo hijo tiene un color asignado,
  # usar el color del hijo (que refleja lo que predomina en su subárbol).
  # ------------------------------------------------------------------
  for (i in seq_len(nrow(edges))) {
    padre <- edges[i, 1]
    hijo  <- edges[i, 2]

    # Solo aplica a aristas entre nodos internos
    if (hijo <= n_tips) next

    # Si ya tiene un color real, no tocar
    if (edge_colors[i] != "gray70") next

    # Usar el color del nodo hijo si está disponible
    color_hijo <- node_color[hijo]
    if (nchar(color_hijo) > 0 && color_hijo != "gray70") {
      edge_colors[i] <- color_hijo
    }
  }

  # ------------------------------------------------------------------
  # Paso 4: colorear las aristas que salen del nodo raíz.
  #
  # La raíz no tiene arista entrante, por lo que node_color[raiz]
  # queda vacío y ape dibuja el segmento vertical inicial en gris.
  # Asignamos a cada arista raíz→hijo que siga en gris el color
  # del nodo hijo correspondiente (ya calculado en pasos anteriores).
  # ------------------------------------------------------------------
  raiz <- n_tips + 1L
  hijos_raiz <- edges[edges[, 1] == raiz, 2]

  for (hijo in hijos_raiz) {
    idx_arista <- which(edges[, 1] == raiz & edges[, 2] == hijo)
    if (edge_colors[idx_arista] == "gray70") {
      color_hijo <- node_color[hijo]
      if (nchar(color_hijo) > 0 && color_hijo != "gray70")
        edge_colors[idx_arista] <- color_hijo
    }
  }

  edge_colors
}
