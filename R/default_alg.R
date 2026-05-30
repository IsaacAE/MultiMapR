################################################################################
# default_alg.R
#
# Default ancestral reconstruction algorithm for MultiMapR.
# Assigns colors to internal branches by depth-weighted majority vote.
#
# Exported function:
#   default_algorithm(tree, tip_colors, edge_colors) -> updated edge_colors
#
# Loaded automatically from main_MultiMapR.R when the user
# selects algorithm 1 (default).
################################################################################


#' Gets all descendants of a node (iterative)
#'
#' Iterative version of getDescendants() that avoids recursion problems
#' in large trees and does not use <<-.
#'
#' @param tree  phylo object.
#' @param node  Index of the starting node.
#' @return Integer vector of descendant node indices.
.get_descendants <- function(tree, node) {
  edges   <- tree$edge
  pending <- node
  result  <- integer(0)

  while (length(pending) > 0) {
    nd      <- pending[1]
    pending <- pending[-1]
    children <- edges[edges[, 1] == nd, 2]
    result  <- c(result, children)
    pending <- c(pending, children)
  }
  result
}


#' Default algorithm: colors branches by depth-weighted majority vote
#'
#' Steps:
#'   1. Colors the edges leading to each tip with that tip's color.
#'   2. Traverses internal nodes in postorder (deepest first).
#'   3. For each internal node, determines the predominant color among
#'      its descendants weighted by depth; that color is assigned to
#'      the edge ENTERING the node.
#'   4. After assigning all nodes, propagates each internal node's color
#'      to the vertical edges connecting it to its parent, preventing
#'      segments from remaining gray when the parent had no color yet.
#'
#' @param tree        phylo object.
#' @param tip_colors  Named character vector: tip_colors[i] = color of tip i.
#' @param edge_colors Named character vector initialised with gray
#'                    (names: "parent-child").
#' @return Named character vector of length nrow(tree$edge) with updated
#'         colors for each edge.
default_algorithm <- function(tree, tip_colors, edge_colors) {
  n_tips  <- Ntip(tree)
  n_nodes <- Nnode(tree)
  edges   <- tree$edge

  # ------------------------------------------------------------------
  # Step 1: edges leading to tips
  # ------------------------------------------------------------------
  for (i in seq_len(n_tips)) {
    idx              <- which(edges[, 2] == i)
    edge_colors[idx] <- tip_colors[i]
  }

  # ------------------------------------------------------------------
  # Step 2: postorder traversal of internal nodes (deepest first)
  # ------------------------------------------------------------------
  all_nodes    <- (n_tips + 1):(n_tips + n_nodes)
  node_depths  <- node.depth.edgelength(tree)
  sorted_nodes <- all_nodes[order(node_depths[all_nodes], decreasing = TRUE)]

  # Store the color assigned to each node for use in step 4
  node_color <- character(n_tips + n_nodes)
  node_color[seq_len(n_tips)] <- tip_colors

  for (nd in sorted_nodes) {
    desc <- .get_descendants(tree, nd)
    if (length(desc) == 0) next

    desc_colors     <- character(0)
    internal_nodes  <- integer(0)

    for (d in desc) {
      if (d <= n_tips) {
        desc_colors <- c(desc_colors, tip_colors[d])
      } else {
        idx            <- which(edges[, 2] == d)
        desc_colors    <- c(desc_colors, edge_colors[idx])
        internal_nodes <- c(internal_nodes, d)
      }
    }

    # Priority color: the one from the descendant internal node with the
    # fewest descendants of its own (i.e. the closest in the hierarchy)
    if (length(internal_nodes) > 0) {
      n_desc     <- sapply(internal_nodes,
                           function(x) length(.get_descendants(tree, x)))
      prio_node  <- internal_nodes[which.min(n_desc)]
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
  # Step 3: propagate colors to edges between internal nodes that
  # remained gray.
  #
  # Problem: in a phylogram, the parent->child edge has a vertical
  # segment (the "elbow") that ape draws with that edge's color.
  # If the edge remained gray because the child node was still gray
  # when the parent was processed, the vertical stays uncolored.
  #
  # Solution: for each parent->child edge where BOTH are internal nodes,
  # if the color is gray but the child node has an assigned color,
  # use the child's color (which reflects what predominates in its subtree).
  # ------------------------------------------------------------------
  for (i in seq_len(nrow(edges))) {
    parent <- edges[i, 1]
    child  <- edges[i, 2]

    # Only applies to edges between internal nodes
    if (child <= n_tips) next

    # If it already has a real color, leave it
    if (edge_colors[i] != "gray70") next

    # Use the child node's color if available
    child_color <- node_color[child]
    if (nchar(child_color) > 0 && child_color != "gray70") {
      edge_colors[i] <- child_color
    }
  }

  # ------------------------------------------------------------------
  # Step 4: color the edges leaving the root node.
  #
  # The root has no incoming edge, so node_color[root] stays empty
  # and ape draws the initial vertical segment in gray.
  # We assign to each root->child edge that is still gray the color
  # of the corresponding child node (already computed in previous steps).
  # ------------------------------------------------------------------
  root       <- n_tips + 1L
  root_children <- edges[edges[, 1] == root, 2]

  for (child in root_children) {
    edge_idx <- which(edges[, 1] == root & edges[, 2] == child)
    if (edge_colors[edge_idx] == "gray70") {
      child_color <- node_color[child]
      if (nchar(child_color) > 0 && child_color != "gray70")
        edge_colors[edge_idx] <- child_color
    }
  }

  edge_colors
}