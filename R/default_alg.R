# Default ancestral reconstruction: colors branches by depth-weighted majority vote
# of their descendants.
#   default_algorithm(tree, tip_colors, edge_colors) -> updated edge_colors


#' Gets all descendants of a node (iterative)
#' @noRd
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
#' @noRd
default_algorithm <- function(tree, tip_colors, edge_colors) {
  n_tips  <- Ntip(tree)
  n_nodes <- Nnode(tree)
  edges   <- tree$edge

  # Step 1: edges leading to tips
  for (i in seq_len(n_tips)) {
    idx              <- which(edges[, 2] == i)
    edge_colors[idx] <- tip_colors[i]
  }

  # Step 2: postorder traversal of internal nodes (deepest first)
  all_nodes    <- (n_tips + 1):(n_tips + n_nodes)
  # Without edge lengths ape reads uninitialised memory, so use unit lengths.
  depth_tree <- tree
  if (is.null(depth_tree$edge.length))
    depth_tree$edge.length <- rep(1, nrow(depth_tree$edge))
  node_depths  <- node.depth.edgelength(depth_tree)
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

  # Step 3: propagate colors to edges between internal nodes that stayed gray. ape
  # draws a phylogram's vertical elbow with the edge's color, so if the child was
  # still gray when the parent was processed the elbow would stay uncolored; use
  # the child's color (what predominates in its subtree).
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

  # Step 4: color the root->child edges that are still gray with the child's color
  # (the root has no incoming edge, so ape would draw its vertical in gray).
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
