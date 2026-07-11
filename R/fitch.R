################################################################################
# fitch.R
#
# Adapter for the Fitch algorithm for use as 'external_algorithm'
# within the MultiMapR package.
#
# AMBIGUITY COLOR
#   All ambiguous / missing-data states (?, -, and unresolved internal nodes)
#   are ALWAYS rendered in FITCH_AMBIG_COLOR regardless of any color the user
#   may have assigned to those states in the palette.  This is intentional:
#   Fitch unresolved states carry a specific analytical meaning that must be
#   visually distinguishable from any real character state.

#' Fixed color for Fitch ambiguity / missing data.
#' Used by external_algorithm() and exported so core_render can add it to the legend.
FITCH_AMBIG_COLOR <- "#FF1493"   # deep pink / fuchsia
#
# When MultiMapR asks for the algorithm, select option 2 (Fitch).
# The optimization mode is chosen interactively from the terminal:
#   1 = ACCTRAN     (accelerated transformations -- toward the tips)
#   2 = DELTRAN     (delayed transformations -- toward the root)
#   3 = Unambiguous (only nodes with a unique state after both Fitch passes)
#
# SIGNATURE expected by MultiMapR:
#   external_algorithm(tree, tip_colors, edge_colors, config) -> edge_colors
#
#   - tree        : phylo object
#   - tip_colors  : named vector (color per tip, in tip.label order)
#   - edge_colors : vector initialized with "gray70" (length == nrow(tree$edge))
#   - config      : MultiMapR configuration list; config$fitch_mode is used
#                   (string: "acctran", "deltran" or "unambiguous") to
#                   determine the optimization method. Default: "deltran".
#
# NOTE: tip_colors already contains the color assigned by the user to each state;
#       this function inverts that color->state mapping to run Fitch and then
#       returns the resulting branch colors in the same format.
################################################################################

# ============================================================================ #
# DOWNPASS (tips -> root)
# ============================================================================ #
.fitch_downpass <- function(phylo, tip_states) {
  n_tips      <- Ntip(phylo)
  n_nodes     <- phylo$Nnode
  total_nodes <- n_tips + n_nodes

  all_states <- sort(unique(tip_states[!tip_states %in% c("?", "-")]))
  sets       <- vector("list", total_nodes)
  operations <- character(total_nodes)

  for (i in seq_len(n_tips)) {
    st <- tip_states[i]
    if (is.na(st) || st == "?" || st == "-" || trimws(st) == "") {
      sets[[i]] <- all_states
    } else {
      if (grepl("^[0-9]+$", st) && nchar(st) > 1) {
        sets[[i]] <- unique(strsplit(st, "")[[1]])
      } else if (grepl("[,/&|]", st)) {
        sets[[i]] <- unique(trimws(strsplit(st, "[,/&|]")[[1]]))
      } else {
        sets[[i]] <- st
      }
    }
    operations[i] <- "terminal"
  }

  internal_nodes <- (n_tips + 1):total_nodes
  count_desc <- function(nd) {
    if (nd <= n_tips) return(1L)
    children <- phylo$edge[phylo$edge[, 1] == nd, 2]
    sum(sapply(children, count_desc))
  }
  desc_counts    <- sapply(internal_nodes, count_desc)
  postorder_seq  <- internal_nodes[order(desc_counts)]

  for (nd in postorder_seq) {
    children <- phylo$edge[phylo$edge[, 1] == nd, 2]

    # Swofford & Maddison rule: frequency count instead of binary intersect/union
    child_states <- unlist(lapply(children, function(h) sets[[h]]))

    if (length(child_states) > 0) {
      state_counts <- table(child_states)
      max_freq     <- max(state_counts)
      sets[[nd]]   <- names(state_counts)[state_counts == max_freq]

      # If the max frequency equals the number of children, it was a perfect intersection
      operations[nd] <- if (max_freq == length(children)) "inter" else "union"
    } else {
      # Escape hatch: avoid undefined sets when children resolve to no states
      sets[[nd]] <- character(0)
      operations[nd] <- "union"
    }
  }

  list(sets       = sets,
       operations = operations,
       n_tips     = n_tips,
       n_nodes    = n_nodes)
}


# ============================================================================ #
# UPPASS (root -> tips)
#
# Implements the Swofford & Maddison (1987) rule to compute the MPR sets
# for each node:
#   - If the node was a UNION in the downpass:
#       MPR(nd) = union(s_down[nd], MPR(parent))
#   - If the node was an INTERSECTION in the downpass:
#       MPR(nd) = union(s_down[nd], intersect(MPR(parent), child_states))
#       where child_states = states of all children in s_down
# ============================================================================ #
.fitch_uppass <- function(phylo, down_result) {
  s_down     <- down_result$sets
  n_tips     <- down_result$n_tips
  root       <- n_tips + 1L
  s_up       <- s_down

  preorder_int <- function(nd) {
    res      <- nd
    children <- phylo$edge[phylo$edge[, 1] == nd, 2]
    for (h in children) if (h > n_tips) res <- c(res, preorder_int(h))
    res
  }
  preorder_seq <- preorder_int(root)

  for (nd in preorder_seq[-1]) {
    parent <- phylo$edge[phylo$edge[, 2] == nd, 1][1]

    c_node   <- s_down[[nd]]
    c_parent <- s_up[[parent]]

    children     <- phylo$edge[phylo$edge[, 1] == nd, 2]
    child_states <- unlist(lapply(children, function(h) s_down[[h]]))

    # Swofford & Maddison rule for the uppass
    if (length(child_states) > 0) {
      state_counts <- table(child_states)
      max_freq     <- max(state_counts)

      valid_parent_states <- c_parent[c_parent %in% names(state_counts)[state_counts >= (max_freq - 1)]]
      s_up[[nd]] <- union(c_node, valid_parent_states)
    } else {
      s_up[[nd]] <- union(c_node, c_parent)
    }
  }

  down_result$sets <- s_up
  down_result
}


# ============================================================================ #
# OPTIMIZATION (ACCTRAN / DELTRAN / Unambiguous) -> one unique state per node
#
# current_mode (derived from config$fitch_mode inside the function):
#   "acctran"     -- on ambiguities, uses DOWNPASS sets;
#                   if the parent state is in the child's set,
#                   it is inherited (accelerates changes toward the tips).
#   "deltran"     -- on ambiguities, uses MPR UPPASS sets;
#                   if the parent state is in the child's set,
#                   it is inherited (delays changes toward the root).
#   "unambiguous" -- a node is unambiguous EXCLUSIVELY when ACCTRAN and
#                   DELTRAN assign exactly the same state; the rest
#                   remain as NA. The root is unambiguous only if its MPR
#                   set has exactly one element.
#
# STRICT BUSINESS RULE (unambiguous):
#   A node is unambiguous if and only if acc[nd] == del[nd] (same string, no NA).
#   This condition is not relaxed: if they share a state it is unambiguous; if
#   they differ -- even if both are valid -- the node remains NA.
#
# Terminals are always untouched: copied directly from tip_states,
# and the pre-order traversal operates only on internal nodes.
# ============================================================================ #
.fitch_optimize <- function(phylo, down_result, up_result, tip_states, config = NULL) {
  s_down  <- down_result$sets
  s_up    <- up_result$sets
  n_tips  <- down_result$n_tips
  n_nodes <- down_result$n_nodes
  total   <- n_tips + n_nodes
  root    <- n_tips + 1L

  current_mode <- tolower(config$fitch_mode %||% "deltran")

  calc_mode <- function(m) {
    opt <- character(total)
    opt[seq_len(n_tips)] <- tip_states

    # Strict tie-break with sort() at the root
    opt[root] <- sort(s_up[[root]])[1]

    preorder_int <- function(nd) {
      res      <- nd
      children <- phylo$edge[phylo$edge[, 1] == nd, 2]
      for (h in children) if (h > n_tips) res <- c(res, preorder_int(h))
      res
    }
    preorder_seq <- preorder_int(root)

    for (nd in preorder_seq[-1]) {
      parent    <- phylo$edge[phylo$edge[, 2] == nd, 1][1]
      st_parent <- opt[parent]

      c_nd_up   <- s_up[[nd]]
      c_nd_down <- s_down[[nd]]

      # sort() applied to guarantee deterministic alphabetic order on ties
      if (length(c_nd_up) == 1L) {
        opt[nd] <- c_nd_up[1]
      } else if (m == "acctran") {
        if (!is.na(st_parent) && st_parent %in% c_nd_down) opt[nd] <- st_parent
        else opt[nd] <- sort(c_nd_down)[1]
      } else {
        if (!is.na(st_parent) && st_parent %in% c_nd_up) opt[nd] <- st_parent
        else opt[nd] <- sort(c_nd_up)[1]
      }
    }
    opt
  }

  if (current_mode %in% c("acctran", "deltran")) {
    return(calc_mode(current_mode))
  }

  acc       <- calc_mode("acctran")
  del       <- calc_mode("deltran")
  opt_unamb <- character(total)

  opt_unamb[seq_len(n_tips)] <- tip_states

  for (nd in (n_tips + 1L):total) {
    if (nd == root) {
      # Strict rigor at the root for Unambiguous
      if (length(s_up[[root]]) == 1L) opt_unamb[root] <- s_up[[root]][1]
      else                            opt_unamb[root] <- NA_character_
    } else {
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
# PUBLIC FUNCTION -- signature compatible with MultiMapR
#
#   external_algorithm(tree, tip_colors, edge_colors, config) -> edge_colors
#
# AMBIGUITY POLICY (Fitch-specific):
#   Any state that is ?, -, or that cannot be unambiguously resolved
#   (NA after .fitch_optimize) is painted with FITCH_AMBIG_COLOR (deep pink).
#   This overrides whatever color the user assigned to those states so that
#   Fitch ambiguity is always visually distinct from real character states.
#   This behavior is EXCLUSIVE to this algorithm; other algorithms respect
#   the user palette in full.
# ============================================================================ #
external_algorithm <- function(tree, tip_colors, edge_colors, config = NULL) {
  n_tips <- Ntip(tree)

  tip_states <- names(tip_colors)

  if (!is.null(tip_states)) {
    # 1. Unify missing-data tokens: treat inapplicable "-" as unknown "?"
    tip_states[tip_states == "-"] <- "?"
    names(tip_colors) <- tip_states

    # 2. Build palette from vector names.
    #    Strip "?" so it never appears as a resolved state color.
    palette <- setNames(as.character(tip_colors), tip_states)
    palette <- palette[!duplicated(names(palette))]
    palette <- palette[names(palette) != "?"]

    # 3. Ambiguity / missing data: always FITCH_AMBIG_COLOR (user choice ignored)
    ambiguous_color <- FITCH_AMBIG_COLOR

  } else {
    # Safety fallback when tip_colors has no names
    unique_colors  <- unique(tip_colors[tip_colors != "gray70"])
    state_labels   <- as.character(seq_len(length(unique_colors)) - 1L)
    palette        <- setNames(unique_colors, state_labels)
    color_to_state <- setNames(state_labels, unique_colors)

    tip_states <- vapply(tip_colors, function(col) {
      if (col == "gray70" || is.na(col)) "?"
      else color_to_state[col]
    }, character(1))
    names(tip_states) <- NULL
    ambiguous_color   <- FITCH_AMBIG_COLOR   # consistent even in fallback
  }

  down_result <- .fitch_downpass(tree, tip_states)
  up_result   <- .fitch_uppass(tree, down_result)
  states      <- .fitch_optimize(tree, down_result, up_result, tip_states, config = config)

  for (i in seq_len(nrow(tree$edge))) {
    child       <- tree$edge[i, 2]
    child_state <- states[child]

    if (!is.na(child_state) && child_state %in% names(palette)) {
      edge_colors[i] <- palette[child_state]
    } else {
      # NA (unresolved after optimize), "?" or any unrecognised state -> fuchsia
      edge_colors[i] <- ambiguous_color
    }
  }

  edge_colors
}
