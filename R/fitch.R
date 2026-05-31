################################################################################
# fitch.R
#
# Adapter for the Fitch algorithm for use as 'external_algorithm'
# within the MultiMapR package.
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

  # Unique real states (excludes "?" and "-")
  all_states <- sort(unique(tip_states[!tip_states %in% c("?", "-")]))

  sets       <- vector("list", total_nodes)
  operations <- character(total_nodes)

  # Initialize terminals
  for (i in seq_len(n_tips)) {
    st <- tip_states[i]
    if (is.na(st) || st == "?" || st == "-" || trimws(st) == "") {
      sets[[i]] <- all_states          # unknown = all states
    } else {
      if (grepl("^[0-9]+$", st) && nchar(st) > 1) {
        # Numeric polymorphisms without separator (e.g. "01" -> "0", "1")
        sets[[i]] <- unique(strsplit(st, "")[[1]])
      } else if (grepl("[,/&|]", st)) {
        # Polymorphisms with explicit separator (e.g. "0/1" or "urban,forest")
        sets[[i]] <- unique(trimws(strsplit(st, "[,/&|]")[[1]]))
      } else {
        # Whole words or single-character states ("forest", "granivore", "0", "A")
        sets[[i]] <- st
      }
    }
    operations[i] <- "terminal"
  }

  # Sort internal nodes from tips to root (post-order by descendant count)
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
    # Generalized Fitch: works with polytomies (2 or more children)
    inter <- Reduce(intersect, lapply(children, function(h) sets[[h]]))

    if (length(inter) > 0) {
      sets[[nd]]       <- inter
      operations[nd]   <- "inter"
    } else {
      sets[[nd]]       <- Reduce(union, lapply(children, function(h) sets[[h]]))
      operations[nd]   <- "union"
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
  operations <- down_result$operations
  n_tips     <- down_result$n_tips

  root  <- n_tips + 1L
  s_up  <- s_down  # Initialize MPR as a copy of the downpass

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
    c_parent <- s_up[[parent]]  # Definitive MPR of the parent

    if (operations[nd] == "union") {
      s_up[[nd]] <- union(c_node, c_parent)
    } else {
      children      <- phylo$edge[phylo$edge[, 1] == nd, 2]
      child_states  <- unlist(lapply(children, function(h) s_down[[h]]))
      valid         <- intersect(c_parent, child_states)
      s_up[[nd]]    <- union(c_node, valid)
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
.fitch_optimize <- function(phylo, down_result, up_result, tip_states,
                            config = NULL) {
  s_down  <- down_result$sets
  s_up    <- up_result$sets
  n_tips  <- down_result$n_tips
  n_nodes <- down_result$n_nodes
  total   <- n_tips + n_nodes
  root    <- n_tips + 1L

  # ------------------------------------------------------------------
  # LOCAL OPTIONS: extract mode from config (never from a global)
  # ------------------------------------------------------------------
  current_mode <- tolower(config$fitch_mode %||% "deltran")

  # ------------------------------------------------------------------
  # Helper: computes the optimization for ACCTRAN or DELTRAN
  # Argument `m`: "acctran" or "deltran"
  # ------------------------------------------------------------------
  calc_mode <- function(m) {
    opt <- character(total)
    opt[seq_len(n_tips)] <- tip_states  # Terminals: always untouched

    # Root takes the first element of the MPR set (single starting point)
    opt[root] <- s_up[[root]][1]

    # Pre-order traversal exclusively over internal nodes
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

      if (length(c_nd_up) == 1L) {
        # Unambiguous node in MPR: direct assignment
        opt[nd] <- c_nd_up[1]
      } else if (m == "acctran") {
        # ACCTRAN: inherit parent if it is in the DOWNPASS set
        if (!is.na(st_parent) && st_parent %in% c_nd_down) opt[nd] <- st_parent
        else opt[nd] <- c_nd_down[1]
      } else {
        # DELTRAN: inherit parent if it is in the MPR UPPASS set
        if (!is.na(st_parent) && st_parent %in% c_nd_up) opt[nd] <- st_parent
        else opt[nd] <- c_nd_up[1]
      }
    }
    opt
  }

  # ------------------------------------------------------------------
  # Dispatch according to current_mode
  # ------------------------------------------------------------------
  if (current_mode %in% c("acctran", "deltran")) {
    return(calc_mode(current_mode))
  }

  # Unambiguous: unambiguous EXCLUSIVELY if ACCTRAN and DELTRAN agree
  acc      <- calc_mode("acctran")
  del      <- calc_mode("deltran")
  opt_unamb <- character(total)

  # Terminals: always untouched
  opt_unamb[seq_len(n_tips)] <- tip_states

  for (nd in (n_tips + 1L):total) {
    if (nd == root) {
      # Root is unambiguous only if its MPR set has exactly one element
      if (length(s_up[[root]]) == 1L) opt_unamb[root] <- s_up[[root]][1]
      else                             opt_unamb[root] <- NA_character_
    } else {
      # Internal node: unambiguous if and only if ACCTRAN == DELTRAN (no NA)
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
# tip_colors : named vector (name = tip state, value = color).
#              If names(tip_colors) is available, it is used directly as the
#              state->color palette, including the color assigned to "?" or "-".
# edge_colors: initialized color vector (completely overwritten).
# ============================================================================ #
external_algorithm <- function(tree, tip_colors, edge_colors, config = NULL) {
  n_tips <- Ntip(tree)

  tip_states <- names(tip_colors)

  if (!is.null(tip_states)) {
    # Extract palette directly from vector names
    palette <- setNames(as.character(tip_colors), tip_states)
    palette <- palette[!duplicated(names(palette))]

    # Preserve the color assigned to ambiguity / missing data
    ambiguous_color <- "gray70"
    for (m in c("?", "-")) {
      if (m %in% names(palette) && palette[[m]] != "gray70") {
        ambiguous_color <- palette[[m]]
        break
      }
    }
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
    ambiguous_color   <- "gray70"
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
      # Use the preserved color for "?" instead of always forcing gray
      edge_colors[i] <- ambiguous_color
    }
  }

  edge_colors
}