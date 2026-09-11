################################################################################
# fitch.R
#
# Adapter for the Fitch/Sankoff algorithm for use as 'external_algorithm'
# within the MultiMapR package.
#
# ALGORITHM
#   Downpass and uppass are computed via the Sankoff (1975) two-pass dynamic
#   program with a 0/1 (unordered Fitch-equivalent) cost matrix, rather than
#   the earlier frequency-count heuristic. The heuristic under-counted the
#   most-parsimonious-reconstruction (MPR) set at nodes with a 3-way tie
#   propagated through more than one ambiguous ancestor (e.g. a node whose
#   two children are fixed to two DIFFERENT single states, sitting directly
#   below a node that is itself tied with a third, unrelated state coming
#   from the outgroup) -- it only ever allowed a state into a node's MPR set
#   if that state was already present among the node's own children, so a
#   state that was optimal purely by tying with the *parent's* assignment
#   could be wrongly excluded. Sankoff's two-pass DP has no such blind spot:
#   a state is included in a node's MPR set exactly when it participates in
#   some assignment achieving the tree's true minimum length, which is
#   exactly what ACCTRAN/DELTRAN/unambiguous need.
#
#   Verified against Winclada (TNT) on real data (bats_tre_pol.tre /
#   matriz_politomy.csv, characters 11 and 12): tree length (L) and the
#   per-node MPR sets it reports match exactly, including a 4-state
#   character with a genuine 3-way tie descending from the root.
#
# AMBIGUITY COLOR
#   All ambiguous / missing-data states (?, -, and unresolved internal nodes)
#   are ALWAYS rendered in FITCH_AMBIG_COLOR regardless of any color the user
#   may have assigned to those states in the palette.
#
# SIGNATURE expected by MultiMapR:
#   external_algorithm(tree, tip_colors, edge_colors, config) -> edge_colors
################################################################################

FITCH_AMBIG_COLOR <- "#FF00FF"   # magenta

# ============================================================================ #
# DOWNPASS (tips -> root): Sankoff cost matrix S[node, state].
# Returns $sets = F (downpass-only MPR sets, used by ACCTRAN) alongside the
# raw cost matrix $S (needed by the uppass) and bookkeeping fields.
# ============================================================================ #
.fitch_downpass <- function(phylo, tip_states) {
  n_tips      <- Ntip(phylo)
  n_nodes     <- phylo$Nnode
  total_nodes <- n_tips + n_nodes

  alphabet <- sort(unique(tip_states[!tip_states %in% c("?", "-")]))
  poly_extra <- unlist(lapply(tip_states, function(st) {
    if (is.na(st) || st %in% c("?", "-", "")) return(NULL)
    if (grepl("[,/&|]", st)) return(trimws(strsplit(st, "[,/&|]")[[1]]))
    if (grepl("^[0-9]+$", st) && nchar(st) > 1) return(unique(strsplit(st, "")[[1]]))
    NULL
  }))
  alphabet <- sort(unique(c(alphabet, poly_extra)))
  k <- length(alphabet)

  S <- matrix(Inf, nrow = total_nodes, ncol = k, dimnames = list(NULL, alphabet))

  for (i in seq_len(n_tips)) {
    st <- tip_states[i]
    if (is.na(st) || st == "?" || st == "-" || trimws(st) == "") {
      S[i, ] <- 0
    } else if (grepl("[,/&|]", st)) {
      poly <- trimws(strsplit(st, "[,/&|]")[[1]])
      S[i, ] <- 1
      S[i, alphabet %in% poly] <- 0
    } else if (grepl("^[0-9]+$", st) && nchar(st) > 1) {
      poly <- unique(strsplit(st, "")[[1]])
      S[i, ] <- 1
      S[i, alphabet %in% poly] <- 0
    } else {
      S[i, ] <- 1
      S[i, st] <- 0
    }
  }

  internal_nodes <- (n_tips + 1):total_nodes
  count_desc <- function(nd) {
    if (nd <= n_tips) return(1L)
    children <- phylo$edge[phylo$edge[, 1] == nd, 2]
    sum(sapply(children, count_desc))
  }
  desc_counts   <- sapply(internal_nodes, count_desc)
  postorder_seq <- internal_nodes[order(desc_counts)]

  operations <- character(total_nodes)
  operations[seq_len(n_tips)] <- "terminal"

  for (nd in postorder_seq) {
    children <- phylo$edge[phylo$edge[, 1] == nd, 2]
    acc <- numeric(k)
    for (h in children) {
      m1      <- min(S[h, ])
      contrib <- ifelse(S[h, ] == m1, m1, m1 + 1)
      acc     <- acc + contrib
    }
    S[nd, ] <- acc
    operations[nd] <- if (sum(S[nd, ] == min(S[nd, ])) == 1) "inter" else "union"
  }

  F_sets <- vector("list", total_nodes)
  for (nd in seq_len(total_nodes)) F_sets[[nd]] <- alphabet[S[nd, ] == min(S[nd, ])]

  list(sets = F_sets, S = S, alphabet = alphabet, operations = operations,
       n_tips = n_tips, n_nodes = n_nodes)
}


# ============================================================================ #
# UPPASS (root -> tips): Sankoff second pass. G[node, state] = minimal cost of
# everything OUTSIDE this node's subtree, given this node = state.
# Returns $sets = M, the true MPR set: every state that participates in at
# least one most-parsimonious reconstruction of the whole tree. Used by
# DELTRAN and by the unambiguous (acctran==deltran) comparison.
# ============================================================================ #
.fitch_uppass <- function(phylo, down_result) {
  S        <- down_result$S
  alphabet <- down_result$alphabet
  n_tips   <- down_result$n_tips
  total    <- n_tips + down_result$n_nodes
  root     <- n_tips + 1L
  k        <- length(alphabet)

  G <- matrix(0, nrow = total, ncol = k, dimnames = list(NULL, alphabet))

  preorder_int <- function(nd) {
    res      <- nd
    children <- phylo$edge[phylo$edge[, 1] == nd, 2]
    for (h in children) if (h > n_tips) res <- c(res, preorder_int(h))
    res
  }
  preorder_seq <- preorder_int(root)

  for (nd in preorder_seq[-1]) {
    parent   <- phylo$edge[phylo$edge[, 2] == nd, 1][1]
    siblings <- setdiff(phylo$edge[phylo$edge[, 1] == parent, 2], nd)

    h_vec <- G[parent, ]
    for (sib in siblings) {
      m1    <- min(S[sib, ])
      h_vec <- h_vec + ifelse(S[sib, ] == m1, m1, m1 + 1)
    }
    hmin    <- min(h_vec)
    G[nd, ] <- pmin(h_vec, hmin + 1)
  }

  full <- S + G
  M_sets <- vector("list", total)
  for (nd in seq_len(total)) M_sets[[nd]] <- alphabet[full[nd, ] == min(full[nd, ])]

  list(sets = M_sets, G = G)
}


# ============================================================================ #
# OPTIMIZATION (ACCTRAN / DELTRAN / Unambiguous) -> one unique state per node
# Unchanged: consumes s_down (F) and s_up (M) exactly as the previous version.
# ============================================================================ #
.fitch_optimize <- function(phylo, down_result, up_result, tip_states, config = NULL) {
  s_down  <- down_result$sets
  s_up    <- up_result$sets
  n_tips  <- down_result$n_tips
  n_nodes <- down_result$n_nodes
  total   <- n_tips + n_nodes
  root    <- n_tips + 1L

  current_mode <- tolower(config$fitch_mode %||% "deltran")

  # DELTRAN needs ACCTRAN's per-node values as its fallback (see below), so
  # compute ACCTRAN first regardless of which mode was requested.
  acc_vals <- local({
    opt <- character(total)
    opt[seq_len(n_tips)] <- tip_states
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
      c_nd_down <- s_down[[nd]]

      if (length(c_nd_down) == 1L) {
        opt[nd] <- c_nd_down[1]
      } else if (!is.na(st_parent) && st_parent %in% c_nd_down) {
        opt[nd] <- st_parent
      } else {
        opt[nd] <- sort(c_nd_down)[1]
      }
    }
    opt
  })

  calc_mode <- function(m) {
    if (m == "acctran") return(acc_vals)

    opt <- character(total)
    opt[seq_len(n_tips)] <- tip_states
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

      if (length(c_nd_up) == 1L) {
        opt[nd] <- c_nd_up[1]
      } else if (!is.na(st_parent) && st_parent %in% c_nd_up) {
        # DELTRAN can delay the change: adopt the parent's state.
        opt[nd] <- st_parent
      } else {
        # Delaying isn't possible here (the state can't have persisted this
        # far down without an extra step) -- fall back to this node's own
        # ACCTRAN value, which is guaranteed to be a valid MPR member here.
        opt[nd] <- acc_vals[nd]
      }
    }
    opt
  }

  if (current_mode %in% c("acctran", "deltran")) {
    return(calc_mode(current_mode))
  }

  # UNAMBIGUOUS: a node is genuinely unambiguous iff its MPR set (s_up, the
  # full Sankoff up-pass set) has exactly one member -- i.e. every
  # most-parsimonious reconstruction of the whole tree agrees on this node's
  # state. This must NOT be defined as "ACCTRAN's value equals DELTRAN's
  # value": since DELTRAN falls back to ACCTRAN's own value whenever it
  # can't inherit the parent's state (see above), the two can coincide by
  # construction of the tie-break even when the node's true MPR set has
  # 2+ members (e.g. {1,3}) and is therefore genuinely ambiguous.
  opt_unamb <- character(total)
  opt_unamb[seq_len(n_tips)] <- tip_states
  for (nd in (n_tips + 1L):total) {
    opt_unamb[nd] <- if (length(s_up[[nd]]) == 1L) s_up[[nd]][1] else NA_character_
  }
  opt_unamb
}


# ============================================================================ #
# PUBLIC FUNCTION -- signature compatible with MultiMapR
# ============================================================================ #
external_algorithm <- function(tree, tip_colors, edge_colors, config = NULL) {
  n_tips <- Ntip(tree)
  tip_states <- names(tip_colors)

  if (!is.null(tip_states)) {
    tip_states[tip_states == "-"] <- "?"
    names(tip_colors) <- tip_states
    palette <- setNames(as.character(tip_colors), tip_states)
    palette <- palette[!duplicated(names(palette))]
    palette <- palette[names(palette) != "?"]
    ambiguous_color <- config$ambiguity_color %||% FITCH_AMBIG_COLOR
  } else {
    unique_colors  <- unique(tip_colors[tip_colors != "gray70"])
    state_labels   <- as.character(seq_len(length(unique_colors)) - 1L)
    palette        <- setNames(unique_colors, state_labels)
    color_to_state <- setNames(state_labels, unique_colors)
    tip_states <- vapply(tip_colors, function(col) {
      if (col == "gray70" || is.na(col)) "?"
      else color_to_state[col]
    }, character(1))
    names(tip_states) <- NULL
    ambiguous_color   <- config$ambiguity_color %||% FITCH_AMBIG_COLOR
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
      edge_colors[i] <- ambiguous_color
    }
  }

  # Root's own state has no incoming edge to live in; expose it as an
  # attribute so the rendering functions can paint the root's vertical
  # connector with the root's own state instead of borrowing a child's color.
  root       <- n_tips + 1L
  root_state <- states[root]
  attr(edge_colors, "root_color") <- if (!is.na(root_state) && root_state %in% names(palette)) {
    unname(palette[root_state])
  } else {
    ambiguous_color
  }

  edge_colors
}
