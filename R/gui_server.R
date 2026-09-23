################################################################################
# gui_server.R
#
# Server logic of the MultiMapR interface (see gui.R / gui_ui.R).
#   Section 1  helpers (log entries, data loading)
#   Section 2  server: language, data page
#   Section 3  server: map page (characters, colors, tree, preview, export,
#              ACCTRAN / DELTRAN comparison)
################################################################################


# ==============================================================================
# SECTION 1 -- HELPERS
# ==============================================================================

#' Creates a log entry
#'
#' Entries keep the dictionary key and its arguments (not the text), so the
#' log is re-translated when the language changes. Raw package messages are
#' stored in \code{text}.
#' @param level "info", "warn", "error" or "ok".
#' @noRd
.gui_log_entry <- function(level, key = NULL, ..., text = NULL) {
  list(level = level, key = key, args = list(...), text = text,
       time = format(Sys.time(), "%H:%M:%S"))
}

#' Converts captured package output into log entries
#' @noRd
.gui_log_from_capture <- function(lines, rename = NULL) {
  lines <- unique(lines[nzchar(trimws(lines)) & !grepl("^\\s*===", lines)])
  if (!is.null(rename))
    for (k in seq_along(rename)) lines <- gsub(names(rename)[k], rename[[k]], lines, fixed = TRUE)
  lapply(lines, function(l) {
    lvl <- if (startsWith(l, "Warning:")) "warn" else "info"
    .gui_log_entry(lvl, text = sub("^\\[MultiMapR\\]\\s*", "", sub("^Warning:\\s*", "", l)))
  })
}

#' Renders a list of log entries
#' @noRd
.gui_render_log <- function(entries, tr) {
  tags <- shiny::tags
  if (length(entries) == 0)
    return(tags$div(class = "mm-log-empty", tr("log_empty")))
  lapply(entries, function(e) {
    txt <- if (!is.null(e$key)) do.call(tr, c(list(e$key), e$args)) else e$text
    tags$div(class = paste0("mm-log-line mm-", e$level), title = e$time,
             tags$span(class = "mm-lvl", tr(paste0("lvl_", e$level))),
             tags$span(txt))
  })
}

#' Number of warnings in a log, formatted
#' @noRd
.gui_log_meta <- function(entries, tr) {
  n <- sum(vapply(entries, function(e) identical(e$level, "warn"), logical(1)))
  if (n == 1) tr("n_warning") else tr("n_warnings", n)
}

#' Reads tree and matrix step by step so errors can say which file failed
#'
#' @param tree_path,matrix_path Paths on disk.
#' @param names  Display names (the uploaded files live under temporary names).
#' @param opts   list(tree_format, matrix_format, sep, species_col, header,
#'               normalize_spaces, prune).
#' @return list(ok, tree, characters, log, error = list(title_key, message, ntip)).
#' @noRd
.gui_read_inputs <- function(tree_path, matrix_path, names, opts) {
  log <- list()
  rename <- stats::setNames(as.list(names), c(tree_path, matrix_path))
  fail <- function(title_key, e, ntip = NULL) {
    msg <- conditionMessage(e)
    for (k in seq_along(rename)) msg <- gsub(names(rename)[k], rename[[k]], msg, fixed = TRUE)
    log[[length(log) + 1]] <<- .gui_log_entry("error", text = msg)
    log[[length(log) + 1]] <<- .gui_log_entry("error", "log_aborted")
    list(ok = FALSE, log = log,
         error = list(title_key = title_key, message = msg, ntip = ntip))
  }

  # -- 1. tree -----------------------------------------------------------------
  tree_fmt <- if (identical(opts$tree_format, "tnt")) "newick" else opts$tree_format
  res <- tryCatch(.gui_capture(read_tree(tree_path, tree_format = tree_fmt)),
                  error = function(e) e)
  if (inherits(res, "error")) return(fail("err_tree_title", res))
  tree <- res$value
  log  <- c(log, .gui_log_from_capture(res$log, rename))
  has_bl <- !is.null(tree$edge.length)
  log[[length(log) + 1]] <- .gui_log_entry(
    "info", if (has_bl) "log_tree_read_bl" else "log_tree_read", names[1], Ntip(tree))

  # -- 2. matrix ---------------------------------------------------------------
  sep <- switch(opts$sep, tab = "\t", opts$sep)
  mfmt <- if (identical(opts$matrix_format, "auto")) suppressMessages(.detect_matrix_format(matrix_path))
          else opts$matrix_format
  res <- tryCatch(.gui_capture(read_characters(matrix_path,
                                               format           = opts$matrix_format,
                                               sep              = sep,
                                               species_col      = opts$species_col,
                                               normalize_spaces = opts$normalize_spaces,
                                               header           = opts$header)),
                  error = function(e) e)
  if (inherits(res, "error")) return(fail("err_matrix_title", res, Ntip(tree)))
  characters <- res$value
  log <- c(log, .gui_log_from_capture(res$log, rename))
  log[[length(log) + 1]] <- .gui_log_entry("info", "log_matrix_read", names[2], toupper(mfmt),
                                           nrow(characters), ncol(characters) - 1L)

  # -- 3. names ----------------------------------------------------------------
  res <- .gui_capture(.harmonize_species(tree, characters))
  characters <- res$value
  log <- c(log, .gui_log_from_capture(res$log, rename))

  list(ok = TRUE, tree = tree, characters = characters, log = log)
}

#' Reconciles tree and matrix (optionally pruning tips absent from the matrix)
#'
#' @return list(ok, tree, characters, aligned, log) or list(ok = FALSE, log, error).
#' @noRd
.gui_reconcile <- function(tree, characters, prune) {
  log <- list()
  species <- as.character(characters$Species)
  missing_tips <- setdiff(tree$tip.label, species)
  extra_rows   <- setdiff(species, tree$tip.label)
  shown <- function(x) paste0(paste(utils::head(x, 6), collapse = ", "),
                              if (length(x) > 6) ", \u2026" else "")

  if (length(extra_rows) > 0)
    log[[length(log) + 1]] <- .gui_log_entry("info", "log_rows_extra", length(extra_rows),
                                             shown(extra_rows))
  if (length(missing_tips) > 0) {
    log[[length(log) + 1]] <- .gui_log_entry("warn", "log_tips_missing", length(missing_tips),
                                             shown(missing_tips))
    if (!prune || length(missing_tips) > Ntip(tree) - 3L) {
      log[[length(log) + 1]] <- .gui_log_entry("error", "log_aborted")
      return(list(ok = FALSE, log = log,
                  error = list(title_key = "err_match_title",
                               message_key = if (prune) "err_match_few" else "err_match_body",
                               args = list(length(missing_tips), shown(missing_tips)),
                               ntip = Ntip(tree))))
    }
    tree <- drop.tip(tree, missing_tips)
    log[[length(log) + 1]] <- .gui_log_entry("warn", "log_pruned", length(missing_tips))
  }

  aligned <- align_tree_data(tree, characters)
  rownames(aligned) <- NULL

  chars <- setdiff(colnames(aligned), "Species")
  poly <- chars[vapply(chars, function(ch) any(grepl("/", aligned[[ch]], fixed = TRUE)),
                       logical(1))]
  if (length(poly) > 0) {
    ex_val <- as.character(aligned[[poly[1]]][grepl("/", aligned[[poly[1]]], fixed = TRUE)][1])
    log[[length(log) + 1]] <- .gui_log_entry("info", "log_poly", length(poly), ex_val, poly[1])
  }
  log[[length(log) + 1]] <- .gui_log_entry("ok", "log_ready")
  list(ok = TRUE, tree = tree, characters = characters, aligned = aligned, log = log)
}


# ==============================================================================
# SECTION 2 -- SERVER
# ==============================================================================

#' Builds the server function
#' @param lang0    Initial language.
#' @param initial  NULL or list(tree, characters, kind, names, log).
#' @noRd
.gui_server <- function(lang0, initial) {
  function(input, output, session) {
    tags <- shiny::tags
    isolate <- shiny::isolate

    # ---- Language -------------------------------------------------------------
    dict <- .gui_dict()
    session$sendCustomMessage("mm-dict", list(dict = lapply(dict, as.list), lang = lang0))
    lang <- shiny::reactiveVal(lang0)
    shiny::observeEvent(input$mmr_lang, {
      if (input$mmr_lang %in% c("en", "es")) lang(input$mmr_lang)
    })
    trr <- shiny::reactive(.gui_translator(lang()))

    shiny::observeEvent(lang(), {
      tr <- trr()
      for (id in c("tree_format", "matrix_format", "csv_header", "ladderize", "stats_sort"))
        shiny::updateSelectInput(session, id, choices = .gui_select_choices(tr, id),
                                 selected = isolate(input[[id]]))
    }, ignoreInit = TRUE)

    # ---- State ------------------------------------------------------------------
    dat        <- shiny::reactiveVal(NULL)   # loaded data set
    data_state <- shiny::reactiveVal("empty")
    data_err   <- shiny::reactiveVal(NULL)
    data_log   <- shiny::reactiveVal(list())
    map_log    <- shiny::reactiveVal(list())
    src        <- shiny::reactiveValues(tree = NULL, matrix = NULL, kind = "files")
    data_id    <- 0L

    push_log <- function(rv, entries) {
      if (length(entries) == 0) return(invisible())
      rv(utils::tail(c(isolate(rv()), entries), 250))
    }

    output$data_state <- shiny::renderText(data_state())
    shiny::outputOptions(output, "data_state", suspendWhenHidden = FALSE)

    # ---- Data loading ---------------------------------------------------------------
    set_data <- function(tree, characters, prune, entries, kind, names) {
      rec <- .gui_reconcile(tree, characters, prune)
      entries <- c(entries, rec$log)
      if (!rec$ok) {
        data_log(entries)
        data_err(rec$error)
        dat(NULL)
        data_state("error")
        return(invisible(FALSE))
      }
      data_id <<- data_id + 1L
      chars  <- setdiff(colnames(rec$aligned), "Species")
      states <- lapply(chars, function(ch) sort_states(as.character(unique(rec$aligned[[ch]]))))
      names(states) <- chars
      data_log(entries)
      data_err(NULL)
      dat(list(tree = rec$tree, characters = rec$characters, aligned = rec$aligned,
               chars = chars, states = states, stats = character_stats(rec$aligned),
               id = data_id, has_lengths = !is.null(rec$tree$edge.length),
               kind = kind, names = names))
      data_state("loaded")
      invisible(TRUE)
    }

    read_opts <- function() {
      sp_col <- trimws(input$species_col %||% "1")
      if (!nzchar(sp_col)) sp_col <- "1"
      list(tree_format      = input$tree_format %||% "auto",
           matrix_format    = input$matrix_format %||% "auto",
           sep              = input$csv_sep %||% "auto",
           species_col      = if (grepl("^[0-9]+$", sp_col)) as.integer(sp_col) else sp_col,
           header           = switch(input$csv_header %||% "auto", yes = TRUE, no = FALSE, NA),
           normalize_spaces = isTRUE(input$normalize_spaces),
           prune            = isTRUE(input$prune_tips))
    }

    load_from_files <- function() {
      if (is.null(src$tree) || is.null(src$matrix)) {
        data_log(list(.gui_log_entry("warn", "need_files")))
        if (is.null(isolate(dat()))) data_state("empty")
        shiny::showNotification(isolate(trr())("need_files"), type = "warning")
        return(invisible())
      }
      opts <- read_opts()
      res <- .gui_read_inputs(src$tree$path, src$matrix$path,
                              c(src$tree$name, src$matrix$name), opts)
      if (!res$ok) {
        data_log(res$log)
        data_err(res$error)
        dat(NULL)
        data_state("error")
        return(invisible())
      }
      set_data(res$tree, res$characters, opts$prune, res$log, src$kind,
               c(src$tree$name, src$matrix$name))
    }

    label_file <- function(id, text) session$sendCustomMessage("mm-file-label", list(id = id, text = text))

    shiny::observeEvent(input$tree_file, {
      src$tree <- list(path = input$tree_file$datapath, name = input$tree_file$name)
      src$kind <- "files"
    })
    shiny::observeEvent(input$matrix_file, {
      src$matrix <- list(path = input$matrix_file$datapath, name = input$matrix_file$name)
      src$kind <- "files"
    })
    shiny::observeEvent(input$load_btn, load_from_files())
    shiny::observeEvent(input$retry_btn, load_from_files())

    shiny::observeEvent(input$example_btn, {
      tp <- system.file("extdata", "bats_tre.tre", package = "MultiMapR")
      mp <- system.file("extdata", "bats_matrix.csv", package = "MultiMapR")
      if (!nzchar(tp) || !nzchar(mp)) {
        shiny::showNotification(trr()("example_missing"), type = "error")
        return()
      }
      src$tree   <- list(path = tp, name = basename(tp))
      src$matrix <- list(path = mp, name = basename(mp))
      src$kind   <- "example"
      label_file("tree_file", basename(tp))
      label_file("matrix_file", basename(mp))
      load_from_files()
    })

    shiny::observeEvent(input$clear_btn, {
      src$tree <- NULL; src$matrix <- NULL; src$kind <- "files"
      dat(NULL); data_err(NULL); data_log(list()); map_log(list())
      data_state("empty")
      label_file("tree_file", "")
      label_file("matrix_file", "")
    })

    if (!is.null(initial)) {
      entries <- .gui_log_from_capture(initial$log)
      entries[[length(entries) + 1]] <- .gui_log_entry(
        "info", if (!is.null(initial$tree$edge.length)) "log_tree_read_bl" else "log_tree_read",
        initial$names[1], Ntip(initial$tree))
      entries[[length(entries) + 1]] <- .gui_log_entry(
        "info", "log_matrix_read", initial$names[2], "R",
        nrow(initial$characters), ncol(initial$characters) - 1L)
      isolate(set_data(initial$tree, initial$characters, prune = FALSE, entries,
                       kind = "preloaded", names = initial$names))
      shiny::observe({
        tr <- trr()
        d <- isolate(dat())
        if (is.null(d) || !identical(d$kind, "preloaded")) return()
        if (identical(initial$kind, "objects")) {
          label_file("tree_file", tr("obj_phylo"))
          label_file("matrix_file", tr("obj_df"))
        } else {
          label_file("tree_file", initial$names[1])
          label_file("matrix_file", initial$names[2])
        }
      })
    }

    # ---- Navigation -------------------------------------------------------------
    shiny::observeEvent(input$goto_map,    bslib::nav_select("nav", "map"))
    shiny::observeEvent(input$formats_btn, bslib::nav_select("nav", "help"))
    shiny::observeEvent(input$goto_data,   bslib::nav_select("nav", "data"))
    shiny::observeEvent(input$open_opts,   bslib::accordion_panel_open("read_opts_acc", "opts"))

    # ---- Data page outputs ---------------------------------------------------------
    output$src_badge <- shiny::renderUI({
      tr <- trr()
      d <- dat()
      key <- if (!is.null(d) && identical(d$kind, "preloaded")) "badge_preloaded"
             else if (!is.null(d) && identical(d$kind, "example")) "badge_example"
             else "badge_local"
      tags$span(class = "mm-badge", tr(key))
    })

    output$data_error <- shiny::renderUI({
      tr <- trr()
      e <- data_err()
      shiny::req(e)
      body <- if (!is.null(e$message_key)) do.call(tr, c(list(e$message_key), e$args)) else e$message
      tags$div(class = "mm-alert", role = "alert",
               tags$span(class = "mm-alert-icon", `aria-hidden` = "true", "!"),
               tags$div(
                 tags$h4(tr(e$title_key)),
                 tags$p(body, " ", tr("err_hint"),
                        if (!is.null(e$ntip)) paste0(" ", tr("err_tree_ok", e$ntip))),
                 tags$div(class = "mm-actions",
                          shiny::actionButton("open_opts", tr("open_opts"), class = "btn-mm-secondary"),
                          shiny::actionButton("retry_btn", tr("retry"), class = "btn-mm-primary"))))
    })

    output$data_log_err  <- shiny::renderUI(.gui_render_log(data_log(), trr()))
    output$data_log_ok   <- shiny::renderUI(.gui_render_log(data_log(), trr()))
    output$data_log_err_meta <- shiny::renderText(.gui_log_meta(data_log(), trr()))
    output$data_log_ok_meta  <- shiny::renderText(.gui_log_meta(data_log(), trr()))

    output$kpis <- shiny::renderUI({
      tr <- trr()
      d <- dat(); shiny::req(d)
      tot  <- d$stats[d$stats$character == "TOTAL", ]
      body <- d$stats[d$stats$character != "TOTAL", ]
      st <- unique(unlist(lapply(d$states, function(s) setdiff(s, c("?", "-")))))
      st <- unique(unlist(strsplit(st, "/", fixed = TRUE)))
      chars_note <- if (length(st) && all(grepl("^[0-9]+$", st))) {
        rng <- range(as.integer(st))
        if (max(body$n_states, na.rm = TRUE) <= 2) tr("kpi_binary", rng[1], rng[2])
        else tr("kpi_multistate", rng[1], rng[2])
      } else tr("kpi_states_max", max(body$n_states, na.rm = TRUE))
      kpi <- function(cls, label, value, note) {
        tags$div(class = paste("mm-kpi", cls),
                 tags$div(class = "mm-kpi-label", label),
                 tags$div(class = "mm-kpi-value", value),
                 tags$div(class = "mm-kpi-note", note))
      }
      tags$div(class = "mm-kpis",
               kpi("", tr("vb_tips"), Ntip(d$tree),
                   if (d$has_lengths) tr("vb_lengths") else tr("vb_nolengths")),
               kpi("", tr("vb_chars"), length(d$chars), chars_note),
               kpi("mm-kpi-missing", tr("vb_missing"), sprintf("%.1f%%", tot$pct_missing),
                   tr("n_cells", tot$n_missing)),
               kpi("mm-kpi-inapp", tr("vb_inapp"), sprintf("%.1f%%", tot$pct_inapplicable),
                   tr("n_cells", tot$n_inapplicable)))
    })

    output$stats_rows <- shiny::renderText({
      d <- dat(); shiny::req(d)
      trr()("n_rows", length(d$chars))
    })

    output$stats_table <- shiny::renderUI({
      tr <- trr()
      d <- dat(); shiny::req(d)
      st   <- d$stats
      body <- st[st$character != "TOTAL", ]
      tot  <- st[st$character == "TOTAL", ]
      cnt <- function(n, p) sprintf("%d \u00B7 %.0f%%", n, p)
      th <- function(i, key, num = FALSE, dot = NULL) {
        tags$th(scope = "col", `data-col` = i - 1L, `data-type` = if (num) "num" else "txt",
                class = if (num) "mm-num", `aria-sort` = "none",
                tags$button(type = "button",
                            if (!is.null(dot)) tags$span(class = paste("mm-dot", dot), `aria-hidden` = "true"),
                            tr(key), tags$span(class = "mm-sort", "\u25B2")))
      }
      # Stacked scored / missing / inapplicable bar, colored by the user's picks (CSS vars)
      dist <- function(r) {
        seg <- function(cls, n) {
          if (n > 0) tags$span(class = cls, style = sprintf("width:%.2f%%", 100 * n / r$n_taxa))
        }
        tags$span(class = "mm-dist", `aria-hidden` = "true",
                  seg("mm-c-s", r$n_scored), seg("mm-c-m", r$n_missing),
                  seg("mm-c-i", r$n_inapplicable))
      }
      rows <- lapply(seq_len(nrow(body)), function(i) {
        r <- body[i, ]
        tags$tr(tags$td(class = "mm-id", `data-v` = r$character, r$character),
                tags$td(class = "mm-num", `data-v` = r$n_states, r$n_states),
                tags$td(class = "mm-num", `data-v` = r$n_scored, r$n_scored),
                tags$td(class = "mm-num", `data-v` = r$n_missing, cnt(r$n_missing, r$pct_missing)),
                tags$td(class = "mm-num", `data-v` = r$n_inapplicable,
                        cnt(r$n_inapplicable, r$pct_inapplicable)),
                tags$td(class = "mm-dist-cell", dist(r)))
      })
      tags$table(class = "mm-table",
                 tags$thead(tags$tr(th(1, "col_character"), th(2, "col_states", TRUE),
                                    th(3, "col_scored", TRUE, "mm-c-s"),
                                    th(4, "col_missing", TRUE, "mm-c-m"),
                                    th(5, "col_inapp", TRUE, "mm-c-i"),
                                    tags$th(scope = "col", tr("col_dist")))),
                 tags$tbody(rows),
                 tags$tfoot(tags$tr(tags$td("TOTAL"), tags$td(class = "mm-num", ""),
                                    tags$td(class = "mm-num", tot$n_scored),
                                    tags$td(class = "mm-num", cnt(tot$n_missing, tot$pct_missing)),
                                    tags$td(class = "mm-num",
                                            cnt(tot$n_inapplicable, tot$pct_inapplicable)),
                                    tags$td(class = "mm-dist-cell", dist(tot)))))
    })

    # ---- Data-category colors (scored / missing / inapplicable) and borders -------
    stats_colors <- shiny::reactive(.resolve_stats_colors(c(
      scored       = input$stats_col_scored,
      missing      = input$stats_col_missing,
      inapplicable = input$stats_col_inapp)))
    stats_border <- shiny::reactive(if (isFALSE(input$stats_borders)) NA else "grey30")
    heat_values  <- shiny::reactive(if (isFALSE(input$heat_values)) FALSE else NULL)

    # The HTML views (matrix, table bars, KPIs) take the colors from CSS variables,
    # so changing a color does not re-render them.
    shiny::observe({
      cols <- .gui_hex(stats_colors())
      ink  <- .contrast_ink(cols)
      session$sendCustomMessage("mm-stats-colors", list(
        scored = cols[[1]], missing = cols[[2]], inapplicable = cols[[3]],
        scored_ink = ink[[1]], missing_ink = ink[[2]], inapplicable_ink = ink[[3]],
        borders = !isFALSE(input$stats_borders)))
    })
    shiny::observeEvent(input$stats_col_reset, {
      ids <- c(scored = "stats_col_scored", missing = "stats_col_missing",
               inapplicable = "stats_col_inapp")
      for (k in names(ids))
        session$sendInputMessage(ids[[k]], list(value = .gui_hex(CHAR_STATS_COLORS[[k]])))
    })

    stats_height <- function() {
      d <- dat(); if (is.null(d)) return(300)
      max(320, round((0.28 * length(d$chars) + 1.5) * 80))
    }
    heat_layout <- shiny::reactive({
      d <- dat(); shiny::req(d)
      .completeness_layout(as.character(d$aligned$Species), d$chars, cell_in = 0.22)
    })
    # Wide matrices get a plot wider than the box (it scrolls) so cells stay legible;
    # a fixed-width container stops Shiny from squeezing the image to 100%.
    output$heat_box <- shiny::renderUI({
      shiny::plotOutput("heat_plot", height = "auto",
                        width = paste0(max(600, round(heat_layout()$width * 96)), "px"))
    })
    heat_height <- function() {
      if (is.null(dat())) return(300)
      max(360, round(heat_layout()$height * 96))
    }
    output$stats_plot <- shiny::renderPlot({
      d <- dat(); shiny::req(d)
      plot_character_stats(d$aligned, sort_by = input$stats_sort %||% "none",
                           colors = stats_colors(), border = stats_border())
    }, height = stats_height, res = 96, bg = "white")
    output$heat_plot <- shiny::renderPlot({
      d <- dat(); shiny::req(d)
      plot_character_completeness(d$aligned, colors = stats_colors(),
                                  border = stats_border(), show_values = heat_values())
    }, height = heat_height, res = 96, bg = "white")

    # Interactive taxon x character matrix (hover highlights row + column, see multimapr.js)
    output$matrix_table <- shiny::renderUI({
      tr <- trr()
      d <- dat(); shiny::req(d)
      esc     <- htmltools::htmlEscape
      species <- as.character(d$aligned$Species)
      n_show  <- min(length(d$chars), max(1L, floor(60000 / max(1L, length(species)))))
      chars   <- d$chars[seq_len(n_show)]
      cells <- vapply(chars, function(ch) {
        v <- trimws(as.character(d$aligned[[ch]]))
        v[is.na(v) | v == ""] <- "?"
        cls <- c(scored = "mm-c-s", missing = "mm-c-m",
                 inapplicable = "mm-c-i")[.classify_character_values(v)]
        paste0('<td class="', cls, '">', esc(v), "</td>")
      }, character(length(species)))
      cells <- matrix(cells, nrow = length(species))
      body <- paste0('<tr><th scope="row">', esc(species), "</th>",
                     apply(cells, 1, paste, collapse = ""), "</tr>", collapse = "")
      head <- paste0('<tr><th scope="col" class="mm-mx-corner">', esc(tr("col_taxon")), "</th>",
                     paste0('<th scope="col"><span>', esc(chars), "</span></th>", collapse = ""),
                     "</tr>")
      note <- if (n_show < length(d$chars))
        tags$div(class = "mm-note mm-mx-note", tr("mx_trunc", n_show, length(d$chars)))
      htmltools::tagList(note,
              shiny::HTML(paste0('<table class="mm-matrix"><thead>', head, "</thead><tbody>",
                                 body, "</tbody></table>")))
    })

    output$dl_stats_csv <- shiny::downloadHandler(
      filename = function() "MultiMapR_character_stats.csv",
      content  = function(file) {
        d <- dat(); shiny::req(d)
        utils::write.csv(d$stats, file, row.names = FALSE)
      })

    completeness_download <- function(format) {
      shiny::downloadHandler(
        filename = function() {
          kind <- if (identical(input$cmp_tabs, "per")) "character_stats" else "completeness"
          paste0("MultiMapR_", kind, ".", format)
        },
        content = function(file) {
          d <- dat(); shiny::req(d)
          base <- tempfile("mmr_stats_")
          .gui_with_sink_device(utils::capture.output(
            if (identical(input$cmp_tabs, "per")) {
              plot_character_stats(d$aligned, sort_by = input$stats_sort %||% "none",
                                   colors = stats_colors(), border = stats_border(),
                                   export_filename = base, export_format = format)
            } else {
              # "matrix" and "heat" tabs both export the heatmap (with states in cells
              # unless switched off on the heatmap tab)
              plot_character_completeness(d$aligned, colors = stats_colors(),
                                          border = stats_border(), show_values = heat_values(),
                                          export_filename = base, export_format = format)
            }))
          file.copy(paste0(base, ".", format), file, overwrite = TRUE)
        })
    }
    output$dl_cmp_png <- completeness_download("png")
    output$dl_cmp_pdf <- completeness_download("pdf")

    map <- .gui_server_map(input, output, session, dat, trr, map_log, push_log)
  }
}


# ==============================================================================
# SECTION 3 -- MAP PAGE
# ==============================================================================

#' Server logic of "2 - Map characters"
#' @noRd
.gui_server_map <- function(input, output, session, dat, trr, map_log, push_log) {
  tags <- shiny::tags
  isolate <- shiny::isolate
  session_start <- format(Sys.time(), "%H:%M")

  mode_num <- shiny::reactive(if (identical(input$mode, "ancestral")) 2L else 1L)

  # ---- Characters (selection order = drawing order) -----------------------------
  chars_sel <- shiny::reactiveVal(character(0))
  shiny::observeEvent(dat(), {
    d <- dat()
    chars_sel(if (is.null(d)) character(0) else d$chars[1])
  }, ignoreNULL = FALSE)

  shiny::observeEvent(input$char_toggle, {
    ch <- input$char_toggle
    d <- dat()
    if (is.null(d) || !ch %in% d$chars) return()
    cur <- chars_sel()
    chars_sel(if (ch %in% cur) setdiff(cur, ch) else c(cur, ch))
  })
  shiny::observe({
    session$sendCustomMessage("mm-chars", list(sel = I(chars_sel())))
  })

  chars_used <- shiny::reactive({
    sel <- chars_sel()
    if (mode_num() == 2L) utils::head(sel, 3L) else sel
  })

  output$char_list <- shiny::renderUI({
    tr <- trr()
    d <- dat()
    if (is.null(d)) return(tags$div(class = "mm-note", tr("need_data")))
    body <- d$stats[d$stats$character != "TOTAL", ]
    sel <- isolate(chars_sel())
    tags$div(class = "mm-charlist", role = "group", `aria-label` = tr("step_chars"),
             tags$div(class = "mm-charlist-head", `aria-hidden` = "true",
                      tags$span(tr("col_character")), tags$span(tr("col_states")),
                      tags$span("? %")),
             lapply(seq_len(nrow(body)), function(i) {
               ch <- body$character[i]
               on <- ch %in% sel
               tags$button(type = "button",
                           class = paste("mm-charrow", if (on) "is-selected"),
                           `data-char` = ch, `aria-pressed` = if (on) "true" else "false",
                           tags$span(class = "mm-cname",
                                     tags$span(class = "mm-box", `aria-hidden` = "true"), ch),
                           tags$span(class = "mm-cnum", body$n_states[i]),
                           tags$span(class = "mm-cnum", sprintf("%.0f%%", body$pct_missing[i])))
             }))
  })

  output$chips <- shiny::renderUI({
    tr <- trr()
    sel <- chars_sel()
    if (length(sel) == 0) return(tags$div(class = "mm-note mb-2", tr("no_chars")))
    used <- chars_used()
    tags$div(class = "mm-chips",
             lapply(sel, function(ch) {
               tags$span(class = paste("mm-chip", if (!ch %in% used) "mm-chip-off"),
                         ch,
                         tags$button(type = "button", class = "mm-chip-x", `data-char` = ch,
                                     `aria-label` = tr("remove_char", ch), title = tr("remove_char", ch),
                                     "\u00D7"))
             }))
  })

  output$trim_warn <- shiny::renderUI({
    tr <- trr()
    sel <- chars_sel()
    if (mode_num() != 2L || length(sel) <= 3L) return(NULL)
    tags$div(class = "mm-warnbox", role = "status",
             tags$b(tr("lvl_warn")), " ",
             tr("trim_warn", length(sel), paste(utils::head(sel, 3L), collapse = ", ")))
  })

  # ---- States and colors ------------------------------------------------------------
  card_prefix <- function(d, ch) paste0("d", d$id, "_c", match(ch, d$chars))
  card_states <- function(d, ch) setdiff(d$states[[ch]], c("?", "-"))
  card_choice <- function(d, ch) {
    pal <- input[[paste0(card_prefix(d, ch), "_pal")]] %||% "global"
    if (identical(pal, "global")) input$pal_global %||% "auto" else pal
  }

  output$char_cards <- shiny::renderUI({
    tr <- trr()
    d <- dat()
    if (is.null(d)) return(NULL)
    sel <- chars_used()
    if (length(sel) == 0) return(NULL)
    body <- d$stats[d$stats$character != "TOTAL", ]
    mt <- isolate(mode_num())

    isolate(lapply(seq_along(sel), function(slot) {
      ch     <- sel[slot]
      pfx    <- card_prefix(d, ch)
      states <- card_states(d, ch)
      vals   <- as.character(d$aligned[[ch]])
      srow   <- body[body$character == ch, ]
      defaults <- .gui_palette_colors(card_choice(d, ch), length(states), slot, mt, length(sel))

      rows <- lapply(seq_along(states), function(si) {
        inc_id <- paste0(pfx, "_inc_", si)
        col_id <- paste0(pfx, "_col_", si)
        inc <- input[[inc_id]] %||% TRUE
        tags$div(class = "mm-state",
                 tags$input(type = "checkbox", id = inc_id, class = "mm-state-inc",
                            checked = if (isTRUE(inc)) NA, `aria-label` = tr("include_state", states[si])),
                 .gui_color_input(col_id, input[[col_id]] %||% defaults[si],
                                  label = tr("state_color", states[si])),
                 tags$span(class = "mm-state-label", states[si]),
                 tags$span(class = "mm-state-count", tr("n_tax", sum(vals == states[si], na.rm = TRUE))))
      })
      special <- lapply(c("?", "-"), function(s) {
        tags$div(class = "mm-state mm-special",
                 tags$input(type = "checkbox", disabled = NA, checked = NA,
                            `aria-label` = tr(if (s == "?") "state_missing" else "state_inapp")),
                 tags$span(class = "mm-state-swatch", style = "background:#b3b3b3", `aria-hidden` = "true"),
                 tags$span(class = "mm-state-label",
                           paste(s, tr(if (s == "?") "state_missing" else "state_inapp"))),
                 tags$span(class = "mm-state-count", tr("n_tax", sum(vals == s, na.rm = TRUE))))
      })

      tags$div(class = "mm-ccard",
               tags$div(class = "mm-ccard-head",
                        tags$span(class = "mm-ccard-id", ch),
                        tags$span(class = "mm-ccard-meta",
                                  tr("card_meta", srow$n_scored, length(states)))),
               tags$div(class = "mm-ccard-pal",
                        tags$span(class = "mm-upper", tr("palette")),
                        shiny::selectInput(paste0(pfx, "_pal"), NULL,
                                           .gui_select_choices(tr, "card_palette"),
                                           selected = input[[paste0(pfx, "_pal")]] %||% "global",
                                           selectize = FALSE, width = "100%")),
               tags$div(class = "mm-ccard-body", rows, special),
               tags$div(class = "mm-ccard-err", role = "alert", tr("no_states", ch)))
    }))
  })

  recolor <- function(d, ch, choice) {
    sel <- chars_used()
    slot <- match(ch, sel, nomatch = 1L)
    states <- card_states(d, ch)
    cols <- .gui_palette_colors(choice, length(states), slot, mode_num(), length(sel))
    pfx <- card_prefix(d, ch)
    for (si in seq_along(states))
      session$sendInputMessage(paste0(pfx, "_col_", si), list(value = cols[si]))
  }

  shiny::observeEvent(input$pal_global, {
    d <- dat(); if (is.null(d)) return()
    for (ch in chars_used()) {
      if (identical(input[[paste0(card_prefix(d, ch), "_pal")]] %||% "global", "global"))
        recolor(d, ch, input$pal_global)
    }
  }, ignoreInit = TRUE)

  applied_palette <- new.env(parent = emptyenv())
  palette_observers <- list()
  shiny::observeEvent(dat(), {
    for (o in palette_observers) o$destroy()
    palette_observers <<- list()
    d <- dat(); if (is.null(d)) return()
    palette_observers <<- lapply(d$chars, function(ch) {
      pal_id <- paste0(card_prefix(d, ch), "_pal")
      shiny::observeEvent(input[[pal_id]], {
        choice <- input[[pal_id]]
        if (identical(applied_palette[[pal_id]], choice)) return()
        first_time <- is.null(applied_palette[[pal_id]])
        assign(pal_id, choice, envir = applied_palette)
        if (first_time) return()   # colors already initialised by the card
        recolor(d, ch, if (identical(choice, "global")) input$pal_global %||% "auto" else choice)
      }, ignoreInit = TRUE)
    })
  }, ignoreNULL = FALSE)

  # ---- Tree options -------------------------------------------------------------------
  bl_pref <- shiny::reactiveVal(TRUE)
  bl_disabled <- shiny::reactive({
    d <- dat()
    is.null(d) || !d$has_lengths || identical(input$topology, "cladogram")
  })
  shiny::observeEvent(input$use_bl, {
    if (!isolate(bl_disabled())) bl_pref(isTRUE(input$use_bl))
  })
  output$use_bl_ui <- shiny::renderUI({
    tr <- trr()
    d <- dat()
    off <- bl_disabled()
    note <- if (is.null(d) || !d$has_lengths) tr("bl_none")
            else if (identical(input$topology, "cladogram")) tr("bl_clado")
            else tr("bl_has")
    cb <- shiny::checkboxInput("use_bl", tr("use_bl"), value = !off && isolate(bl_pref()))
    if (off) cb <- htmltools::tagQuery(cb)$find("input")$addAttrs(disabled = NA)$allTags()
    tags$div(class = paste("mb-3", if (off) "mm-disabled"),
             cb, tags$div(class = "mm-note", style = "margin-top:-6px", note))
  })
  uses_lengths <- shiny::reactive(!bl_disabled() && isTRUE(input$use_bl %||% bl_pref()))

  output$tip_note <- shiny::renderUI({
    tr <- trr()
    if (mode_num() == 1L)
      return(tags$div(class = "mm-note", tr("tip_simple")))
    if (uses_lengths())
      return(tags$div(class = "mm-note mm-note-warn", tr("tip_bl")))
    NULL
  })

  # ---- Configuration (same structure as setup_mapping_config()) ------------------------
  spec <- shiny::reactive({
    d <- dat()
    if (is.null(d)) return(list(status = "empty"))
    sel <- chars_used()
    if (length(sel) == 0) return(list(status = "nochars"))
    mt <- mode_num()

    colors_by_char <- list()
    empty_chars    <- character(0)
    all_selected   <- TRUE
    for (slot in seq_along(sel)) {
      ch     <- sel[slot]
      pfx    <- card_prefix(d, ch)
      states <- card_states(d, ch)
      defaults <- .gui_palette_colors(isolate(card_choice(d, ch)), length(states), slot, mt,
                                      length(sel))
      cols <- vapply(seq_along(states), function(si)
        input[[paste0(pfx, "_col_", si)]] %||% defaults[si], character(1))
      inc  <- vapply(seq_along(states), function(si)
        isTRUE(input[[paste0(pfx, "_inc_", si)]] %||% TRUE), logical(1))
      if (length(states) > 0 && !any(inc)) {
        empty_chars <- c(empty_chars, ch)
        next
      }
      all_selected <- all_selected && all(inc)
      colors_by_char[[ch]] <- stats::setNames(cols[inc], states[inc])
    }
    draw <- setdiff(sel, empty_chars)
    if (length(draw) == 0) return(list(status = "nostates", empty = empty_chars))

    tree_type <- input$topology %||% "phylogram"
    pch <- as.integer(input$shape %||% "21")
    size <- input$shape_size %||% 1.4
    config <- list(mapping_type         = mt,
                   aligned_data         = d$aligned,
                   datos_ord            = d$aligned,
                   ambiguity_color      = NULL,
                   caracteres           = draw,
                   colores_por_caracter = colors_by_char,
                   tipo_arbol           = tree_type,
                   grosor               = input$branch_width %||% 2,
                   use_edge_length      = uses_lengths(),
                   legend_corner        = input$legend_corner %||% "bottomleft",
                   export_filename      = NULL,
                   export_format        = NULL,
                   width                = NULL,
                   height               = NULL)
    if (mt == 1L) {
      shapes <- !identical(input$presentation, "labels")
      config$simple_mode <- if (shapes) "figures" else "tip_color"
      config$pch_figura  <- if (shapes) pch else 21L
      config$tam_figura  <- if (shapes) size else 1
    } else {
      config$funcion_multi <- if (identical(input$viz, "tree_tips")) 2L else 1L
      config$mapear_todos  <- all_selected
      config$algoritmo     <- if (identical(input$algo, "fitch")) 2L else 1L
      config$fitch_mode    <- if (config$algoritmo == 2L) input$optim %||% "acctran" else NULL
      if (config$algoritmo == 2L) config$ambiguity_color <- input$ambiguity_color %||% "#FF00FF"
      if (config$funcion_multi == 2L) {
        config$pch_figura <- pch
        config$tam_figura <- size
      }
    }
    ladder <- switch(input$ladderize %||% "TRUE", "TRUE" = TRUE, "right" = "right", FALSE)
    list(status = "ok",
         tree   = d$tree,
         config = config,
         empty  = empty_chars,
         trimmed = if (mt == 2L && length(chars_sel()) > 3L) utils::head(chars_sel(), 3L),
         render = list(branch_width     = input$branch_width %||% 2,
                       use_edge_length  = uses_lengths(),
                       ladderize        = ladder,
                       terminal_stretch = if (mt == 2L && !uses_lengths())
                                            input$tip_mult %||% 1 else 1))
  })
  spec_d <- shiny::debounce(spec, 400)

  custom_dims <- shiny::reactive({
    if (!identical(input$dims, "custom")) return(list(w = NULL, h = NULL))
    w <- input$exp_w; h <- input$exp_h
    ok <- function(x) is.numeric(x) && length(x) == 1 && !is.na(x) && x > 0
    if (ok(w) && ok(h)) list(w = w, h = h) else list(w = NULL, h = NULL)
  })

  # ---- Summaries ------------------------------------------------------------------------
  topo_label <- function(tr, x) tr(x %||% "phylogram")
  n_chars_label <- function(tr, n) if (n == 1) tr("n_char_1") else tr("n_chars", n)
  optim_label <- function(tr) switch(input$optim %||% "acctran",
                                     acctran = "ACCTRAN", deltran = "DELTRAN", tr("unambiguous"))

  output$sum_A <- shiny::renderText({
    tr <- trr()
    if (mode_num() == 1L) return(tr("sum_simple"))
    if (identical(input$algo, "fitch")) paste0(tr("sum_ancestral"), " \u00B7 ", optim_label(tr))
    else paste0(tr("sum_ancestral"), " \u00B7 ", tr("sum_majority"))
  })
  output$sum_B <- shiny::renderText({
    n <- length(chars_used())
    if (n == 0) "\u2014" else n_chars_label(trr(), n)
  })
  output$sum_C <- shiny::renderText(trr()(.GUI_PALETTE_KEYS[[input$pal_global %||% "auto"]]))
  output$sum_D <- shiny::renderText(topo_label(trr(), input$topology))
  output$sum_E <- shiny::renderText({
    arrow <- switch(input$legend_corner %||% "bottomleft",
                    topleft = "\u2196", topright = "\u2197",
                    bottomleft = "\u2199", bottomright = "\u2198")
    paste0(arrow, " \u00B7 ", input$preview_height %||% 640, "px")
  })
  output$sum_F <- shiny::renderText(toupper(input$fmt %||% "png"))

  output$side_counter <- shiny::renderText({
    tr <- trr()
    d <- dat()
    if (is.null(d)) tr("no_data_short") else tr("side_counter", Ntip(d$tree), length(d$chars))
  })

  output$corner_label <- shiny::renderText(trr()(paste0("corner_", input$legend_corner %||% "bottomleft")))

  output$cfg_summary <- shiny::renderText({
    tr <- trr()
    n <- length(chars_used())
    parts <- if (mode_num() == 1L) {
      c(tr("sum_simple"),
        if (identical(input$presentation, "labels")) tr("sum_labels") else tr("sum_shapes"))
    } else if (identical(input$algo, "fitch")) {
      c("Fitch", optim_label(tr))
    } else {
      tr("sum_majority_long")
    }
    paste(c(parts, n_chars_label(tr, n), topo_label(tr, input$topology)), collapse = " \u00B7 ")
  })

  output$computing_detail <- shiny::renderText({
    tr <- trr()
    d <- dat(); shiny::req(d)
    algo <- if (mode_num() == 1L) tr("sum_simple")
            else if (identical(input$algo, "fitch")) "Fitch" else tr("sum_majority")
    paste(algo, n_chars_label(tr, length(chars_used())),
          tr("n_tips", Ntip(d$tree)), tr("n_nodes", d$tree$Nnode), sep = " \u00B7 ")
  })

  output$export_btn_label <- shiny::renderText(trr()("export_fmt", toupper(input$fmt %||% "png")))

  output$session_stamp <- shiny::renderText(trr()("session_at", session_start))

  # ---- Preview ---------------------------------------------------------------------------
  last_png <- NULL
  preview <- shiny::reactive({
    input$recalc
    input$recalc_err
    s <- spec_d()
    if (!identical(s$status, "ok")) return(s)
    wcss <- session$clientData$output_preview_img_width %||% 0
    if (!is.numeric(wcss) || wcss < 240) wcss <- 900
    wcss <- wcss - 4
    ratio <- min(max(session$clientData$pixelratio %||% 1, 1), 2)
    cd <- isolate(custom_dims())
    dims <- .gui_figure_dims(Ntip(s$tree), s$config$tipo_arbol, cd$w, cd$h)
    file <- tempfile("mmr_preview_", fileext = ".png")
    res <- tryCatch(.gui_render_png(s, file, wcss * ratio, dims), error = function(e) e)
    if (!is.null(last_png) && file.exists(last_png)) unlink(last_png)
    if (inherits(res, "error")) {
      last_png <<- NULL
      return(list(status = "error", message = conditionMessage(res), spec = s))
    }
    last_png <<- file
    list(status = "ok", path = file, spec = s, log = res$log, seconds = res$seconds,
         width = round(wcss), height = round(res$height_px / ratio))
  })

  output$preview_img <- shiny::renderImage({
    p <- preview()
    shiny::req(identical(p$status, "ok"))
    list(src = p$path, contentType = "image/png", width = p$width, height = p$height,
         alt = isolate(trr())("preview_alt"))
  }, deleteFile = FALSE)

  output$preview_state <- shiny::renderUI({
    tr <- trr()
    p <- preview()
    switch(p$status,
      empty = .gui_empty_state(tr, "empty_title", tr("map_empty_body"),
                               actions = shiny::actionButton("goto_data", tr("goto_data"),
                                                             class = "btn-mm-accent-outline")),
      nochars = .gui_empty_state(tr, "map_nochars_title", tr("map_nochars_body")),
      nostates = .gui_empty_state(tr, "map_nostates_title",
                                  tr("map_nostates_body", paste(p$empty, collapse = ", "))),
      error = tags$div(class = "mm-alert", role = "alert",
                       tags$span(class = "mm-alert-icon", `aria-hidden` = "true", "!"),
                       tags$div(tags$h4(tr("render_error_title")),
                                tags$p(p$message),
                                shiny::actionButton("recalc_err", tr("recalc"), class = "btn-mm-primary"))),
      NULL)
  })

  # Map log: one block of lines per regenerated preview
  shiny::observeEvent(preview(), {
    p <- preview()
    s <- p$spec
    if (!p$status %in% c("ok", "error", "nostates")) return()
    entries <- list()
    add <- function(...) entries[[length(entries) + 1]] <<- .gui_log_entry(...)
    if (identical(p$status, "nostates")) {
      for (ch in p$empty) add("error", "log_no_states", ch)
      push_log(map_log, entries)
      return()
    }
    cfg <- s$config
    add("info", "log_active", paste(cfg$caracteres, collapse = ", "))
    if (length(s$trimmed)) add("warn", "log_trimmed", length(isolate(chars_sel())),
                               paste(s$trimmed, collapse = ", "))
    for (ch in s$empty) add("error", "log_no_states", ch)
    if (identical(cfg$tipo_arbol, "cladogram") && isTRUE(isolate(dat())$has_lengths))
      add("warn", "log_clado_bl")
    if (identical(p$status, "error")) {
      add("error", text = p$message)
    } else {
      entries <- c(entries, .gui_log_from_capture(p$log))
      if (isTRUE(cfg$algoritmo == 2L)) {
        n <- .gui_count_ambiguous(s, cfg$fitch_mode)
        if (!is.na(n)) add("info", "log_ambig", toupper(cfg$fitch_mode), n,
                           toupper(.gui_hex(.get_fitch_ambig_color(cfg))))
      }
      add("ok", "log_rendered", sprintf("%.2f", p$seconds))
    }
    push_log(map_log, entries)
  })

  output$map_log <- shiny::renderUI(.gui_render_log(map_log(), trr()))

  # ---- ACCTRAN / DELTRAN comparison --------------------------------------------------------
  cmp_open <- shiny::reactiveVal(FALSE)
  shiny::observeEvent(input$cmp_toggle, cmp_open(!cmp_open()))
  shiny::observeEvent(input$cmp_close, cmp_open(FALSE))
  cmp_active <- shiny::reactive(cmp_open() && mode_num() == 2L && identical(input$algo, "fitch"))
  output$cmp_visible <- shiny::renderText(if (cmp_active()) "yes" else "no")
  shiny::outputOptions(output, "cmp_visible", suspendWhenHidden = FALSE)
  shiny::observe(session$sendCustomMessage("mm-class", list(id = "cmp_toggle", cls = "is-on",
                                                            on = cmp_open())))
  shiny::observeEvent(input$use_optim, {
    shiny::updateRadioButtons(session, "optim", selected = input$use_optim)
  })

  lapply(c("acctran", "deltran"), function(m) {
    output[[paste0("cmp_img_", m)]] <- shiny::renderImage({
      shiny::req(cmp_active())
      s <- spec_d()
      shiny::req(identical(s$status, "ok"), isTRUE(s$config$algoritmo == 2L))
      s$config$fitch_mode <- m
      dims <- .gui_figure_dims(Ntip(s$tree), s$config$tipo_arbol)
      file <- tempfile(paste0("mmr_cmp_", m, "_"), fileext = ".png")
      .gui_render_png(s, file, 900, dims)
      list(src = file, contentType = "image/png", alt = paste("Fitch", toupper(m)))
    }, deleteFile = TRUE)

    output[[paste0("cmp_note_", m)]] <- shiny::renderText({
      shiny::req(cmp_active())
      s <- spec_d()
      shiny::req(identical(s$status, "ok"), isTRUE(s$config$algoritmo == 2L))
      tr <- trr()
      n <- .gui_count_ambiguous(s, m)
      paste0(tr(paste0("cmp_note_", m)), " \u00B7 ",
             if (is.na(n)) "\u2014" else tr("n_ambig", n))
    })

    output[[paste0("cmp_use_", m)]] <- shiny::renderUI({
      tr <- trr()
      active <- identical(input$optim, m)
      tags$button(type = "button", class = paste("mm-cmp-use", if (active) "is-active"),
                  disabled = if (active) NA,
                  onclick = sprintf("Shiny.setInputValue('use_optim','%s',{priority:'event'})", m),
                  if (active) tr("in_use") else tr("use_x", toupper(m)))
    })
  })

  # ---- Export --------------------------------------------------------------------------------
  auto_name <- shiny::reactive({
    sel <- chars_used()
    topo <- tolower(topo_label(trr(), input$topology))
    if (length(sel) == 0) paste0("MultiMapR_", topo)
    else paste0("MultiMapR_", paste(sel, collapse = "-"), "_", topo)
  })
  name_state <- shiny::reactiveValues(edited = FALSE, last = NULL)
  shiny::observe({
    n <- auto_name()
    if (!isolate(name_state$edited)) {
      name_state$last <- n
      shiny::updateTextInput(session, "export_name", value = n)
    }
  })
  shiny::observeEvent(input$export_name, {
    v <- input$export_name
    if (!nzchar(trimws(v))) {
      name_state$edited <- FALSE
      name_state$last <- auto_name()
      shiny::updateTextInput(session, "export_name", value = name_state$last)
    } else if (!identical(v, name_state$last)) {
      name_state$edited <- TRUE
    }
  })

  output$dims_note <- shiny::renderText({
    tr <- trr()
    d <- dat()
    n <- if (is.null(d)) 0L else Ntip(d$tree)
    dims <- .gui_figure_dims(n, input$topology)
    if (identical(input$topology, "fan")) tr("dims_note_fan", dims$width, dims$height)
    else tr("dims_note", n, dims$height)
  })
  shiny::observeEvent(list(dat(), input$topology), {
    d <- dat(); if (is.null(d)) return()
    dims <- .gui_figure_dims(Ntip(d$tree), input$topology)
    shiny::updateNumericInput(session, "exp_w", value = dims$width)
    shiny::updateNumericInput(session, "exp_h", value = round(dims$height * 2) / 2)
  })

  exports_dir <- file.path(normalizePath(getwd(), winslash = "/"), "Exports")
  output$exports_path <- shiny::renderText(paste0(exports_dir, "/"))

  export_args <- function() {
    s <- spec()
    if (!identical(s$status, "ok")) stop(trr()("export_nothing"), call. = FALSE)
    cd <- custom_dims()
    list(spec = s, format = input$fmt %||% "png",
         filename = .gui_clean_filename(input$export_name, isolate(auto_name())),
         width = cd$w, height = cd$h)
  }

  map_download <- function() {
    shiny::downloadHandler(
      filename = function() {
        a <- tryCatch(export_args(), error = function(e) NULL)
        if (is.null(a)) "MultiMapR.png" else paste0(a$filename, ".", a$format)
      },
      content = function(file) {
        a <- export_args()
        tmp <- tempfile("mmr_export_")
        dir.create(tmp)
        res <- tryCatch(.gui_export_mapping(a$spec, a$format, a$filename, a$width, a$height,
                                            out_dir = tmp),
                        error = function(e) e)
        if (inherits(res, "error")) {
          push_log(map_log, list(.gui_log_entry("error", "log_export_error", conditionMessage(res))))
          stop(conditionMessage(res), call. = FALSE)
        }
        push_log(map_log, list(.gui_log_entry("ok", "log_downloaded",
                                              paste0(a$filename, ".", a$format))))
        file.copy(res$path, file, overwrite = TRUE)
      })
  }
  output$dl_map     <- map_download()
  output$dl_map_top <- map_download()

  saved <- shiny::reactiveVal(NULL)
  shiny::observeEvent(input$save_map, {
    a <- tryCatch(export_args(), error = function(e) e)
    if (inherits(a, "error")) {
      push_log(map_log, list(.gui_log_entry("error", text = conditionMessage(a))))
      return()
    }
    res <- tryCatch(.gui_export_mapping(a$spec, a$format, a$filename, a$width, a$height,
                                        out_dir = getwd()),
                    error = function(e) e)
    if (inherits(res, "error")) {
      push_log(map_log, list(.gui_log_entry("error", "log_export_error", conditionMessage(res))))
      saved(NULL)
      return()
    }
    rel <- paste0("Exports/", basename(res$path))
    push_log(map_log, list(.gui_log_entry("ok", "saved_to", rel)))
    saved(rel)
  })
  output$save_ok <- shiny::renderUI({
    rel <- saved(); shiny::req(rel)
    tags$div(class = "mm-okbox", role = "status", tags$b(trr()("lvl_ok")), " ",
             trr()("saved_to", rel))
  })

  # Returned for tests (shiny::testServer)
  list(preview = preview, spec = spec, chars_used = chars_used, map_log = map_log,
       export_args = export_args)
}
