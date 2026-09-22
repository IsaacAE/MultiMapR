################################################################################
# gui_ui.R
#
# Layout and components of the MultiMapR interface (see gui.R).
# Every visible string comes from the dictionary (inst/app/i18n.txt): static
# text is emitted as <span data-i18n="key"> so the language can be switched
# in the browser without reloading; server-rendered fragments are translated
# in R with the current language.
################################################################################


# ==============================================================================
# SECTION 1 -- SMALL COMPONENTS
# ==============================================================================

#' Translatable text node
#' @param key Dictionary key.
#' @param tr  Translator.
#' @param tag Tag function (default span).
#' @noRd
.gui_t <- function(key, tr, tag = shiny::tags$span, ...) {
  tag(`data-i18n` = key, ..., tr(key))
}

#' Title + one-line explanation used as a radio "card" choice
#' @noRd
.gui_rcard <- function(tr, title_key, hint_key = NULL) {
  tags <- shiny::tags
  tags$span(class = "mm-rcard",
            .gui_t(title_key, tr, class = "mm-rcard-title"),
            if (!is.null(hint_key)) .gui_t(hint_key, tr, class = "mm-rcard-hint"))
}

#' Group label
#' @noRd
.gui_label <- function(key, tr) .gui_t(key, tr, class = "mm-label")

#' Inline SVG icons (topologies, empty state, logo)
#' @noRd
.gui_icon <- function(name) {
  svg <- switch(name,
    phylogram = '<svg width="20" height="18" viewBox="0 0 20 18" aria-hidden="true"><path d="M2 9 H6 M6 3 V15 M6 3 H13 M6 9 H11 M6 15 H16" stroke-width="1.6" fill="none" stroke-linecap="round"/></svg>',
    cladogram = '<svg width="20" height="18" viewBox="0 0 20 18" aria-hidden="true"><path d="M2 9 L10 3 H17 M2 9 H17 M2 9 L10 15 H17" stroke-width="1.6" fill="none" stroke-linecap="round" stroke-linejoin="round"/></svg>',
    fan       = '<svg width="20" height="18" viewBox="0 0 20 18" aria-hidden="true"><circle cx="10" cy="9" r="2" fill="none" stroke-width="1.4"/><path d="M10 7 V2 M12 9 H18 M10 11 V16 M8 9 H2" stroke-width="1.6" fill="none" stroke-linecap="round"/></svg>',
    empty     = '<svg width="54" height="54" viewBox="0 0 54 54" aria-hidden="true"><rect x="1" y="1" width="52" height="52" rx="8" fill="none" stroke-width="1.4"/><path d="M14 27 H21 M21 13 V41 M21 13 H33 M21 27 H33 M21 41 H33" stroke-width="1.6" fill="none" stroke-linecap="round"/><circle cx="37" cy="13" r="3"/><circle cx="37" cy="27" r="3"/><circle cx="37" cy="41" r="3"/></svg>',
    logo      = '<svg width="26" height="26" viewBox="0 0 26 26" aria-hidden="true"><rect x="0" y="0" width="26" height="26" rx="5" class="mm-logo-bg"/><path d="M6 13 H10 M10 6 V20 M10 6 H15 M10 13 H15 M10 20 H15" class="mm-logo-tree" stroke-width="1.6" fill="none" stroke-linecap="round"/><circle cx="17.5" cy="6" r="2.1" fill="#E69F00"/><circle cx="17.5" cy="13" r="2.1" fill="#56B4E9"/><circle cx="17.5" cy="20" r="2.1" fill="#009E73"/></svg>')
  shiny::HTML(svg)
}

#' Palette row: name + color strip
#' @noRd
.gui_palette_row <- function(tr, key) {
  tags <- shiny::tags
  cols <- if (identical(key, "auto")) PALETAS_ACCESIBLES$okabe else .gui_palettes()[[key]]
  cols <- .gui_hex(cols)
  tags$span(class = "mm-pal",
            .gui_t(.GUI_PALETTE_KEYS[[key]], tr, class = "mm-pal-name"),
            tags$span(class = "mm-strip", `aria-hidden` = "true",
                      lapply(cols, function(cl) tags$span(style = paste0("background:", cl)))))
}

#' Color input (HTML5 <input type='color'>, bound by multimapr.js)
#' @noRd
.gui_color_input <- function(id, value, label = NULL, hex_target = NULL, class = NULL) {
  shiny::tags$input(id = id, type = "color", class = paste("mm-color", class),
                    value = .gui_hex(value), `aria-label` = label, title = label,
                    `data-hex-target` = hex_target)
}

#' File input styled as in the design (label above, help below)
#' @noRd
.gui_file <- function(tr, id, label_key, help_key) {
  tags <- shiny::tags
  fi <- shiny::fileInput(id, label = NULL, buttonLabel = .gui_t("browse", tr),
                         placeholder = tr("no_file"))
  fi <- htmltools::tagQuery(fi)$find("input.form-control")$
    addAttrs(`data-i18n-ph` = "no_file", `aria-labelledby` = paste0(id, "-lbl"))$allTags()
  tags$div(class = "mm-file",
           .gui_t(label_key, tr, class = "mm-upper", id = paste0(id, "-lbl")),
           fi,
           .gui_t(help_key, tr, tag = tags$div, class = "mm-help"))
}

#' Choices of the select inputs that must be re-labelled on language change
#' @noRd
.gui_select_choices <- function(tr, which) {
  switch(which,
    tree_format   = stats::setNames(c("auto", "newick", "nexus", "tnt"),
                                    c(tr("auto"), "Newick", "NEXUS", "TNT")),
    matrix_format = stats::setNames(c("auto", "csv", "tnt", "nexus"),
                                    c(tr("auto"), "CSV", "TNT", "NEXUS")),
    csv_header    = stats::setNames(c("auto", "yes", "no"),
                                    c(tr("auto"), tr("yes"), tr("no"))),
    ladderize     = stats::setNames(c("TRUE", "right", "FALSE"),
                                    c(tr("ladder_bottom"), tr("ladder_top"), tr("ladder_none"))),
    stats_sort    = stats::setNames(c("none", "pct_missing", "pct_inapplicable", "character"),
                                    c(tr("sort_none"), tr("sort_missing"), tr("sort_inapp"),
                                      tr("sort_name"))),
    card_palette  = c(stats::setNames("global", tr("pal_use_global")),
                      stats::setNames(names(.GUI_PALETTE_KEYS),
                                      vapply(.GUI_PALETTE_KEYS, tr, character(1)))))
}

#' Card with a message log (see .gui_render_log in gui_server.R)
#' @noRd
.gui_log_card <- function(tr, id, title_key, kind = c("data", "map"), meta_id = NULL,
                          footer = NULL) {
  kind <- match.arg(kind)
  tags <- shiny::tags
  tags$div(class = "mm-card",
           tags$div(class = "mm-card-head",
                    .gui_t(title_key, tr, class = "mm-title"),
                    if (!is.null(meta_id))
                      tags$span(class = "mm-meta", shiny::textOutput(meta_id, inline = TRUE))),
           tags$div(class = "mm-card-body",
                    shiny::uiOutput(id, class = paste0("mm-log mm-log-", kind),
                                    role = "log", `aria-live` = "polite")),
           footer)
}

#' Empty state (icon + title + text + actions)
#' @noRd
.gui_empty_state <- function(tr, title_key, body, actions = NULL, class = NULL) {
  tags <- shiny::tags
  tags$div(class = paste("mm-empty", class),
           .gui_icon("empty"),
           .gui_t(title_key, tr, tag = tags$h3),
           tags$p(body),
           if (!is.null(actions)) tags$div(class = "mm-actions", actions))
}


# ==============================================================================
# SECTION 2 -- PAGE LAYOUT
# ==============================================================================

#' Builds the interface layout
#' @param tr Translator (initial language).
#' @noRd
.gui_ui <- function(tr) {
  tags <- shiny::tags
  lang <- attr(tr, "lang")

  brand <- tags$span(
    class = "mm-brand d-flex align-items-center gap-2",
    .gui_icon("logo"),
    tags$span(class = "mm-brand-text",
              tags$span(class = "mm-brand-name", "MultiMapR"),
              tags$span(class = "mm-brand-sub",
                        paste0("v", utils::packageVersion("MultiMapR"),
                               " \u00B7 R \u2265 4.1"))))

  tools <- tags$div(
    class = "mm-navtools d-flex align-items-center gap-2",
    tags$div(class = "mm-lang", role = "group", `aria-label` = "Language / Idioma",
             tags$button(type = "button", `data-lang` = "es",
                         class = if (lang == "es") "is-active", "ES"),
             tags$button(type = "button", `data-lang` = "en",
                         class = if (lang == "en") "is-active", "EN")),
    tags$button(type = "button", class = "mm-theme-btn", `aria-pressed` = "false",
                tags$span(class = "mm-theme-dot", `aria-hidden` = "true"),
                .gui_t("theme_light", tr, class = "mm-theme-label")))

  bslib::page_navbar(
    title        = brand,
    window_title = "MultiMapR",
    id           = "nav",
    lang         = lang,
    fillable     = FALSE,
    theme        = .gui_theme(),
    header       = tags$head(
      tags$link(rel = "icon", type = "image/svg+xml", href = "mmr-assets/logo.svg"),
      tags$link(rel = "stylesheet", href = "mmr-assets/multimapr.css"),
      tags$script(src = "mmr-assets/multimapr.js")),
    bslib::nav_panel(.gui_t("tab_data", tr), value = "data", .gui_page_data(tr)),
    bslib::nav_panel(.gui_t("tab_map", tr),  value = "map",  .gui_page_map(tr)),
    bslib::nav_panel(.gui_t("tab_help", tr), value = "help", .gui_page_help(tr)),
    bslib::nav_spacer(),
    bslib::nav_item(tools)
  )
}


# ==============================================================================
# SECTION 3 -- 1 - DATA
# ==============================================================================

#' @noRd
.gui_page_data <- function(tr) {
  tags <- shiny::tags
  t <- function(key, ...) .gui_t(key, tr, ...)

  input_card <- tags$div(
    class = "mm-card",
    tags$div(class = "mm-card-head",
             t("data_input", class = "mm-title"),
             shiny::uiOutput("src_badge", inline = TRUE)),
    tags$div(
      class = "mm-card-body",
      .gui_file(tr, "tree_file", "tree_file", "tree_file_help"),
      .gui_file(tr, "matrix_file", "matrix_file", "matrix_file_help"),
      bslib::accordion(
        id = "read_opts_acc", open = FALSE, class = "mm-adv",
        bslib::accordion_panel(
          title = t("read_opts"), value = "opts",
          tags$div(
            class = "mm-adv-grid",
            shiny::selectInput("tree_format", t("tree_format"),
                               .gui_select_choices(tr, "tree_format"), selectize = FALSE),
            shiny::selectInput("matrix_format", t("matrix_format"),
                               .gui_select_choices(tr, "matrix_format"), selectize = FALSE),
            tags$div(class = "mm-seg",
                     shiny::radioButtons("csv_sep", t("csv_sep"), inline = TRUE,
                                         choiceNames = list(t("auto_short"), ",", ";",
                                                            tags$span(title = "Tab", "\u21E5")),
                                         choiceValues = c("auto", ",", ";", "tab"))),
            shiny::textInput("species_col", t("species_col"), value = "1"),
            shiny::selectInput("csv_header", t("header"),
                               .gui_select_choices(tr, "csv_header"), selectize = FALSE),
            tags$div(class = "mm-full",
                     shiny::checkboxInput("normalize_spaces", t("normalize_spaces"), TRUE),
                     shiny::checkboxInput("prune_tips", t("prune"), FALSE))))),
      tags$div(class = "mm-actions",
               shiny::actionButton("load_btn", t("load"), class = "btn-mm-primary"),
               shiny::actionButton("clear_btn", t("clear"), class = "btn-mm-secondary"))))

  empty_card <- .gui_empty_state(
    tr, "empty_title",
    tags$span(t("empty_body_1", .noWS = "outside"), tags$b(t("load"), .noWS = c("outside", "inside")),
              t("empty_body_2", .noWS = "outside"),
              tags$b(t("tab_map"), .noWS = c("outside", "inside")), ".", .noWS = "inside"),
    actions = list(
      shiny::actionButton("example_btn", t("load_example"), class = "btn-mm-accent-outline"),
      shiny::actionButton("formats_btn", t("see_formats"), class = "btn-mm-secondary")))

  stats_card <- tags$div(
    class = "mm-card",
    tags$div(class = "mm-card-head",
             tags$div(class = "d-flex align-items-baseline gap-2",
                      t("stats_title", class = "mm-title"),
                      tags$span(class = "mm-meta", shiny::textOutput("stats_rows", inline = TRUE))),
             shiny::downloadButton("dl_stats_csv", t("download_csv"), icon = NULL,
                                   class = "btn-mm-secondary btn-mm-sm")),
    tags$div(class = "mm-table-wrap", shiny::uiOutput("stats_table")))

  completeness <- bslib::navset_card_tab(
    id = "cmp_tabs",
    bslib::nav_panel(t("per_char"), value = "per",
                     tags$div(class = "mm-plotbox",
                              shiny::plotOutput("stats_plot", height = "auto"))),
    bslib::nav_panel(t("heatmap"), value = "heat",
                     tags$div(class = "mm-plotbox",
                              shiny::plotOutput("heat_plot", height = "auto"))),
    bslib::nav_spacer(),
    bslib::nav_item(tags$div(
      class = "mm-tabtools d-flex align-items-center gap-2 py-1",
      shiny::conditionalPanel("input.cmp_tabs == 'per'",
                              shiny::selectInput("stats_sort", t("sort_by"),
                                                 .gui_select_choices(tr, "stats_sort"),
                                                 selectize = FALSE)),
      shiny::downloadButton("dl_cmp_png", "PNG", icon = NULL, class = "btn-mm-secondary btn-mm-sm"),
      shiny::downloadButton("dl_cmp_pdf", "PDF", icon = NULL, class = "btn-mm-secondary btn-mm-sm")))
  )
  completeness <- htmltools::tagAppendAttributes(completeness, class = "mm-tabcard")

  goto_map <- tags$div(class = "mm-card-foot",
                       shiny::actionButton("goto_map", t("goto_map"), class = "btn-mm-primary"))

  tags$div(
    class = "mm-page mm-page-data",
    tags$div(
      class = "mm-data-grid",
      input_card,
      tags$div(
        shiny::conditionalPanel("output.data_state == 'empty'", empty_card),
        shiny::conditionalPanel("output.data_state == 'error'",
                                shiny::uiOutput("data_error"),
                                .gui_log_card(tr, "data_log_err", "log_title", "data",
                                              meta_id = "data_log_err_meta")),
        shiny::conditionalPanel("output.data_state == 'loaded'",
                                shiny::uiOutput("kpis"),
                                stats_card))),
    shiny::conditionalPanel(
      "output.data_state == 'loaded'",
      tags$div(class = "mm-data-grid2",
               completeness,
               .gui_log_card(tr, "data_log_ok", "log_title", "data",
                             meta_id = "data_log_ok_meta", footer = goto_map))))
}


# ==============================================================================
# SECTION 4 -- 2 - MAP CHARACTERS
# ==============================================================================

#' @noRd
.gui_page_map <- function(tr) {
  tags <- shiny::tags
  t  <- function(key, ...) .gui_t(key, tr, ...)
  rc <- function(title, hint = NULL) .gui_rcard(tr, title, hint)

  step <- function(letter, key, ...) {
    bslib::accordion_panel(
      title = tags$span(class = "mm-step-title",
                        tags$span(class = "mm-step-badge", `aria-hidden` = "true", letter),
                        t(key, class = "mm-step-name"),
                        tags$span(class = "mm-step-sum",
                                  shiny::textOutput(paste0("sum_", letter), inline = TRUE))),
      value = letter, ...)
  }

  shape_choice <- function(glyph, key) {
    tags$span(title = tr(key), `data-i18n-title` = key,
              tags$span(`aria-hidden` = "true", glyph),
              t(key, class = "mm-sr"))
  }
  topo_choice <- function(key) tags$span(.gui_icon(key), t(key))
  corner_choice <- function(pos) {
    key <- paste0("corner_", pos)
    tags$span(title = tr(key), `data-i18n-title` = key,
              tags$span(class = paste0("mm-corner-box mm-c-", pos), `aria-hidden` = "true"),
              t(key, class = "mm-sr"))
  }

  # ---- Step A -- mapping type ------------------------------------------------
  step_a <- step(
    "A", "step_mapping",
    tags$div(class = "mm-rcards",
             shiny::radioButtons("mode", NULL,
                                 choiceNames  = list(rc("simple", "simple_hint"),
                                                     rc("ancestral", "ancestral_hint")),
                                 choiceValues = c("simple", "ancestral"))),
    tags$div(
      class = "mm-cond",
      shiny::conditionalPanel(
        "input.mode == 'simple'",
        .gui_label("presentation", tr),
        tags$div(class = "mm-rcards mm-compact",
                 shiny::radioButtons("presentation", NULL,
                                     choiceNames  = list(rc("pres_shapes", "pres_shapes_hint"),
                                                         rc("pres_labels", "pres_labels_hint")),
                                     choiceValues = c("shapes", "labels")))),
      shiny::conditionalPanel(
        "input.mode == 'ancestral'",
        .gui_label("viz", tr),
        tags$div(class = "mm-rcards mm-compact",
                 shiny::radioButtons("viz", NULL,
                                     choiceNames  = list(rc("viz_branches", "viz_branches_hint"),
                                                         rc("viz_tree_tips", "viz_tree_tips_hint")),
                                     choiceValues = c("branches", "tree_tips")))),
      shiny::conditionalPanel(
        "(input.mode == 'simple' && input.presentation == 'shapes') || (input.mode == 'ancestral' && input.viz == 'tree_tips')",
        tags$div(class = "mm-subblock",
                 .gui_label("figure", tr),
                 tags$div(class = "mm-btns mm-shapes",
                          shiny::radioButtons("shape", NULL, inline = TRUE,
                                              choiceNames = list(shape_choice("\u25CF", "circle"),
                                                                 shape_choice("\u25A0", "square"),
                                                                 shape_choice("\u25B2", "triangle"),
                                                                 shape_choice("\u25C6", "diamond")),
                                              choiceValues = c("21", "22", "24", "23"))),
                 shiny::sliderInput("shape_size", t("size"), min = 0.6, max = 3,
                                    value = 1.4, step = 0.1, ticks = FALSE))),
      shiny::conditionalPanel(
        "input.mode == 'ancestral'",
        tags$div(
          class = "mm-subblock",
          .gui_label("algorithm", tr),
          tags$div(class = "mm-rcards mm-compact",
                   shiny::radioButtons("algo", NULL,
                                       choiceNames  = list(rc("algo_majority", "algo_majority_hint"),
                                                           rc("algo_fitch", "algo_fitch_hint")),
                                       choiceValues = c("majority", "fitch"))),
          shiny::conditionalPanel(
            "input.algo == 'fitch'",
            tags$div(
              class = "mm-inner",
              .gui_label("resolution", tr),
              tags$div(class = "mm-rcards mm-compact",
                       shiny::radioButtons("optim", NULL,
                                           choiceNames = list(rc("acctran", "acctran_hint"),
                                                              rc("deltran", "deltran_hint"),
                                                              rc("unambiguous", "unambiguous_hint")),
                                           choiceValues = c("acctran", "deltran", "unambiguous"))),
              tags$div(class = "mm-ambig-row",
                       t("ambig_color", class = "mm-grow", id = "ambig_lbl"),
                       .gui_color_input("ambiguity_color", "#FF00FF", label = tr("ambig_color"),
                                        hex_target = "ambig_hex", class = "mm-color-lg"),
                       tags$span(id = "ambig_hex", class = "mm-hex", "#FF00FF")))))))
  )

  # ---- Step B -- characters -------------------------------------------------
  step_b <- step(
    "B", "step_chars",
    shiny::conditionalPanel("input.mode == 'simple'",
                            t("chars_hint_simple", tag = tags$div, class = "mm-note")),
    shiny::conditionalPanel("input.mode == 'ancestral'",
                            t("chars_hint_anc", tag = tags$div, class = "mm-note")),
    tags$input(id = "char_search", type = "search", class = "form-control mm-search",
               placeholder = tr("search_ph"), `data-i18n-ph` = "search_ph",
               `aria-label` = tr("search_ph"), autocomplete = "off"),
    shiny::uiOutput("chips"),
    shiny::uiOutput("trim_warn"),
    shiny::uiOutput("char_list")
  )

  # ---- Step C -- states and colors -----------------------------------------
  step_c <- step(
    "C", "step_colors",
    .gui_label("pal_global", tr),
    tags$div(class = "mm-palettes",
             shiny::radioButtons("pal_global", NULL,
                                 choiceNames  = lapply(names(.GUI_PALETTE_KEYS),
                                                       function(k) .gui_palette_row(tr, k)),
                                 choiceValues = names(.GUI_PALETTE_KEYS))),
    t("colors_hint", tag = tags$div, class = "mm-note mb-1"),
    shiny::uiOutput("char_cards")
  )

  # ---- Step D -- tree --------------------------------------------------------
  step_d <- step(
    "D", "step_tree",
    .gui_label("topology", tr),
    tags$div(class = "mm-btns mm-topo mb-3",
             shiny::radioButtons("topology", NULL, inline = TRUE,
                                 choiceNames  = list(topo_choice("phylogram"),
                                                     topo_choice("cladogram"),
                                                     topo_choice("fan")),
                                 choiceValues = c("phylogram", "cladogram", "fan"))),
    shiny::uiOutput("use_bl_ui"),
    shiny::selectInput("ladderize", t("ladderize"), .gui_select_choices(tr, "ladderize"),
                       selectize = FALSE),
    shiny::sliderInput("branch_width", t("branch_width"), min = 0.5, max = 8,
                       value = 2, step = 0.5, ticks = FALSE),
    tags$div(class = "mm-scale mb-3", t("thin_1"), t("normal_2"), t("thick_4"), tags$span("8")),
    shiny::sliderInput("tip_mult", t("terminal_stretch"), min = 1, max = 4,
                       value = 1, step = 0.25, pre = "\u00D7", ticks = FALSE),
    shiny::uiOutput("tip_note")
  )

  # ---- Step E -- legend and view ---------------------------------------------
  step_e <- step(
    "E", "step_view",
    .gui_label("legend_corner", tr),
    tags$div(class = "mm-btns mm-corner",
             shiny::radioButtons("legend_corner", NULL, inline = TRUE,
                                 choiceNames  = lapply(c("topleft", "topright", "bottomleft",
                                                         "bottomright"), corner_choice),
                                 choiceValues = c("topleft", "topright", "bottomleft", "bottomright"),
                                 selected = "bottomleft")),
    tags$div(class = "mm-corner-label mb-3", shiny::textOutput("corner_label", inline = TRUE)),
    shiny::sliderInput("preview_height", t("preview_height"), min = 360, max = 1000,
                       value = 640, step = 20, post = " px", ticks = FALSE)
  )

  # ---- Step F -- export ------------------------------------------------------
  step_f <- step(
    "F", "step_export",
    .gui_label("format", tr),
    tags$div(class = "mm-rcards mm-compact mb-3",
             shiny::radioButtons("fmt", NULL,
                                 choiceNames  = list(rc("fmt_png", "fmt_png_hint"),
                                                     rc("fmt_pdf", "fmt_pdf_hint")),
                                 choiceValues = c("png", "pdf"))),
    shiny::textInput("export_name", t("filename"), value = ""),
    .gui_label("dims", tr),
    tags$div(class = "mm-rcards mm-compact mm-row",
             shiny::radioButtons("dims", NULL,
                                 choiceNames  = list(rc("dims_auto"), rc("dims_custom")),
                                 choiceValues = c("auto", "custom"))),
    shiny::conditionalPanel("input.dims == 'auto'",
                            tags$div(class = "mm-dims-note",
                                     shiny::textOutput("dims_note", inline = TRUE))),
    shiny::conditionalPanel("input.dims == 'custom'",
                            tags$div(class = "mm-custom-dims",
                                     shiny::numericInput("exp_w", t("width_in"), value = 12,
                                                         min = 1, step = 0.5),
                                     shiny::numericInput("exp_h", t("height_in"), value = 10,
                                                         min = 1, step = 0.5))),
    tags$div(class = "mm-exp-actions",
             shiny::downloadButton("dl_map", t("download"), icon = NULL, class = "btn-mm-primary"),
             shiny::actionButton("save_map", t("save_exports"), class = "btn-mm-secondary")),
    tags$div(class = "mm-path", shiny::textOutput("exports_path", inline = TRUE)),
    shiny::uiOutput("save_ok")
  )

  side <- tags$aside(
    class = "mm-side", `aria-label` = tr("config"),
    tags$div(class = "mm-side-head",
             t("config", tag = tags$h2),
             tags$span(class = "mm-meta", shiny::textOutput("side_counter", inline = TRUE))),
    bslib::accordion(id = "steps", multiple = FALSE, open = "A", class = "mm-steps",
                     step_a, step_b, step_c, step_d, step_e, step_f))

  # ---- Preview column ---------------------------------------------------------
  topbar <- tags$div(
    class = "mm-card mm-topbar",
    t("preview", class = "mm-upper"),
    tags$span(class = "mm-pill", shiny::textOutput("cfg_summary", inline = TRUE)),
    tags$span(class = "mm-spacer"),
    shiny::conditionalPanel("input.mode == 'ancestral' && input.algo == 'fitch'",
                            shiny::actionButton("cmp_toggle", t("compare"),
                                                class = "btn-mm-secondary btn-mm-toggle")),
    shiny::actionButton("recalc", t("recalc"), class = "btn-mm-secondary"),
    shiny::downloadButton("dl_map_top", shiny::textOutput("export_btn_label", inline = TRUE),
                          icon = NULL, class = "btn-mm-primary"))

  cmp_col <- function(m) {
    tags$div(class = "mm-cmp-col",
             tags$div(class = "mm-cmp-colhead",
                      tags$b(paste0("Fitch \u00B7 ", toupper(m))),
                      shiny::uiOutput(paste0("cmp_use_", m), inline = TRUE)),
             tags$div(class = "mm-cmp-thumb",
                      shiny::imageOutput(paste0("cmp_img_", m), height = "230px", width = "100%")),
             tags$div(class = "mm-cmp-note", shiny::textOutput(paste0("cmp_note_", m), inline = TRUE)))
  }
  compare <- shiny::conditionalPanel(
    "output.cmp_visible == 'yes'",
    tags$div(class = "mm-card mm-cmp",
             tags$div(class = "mm-card-head",
                      t("cmp_title", class = "mm-title"),
                      tags$button(id = "cmp_close", type = "button", class = "mm-x action-button",
                                  title = tr("close"), `data-i18n-title` = "close",
                                  `aria-label` = tr("close"), "\u00D7")),
             tags$div(class = "mm-cmp-grid", cmp_col("acctran"), cmp_col("deltran"))))

  canvas <- tags$div(
    class = "mm-canvas-card",
    tags$div(id = "mm-canvas", class = "mm-canvas",
             shiny::uiOutput("preview_state"),
             shiny::imageOutput("preview_img", width = "100%", height = "auto")),
    tags$div(class = "mm-computing", role = "status",
             tags$div(class = "mm-spinner", `aria-hidden` = "true"),
             t("computing", tag = tags$div, class = "mm-computing-title"),
             tags$div(class = "mm-computing-detail",
                      shiny::textOutput("computing_detail", inline = TRUE))),
    t("canvas_foot", tag = tags$div, class = "mm-canvas-foot"))

  main <- tags$div(
    class = "mm-main", style = "min-width:0",
    topbar, compare, canvas,
    .gui_log_card(tr, "map_log", "map_log_title", "map", meta_id = "session_stamp"))

  tags$div(class = "mm-page mm-page-map",
           tags$div(class = "mm-map-grid", side, main))
}


# ==============================================================================
# SECTION 5 -- HELP
# ==============================================================================

#' @noRd
.gui_page_help <- function(tr) {
  tags <- shiny::tags
  t <- function(key, ...) .gui_t(key, tr, ...)
  card <- function(title_key, ...) {
    tags$div(class = "mm-card",
             tags$div(class = "mm-card-head", t(title_key, class = "mm-title")),
             tags$div(class = "mm-card-body", ...))
  }
  code <- function(x) tags$code(x)

  csv_rows <- paste("Saccopteryx_bilineata,0,1,?,0",
                    "Noctilio_albiventris,1,1,0,1",
                    "Noctilio_leporinus,1,1,0,1",
                    "Mystacina_tuberculata,0,-,2,0", sep = "\n")

  modes <- c("simple", "ancestral", "acctran", "deltran", "unambiguous", "algo_majority")

  tags$div(
    class = "mm-page mm-page-help",
    t("help_title", tag = tags$h1, class = "mm-help-title"),
    t("help_intro", tag = tags$p, class = "mm-help-intro"),

    card("help_formats",
         tags$table(class = "mm-help-table",
                    tags$tr(t("help_tree", tag = tags$th, scope = "row"),
                            tags$td("Newick ", code(".tre .tree .nwk"), " \u00B7 NEXUS ",
                                    code(".nex .nexus"), " \u00B7 TNT. ", t("help_tree_note"))),
                    tags$tr(t("help_matrix", tag = tags$th, scope = "row"),
                            tags$td(t("help_matrix_csv"), " \u00B7 TNT ", code("xread"),
                                    " \u00B7 NEXUS ", code("MATRIX"))),
                    tags$tr(t("help_names", tag = tags$th, scope = "row"),
                            tags$td(t("help_names_note"))))),

    tags$div(class = "mm-help-grid2",
             card("help_csv_header",
                  tags$pre(class = "mm-pre", .noWS = "inside", paste0("Species,C1,C2,C3,diet\n", csv_rows))),
             card("help_csv_noheader",
                  tags$pre(class = "mm-pre", .noWS = "inside", paste0(csv_rows, "\n"),
                           t("help_csv_noheader_note", class = "mm-cmt", .noWS = "outside")))),

    card("help_poly",
         tags$div(class = "mm-codes",
                  code("0/1"),          t("help_poly_slash"),
                  code("urban,forest"), t("help_poly_comma"),
                  code("01"),           t("help_poly_concat"),
                  code("?"),            t("help_missing"),
                  code("-"),            t("help_inapp"))),

    card("help_modes",
         tags$div(class = "mm-modes",
                  lapply(modes, function(m) {
                    tags$div(class = "mm-mode",
                             t(m, tag = tags$h4),
                             t(paste0("help_", m), tag = tags$p))
                  }))),

    card("help_console",
         t("help_console_text", tag = tags$p, class = "mb-0")),

    card("help_cite",
         tags$div(class = "mm-cite",
                  tags$p(paste0("Alc\u00E1ntara Estrada, K. I., Linares R\u00EDos, J. & D\u00EDaz Cruz, J. A. (",
                                format(Sys.Date(), "%Y"), "). "),
                         tags$i("MultiMapR: Mapping and Visualization of Discrete Characters on Phylogenies"),
                         paste0(". R package version ", utils::packageVersion("MultiMapR"), ".")),
                  tags$p("Paradis, E. & Schliep, K. (2019). ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. ",
                         tags$i("Bioinformatics"), " 35, 526\u2013528."),
                  tags$p("Fitch, W. M. (1971). Toward defining the course of evolution: minimum change for a specific tree topology. ",
                         tags$i("Systematic Zoology"), " 20, 406\u2013416."),
                  tags$p("Swofford, D. L. & Maddison, W. P. (1987). Reconstructing ancestral character states under Wagner parsimony. ",
                         tags$i("Mathematical Biosciences"), " 87, 199\u2013229."),
                  tags$p("Okabe, M. & Ito, K. (2008). Color Universal Design (CUD): how to make figures and presentations that are friendly to colorblind people.")))
  )
}
