# Shiny + bslib interface (run_multimapr_app()). It builds the same config list as
# setup_mapping_config() and drives .render_configured_mapping(), so its figures
# match those of execute_phylogeny(). Layout lives in gui_ui.R, server logic in
# gui_server.R, assets (css, js, i18n.txt) in inst/app/. shiny and bslib are only
# suggested dependencies, checked when the app is launched.


# Public entry point

#' Launch the MultiMapR graphical interface
#'
#' Opens an interactive application (Shiny) that gives access to every
#' MultiMapR feature without typing answers in the console: loading a tree
#' and a character matrix, reviewing missing / inapplicable data, simple
#' mapping on terminals, ancestral reconstruction (depth-weighted majority or
#' Fitch ACCTRAN / DELTRAN / Unambiguous, with a side-by-side ACCTRAN / DELTRAN
#' comparison), one to three superimposed characters, the three tree
#' topologies, and PNG / PDF export.
#'
#' Plots are drawn with the same rendering engine used by
#' \code{\link{execute_phylogeny}}, so a figure configured in the interface
#' and one configured in the console menus are identical. The interface is
#' bilingual (English / Spanish, switchable at any time from the navigation
#' bar) and has a light and a dark theme; the figure preview is always drawn
#' on white, as in the exported file.
#'
#' @param tree        Optional tree to preload: a \code{phylo} object or a
#'                    path to a tree file. Must be given together with
#'                    \code{characters}.
#' @param characters  Optional character matrix to preload: a data.frame with
#'                    a \code{"Species"} column or a path to a matrix file
#'                    (CSV, TNT or NEXUS).
#' @param lang        Initial interface language: \code{"en"} (English,
#'                    default) or \code{"es"} (Spanish).
#' @param launch.browser Open the interface in the web browser (default
#'                    \code{TRUE}). Passed to \code{shiny::runApp()}.
#' @param port        Port for the local server. \code{NULL} (default) picks a
#'                    free one.
#' @param ...         Further arguments for \code{load_data()} when
#'                    \code{tree} and \code{characters} are file paths (e.g.
#'                    \code{normalize_spaces = TRUE}).
#' @return Called for its side effect (runs the app until it is closed).
#'   Files saved with the "Save to Exports/" button are written to the
#'   \code{Exports} folder of the current working directory, as in the
#'   console workflow.
#' @examples
#' if (interactive()) {
#'   tree <- system.file("extdata", "bats_tre.tre", package = "MultiMapR")
#'   chars <- system.file("extdata", "bats_matrix.csv", package = "MultiMapR")
#'
#'   # Empty interface: load the files from the "Data" tab
#'   run_multimapr_app()
#'
#'   # Spanish interface with data already loaded
#'   run_multimapr_app(tree, chars, lang = "es")
#' }
#' @export
run_multimapr_app <- function(tree = NULL, characters = NULL,
                              lang = c("en", "es"),
                              launch.browser = TRUE,
                              port = NULL, ...) {
  missing_pkgs <- c("shiny", "bslib")[!vapply(c("shiny", "bslib"), requireNamespace,
                                              logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0)
    stop("The MultiMapR interface needs the package(s): ",
         paste(missing_pkgs, collapse = ", "), ".\n  Install with: install.packages(c(",
         paste0('"', missing_pkgs, '"', collapse = ", "), "))", call. = FALSE)
  lang <- match.arg(lang)

  initial <- NULL
  if (!is.null(tree) || !is.null(characters)) {
    if (is.null(tree) || is.null(characters))
      stop("`tree` and `characters` must be supplied together.", call. = FALSE)
    if (is.character(tree) && is.character(characters)) {
      res <- .gui_capture(load_data(tree, characters, ...))
      initial <- list(tree = res$value$tree, characters = res$value$characters,
                      kind = "files", names = basename(c(tree, characters)),
                      log = res$log)
    } else if (inherits(tree, "phylo") && is.data.frame(characters)) {
      initial <- list(tree = tree, characters = characters,
                      kind = "objects", names = c("phylo", "data.frame"),
                      log = character(0))
    } else {
      stop("`tree` and `characters` must be two file paths, or a 'phylo' object ",
           "and a data.frame.", call. = FALSE)
    }
  }

  app <- .mmr_app(initial = initial, lang = lang)
  shiny::runApp(app, launch.browser = launch.browser, port = port)
}


#' Builds the MultiMapR Shiny application object
#' @noRd
.mmr_app <- function(initial = NULL, lang = "en") {
  www <- .gui_app_file("www")
  if (!nzchar(www))
    stop("MultiMapR interface assets not found (inst/app/www). Reinstall the package.",
         call. = FALSE)
  shiny::addResourcePath("mmr-assets", www)
  tr <- .gui_translator(lang)
  shiny::shinyApp(ui     = .gui_ui(tr),
                  server = .gui_server(lang, initial))
}


# Assets, theme and translations

#' Path to a file shipped in inst/app
#' @noRd
.gui_app_file <- function(...) system.file("app", ..., package = "MultiMapR")

#' bslib theme with the design tokens (dark mode comes from multimapr.css)
#' @noRd
.gui_theme <- function() {
  sans <- bslib::font_collection(
    bslib::font_google("Source Sans 3", wght = c(400, 600, 700), ital = c(0, 1),
                       local = FALSE),
    "system-ui", "-apple-system", "Segoe UI", "sans-serif")
  mono <- bslib::font_collection(
    bslib::font_google("IBM Plex Mono", wght = c(400, 500), local = FALSE),
    "ui-monospace", "Consolas", "monospace")
  bslib::bs_theme(
    version   = 5,
    bg        = "#f6f6f4", fg = "#1b1c1a",
    primary   = "#2a6f97", secondary = "#5d5f59", info = "#2a6f97",
    success   = "#1e6b4a", warning = "#8a5d00", danger = "#9d2b2b",
    base_font = sans, code_font = mono, heading_font = sans,
    "border-color"  = "#dcdcd6",
    "border-radius" = "6px",
    "font-size-base" = "0.84375rem"
  )
}

#' Interface dictionary (inst/app/i18n.txt: key | es | en)
#' @noRd
.gui_dict <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    path <- .gui_app_file("i18n.txt")
    if (!nzchar(path)) stop("MultiMapR interface dictionary not found (inst/app/i18n.txt).")
    lines <- readLines(path, encoding = "UTF-8", warn = FALSE)
    lines <- lines[nzchar(trimws(lines)) & !startsWith(lines, "#")]
    parts <- strsplit(lines, "|", fixed = TRUE)
    parts <- parts[vapply(parts, length, integer(1)) >= 3L]
    dict  <- lapply(parts, function(p) c(es = enc2utf8(p[2]), en = enc2utf8(p[3])))
    names(dict) <- vapply(parts, `[`, character(1), 1L)
    cache <<- dict
    dict
  }
})

#' Returns a translation function for the interface strings
#' @noRd
.gui_translator <- function(lang) {
  dict <- .gui_dict()
  tr <- function(key, ...) {
    txt <- dict[[key]]
    if (is.null(txt)) return(key)
    txt <- unname(txt[[lang]] %||% txt[["en"]])
    if (length(list(...)) > 0) sprintf(txt, ...) else txt
  }
  attr(tr, "lang") <- lang
  tr
}


# Non-reactive helpers

#' Evaluates an expression capturing console output, messages and warnings
#' @noRd
.gui_capture <- function(expr) {
  msgs  <- character(0)
  value <- NULL
  out <- utils::capture.output(
    value <- withCallingHandlers(
      expr,
      message = function(m) {
        msgs <<- c(msgs, sub("\n$", "", conditionMessage(m)))
        invokeRestart("muffleMessage")
      },
      warning = function(w) {
        msgs <<- c(msgs, paste("Warning:", conditionMessage(w)))
        invokeRestart("muffleWarning")
      }),
    type = "output")
  out <- out[nzchar(trimws(out))]
  list(value = value, log = c(out, msgs))
}

#' Converts any R color (name or hex) to a "#RRGGBB" string
#' @noRd
.gui_hex <- function(col) {
  m <- grDevices::col2rgb(col)
  grDevices::rgb(m[1, ], m[2, ], m[3, ], maxColorValue = 255)
}

#' Palettes offered by the interface (all defined in cli_menu.R)
#' @noRd
.gui_palettes <- function() {
  list(okabe  = ACCESSIBLE_PALETTES$okabe,
       tol     = ACCESSIBLE_PALETTES$tol,
       plasma  = ACCESSIBLE_PALETTES$plasma,
       reds    = MULTIMAP_RAMPS$reds,
       greens  = MULTIMAP_RAMPS$greens,
       blues   = MULTIMAP_RAMPS$blues,
       golds   = MULTIMAP_RAMPS$golds,
       purples = MULTIMAP_RAMPS$purples,
       wc      = PRESET_PALETTES$wc,
       MM1     = PRESET_PALETTES$MM1)
}

#' Dictionary key (or literal name) for each palette
#' @noRd
.GUI_PALETTE_KEYS <- c(auto = "pal_auto", okabe = "pal_okabe", tol = "pal_tol",
                       plasma = "pal_plasma", reds = "pal_reds", greens = "pal_greens",
                       blues = "pal_blues", golds = "pal_golds", purples = "pal_purples",
                       wc = "pal_wc", MM1 = "pal_mm1")

#' Resolves a palette choice to n hex colors
#' @noRd
.gui_palette_colors <- function(choice, n, slot = 1L, mapping_type = 1L, n_selected = 1L) {
  pals <- .gui_palettes()
  pal <- if (!identical(choice, "auto") && choice %in% names(pals)) {
    pals[[choice]]
  } else if (mapping_type == 2L && n_selected > 1L) {
    MULTIMAP_RAMPS[[min(slot, length(MULTIMAP_RAMPS))]]
  } else if (mapping_type == 2L) {
    ACCESSIBLE_PALETTES$okabe
  } else {
    ACCESSIBLE_PALETTES[[((slot - 1L) %% length(ACCESSIBLE_PALETTES)) + 1L]]
  }
  if (n <= 0) return(character(0))
  .gui_hex(rep(pal, length.out = n))
}

#' Cleans a user-typed file name (no extension, no path, no forbidden chars)
#' @noRd
.gui_clean_filename <- function(x, fallback) {
  x <- trimws(x %||% "")
  x <- sub("\\.(png|pdf)$", "", x, ignore.case = TRUE)
  x <- gsub("[\\\\/:*?\"<>|]", "_", x)
  if (nzchar(x)) x else gsub("[\\\\/:*?\"<>|]", "_", fallback)
}

#' Default color for missing / inapplicable tips in the interface (gray70)
#' @noRd
.GUI_MISSING_DEFAULT <- "#B3B3B3"

#' Figure dimensions in inches, as computed by .export_device()
#' @noRd
.gui_figure_dims <- function(n_tips, tree_type, width = NULL, height = NULL,
                             branch_angle = NULL) {
  h <- n_tips * 0.25 + 2
  w <- .emtree_auto_width(tree_type, h, branch_angle)
  if (identical(tree_type, "fan")) {
    h <- max(h, w)
    w <- h
  }
  list(width  = if (!is.null(width))  width  else w,
       height = if (!is.null(height)) height else h,
       auto_height = h)
}

#' Runs a drawing expression with a throw-away device as the active device
#' @noRd
.gui_with_sink_device <- function(expr) {
  prev <- dev.cur()
  pdf(NULL)
  sink_dev <- dev.cur()
  on.exit({
    if (sink_dev %in% dev.list()) dev.off(sink_dev)
    if (prev > 1L && prev %in% dev.list()) dev.set(prev)
  }, add = TRUE)
  force(expr)
}

#' Draws a configured mapping into a PNG file
#' @noRd
.gui_render_png <- function(spec, file, width_px, dims) {
  res      <- width_px / dims$width
  height_px <- max(50L, round(dims$height * res))
  prev <- dev.cur()
  grDevices::png(file, width = round(width_px), height = height_px, res = res, bg = "white")
  t0 <- Sys.time()
  out <- tryCatch(
    .gui_capture(do.call(.render_configured_mapping,
                         c(list(phylogeny = spec$tree, config = spec$config), spec$render))),
    finally = .close_device_restore(prev))
  list(log = out$log, seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")),
       height_px = height_px)
}

#' Exports a configured mapping through the package export engine
#' @noRd
.gui_export_mapping <- function(spec, format, filename, width, height, out_dir) {
  old_wd <- setwd(out_dir)
  on.exit(setwd(old_wd), add = TRUE)

  config <- spec$config
  config$export_filename <- filename
  config$export_format   <- format
  config$width           <- width
  config$height          <- height

  res <- .gui_with_sink_device(.gui_capture(
    do.call(.render_configured_mapping,
            c(list(phylogeny = spec$tree, config = config), spec$render))))
  list(path = normalizePath(file.path(out_dir, "Exports", paste0(basename(filename), ".", format)),
                            mustWork = FALSE),
       log  = res$log)
}

#' Counts the branches left ambiguous by a Fitch optimization
#' @noRd
.gui_count_ambiguous <- function(spec, mode) {
  cfg <- spec$config
  cfg$fitch_mode <- mode
  ambig <- toupper(.gui_hex(.get_fitch_ambig_color(cfg)))
  tree  <- spec$tree
  resolve_color <- function(value, state_colors) {
    value <- as.character(value)
    if (is.na(value) || value == "") return("gray70")
    if (value %in% names(state_colors)) return(state_colors[[value]])
    "gray70"
  }
  tryCatch({
    total <- 0L
    for (ch in cfg$characters) {
      tip_colors <- sapply(cfg$aligned_data[[ch]], resolve_color,
                           state_colors = cfg$colors_by_character[[ch]])
      ec <- external_algorithm(tree, tip_colors, init_edge_colors(tree), cfg)
      total <- total + sum(toupper(.gui_hex(unname(ec))) == ambig)
    }
    as.integer(total)
  }, error = function(e) NA_integer_)
}
