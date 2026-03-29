library(shiny)
library(bslib)
library(DT)
library(AMRgen)
library(readxl)

addResourcePath("assets", normalizePath("."))

capture_conditions <- function(expr) {
  warnings <- character()
  messages <- character()

  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )

  list(
    value = value,
    warnings = unique(warnings),
    messages = unique(messages)
  )
}

phenotype_columns <- function(data) {
  cols <- grep("^(pheno|ecoff)", names(data), value = TRUE)
  cols[!vapply(data[cols], is.list, logical(1))]
}

default_export_pheno_col <- function(choices, interpret_eucast = FALSE, interpret_clsi = FALSE) {
  if (!length(choices)) {
    return(character())
  }

  if (isTRUE(interpret_eucast)) {
    eucast_choice <- choices[grepl("eucast", choices, ignore.case = TRUE)][1]
    if (!is.na(eucast_choice)) {
      return(eucast_choice)
    }
  }

  if (isTRUE(interpret_clsi)) {
    clsi_choice <- choices[grepl("clsi", choices, ignore.case = TRUE)][1]
    if (!is.na(clsi_choice)) {
      return(clsi_choice)
    }
  }

  if ("pheno_provided" %in% choices) {
    return("pheno_provided")
  }

  choices[[1]]
}

format_choice_label <- function(x) {
  switch(
    x,
    sensititre = "Sensititre",
    vitek = "Vitek",
    phoenix = "Phoenix",
    microscan = "MicroScan",
    whonet = "WHONET",
    ncbi = "NCBI browser download",
    ncbi_biosample = "NCBI cloud download",
    ebi_web = "EBI browser download",
    ebi_ftp = "EBI ftp download",
    x
  )
}

guess_format_from_filename <- function(filename) {
  filename <- tolower(filename %||% "")

  format_patterns <- c(
    ebi_ftp = "(^|[^a-z])ebi[_. -]?ftp([^a-z]|$)",
    ebi_web = "(^|[^a-z])ebi[_. -]?web([^a-z]|$)",
    ebi = "(^|[^a-z])ebi([^a-z]|$)",
    ncbi = "(^|[^a-z])ncbi([^a-z]|$)",
    whonet = "(^|[^a-z])whonet([^a-z]|$)",
    sensititre = "(^|[^a-z])sensititre([^a-z]|$)",
    microscan = "(^|[^a-z])microscan([^a-z]|$)",
    phoenix = "(^|[^a-z])phoenix([^a-z]|$)",
    vitek = "(^|[^a-z])vitek([^a-z]|$)"
  )

  for (format_name in names(format_patterns)) {
    if (grepl(format_patterns[[format_name]], filename, perl = TRUE)) {
      return(format_name)
    }
  }

  NULL
}

prepare_import_input <- function(file_info, format = NULL) {
  ext <- tolower(tools::file_ext(file_info$name))

  if (identical(format, "phoenix") && ext %in% c("xlsx", "xls")) {
    return(list(
      input = file_info$datapath,
      note = "Passed the Phoenix Excel file directly to AMRgen so the package can apply its own header detection."
    ))
  }

  if (ext %in% c("xlsx", "xls")) {
    sheet_names <- readxl::excel_sheets(file_info$datapath)
    data <- readxl::read_excel(file_info$datapath, sheet = sheet_names[[1]])

    return(list(
      input = as.data.frame(data),
      note = paste("Read Excel workbook sheet", shQuote(sheet_names[[1]]), "before dispatching to import_pheno().")
    ))
  }

  list(
    input = file_info$datapath,
    note = NULL
  )
}

read_tabular_input <- function(file_info) {
  ext <- tolower(tools::file_ext(file_info$name))

  if (ext %in% c("xlsx", "xls")) {
    sheet_names <- readxl::excel_sheets(file_info$datapath)
    return(as.data.frame(readxl::read_excel(file_info$datapath, sheet = sheet_names[[1]])))
  }

  if (ext == "csv") {
    return(utils::read.csv(file_info$datapath, stringsAsFactors = FALSE, check.names = FALSE))
  }

  if (ext %in% c("tsv", "txt")) {
    first_line <- readLines(file_info$datapath, n = 1, warn = FALSE)
    sep <- if (grepl("\t", first_line, fixed = TRUE)) "\t" else ","
    return(utils::read.table(
      file_info$datapath,
      header = TRUE,
      sep = sep,
      quote = "\"",
      stringsAsFactors = FALSE,
      check.names = FALSE,
      comment.char = ""
    ))
  }

  stop("Unsupported BioSample mapping file type. Please use CSV, TSV, TXT, XLS, or XLSX.")
}

normalise_name <- function(x) {
  tolower(gsub("[^a-z0-9]", "", x %||% ""))
}

find_biosample_column <- function(data) {
  cols <- names(data)
  match <- cols[vapply(cols, function(x) normalise_name(x) == "biosample", logical(1))]

  if (!length(match)) {
    return(NULL)
  }

  match[[1]]
}

default_mapping_id_column <- function(data, biosample_col) {
  id_candidates <- setdiff(names(data), biosample_col)

  if (!length(id_candidates)) {
    return(NULL)
  }

  normalised <- stats::setNames(vapply(id_candidates, normalise_name, character(1)), id_candidates)

  sample_match <- names(normalised)[grepl("sample", normalised, fixed = TRUE)]
  if (length(sample_match)) {
    return(sample_match[[1]])
  }

  identifier_match <- names(normalised)[grepl("identifier", normalised, fixed = TRUE)]
  if (length(identifier_match)) {
    return(identifier_match[[1]])
  }

  id_match <- names(normalised)[grepl("id", normalised, fixed = TRUE)]
  if (length(id_match)) {
    return(id_match[[1]])
  }

  id_candidates[[1]]
}

apply_biosample_mapping <- function(data, mapping_df, sample_id_col, biosample_col) {
  if (is.null(mapping_df) || !nrow(mapping_df)) {
    return(data)
  }

  sample_ids <- as.character(mapping_df[[sample_id_col]])
  biosamples <- as.character(mapping_df[[biosample_col]])
  valid <- nzchar(sample_ids) & nzchar(biosamples)

  if (!any(valid)) {
    stop("No usable BioSample mappings were found in the uploaded mapping file.")
  }

  lookup <- stats::setNames(biosamples[valid], sample_ids[valid])
  mapped_ids <- unname(lookup[as.character(data$id)])

  if (any(!is.na(mapped_ids))) {
    data$id <- ifelse(is.na(mapped_ids), as.character(data$id), mapped_ids)
  }

  data
}

has_invalid_biosample_ids <- function(ids, valid_prefixes) {
  ids <- unique(as.character(ids))
  ids <- ids[!is.na(ids) & nzchar(ids)]

  if (!length(ids)) {
    return(FALSE)
  }

  pattern <- paste0("^(", paste(valid_prefixes, collapse = "|"), ")")
  any(!grepl(pattern, ids, ignore.case = TRUE))
}

all_biosample_ids_match <- function(ids, valid_prefixes) {
  ids <- unique(as.character(ids))
  ids <- ids[!is.na(ids) & nzchar(ids)]

  if (!length(ids)) {
    return(FALSE)
  }

  pattern <- paste0("^(", paste(valid_prefixes, collapse = "|"), ")")
  all(grepl(pattern, ids, ignore.case = TRUE))
}

format_importer_name <- function(format) {
  switch(
    format,
    sensititre = "import_sensititre_ast",
    vitek = "import_vitek_ast",
    phoenix = "import_phoenix_ast",
    microscan = "import_microscan_ast",
    whonet = "import_whonet_ast",
    ncbi = "import_ncbi_ast",
    ncbi_biosample = "import_ncbi_biosample",
    ebi = "import_ebi_ast",
    ebi_web = "import_ebi_ast",
    ebi_ftp = "import_ebi_ast_ftp",
    NULL
  )
}

format_specific_args <- function(format) {
  importer_name <- format_importer_name(format)

  if (is.null(importer_name) || !exists(importer_name, where = asNamespace("AMRgen"), inherits = FALSE)) {
    return(list())
  }

  args <- formals(getFromNamespace(importer_name, "AMRgen"))
  common_args <- c(
    "input", "source", "species", "ab",
    "interpret_eucast", "interpret_clsi", "interpret_ecoff"
  )

  args[setdiff(names(args), common_args)]
}

extra_arg_input_id <- function(arg_name) {
  paste0("extra_arg__", arg_name)
}

render_extra_arg_input <- function(arg_name, default_value) {
  input_id <- extra_arg_input_id(arg_name)

  if (is.logical(default_value)) {
    return(
      checkboxInput(
        inputId = input_id,
        label = arg_name,
        value = isTRUE(default_value)
      )
    )
  }

  if (is.numeric(default_value) && length(default_value) == 1) {
    return(
      numericInput(
        inputId = input_id,
        label = arg_name,
        value = default_value
      )
    )
  }

  placeholder <- if (is.null(default_value)) "" else as.character(default_value)

  textInput(
    inputId = input_id,
    label = arg_name,
    value = placeholder
  )
}

collect_extra_import_args <- function(input, format) {
  args <- format_specific_args(format)

  if (!length(args)) {
    return(list())
  }

  values <- lapply(names(args), function(arg_name) {
    input[[extra_arg_input_id(arg_name)]]
  })
  names(values) <- names(args)

  for (arg_name in names(values)) {
    default_value <- args[[arg_name]]

    if (is.null(default_value) && is.character(values[[arg_name]]) && !nzchar(values[[arg_name]])) {
      values[[arg_name]] <- NULL
    }
  }

  values
}

supported_export_args <- function(fn_name) {
  if (!exists(fn_name, where = asNamespace("AMRgen"), inherits = FALSE)) {
    return(character())
  }

  names(formals(getFromNamespace(fn_name, "AMRgen")))
}

filter_supported_args <- function(arg_list, fn_name) {
  supported <- supported_export_args(fn_name)
  arg_list[names(arg_list) %in% supported]
}

null_if_blank <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }

  if (is.character(x) && length(x) == 1 && !nzchar(x)) {
    return(NULL)
  }

  x
}

format_amrgen_error <- function(message) {
  if (grepl("AMRgen[.]rdb.+corrupt", message, ignore.case = TRUE)) {
    return(paste(
      "AMRgen appears to be installed incorrectly or is corrupted.",
      "Please reinstall it with remotes::install_github('AMRverse/AMRgen')",
      "and fully restart the R session before relaunching the app."
    ))
  }

  message
}

export_standard_for_pheno <- function(pheno_col) {
  pheno_col <- tolower(pheno_col %||% "")

  if (grepl("eucast", pheno_col, fixed = TRUE)) {
    return("EUCAST 2025")
  }

  if (grepl("clsi", pheno_col, fixed = TRUE)) {
    return("CLSI 2025")
  }

  NULL
}

guideline_for_pheno <- function(pheno_col) {
  pheno_col <- tolower(pheno_col %||% "")

  if (grepl("eucast", pheno_col, fixed = TRUE)) {
    return("EUCAST")
  }

  if (grepl("clsi", pheno_col, fixed = TRUE)) {
    return("CLSI")
  }

  NULL
}

version_for_pheno <- function(pheno_col, current_value = NULL) {
  guideline <- guideline_for_pheno(pheno_col)

  if (!is.null(guideline)) {
    return("2025")
  }

  if (!is.null(current_value) && nzchar(current_value)) {
    return(current_value)
  }

  ""
}

ncbi_guideline_choices <- function() {
  c("BSAC", "CLSI", "DIN", "EUCAST", "NARMS", "NCCLS", "SFM", "SIR", "WRG", "missing")
}

guideline_from_testing_standard <- function(x) {
  values <- unique(as.character(x))
  values <- values[!is.na(values) & nzchar(values)]

  if (!length(values)) {
    return(NULL)
  }

  valid_choices <- ncbi_guideline_choices()

  for (value in values) {
    upper_value <- toupper(value)
    matched <- valid_choices[toupper(valid_choices) == upper_value]
    if (length(matched)) {
      return(matched[[1]])
    }

    matched <- valid_choices[valid_choices != "missing" & grepl(
      paste0("\\b", valid_choices[valid_choices != "missing"], "\\b"),
      upper_value
    )]
    if (length(matched)) {
      return(matched[[1]])
    }
  }

  NULL
}

default_ncbi_guideline <- function(data, pheno_col) {
  if (!is.null(data) && "testing_standard" %in% names(data)) {
    from_data <- guideline_from_testing_standard(data$testing_standard)
    if (!is.null(from_data)) {
      return(from_data)
    }
  }

  from_pheno <- guideline_for_pheno(pheno_col)
  if (!is.null(from_pheno)) {
    return(from_pheno)
  }

  "missing"
}

breakpoint_version_for_pheno <- function(pheno_col, current_value = NULL) {
  version <- version_for_pheno(pheno_col)

  if (!is.null(version) && nzchar(version)) {
    return(version)
  }

  if (!is.null(current_value) && nzchar(current_value)) {
    return(current_value)
  }

  ""
}

prepare_export_data <- function(data, pheno_col, guideline = NULL) {
  guideline <- null_if_blank(guideline)
  standard <- guideline %||% guideline_for_pheno(pheno_col)

  if (is.null(standard)) {
    return(data)
  }

  data$guideline <- standard
  data
}

apply_ncbi_testing_standard <- function(data, pheno_col, guideline = NULL) {
  guideline <- null_if_blank(guideline)
  standard <- NULL

  if (!is.null(guideline) && nzchar(guideline)) {
    standard <- guideline
  } else {
    standard <- export_standard_for_pheno(pheno_col)
  }

  if (!is.null(standard) && "testing_standard" %in% names(data)) {
    data$testing_standard <- standard
  }

  data
}

apply_ncbi_export_overrides <- function(data, guideline = NULL, vendor = NULL) {
  guideline <- null_if_blank(guideline)
  vendor <- null_if_blank(vendor)

  if (!is.null(guideline) && nzchar(guideline) && "testing_standard" %in% names(data)) {
    data$testing_standard <- guideline
  }

  if (!is.null(vendor) && nzchar(vendor) && "vendor" %in% names(data)) {
    data$vendor <- vendor
  }

  data
}

apply_ebi_testing_standard <- function(data, pheno_col, guideline = NULL) {
  guideline <- null_if_blank(guideline)
  standard <- guideline

  if (is.null(standard) || !nzchar(standard)) {
    standard <- guideline_for_pheno(pheno_col)
  }

  if (!is.null(standard) && "ast_standard" %in% names(data)) {
    data$ast_standard <- standard
  }

  data
}

apply_ebi_breakpoint_version <- function(data, breakpoint_version = NULL) {
  breakpoint_version <- null_if_blank(breakpoint_version)
  if (!is.null(breakpoint_version) && nzchar(breakpoint_version) && "breakpoint_version" %in% names(data)) {
    data$breakpoint_version <- breakpoint_version
  }

  data
}

app_theme <- bs_theme(
  version = 5,
  bootswatch = "minty",
  primary = "#68349A",
  secondary = "#355c7d",
  base_font = font_google("Source Sans 3"),
  heading_font = font_google("Merriweather Sans"),
  code_font = font_google("IBM Plex Mono")
)

ui <- page_sidebar(
  title = NULL,
  theme = app_theme,
  tags$head(
    tags$style(HTML("
      .btn,
      .btn-primary,
      .btn-default,
      .btn-secondary,
      .btn-outline-secondary {
        background-color: #68349A;
        border-color: #68349A;
        color: #FFFFFF;
      }
      .btn:hover,
      .btn:focus,
      .btn:active,
      .btn-primary:hover,
      .btn-primary:focus,
      .btn-primary:active,
      .btn-default:hover,
      .btn-default:focus,
      .btn-default:active,
      .btn-secondary:hover,
      .btn-secondary:focus,
      .btn-secondary:active,
      .btn-outline-secondary:hover,
      .btn-outline-secondary:focus,
      .btn-outline-secondary:active {
        background-color: #68349A;
        border-color: #68349A;
        color: #FFFFFF;
        filter: brightness(0.95);
      }
      .nav-pills .nav-link,
      .nav-tabs .nav-link {
        color: #68349A;
      }
      .nav-pills .nav-link.active,
      .nav-tabs .nav-link.active {
        background-color: #68349A;
        border-color: #68349A;
        color: #FFFFFF;
      }
      .progress-bar,
      .shiny-file-input-progress .progress-bar {
        background-color: #A485C4 !important;
      }
      .sidebar,
      .bslib-sidebar-layout .sidebar,
      .bslib-sidebar-layout .sidebar-content,
      .bslib-sidebar-layout .sidebar > .sidebar-content,
      .bslib-sidebar-layout .sidebar .shiny-input-container,
      .bslib-sidebar-layout .sidebar .form-group,
      .bslib-sidebar-layout .sidebar .selectize-control,
      .bslib-sidebar-layout .sidebar .selectize-input,
      .bslib-sidebar-layout .sidebar .selectize-dropdown,
      .bslib-sidebar-layout .sidebar .form-control,
      .bslib-sidebar-layout .sidebar .well {
        background-color: #FFFFFF !important;
      }
    "))
  ),
  sidebar = sidebar(
    width = 360,
    style = "background-color: #FFFFFF; display: flex; flex-direction: column; height: 100%;",
    div(
      style = "background-color: #FFFFFF; padding: 10px 16px 16px 16px; margin: -1rem -1rem 1rem -1rem;",
      div(
        style = "display: flex; align-items: center; justify-content: center;",
        tags$a(
          href = "https://amrgen.org",
          target = "_blank",
          rel = "noopener noreferrer",
          tags$img(
            src = "assets/logo.png",
            alt = "AMRgen logo",
            style = "max-width: 100%; max-height: 72px; height: auto;"
          )
        )
      )
    ),
    div(
      style = "display: flex; flex-direction: column; gap: 0.5rem; flex: 1 1 auto; background-color: #FFFFFF;",
      fileInput(
        "input_file",
        "Upload antibiogram file",
        accept = c(".csv", ".tsv", ".txt", ".gz", ".xls", ".xlsx")
      ),
      selectInput(
        "input_format",
        "Input format",
        choices = c(
          "Sensititre" = "sensititre",
          "Vitek" = "vitek",
          "Phoenix" = "phoenix",
          "MicroScan" = "microscan",
          "WHONET" = "whonet",
          "NCBI browser download" = "ncbi",
          "NCBI cloud download" = "ncbi_biosample",
          "EBI browser download" = "ebi_web",
          "EBI ftp download" = "ebi_ftp"
        ),
        selected = "sensititre"
      ),
      tags$h5("Interpret resistance categories?", style = "margin-top: 0.5rem; color: #3E0B5C;"),
      checkboxInput("interpret_eucast", "Interpret using EUCAST breakpoints", FALSE),
      checkboxInput("interpret_clsi", "Interpret using CLSI breakpoints", FALSE),
      checkboxInput("interpret_ecoff", "Interpret WT/NWT using ECOFF", FALSE),
      actionButton("run_import", "Import and Summarise", class = "btn-primary"),
      uiOutput("extra_import_args_ui"),
      tags$h5("Other options", style = "margin-top: 1rem; color: #3E0B5C;"),
      tags$p(
        "Input files will normally specify species and antibiotic for each data point. If not, and your file includes data for a single species or antibiotic, you can provide the value here instead.",
        style = "font-size: 0.9rem; color: #4a4a4a; margin-bottom: 0.75rem;"
      ),
      textInput("species_override", "Species override"),
      textInput("antibiotic_override", "Antibiotic override")
    ),
    div(
      style = "margin-top: auto; padding-top: 1rem; background-color: #FFFFFF; text-align: center;",
      tags$a(
        href = "https://github.com/AMRverse/AMRgenconverter",
        target = "_blank",
        rel = "noopener noreferrer",
        style = "color: #3E0B5C; font-weight: 600; text-decoration: none;",
        "Source code"
      )
    )
  ),
  navset_card_tab(
    id = "main_tabs",
    selected = "instructions",
    nav_panel(
      "Instructions",
      value = "instructions",
      div(
        style = "padding: 0.5rem 0.25rem 0.25rem 0.25rem;",
        tags$h4("How to use this app", style = "color: #3E0B5C; margin-bottom: 1rem;"),
        tags$p(HTML("This app imports antibiogram data from a variety of input formats, and exports submission-ready files suitable for submission to <a href=\"https://www.ncbi.nlm.nih.gov/biosample/docs/antibiogram/\" target=\"_blank\" rel=\"noopener noreferrer\">NCBI</a> or <a href=\"https://www.ebi.ac.uk/amr/amr_submission_guide/\" target=\"_blank\" rel=\"noopener noreferrer\">ENA (via EBI AMR portal)</a> as BioSample data. It is powered by the <a href=\"https://amrgen.org/\" target=\"_blank\" rel=\"noopener noreferrer\">AMRgen</a> and <a href=\"https://amr-for-r.org/\" target=\"_blank\" rel=\"noopener noreferrer\">AMR for R</a> packages.")),
        tags$ol(
          tags$li("Upload your antibiogram file in a supported format. Click the 'Input format' selector to see the list of supported input formats."),
          tags$li('If you would like to re-interpret the MIC or disk measurements in the input file using latest breakpoints, check the relevant boxes under "Interpret resistance categories?"'),
          tags$li(HTML("Optionally, you can specify any additional import arguments using the fields below the Import button. See the AMRgen <a href=\"https://amrgen.org/reference/import_pheno.html\" target=\"_blank\" rel=\"noopener noreferrer\">import_pheno()</a> documentation for further guidance.")),
          tags$li("Click 'Import and Summarise' to load the data."),
          tags$li(HTML("Review the imported data and summary tables, using the tabs above to navigate. See the AMRgen <a href=\"https://amrgen.org/reference/summarise_pheno.html\" target=\"_blank\" rel=\"noopener noreferrer\">summarise_pheno()</a> documentation for help with interpreting the summary tables.")),
          tags$li(
            HTML("Open the Exports tab to convert the imported data to submission-ready formats. For further details of the exported file formats and options, see the AMRgen functions:"),
            tags$ul(
              style = "margin-top: 0.5rem;",
              tags$li(HTML("<a href=\"https://amrgen.org/reference/export_ncbi_ast.html\" target=\"_blank\" rel=\"noopener noreferrer\">export_ncbi_ast()</a>")),
              tags$li(HTML("<a href=\"https://amrgen.org/reference/export_ebi_ast.html\" target=\"_blank\" rel=\"noopener noreferrer\">export_ebi_ast()</a>"))
            )
          ),
          tags$li("Note that NCBI and ENA require valid BioSamples be included in the submission files. If needed you can upload a second file mapping sample identifiers from your input antibiogram file, to BioSample accessions, under the Exports tab before downloading. This file can be in any tabular format, you just need to indicate which column contains sample identifiers that match your input file, and which column contains the BioSample accessions to replace these with in the exported files.")
        )
      )
    ),
    nav_panel(
      "Imported data",
      value = "imported_data",
      DTOutput("imported_table")
    ),
    nav_panel(
      "Summary",
      value = "summary",
      navset_card_tab(
        id = "summary_tabs",
        nav_panel("Overview", DTOutput("summary_uniques")),
        nav_panel("Drugs", DTOutput("summary_drugs")),
        nav_panel("Details", DTOutput("summary_details")),
        nav_panel("Phenotypes", uiOutput("pheno_summary_ui"))
      )
    ),
    nav_panel(
      "Exports",
      value = "exports",
      layout_sidebar(
        sidebar = sidebar(
          width = 320,
          position = "right",
          title = "Optional BioSample mapping",
          uiOutput("biosample_mapping_help_ui"),
          fileInput(
            "biosample_file",
            "Upload BioSample mapping file",
            accept = c(".csv", ".tsv", ".txt", ".xls", ".xlsx")
          ),
          uiOutput("biosample_sample_col_ui")
        ),
        navset_card_tab(
          id = "export_tabs",
          nav_panel(
            "NCBI",
            tags$p(
              HTML('Export a TSV file suitable for submission to NCBI. See <a href="https://www.ncbi.nlm.nih.gov/biosample/docs/antibiogram/" target="_blank" rel="noopener noreferrer">NCBI BioSample Antibiograms</a> for details.'),
              style = "margin-bottom: 1rem;"
            ),
            div(
              style = "margin-bottom: 24px;",
              layout_columns(
                col_widths = c(4, 8),
                div(style = "padding-top: 31px;", downloadButton("download_ncbi", "Download NCBI TSV")),
                uiOutput("export_pheno_col_ncbi_ui")
              ),
              div(
                style = "margin-top: 16px;",
                layout_columns(
                  col_widths = c(6, 6),
                  selectInput(
                    "ncbi_guideline",
                    "Guideline",
                    choices = ncbi_guideline_choices(),
                    selected = "missing"
                  ),
                  textInput("ncbi_vendor", "Vendor (optional)")
                )
              )
            ),
            hr(style = "margin-top: 0; margin-bottom: 24px;"),
            uiOutput("ncbi_biosample_warning_ui"),
            uiOutput("ncbi_preview_ui")
          ),
          nav_panel(
            "EBI table",
            value = "ebi_table",
            tags$p(
              HTML('Export a TSV file with the fields required for ENA (for review purposes only). To prepare submission-ready JSON files for ENA use the \'EBI JSON\' tab. See <a href="https://www.ebi.ac.uk/amr/amr_submission_guide/" target="_blank" rel="noopener noreferrer">AMR Submission Guide</a> for details.'),
              style = "margin-bottom: 1rem;"
            ),
            div(
              style = "margin-bottom: 24px;",
              layout_columns(
                col_widths = c(4, 8),
                div(style = "padding-top: 31px;", downloadButton("download_ebi_tsv", "Download EBI TSV")),
                uiOutput("export_pheno_col_ebi_table_ui")
              ),
              div(
                style = "margin-top: 16px;",
                layout_columns(
                  col_widths = c(6, 6),
                  textInput("ebi_guideline", "AST standard"),
                  textInput("ebi_breakpoint_version", "Breakpoint version", value = "")
                )
              )
            ),
            hr(style = "margin-top: 0; margin-bottom: 24px;"),
            uiOutput("ebi_biosample_warning_ui"),
            uiOutput("ebi_preview_ui")
          ),
          nav_panel(
            "EBI JSON",
            value = "ebi_json",
            tags$p(
              HTML('Export a zipped directory of submission-ready JSON files for ENA. See <a href="https://www.ebi.ac.uk/amr/amr_submission_guide/" target="_blank" rel="noopener noreferrer">AMR Submission Guide</a> for details.'),
              style = "margin-bottom: 1rem;"
            ),
            div(
              style = "margin-bottom: 24px;",
              layout_columns(
                col_widths = c(4, 8),
                div(style = "padding-top: 31px;", uiOutput("download_ebi_json_top_ui")),
                uiOutput("export_pheno_col_ebi_json_ui")
              ),
              div(
                style = "margin-top: 8px; margin-bottom: 8px;",
                layout_columns(
                  col_widths = c(6, 6),
                  textInput("ebi_guideline_json", "AST standard"),
                  textInput("ebi_breakpoint_version_json", "Breakpoint version", value = "")
                )
              )
            ),
            div(
              style = "margin-top: 0; margin-bottom: 8px;",
              layout_columns(
                col_widths = c(6, 6),
                textInput("ebi_domain", "EBI domain", value = "self.ExampleDomain"),
                textInput("ebi_submission_account", "EBI submission account (optional)")
              )
            ),
            uiOutput("ebi_json_biosample_warning_ui"),
            uiOutput("ebi_json_info_ui")
          )
        )
      )
    )
  ),
  card(
    full_screen = TRUE,
    card_header("Workflow status"),
    verbatimTextOutput("status_log")
  )
)

server <- function(input, output, session) {
  rv <- reactiveValues(
    imported = NULL,
    summary = NULL,
    ncbi = NULL,
    ebi = NULL,
    status = "Upload a supported AST file, choose the source format, and run the import.",
    export_pheno_col = NULL,
    ebi_guideline = "",
    ebi_breakpoint_version = ""
  )

  biosample_mapping_info <- reactive({
    if (is.null(input$biosample_file)) {
      return(NULL)
    }

    mapping_df <- tryCatch(
      read_tabular_input(input$biosample_file),
      error = function(e) e
    )

    if (inherits(mapping_df, "error")) {
      stop(conditionMessage(mapping_df))
    }

    biosample_col <- find_biosample_column(mapping_df)
    if (is.null(biosample_col)) {
      biosample_col <- names(mapping_df)[[1]]
    }

    list(
      data = mapping_df,
      biosample_default = biosample_col
    )
  })

  output$biosample_sample_col_ui <- renderUI({
    info <- tryCatch(biosample_mapping_info(), error = function(e) e)

    if (is.null(info)) {
      return(NULL)
    }

    if (inherits(info, "error")) {
      return(tags$p(conditionMessage(info), style = "color: #b22222;"))
    }

    biosample_choices <- names(info$data)
    raw_id_default <- default_mapping_id_column(info$data, info$biosample_default)
    raw_id_choices <- names(info$data)

    tagList(
      selectInput(
        "biosample_biosample_col",
        "BioSample column",
        choices = biosample_choices,
        selected = info$biosample_default
      ),
      selectInput(
        "biosample_sample_col",
        "Raw sample ID column",
        choices = raw_id_choices,
        selected = raw_id_default %||% raw_id_choices[[1]]
      )
    )
  })

  observeEvent(input$input_file, {
    req(input$input_file$name)

    guessed_format <- guess_format_from_filename(input$input_file$name)

    if (!is.null(guessed_format)) {
      updateSelectInput(
        session,
        "input_format",
        selected = guessed_format
      )
    }
  })

  output$extra_import_args_ui <- renderUI({
    args <- format_specific_args(input$input_format)

    if (!length(args)) {
      return(NULL)
    }

    tagList(
      tags$h5("Format-specific options", style = "margin-top: 0.5rem; color: #3E0B5C;"),
      lapply(names(args), function(arg_name) {
        render_extra_arg_input(arg_name, args[[arg_name]])
      })
    )
  })

  output$export_pheno_col_ncbi_ui <- renderUI({
    choices <- phenotype_columns(rv$imported %||% data.frame())
    if (!length(choices)) {
      return(helpText("Import data to choose a phenotype column for export."))
    }

    selectInput("export_pheno_col_ncbi", "Phenotype column for exports", choices = choices, selected = rv$export_pheno_col %||% choices[[1]])
  })

  output$export_pheno_col_ebi_table_ui <- renderUI({
    choices <- phenotype_columns(rv$imported %||% data.frame())
    if (!length(choices)) {
      return(helpText("Import data to choose a phenotype column for export."))
    }

    selectInput("export_pheno_col_ebi_table", "Phenotype column for exports", choices = choices, selected = rv$export_pheno_col %||% choices[[1]])
  })

  output$export_pheno_col_ebi_json_ui <- renderUI({
    choices <- phenotype_columns(rv$imported %||% data.frame())
    if (!length(choices)) {
      return(helpText("Import data to choose a phenotype column for export."))
    }

    selectInput("export_pheno_col_ebi_json", "Phenotype column for exports", choices = choices, selected = rv$export_pheno_col %||% choices[[1]])
  })

  observeEvent(input$export_pheno_col_ncbi, {
    req(input$export_pheno_col_ncbi)
    rv$export_pheno_col <- input$export_pheno_col_ncbi
  }, ignoreInit = TRUE)

  observeEvent(input$export_pheno_col_ebi_table, {
    req(input$export_pheno_col_ebi_table)
    rv$export_pheno_col <- input$export_pheno_col_ebi_table
  }, ignoreInit = TRUE)

  observeEvent(input$export_pheno_col_ebi_json, {
    req(input$export_pheno_col_ebi_json)
    rv$export_pheno_col <- input$export_pheno_col_ebi_json
  }, ignoreInit = TRUE)

  observe({
    req(rv$export_pheno_col)
    choices <- phenotype_columns(rv$imported %||% data.frame())
    req(length(choices))

    updateSelectInput(session, "export_pheno_col_ncbi", choices = choices, selected = rv$export_pheno_col)
    updateSelectInput(session, "export_pheno_col_ebi_table", choices = choices, selected = rv$export_pheno_col)
    updateSelectInput(session, "export_pheno_col_ebi_json", choices = choices, selected = rv$export_pheno_col)
  })

  observeEvent(rv$export_pheno_col, {
    req(rv$export_pheno_col)

    default_guideline <- guideline_for_pheno(rv$export_pheno_col)
    default_version <- breakpoint_version_for_pheno(rv$export_pheno_col)
    default_ncbi_guideline_value <- default_ncbi_guideline(rv$imported, rv$export_pheno_col)

    rv$ebi_guideline <- default_guideline %||% ""
    rv$ebi_breakpoint_version <- default_version %||% ""

    updateTextInput(
      session,
      "ebi_breakpoint_version",
      value = rv$ebi_breakpoint_version
    )
    updateTextInput(
      session,
      "ebi_guideline",
      value = rv$ebi_guideline
    )

    updateTextInput(
      session,
      "ebi_breakpoint_version_json",
      value = rv$ebi_breakpoint_version
    )
    updateTextInput(
      session,
      "ebi_guideline_json",
      value = rv$ebi_guideline
    )
    updateSelectInput(
      session,
      "ncbi_guideline",
      selected = default_ncbi_guideline_value
    )
  }, ignoreInit = TRUE)

  observeEvent(input$ebi_guideline, {
    req(!is.null(input$ebi_guideline))
    rv$ebi_guideline <- input$ebi_guideline
    if (!identical(input$ebi_guideline_json, rv$ebi_guideline)) {
      updateTextInput(session, "ebi_guideline_json", value = rv$ebi_guideline)
    }
  }, ignoreInit = TRUE)

  observeEvent(input$ebi_guideline_json, {
    req(!is.null(input$ebi_guideline_json))
    rv$ebi_guideline <- input$ebi_guideline_json
    if (!identical(input$ebi_guideline, rv$ebi_guideline)) {
      updateTextInput(session, "ebi_guideline", value = rv$ebi_guideline)
    }
  }, ignoreInit = TRUE)

  observeEvent(input$ebi_breakpoint_version, {
    req(!is.null(input$ebi_breakpoint_version))
    rv$ebi_breakpoint_version <- input$ebi_breakpoint_version
    if (!identical(input$ebi_breakpoint_version_json, rv$ebi_breakpoint_version)) {
      updateTextInput(session, "ebi_breakpoint_version_json", value = rv$ebi_breakpoint_version)
    }
  }, ignoreInit = TRUE)

  observeEvent(input$ebi_breakpoint_version_json, {
    req(!is.null(input$ebi_breakpoint_version_json))
    rv$ebi_breakpoint_version <- input$ebi_breakpoint_version_json
    if (!identical(input$ebi_breakpoint_version, rv$ebi_breakpoint_version)) {
      updateTextInput(session, "ebi_breakpoint_version", value = rv$ebi_breakpoint_version)
    }
  }, ignoreInit = TRUE)

  ebi_guideline_value <- reactive({
    if (!is.null(input$ebi_guideline_json) && nzchar(input$ebi_guideline_json)) {
      return(input$ebi_guideline_json)
    }
    if (!is.null(input$ebi_guideline) && nzchar(input$ebi_guideline)) {
      return(input$ebi_guideline)
    }
    rv$ebi_guideline %||% ""
  })

  ebi_breakpoint_version_value <- reactive({
    if (!is.null(input$ebi_breakpoint_version_json) && nzchar(input$ebi_breakpoint_version_json)) {
      return(input$ebi_breakpoint_version_json)
    }
    if (!is.null(input$ebi_breakpoint_version) && nzchar(input$ebi_breakpoint_version)) {
      return(input$ebi_breakpoint_version)
    }
    rv$ebi_breakpoint_version %||% ""
  })

  output$biosample_mapping_help_ui <- renderUI({
    if (identical(input$export_tabs, "ebi_table") || identical(input$export_tabs, "ebi_json")) {
      return(tags$p(
        'EBI submission files require valid ENA BioSample identifiers for each sample (e.g. "SAMExxxx" or "ERSxxx"). If your antibiogram file contains sample identifiers other than ENA BioSample accessions, you should upload a file mapping those sample identifiers to valid BioSamples.',
        style = "margin-bottom: 0.75rem;"
      ))
    }

    tags$p(
      'NCBI submission files require valid NCBI BioSample identifiers for each sample (e.g. "SAMNxxxx" or "SRSxxx"). If your antibiogram file contains sample identifiers other than NCBI BioSample accessions, you should upload a file mapping those sample identifiers to valid BioSamples.',
      style = "margin-bottom: 0.75rem;"
    )
  })

  output$ncbi_biosample_warning_ui <- renderUI({
    req(rv$imported)

    ids <- export_source_data()$id

    if (all_biosample_ids_match(ids, c("SAME", "ERS"))) {
      return(tags$p(
        "WARNING: Sample identifiers should be valid NCBI BioSample accessions. The BioSamples currently loaded appear to be EBI BioSamples. To export to EBI submission format, click the 'EBI table' or 'EBI JSON' tab instead.",
        style = "color: #68349A; font-weight: 600; margin-bottom: 0.75rem;"
      ))
    }

    if (has_invalid_biosample_ids(ids, c("SAMN", "SRS"))) {
      return(tags$p(
        "WARNING: Sample identifiers should be valid NCBI BioSample accessions. Upload a BioSample mapping file before exporting.",
        style = "color: #68349A; font-weight: 600; margin-bottom: 0.75rem;"
      ))
    }

    NULL
  })

  output$ebi_biosample_warning_ui <- renderUI({
    req(rv$imported)

    ids <- export_source_data()$id

    if (all_biosample_ids_match(ids, c("SAMN", "SRS"))) {
      return(tags$p(
        "WARNING: Sample identifiers should be valid ENA BioSample accessions. The BioSamples currently loaded appear to be NCBI BioSamples. To export to NCBI submission format, click the 'NCBI' tab instead.",
        style = "color: #68349A; font-weight: 600; margin-bottom: 0.75rem;"
      ))
    }

    if (has_invalid_biosample_ids(ids, c("SAME", "ERS"))) {
      return(tags$p(
        "WARNING: Sample identifiers should be valid ENA BioSample accessions. Upload a BioSample mapping file before exporting.",
        style = "color: #68349A; font-weight: 600; margin-bottom: 0.75rem;"
      ))
    }

    NULL
  })

  output$ebi_json_biosample_warning_ui <- renderUI({
    req(rv$imported)

    ids <- export_source_data()$id

    if (all_biosample_ids_match(ids, c("SAMN", "SRS"))) {
      return(tags$p(
        "WARNING: Sample identifiers should be valid ENA BioSample accessions. The BioSamples currently loaded appear to be NCBI BioSamples. To export to NCBI submission format, click the 'NCBI' tab instead.",
        style = "color: #68349A; font-weight: 600; margin-bottom: 0.75rem;"
      ))
    }

    if (has_invalid_biosample_ids(ids, c("SAME", "ERS"))) {
      return(tags$p(
        "WARNING: Sample identifiers should be valid ENA BioSample accessions. Upload a BioSample mapping file before exporting.",
        style = "color: #68349A; font-weight: 600; margin-bottom: 0.75rem;"
      ))
    }

    NULL
  })

  output$ncbi_preview_ui <- renderUI({
    if (is.null(rv$imported)) {
      return(p("Import data to preview and download NCBI exports."))
    }

    DTOutput("ncbi_preview")
  })

  output$ebi_preview_ui <- renderUI({
    if (is.null(rv$imported)) {
      return(p("Import data to preview and download EBI exports."))
    }

    DTOutput("ebi_preview")
  })

  output$ebi_json_info_ui <- renderUI({
    if (is.null(rv$imported)) {
      return(p("Import data to prepare EBI JSON exports."))
    }

    tagList(
      if (!nzchar(input$ebi_domain)) {
        p("Enter the domain to enable the JSON download.")
      }
    )
  })

  output$download_ebi_json_top_ui <- renderUI({
    if (is.null(rv$imported)) {
      return(NULL)
    }

    if (nzchar(input$ebi_domain)) {
      return(downloadButton("download_ebi_json", "Download EBI JSON zip"))
    }

    tags$span("Complete fields below", style = "color: #6c757d;")
  })

  observeEvent(input$run_import, {
    req(input$input_file)
    prepared_input <- prepare_import_input(input$input_file, format = input$input_format)

    import_args <- list(
      input = prepared_input$input,
      format = input$input_format,
      interpret_eucast = isTRUE(input$interpret_eucast),
      interpret_clsi = isTRUE(input$interpret_clsi),
      interpret_ecoff = isTRUE(input$interpret_ecoff),
      species = if (nzchar(input$species_override)) input$species_override else NULL,
      ab = if (nzchar(input$antibiotic_override)) input$antibiotic_override else NULL,
      source = if (nzchar(input$source_override)) input$source_override else NULL
    )
    import_args <- c(import_args, collect_extra_import_args(input, input$input_format))

    result <- tryCatch(
      capture_conditions(do.call(AMRgen::import_pheno, import_args)),
      error = function(e) e
    )

    if (inherits(result, "error")) {
      rv$imported <- NULL
      rv$summary <- NULL
      rv$ncbi <- NULL
      rv$ebi <- NULL
      rv$status <- paste("Import failed:", conditionMessage(result))
      return()
    }

    rv$imported <- result$value
    rv$export_pheno_col <- default_export_pheno_col(
      phenotype_columns(rv$imported),
      interpret_eucast = isTRUE(input$interpret_eucast),
      interpret_clsi = isTRUE(input$interpret_clsi)
    )

    summary_result <- tryCatch(
      capture_conditions(AMRgen::summarise_pheno(rv$imported)),
      error = function(e) e
    )

    if (inherits(summary_result, "error")) {
      rv$summary <- NULL
      rv$status <- paste(
        c(
          paste("Imported", nrow(rv$imported), "rows."),
          if (length(result$warnings)) paste("Import warnings:", paste(result$warnings, collapse = " | ")) else NULL,
          paste("Summary failed:", conditionMessage(summary_result))
        ),
        collapse = "\n"
      )
      return()
    }

    rv$summary <- summary_result$value

    rv$status <- paste(
      c(
        paste(
          "Imported", nrow(rv$imported), "rows from",
          shQuote(input$input_file$name),
          "using", format_choice_label(input$input_format), "format."
        ),
        prepared_input$note,
        if (length(result$warnings)) paste("Import warnings:", paste(result$warnings, collapse = " | ")) else NULL,
        if (length(summary_result$warnings)) paste("Summary warnings:", paste(summary_result$warnings, collapse = " | ")) else NULL
      ),
      collapse = "\n"
    )

    bslib::nav_select("main_tabs", "imported_data")
  })

  export_source_data <- reactive({
    req(rv$imported)

    data <- rv$imported
    info <- tryCatch(biosample_mapping_info(), error = function(e) e)

    if (inherits(info, "error")) {
      stop(conditionMessage(info))
    }

    if (is.null(info)) {
      return(data)
    }

    biosample_col <- input$biosample_biosample_col %||% info$biosample_default
    sample_id_col <- input$biosample_sample_col

    if (is.null(biosample_col) || !nzchar(biosample_col)) {
      stop("Select the BioSample column from the BioSample mapping file.")
    }

    if (is.null(sample_id_col) || !nzchar(sample_id_col)) {
      stop("Select the raw sample ID column from the BioSample mapping file.")
    }

    if (identical(biosample_col, sample_id_col)) {
      stop("The BioSample column and raw sample ID column must be different.")
    }

    apply_biosample_mapping(
      data = data,
      mapping_df = info$data,
      sample_id_col = sample_id_col,
      biosample_col = biosample_col
    )
  })

  export_tables <- reactive({
    req(rv$export_pheno_col)
    export_data <- prepare_export_data(
      export_source_data(),
      rv$export_pheno_col,
      guideline = input$ncbi_guideline
    )

    ncbi_call_args <- filter_supported_args(
      list(
        data = export_data,
        pheno_col = rv$export_pheno_col,
        guideline = input$ncbi_guideline,
        vendor = null_if_blank(input$ncbi_vendor)
      ),
      "export_ncbi_ast"
    )

    ncbi_result <- capture_conditions(
      do.call(AMRgen::export_ncbi_ast, ncbi_call_args)
    )
    ncbi_result$value <- apply_ncbi_testing_standard(
      ncbi_result$value,
      rv$export_pheno_col,
      guideline = input$ncbi_guideline
    )
    ncbi_result$value <- apply_ncbi_export_overrides(
      ncbi_result$value,
      guideline = input$ncbi_guideline,
      vendor = input$ncbi_vendor
    )

    ebi_args <- filter_supported_args(list(
      data = prepare_export_data(
        export_source_data(),
        rv$export_pheno_col,
        guideline = ebi_guideline_value()
      ),
      pheno_col = rv$export_pheno_col,
      guideline = ebi_guideline_value(),
      breakpoint_version = ebi_breakpoint_version_value()
    ), "export_ebi_ast")

    ebi_result <- capture_conditions(do.call(AMRgen::export_ebi_ast, ebi_args))
    ebi_result$value <- apply_ebi_testing_standard(
      ebi_result$value,
      rv$export_pheno_col,
      guideline = ebi_guideline_value()
    )
    ebi_result$value <- apply_ebi_breakpoint_version(
      ebi_result$value,
      breakpoint_version = ebi_breakpoint_version_value()
    )

    rv$ncbi <- ncbi_result$value
    rv$ebi <- ebi_result$value

    list(
      ncbi = ncbi_result,
      ebi = ebi_result
    )
  })

  output$rows_imported <- renderText({
    if (is.null(rv$imported)) "0" else format(nrow(rv$imported), big.mark = ",")
  })

  output$pheno_col_count <- renderText({
    if (is.null(rv$imported)) {
      "0"
    } else {
      as.character(length(phenotype_columns(rv$imported)))
    }
  })

  output$selected_format <- renderText({
    format_choice_label(input$input_format)
  })

  output$status_log <- renderText({
    rv$status
  })

  output$imported_table <- renderDT({
    req(rv$imported)
    datatable(rv$imported, filter = "top", options = list(scrollX = TRUE, pageLength = 10))
  })

  output$summary_uniques <- renderDT({
    req(rv$summary)
    datatable(rv$summary$uniques, options = list(dom = "t"))
  })

  output$summary_drugs <- renderDT({
    req(rv$summary)
    datatable(rv$summary$drugs, options = list(scrollX = TRUE, pageLength = 10))
  })

  output$summary_details <- renderDT({
    req(rv$summary)
    datatable(rv$summary$details, options = list(scrollX = TRUE, pageLength = 10))
  })

  output$pheno_summary_ui <- renderUI({
    req(rv$summary)
    if (!length(rv$summary$pheno_counts_list)) {
      return(p("No phenotype count tables were generated for the imported dataset."))
    }

    tabs <- lapply(names(rv$summary$pheno_counts_list), function(name) {
      nav_panel(name, DTOutput(paste0("pheno_", name)))
    })
    do.call(navset_card_tab, tabs)
  })

  observe({
    req(rv$summary)
    pheno_tables <- rv$summary$pheno_counts_list

    for (name in names(pheno_tables)) {
      local({
        output_id <- paste0("pheno_", name)
        table_data <- pheno_tables[[name]]

        output[[output_id]] <- renderDT({
          datatable(table_data, options = list(scrollX = TRUE, pageLength = 10))
        })
      })
    }
  })

  output$ncbi_preview <- renderDT({
    tables <- export_tables()
    datatable(
      tables$ncbi$value,
      options = list(scrollX = TRUE, pageLength = 10, dom = "tip")
    )
  })

  output$ebi_preview <- renderDT({
    tables <- export_tables()
    datatable(
      tables$ebi$value,
      options = list(scrollX = TRUE, pageLength = 10, dom = "tip")
    )
  })

  output$download_ncbi <- downloadHandler(
    filename = function() {
      paste0(tools::file_path_sans_ext(input$input_file$name), "_ncbi.tsv")
    },
    content = function(file) {
      export_data <- prepare_export_data(
        export_source_data(),
        rv$export_pheno_col,
        guideline = input$ncbi_guideline
      )
      ncbi_table <- do.call(
        AMRgen::export_ncbi_ast,
        filter_supported_args(
          list(
            data = export_data,
            pheno_col = rv$export_pheno_col,
            guideline = input$ncbi_guideline,
            vendor = null_if_blank(input$ncbi_vendor)
          ),
          "export_ncbi_ast"
        )
      )

      ncbi_table <- apply_ncbi_testing_standard(
        ncbi_table,
        rv$export_pheno_col,
        guideline = input$ncbi_guideline
      )
      ncbi_table <- apply_ncbi_export_overrides(
        ncbi_table,
        guideline = input$ncbi_guideline,
        vendor = input$ncbi_vendor
      )

      utils::write.table(
        ncbi_table,
        file = file,
        sep = "\t",
        row.names = FALSE,
        quote = TRUE,
        na = ""
      )
    }
  )

  output$download_ebi_tsv <- downloadHandler(
    filename = function() {
      paste0(tools::file_path_sans_ext(input$input_file$name), "_ebi.tsv")
    },
    content = function(file) {
      export_data <- prepare_export_data(
        export_source_data(),
        rv$export_pheno_col,
        guideline = ebi_guideline_value()
      )
      ebi_table <- do.call(
        AMRgen::export_ebi_ast,
        filter_supported_args(list(
          data = export_data,
          pheno_col = rv$export_pheno_col,
          guideline = ebi_guideline_value(),
          breakpoint_version = ebi_breakpoint_version_value()
        ), "export_ebi_ast")
      )
      ebi_table <- apply_ebi_testing_standard(
        ebi_table,
        rv$export_pheno_col,
        guideline = ebi_guideline_value()
      )
      ebi_table <- apply_ebi_breakpoint_version(
        ebi_table,
        breakpoint_version = ebi_breakpoint_version_value()
      )

      utils::write.table(
        ebi_table,
        file = file,
        sep = "\t",
        row.names = FALSE,
        quote = TRUE,
        na = ""
      )
    }
  )

  output$download_ebi_json <- downloadHandler(
    filename = function() {
      paste0(tools::file_path_sans_ext(input$input_file$name), "_ebi_json.zip")
    },
    contentType = "application/zip",
    content = function(file) {
      tryCatch({
        req(rv$export_pheno_col)
        if (!nzchar(input$ebi_domain)) {
          stop("Please provide an EBI domain before downloading JSON.")
        }

        export_data <- prepare_export_data(
          export_source_data(),
          rv$export_pheno_col,
          guideline = ebi_guideline_value()
        )

        tmp_dir <- tempfile("ebi_json_")
        dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

        do.call(
          AMRgen::export_ebi_ast,
          filter_supported_args(list(
            data = export_data,
            pheno_col = rv$export_pheno_col,
            guideline = null_if_blank(ebi_guideline_value()),
            breakpoint_version = null_if_blank(ebi_breakpoint_version_value()),
            submission_account = if (is.null(input$ebi_submission_account)) "" else input$ebi_submission_account,
            domain = input$ebi_domain,
            output_dir = tmp_dir
          ), "export_ebi_ast")
        )

        json_files <- list.files(tmp_dir, full.names = TRUE)
        if (!length(json_files)) {
          stop("No EBI JSON files were created.")
        }

        tmp_zip <- tempfile(fileext = ".zip")
        zip::zipr(
          zipfile = tmp_zip,
          files = json_files,
          root = tmp_dir
        )
        ok <- file.copy(tmp_zip, file, overwrite = TRUE)
        if (!ok) {
          stop("Failed to prepare the EBI JSON zip file for download.")
        }
      }, error = function(e) {
        friendly_message <- format_amrgen_error(conditionMessage(e))
        rv$status <- paste("EBI JSON export failed:", friendly_message)
        stop(friendly_message)
      })
    }
  )
}

shinyApp(ui, server)
