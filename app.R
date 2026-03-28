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

prepare_import_input <- function(file_info) {
  ext <- tolower(tools::file_ext(file_info$name))

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

breakpoint_version_for_pheno <- function(pheno_col, current_value = NULL) {
  standard <- export_standard_for_pheno(pheno_col)

  if (!is.null(standard)) {
    return(standard)
  }

  if (!is.null(current_value) && nzchar(current_value)) {
    return(current_value)
  }

  "EUCAST 2024"
}

prepare_export_data <- function(data, pheno_col) {
  standard <- export_standard_for_pheno(pheno_col)

  if (is.null(standard)) {
    return(data)
  }

  data$guideline <- standard
  data
}

apply_ncbi_testing_standard <- function(data, pheno_col) {
  standard <- export_standard_for_pheno(pheno_col)

  if (!is.null(standard) && "testing_standard" %in% names(data)) {
    data$testing_standard <- standard
  }

  data
}

apply_ebi_testing_standard <- function(data, pheno_col, breakpoint_version = NULL) {
  standard <- breakpoint_version

  if (is.null(standard) || !nzchar(standard)) {
    standard <- export_standard_for_pheno(pheno_col)
  }

  if (!is.null(standard) && "ast_standard" %in% names(data)) {
    data$ast_standard <- standard
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
      checkboxInput("interpret_eucast", "Interpret using EUCAST breakpoints", FALSE),
      checkboxInput("interpret_clsi", "Interpret using CLSI breakpoints", FALSE),
      checkboxInput("interpret_ecoff", "Interpret WT/NWT using ECOFF", FALSE),
      actionButton("run_import", "Import and Summarise", class = "btn-primary"),
      tags$h5("Other options", style = "margin-top: 1rem; color: #3E0B5C;"),
      textInput("species_override", "Species override"),
      textInput("antibiotic_override", "Antibiotic override"),
      textInput("source_override", "Source override")
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
    nav_panel(
      "Imported data",
      DTOutput("imported_table")
    ),
    nav_panel(
      "Summary",
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
      navset_card_tab(
        id = "export_tabs",
        nav_panel(
          "NCBI",
          div(
            style = "margin-bottom: 20px;",
            layout_columns(
              col_widths = c(4, 8),
              div(style = "padding-top: 31px;", downloadButton("download_ncbi", "Download NCBI TSV")),
              uiOutput("export_pheno_col_ncbi_ui")
            )
          ),
          hr(style = "margin-top: 0; margin-bottom: 24px;"),
          uiOutput("ncbi_preview_ui")
        ),
        nav_panel(
          "EBI table",
          value = "ebi_table",
          div(
            style = "margin-bottom: 20px;",
            layout_columns(
              col_widths = c(3, 4, 5),
              div(style = "padding-top: 31px;", downloadButton("download_ebi_tsv", "Download EBI TSV")),
              uiOutput("export_pheno_col_ebi_table_ui"),
              textInput("ebi_breakpoint_version", "EBI breakpoint version", value = "EUCAST 2024")
            )
          ),
          hr(style = "margin-top: 0; margin-bottom: 24px;"),
          uiOutput("ebi_preview_ui")
        ),
        nav_panel(
          "EBI JSON",
          value = "ebi_json",
          div(
            style = "margin-bottom: 20px;",
            layout_columns(
              col_widths = c(3, 4, 5),
              div(style = "padding-top: 31px;", uiOutput("download_ebi_json_top_ui")),
              uiOutput("export_pheno_col_ebi_json_ui"),
              textInput("ebi_breakpoint_version_json", "EBI breakpoint version", value = "EUCAST 2024")
            )
          ),
          hr(style = "margin-top: 0; margin-bottom: 24px;"),
          textInput("ebi_submission_account", "EBI submission account (optional)"),
          textInput("ebi_domain", "EBI domain", value = "self.ExampleDomain"),
          helpText("EBI JSON export requires the breakpoint version and domain. Submission account is optional."),
          uiOutput("ebi_json_info_ui")
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
    export_pheno_col = NULL
  )

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

    updateTextInput(
      session,
      "ebi_breakpoint_version",
      value = breakpoint_version_for_pheno(
        rv$export_pheno_col,
        current_value = input$ebi_breakpoint_version
      )
    )

    updateTextInput(
      session,
      "ebi_breakpoint_version_json",
      value = breakpoint_version_for_pheno(
        rv$export_pheno_col,
        current_value = input$ebi_breakpoint_version_json
      )
    )
  }, ignoreInit = TRUE)

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
      if (!nzchar(input$ebi_breakpoint_version_json) ||
          !nzchar(input$ebi_domain)) {
        p("Enter the breakpoint version and domain to enable the JSON download.")
      },
      p("The files are generated from the selected phenotype export column and the EBI metadata shown in this tab.")
    )
  })

  output$download_ebi_json_top_ui <- renderUI({
    if (is.null(rv$imported)) {
      return(NULL)
    }

    if (nzchar(input$ebi_breakpoint_version_json) &&
        nzchar(input$ebi_domain)) {
      return(downloadButton("download_ebi_json", "Download EBI JSON zip"))
    }

    tags$span("Complete fields below", style = "color: #6c757d;")
  })

  observeEvent(input$run_import, {
    req(input$input_file)
    prepared_input <- prepare_import_input(input$input_file)

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
  })

  export_tables <- reactive({
    req(rv$imported, rv$export_pheno_col)
    export_data <- prepare_export_data(rv$imported, rv$export_pheno_col)

    ncbi_result <- capture_conditions(
      AMRgen::export_ncbi_ast(
        export_data,
        pheno_col = rv$export_pheno_col
      )
    )
    ncbi_result$value <- apply_ncbi_testing_standard(
      ncbi_result$value,
      rv$export_pheno_col
    )

    ebi_args <- list(
      data = export_data,
      pheno_col = rv$export_pheno_col,
      breakpoint_version = input$ebi_breakpoint_version
    )

    ebi_result <- capture_conditions(do.call(AMRgen::export_ebi_ast, ebi_args))
    ebi_result$value <- apply_ebi_testing_standard(
      ebi_result$value,
      rv$export_pheno_col,
      breakpoint_version = input$ebi_breakpoint_version
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
      export_data <- prepare_export_data(rv$imported, rv$export_pheno_col)
      ncbi_table <- AMRgen::export_ncbi_ast(
        export_data,
        pheno_col = rv$export_pheno_col
      )

      ncbi_table <- apply_ncbi_testing_standard(
        ncbi_table,
        rv$export_pheno_col
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
      export_data <- prepare_export_data(rv$imported, rv$export_pheno_col)
      ebi_table <- do.call(
        AMRgen::export_ebi_ast,
        list(
          data = export_data,
          pheno_col = rv$export_pheno_col,
          breakpoint_version = input$ebi_breakpoint_version
        )
      )
      ebi_table <- apply_ebi_testing_standard(
        ebi_table,
        rv$export_pheno_col,
        breakpoint_version = input$ebi_breakpoint_version
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
      req(rv$imported, rv$export_pheno_col)
      if (!nzchar(input$ebi_breakpoint_version_json)) {
        stop("Please provide an EBI breakpoint version before downloading JSON.")
      }
      if (!nzchar(input$ebi_domain)) {
        stop("Please provide an EBI domain before downloading JSON.")
      }

      export_data <- prepare_export_data(rv$imported, rv$export_pheno_col)

      tmp_dir <- tempfile("ebi_json_")
      dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

      do.call(
        AMRgen::export_ebi_ast,
        list(
          data = export_data,
          pheno_col = rv$export_pheno_col,
          breakpoint_version = input$ebi_breakpoint_version_json,
          submission_account = input$ebi_submission_account,
          domain = input$ebi_domain,
          output_dir = tmp_dir
        )
      )

      json_files <- list.files(tmp_dir, full.names = TRUE)
      if (!length(json_files)) {
        stop("No EBI JSON files were created.")
      }

      zip::zipr(
        zipfile = file,
        files = json_files,
        root = tmp_dir
      )
    }
  )
}

shinyApp(ui, server)
