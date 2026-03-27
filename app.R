library(shiny)
library(bslib)
library(DT)
library(AMRgen)
library(readxl)

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

format_choice_label <- function(x) {
  switch(
    x,
    ebi = "EBI antibiogram",
    ebi_web = "EBI web export",
    ebi_ftp = "EBI FTP export",
    ncbi = "NCBI antibiogram",
    vitek = "Vitek",
    phoenix = "Phoenix",
    microscan = "MicroScan",
    sensititre = "Sensititre",
    whonet = "WHONET",
    x
  )
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
  title = "shinyAMR",
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
    "))
  ),
  sidebar = sidebar(
    width = 360,
    fileInput(
      "input_file",
      "Upload antibiogram file",
      accept = c(".csv", ".tsv", ".txt", ".gz", ".xls", ".xlsx")
    ),
    selectInput(
      "input_format",
      "Input format",
      choices = c(
        "EBI antibiogram" = "ebi",
        "EBI web export" = "ebi_web",
        "EBI FTP export" = "ebi_ftp",
        "NCBI antibiogram" = "ncbi",
        "Vitek" = "vitek",
        "Phoenix" = "phoenix",
        "MicroScan" = "microscan",
        "Sensititre" = "sensititre",
        "WHONET" = "whonet"
      ),
      selected = "ncbi"
    ),
    checkboxInput("interpret_eucast", "Interpret using EUCAST breakpoints", FALSE),
    checkboxInput("interpret_clsi", "Interpret using CLSI breakpoints", FALSE),
    checkboxInput("interpret_ecoff", "Interpret WT/NWT using ECOFF", FALSE),
    textInput("species_override", "Species override (optional)"),
    textInput("antibiotic_override", "Antibiotic override (optional)"),
    textInput("source_override", "Source override (optional)"),
    actionButton("run_import", "Import and Summarise", class = "btn-primary"),
    hr(),
    uiOutput("pheno_col_ui"),
    textInput("ebi_breakpoint_version", "EBI breakpoint version", value = "EUCAST 2024"),
    textInput("ebi_submission_account", "EBI submission account (optional)"),
    textInput("ebi_domain", "EBI domain", value = "self.ExampleDomain"),
    helpText(
      "EBI JSON download is enabled when breakpoint version, submission account, and domain are provided."
    )
  ),
  layout_column_wrap(
    width = 1 / 3,
    value_box(
      title = "Rows imported",
      value = textOutput("rows_imported"),
      theme = "primary"
    ),
    value_box(
      title = "Phenotype columns",
      value = textOutput("pheno_col_count"),
      theme = "secondary"
    ),
    value_box(
      title = "Selected format",
      value = textOutput("selected_format"),
      theme = "success"
    )
  ),
  card(
    full_screen = TRUE,
    card_header("Workflow status"),
    verbatimTextOutput("status_log")
  ),
  navset_card_tab(
    id = "main_tabs",
    title = "Results",
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
      p("Preview the converted outputs and download them below."),
      h5("NCBI preview"),
      DTOutput("ncbi_preview"),
      br(),
      downloadButton("download_ncbi", "Download NCBI TSV"),
      hr(),
      h5("EBI preview"),
      DTOutput("ebi_preview"),
      br(),
      downloadButton("download_ebi_tsv", "Download EBI TSV"),
      br(),
      br(),
      downloadButton("download_ebi_json", "Download EBI JSON zip")
    )
  )
)

server <- function(input, output, session) {
  rv <- reactiveValues(
    imported = NULL,
    summary = NULL,
    ncbi = NULL,
    ebi = NULL,
    status = "Upload a supported AST file, choose the source format, and run the import."
  )

  output$pheno_col_ui <- renderUI({
    choices <- phenotype_columns(rv$imported %||% data.frame())
    if (!length(choices)) {
      return(helpText("Import data to choose a phenotype column for export."))
    }

    selectInput(
      "export_pheno_col",
      "Phenotype column for exports",
      choices = choices,
      selected = if ("pheno_provided" %in% choices) "pheno_provided" else choices[[1]]
    )
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
    req(rv$imported, input$export_pheno_col)

    ncbi_result <- capture_conditions(
      AMRgen::export_ncbi_ast(
        rv$imported,
        pheno_col = input$export_pheno_col
      )
    )

    ebi_args <- list(
      data = rv$imported,
      pheno_col = input$export_pheno_col,
      breakpoint_version = input$ebi_breakpoint_version
    )

    ebi_result <- capture_conditions(do.call(AMRgen::export_ebi_ast, ebi_args))

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
    datatable(tables$ncbi$value, options = list(scrollX = TRUE, pageLength = 10))
  })

  output$ebi_preview <- renderDT({
    tables <- export_tables()
    datatable(tables$ebi$value, options = list(scrollX = TRUE, pageLength = 10))
  })

  output$download_ncbi <- downloadHandler(
    filename = function() {
      paste0(tools::file_path_sans_ext(input$input_file$name), "_ncbi.tsv")
    },
    content = function(file) {
      AMRgen::export_ncbi_ast(
        rv$imported,
        file = file,
        overwrite = TRUE,
        pheno_col = input$export_pheno_col
      )
    }
  )

  output$download_ebi_tsv <- downloadHandler(
    filename = function() {
      paste0(tools::file_path_sans_ext(input$input_file$name), "_ebi.tsv")
    },
    content = function(file) {
      ebi_table <- do.call(
        AMRgen::export_ebi_ast,
        list(
          data = rv$imported,
          pheno_col = input$export_pheno_col,
          breakpoint_version = input$ebi_breakpoint_version
        )
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
    content = function(file) {
      req(
        nzchar(input$ebi_breakpoint_version),
        nzchar(input$ebi_submission_account),
        nzchar(input$ebi_domain)
      )

      tmp_dir <- tempfile("ebi_json_")
      dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

      do.call(
        AMRgen::export_ebi_ast,
        list(
          data = rv$imported,
          pheno_col = input$export_pheno_col,
          breakpoint_version = input$ebi_breakpoint_version,
          submission_account = input$ebi_submission_account,
          domain = input$ebi_domain,
          output_dir = tmp_dir
        )
      )

      json_files <- list.files(tmp_dir, full.names = TRUE)
      if (!length(json_files)) {
        stop("No EBI JSON files were created.")
      }

      old_wd <- getwd()
      on.exit(setwd(old_wd), add = TRUE)
      setwd(tmp_dir)
      utils::zip(
        zipfile = file,
        files = basename(json_files)
      )
    }
  )
}

shinyApp(ui, server)
