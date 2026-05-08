library(shiny)
library(DT)
library(readr)
library(dplyr)

# ---------------------------------------------------------------------------
# Locate the project root by searching upward for the single input table
# (same convention as before, but pointed at the new combined table).
# ---------------------------------------------------------------------------
INPUT_REL <- "actinopteriigy_cne_final_table.tsv"

find_project_root <- function(target = INPUT_REL, max_up = 6) {
  p <- tryCatch(normalizePath(getwd(), winslash = "/"), error = function(e) {
    getwd()
  })
  for (i in 0:max_up) {
    candidate <- file.path(p, target)
    if (file.exists(candidate)) {
      return(p)
    }
    parent <- dirname(p)
    if (identical(parent, p)) {
      break
    }
    p <- parent
  }
  NULL
}

proj_root <- find_project_root()
if (!is.null(proj_root)) {
  if (
    !identical(
      normalizePath(getwd(), winslash = "/"),
      normalizePath(proj_root, winslash = "/")
    )
  ) {
    message("Setting working directory to project root: ", proj_root)
    setwd(proj_root)
  }
} else {
  message(
    "Warning: could not locate '",
    INPUT_REL,
    "' in parent paths; leaving working directory as ",
    getwd()
  )
}

# ---------------------------------------------------------------------------
# Single-table input. Everything the app needs (including the four binary
# upset-membership columns) is already in this file -- no side merges.
# ---------------------------------------------------------------------------
fullData <- tryCatch(
  readr::read_tsv(INPUT_REL, show_col_types = FALSE),
  error = function(e) {
    message(sprintf("Failed to read %s: %s", INPUT_REL, e$message))
    NULL
  }
)

# Binary columns produced by cne_analysis.R. Order = display order in sidebar.
BIN_COLS <- list(
  in_atac_peak = "ATAC peak overlap",
  nearby_gene_active = "Nearby gene active",
  overlaps_chan_enhancer = "Chan enhancer overlap",
  overlaps_yuesong_cne = "YueSong CNE overlap"
)

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Actinopterygii-specific CNE browser"),
  sidebarLayout(
    sidebarPanel(
      uiOutput("chrom_ui"),
      uiOutput("phastcons_ui"),
      uiOutput("width_ui"),
      tags$hr(),
      tags$h4("Set-membership filters"),
      tags$p(tags$em(
        "Each filter is independent; filters are combined with AND."
      )),
      uiOutput("binary_filters_ui"),
      tags$hr(),
      textInput("gene_search", "Gene / annotation contains:", value = ""),
      actionButton("reset", "Reset filters"),
      br(),
      br(),
      downloadButton("download_filtered", "Download filtered TSV")
    ),
    mainPanel(
      DT::dataTableOutput("table"),
      br(),
      verbatimTextOutput("stats")
    )
  )
)

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------
server <- function(input, output, session) {
  # ---- reset everything to "no filter" ----
  observeEvent(input$reset, {
    if (is.null(fullData)) {
      return()
    }
    try(updateSelectInput(session, "chrom", selected = "All"), silent = TRUE)

    if ("phastcons" %in% colnames(fullData)) {
      ph_min <- min(fullData$phastcons, na.rm = TRUE)
      ph_max <- max(fullData$phastcons, na.rm = TRUE)
      updateSliderInput(session, "phastcons", value = c(ph_min, ph_max))
    }
    if ("width" %in% colnames(fullData)) {
      wd_min <- min(fullData$width, na.rm = TRUE)
      wd_max <- max(fullData$width, na.rm = TRUE)
      updateSliderInput(session, "width", value = c(wd_min, wd_max))
    }
    for (col in names(BIN_COLS)) {
      if (col %in% colnames(fullData)) {
        try(
          updateRadioButtons(session, paste0("flt_", col), selected = "all"),
          silent = TRUE
        )
      }
    }
    updateTextInput(session, "gene_search", value = "")
  })

  # ---- dynamic UI bits ----
  output$chrom_ui <- renderUI({
    if (is.null(fullData)) {
      return(tags$div("Data not loaded."))
    }
    chrs <- unique(as.character(fullData$seqnames))
    selectInput(
      "chrom",
      "Chromosome / contig:",
      choices = c("All", sort(chrs)),
      selected = "All"
    )
  })

  output$phastcons_ui <- renderUI({
    if (is.null(fullData) || !"phastcons" %in% colnames(fullData)) {
      return(NULL)
    }
    minv <- min(fullData$phastcons, na.rm = TRUE)
    maxv <- max(fullData$phastcons, na.rm = TRUE)
    sliderInput(
      "phastcons",
      "phastCons score range:",
      min = minv,
      max = maxv,
      value = c(minv, maxv)
    )
  })

  output$width_ui <- renderUI({
    if (is.null(fullData) || !"width" %in% colnames(fullData)) {
      return(NULL)
    }
    minw <- min(fullData$width, na.rm = TRUE)
    maxw <- max(fullData$width, na.rm = TRUE)
    sliderInput(
      "width",
      "CNE width range (bp):",
      min = minw,
      max = maxw,
      value = c(minw, maxw)
    )
  })

  # One independent radio-button group per binary column.
  # Three states: All / With / Without. They combine with AND.
  output$binary_filters_ui <- renderUI({
    if (is.null(fullData)) {
      return(NULL)
    }
    tagList(
      lapply(names(BIN_COLS), function(col) {
        if (!col %in% colnames(fullData)) {
          return(NULL)
        }
        radioButtons(
          inputId = paste0("flt_", col),
          label = BIN_COLS[[col]],
          choices = c("All" = "all", "With" = "with", "Without" = "without"),
          selected = "all",
          inline = TRUE
        )
      })
    )
  })

  # ---- core filtering pipeline ----
  filtered <- reactive({
    df <- fullData
    if (is.null(df)) {
      return(NULL)
    }

    if (!is.null(input$chrom) && input$chrom != "All") {
      df <- df %>% filter(seqnames == input$chrom)
    }
    if (!is.null(input$phastcons) && "phastcons" %in% colnames(df)) {
      df <- df %>%
        filter(phastcons >= input$phastcons[1], phastcons <= input$phastcons[2])
    }
    if (!is.null(input$width) && "width" %in% colnames(df)) {
      df <- df %>% filter(width >= input$width[1], width <= input$width[2])
    }

    # Apply each binary-column filter independently (AND-combined).
    # Users can simultaneously require, exclude, or ignore each set.
    for (col in names(BIN_COLS)) {
      if (!col %in% colnames(df)) {
        next
      }
      val <- input[[paste0("flt_", col)]]
      if (is.null(val) || val == "all") {
        next
      }
      if (val == "with") {
        df <- df %>% filter(.data[[col]] == 1)
      } else if (val == "without") {
        df <- df %>% filter(.data[[col]] == 0)
      }
    }

    if (!is.null(input$gene_search) && nzchar(input$gene_search)) {
      pattern <- tolower(input$gene_search)
      df <- df %>%
        filter(
          if_any(
            everything(),
            ~ grepl(pattern, tolower(as.character(.x)), fixed = TRUE)
          )
        )
    }
    df
  })

  # ---- outputs ----
  output$table <- DT::renderDataTable({
    df <- filtered()
    if (is.null(df)) {
      return(DT::datatable(
        data.frame(message = "No data"),
        options = list(dom = 't')
      ))
    }
    DT::datatable(
      df,
      options = list(pageLength = 25, scrollX = TRUE),
      filter = 'top'
    )
  })

  output$stats <- renderPrint({
    df <- filtered()
    if (is.null(df)) {
      cat("No data\n")
      return()
    }
    cat(sprintf("Rows shown: %d (total: %d)\n", nrow(df), nrow(fullData)))
    cat("Columns:", paste(colnames(fullData), collapse = ", "), "\n")
    cat("\nSet-membership counts in current view:\n")
    for (col in names(BIN_COLS)) {
      if (col %in% colnames(df)) {
        cat(sprintf(
          "  %-25s = 1 : %d\n",
          BIN_COLS[[col]],
          sum(df[[col]] == 1, na.rm = TRUE)
        ))
      }
    }
  })

  output$download_filtered <- downloadHandler(
    filename = function() {
      paste0("actinopteriigy_cne_filtered_", Sys.Date(), ".tsv")
    },
    content = function(file) {
      df <- filtered()
      readr::write_tsv(df, file)
    }
  )
}

shinyApp(ui, server)
