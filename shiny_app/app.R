library(shiny)
library(DT)
library(readr)
library(dplyr)

# Ensure working directory is the project root (so relative paths like "output/..." resolve)
find_project_root <- function(target = file.path("output", "teleost_specific_cne.csv"), max_up = 6) {
  p <- tryCatch(normalizePath(getwd(), winslash = "/"), error = function(e) getwd())
  for (i in 0:max_up) {
    candidate <- file.path(p, target)
    if (file.exists(candidate)) return(p)
    parent <- dirname(p)
    if (identical(parent, p)) break
    p <- parent
  }
  return(NULL)
}

proj_root <- find_project_root()
if (!is.null(proj_root)) {
  if (!identical(normalizePath(getwd(), winslash = "/"), normalizePath(proj_root, winslash = "/"))) {
    message("Setting working directory to project root: ", proj_root)
    setwd(proj_root)
  }
} else {
  message("Warning: could not locate 'output/teleost_specific_cne.csv' in parent paths; leaving working directory as ", getwd())
}

# Path to the CNE CSV (relative to project root)
csv_path <- file.path("output", "teleost_specific_cne.csv")

fullData <- tryCatch(
  readr::read_csv(csv_path, show_col_types = FALSE),
  error = function(e) {
    message(sprintf("Failed to read %s: %s", csv_path, e$message))
    NULL
  }
)

# Try to read ATAC overlap annotations and add integer column `atac_overlap` (1/0)
atac_path <- file.path("output", "teleost_cne_atac_overlaps_annotated.tsv")
atac_overlaps <- tryCatch(
  readr::read_tsv(atac_path, show_col_types = FALSE),
  error = function(e) NULL
)
if (!is.null(fullData)) {
  if (!is.null(atac_overlaps)) {
    atac_marker <- atac_overlaps %>% distinct(seqnames, start, end) %>% mutate(atac_overlap = 1L)
    fullData <- fullData %>% left_join(atac_marker, by = c("seqnames", "start", "end")) %>%
      mutate(atac_overlap = ifelse(is.na(atac_overlap), 0L, as.integer(atac_overlap)))
  } else {
    fullData <- fullData %>% mutate(atac_overlap = 0L)
  }
}

ui <- fluidPage(
  titlePanel("Teleost-specific CNE browser"),
  sidebarLayout(
    sidebarPanel(
      uiOutput("chrom_ui"),
      uiOutput("phastcons_ui"),
      uiOutput("width_ui"),
      uiOutput("atac_filter_ui"),
      textInput("gene_search", "Gene / annotation contains:", value = ""),
      actionButton("reset", "Reset filters"),
      br(),
      downloadButton("download_filtered", "Download filtered CSV")
    ),
    mainPanel(
      DT::dataTableOutput("table"),
      br(),
      verbatimTextOutput("stats")
    )
  )
)

server <- function(input, output, session) {
  data <- reactiveVal(fullData)

  observeEvent(input$reset, {
    if (is.null(fullData)) return()
    try(updateSelectInput(session, "chrom", selected = "All"), silent = TRUE)
    # reset sliders to full-data ranges
    ph_min <- min(fullData$phastcons, na.rm = TRUE)
    ph_max <- max(fullData$phastcons, na.rm = TRUE)
    updateSliderInput(session, "phastcons", value = c(ph_min, ph_max))
    wd_min <- min(fullData$width, na.rm = TRUE)
    wd_max <- max(fullData$width, na.rm = TRUE)
    updateSliderInput(session, "width", value = c(wd_min, wd_max))
    try(updateRadioButtons(session, "atac_filter", selected = "all"), silent = TRUE)
    updateTextInput(session, "gene_search", value = "")
  })

  output$chrom_ui <- renderUI({
    if (is.null(fullData)) return(tags$div("Data not loaded."))
    chrs <- unique(as.character(fullData$seqnames))
    selectInput("chrom", "Chromosome/contig:", choices = c("All", sort(chrs)), selected = "All")
  })

  output$phastcons_ui <- renderUI({
    if (is.null(fullData)) return(NULL)
    minv <- min(fullData$phastcons, na.rm = TRUE)
    maxv <- max(fullData$phastcons, na.rm = TRUE)
    sliderInput("phastcons", "phastCons score range:", min = minv, max = maxv, value = c(minv, maxv))
  })

  output$width_ui <- renderUI({
    if (is.null(fullData)) return(NULL)
    minw <- min(fullData$width, na.rm = TRUE)
    maxw <- max(fullData$width, na.rm = TRUE)
    sliderInput("width", "CNE width range (bp):", min = minw, max = maxw, value = c(minw, maxw))
  })

  output$atac_filter_ui <- renderUI({
    if (is.null(fullData)) return(NULL)
    if (!"atac_overlap" %in% colnames(fullData)) return(NULL)
    radioButtons("atac_filter", "ATACseq overlap:", choices = c("All" = "all", "With ATAC peak" = "with", "Without ATAC peak" = "without"), selected = "all")
  })

  filtered <- reactive({
    df <- fullData
    if (is.null(df)) return(NULL)

    if (!is.null(input$chrom) && input$chrom != "All") {
      df <- df %>% filter(seqnames == input$chrom)
    }

    if (!is.null(input$phastcons)) {
      df <- df %>% filter(phastcons >= input$phastcons[1], phastcons <= input$phastcons[2])
    }

    if (!is.null(input$width)) {
      df <- df %>% filter(width >= input$width[1], width <= input$width[2])
    }

    if (!is.null(input$atac_filter) && input$atac_filter != "all") {
      if (input$atac_filter == "with") {
        df <- df %>% filter(atac_overlap == 1)
      } else if (input$atac_filter == "without") {
        df <- df %>% filter(atac_overlap == 0)
      }
    }

    if (nzchar(input$gene_search)) {
      pattern <- tolower(input$gene_search)
      df <- df %>% filter(if_any(everything(), ~ grepl(pattern, tolower(as.character(.x)), fixed = TRUE)))
    }

    df
  })

  output$table <- DT::renderDataTable({
    df <- filtered()
    if (is.null(df)) return(DT::datatable(data.frame(message = "No data"), options = list(dom = 't')))
    DT::datatable(df, options = list(pageLength = 25, scrollX = TRUE), filter = 'top')
  })

  output$stats <- renderPrint({
    df <- filtered()
    if (is.null(df)) return("No data")
    cat(sprintf("Rows shown: %d (total: %d)\n", nrow(df), nrow(fullData)))
    cat("Columns:", paste(colnames(fullData), collapse = ", "), "\n")
  })

  output$download_filtered <- downloadHandler(
    filename = function() {
      paste0("teleost_specific_cne_filtered_", Sys.Date(), ".csv")
    },
    content = function(file) {
      df <- filtered()
      readr::write_csv(df, file)
    }
  )
}

shinyApp(ui, server)
