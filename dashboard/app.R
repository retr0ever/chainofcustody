library(shiny)
library(reticulate)
library(DT)

app_dir <- normalizePath(getwd(), mustWork = TRUE)
repo_root <- normalizePath(file.path(app_dir, ".."), mustWork = TRUE)

py <- file.path(repo_root, ".venv", "bin", "python")
if (!file.exists(py)) py <- file.path(repo_root, ".venv", "Scripts", "python.exe")

if (!file.exists(py)) {
  stop("Could not find uv python at: ", py,
       "\nRun from repo root:\n  uv venv .venv\n  uv sync\n")
}

use_python(py, required = TRUE)
py_run_string(sprintf('import sys; sys.path.insert(0, "%s")', repo_root))

api <- import("chainofcustody.dashboard_api.api", delay_load = FALSE)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

ui <- fluidPage(
  titlePanel("ChainOfCustody — UTR Optimizer"),

  sidebarLayout(
    sidebarPanel(
      textInput("gene", "Gene", value = "POU5F1"),
      textInput("cell", "Target cell type", value = "Dendritic_cell"),

      numericInput("utr5_min", "5'UTR min", value = 10, min = 1),
      numericInput("utr5_max", "5'UTR max", value = 100, min = 1),
      numericInput("pop", "Population", value = 128, min = 2),
      numericInput("ngen", "Generations", value = 10, min = 1),
      textInput("run_name", "Run name (optional)", value = ""),

      actionButton("run", "Run", class = "btn-primary"),
      tags$hr(),
      tags$small("Python: "),
      verbatimTextOutput("py_path", placeholder = TRUE)
    ),

    mainPanel(
      uiOutput("status"),
      tags$hr(),

      tags$h4("Best summary"),
      verbatimTextOutput("summary"),

      tags$hr(),
      tags$h4("miRNA binding sites used"),
      DTOutput("mirna_tbl"),

      tags$hr(),
      tags$h4("Secondary structure plots (best)"),
      fluidRow(
        column(4, imageOutput("plot_full")),
        column(4, imageOutput("plot_utr5")),
        column(4, imageOutput("plot_utr3"))
      ),

      tags$hr(),
      tags$h4("History"),
      DTOutput("history_tbl"),

      tags$hr(),
      tags$h4("Python traceback (only if error)"),
      verbatimTextOutput("traceback")
    )
  )
)

server <- function(input, output, session) {

  output$py_path <- renderPrint({
    cfg <- py_config()
    cat(cfg$python, "\n")
  })

  res <- eventReactive(input$run, {
    withProgress(message = "Running Python pipeline", value = 0, {
      incProgress(0.05)

      run_name <- if (nchar(trimws(input$run_name)) == 0) NULL else input$run_name

      out <- api$optimize_and_plot(
        gene = input$gene,
        target_cell_type = input$cell,
        utr5_min = as.integer(input$utr5_min),
        utr5_max = as.integer(input$utr5_max),
        pop_size = as.integer(input$pop),
        n_gen = as.integer(input$ngen),
        out_dir = file.path(repo_root, "runs"),
        run_name = run_name
      )

      incProgress(0.95)
      out
    })
  })

  output$status <- renderUI({
    req(res())
    if (isTRUE(res()$ok)) {
      div(style="padding:10px; background:#e8f5e9; border:1px solid #c8e6c9;",
          strong("OK"), " — results written to: ", code(res()$run_dir))
    } else {
      div(style="padding:10px; background:#ffebee; border:1px solid #ffcdd2;",
          strong("ERROR"), " — ", code(res()$error))
    }
  })

  output$summary <- renderPrint({
    req(res())
    if (!isTRUE(res()$ok)) {
      cat("No summary. Fix the error above.\n")
      return()
    }
    best <- res()$best
    lens <- res()$lengths
    cat("Run dir:", res()$run_dir, "\n")
    cat("CDS length:", lens$cds, "nt\n")
    cat("3'UTR length:", lens$utr3, "nt\n\n")
    cat("Best label:", best$label, "\n")
    cat("Overall:", best$overall, "\n")
  })

  output$mirna_tbl <- renderDT({
    req(res())
    dat <- tryCatch(as.data.frame(res()$mirna_sites), error = function(e) data.frame())
    datatable(dat, options = list(pageLength = 10, scrollX = TRUE))
  })

  output$history_tbl <- renderDT({
    req(res())
    dat <- tryCatch(as.data.frame(res()$history), error = function(e) data.frame())
    datatable(dat, options = list(pageLength = 10, scrollX = TRUE))
  })

  # This is the second half of "display the plot":
  # Shiny serves the PNG file path returned by Python.
  renderPlotImg <- function(key) {
    renderImage({
      req(res())
      path <- ""
      try({ path <- res()$plots[[key]]$files$png }, silent = TRUE)

      if (is.null(path) || nchar(path) == 0 || !file.exists(path)) {
        tmp <- tempfile(fileext = ".gif")
        writeBin(as.raw(c(
          0x47,0x49,0x46,0x38,0x39,0x61,0x01,0x00,0x01,0x00,0x80,0x00,0x00,
          0x00,0x00,0x00,0xff,0xff,0xff,0x21,0xf9,0x04,0x01,0x00,0x00,0x00,
          0x00,0x2c,0x00,0x00,0x00,0x00,0x01,0x00,0x01,0x00,0x00,0x02,0x02,
          0x44,0x01,0x00,0x3b
        )), tmp)
        list(src = tmp, contentType = "image/gif", width = "100%")
      } else {
        list(src = path, contentType = "image/png", width = "100%")
      }
    }, deleteFile = FALSE)
  }

  output$plot_full <- renderPlotImg("full")
  output$plot_utr5 <- renderPlotImg("utr5")
  output$plot_utr3 <- renderPlotImg("utr3")

  output$traceback <- renderPrint({
    req(res())
    if (isTRUE(res()$ok)) "" else cat(res()$traceback %||% "")
  })
}

shinyApp(ui, server)