#' Launch the DEGinsight Shiny interface
#'
#' Prepare enrichment and preview the evidence sent to the configured model.
#' The app only submits a model request when you click Run LLM. It reads the
#' model endpoint and credentials from LLM_BASE_URL and LLM_API_KEY.
#'
#' @param launch.browser Whether to open a browser when the app starts.
#' @param port Local port. NULL lets Shiny choose an available port.
#' @return The value returned by `shiny::runApp()` when the app stops.
#' @examples
#' \dontrun{
#' run_shiny()
#' }
#' @export
run_shiny <- function(launch.browser = interactive(), port = NULL) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Install shiny with install.packages('shiny').", call. = FALSE)
  }
  shiny::runApp(shiny::shinyApp(shiny_ui(), shiny_server),
                host = "127.0.0.1", port = port, launch.browser = launch.browser)
}

#' @noRd
shiny_ui <- function() {
  shiny::fluidPage(
    title = "DEGinsight",
    shiny::titlePanel("DEGinsight"),
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        shiny::selectInput("mode", "Input",
          c("Bundled example" = "example", "Enrichment table" = "enrichment",
            "Differential expression for ORA" = "ora", "Signed ranks for GSEA" = "gsea")),
        shiny::conditionalPanel("input.mode !== 'example'",
          shiny::fileInput("table", "Input table", accept = c(".csv", ".tsv")),
          shiny::textOutput("schema")),
        shiny::conditionalPanel("input.mode === 'ora' || input.mode === 'gsea'",
          shiny::fileInput("gmt", "Gene sets", accept = ".gmt")),
        shiny::conditionalPanel("input.mode === 'ora'",
          shiny::selectInput("direction", "Direction", c("mixed", "up", "down")),
          shiny::numericInput("padj", "Adjusted P cutoff", 0.05, min = 0, max = 1),
          shiny::numericInput("logfc", "Absolute log2 fold change cutoff", 1, min = 0)),
        shiny::textInput("contrast", "Contrast", "treated versus control"),
        shiny::textInput("species", "Species", "Homo sapiens"),
        shiny::textInput("tissue", "Tissue or cell type (optional)"),
        shiny::numericInput("top_n", "Maximum pathways in prompt", 20, min = 1, step = 1),
        shiny::textInput("model", "Model ID", "deepseek-flash"),
        shiny::actionButton("prepare", "Prepare evidence", class = "btn-primary"),
        shiny::actionButton("interpret", "Run LLM"),
        shiny::helpText("Run LLM sends the last prepared prompt. Prepare again after changing inputs."),
        shiny::textOutput("status")
      ),
      shiny::mainPanel(
        shiny::tabsetPanel(
          shiny::tabPanel("Enrichment", shiny::tableOutput("enrichment"),
                          shiny::downloadButton("table_csv", "Download enrichment CSV")),
          shiny::tabPanel("Prompt", shiny::verbatimTextOutput("prompt")),
          shiny::tabPanel("Interpretation", shiny::uiOutput("hypotheses")),
          shiny::tabPanel("Citation audit", shiny::verbatimTextOutput("audit"),
                          shiny::downloadButton("result_json", "Download interpretation JSON"))
        )
      )
    )
  )
}

#' @noRd
shiny_read_table <- function(upload, required) {
  if (is.null(upload)) stop("Upload an input table.", call. = FALSE)
  extension <- tolower(tools::file_ext(upload$name))
  if (!extension %in% c("csv", "tsv")) stop("Use a CSV or TSV input table.", call. = FALSE)
  read <- if (extension == "csv") utils::read.csv else utils::read.delim
  tab <- read(upload$datapath, check.names = FALSE, stringsAsFactors = FALSE)
  if (anyDuplicated(names(tab)) || !all(required %in% names(tab))) {
    stop("Required columns: ", paste(required, collapse = ", "), call. = FALSE)
  }
  if (!nrow(tab)) stop("The input table is empty.", call. = FALSE)
  tab
}

#' @noRd
shiny_prepare_enrichment <- function(input) {
  mode <- match.arg(input$mode, c("example", "enrichment", "ora", "gsea"))
  if (mode == "example") {
    example <- new.env()
    utils::data("example_data", package = "DEGinsight", envir = example)
    return(example$example_ora)
  }
  columns <- switch(mode, enrichment = c("pathway", "p_adjust"),
                    ora = c("gene", "padj", "log2FoldChange"), gsea = c("gene", "stat"))
  tab <- shiny_read_table(input$table, columns)
  if (mode == "enrichment") return(tab)
  if (is.null(input$gmt)) stop("Upload a GMT file.", call. = FALSE)
  sets <- read_gmt(input$gmt$datapath)
  if (!length(sets)) stop("No gene sets could be read from the GMT file.", call. = FALSE)
  if (mode == "gsea") return(run_gsea(stats::setNames(tab$stat, tab$gene), sets))
  if (!is.finite(input$padj) || input$padj < 0 || input$padj > 1 ||
      !is.finite(input$logfc) || input$logfc < 0) {
    stop("Use an adjusted P cutoff in [0, 1] and a nonnegative fold change cutoff.")
  }
  selected <- select_ora_genes(tab, "gene", "padj", "log2FoldChange",
                               input$padj, input$logfc, NULL, input$direction)
  run_ora(selected$genes, selected$universe, sets)
}

#' @noRd
shiny_server <- function(input, output, session) {
  prepared <- shiny::reactiveVal(NULL)
  result <- shiny::reactiveVal(NULL)
  status <- shiny::reactiveVal("Choose an input and prepare evidence.")
  safe <- function(expr) {
    tryCatch(expr, error = function(e) {
      status(conditionMessage(e))
      NULL
    })
  }
  output$status <- shiny::renderText(status())
  output$schema <- shiny::renderText(switch(input$mode,
    enrichment = "Columns: pathway, p_adjust. Optional: genes (comma separated).",
    ora = "Columns: gene, padj, log2FoldChange. Include all tested genes.",
    gsea = "Columns: gene, stat. Use unique genes and signed ranks."))
  shiny::observeEvent(input$prepare, {
    prepared(NULL)
    result(NULL)
    current <- safe({
      enrichment <- shiny_prepare_enrichment(input)
      model <- trimws(input$model)
      context <- list(contrast = trimws(input$contrast), species = trimws(input$species),
                      tissue = if (nzchar(trimws(input$tissue))) trimws(input$tissue) else NULL,
                      direction = if (input$mode == "ora") input$direction else "mixed",
                      top_n = input$top_n)
      if (!nzchar(model) || !nzchar(context$contrast) || !nzchar(context$species)) {
        stop("Enter the model, contrast and species.", call. = FALSE)
      }
      prompt <- do.call(build_interpretation_prompt, c(list(enrichment = enrichment), context))
      list(table = enrichment, context = context, prompt = prompt, model = model)
    })
    if (!is.null(current)) {
      prepared(current)
      n <- length(current$prompt$input_pathways)
      status(if (n) sprintf("%d pathways prepared for %s. Review the prompt before running the model.",
                            n, current$model) else
            "No pathways available for interpretation.")
    }
  })
  shiny::observeEvent(input$interpret, {
    result(NULL)
    current <- prepared()
    if (is.null(current) || !length(current$prompt$input_pathways)) {
      status("Prepare evidence containing at least one pathway first.")
      return()
    }
    res <- safe(shiny::withProgress(message = "Requesting model interpretation", value = 0.5, {
      do.call(interpret_enrichment, c(list(enrichment = current$table, model = current$model),
                                      current$context))
    }))
    if (!is.null(res)) {
      result(res)
      status(paste(length(res$hypotheses), "hypotheses returned."))
    }
  })
  output$enrichment <- shiny::renderTable({
    shiny::req(prepared())
    utils::head(prepared()$table, 30)
  })
  output$prompt <- shiny::renderText({
    shiny::req(prepared())
    prepared()$prompt$prompt
  })
  output$hypotheses <- shiny::renderUI({
    x <- result()
    shiny::req(x)
    shiny::tagList(
      shiny::tags$p(x$overall_summary),
      lapply(x$hypotheses, function(h) shiny::tagList(
        shiny::tags$h4(h$name),
        shiny::tags$p(h$mechanism),
        shiny::tags$p("Supporting pathways:", paste(unlist(h$supporting_pathways), collapse = ", ")),
        shiny::tags$p("Driver genes:", paste(unlist(h$driver_genes), collapse = ", ")),
        shiny::tags$p("Model evidence label:", h$evidence_grade),
        shiny::tags$p("Proposed test:", h$follow_up_experiment)
      ))
    )
  })
  output$audit <- shiny::renderPrint({
    x <- result()
    shiny::req(x)
    list(grounded_pathways = x$validation$grounded,
         unsupported_pathways = x$validation$unsupported_pathways,
         unsupported_driver_genes = x$validation$unsupported_driver_genes)
  })
  output$table_csv <- shiny::downloadHandler(
    filename = "deginsight-enrichment.csv",
    content = function(path) utils::write.csv(shiny::req(prepared())$table, path, row.names = FALSE)
  )
  output$result_json <- shiny::downloadHandler(
    filename = "deginsight-interpretation.json",
    content = function(path) jsonlite::write_json(unclass(shiny::req(result())), path,
                                                  auto_unbox = TRUE, pretty = TRUE)
  )
}
