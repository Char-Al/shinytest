library(shiny)
library(bslib)
library(plotly)
library(DT)

phenotype_colors <- c(USH = "#e74c3c", RP = "#00b894")
user_color <- "#74b9ff"
protein_length <- 5202
chart_title <- "USH2A contact-map"

# Domain reference: UniProt O75445 (usherin), isoform 1 = 5202 aa -----
genotypes_All <- read.csv("dataSetAll.csv", sep = ";")
domains_all <- read.csv("domains.csv")

band_color_for <- function(name) {
  palette <- c(
    "Fibronectin" = "#5b8db8",
    "Laminin" = "#7d6bbf",
    "Transmembrane" = "#e08a3c",
    "Cytoplasmic" = "#9aa5ad"
  )
  hits <- palette[vapply(names(palette), function(k) grepl(k, name, fixed = TRUE), logical(1))]
  if (length(hits) > 0) unname(hits[1]) else "#c8d6e5"
}

# Helper: mirror each genotype pair across the diagonal so that USH sits above the #
# diagonal and RP below it (a pair (i,j) and (j,i) is the same edge on a contact-map) #
map_region <- function(r1, r2, region = c("USH", "RP"), eps = 0.004) {
  region <- match.arg(region)
  x <- r1 / protein_length
  y <- r2 / protein_length
  swap <- (region == "USH" & y < x) | (region == "RP" & y > x)
  xs <- ifelse(swap, y, x)
  ys <- ifelse(swap, x, y)
  on_diagonal <- xs == ys
  if (region == "USH") {
    ys[on_diagonal] <- ys[on_diagonal] + eps
  } else {
    xs[on_diagonal] <- xs[on_diagonal] + eps
  }
  list(x = xs, y = ys)
}

axis_strips <- list()
for (i in seq_len(nrow(domains_all))) {
  dn <- domains_all$name[i]
  ds <- domains_all$start[i] / protein_length
  de <- domains_all$end[i] / protein_length
  col <- band_color_for(dn)
  axis_strips[[length(axis_strips) + 1]] <- list(
    type = "rect", xref = "x", yref = "paper",
    x0 = ds, x1 = de, y0 = 0, y1 = 1,
    fillcolor = col, opacity = 0.05, line = list(width = 0), layer = "below"
  )
  axis_strips[[length(axis_strips) + 1]] <- list(
    type = "rect", xref = "x", yref = "paper",
    x0 = ds, x1 = de, y0 = 0, y1 = 0.05,
    fillcolor = col, opacity = 0.55, line = list(width = 0), layer = "below"
  )
  axis_strips[[length(axis_strips) + 1]] <- list(
    type = "rect", xref = "paper", yref = "y",
    x0 = 0, x1 = 0.05, y0 = ds, y1 = de,
    fillcolor = col, opacity = 0.55, line = list(width = 0), layer = "below"
  )
}

axis_labels <- Filter(Negate(is.null), lapply(seq_len(nrow(domains_all)), function(i) {
  dn <- domains_all$name[i]
  ds <- domains_all$start[i] / protein_length
  de <- domains_all$end[i] / protein_length
  if (de - ds < 0.045) return(NULL)
  list(
    x = (ds + de) / 2, y = 0.025,
    xref = "x", yref = "paper",
    text = dn, textangle = 90, showarrow = FALSE,
    font = list(size = 8, color = "#f5f8fa"), opacity = 0.8,
    xanchor = "center", yanchor = "middle"
  )
}))

# Define UI for USH2A_App -----
ui <- fluidPage(
  theme = bs_theme(bootswatch = "flatly"),
  titlePanel("USH2A_App"),
  sidebarLayout(
    position = "right",
    sidebarPanel(
      "Plot your own genotype:",
      selectInput("genotype", "Select phenotype, note that the first variant has to be the 5' one", choices = c("USH", "RP"), selected = "USH"),
      br(""),
      "First residue position:",
      numericInput("residue1", "for p.(Cys759Phe), enter 759", value = 1, min = 1, max = protein_length),
      "Second residue position:",
      numericInput("residue2", "for Cys934Trp enter 934", value = 5202, min = 1, max = protein_length),
      br(""),
      actionButton("plotButton", "Plot genotype", class = "btn-primary"),
      br(""),
      br(""),
      "Tip: pre-fill the form with a link like ?pheno=USH&r1=759&r2=934",
      br(""),
      "USH2A_App is a shiny app, available at https://char-al.github.io/shinytest/. The app runs entirely in your web browser (WebR), no installation is needed.",
      br(""),
      code("library(shiny); library(plotly)")
    ),
    mainPanel(
      "This app can be used to create a plot, inspired by the contact-map used with Hi-C. The created map allows direct observation of genotype's distribution. For example, you can use it to see if some regions are regularly associated in a phenotype as they will create a cluster of dots",
      plotlyOutput("scatterPlot", height = "500px"),
      hr(),
      DTOutput("genoTable")
    )
  )
)

# Define server logic required to draw the contact-map plot ----
server <- function(input, output, session) {

  # Pre-fill inputs from URL query string (?pheno=USH&r1=759&r2=934) -----
  qs <- parseQueryString(session$clientData$url_search)
  if (!is.null(qs$pheno) && toupper(qs$pheno) %in% names(phenotype_colors)) {
    updateSelectInput(session, "genotype", selected = toupper(qs$pheno))
  }
  for (par in c("r1", "r2")) {
    val <- suppressWarnings(as.numeric(qs[[par]]))
    if (length(val) == 1L && is.finite(val)) {
      updateNumericInput(session, if (par == "r1") "residue1" else "residue2",
                         value = max(1, min(protein_length, val)))
    }
  }

  # Determine coordinates based on the phenotype in the genotypes_All -----
  n_ush <- sum(genotypes_All$phenotype == "USH")
  n_rp <- sum(genotypes_All$phenotype == "RP")
  ush <- map_region(genotypes_All$residue1_position[genotypes_All$phenotype == "USH"],
                    genotypes_All$residue2_position[genotypes_All$phenotype == "USH"], "USH")
  rp <- map_region(genotypes_All$residue1_position[genotypes_All$phenotype == "RP"],
                   genotypes_All$residue2_position[genotypes_All$phenotype == "RP"], "RP")
  x_values <- ush$x
  y_values <- ush$y
  z_values <- rp$x
  a_values <- rp$y

  # Initialize the plot with the line representing the 5202 residue protein ------
  p <- plot_ly() %>%
    add_trace(type = "scatter", mode = "lines", x = c(0, 1), y = c(0, 1), line = list(color = "black"), name = "USH2A protein") %>%
    add_trace(type = "scatter", mode = "markers", x = x_values, y = y_values, marker = list(color = "#e74c3c", size = 4), hoverinfo = "text", text = paste("variant 1:", genotypes_All$residue1_nucleotide[genotypes_All$phenotype =="USH"], " - ", genotypes_All$residue1_protein[genotypes_All$phenotype == "USH"], ", ", genotypes_All$residue1_domain[genotypes_All$phenotype == "USH"],", (GN005) :", genotypes_All$residue1_GN005[genotypes_All$phenotype == "USH"], "<br>variant 2:", genotypes_All$residue2_nucleotide[genotypes_All$phenotype == "USH"], " - ", genotypes_All$residue2_protein[genotypes_All$phenotype == "USH"],", ", genotypes_All$residue2_domain[genotypes_All$phenotype == "USH"],", (GN005) :", genotypes_All$residue2_GN005[genotypes_All$phenotype == "USH"]) ,name = "USH all genotypes") %>%
    add_trace(type = "scatter", mode = "markers", x = z_values, y = a_values, marker = list(color = "#00b894", size = 3), hoverinfo = "text", text = paste("variant 1:", genotypes_All$residue1_nucleotide[genotypes_All$phenotype == "RP"], " - ", genotypes_All$residue1_protein[genotypes_All$phenotype == "RP"], ", ",genotypes_All$residue1_domain[genotypes_All$phenotype == "RP"],", (GN005) :", genotypes_All$residue1_GN005[genotypes_All$phenotype == "RP"], "<br>variant 2:", genotypes_All$residue2_nucleotide[genotypes_All$phenotype == "RP"], " - ", genotypes_All$residue2_protein[genotypes_All$phenotype == "RP"],", ", genotypes_All$residue2_domain[genotypes_All$phenotype == "RP"],", (GN005) :", genotypes_All$residue2_GN005[genotypes_All$phenotype == "RP"]), name = "RP all genotypes") %>%
    layout(
      title = chart_title,
      shapes = axis_strips,
      annotations = axis_labels,
      xaxis = list(range = c(0, 1), showticklabels = FALSE, showline = FALSE, linewidth = 2, mirror = "ticks", showgrid = FALSE, zeroline = FALSE, rangeslider = list(visible = FALSE)),
      yaxis = list(range = c(0, 1), showticklabels = FALSE, showline = FALSE, linewidth = 2, mirror = "ticks", showgrid = FALSE, zeroline = FALSE, rangeslider = list(visible = FALSE)),
      showlegend = TRUE
    )

  current_plot <- reactiveVal(p)

  # Show the plot with pre-calculted points before eventual update by the user -----
  output$scatterPlot <- renderPlotly({
    config(current_plot(),
      displaylogo = FALSE,
      toImageButtonOptions = list(format = "png", filename = "USH2A_contact_map", scale = 2)
    )
  })

  # Genotype reference table, extended with user-added genotypes -----
  user_rows <- reactiveVal(NULL)
  output$genoTable <- renderDT({
    base_table <- data.frame(
      Source = "Dataset",
      Phenotype = genotypes_All$phenotype,
      Residue1 = genotypes_All$residue1_position,
      Variant1 = genotypes_All$residue1_protein,
      Residue2 = genotypes_All$residue2_position,
      Variant2 = genotypes_All$residue2_protein,
      stringsAsFactors = FALSE
    )
    DT::datatable(
      rbind(base_table, user_rows()),
      rownames = FALSE,
      options = list(pageLength = 25, scrollX = TRUE)
    )
  })

  # Add points to the plot when the button is clicked -----
  observeEvent(input$plotButton, {
    user_point <- map_region(input$residue1, input$residue2, input$genotype)
    hover <- paste0("User genotype (", input$genotype, ")<br>variant 1: residue ", input$residue1, "<br>variant 2: residue ", input$residue2)

    current_plot(current_plot() %>%
      add_trace(type = "scatter", mode = "markers", x = user_point$x, y = user_point$y,
                marker = list(color = user_color, size = 10),
                hoverinfo = "text",
                text = hover,
                showlegend = FALSE))

    user_rows(rbind(user_rows(), data.frame(
      Source = "User",
      Phenotype = input$genotype,
      Residue1 = input$residue1,
      Variant1 = paste0("residue ", input$residue1),
      Residue2 = input$residue2,
      Variant2 = paste0("residue ", input$residue2),
      stringsAsFactors = FALSE
    )))
  })
}

# Run the USH2A_app ----
shinyApp(ui = ui, server = server)
