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

domain_options <- c("All", sort(unique(stats::na.omit(c(genotypes_All$residue1_domain, genotypes_All$residue2_domain)))))
impact_options <- c("All", sort(unique(stats::na.omit(c(genotypes_All$residue1_GN005, genotypes_All$residue2_GN005)))))
class_options <- c("All", sort(unique(stats::na.omit(c(genotypes_All$residue1_class, genotypes_All$residue2_class)))))

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

axis_shapes <- list()
for (i in seq_len(nrow(domains_all))) {
  dn <- domains_all$name[i]
  ds <- domains_all$start[i] / protein_length
  de <- domains_all$end[i] / protein_length
  col <- band_color_for(dn)
  axis_shapes[[length(axis_shapes) + 1]] <- list(
    type = "rect", xref = "x", yref = "paper",
    x0 = ds, x1 = de, y0 = 0, y1 = 1,
    fillcolor = col, opacity = 0.15, line = list(width = 0), layer = "below"
  )
  axis_shapes[[length(axis_shapes) + 1]] <- list(
    type = "rect", xref = "paper", yref = "y",
    x0 = 0, x1 = 1, y0 = ds, y1 = de,
    fillcolor = col, opacity = 0.15, line = list(width = 0), layer = "below"
  )
  axis_shapes[[length(axis_shapes) + 1]] <- list(
    type = "rect", xref = "x", yref = "paper",
    x0 = ds, x1 = de, y0 = 0, y1 = 0.05,
    fillcolor = col, opacity = 0.5, line = list(width = 0), layer = "below"
  )
  axis_shapes[[length(axis_shapes) + 1]] <- list(
    type = "rect", xref = "paper", yref = "y",
    x0 = 0, x1 = 0.05, y0 = ds, y1 = de,
    fillcolor = col, opacity = 0.5, line = list(width = 0), layer = "below"
  )
}

axis_labels <- Filter(Negate(is.null), lapply(seq_len(nrow(domains_all)), function(i) {
  dn <- domains_all$name[i]
  ds <- domains_all$start[i] / protein_length
  de <- domains_all$end[i] / protein_length
  if (de - ds < 0.045) return(NULL)
  list(
    list(
      x = (ds + de) / 2, y = 0.025,
      xref = "x", yref = "paper",
      text = dn, textangle = 90, showarrow = FALSE,
      font = list(size = 8, color = "#f5f8fa"), opacity = 0.8,
      xanchor = "center", yanchor = "middle"
    ),
    list(
      x = 0.025, y = (ds + de) / 2,
      xref = "paper", yref = "y",
      text = dn, textangle = -90, showarrow = FALSE,
      font = list(size = 8, color = "#f5f8fa"), opacity = 0.8,
      xanchor = "center", yanchor = "middle"
    )
  )
}))
axis_labels <- unlist(axis_labels, recursive = FALSE)

# Build the full contact-map plot from a (filtered) dataset frame -----
make_plot <- function(df) {
  df_ush <- df[df$phenotype == "USH" & !is.na(df$residue1_position) & !is.na(df$residue2_position), , drop = FALSE]
  df_rp <- df[df$phenotype == "RP" & !is.na(df$residue1_position) & !is.na(df$residue2_position), , drop = FALSE]
  ush <- map_region(df_ush$residue1_position, df_ush$residue2_position, "USH")
  rp <- map_region(df_rp$residue1_position, df_rp$residue2_position, "RP")

  plot_ly() %>%
    add_trace(type = "scatter", mode = "lines", x = c(0, 1), y = c(0, 1), line = list(color = "black"), name = "USH2A protein") %>%
    add_trace(type = "scatter", mode = "markers", x = ush$x, y = ush$y, marker = list(color = "#e74c3c", size = 4), hoverinfo = "text", text = paste("variant 1:", df_ush$residue1_nucleotide, " - ", df_ush$residue1_protein, ", ", df_ush$residue1_domain, ", (GN005) : ", df_ush$residue1_GN005, "<br>variant 2:", df_ush$residue2_nucleotide, " - ", df_ush$residue2_protein, ", ", df_ush$residue2_domain, ", (GN005) : ", df_ush$residue2_GN005), name = paste0("USH all genotypes (n = ", nrow(df_ush), ")")) %>%
    add_trace(type = "scatter", mode = "markers", x = rp$x, y = rp$y, marker = list(color = "#00b894", size = 3), hoverinfo = "text", text = paste("variant 1:", df_rp$residue1_nucleotide, " - ", df_rp$residue1_protein, ", ", df_rp$residue1_domain, ", (GN005) : ", df_rp$residue1_GN005, "<br>variant 2:", df_rp$residue2_nucleotide, " - ", df_rp$residue2_protein, ", ", df_rp$residue2_domain, ", (GN005) : ", df_rp$residue2_GN005), name = paste0("RP all genotypes (n = ", nrow(df_rp), ")")) %>%
    layout(
      title = chart_title,
      shapes = axis_shapes,
      annotations = axis_labels,
      xaxis = list(autorange = TRUE, showticklabels = FALSE, showline = FALSE, linewidth = 2, mirror = "ticks", showgrid = FALSE, zeroline = FALSE, rangeslider = list(visible = FALSE)),
      yaxis = list(autorange = TRUE, showticklabels = FALSE, showline = FALSE, linewidth = 2, mirror = "ticks", showgrid = FALSE, zeroline = FALSE, rangeslider = list(visible = FALSE)),
      showlegend = TRUE
    )
}

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
      fluidRow(
        column(6, selectInput("cls1", "Class variant 1:", class_options, selected = "All")),
        column(6, selectInput("cls2", "Class variant 2:", class_options, selected = "All"))
      ),
      fluidRow(
        column(6, selectInput("dom1", "Domain variant 1:", domain_options, selected = "All")),
        column(6, selectInput("dom2", "Domain variant 2:", domain_options, selected = "All"))
      ),
      fluidRow(
        column(6, selectInput("imp1", "Impact variant 1 (GN005):", impact_options, selected = "All")),
        column(6, selectInput("imp2", "Impact variant 2 (GN005):", impact_options, selected = "All"))
      ),
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
  # 1) native shiny (RStudio / runApp): session$clientData works -----
  observeEvent(session$clientData$url_search, {
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
  }, once = TRUE)

  # 2) shinylive export: the page query string is not forwarded to the session #
  # (posit-dev/shinylive#79), so a head script pushes it via Shiny.setInputValue #
  observeEvent(input$url_query, {
    q <- tryCatch(jsonlite::fromJSON(input$url_query), error = function(e) NULL)
    if (is.null(q) || !is.list(q)) return()
    pheno <- q$pheno
    if (!is.null(pheno) && toupper(pheno) %in% names(phenotype_colors)) {
      updateSelectInput(session, "genotype", selected = toupper(pheno))
    }
    for (par in c("r1", "r2")) {
      val <- suppressWarnings(as.numeric(q[[par]]))
      if (length(val) == 1L && is.finite(val)) {
        updateNumericInput(session, if (par == "r1") "residue1" else "residue2",
                           value = max(1, min(protein_length, val)))
      }
    }
  })

  # Cross-filtered dataset (All = no filter; NA rows are kept out of active filters) -----
  df_filtered <- reactive({
    c1 <- input$cls1 == "All" | genotypes_All$residue1_class %in% input$cls1
    c2 <- input$cls2 == "All" | genotypes_All$residue2_class %in% input$cls2
    d1 <- input$dom1 == "All" | genotypes_All$residue1_domain %in% input$dom1
    d2 <- input$dom2 == "All" | genotypes_All$residue2_domain %in% input$dom2
    i1 <- input$imp1 == "All" | genotypes_All$residue1_GN005 %in% input$imp1
    i2 <- input$imp2 == "All" | genotypes_All$residue2_GN005 %in% input$imp2
    genotypes_All[c1 & c2 & d1 & d2 & i1 & i2, , drop = FALSE]
  })

  # User-added genotypes, re-applied on top of every filtered rebuild ----
  user_pts <- reactiveVal(data.frame(x = numeric(0), y = numeric(0), text = character(0)))

  current_plot <- reactive({
    pl <- make_plot(df_filtered())
    up <- user_pts()
    if (nrow(up) > 0) {
      pl <- pl %>%
        add_trace(type = "scatter", mode = "markers", x = up$x, y = up$y,
                  marker = list(color = user_color, size = 10),
                  hoverinfo = "text",
                  text = up$text,
                  showlegend = FALSE)
    }
    pl
  })

  # Show the plot with pre-calculted points before eventual update by the user -----
  output$scatterPlot <- renderPlotly({
    config(current_plot(),
      displaylogo = FALSE,
      toImageButtonOptions = list(format = "png", filename = "USH2A_contact_map", scale = 2)
    )
  })

  # Genotype reference table, synced with the filters and user additions -----
  output$genoTable <- renderDT({
    df <- df_filtered()
    base_table <- data.frame(
      Source = "Dataset",
      Phenotype = df$phenotype,
      Residue1 = df$residue1_position,
      Variant1 = df$residue1_protein,
      Domain1 = df$residue1_domain,
      Impact1 = df$residue1_GN005,
      Residue2 = df$residue2_position,
      Variant2 = df$residue2_protein,
      Domain2 = df$residue2_domain,
      Impact2 = df$residue2_GN005,
      stringsAsFactors = FALSE
    )
    user_rows <- data.frame(
      Source = "User",
      Phenotype = rep(NA_character_, nrow(user_pts())),
      Residue1 = NA_real_,
      Variant1 = NA_character_,
      Domain1 = NA_character_,
      Impact1 = NA_character_,
      Residue2 = NA_real_,
      Variant2 = NA_character_,
      Domain2 = NA_character_,
      Impact2 = NA_character_,
      stringsAsFactors = FALSE
    )
    if (nrow(user_pts()) > 0) {
      user_rows$Phenotype <- sub("User genotype \\((\\w+)\\).*", "\\1", user_pts()$text)
      user_rows$Residue1 <- as.numeric(sub(".*variant 1: residue (\\d+).*", "\\1", user_pts()$text))
      user_rows$Residue2 <- as.numeric(sub(".*variant 2: residue (\\d+).*", "\\1", user_pts()$text))
      user_rows$Variant1 <- paste0("residue ", user_rows$Residue1)
      user_rows$Variant2 <- paste0("residue ", user_rows$Residue2)
    }
    DT::datatable(
      rbind(base_table, user_rows),
      rownames = FALSE,
      options = list(pageLength = 25, scrollX = TRUE)
    )
  })

  # Add points to the plot when the button is clicked -----
  observeEvent(input$plotButton, {
    user_point <- map_region(input$residue1, input$residue2, input$genotype)
    hover <- paste0("User genotype (", input$genotype, ")<br>variant 1: residue ", input$residue1, "<br>variant 2: residue ", input$residue2)

    user_pts(rbind(user_pts(), data.frame(
      x = user_point$x,
      y = user_point$y,
      text = hover,
      stringsAsFactors = FALSE
    )))
  })
}

# Run the USH2A_app ----
shinyApp(ui = ui, server = server)
