library(shiny)
library(plotly)

# Load cured genotypes from CSV file (appdir is the working directory) -----
genotypes_All <- read.csv("dataSetAll.csv", sep = ";")

# Define UI for USH2A_App -----
ui <- fluidPage(
  titlePanel("USH2A_App"),
  sidebarLayout(
    position = "right",
    sidebarPanel(
      "Plot your own genotype:",
      selectInput("genotype", "Select phenotype, note that the first variant has to be the 5' one", choices = c("USH", "RP"), selected = "USH"),
      br(""),
      "First residue position:",
      numericInput("residue1", "for p.(Cys759Phe), enter 759", value = 1),
      "Second residue position:",
      numericInput("residue2", "for Cys934Trp enter 934", value = 5202),
      br(""),
      actionButton("plotButton", "Plot genotype"),
      br(""),
      "USH2A_App is a shiny app, available at https://char-al.github.io/shinytest/. The app runs entirely in your web browser (WebR), no installation is needed.",
      br(""),
      code("library(shiny); library(plotly)")
    ),
    mainPanel("This app can be used to create a plot, inspired by the contact-map used with Hi-C. The created map allows direct observation of genotype's distribution. For example, you can use it to see if some regions are regularly associated in a phenotype as they will create a cluster of dots", plotlyOutput("scatterPlot"))
  )
)

# Define server logic required to draw the contact-map plot ----
server <- function(input, output) {

  # Determine coordinates based on the phenotype in the genotypes_All -----
  x_values <- genotypes_All$residue1_position[genotypes_All$phenotype == "USH"] / 5202
  y_values <- genotypes_All$residue2_position[genotypes_All$phenotype == "USH"] / 5202
  z_values <- genotypes_All$residue1_position[genotypes_All$phenotype == "RP"] / 5202
  a_values <- genotypes_All$residue2_position[genotypes_All$phenotype == "RP"] / 5202

  # Initialize the plot with the line representing the 5202 residue protein ------
  p <- plot_ly() %>%
    add_trace(type = "scatter", mode = "lines", x = c(0, 1), y = c(0, 1), line = list(color = "black"), name = "USH2A protein") %>%
    add_trace(type = "scatter", mode = "markers", x = x_values, y = y_values, marker = list(color = "#e74c3c", size = 4), hoverinfo = "text", text = paste("variant 1:", genotypes_All$residue1_nucleotide[genotypes_All$phenotype =="USH"], " - ", genotypes_All$residue1_protein[genotypes_All$phenotype == "USH"], ", ", genotypes_All$residue1_domain[genotypes_All$phenotype == "USH"],", (GN005) :", genotypes_All$residue1_GN005[genotypes_All$phenotype == "USH"], "<br>variant 2:", genotypes_All$residue2_nucleotide[genotypes_All$phenotype == "USH"], " - ", genotypes_All$residue2_protein[genotypes_All$phenotype == "USH"],", ", genotypes_All$residue2_domain[genotypes_All$phenotype == "USH"],", (GN005) :", genotypes_All$residue2_GN005[genotypes_All$phenotype == "USH"]) ,name = "USH all genotypes") %>%
    add_trace(type = "scatter", mode = "markers", x = z_values, y = a_values, marker = list(color = "#00b894", size = 3), hoverinfo = "text", text = paste("variant 1:", genotypes_All$residue1_nucleotide[genotypes_All$phenotype == "RP"], " - ", genotypes_All$residue1_protein[genotypes_All$phenotype == "RP"], ", ",genotypes_All$residue1_domain[genotypes_All$phenotype == "RP"],", (GN005) :", genotypes_All$residue1_GN005[genotypes_All$phenotype == "RP"], "<br>variant 2:", genotypes_All$residue2_nucleotide[genotypes_All$phenotype == "RP"], " - ", genotypes_All$residue2_protein[genotypes_All$phenotype == "RP"],", ", genotypes_All$residue2_domain[genotypes_All$phenotype == "RP"],", (GN005) :", genotypes_All$residue2_GN005[genotypes_All$phenotype == "RP"]), name = "RP all genotypes") %>%
    layout(
      xaxis = list(range = c(0, 1), showticklabels = FALSE, showline = FALSE, linewidth = 2, mirror = "ticks", showgrid = FALSE, zeroline = FALSE, rangeslider = list(visible = FALSE)),
      yaxis = list(range = c(0, 1), showticklabels = FALSE, showline = FALSE, linewidth = 2, mirror = "ticks", showgrid = FALSE, zeroline = FALSE, rangeslider = list(visible = FALSE)),
      showlegend = TRUE
    )

  # Show the plot with pre-calculted points before eventual update by the user -----
  output$scatterPlot <- renderPlotly({
    p
  })

  # Add points to the plot when the button is clicked -----
  observeEvent(input$plotButton, {
    if (input$genotype == "USH") {
      f_value <- input$residue2 / 5202
      g_value <- input$residue1 / 5202
    } else if (input$genotype == "RP") {
      f_value <- input$residue1 / 5202
      g_value <- input$residue2 / 5202
    }

    p <- p %>%
      add_trace(type = "scatter", mode = "markers", x = f_value, y = g_value,
                marker = list(color = "#74b9ff", size = 10),
                hoverinfo = "text",
                text = paste("User variant1:", f_value * 5202, "<br>User variant2:", g_value * 5202))

    # Update the plot -----
    output$scatterPlot <- renderPlotly({
      p
    })
  })
}

# Run the USH2A_app ----
shinyApp(ui = ui, server = server)
