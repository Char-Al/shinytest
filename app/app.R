library(shiny)

ui <- fluidPage(
  titlePanel("Mon app Shiny sur GitHub Pages"),
  sidebarLayout(
    sidebarPanel(
      sliderInput(
        inputId = "n_bins",
        label = "Nombre de barres :",
        min = 10,
        max = 60,
        value = 30
      ),
      selectInput(
        inputId = "variable",
        label = "Variable :",
        choices = c(
          "Longueur du pétale" = "Petal.Length",
          "Largeur du pétale" = "Petal.Width",
          "Longueur du sépale" = "Sepal.Length",
          "Largeur du sépale" = "Sepal.Width"
        )
      )
    ),
    mainPanel(
      plotOutput(outputId = "dist_plot")
    )
  )
)

server <- function(input, output, session) {
  output$dist_plot <- renderPlot({
    hist(
      iris[[input$variable]],
      breaks = seq(0, 8, length.out = input$n_bins + 1),
      col = "#447099",
      border = "white",
      main = input$variable,
      xlab = input$variable
    )
  })
}

shinyApp(ui, server)
