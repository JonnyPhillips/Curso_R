#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(tidyverse)
library(shiny)

data <- read_csv("C:\\Users\\Jonny\\Google Drive\\Academic\\FGV-SP\\Curso_R_Slides\\Class_9\\data.csv")
data <- data %>% mutate(NOME_PARTIDO=iconv(NOME_PARTIDO, "ASCII", "UTF-8", sub=""))


# Define UI for application that draws a histogram
ui <- fluidPage(
  
  selectInput("Estado",
              label="Estado",
              choices = c("SP","RJ"),
              selected = "SP"),
  
  plotOutput("Graph")
)

# Define server logic required to draw a histogram
server <- function(input, output) {
   
  new_data <- reactive({
    data %>% filter(UF==input$Estado)
  })
  
  output$Graph <- renderPlot({
    new_data() %>% ggplot() +
      geom_col(aes(x=NOME_PARTIDO,y=QTDE_VOTOS)) +
      theme_classic() +
      ylab("Quantidade de Votos")
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)

