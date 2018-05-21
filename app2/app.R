library(tidyverse)
library(shiny)

#data <- read_csv("C:\\Users\\Jonny\\Google Drive\\Academic\\FGV-SP\\Curso_R_Slides\\Class_9\\data.csv")
#data <- data %>% mutate(NOME_PARTIDO=iconv(NOME_PARTIDO, "ASCII", "UTF-8", sub=""))

ui <- fluidPage( 
  
  selectInput("Estado",
              label="Estado",
              choices = c("SP","RJ"),
              selected = "SP"),
  
  plotOutput("Graph")
)

server <- function(input, output, session) { 
  
  new_data <- reactive({
    data %>% filter(UF==input$Estado)
  })
  
  output$Graph <- renderPlot({
    new_data() %>% ggplot() +
      geom_col(aes(x=NOME_PARTIDO,y=QTDE_VOTOS))
  })
}
shinyApp(ui = ui, server = server)