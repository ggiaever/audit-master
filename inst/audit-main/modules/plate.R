
library(shiny)
library(ggplot2)
library(dplyr)
library(readgp1)
library(reshape2)

#x = read.csv("/home/ggiaever/RProjects/ggAUDIT/2019-03-13_mtx-dox-hiphop_2ndexp-Results (1).txt",stringsAsFactors = F)



#names(x)[1:3]=c("plate","date","runtime")

#x$plate = as.numeric(substr(x$plate,7,7))

#if(ncol(x)== 51) names(x)[4:51] = paste0(rep(LETTERS[1:6],each=8),rep(1:8,6))






#plateUI <- function(id){
 # ns = NS(id)

plateUI <-  function(id, label = NULL) {
  ns <- NS(id)

  fluidPage(
  titlePanel("plate plot"),

  fluidRow(
    column(3,
      numericInput(ns("plate"),
        h3("plate number"),
        value = 1)),
    column(3,
      numericInput(ns("num"),
        h3("time limit (hrs)"),
        value = 15)),
  plotOutput(ns("plateplot"))

  ))
}

# Define server logic ----
plateSERV <- function(input, output,session,datafile) {
  outputData <- reactiveVal(value = NULL)



  data = reactive({
    #dplyr::distinct(datafile(), run, plate, well)
    subdata = datafile() %>% filter(plate == as.numeric(input$plate))


    # n = 0:(nrow( subdata)-1)
    # n2 = n*0.2189
    #
    subdata$runtime = round(subdata$runtime/3600,2)
    #
    # m = melt(subdata[,3:51],id.vars = "runtime",measure.vars = names(subdata)[4:51],variable.name = "well",value.name = "measure")
    subdata$row = substring(subdata$well,1,1)
    subdata$coln = substring(subdata$well,2,2)
    return(subdata)

  })

  #datafile = x %>% filter(plate == as.numeric(input$plate))
  output$plateplot <- renderPlot({
  g=ggplot(subset(data(), runtime <= as.numeric(input$num)),aes(x = runtime,y=measure, col = well))+ geom_point()+theme_bw()
  #g2 = g + geom_line(aes(x=runtime,y = ref,col = "black"))
  g1 = g+ facet_grid(row~coln) +  theme(legend.position="none")
  g1=g1+theme(strip.background =element_rect(fill="dodgerblue")) + theme(strip.text = element_text(colour = 'white'))

  g1

})
}


# Run the app ----
shinyApp(ui = plateUI , server = plateSERV)



















