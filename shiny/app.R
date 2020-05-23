### Process for saving input data
#rm(list=ls())
#germ <- read.csv('C:/Users/Maggie/Documents/WetlandEcology/RevegModel/ShinyApp/GerminationSeedlingExampleMatrix.csv', header=TRUE)
#early <- read.csv('C:/Users/Maggie/Documents/WetlandEcology/RevegModel/ShinyApp/EarlySeedlingExampleMatrix.csv', header=TRUE)
#late <- read.csv('C:/Users/Maggie/Documents/WetlandEcology/RevegModel/ShinyApp/LateSeedlingExampleMatrix.csv', header=TRUE)
## then save workspace as "global.RData" in same folder as app.R script, or run
#save.image("~/WetlandEcology/RevegModel/ShinyApp/ModelExample/global.RData")

### Load dependencies
library(shiny)
library(DT)
library(ggplot2)
library(RColorBrewer)

### Example option 4
### input environmental variables, output table of options
### based on standard number of seeds sown (1000) per source
ui <- fluidPage(
   # Web Application title
   titlePanel("Wetland Seed Restoration Webapp Example"),
   # Sidebar with a slider input for number of bins 
    sidebarLayout(
      sidebarPanel(
         selectInput(inputId="water",
                  label="Water Availability",
                  choices=c("Wet", "Dry")),
         selectInput("temp",
                  label="Seeding Temperature",
                  c("32-15", "28-10", "36-20"))),
      # Show a plot of the generated distribution
      mainPanel(h1("Top 10 Sources"),
                h2(id="warning", "WARNING: Outputs are not based on real data - yet!"),
                tags$style(HTML("#warning{color: red;}")),
                dataTableOutput("BestSources"),
                h1("Outputs for All Sources"),
                h2(id="warning", "WARNING: Outputs are not based on real data - yet!"),
                plotOutput("SourcePlot"),
                h1("Life Stage Abundance Over Time"),
                h2(id="warning", "WARNING: Outputs are not based on real data - yet!"),
                plotOutput("LifeStage1"),
                plotOutput("LifeStage2"),
                plotOutput("LifeStage3"),
                h2(id="warning", "WARNING: Outputs are not based on real data - yet!"))
   )#end sidebar
)#end fluidPage

server <- function(input, output) {
  ### Load output rates from model (EXAMPLES FOR NOW)
  load('global.RData')

  ### Calculate number of plants in each pop at each timestep
  NSown <- 1000
  germ_ct <- germ
  germ_ct[,5] <- round(NSown*germ[,5], digits=0)
  for (t in 2:(ncol(germ)-4)){
    for (i in 1:nrow(germ_ct)){
      germ_ct[i,t+4] <- round(germ_ct[i,t+3] * germ[i,t+4], digits=0)
    }#for i
  }#for t
  germ_ct[,ncol(germ_ct)] <- germ_ct[,ncol(germ_ct)-1]
  for (i in 1:nrow(germ_ct)){
    for (t in 5:41){
      germ_ct[i,t] <- NSown-germ_ct[i,t]
    }#for i 
  }#for t
  early_ct <- early
  early_ct[,5] <- 0
  for (t in 2:(ncol(germ)-4)){
    for (i in 1:nrow(germ_ct)){
      if ((1-early[i,t+4])==0){
        early_ct[i,t+4] <- early_ct[i,t+3]
      } else {
        early_ct[i,t+4] <- round(germ_ct[i,t+4] * (1-early[i,t+4]), digits=0)
      }#else
    }# for i
  }#for t
  late_ct <- late
  late_ct[,5] <- 0
  for (t in 2:(ncol(germ)-4)){
    for (i in 1:nrow(germ_ct)){
      if ((1-late[i,t+4])==0){
        late_ct[i,t+4] <- late_ct[i,t+3]
      } else {
        late_ct[i,t+4] <- round(early_ct[i,t+4] * (1-late[i,t+4]), digits=0)
      }#else
    }#for i
  }# for t
  
  # set up reactive context
  #inputWP <- reactive({
  #  switch(input$water,
  #         "wet" = 1, 
  #         "dry" = 2)
  #})
  #inputTemp <- reactive({
  #  switch(input$temp,
  #         "28-10" = 1,
  #         "32-15" = 2,
  #         "36-20" = )
  #})
  
  ## return subset table matching conditions
  ###NOTE: These need to be adjusted to be able to accept real input values
  outputData <- late_ct[which(late_ct$WP==1 & late_ct$Temp==2),]
  ## get max 10 results for input conditions
  best_options <- outputData[order(outputData[,ncol(outputData)], decreasing=TRUE),]
  best_options_truncated <- best_options[c(1:4, ncol(late_ct))]
  names(best_options_truncated)[5] <- "Final # Plants"  
  output$BestSources <- DT::renderDataTable({best_options_truncated})
  
  ## plot # adult plants for each source over time for subset
  output$SourcePlot <- renderPlot({
    ggplot(data=outputData)+
      geom_bar(aes(x=Source, y=T41, col=Species, fill=Species), stat="identity")+
      xlab("Source population")+
      ylab("Number final plants")+
      theme_classic()+
      theme(axis.text.x = element_text(angle=90))+
      ggtitle("End-of-Season Plant Abundance per Source",
              subtitle = paste("Environmental Conditions:", input$water, "and", input$temp))
  })#end SourcePlot

  ## plot germination per source
  output$LifeStage1 <- renderPlot({
    # convert data to long format
    ## NOTE: Probably worth saving this in global.RData file for efficiency of webapp
    long_format <- as.data.frame(matrix(NA, nrow=nrow(germ_ct)*ncol(germ_ct), ncol=6))
    row <- 0#set counter
    names(long_format) <- c("Species", "Source", "WP", "Temp", "Timestep", "Count")
    for (i in 1:nrow(germ_ct)){
      for (t in 1:(ncol(germ_ct)-4)){
        row = row + 1
        long_format$Species[row] <- as.character(germ_ct$Species[i])
        long_format$Source[row] <- as.character(germ_ct$Source[i])
        long_format$WP[row] <- germ_ct$WP[i]
        long_format$Temp[row] <- germ_ct$Temp[i]
        long_format$Timestep[row] <- t
        long_format$Count[row] <- germ_ct[i,t+4]
      }#t
    }#i
    # take subset of long format based on inputs
    germ_subset <- long_format[which(long_format$WP==1 & long_format$Temp==2),]     # run plot
    ggplot(germ_subset)+
      geom_line(aes(x=Timestep, y=Count, group=Source, col=Source))+
      theme_bw()+
      facet_wrap(~Species)+
      ggtitle("Germination Counts Over Time", subtitle="per species and source")+
      xlab("Timestep")+
      ylab("Number germinations")
  })#end lifestage1 plot

  ## plot early seedlings per source
  output$LifeStage2 <- renderPlot({
    # convert data to long format
    ## NOTE: Probably worth saving this in global.RData file for efficiency of webapp
    long_format <- as.data.frame(matrix(NA, nrow=nrow(early_ct)*ncol(early_ct), ncol=6))
    row <- 0#set counter
    names(long_format) <- c("Species", "Source", "WP", "Temp", "Timestep", "Count")
    for (i in 1:nrow(early_ct)){
      for (t in 1:(ncol(early_ct)-4)){
        row = row + 1
        long_format$Species[row] <- as.character(early_ct$Species[i])
        long_format$Source[row] <- as.character(early_ct$Source[i])
        long_format$WP[row] <- early_ct$WP[i]
        long_format$Temp[row] <- early_ct$Temp[i]
        long_format$Timestep[row] <- t
        long_format$Count[row] <- early_ct[i,t+4]
      }#t
    }#i
    # take subset of long format based on inputs
    early_subset <- long_format[which(long_format$WP==1 & long_format$Temp==2),]     # run plot
    ggplot(early_subset)+
      geom_line(aes(x=Timestep, y=Count, group=Source, col=Source))+
      theme_bw()+
      facet_wrap(~Species)+
      ggtitle("Early (7-day) Seedlings Counts Over Time", subtitle="per species and source")+
      xlab("Timestep")+
      ylab("Number early seedlings")
  })#end lifestage2 plot
  
  ## plot late seedling per source
  output$LifeStage3 <- renderPlot({
    # convert data to long format
    ## NOTE: Probably worth saving this in global.RData file for efficiency of webapp
    long_format <- as.data.frame(matrix(NA, nrow=nrow(late_ct)*ncol(late_ct), ncol=6))
    row <- 0#set counter
    names(long_format) <- c("Species", "Source", "WP", "Temp", "Timestep", "Count")
    for (i in 1:nrow(late_ct)){
      for (t in 1:(ncol(late_ct)-4)){
        row = row + 1
        long_format$Species[row] <- as.character(late_ct$Species[i])
        long_format$Source[row] <- as.character(late_ct$Source[i])
        long_format$WP[row] <- late_ct$WP[i]
        long_format$Temp[row] <- late_ct$Temp[i]
        long_format$Timestep[row] <- t
        long_format$Count[row] <- late_ct[i,t+4]
      }#t
    }#i
    # take subset of long format based on inputs
    late_subset <- long_format[which(long_format$WP==1 & long_format$Temp==2),]     # run plot
    ggplot(late_subset)+
      geom_line(aes(x=Timestep, y=Count, group=Source, col=Source))+
      theme_bw()+
      facet_wrap(~Species)+
      ggtitle("Late (21-day) Seedling Counts Over Time", subtitle="per species and source")+
      xlab("Timestep")+
      ylab("Number germinations")
  })#end lifestage3 plot
    
}#end server function

# Run the application 
shinyApp(ui=ui, server=server)
