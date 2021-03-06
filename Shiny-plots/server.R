library(shiny) # Load libraries we will be using
library(dplyr)

shinyServer(function(input, output, session) {

  
  observe({
    if (input$selectall > 0) {
      if (input$selectall %% 2 == 0){
        updateCheckboxGroupInput(session  = session, 
                                 inputId  = "country",
                                 choices  = unique(as.character(anim.data$Country)),
                                 selected = unique(as.character(anim.data$Country)))
        
                                # choices  = subset(unique(as.character(anim.data$Country)), unique(as.character(anim.data$Country)) != "None"),
                                # selected = subset(unique(as.character(anim.data$Country)), unique(as.character(anim.data$Country)) != "None"))
        
      } else {
        updateCheckboxGroupInput(session  = session, 
                                 inputId  = "country",
                                 choices  = unique(as.character(anim.data$Country)),
                                 selected = "")
        
      }
  }
  })
  
  yearData <- reactive({
    
#    if(input$radio == "All"){
#      df <- anim.data %>%
#      filter(Year == input$year  ) %>%
#       select(Country, Coefficient.of.Variation, Incidence, Mean.vaccination, Mean.birth.rate) 
#    } else if (input$radio == "Manual Select"){
#      df <- anim.data %>%
#      filter(Year == input$year & Country %in%  input$country) %>%
#      select(Country, Coefficient.of.Variation, Incidence, Mean.vaccination, Mean.birth.rate) 
#    }
#    else{
#      df <- anim.data %>%
#      filter(Year == input$year & Country %in%  input$country) %>%
#      select(Country, Coefficient.of.Variation, Incidence, Mean.vaccination, Mean.birth.rate) 
#    }
    if(length(input$country) > 1){
      
    }
    
    #df <- subset(anim.data, unique(as.character(anim.data$Country)) != "None") %>%
    p = anim.data
    p$Incidence  =  round(sqrt(p$Incidence),2)
      df<- p %>%
      filter(Year == input$year & Country %in%  input$country) %>%
      
      select(Country, Coefficient.of.Variation, Incidence, Mean.vaccination, Mean.birth.rate) 
    
  })
  
#  yearData <- reactive({
    # Filter to the desired year, and put the columns
    # in the order that Google's Bubble Chart expects
    # them (name, x, y, color, size). Also sort by region
    # so that Google Charts orders and colors the regions
    # consistently.
 #   df <- anim.data %>%
        
      #df <- anim.data[anim.data$Country %in% input$country, ]
      #if(input$radio == "All"){
      #  filter(Year == input$year & Country ) %>%
      #  select(Country, Coefficient.of.Variation, Incidence, Mean.vaccination, Mean.birth.rate) 
      #} else {
#      print(input$year)
 #       filter(Year == input$year & Country  %in%  input$country) %>%
#        select(Country, Coefficient.of.Variation, Incidence, Mean.vaccination, Mean.birth.rate) 
      #}
    
  
      
    #arrange(Region)
#  })
  
  output$chart <- reactive({
    # Return the data and options
    list(
      data = googleDataTable(yearData()),
      options = list(
        title = sprintf(
          "Incidence vs coefficient of variation %s-%s ",
          floor(input$year), floor(input$year) + window.length - 1)#,
        #series = series
      )
    )
  })
})