library(shiny) # Load libraries we will be using
library(dplyr)

shinyServer(function(input, output, session) {

  
  observe({
    if (input$selectall > 0) {
      if (input$selectall %% 2 == 0){
        updateCheckboxGroupInput(session  = session, 
                                 inputId  = "country",
                                 choices  =  unique(as.character(anim.data$Country)),
                                 selected =  unique(as.character(anim.data$Country)))
        updateCheckboxGroupInput(session  = session, 
                                 inputId  = "A",
                                 choices  = unique(subset(anim.data, anim.data$WHO_REGION == "EMR")$Country),
                                 selected = unique(subset(anim.data, anim.data$WHO_REGION == "EMR")$Country))
        updateCheckboxGroupInput(session  = session, 
                                 inputId  = "B",
                                 choices  = unique(subset(anim.data, anim.data$WHO_REGION == "EUR")$Country),
                                 selected = unique(subset(anim.data, anim.data$WHO_REGION == "EUR")$Country))
        updateCheckboxGroupInput(session  = session, 
                                 inputId  = "C",
                                 choices  =  unique(subset(anim.data, anim.data$WHO_REGION == "AFR")$Country),
                                 selected = unique(subset(anim.data, anim.data$WHO_REGION == "AFR")$Country))
        updateCheckboxGroupInput(session  = session, 
                                 inputId  = "D",
                                 choices  =  unique(subset(anim.data, anim.data$WHO_REGION == "AMR")$Country),
                                 selected = unique(subset(anim.data, anim.data$WHO_REGION == "AMR")$Country))
        updateCheckboxGroupInput(session  = session, 
                                 inputId  = "E",
                                 choices  =  unique(subset(anim.data, anim.data$WHO_REGION == "WPR")$Country),
                                 selected = unique(subset(anim.data, anim.data$WHO_REGION == "WPR")$Country))
        updateCheckboxGroupInput(session  = session, 
                                 inputId  = "F",
                                 choices  =  unique(subset(anim.data, anim.data$WHO_REGION == "SEAR")$Country),
                                 selected = unique(subset(anim.data, anim.data$WHO_REGION == "SEAR")$Country))
        
                                # choices  = subset(unique(as.character(anim.data$Country)), unique(as.character(anim.data$Country)) != "None"),
                                # selected = subset(unique(as.character(anim.data$Country)), unique(as.character(anim.data$Country)) != "None"))
        
      } else {
        updateCheckboxGroupInput(session  = session, 
                                 inputId  = "country",
                                 choices  = unique(as.character(anim.data$Country)),
                                 selected = "")
        updateCheckboxGroupInput(session  = session, 
                                 inputId  = "A",
                                 choices  = unique(subset(anim.data, anim.data$WHO_REGION == "EMR")$Country),
                                 selected = "")
        updateCheckboxGroupInput(session  = session, 
                                 inputId  = "B",
                                 choices  = unique(subset(anim.data, anim.data$WHO_REGION == "EUR")$Country),
                                 selected = "")
        updateCheckboxGroupInput(session  = session, 
                                 inputId  = "C",
                                 choices  =  unique(subset(anim.data, anim.data$WHO_REGION == "AFR")$Country),
                                 selected = "")
        updateCheckboxGroupInput(session  = session, 
                                 inputId  = "D",
                                 choices  =  unique(subset(anim.data, anim.data$WHO_REGION == "AMR")$Country),
                                 selected = "")
        updateCheckboxGroupInput(session  = session, 
                                 inputId  = "E",
                                 choices  =  unique(subset(anim.data, anim.data$WHO_REGION == "WPR")$Country),
                                 selected = "")
        updateCheckboxGroupInput(session  = session, 
                                 inputId  = "F",
                                 choices  =  unique(subset(anim.data, anim.data$WHO_REGION == "SEAR")$Country),
                                 selected = "")
      }
  }
  
  })
  
  yearData <- reactive({
    
    if(input$radio == "AFR"){
      #j = which(regions == input$radio)
      p = anim.data
      p$Incidence  =  round(sqrt(p$Incidence),2)
      df<- p  %>%
        filter(Year == input$year & WHO_REGION == "AFR" & Country %in% input$C) %>%
        select(Country, Coefficient.of.Variation, Incidence, Mean.vaccination, Mean.birth.rate) 
    }  else if (input$radio == "EMR"){
      p = anim.data
      p$Incidence  =  round(sqrt(p$Incidence),2)
      df<- p %>%
        filter(Year == input$year & WHO_REGION == "EMR" & Country %in% input$A) %>%
        select(Country, Coefficient.of.Variation, Incidence, Mean.vaccination, Mean.birth.rate)  
    }else if (input$radio == "EUR"){
      p = anim.data
      p$Incidence  =  round(sqrt(p$Incidence),2)
      df<- p %>%
        filter(Year == input$year & WHO_REGION == "EUR" & Country %in% input$B) %>%
        select(Country, Coefficient.of.Variation, Incidence, Mean.vaccination, Mean.birth.rate)  
    } else if (input$radio == "AMR"){
      p = anim.data
      p$Incidence  =  round(sqrt(p$Incidence),2)
      df<- p %>%
        filter(Year == input$year & WHO_REGION == "AMR" & Country %in% input$D) %>%
        select(Country, Coefficient.of.Variation, Incidence, Mean.vaccination, Mean.birth.rate)  
    } else if (input$radio == "WPR"){
      p = anim.data
      p$Incidence  =  round(sqrt(p$Incidence),2)
      df<- p %>%
        filter(Year == input$year & WHO_REGION == "WPR" & Country %in% input$E) %>%
        select(Country, Coefficient.of.Variation, Incidence, Mean.vaccination, Mean.birth.rate)  
    } else if (input$radio == "SEAR"){    
      
      p = anim.data
      p$Incidence  =  round(sqrt(p$Incidence),2)
      df<- p %>%
        filter(Year == input$year & WHO_REGION == "SEAR" & Country %in% input$F) %>%
        select(Country, Coefficient.of.Variation, Incidence, Mean.vaccination, Mean.birth.rate)  
    }  else {
      p = anim.data
      p$Incidence  =  round(sqrt(p$Incidence),2)
      df<- p %>%
        filter(Year == input$year & Country %in% input$country) %>%
        
        select(Country, Coefficient.of.Variation, Incidence, Mean.vaccination, Mean.birth.rate)
      
    }
  #  if(input$radio == "AFR"){
  #    p = anim.data
  #    p$Incidence  =  round(sqrt(p$Incidence),2)
  #    df<- p %>%
  #    filter(Year == input$year & WHO_REGION == "AFR") %>%
  #      
  #     select(Country, Coefficient.of.Variation, Incidence, Mean.vaccination, Mean.birth.rate) 
        
    #  updateCheckboxGroupInput(session  = session, 
    #                           inputId  = "country",
    #                           choices  = subset(anim.data$Country, anim.data$WHO_REGION == "AFR"),
    #                           selected = subset(anim.data$Country, anim.data$WHO_REGION == "AFR"))
    #} else if (input$radio == "SEAR"){
    #  p = anim.data
    #  p$Incidence  =  round(sqrt(p$Incidence),2)
    #  df<- p %>%
    #    filter(Year == input$year & WHO_REGION == "SEAR") %>%
        
    #    select(Country, Coefficient.of.Variation, Incidence, Mean.vaccination, Mean.birth.rate)
    #}
    
    
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