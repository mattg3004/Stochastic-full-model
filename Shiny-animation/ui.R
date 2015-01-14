  library(googleCharts)
  
  # Use global max/min for axes so the view window stays
  # constant as the user moves between years
  xlim <- list(
    min = 0, #min(anim.data$Health.Expenditure) - 500,
    max = max(anim.data$Coefficient.of.Variation, na.rm = TRUE) + 0.05
  )
  ylim <- list(
    min = 0,  #min(anim.data$Life.Expectancy),
    max = max(anim.data$Incidence, na.rm = TRUE) + 2
  )
  
  
  
  
  
  shinyUI(fluidPage(
    # This line loads the Google Charts JS library
    googleChartsInit(),
    
    # Use the Google webfont "Source Sans Pro"
    tags$link(
      href=paste0("http://fonts.googleapis.com/css?",
                  "family=Source+Sans+Pro:300,600,300italic"),
      rel="stylesheet", type="text/css"),
    tags$style(type="text/css",
               "body {font-family: 'Source Sans Pro'}"
    ),
    
    h2("Google Charts demo"),
    
    
    googleBubbleChart("chart",
                      width="100%", height = "700px",
                      # Set the default options for this chart; they can be
                      # overridden in server.R on a per-update basis. See
                      # https://developers.google.com/chart/interactive/docs/gallery/bubblechart
                      # for option documentation.
                      options = list(
                        fontName = "Source Sans Pro",
                        fontSize = 13,
                        # Set axis labels and ranges
                        hAxis = list(
                          title = "Coefficient of variation",
                          viewWindow = xlim
                        ),
                        vAxis = list(
                          title = "Incidence per 1000",
                          viewWindow = ylim
                        ),
                        
                          colorAxis = list(
                            colors = c('purple', 'green')
                        ),
                        
                        # The default padding is a little too spaced out
                        chartArea = list(
                          top = 50, left = 75,
                          height = "75%", width = "75%"
                        ),
                        # Allow pan/zoom
                        explorer = list(),
                        # Set bubble visual props
                        bubble = list(
                          opacity = 0.4, stroke = "none",
                          # Hide bubble label
                          textStyle = list(
                           
                          )
                        ),
                        # Set fonts
                        titleTextStyle = list(
                          fontSize = 16
                        ),
                        tooltip = list(
                          textStyle = list(
                            fontSize = 12
                          )
                        )
                      )
    ),
    fluidRow(
      shiny::column(2, offset = 4,
                    sliderInput("year", "Year",
                                min = min(anim.data$Year), max = max(anim.data$Year),
                                value = min(anim.data$Year), animate = TRUE)
                   # checkboxGroupInput("variable", "Variable:", choices = names(anim.data$Country)),
                 #  selectInput("country", 
                  #              "Country:", 
                  #                unique(as.character(anim.data$Country)), multiple=TRUE, selectize = TRUE)
                  # server = function(input, output, session) {
                  #   output$content <- renderTable({
                  #     if(is.null(input$variable))
                  #       return()
                  #     
                  #     anim.data[input$variable]
                  #   })
                  # }
      )
    )
  ))