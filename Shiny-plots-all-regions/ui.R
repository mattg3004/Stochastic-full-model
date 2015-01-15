

library(googleCharts)

# Use global max/min for axes so the view window stays
# constant as the user moves between years
xlim <- list(
  min = 0, #min(anim.data$Health.Expenditure) - 500,
  max = max(anim.data$Coefficient.of.Variation, na.rm = TRUE) + 0.2
)
ylim <- list(
  #min = 0,  #min(anim.data$Life.Expectancy),
  min = sqrt(min(anim.data$Incidence, na.rm = TRUE)) -0.2,
  max = sqrt(max(anim.data$Incidence, na.rm = TRUE)) + 0.2
  #max = 50
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
  
  h2("Measles data 1980 - 2014"),
  
  sidebarLayout(
  sidebarPanel(
    
    actionButton("selectall", label="Select/Deselect all"),
    
    
    radioButtons(
      inputId = "radio",
      label   = "Variable Selection Type:",
      choices =  c("All", subset(unique(as.character(anim.data$WHO_REGION)), unique(as.character(anim.data$WHO_REGION)) !="") ),# subset(unique(as.character(anim.data$WHO_REGION)), unique(as.character(anim.data$WHO_REGION)) !="")),#,unique(as.character(anim.data$WHO_REGION))),
      selected="All"),

 #   conditionalPanel(
#      condition = "input$radio == regions[1]",
#        checkboxGroupInput(
#        inputId = "country",  # ID to be used in server.R
#        label="Select Region:",
#        choices= unique(as.character(anim.data$WHO_REGION)), 
#        selected = countries.to.plot
#      ),
#    ),
    
  # conditionalPanel(
  #   condition = "input.radio == 'Manual Select'",
  #    checkboxGroupInput(
  #      inputId = "country",  # ID to be used in server.R
  #      label="Select Countries:",
  #      choices= unique(as.character(anim.data$Country)), 
  #      selected = c()
  #    )
  #  )#,
    
    checkboxGroupInput(inputId  =  "country",  # ID to be used in server.R
                      label    =  "Select Countries:",
                       choices  =  unique(as.character(anim.data$Country)),
                       #choices  =  subset(unique(as.character(anim.data$Country)), unique(as.character(anim.data$Country)) != "None"),
                       selected =  unique(as.character(anim.data$Country))
     )
    ),
  mainPanel(googleBubbleChart("chart",
                              width="100%", height = "600px",
                              # Set the default options for this chart; they can be
                              # overridden in server.R on a per-update basis. See
                              # https://developers.google.com/chart/interactive/docs/gallery/bubblechart
                              # for option documentation.
                              options = list(
                                fontName = "Source Sans Pro",
                                fontSize = 12,
                                # Set axis labels and ranges
                                hAxis = list(
                                  title = "Coefficient of variation",
                                  viewWindow = xlim
                                ),
                                vAxis = list(
                                  title = "Sqrt of incidence per 1000",
                                  viewWindow = ylim,
                                  #ticks = c(0,4,16,36,64, 100, 164),
                                  #maxValue = max(anim.data$Incidence, na.rm = TRUE) + 2,
                                  #minValue = min(anim.data$Incidence, na.rm = TRUE),
                                 # viewWindow.min  =  min,
                                 # viewWindow.max  =  max,
                                  logScale = F
                                ),
                                
                                colorAxis = list(
                                  colors = c('firebrick', 'red', 'yellow', 'palegreen', 'darkgreen'),
                                  colors = c('firebrick', 'red', 'yellow', 'skyblue', 'blue')
                                ),
                                sizeAxis = list(minValue = 2,  maxSize = 25),
                                
                                # The default padding is a little too spaced out
                                chartArea = list(
                                  top = 50, left = 75,
                                  height = "75%", width = "75%"
                                ),
                                # Allow pan/zoom
                                explorer = list(),
                                # Set bubble visual props
                                bubble = list(
                                  opacity = 0.55, stroke = "none",
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
  shiny::column(2, offset = 4,
                sliderInput(inputId = "year", 
                            label = sprintf("Beginning year of %s period", window.length),
                            min = min(anim.data$Year), max = max(anim.data$Year), step = 0.1, 
                            format = "####", value = min(anim.data$Year), width = '200%',
                            animate = animationOptions(interval = 25,  loop = FALSE, 
                                                       playButton = 'Play', 
                                                       pauseButton = 'Pause'))
  )
  
    # Outputs excluded for brevity 
  )
  
)#,
  
  
#  fluidRow(
#    shiny::column(2, offset = 4,
#                  sliderInput("year", "Year",
#                              min = min(anim.data$Year), max = max(anim.data$Year), step = 0.1,
#                              value = min(anim.data$Year), animate = TRUE),
                  

 #                   checkboxGroupInput(inputId="country",  # ID to be used in server.R
#                                       label="Select Countries:",
#                                       choices=unique(as.character(anim.data$Country)),
#                                       selected = countries.to.plot 
#                  )
                  
                  #mainPanel(
                  #  htmlOutput("plot")  # Argument name from server.R
                  #)
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
  #  )
  #)
))


