#' Offshore wind scenario
#'
#' Launches an interactive Shiny application to define custom scenarios for use in the bioenergetic model.
#' 
#' @import shiny
#' @rawNamespace import(shinyjs, except = runExample)
#' @import shinythemes
#' @import shinydashboard
#' @import shinyFeedback
#' @import shinyvalidate
#' @import leaflet
#' @rawNamespace import(dplyr, except = c(first, last, between))
#' @export
#' @author Rob Schick, Phil J. Bouchet, Enrico Pirotta
#' @return An object of class \code{narwscenario}.

scenario <- function() {

  # UI function ----
  
  ui <- function() {
    
    # Header 
    header <-  shinydashboard::dashboardHeader(title = span("BOEM Offshore Wind Scenario",
                                                            style = "color: #fff;"), # Customize title style as needed
                                               tags$li(class = "dropdown",
                                                       # tags$img(src = 'BOEM_logo.png', height = '40px', 
                                                       #          style = "position: absolute; right: 20px; top: 10px;")
                                               ))

    # Side bar
    sidebar <- shinydashboard::dashboardSidebar(
      tags$head(
        tags$style(HTML("
        /* Reduce space between checkbox inputs */
        .checkbox { 
          margin-top: 1px; 
          margin-bottom: 1px; 
        }
        
        /* Change color of checkbox labels */
        .checkbox label { 
          color: #b8c7ce; 
        }
      "))
      ),
      sidebarMenuOutput("sidebar"),
      hr(), 
      tags$head(
        tags$style(HTML("
        .btn-download {
          color: #fff;
          background-color: #5cb85c;
          border-color: #4cae4c;
        }

        .btn-download:hover {
          color: #fff;
          background-color: #449d44;
          border-color: #398439;
        }
        
        .element-padding {
          padding: 0px 10px; 
        }
      "))
      ),
      div(class = "element-padding",
          helpText("Use these checkboxes to plot different geographic features on the map.")
      ),
      checkboxInput("windfarm1", label = 'Wind Farm Boundary', value = F),
      checkboxInput("wearoute1", label = 'Route(s) to WEA', value = F),
      checkboxInput("piles", label = 'Turbine Locations', value = F),
      checkboxInput('tss', label = 'Traffic Separation Schemes', value = F),
      checkboxInput('atba', label = 'Areas to Be Avoided', value = F),
      checkboxInput('speed_zones', label = 'NARW Speed Zones', value = F),
      div(class = "element-padding",
          downloadButton("download", "Download Scenario", class = "btn-download")
      )
    )
    
    # Body 
    body <- shinydashboard::dashboardBody(
      shinyjs::useShinyjs(),
      shinyFeedback::useShinyFeedback(),
      fluidRow(
        
        
        # Sidebar 
        column(width = 5,
               
               # Map editor
               tabItems(
                 
                 # User Guide
                 tabItem(tabName = 'about',
                         h3('How to use the App'),
                         helpText("The BOEM offshore wind scenario tool allows developers, regulators, and managers to simulate realistic scenarios of wind farm construction in support of a bioenergetic model of critically endangered North Atlantic right whales (Eubalaena glacialis). The model is available as a software package implemented in the R programming language; further information and a detailed tutorial can be found on the package website at https://pjbouchet.github.io/narwind/."),
                         h3("General User Interface (GUI)"),
                         helpText("The tool’s GUI consists of a tabset panel for scenario parameterization (left) and an interactive map for on-the-fly visualization (right). The “Map” tab allows users to show/hide spatial features of interest, including wind energy areas (WEAs), turbine locations, and vessel transit routes. Zoom/panning controls can be used to navigate through the map interactively."),
                         h3("Using the tool"),
                         helpText("An example scenario involving three sites (2 x Southern New England, 1 x Virginia) is pre-loaded and shown by default. Bespoke scenarios can be defined by uploading custom data  and modifying scenario parameters using the controls found in the Piling Noise and Vessel Traffic tabs."),
                         hr(),
                         h3("Piling Noise"),
                         helpText("The acoustic footprint of pile-driving can be visualized as a circular buffer centered on a given piling location. Buffer boundaries mark the distance at which noise levels have decreased to a target value, under the assumption that transmission loss (TL) depends on log-range (log10R) and frequency-specific absorption. Ticking the “Show noise footprint” checkbox will add buffers to the map and activate additional controls, including: "),
                         h4("Target Noise Level"),
                         helpText("This slider controls the size of the footprint shown; a larger value means a smaller buffer, as noise levels are highest at the piling site and decline with increasing distance from the source."),
                         h4("Piling Date"),
                         helpText("This slider can be used to cycle through the calendar year and display wind farm scenario parameters for a chosen day. Use the play button to launch a day-by-day animation."),
                         h4("Log-range Coefficient"),
                         helpText("As a sound wave propagates from a localized source, its energy spreads over an increasingly larger area and its intensity therefore declines. This parameter controls the rate at which spreading loss occurs. Higher values capture stronger range attenuation."),
                         h4("Absorption Coefficient"),
                         helpText("This parameter determines the rate at which sound energy is transformed into heat as the sound wave propagates away from the source, inducing friction between water molecules. The value of the absorption coefficient (in dB/km) depends on the frequency of the sound, and likely varies with water properties such as temperature."),
                         h4("Noise Mitigation"),
                         helpText("This parameter controls the reduction in source level achieved using noise abatement systems such as bubble curtains"),
                         hr(),
                         h3('Vessel Traffic'),
                         helpText("TBD")
                         
                         
                         
                         
                 ),
                 
                 # Data tab 
                 
                 # Select Data    
                 tabItem(tabName = 'map',
                         
                         selectInput("scenario", label = "Offshore Wind Scenario", 
                                     choices = c("Scenario 1", 
                                                 "Scenario 2",
                                                 "Scenario 3", 
                                                 "Custom Scenario"),
                                     selected = "Scenario 1"),
                         conditionalPanel(condition = "input.scenario == 'Custom Scenario'",
                                          fileInput("custom_shp", "Upload Your WEA Boundary Shape File", 
                                                    accept = c(".shp", ".dbf", ".sbn", ".sbx", ".shx", ".prj"), 
                                                    multiple = TRUE),
                                          fileInput("custom_piles", "Upload Your Piling Locations & Timing Shape File", 
                                                    accept = c(".shp", ".dbf", ".sbn", ".sbx", ".shx", ".prj"), 
                                                    multiple = TRUE),
                                          fileInput("custom_route", "Upload Your Vessel Route(s) Shape File", 
                                                    accept = c(".shp", ".dbf", ".sbn", ".sbx", ".shx", ".prj"), 
                                                    multiple = TRUE),
                                          selectInput("custom_phase", label = 'Development Phase',
                                                      choices = c('Construction', 'Operation and Management'))
                         )
                         
                         
                 ),
                 
                 tabItem(tabName = "piling",
                         uiOutput("dynamicTab")
                 ),
                 
                 # Mitigation of Sound - Table
                 tabItem(tabName = 'vessel_traffic',
                         
                         # helpText("Plot vessel traffic density?"),
                         checkboxInput("vessel_rast_yesno", label = 'Vessel Density', value = FALSE),
                         conditionalPanel(
                           condition = "input.vessel_rast_yesno == true",
                           radioButtons("vessel_scenario",
                                        label = 'Choose Which Vessel Density Layer',
                                        choices = c('Baseline' = 'baseline',
                                                    'Scenario' = 'scenario'),
                                        selected = character(0)),
                           # checkboxGroupInput("vessel_scenario",
                           #                    label = 'Choose Which Vessel Density Layer',
                           #                    choices = c('Baseline' = 'baseline',
                           #                                'Altered Scenario' = 'scenario')),
                         ), # end conditional panel
                         selectInput("month", label = 'Choose a Month to Plot',
                                     choices = c("Jan", "Feb", "Mar", "Apr",
                                                 "May", "Jun", "Jul", "Aug",
                                                 "Sep", "Oct", "Nov", "Dec")),
                         useShinyjs(),  # Initialize shinyjs
                         tags$head(
                           tags$script(
                             HTML("
                                   $(document).on('click', '#saveBtn', function() {
                                   $(this).css('background-color', '#2ca25f');  // Change button color to a green from color brewer
                                   var btn = $(this);
                                   setTimeout(function() {
                                     btn.css('background-color', '');  // Reset to default color after delay
                                   }, 4000);  // Delay in milliseconds
                                 });
                               ")
                           )
                         ),
                         DT::DTOutput("editableTable"),
                         #   uiOutput("tabContent"),
                         actionButton("saveBtn", "Save Changes")
                         
                         
                         
                         
                 ),
                 
                 # Marine Mammals
                 tabItem(tabName = "whale_density",
                         
                         helpText("Plot whale density?"),
                         checkboxInput("mm_rast", label = 'Whale Density', value = F),
                         selectInput("mm_month", label = 'Choose a Month to Plot',
                                     choices = c("Jan", "Feb", "Mar", "Apr",
                                                 "May", "Jun", "Jul", "Aug",
                                                 "Sep", "Oct", "Nov", "Dec"),
                                     selected = "Apr")

                 )

               )

        ),
        
        # Main Display
        column(width = 6,
               
               # Map
               box(width = NULL, solidHeader = T, collapsible = T, 
                   status = 'primary', title = 'Map', 
                   
                   
                   leaflet::leafletOutput("map", height = 600),
                   
                   helpText("These data are preliminary data, subject to change,
                        and not to be used without permission from the contributor(s)")
                   
               )
               
        ) # end column main display
        
      ),
      
    )
    
    # Construct UI 
    shinydashboard::dashboardPage(
      header,
      sidebar,
      body
    )
  }
 
  # Server function----
  
  server <- function(input, output, session) {
    
    # Load required scripts and data objects
    
    # Coordinate Reference System 
    brs_crs <- paste("+proj=aea +lat_1=27.33333333333333",
                     "+lat_2=40.66666666666666 +lat_0=34 +lon_0=-78",
                     "+x_0=0 +y_0=0",
                     "+ellps=WGS84 +datum=WGS84 +units=m +no_defs")
    
    output$sidebar <- renderMenu({
      
      sidebarMenu(
        
        id = "tabs",
        menuItem("About", tabName = "about", icon = icon("circle-info")),
        menuItem("Scenario Selection", tabName = "map", icon = icon("map-marker")),
        menuItem("Piling Noise", tabName = "piling", icon = icon("volume-up")),
        menuItem("Vessel Traffic", tabName = "vessel_traffic", icon = icon("ship")),
        menuItem("Whale Density", tabName = "whale_density", icon = icon("signal"))
        
      )
      
    })
    
    # Reactive expression for vesselmaps
    vesselmaps <- reactiveVal(map_vessels(obj = scenario_01,
                                          z = "dist",
                                          strike_scalar = 1e-7,
                                          which.month = 1:12,
                                          vessel.speed = NULL,
                                          speed.limit = c(10,10,10,10,rep(NA,8)),
                                          baseline = TRUE,
                                          do.plot = FALSE,
                                          spgdf = FALSE))
    
    # Output Tab Content for Vessel Table 
    output$tabContent <- renderUI({
      if(input$tabs == "vessel_traffic") {
        DT::DTOutput("editableTable")
      } 
    })
    
    # Reset Values 
    observeEvent(input$reset_input, {
      shinyjs::reset(id = "")
    })
    
    # Shiny Feedback 
    iv <- InputValidator$new()
    iv$add_rule("alpha", sv_between(0, 30))
    iv$add_rule("abs_coef", sv_between(0, 5))
    iv$enable()
    
    # Transmission Loss Functions 
    dB2km <- function(target.dB, SL = 200, logfac = 15, a = 1.175){
      opt.fun <- function(r, SL, target.L, logfac){
        return((SL-TL(r, logfac = logfac, a = a)-target.L)^2)}
      out <- optimize(f = opt.fun, interval = c(0,30000), 
                      SL = SL, 
                      target.L = target.dB, 
                      logfac = logfac)
      return(out$minimum)
    }
    
    TL <- function(r, logfac = 15, a = 1.175){
      # r in km
      # alpha is in dB/km
      loss <- logfac * log10(r*1000)
      loss[loss < 0] <- 0
      loss <- loss + a * r
      return(loss)
    }
    
    
    # Scenario Selector
    uploadedDataVal <- reactiveVal()
    
    selectedList <- reactive({
      
      if (input$scenario == "Scenario 1") {
        return(scenario_01)  
      } else if (input$scenario == "Scenario 2") {
        return(scenario_02)  
      } else if (input$scenario == "Scenario 3") {
        return(scenario_03)  # for the third scenario
      } else {
        
        validate(need(tools::file_ext(input$file$name) == "rda", 
                      "Please upload a valid .rda file"))
        
        # Create a new environment for loading
        load_env <- new.env()
        
        # Load the .rda file into the new environment
        load(input$file$datapath, envir = load_env)
        
        # Get the names of all objects in the new environment
        loaded_objs <- ls(load_env)
        
        # Ensure only one object is loaded
        validate(need(length(loaded_objs) == 1, "The .rda file should contain exactly one object. Please remove extra objects and try again."))
        
        # Fetch the actual object from the new environment
        uploadedData <- get(loaded_objs[[1]], envir = load_env)
        
        # Update the reactiveVal with the new data
        uploadedDataVal(uploadedData)
        
        # Further checks, e.g. validate the structure of uploadedData
        validate(need(is.list(uploadedData) && "locs" %in% names(uploadedData), "The uploaded object must be a list with a 'locs' element."))
        
        return(uploadedData)
        
      }
    })
    
    # Construction & Operation Phase
    currentPhase <- reactive({
      selectedList()[["phase"]]
    })
    
    # Dynamic Sidebar Tab
    output$dynamicTab <- renderUI({
      if (currentPhase() == 1) {
        tabItem(tabName = 'piling',
                h4('Acoustic Footprints'),
                helpText("The acoustic",
                         " footprint of pile-driving can be visualized",
                         " as a circular buffer on a given piling location. ",
                         "Buffer boundaries mark the distance at which noise levels ",
                         "have decreased to a target value, based on a simple ",
                         "transmission loss model."),
                # br(),
                h4('Visualize Footprints'),
                checkboxInput("show_buffer", "Show Noise Footprint", value = FALSE),
                sliderInput("bufferRadius", "Target Noise Level (dB)", min = 1, max = 200,
                            value = 100),
                sliderInput(
                  "dates",
                  "Piling date",
                  min = as.Date("2024-01-01"),
                  max = as.Date("2024-04-24"),
                  step = 1,
                  value = as.Date("2024-01-01"),
                  animate =
                    animationOptions(interval = 350, loop = FALSE)
                ),
                
                hr(),
                h4('Change Source Level Values'),
                numericInput("sl", 'Source Level (dB)', min = 10, max = 265, value = 220),
                hr(),
                h4('Transmission Loss'),
                helpText("We use a simple model, which assumes that transmission ",
                         "loss (TL) depends on log-range (log10R) and  ",
                         "frequency-specific absorption (a).  ",
                         "The model is of the form: TL = b x log10R + a x R. "),
                numericInput("alpha", "Log-range coefficient (b)",
                             min = 10, max = 20, value = 15),
                numericInput("abs_coef", "Absorption coefficient (a)",
                             min = 0, max = 5, value = 1.175),
                hr(),
                h4('Noise mitigation'),
                helpText("Noise abatement systems (e.g., bubble curtains)  ",
                         "may be in place to limit the exposure of ",
                         "protected marine species to pile-driving noise."),
                numericInput("bubble", "Magnitude of noise attenuation (dB)",
                             value = 0,
                             min = 0,
                             max = 10,
                             step = 1),
                
                actionButton("reset_input", "Reset inputs")
        )
      } else if (currentPhase() == 2) {
        tabItem(tabName = 'piling',
                h4('Acoustic Footprints'),
                helpText("You've chosen an operation phase, so there are no changes to be made on this panel")
        )
      } else {
        
        tabItem(tabName = 'piling', "No details available")
      }
    })
    
    # Turbines 
    turbine_df <- reactive({
      # Check if 'locs' exists and is not empty
      if ("locs" %in% names(selectedList()) && 
          nrow(selectedList()[["locs"]]) > 0 &&
          !is.null(selectedList()[["locs"]]$date)) {
        data <- selectedList()[["locs"]] %>%
          mutate(date = as.Date(date, format = '%Y-%m-%d'),
                 month = lubridate::month(date, label = TRUE, abbr = TRUE),
                 buffer = 80) %>%
          sf::st_as_sf(coords = c('longitude', 'latitude')) %>%
          sf::st_set_crs(4326)
        
        return(data)
      } else {
        # Return NULL or an alternative value when 'locs' is not suitable
        return(NULL)
      }
    })
    
    
    # Update Dates based on Selector 
    observe({
      turbine_data <- req(turbine_df())
      
      # Calculate min and max dates
      date_range <- range(turbine_data$date)
      min_date <- date_range[1]
      max_date <- date_range[2]
      
      # Update the slider input
      updateSliderInput(session, 'dates', value = min_date,
                        min = min_date, max = max_date)
      
    })
    
    
    # Base Map  
    initial_lat = 38.2
    initial_lng = -73.7
    initial_zoom = 6
    
    output$map <- renderLeaflet({
      leaflet() %>%
        addWMSTiles(baseUrl = "https://server.arcgisonline.com/ArcGIS/rest/services/Ocean/World_Ocean_Base/MapServer/tile/{z}/{y}/{x}", layers = "Ocean_World_Ocean_Base") %>%
        addWMSTiles(baseUrl = "https://server.arcgisonline.com/ArcGIS/rest/services/Ocean/World_Ocean_Reference/MapServer/tile/{z}/{y}/{x}",
                    layers = "Ocean_World_Ocean_Reference", attribution = 'Tiles &copy; Esri &mdash; Sources: GEBCO, NOAA, CHS, OSU, UNH, CSUMB, National Geographic, DeLorme, NAVTEQ, and Esri') %>%
        setView(lat = initial_lat, lng = initial_lng, zoom = initial_zoom) %>%
        addMeasure()
    })
    
    
    
    # Reset Zoom Observer 
    observe({
      input$reset_button
      leafletProxy("map") %>%
        setView(lat = initial_lat, lng = initial_lng, zoom = initial_zoom)
    })
    
    
    # WEA Zone 1 Observer
    
    observe({
      
      # Define proxy and remove effort group
      proxy <- leafletProxy("map")
      proxy %>% clearGroup("windfarm1")
      
      # Add Wind Energy Areas
      if(input$windfarm1){
        
        proxy %>%
          addPolygons(data = windfarm1,
                      fill = FALSE,
                      weight = 2,
                      color = "#006d2c",
                      popup = 'Wind Farm 1',
                      group = 'windfarm1') %>%
          addPolygons(data = windfarm2,
                      fill = FALSE,
                      weight = 2,
                      color = "#006d2c",
                      popup = 'Wind Farm 2',
                      group = 'windfarm1')  %>%
          addPolygons(data = windfarm3,
                      fill = FALSE,
                      weight = 2,
                      color = "#006d2c",
                      popup = 'Wind Farm 3',
                      group = 'windfarm1')
        
        # Switch to show/hide
        ifelse(input$windfarm1, showGroup(proxy, 'windfarm1'),
               hideGroup(proxy, 'windfarm1'))
        
      }
      
    })
    
    # User Guide Text
    output$help_text <- renderUI({
      if (input$help_button > 0) {
        tagList(
          tags$script(HTML('$(document).ready(function() { $("#help_button").popover(); })')),
          "The BOEM offshore wind scenario tool allows developers, regulators, and managers
        to simulate realistic scenarios of wind farm construction in support of a
        bioenergetic model of critically endangered North Atlantic right whales (Eubalaena glacialis).
        The model is available as a software package implemented in the R programming language;
        further information and a detailed tutorial can be found on the package website at
        https://pjbouchet.github.io/narwind/."
        )
      }
    })
    
    
    # WEA Piling Observer
    
    observe({
      
      req(turbine_df())
      
      # Define proxy and remove effort group
      proxy <- leafletProxy("map")
      proxy %>% clearGroup("piles")
      
      # Add Piling Locations within Wind Energy Areas
      if(input$piles){
        
        proxy %>%
          addCircleMarkers(data = turbine_df(),
                           radius = 2,
                           fill = FALSE,
                           weight = 2,
                           stroke = TRUE,
                           color = "#08306b",
                           group = 'piles')
        
      }
      
    })
    
    # Make buffers on the fly; use addCircles with radius in meters
    observe({
      
      req(turbine_df(), input$show_buffer, input$bufferRadius, input$sl, input$bubble, input$alpha, input$abs_coef)
      
      # Since req() ensures that inputs are not NULL, we can safely calculate buffer_radius
      buffer_radius <- if (input$show_buffer) {
        dB2km(input$bufferRadius, SL = input$sl - input$bubble, logfac = input$alpha, a = input$abs_coef) * 1000
      } else {
        0  # Set to zero if show_buffer is FALSE
      }
      
      print(paste("Buffer radius:", buffer_radius))
      
      pt_opacity <- ifelse(input$show_buffer, 0.5, 0)
      stroke_val <- ifelse(input$show_buffer, TRUE, FALSE)
      
      # Applying transformations and optional buffering
      scn1_poly <- turbine_df() %>%
        sf::st_transform(brs_crs) 
      
      if (buffer_radius > 0) {
        scn1_poly <- sf::st_buffer(scn1_poly, dist = buffer_radius)
      }
      
      scn1_poly <- scn1_poly %>%
        sf::st_transform(4326) %>%
        dplyr::filter(date == as.Date(input$dates))  # Ensure date formats are consistent
      
      
      scn1_pt <- turbine_df() %>%
        filter(date == input$dates)
      
      leafletProxy("map") %>%
        clearShapes() %>%
        addPolygons(data = scn1_poly,
                    fillColor = "#fdbe85",
                    color = "#a63603",
                    weight = 1,
                    opacity = 1,
                    fillOpacity = 0.2
        ) %>%
        clearGroup('pile_centroid') %>%
        addCircleMarkers(data = scn1_pt,
                         radius = 0,
                         color = "white",
                         stroke = stroke_val,
                         fillOpacity = pt_opacity,
                         group = 'pile_centroid'
        )
      
    })
    
    
    
    # WEA Ship Routes Observer
    observe({
      
      # Define proxy and remove effort group
      proxy <- leafletProxy("map")
      proxy %>% clearGroup("wearoute1")
      
      # Add Routes to Wind Energy Areas
      if(input$wearoute1){
        
        proxy %>%
          addPolylines(data = route_wf1,
                       fill = FALSE,
                       weight = 2,
                       color = "#636363",
                       group = 'wearoute1',
                       label = "WEA - Route 1") %>%
          addPolylines(data = route_wf2,
                       fill = FALSE,
                       weight = 2,
                       color = "#636363",
                       group = 'wearoute1',
                       label = "WEA - Route 2") %>%
          addPolylines(data = route_wf3,
                       fill = FALSE,
                       weight = 2,
                       color = "#636363",
                       group = 'wearoute1',
                       label = "WEA - Route 3")
        
      }
      
    })
    
    # Marine Mammal Density Observer
    observe({
      
      # Define proxy and remove effort group
      proxy <- leafletProxy("map")
      proxy %>%
        clearGroup("mm")
      
      if(input$mm_rast){
        
        my_raster <- density_narw[[input$mm_month]]
        r <- raster::raster(my_raster)
        
        r_cut <- format_raster(r)
        minValue <- min(terra::values(r_cut), na.rm = TRUE)
        maxValue <- max(terra::values(r_cut), na.rm = TRUE)
        
        proxy %>%
          clearGroup("mm")%>%
          addRasterImage(r_cut,
                         colors = pals::viridis(n = 10),
                         opacity = 1,
                         group = 'mm')
        
        # Switch to show/hide
        ifelse(input$mm_rast, showGroup(proxy, 'mm_rast'),
               hideGroup(proxy, 'mm_rast'))
        
      } # end marine mammal
      
    })
    
    # Vessel Transit - Tabular
    vessel.tbl <- data.frame(
      VesselClass = c("Cable Lay", "Construction/Crane", "Crew Transfer",
                      "Heavy Cargo", "Support Vessels", "Survey", "Tugs"),
      NumberOfVessels = 5,
      # SpeedOfVessels = c(rep(99,7)),
      SpeedOfVessels = c(15, 16, 25,
                         15, 15, 30, 14),
      NumberOfTrips = 10
    )
    
    # Reactive value to store the current state of the table
    currentTable <- reactiveVal(vessel.tbl)
    
    # Update currentTable on each edit
    observeEvent(input$editableTable_cell_edit, {
      info <- input$editableTable_cell_edit

      # Get current table data
      tableData <- currentTable()
      
      # Check if tableData is not NULL
      if (!is.null(tableData)) {
        # Update the cell value
        tableData[info$row, info$col] <- info$value
        
        # Update the reactive value
        currentTable(tableData)
      }
    })
    
    # Render the editable table
    output$editableTable <- DT::renderDT({
      currentTable()
    }, editable = TRUE)
    
    
    # When the user clicks the save button
    # the vessels component of the scenario list object is updated
    updatedScenario <- eventReactive(input$saveBtn, {
      
      update.narwscenario(selectedList(), currentTable())
      
    })
    
    observeEvent(input$saveBtn, {
      updatedMaps <- map_vessels(obj = updatedScenario(),
                                 z = "dist",
                                 strike_scalar = 1e-7,
                                 which.month = 1:12,
                                 vessel.speed = NULL,
                                 speed.limit = c(10,10,10,10,rep(NA,8)),
                                 baseline = TRUE,
                                 do.plot = FALSE,
                                 spgdf = FALSE)
      vesselmaps(updatedMaps)
    })
    
    
    # Vessel Traffic Density Observer
    observe({
      
      # Define proxy and remove effort group
      proxy <- leafletProxy("map")
      proxy %>%
        clearGroup("vessel")
      
      if(input$vessel_rast_yesno && !is.null(input$vessel_scenario) && length(input$vessel_scenario) > 0) {
        
        sel_month <- which(names(traffic)== input$month)
        
        if(input$vessel_scenario == 'baseline') {
          my_raster <- vesselmaps()[[sel_month]][[2]]
        } else if (input$vessel_scenario == 'scenario') {
          my_raster <- vesselmaps()[[sel_month]][[1]]
        }
        
        # Old plotting code
        r <- my_raster
        prange <- c(0, 1000)
        
        proxy %>%
          addRasterImage(r,
                         colors = rev(pals::parula(n = 10)),
                         opacity = 1,
                         group = 'vessel')
        
        showGroup(proxy, 'vessel')
        
      } else {
        hideGroup(proxy, 'vessel')
      } # end vessel traffic
      
    })
    
    # tss observer 
    observe(priority = 4, {
      
      # Define proxy
      proxy <- leafletProxy("map")
      proxy %>% clearGroup('tss')
      
      if(input$tss){
        
        # Plot shipping lanes
        proxy %>%
          addPolygons(data = tss,
                      weight = .5,
                      color = 'grey',
                      fillColor = 'grey',
                      options = pathOptions(clickable = F),
                      group = 'tss')
        
        # Switch to show/hide
        ifelse(input$tss, showGroup(proxy, 'tss'), hideGroup(proxy, 'tss'))
      }
      
    })
    
    # Area to be avoided observer 
    observe(priority = 4, {
      
      # Define proxy
      proxy <- leafletProxy("map")
      proxy %>% clearGroup('atba')
      
      if(input$atba){
        
        # Plot Areas to Be Avoided
        proxy %>%
          addPolygons(data = atba,
                      weight = .5,
                      color = '#2ca25f',
                      fillColor = '#2ca25f',
                      options = pathOptions(clickable = F),
                      group = 'atba')
        
        # Switch to show/hide
        ifelse(input$atba, showGroup(proxy, 'atba'), hideGroup(proxy, 'atba'))
      }
      
    })
    
    
    # NARW Speed Zones observer 
    observe(priority = 4, {
      
      # Define proxy
      proxy <- leafletProxy("map")
      proxy %>% clearGroup('speed_zones')
      
      if(input$speed_zones){
        
        # Plot Areas to Be Avoided
        proxy %>%
          addPolygons(data = speed_zones,
                      weight = .5,
                      color = '#a63603',
                      fillColor = '#a63603',
                      options = pathOptions(clickable = F),
                      group = 'speed_zones')
        
        # Switch to show/hide
        ifelse(input$speed_zones, showGroup(proxy, 'speed_zones'), hideGroup(proxy, 'speed_zones'))
      }
      
    })
    
    # Uncheck Map Features
    observeEvent(input$scenario, {
      
      updateCheckboxInput(session, "windfarm1", value = FALSE)
      updateCheckboxInput(session, "windfarm2", value = FALSE)
      updateCheckboxInput(session, "windfarm3", value = FALSE)
      updateCheckboxInput(session, "wearoute1", value = FALSE)
      updateCheckboxInput(session, "wearoute2", value = FALSE)
      updateCheckboxInput(session, "wearoute3", value = FALSE)
      updateCheckboxInput(session, "piles", value = FALSE)
      updateCheckboxInput(session, "tss", value = FALSE)
      updateCheckboxInput(session, "atba", value = FALSE)
      updateCheckboxInput(session, "speed_zones", value = FALSE)
      
      
    })
    
    # Assemble Parameters 
    my_params <- reactive({
      
      data.frame(scenario = input$scenario,
                 ambient = input$bufferRadius,
                 sourceLvL = input$sl,
                 logrange = input$alpha,
                 absorb = input$abs_coef,
                 lowerdB = input$bubble,
                 phase = input$phase
      )
    })
    
    output$download <- downloadHandler(
      filename = function() {
        paste0("BOEM-NARWIND_parameters_", gsub("[[:blank:]]", "", input$scenario, "_", Sys.time()), ".csv")
      },
      content = function(file) {
        vroom::vroom_write(my_params(), file, ",")
      }
    )
    
    user_shp <- reactive({
      req(input$custom_shp)
      Read_Shapefile(input$custom_shp) %>%
        st_transform(crs = 4326)
    })
    
    # Stop app when closing window
    session$onSessionEnded(function() {
      stopApp()
    })
    
  }
  
  # Build Shiny app----
  shinyApp(ui = ui,
           server = server,
           options = list(launch.browser = TRUE))

} # End function
