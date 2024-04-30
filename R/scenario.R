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
    # ui.R
    # header ----------------------------------------------------------------------
    header <-  shinydashboard::dashboardHeader(title = span("BOEM Offshore Wind Scenario", 
                                                            style = "color: #fff;"), 
                                               tags$li(class = "dropdown")) #,
                                                       # tags$img(src = 'www/BOEM_logo.png', height = '40px', 
                                                       #          style = "position: absolute; right: 20px; top: 10px;")))
                                                        
    # Side bar -----------------
    sidebar <- dashboardSidebar(
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
      uiOutput("dynamicCheckbox"),
      checkboxInput("wearoute1", label = 'Route(s) to WEA', value = F),
      uiOutput("dynamicCheckbox_route"),
      checkboxInput("piles", label = 'Turbine Locations', value = F),
      checkboxInput('tss', label = 'Traffic Separation Schemes', value = F),
      checkboxInput('atba', label = 'Areas to Be Avoided', value = F),
      checkboxInput('speed_zones', label = 'NARW Speed Zones', value = F),
      div(class = "element-padding",
          downloadButton("download", "Download Scenario", class = "btn-download")
      )
    )
    
    # body --------------------------------------------------------------------
    
    body <- shinydashboard::dashboardBody(
      shinyjs::useShinyjs(),
      shinyFeedback::useShinyFeedback(),
      fluidRow(
        
        
        # Sidebar ---
        column(width = 5,
               # map editor
               tabItems(
                 
                 # User Guide
                 tabItem(tabName = 'about',
                         h3('How to use the App'),
                         helpText("The BOEM offshore wind scenario Shiny app allows developers, regulators, and managers to generate construction and operational wind farm scenarios to be used as inputs to a bioenergetic model of critically endangered North Atlantic right whales (Eubalaena glacialis, NARW). It also allows visualization of a range of geographical features relevant to the scenarios and to NARW distribution and density. Users can select from three pre-defined scenarios, which can be modified as described below, or upload files to create a custom scenario. The parameters relating to scenarios generated or modified within the app can be downloaded for input into the bioenergetic model. The app is integrated in the narwind package, and the scenarios generated using the app can be used directly in the functions of the package. A detailed tutorial and vignette for the model and associated software package is available at https://offshore-wind.github.io/narwind/index.html, with specific information pertaining to the Shiny app being found in Tutorial #2: https://offshore-wind.github.io/narwind/articles/scenarios.html"),
                         h3("General User Interface (GUI)"),
                         helpText("The app consists of three main structural elements, namely: 1) A sidebar (left) – which is used for navigation, visualization, and scenario set up. 2) A text area (middle) – which displays information about the app, model parameters, and required user inputs. 3) An interactive map (right) – where scenario parameters and input layers can be visualized on the fly. The map and sidebar always stay in place; however, the contents of the main text area change dynamically in response to user selections. The sidebar is split into a series of tabs (top) and checkboxes (bottom). Tabs give access to the app’s features and, with the exception of the Whale Density tab (which is only used for visualization purposes), allow users to modify scenario parameters. By contrast, checkboxes can only be used to toggle map layers on and off, i.e., they have no role in scenario building. The map is interactive and can be zoomed (using either the +/- buttons in the top left or a scroll of the mouse), panned (by clicking and dragging) and minimized (using the ‘-’ symbol on the top right). Distance measurements can also be taken by clicking on the ruler icon in the top right and following the on-screen instructions."),
                         h3("Scenario Selection"),
                         helpText("Three pre-defined scenarios are already included in the package and app and can be selected from the drop-down menu on the Scenario Selection tab.  Full details of the scenarios can be found in the accompanying project report. All three scenarios relate to three windfarm sites (2 x Southern New England (Farms 1 and 2), 1 x Virginia (Farm 3)):"),
                         h4("Scenario 1: Installation – “unmitigated” case"),
                         helpText("This scenario entails the synchronous construction of all three wind farms at times coinciding with expected peaks in right whale abundance within each respective area (Dec–May in SNE, Feb–March in VA). The objective of this scenario is to explore the potential for cumulative effects of multiple installation activities under unmitigated conditions."),
                         h4("Scenario 2: Installation – “mitigated” case"),
                         helpText("This scenario involves shutdowns of installation operations during the main right whale foraging/calving seasons (Nov–Mar). Construction resumes in May at the southward site (Farm 3) and during the late summer to early fall at Farms 1 and 2. In addition, activities at each of the two SNE farms are asynchronous, with noise abatement systems in place to limit noise impacts. Note that a subset of only 60 monopiles are driven at Farm 1 in order to align with previous BOEM studies."),
                         h4("Scenario 3: Operation"),
                         helpText("In this scenario, we assume that all three farms are in simultaneous and continuous operation. The primary footprint for operations is taken to be vessel traffic to and from wind farm sites. While there is limited knowledge of servicing patterns associated with large-scale windfarm operations off the U.S. east coast, information sourced from BOEM about the types, numbers, routes, and speeds of vessels involved in maintenance and operational activities allows simplifying assumptions to be made concerning servicing scope."),
                         helpText("Once a scenario has been selected, the parameters associated with the selected scenario will be populated throughout the rest of the app. The scenario can then be viewed, and the default scenario parameters can be modified using the controls found in the Piling Noise and Vessel Traffic tabs. Modifications to scenario parameters can be saved into an R object that can be used directly in the bioenergetic model. Custom scenarios can also be uploaded and modified on this tab by selecting “Custom Scenario” in the drop-down menu. This prompts the user to upload files for the piling locations and timings, the vessel routes, and the vessel information table that underpin the scenario they wish to explore. An additional drop-down menu will also become available to indicate whether the custom scenario relates to construction or operation activities. This is important, as it ensures that the scenario object is correctly labeled. The custom scenario can then be viewed in the same way as the pre-defined scenarios, and scenario parameters can be modified using the controls found in the Piling Noise and Vessel Traffic tabs. As with the pre-defined scenarios, modifications made to scenario parameters within the app tabs can be saved into an R object that will be used in the bioenergetic model."),
                         h3("Piling Noise"),
                         helpText("The acoustic footprint of pile-driving can be visualized as a circular buffer centered on a given piling location. Buffer boundaries mark the distance at which noise levels have decreased to a target value, under the assumption that transmission loss (TL) depends on log-range (log10R) and frequency-specific absorption. If the selected scenario relates to the construction phase, then ticking the “Show noise footprint” checkbox will add buffers to the map and activate additional controls. If the scenario relates to the operational phase of a wind farm site, then there are no changes to be made on this tab. When the scenario involves a construction phase, then the following can be viewed and modified by the user: "),
                         h4("Target Noise Level"),
                         helpText("This slider controls the size of the footprint shown; a larger value corresponds to a smaller buffer, as noise levels are highest at the piling site and decline with increasing distance from the source. This slider is for visualization only, i.e., it has no effect on any scenario parameters that feed into the bioenergetic model."),
                         h4("Piling Date"),
                         helpText("This slider can be used to cycle through the calendar year and display wind farm scenario parameters for a chosen day. The play button can be used to launch a day-by-day animation. This slider is for visualization only, i.e., it has no effect on any scenario parameters that feed into the bioenergetic model."),
                         h4("Source level value"),
                         helpText("By default, the source level is set to the value associated with the selected scenario or uploaded scenario but can be modified by the user on this tab. Any changes to this number will be included as a scenario parameter value."),
                         h4("Log-range Coefficient"),
                         helpText("As a sound wave propagates from a localized source, its energy spreads over an increasingly larger area and its intensity therefore declines. This parameter controls the rate at which spreading loss occurs. Higher values capture stronger range attenuation. By default, the Log-range coefficient is set to the value associated with the selected scenario or uploaded scenario but can be modified by the user on this tab. Any changes to this number will be included as a scenario parameter value."),
                         h4("Absorption Coefficient"),
                         helpText("This parameter determines the rate at which sound energy is transformed into heat as the sound wave propagates away from the source, inducing friction between water molecules. The value of the absorption coefficient (in dB/km) depends on the frequency of the sound, and likely varies with water properties such as temperature. By default, the Absorption coefficient is set to the value associated with the selected scenario or uploaded scenario but can be modified by the user on this tab. Any changes to this number will be included as a scenario parameter value."),
                         h4("Noise Mitigation"),
                         helpText("This parameter controls the reduction in source level achieved using noise abatement systems such as bubble curtains. By default, the magnitude of noise attenuation is set to the value associated with the selected scenario or uploaded scenario but can be modified by the user on this tab. Any changes to this number will be included as a scenario parameter value."),
                         helpText("At the bottom of this tab there is a button allowing for values to be reset to the default for the scenario. Any user-modified values will be reflected in the final R object that will be used in the bioenergetic model."),
                         h3('Vessel Traffic'),
                         helpText("The vessel traffic tab allows the user to visualize and modify the vessel strike risk surface that is used in the bioenergetic model. Vessel strike risk is evaluated using a simple metric, Total PLETH, which accounts for both the effect of vessel speed on lethality and for the cumulative exposure of whales to all vessels traversing each grid cell. Total PLETH is scaled to obtain daily strike probabilities for use in the agent-based bioenergetic model. Further details on the calculation of this metric can be found in the accompanying project report. Ticking the “Vessel strike risk” checkbox prompts the user to choose whether they wish to display a baseline vessel strike risk layer or a scenario vessel strike risk layer (relating to the scenario selected on the Scenario Selection tab). The user should also select which month they wish to visualize and be aware that, if a month is chosen that doesn’t coincide with the timing of the chosen scenario, then only baseline vessel strike risk will be displayed.  The “Update Map” button should be clicked to display the selected vessel strike risk surface. The table on this tab allows the user to specify the number of vessels (Nvessels), speed in knots (speed_knt) and number of round trips per foundation (roundtrips_foundation) for each class of vessels (vesselclass) associated with each windfarm site and route; please note that other columns in this table should not be modified. Values in the table can be changed by clicking on the relevant cell. Changes can be visualized on the map by clicking the “Update Map” button."),
                         h3('Whale Density'),
                         helpText("Ticking the “Whale density” checkbox overlays whale density on the map according to the month that is selected."),
                         h3('Download'),
                         helpText("Once all changes to scenario parameters have been made, the scenario object can be downloaded using the green Download Scenario button located on the left. The saved scenario will also be automatically loaded into the open R session, for use in the bioenergetic model once the Shiny app window has been closed. The new scenario object loaded in R will be called “scenario_custom”."),
                         h3('Geographical features'),
                         helpText("The bottom of the left-hand sidebar consists of a list of check boxes that allow different geographical features to be visualized on the map. This includes the boundaries of the pre-defined wind farm sites (not available for custom scenarios), the vessel routes to these wind energy areas (also not available for custom scenarios), the turbine locations, as well as the location of traffic separation schemes, areas to be avoided (for NARW), and NARW speed zones.")

                 ),
                 # data tab ----------------------------------------------------------
                 
                 # Select Data    
                 tabItem(tabName = 'map',
                         
                         selectInput("scenario", label = "Offshore Wind Scenario", 
                                     choices = c("Scenario 1", 
                                                 "Scenario 2",
                                                 "Scenario 3", 
                                                 "Custom Scenario"),
                                     selected = "Scenario 1"),
                         conditionalPanel(condition = "input.scenario == 'Custom Scenario'",
                                          fileInput("custom_vessels", "Upload Your Vessel Info Table (CSV File)", 
                                                    accept = c('.csv'), 
                                                    multiple = TRUE),
                                          fileInput("custom_piles", "Upload Your Piling Locations & Timing Data (CSV File)", 
                                                    accept = c(".csv"), 
                                                    multiple = FALSE),
                                          fileInput("custom_route", "Upload Your Vessel Route(s) Shape File", 
                                                    accept = c(".shp", ".dbf", ".sbn", ".sbx", ".shx", ".prj"), 
                                                    multiple = TRUE),
                                          selectInput("custom_phase", label = 'Development Phase',
                                                      choices = c('Construction' = 1, 
                                                                  'Operation and Management' = 2))
                         )
                         
                         
                 ),
                 
                 tabItem(tabName = "piling",
                         uiOutput("dynamicTab")
                 ),
                 
                 # Mitigation of Sound - Table
                 tabItem(tabName = 'vessel_traffic',
                         h3('How to use this Tab'),
                         helpText("This tab allows the visualization and modification of vessel strike risk surfaces. To view a surface please tick the “Vessel Strike Risk” checkbox, then select “Baseline” or “Scenario” from the drop down menu. After selecting which month to view, please click the “Update Map” button to display the selection. To view and edit the table, please minimize the map by clicking “-“ at the top right of the map. The map can be restored by clicking the “+”. To change a value in the table, please double click in the cell and use the up/down arrows, or enter a value. If table values are changed then please click the “Update Map” button to see those changes reflected in the map."),
                         hr(),
                         # helpText("Plot vessel traffic density?"),
                         checkboxInput("vessel_rast_yesno", label = 'Vessel Strike Risk', value = FALSE),
                         conditionalPanel(
                           condition = "input.vessel_rast_yesno == true",
                           radioButtons("vessel_scenario",
                                        label = 'Choose Which Vessel Strike Risk Layer',
                                        choices = c('Baseline' = 'baseline',
                                                    'Scenario' = 'scenario'),
                                        selected = character(0)),
                         ), # end conditional panel
                         selectInput("month", label = 'Choose a Month to Plot',
                                     choices = c("Jan", "Feb", "Mar", "Apr",
                                                 "May", "Jun", "Jul", "Aug",
                                                 "Sep", "Oct", "Nov", "Dec")),
                         actionButton("saveBtn", "Update Map"),
                         br(), br(),
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
                         DT::DTOutput("editableTable") #,
                         
                         
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
        
        # main Display
        column(width = 6,
               
               # Map
               box(width = NULL, solidHeader = T,collapsible = T, 
                   status = 'primary', title = 'Map', 
                   
                   
                   leafletOutput("map", height = 600),
                   
                   helpText("These data are preliminary data, subject to change,
                        and not to be used without permission from the contributor(s)")
                   
               ) 
               
        ) # end column main display
        
      ),
      
      
    )
    
    
    # construct ui ----------------------------------------------------------
    
    shinydashboard::dashboardPage(
      header,
      sidebar,
      body
    )
    
  }
  
  # Server function -----------
  server <- function(input, output, session) {
    
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
    
    # Reactive expression for vesselmaps --------
    vesselmaps <- reactive({
      # Ensure selectedList() is re-evaluated when it changes
      req(selectedList())  # Use req to handle NULL or not yet available cases
      map_vessels(obj = selectedList(),
                  z = "dist",
                  strike_scalar = 1e-7,
                  which.month = 1:12,
                  vessel.speed = NULL,
                  speed.limit = c(10,10,10,10,rep(NA,8)),
                  baseline = TRUE,
                  do.plot = FALSE,
                  spgdf = FALSE)
    })
    
    # Output Tab Content for Vessel Table -----------------------
    output$tabContent <- renderUI({
      if(input$tabs == "vessel_traffic") {
        DT::DTOutput("editableTable")
      } 
    })
    
    # Reset Values -------- 
    observeEvent(input$reset_input, {
      shinyjs::reset(id = "")
    })
    
    # Shiny Feedback --------  
    iv <- InputValidator$new()
    iv$add_rule("alpha", sv_between(0, 30))
    iv$add_rule("abs_coef", sv_between(0, 5))
    iv$enable()
    
    # Transmission Loss Functions ---------
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
    
    # Scenario Selector -------------
    uploadedDataVal <- reactiveVal()
    
    selectedList <- reactive({
      
      if (input$scenario %in% c("Scenario 1", "Scenario 2", "Scenario 3")) {
        # Return the predefined scenario data
        return(get(paste0("scenario_0", stringr::str_sub(input$scenario, -1))))
      } else {
        
        req(input$custom_vessels, input$custom_piles, input$custom_route)
        file_names <- input$custom_route$name
        file_paths <- input$custom_route$datapath
        file_extensions <- tools::file_ext(file_names)
        paths <- setNames(file_paths, file_extensions)
        print(paths[["shp"]])
        # Validate file types
        validate(
          need(tools::file_ext(input$custom_vessels$name) == "csv", "Please upload a CSV file for vessels."),
          need(tools::file_ext(input$custom_piles$name) == "csv", "Please upload a CSV file for piles."),
          need(tools::file_ext(input$custom_route$name) == "shp", "Please upload a .shp file for the route.")
        )
        
        uploaded_files <- tools::file_ext(input$custom_route$name)
        validate(
          need(all(c("shp", "shx", "dbf") %in% uploaded_files), "Please upload all shapefile components: .shp, .shx, and .dbf files.")
        )
        
        # Read shapefile 
        route_data <- tryCatch({
          Read_Shapefile(input$custom_route)
        }, error = function(e) {
          stop("Failed to read shapefile: ", e$message)
        })
        if (any(sf::st_geometry_type(route_data) %in% c("LINESTRING", "MULTILINESTRING"))) {
          # Convert sf object to SpatialLinesDataFrame
          route_data_sp <- as(route_data, "Spatial")
        } else {
          stop("The geometry type of the route data is not suitable for conversion to SpatialLinesDataFrame.")
        }
        # 
        # Read CSV files: 1) vessel speed table, 2) piling locations
        vessels_data <- read_csv(input$custom_vessels$datapath)
        piles_data <- read_csv(input$custom_piles$datapath)
        piles_data$windfarm <- as.factor(piles_data$windfarm)
        
        # Compile data into a list
        data_list <- list(phase = as.numeric(input$custom_phase),
                          locs = data.table::data.table(piles_data),
                          routes = list(route_data_sp),
                          vessels = data.table::data.table(vessels_data),
                          start.month = c(1, 2, 1),
                          start.day = c(15, 1, 1),
                          piles.per.day = 1,
                          ambient = 80,
                          sourceLvL = 220,
                          lowerdB = 0 ,
                          logrange = 15,
                          absorb = 1.175)
        
        return(data_list)
        
      }
    })
    
    ### Construction & Operation Phase
    currentPhase <- reactive({
      as.numeric(selectedList()[["phase"]])
    })
    
    # Dynamic Sidebar Tab --------------  
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
                            value = selectedList()[['ambient']]),
                sliderInput(
                  "dates",
                  "Piling date",
                  min = range(turbine_df()$date)[1],
                  max = range(turbine_df()$date)[2],
                  step = 1,
                  value = range(turbine_df()$date)[1],
                  animate =
                    animationOptions(interval = 350, loop = FALSE)
                ),
                
                ##
                # turbine_data <- req(turbine_df())
                # 
                # # Calculate min and max dates
                # date_range <- range(turbine_data$date)
                # min_date <- date_range[1]
                # max_date <- date_range[2]
                ##
                
                hr(),
                h4('Change Source Level Values'),
                numericInput("sl", 'Source Level (dB)', min = 10, max = 265, value = selectedList()[['sourceLvL']]),
                hr(),
                h4('Transmission Loss'),
                helpText("We use a simple model, which assumes that transmission ",
                         "loss (TL) depends on log-range (log10R) and  ",
                         "frequency-specific absorption (a).  ",
                         "The model is of the form: TL = b x log10R + a x R. "),
                numericInput("alpha", "Log-range coefficient (b)",
                             min = 10, max = 20, value = selectedList()[['logrange']]),
                numericInput("abs_coef", "Absorption coefficient (a)",
                             min = 0, max = 5, value = selectedList()[['absorb']]),
                hr(),
                h4('Noise mitigation'),
                helpText("Noise abatement systems (e.g., bubble curtains)  ",
                         "may be in place to limit the exposure of ",
                         "protected marine species to pile-driving noise."),
                numericInput("bubble", "Magnitude of noise attenuation (dB)",
                             value = selectedList()[['lowerdB']],
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
    
    ### Turbines  ------------
    turbine_df <- reactive({
      
      if ("locs" %in% names(selectedList()) && nrow(selectedList()[["locs"]]) > 0) {
        data <- selectedList()[["locs"]]
        
        # Check for custom scenario using an input like input$scenario
        if (input$scenario %in% c("scenario_03", "Custom Scenario")) {
          # For custom scenario, just return the locations without additional processing
          return(data)
        } else {
          # Process data for non-custom scenarios
          data <- data %>%
            mutate(buffer = 80) %>%
            sf::st_as_sf(coords = c('longitude', 'latitude')) %>%
            sf::st_set_crs(4326)
          
          return(data)
        }
      } else {
        # Return NULL or an alternative when 'locs' is not suitable
        return(NULL)
      }
      
    })
    
    
    # Update Dates based on Selector --------
    observe({
      turbine_data <- req(turbine_df())
      
      if(!is.null(turbine_data$date)){
        date_range <- range(turbine_data$date)
        min_date <- date_range[1]
        max_date <- date_range[2]
        
        # Update the slider input
        updateSliderInput(session, 'dates', value = min_date,
                          min = min_date, max = max_date)  
      }
      
      
    })
    
    
    # Base Map  ---------------------------------------------------------
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
    
    
    
    # Reset Zoom Observer ------------
    observe({
      input$reset_button
      leafletProxy("map") %>%
        setView(lat = initial_lat, lng = initial_lng, zoom = initial_zoom)
    })
    
    # WEA Zone 1 Observer-----------------
    
    observe({
      
      # define proxy and remove effort group
      proxy <- leafletProxy("map")
      proxy %>% clearGroup("windfarm1")
      
      # # Add Wind Energy Areas
      req(input$windfarm1)
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
        
        # switch to show/hide
        ifelse(input$windfarm1, showGroup(proxy, 'windfarm1'),
               hideGroup(proxy, 'windfarm1'))
        
      }
      
    })
    
    # User Guide Text ------------------
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
    
    # WEA Piling Observer-----------------
    
    observe({
      
      req(turbine_df())
      
      # define proxy and remove effort group
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
    
    # make buffers on the fly - Version 2
    # use addCircles with radius in meters
    observe({
      req(turbine_df(), input$bufferRadius, input$sl, input$bubble, input$alpha, input$abs_coef)
      
      # Clear existing layers upfront to handle toggling off
      leafletProxy("map") %>%
        clearGroup('pile_buffer') %>%
        clearGroup('pile_centroid')
      
      # Proceed only if show_buffer is TRUE
      if (input$show_buffer) {
        # Calculate buffer radius
        buffer_radius <- dB2km(input$bufferRadius, SL = input$sl - input$bubble, logfac = input$alpha, a = input$abs_coef) * 1000
        
        # Applying transformations and optional buffering
        scn1_poly <- turbine_df() %>%
          sf::st_transform(brs_crs) 
        
        if (buffer_radius > 0) {
          scn1_poly <- sf::st_buffer(scn1_poly, dist = buffer_radius)
        }
        
        scn1_poly <- scn1_poly %>%
          sf::st_transform(4326) %>%
          dplyr::filter(date == as.Date(input$dates))  # Ensure date formats are consistent
        
        # Add polygons to the map
        leafletProxy("map") %>%
          addPolygons(
            data = scn1_poly,
            fillColor = "#fdbe85",
            color = "#a63603",
            weight = 1,
            opacity = 1,
            fillOpacity = 0.2,
            group = 'pile_buffer'
          )
      }
      
      # Handle circle markers separately if needed
      if (input$show_buffer) {
        scn1_pt <- turbine_df() %>%
          dplyr::filter(date == as.Date(input$dates))
        
        leafletProxy("map") %>%
          addCircleMarkers(
            data = scn1_pt,
            radius = 0,
            color = "white",
            stroke = ifelse(input$show_buffer, TRUE, FALSE),
            fillOpacity = ifelse(input$show_buffer, 0.5, 0),
            group = 'pile_centroid'
          )
      }
    })
    
    
    # WEA Ship Routes Observer-----------------
    
    observe({
      req(input$wearoute1)
      # define proxy and remove effort group
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
    
    # Marine Mammal Density Observer ---------
    observe({
      
      # define proxy and remove effort group
      proxy <- leafletProxy("map")
      proxy %>%
        clearGroup("mm") %>%
        removeControl("whaleLegend")
      
      if(input$mm_rast){
        
        my_raster <- density_narw[[input$mm_month]]
        r <- raster::raster(my_raster)
        
        r_cut <- format_raster(r)
        # minValue <- min(values(r_cut), na.rm = TRUE)
        # maxValue <- max(values(r_cut), na.rm = TRUE)
        valid_range <- range(raster::values(r_cut), na.rm = TRUE)
        pal <- colorNumeric(palette = pals::viridis(n = 10), 
                            domain = valid_range,
                            na.color = NA)
        
        proxy %>%
          addRasterImage(r_cut,
                         colors = pal,
                         opacity = 1,
                         group = 'mm')%>%
          addLegend(pal = pal,
                    values = raster::values(r_cut),  
                    na.label = 'NA',
                    title = "Whale Density",
                    opacity = 1,
                    position = "bottomright",
                    group = 'mm',
                    layerId = "whaleLegend")  
        
        # switch to show/hide
        ifelse(input$mm_rast, showGroup(proxy, 'mm_rast'),
               hideGroup(proxy, 'mm_rast'))
        
      } # end marine mammal
      
    })
    
    # Initialize currentTable as a reactiveVal
    currentTable <- reactiveVal()
    
    # On start-up or when selectedList() changes, update currentTable
    observe({
      currentTable(selectedList()[["vessels"]])
    })
    
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
    
    # # TODO: Run/update the vessel maps.
    values <- reactiveValues(currentScenario = NULL)
    
    # Then your observeEvent might look like this:
    observeEvent(input$saveBtn, {
      # Update the reactive value to the new scenario
      values$currentScenario <- map_vessels(
        obj = updatedScenario(),
        z = "dist",
        strike_scalar = 1e-7,
        which.month = 1:12,
        vessel.speed = NULL,
        speed.limit = c(10, 10, 10, 10, rep(NA, 8)),
        baseline = TRUE,
        do.plot = FALSE,
        spgdf = FALSE
      )
      
    })
    
    # Adjust vesselmaps to depend on values$currentScenario:
    vesselmaps <- reactive({
      req(values$currentScenario)
      values$currentScenario
    })
    
    # Vessel Traffic Density Observer ---------
    observe({
      # Define the map proxy
      proxy <- leafletProxy("map")
      
      # Clear the existing legend
      proxy %>% clearGroup("vessel")
      proxy %>% removeControl("vesselLegend")
      
      # Check conditions for adding the raster image
      if (input$vessel_rast_yesno && !is.null(input$vessel_scenario) && length(input$vessel_scenario) > 0) {
        sel_month <- which(names(traffic) == input$month)
        
        # Select the raster based on the scenario
        if (input$vessel_scenario == 'baseline') {
          my_raster <- vesselmaps()[[sel_month]][[2]]
        } else if (input$vessel_scenario == 'scenario') {
          my_raster <- vesselmaps()[[sel_month]][[1]]
        }
        
        # Set parameters for the raster image
        r <- my_raster
        prange <- c(0, 1000)  # Adjust this range based on your data
        valid_range <- range(raster::values(r), na.rm = TRUE)
        pal <- colorNumeric(palette = rev(pals::magma(n = 10)), 
                            domain = valid_range,
                            na.color = NA)
        
        # Add raster image and update legend
        proxy %>%
          addRasterImage(r,
                         colors = pal,
                         opacity = 1,
                         group = 'vessel') %>%
          addLegend(pal = pal,
                    values = raster::values(r),  
                    na.label = 'NA',
                    title = "Vessel Strike Risk",
                    opacity = 1,
                    position = "bottomright",
                    group = 'vessel',
                    layerId = "vesselLegend")  
        
        # Show the group containing the raster
        showGroup(proxy, 'vessel')
        
      } else {
        # Hide the vessel group if conditions are not met
        hideGroup(proxy, 'vessel')
      }
    })
    
    
    # tss observer ------------------------------------------------------  
    
    observe(priority = 4, {
      
      # define proxy
      proxy <- leafletProxy("map")
      proxy %>% clearGroup('tss')
      
      if(input$tss){
        
        # plot shipping lanes
        
        proxy %>%
          addPolygons(data = tss,
                      weight = .5,
                      color = 'grey',
                      fillColor = 'grey',
                      options = pathOptions(clickable = F),
                      group = 'tss')
        
        # switch to show/hide
        ifelse(input$tss, showGroup(proxy, 'tss'), hideGroup(proxy, 'tss'))
      }
      
    })
    
    # area to be avoided observer ------------------------------------------------------  
    
    observe(priority = 4, {
      
      # define proxy
      proxy <- leafletProxy("map")
      proxy %>% clearGroup('atba')
      
      if(input$atba){
        
        # plot Areas to Be Avoided
        
        proxy %>%
          addPolygons(data = atba,
                      weight = .5,
                      color = '#2ca25f',
                      fillColor = '#2ca25f',
                      options = pathOptions(clickable = F),
                      group = 'atba')
        
        # switch to show/hide
        ifelse(input$atba, showGroup(proxy, 'atba'), hideGroup(proxy, 'atba'))
      }
      
    })
    
    # NARW Speed Zones observer ------------------------------------------------------  
    
    observe(priority = 4, {
      
      # define proxy
      proxy <- leafletProxy("map")
      proxy %>% clearGroup('speed_zones')
      
      if(input$speed_zones){
        
        # plot Areas to Be Avoided
        
        proxy %>%
          addPolygons(data = speed_zones,
                      weight = .5,
                      color = '#a63603',
                      fillColor = '#a63603',
                      options = pathOptions(clickable = F),
                      group = 'speed_zones')
        
        # switch to show/hide
        ifelse(input$speed_zones, showGroup(proxy, 'speed_zones'), hideGroup(proxy, 'speed_zones'))
      }
      
    })
    
    # Uncheck Map Features
    observeEvent(input$scenario, {
      
      if (!is.null(input$scenario)) {
        if (input$scenario != "Custom Scenario") {
          updateCheckboxInput(session, "windfarm1", value = FALSE)
          updateCheckboxInput(session, "wearoute1", value = FALSE)
        }
        updateCheckboxInput(session, "windfarm2", value = FALSE)
        updateCheckboxInput(session, "windfarm3", value = FALSE)
        updateCheckboxInput(session, "wearoute2", value = FALSE)
        updateCheckboxInput(session, "wearoute3", value = FALSE)
        updateCheckboxInput(session, "piles", value = FALSE)
        updateCheckboxInput(session, "tss", value = FALSE)
        updateCheckboxInput(session, "atba", value = FALSE)
        updateCheckboxInput(session, "speed_zones", value = FALSE)
      }
      
    })
    
    # Assemble Parameters -------
    my_params <- reactive({
      req(currentTable(), input$bufferRadius, input$sl, input$alpha, 
          input$abs_coef, input$bubble, input$custom_phase)
      
      list(vessels = currentTable(),
           ambient = input$bufferRadius,
           sourceLvL = input$sl,
           logrange = input$alpha,
           absorb = input$abs_coef,
           lowerdB = input$bubble,
           phase = as.numeric(input$custom_phase))
      
    })
    
    # # Update Scenario List Object ----- 
    # updatedScenario <- reactive({
    #   # Force evaluation of input$saveBtn to maintain original trigger behavior
    #   input$saveBtn
    #   
    #   # Return the updated list based on the most recent table data
    #   new_scenario <- update.narwscenario(selectedList(), currentTable())
    #   class(new_scenario) <- c('narwscenario', 'list')
    #   return(new_scenario)
    #   
    # })
    
    updatedScenario <- reactive({
      # Log for debugging
      cat("Button Pressed at:", Sys.time(), "\n")
      
      # Force evaluation of input$saveBtn to maintain original trigger behavior
      input$saveBtn
      
      # Return the updated list based on the most recent table data
      new_scenario <- update.narwscenario(selectedList(), currentTable())
      
      # Debug output
      # if (is.list(new_scenario)) {
      #   cat("new_scenario is a list.\n")
      # } else {
      #   cat("new_scenario is not a list. It is a", class(new_scenario), "\n")
      # }
      
      # Set class
      class(new_scenario) <- c('narwscenario', 'list')
      # cat("new_scenario is a", class(new_scenario), "\n")
      
      # Return with the new class
      return(new_scenario)
    })
    
    updatedList <- reactive({
      # print("updatedList executed")
      req(updatedScenario())
      req(my_params())
      modifyList(updatedScenario(), my_params())
      
    })
    
    
    output$download <- downloadHandler(
      filename = function() {
        req(input$scenario)
        paste0("BOEM-NARWIND_parameters_", gsub("[[:blank:]]", "", input$scenario), "_", 
               format(Sys.time(), "%Y%m%d_%H%M%S"), ".RData")
      },
      content = function(file) {
        req(updatedList())
        scenario_custom <- updatedList()
        # Save the object properly
        save(scenario_custom, file = file)
      }
    )
    
    user_shp <- reactive({
      req(input$custom_shp)
      Read_Shapefile(input$custom_shp) %>%
        sf::st_transform(crs = 4326)
    })
    
    # Store a non-reactive copy of the data, updated every time the reactiveValues change
    latestData <- reactiveVal()
    observe({
      req(updatedList())
      latestData(updatedList())
      cat("Non-Reactive latestData Object captured at: ", Sys.time(), "\n") 
    })
    
    # Exit Handling
    session$onSessionEnded(function() {
      print("Session is ending; saving parameters...")
      
      # Completely isolate the retrieval of latestData
      latest_data <- isolate(latestData())
      
      # Check if latest_data is non-null
      if (!is.null(latest_data)) {
        # cat("Structure of latestData:\n")
        # str(latest_data)
        saveParameters(latest_data)
      }
      
      stopApp()
    })
    
    onStop(function() {
      cat("App is stopping now.\n")
    })
    
  }
  
  # Launch Shiny app----
  shinyApp(ui = ui,
           server = server,
           options = list(launch.browser = TRUE))
  
} # End function
