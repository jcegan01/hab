#
# You can run the application by clicking
# the 'Run App' button above.
#
# Created by Alex Ligo, 2021-05-28
#

library(dplyr)
library(tidyr)
library(ggplot2)
library(shiny)
library(plotly)

# load data pre calculated in R notebook 01_wrangling
load( 'data/wqdata_cleaned.rda' )

# load trophic states
load( 'data/reservoir_summary.rda' )

# list of reservoirs
reservoirs <- unique(dat4$reservoir)
# list of parameters
params <- list( "Algae Chl-A indicator" = "chlA"
              , "Ammonia" = "ammonia"
              , "Dissolved Oxygen (mg/l)" = "DO_mgl"
              , "Dissolved Oxygen (% saturation)" = "DO_perc"
              , "Nitrate" = "nitrate"
              , "Outflow/Discharge" = "flow"
              , "Phosphate" = "phosphate"
              , "Secchi Depth" = "secchi_in" 
              , "Total Dissolved Solids" = "TDS"
              , "Total Nitrogen" = "N_total"
              , "Total Organic Carbon" = "C_total"
              , "Total Phosphorous" = "P_total"
              , "Total Suspended Solids/Fixed Suspended Solids" = "TSS"
              , "Water temperature" = "temp_c"
            )
data <- dat4 
par_names <- as.list(unlist(names(params)))
names(par_names) <- unlist(params,use.names=F)
data$par_name <- recode( data$param, !!!par_names )

# number of data points
data_cnt <- data %>% 
    group_by(reservoir,param,par_name) %>%
    summarise(cnt = n(), .groups = "drop_last") %>%
    mutate( label = paste0(par_name,' (n=',cnt,')')
            , across(where(is.factor), as.character) ) 
data_cnt <- data_cnt[with(data_cnt, order(reservoir,par_name)),]

# time range in data
yrrange <- as.numeric(range(dat4$yr))

# UI

ui <- fluidPage(
    # Application title
    titlePanel('HAB parameter visualization'),
    sidebarLayout(
        # Grey Sidebar with input controls
        sidebarPanel(
            selectInput( "reservoir", 'Select reservoir:', reservoirs, selected = 'TAR' )
            , fluidRow( checkboxGroupInput('params', 'Select parameters to display:', params, selected=params[c(1,3,14)]) )
            , fluidRow( radioButtons( 'avg', 'Averaging:',
                                      choices = list('No averaging' = 'raw'
                                                     , 'Average over locations in each depth of reservoir' = 'avg_loc'),
                                      selected = 'avg_loc' ) )
        ), # sidebarPanel
        # date slider and graph for selected reservoir
        mainPanel(
                wellPanel( fluidRow(
                    column(8, sliderInput('years', 'Interval to display:', min=yrrange[1], max=yrrange[2]
                                                , value=c(yrrange[2]-5, yrrange[2]), sep=''))
                    , column(4, sliderInput('depths', 'Select depths:', min=0, max=100, value=c(0, 0))))
                    , fluidRow(
                        tabsetPanel(
                            tabPanel('Time series',
                                 plotlyOutput("distPlot", height = "550px")
                        ), # tabPanel
                            tabPanel('Correlations',
                                plotlyOutput("corPlot", height = "550px")
                            ) # tabPanel('Correlations'
                        ) # tabsetPanel()
                    )
                    ) # wellPanel
        ) # mainPanel
    ) # sidebarLayout
) # fluidPage()

# Define server logic required to draw chart based on input
server <- function(input, output, session) {
    updatedInput <- reactive({
        # This reactive function is executed every time the reservoir or parameter selection change
        
        # the goal is to update the display of number of data points available for each param.
        cnts <- data_cnt %>% filter( reservoir==input$reservoir )
        updateCheckboxGroupInput(session, 'params'
                                  , choiceNames = cnts$label, choiceValues = cnts$param
                                 , selected = input$params)
        
        # the goal is to update the date range in the slider
        df <- data %>% filter( reservoir==input$reservoir & param %in% input$params )
        yrrange <- as.numeric(range(df$yr))
        
        drange <- range(df$depth)

        if ( !anyNA(yrrange) ){
            updateSliderInput(session, 'years', min=yrrange[1], max=yrrange[2]
                              , value=c(yrrange[2]-5, yrrange[2]) )
            updateSliderInput(session, 'depths', min=drange[1], max=drange[2]
                              , value=c(drange[1], drange[1]) )
        }
        
        # reservoir's trophic state to be displayed in chart title 
        trophic <- reservoir.summary3$trophic_state[reservoir.summary3$reservoir == input$reservoir]
        if ( !is.na(trophic) )
            ggtit <- paste0(input$reservoir, ': ', trophic, ' (based on latest 5-year data)')
        else
            ggtit <- paste0(input$reservoir, ': latest 5-year data insufficient to determine trophic state')
        
        return( list(df, ggtit) )
    })
    
    output$distPlot <- renderPlotly({
        # time series plots
        u <- updatedInput()
        df <- u[[1]]
        ggtit <- u[[2]]

        dt <- df %>% filter( yr >= input$years[1] & yr <= input$years[2] 
                               & depth >= input$depths[1] & depth <= input$depths[2] )

        if ( input$avg == 'avg_loc')
            dt <- dt %>% group_by( param, par_name, sample_date, depth ) %>%
                     summarise( value = mean(value), loc_id = 'all', units = first(units)
                                , .groups = "drop_last") %>%
                     ungroup()
        
        p <- ggplot( dt ) +
            geom_point(aes(units=units, color=par_name, loc_id=loc_id, depth=depth, x=sample_date, y=value), size=1.5 ) +
            ggtitle(ggtit) +
            ylab('Parameter Value') +
            theme_bw() +
            theme(plot.title = element_text(hjust = 0.5), 
                  text = element_text(size=12), axis.title.x = element_blank()
                  , legend.title=element_blank()) 

        ggplotly(p, tooltip = c("x","y","units","par_name","loc_id","depth")
                 , dynamicTicks = TRUE) %>%
            layout(legend = list(orientation = "h", y = -0.05))
    }) # closing of distPlot()

    output$corPlot <- renderPlotly({
        # Correlation plot
        u <- updatedInput()
        df <- u[[1]]
        ggtit <- u[[2]]
        
        dt <- df %>% filter( yr >= input$years[1] & yr <= input$years[2] 
                             & depth >= input$depths[1] & depth <= input$depths[2] )
        dp <- dt %>% 
                filter( param %in% input$params[1:2] )  %>% # consider only first two selections for correlations
                group_by( sample_date, depth, param ) %>%
                summarise( value = mean(value), .groups = "drop_last" ) %>%
                ungroup()
        
        dp <- dp %>%
                    pivot_wider( names_from = param, values_from = value ) 
        dp <- dp[complete.cases(dp),]
        
        if ( nrow(dp) & input$params[1] %in% names(dp) & input$params[2] %in% names(dp) ){
            # correlation
            corr = sprintf( "%.3f", cor(dp[,input$params[1]],dp[,input$params[2]]) )
            
            tit <- paste0('Correlation between ', input$params[1], ' and\n', y=input$params[2]
                         , ' = ', corr, ' (n = ', nrow(dp), ')')
            p <- ggplot( dp, aes_string(depth='depth', x=input$params[2], y=input$params[1]) ) +
                geom_point( shape=1 ) +
                geom_smooth(method=lm) +  # Add linear regression line, add shaded confidence region
                theme_bw() +
                ggtitle(ggtit) +
                theme(plot.title = element_text(hjust = 0.5), text = element_text(size=12)) 
            xr <- ggplot_build(p)$layout$panel_scales_x[[1]]$range$range
            yr <- ggplot_build(p)$layout$panel_scales_y[[1]]$range$range
            titx <- xr[1] + (xr[2] - xr[1]) * 0.3
            tity <- yr[1] # + (yr[2] - yr[1]) * 0.1
            p <- p + geom_text( x=titx, y=tity, label=tit, color='blue' )
            
            ggplotly(p, tooltip = c("date","x","y","depth"), dynamicTicks = TRUE) 
        }
    }) # closing of corPlot()
}

# Run the application 
shinyApp(ui = ui, server = server)
