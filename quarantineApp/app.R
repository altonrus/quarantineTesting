#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyjs)
library(ggplot2)
library(data.table)
source("quarantine-functions.R")



##
#TO DO
#
# Debug run_sim, numbers aren't correct.
# Get dt_analysis to update after run_sim is run
#
###



#Incubation time from Lauer 2020 bootstrapped posteriors
incubation_dist_fit_lnorm <- readRDS("ncov_inc_fit_boot.rds")
dt_incubation_dists_lnorm <- data.table(incubation_dist_fit_lnorm@samples)


# Define UI
ui <- tagList(
    useShinyjs(),
    #Navbar with 3 tabs: "Main", "Sim Inputs", "About"
    navbarPage("COVID traveler quarantine policy assessment tool",
               
                 tabPanel("Main", ##########
                          #Max-width to keep it from spreading out too much in wide browzers
                          tags$head(tags$style(type="text/css", ".container-fluid {  max-width: 800px}")),
                          #Controls on top
                          fluidRow(
                              column(
                                  width = 12,
                                  h1("Analysis parameters"),
                                  "Additional parameters can be customized in the *Sim inputs* tab.",
                                  br(),
                                  br(),
                                  fluidRow(
                                      column(
                                          width = 6,
                                          radioButtons("metric", "Metric:",
                                                       c("Days at-risk per infected traveler" = "days",
                                                         "Person-days per 10,000 travelers" = "person-days",
                                                         "Secondary cases" = "sec-cases")),
                                          checkboxInput("include_test",
                                                        "Include testing?",
                                                        value = TRUE)
                                      ),
                                      column(
                                          width = 6,
                                          numericInput("prev",
                                                       "Prevalence of active infection",
                                                       value = 0.001,
                                                       min = 0,
                                                       max = 1),
                                          numericInput("sec_cases_per_day",
                                                       "Secondary infections per person-day infectious in community",
                                                       min = 0,
                                                       value = .2)
                                      )
                                  ),
                                  h1("Results")
                              )
                          ),
                          #Results below
                          fluidRow(
                              column(width = 12, align = "center",
                                     #Plot
                                     h2("Plot"),
                                     plotOutput("sim_plot", height = "400px", width = "550px"),
                                     #table
                                     h2("Table"),
                                     tableOutput("sim_data")
                                     )
                              )
                          ),
                 tabPanel("Sim inputs", ########
                          #Instructions on top
                          fluidRow(
                              column(12,
                                     align = "center",
                                     h4("Customize simulation parameters and then press \"re-run simulation\" below."),
                                     br()
                                     )
                              ),
                          #Controls in 2 cols
                          fluidRow(
                              column(6,
                                     h4("Options"),
                                     radioButtons("n_iters", "Number of iterations",
                                                  c("1000 (slow but full accuracy)" = 1000,
                                                    "100 (faster update)" = 100)),
                                     selectizeInput(
                                         "dur_quarantine", 
                                         "Quarantine durations to compare (days)", 
                                         choices = 0:20,
                                         multiple = TRUE,
                                         options = list(create = TRUE),
                                         selected = c(0, 2, 5, 7, 14)
                                     ),
                                     numericInput("RNseed",
                                                  "Seed for random number generator",
                                                  min = 0,
                                                  value = 91,
                                                  step = 1),
                                     br(),
                                     h4("Test sensitivity"),
                                     sliderInput("sn_presympt",
                                                 "Sensitivity, pre-symptomatic phase",
                                                 min = 0,
                                                 max = 1,
                                                 value = .70),
                                     sliderInput("sn_sympt",
                                                 "Sensitivity, Symptomatic phase",
                                                 min = 0,
                                                 max = 1,
                                                 value = .70),
                                     sliderInput("sn_asympt",
                                                 "Sensitivity, asymptomatic phase",
                                                 min = 0,
                                                 max = 1,
                                                 value = .60)
                              ),
                              column(6,
                                     h4("Probabilities"),
                                     sliderInput("prob_asympt",
                                                 "Asymptomatic infection",
                                                 min = 0,
                                                 max = 1,
                                                 value = .24),
                                     sliderInput("prob_quarantine_compliance",
                                                 "Quarantine compliance",
                                                 min = 0,
                                                 max = 1,
                                                 value = .80),
                                     sliderInput("prob_isolate_sympt",
                                                 "Isolate when symptomatic",
                                                 min = 0,
                                                 max = 1,
                                                 value = .80),
                                     sliderInput("prob_isolate_test",
                                                 "Isolate after positive test",
                                                 min = 0,
                                                 max = 1,
                                                 value = .90),
                                     sliderInput("prob_isolate_both",
                                                 "Isolate with symptoms + positive test",
                                                 min = 0,
                                                 max = 1,
                                                 value = 1.0),
                                     radioButtons("rand_u",
                                                  "How long are travelers infected before arriving?",
                                                  c("Infection has progressed a random percent 0% - 100%" = TRUE,
                                                  "Infected immediately before arriving" = FALSE))
                              )
                              
                          ),
                          #Button to run sim below
                          fluidRow(
                              column(12,
                                     br(),
                                     align = "center",
                                     actionButton("run_sim",
                                                  "Re-run simulation",
                                                  width = '30%',
                                                  style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                                     br(),
                                     br()
                              )
                          )
                          ),
                 tabPanel("About", #########
                          h1("Tool summary"),
                          "This tool is intended to estimate the impact of policies to manage the risk of 
                          COVID-19 introduced by arriving travelers. The policies compared are mandatory 
                          quarantine of varying length with or without testing. The user can customize 
                          several simulation parameters related to compliance with quarantine and isolation, test 
                          sensitivity, the relative prevalence of asymptomatic vs. infection. Uncertainty in outcomes 
                          reflect uncertainty in the distribution of the duration of the 
                          pre-symptomatic-infectious and symptomatic-infectious periods for those with symptomatic 
                          infections, the duration of the asymptomatic-infectious period for those with asymptomatic 
                          infections, and the incubation period for both infection types",
                          br(),
                          br(),
                          "Users can toggle between three metrics and choose whether to display policies with testing 
                          directly on the main tab without re-running the simulation. The metric \"Days at-risk per 
                          infected traveler\" refers to the average time an infected traveler would be at risk of 
                          infecting members of the community because they are infectious and not in quarantine or 
                          isolation. This metric requires only the parameters from the simulation. The outcome 
                          \"Person-days at risk per 10,000 travelers\" has both infected and uninfected travelers 
                          in the denominator. For this outcome, the user must indicate an estimated or assumed 
                          prevalence of active infection among travelers. The outcome \"Secondary cases\" estimates how
                          many community members you would expect to be infected by travelers. Calculating this outcome 
                          requires an estimate of the rate of econdary infections per person-day that an traveler is infectious and 
                          at-risk in the community.",
                          br(),
                          h1("About the authors"),
                          "This tool was created by Alton Russell, PhD candidate at Stanford University and visitor
                          at  McGill Clinical & Health Informatics (MCHI), in collaboration
                          with David Buckeridge, Professor of Epidemiology and Biostatistics at McGill University
                          and director of the surveillance lab at MCHI.",
                          br(),
                          h2("Links"),
                          tags$a(href="http://mchi.mcgill.ca/", "McGill Clinical Health and Informatics (MCHI)"),
                          br(),
                          tags$a(href="http://altonrussell.com/", "Alton's website"),
                          ) #########

               ) ##Close navbar
) ##Close taglist

Values <- reactiveValues(dt_raw = fread("dt_raw.csv"))

# Define server logic required to draw a histogram
server <- function(input, output, session) {

    

    
    #dt_analysis is conductor between input parameters and output (plot and displayed table)
    dt_analysis <- reactive({
        make_dt_analysis(
                Values$dt_raw,
                input$metric,
                input$include_test,
                input$prev,
                input$sec_cases_per_day)
    })
    
    #PLOT
    output$sim_plot <- renderPlot(
        # dt_analysis <- make_dt_analysis(
        #     Values$dt_raw,
        #     input$metric,
        #     input$include_test,
        #     input$prev,
        #     input$sec_cases_per_day),
        if (input$include_test ){
            ggplot(dt_analysis())+
                geom_pointrange(aes(x = testing, y = q0.5, ymin = q0.01, ymax = q0.99, 
                                    color = factor(quarantine_length), group = factor(quarantine_length)),
                                position = position_dodge(width = 0.2), size=1)+
                xlab("")+
                ylab( switch(input$metric,
                             "days" = "Days at risk per infected traveler",
                             "person-days" = "Person-days at risk\nper 10,000 travelers",
                             "sec-cases" = "Secondary cases\nper 10,000 travelers"))+
                scale_color_discrete(name="Quarantine length (days)")+
                scale_x_discrete(labels = c("No test", "Test"))+
                geom_line(aes(x = testing, y = q0.5, group = factor(quarantine_length), 
                              color = factor(quarantine_length)), linetype="dashed",
                          position = position_dodge(width = 0.2))+
                theme_minimal()+
                theme(legend.position = "bottom",
                      text = element_text(size = 18))+
                scale_y_continuous(limits=c(0,NA))
        } else{
            ggplot(dt_analysis())+
                geom_pointrange(aes(x = factor(quarantine_length), y = q0.5, ymin = q0.01, ymax = q0.99, 
                                    color = factor(quarantine_length)),
                                size=1)+
                xlab("Quarantine length (days)")+
                ylab( switch(input$metric,
                             "days" = "Days at risk per infected traveler",
                             "person-days" = "Person-days at risk\nper 10,000 travelers",
                             "sec-cases" = "Secondary cases\nper 10,000 travelers"))+
                scale_color_discrete(name="Quarantine length (days)")+
                theme_minimal()+
                theme(legend.position = "none",
                      text = element_text(size = 18))+
                scale_y_continuous(limits=c(0,NA))
        }
    
    )
    
    #TABLE
    output$sim_data <- renderTable(
        dt_analysis()
        )
    

    
    # observe({ ## Hide sensitivity options if testing not selected
    #     if (input$include_test == FALSE) {
    #         shinyjs::hide("sn_sympt", anim = TRUE)
    #         shinyjs::hide("sn_asympt", anim = TRUE)
    #         shinyjs::hide("sn_presympt", anim = TRUE)
    #     } else {
    #         shinyjs::show("sn_sympt", anim = TRUE)
    #         shinyjs::show("sn_asympt", anim = TRUE)
    #         shinyjs::show("sn_presympt", anim = TRUE)
    #     }
    # })
    observe({ ## Hide the prev and secondary cases if outcome doesn't require them
        if (input$metric == "days") {
            shinyjs::disable("prev")
            shinyjs::disable("sec_cases_per_day")
        } else if (input$metric == "person-days"){
            shinyjs::enable("prev")
            shinyjs::disable("sec_cases_per_day")
        } else { #"Secondary cases
            shinyjs::enable("prev")
            shinyjs::enable("sec_cases_per_day")
        }
    })
    
    #SiMULATION IS RUN
    observeEvent(input$run_sim, {
        #Collect parameters
        sim_params <- list(
            prob_asympt = input$prob_asympt,
            prob_isolate_test = input$prob_isolate_test,
            prob_isolate_sympt = input$prob_isolate_sympt,
            prob_isolate_both = 1.0,
            sn_presympt = input$sn_presympt,
            sn_sympt = input$sn_sympt,
            sn_asympt = input$sn_asympt,
            prob_quarantine_compliance = input$prob_quarantine_compliance,
            dur_presympt_mean_lb = 1.8,
            dur_presympt_mean_ub = 2.8,
            dur_presympt_var_lb = 4.0,
            dur_presympt_var_ub = 6.0,
            dur_sympt_mean_lb = 2.6,
            dur_sympt_mean_ub = 3.9,
            dur_sympt_var_lb = 3.0,
            dur_sympt_var_ub = 4.5,
            dur_asympt_mean_lb = 4.0,
            dur_asympt_mean_ub = 6.0,
            dur_asympt_var_lb = 4.0,
            dur_asympt_var_ub = 6.0,
            n_iters = as.numeric(input$n_iters),
            dur_quarantine = as.numeric(input$dur_quarantine),
            seed = input$RNseed,
            rand_u = input$rand_u
        )
        print(sim_params)
        #Run simulation
        Values$dt_raw <<- run_sim(sim_params)
        })
    
    
}

# Run the application 
shinyApp(ui = ui, server = server)
