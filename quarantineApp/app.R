#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# library(shiny)
# library(shinyjs)
# library(ggplot2)
# library(data.table)
# library(trunc)
library(pacman)
p_load(shiny)
p_load(shinyjs)
p_load(ggplot2)
p_load(data.table)
p_load(truncdist)
p_load(scales)
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
                                                       c("Percent risk reduced" = "risk",
                                                         "Days at-risk per infected traveler" = "days",
                                                         "Adjusted days at-risk per infected traveler"="adj-days",
                                                         "Person-days per 10,000 travelers" = "person-days",
                                                         "Secondary cases" = "sec-cases"))
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
                                     plotOutput("sim_plot", height = "400px", width = "700px"),
                                     br(),
                                     actionButton("reset",
                                                  "Reset to base scenario",
                                                  width = '30%',
                                                  style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                                     br(),
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
                                     h5("Quarantine durations (days)"),
                                     selectizeInput(
                                         "dur_quarantine_alone", 
                                         "With no test", 
                                         choices = 0:20,
                                         multiple = TRUE,
                                         options = list(create = TRUE),
                                         selected = c(0, 2, 5, 7, 14)
                                     ),
                                     selectizeInput(
                                         "dur_quarantine_endtest", 
                                         "With test 24h before quarantine end", 
                                         choices = 1:20,
                                         multiple = TRUE,
                                         options = list(create = TRUE),
                                         selected = c(0, 2, 5, 7, 14)
                                     ),
                                     selectizeInput(
                                         "dur_quarantine_arrivetest", 
                                         "With test on arrival", 
                                         choices = 0:20,
                                         multiple = TRUE,
                                         options = list(create = TRUE),
                                         selected = c(0, 2, 5, 7, 14)
                                     ),
                                     selectizeInput(
                                         "dur_quarantine_pretest", 
                                         "With test 72h before arrival", 
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
                                     radioButtons("infection_timing",
                                                  "How long are travelers infected before arriving?",
                                                  c("Infection progression random but travelers with symptomatic infections are no more than 24h into the symptomatic phase"="rand_nosympt_24hr_before_arrival",
                                                    "Infection progresssion random and does not include symptomatic phase (i.e., no one is traveling with symptoms)" = "rand_presympt",
                                                    "Infection progression random including symptomatic phase (i.e., people are traveling with symptoms)" = "rand_incl_sympt",
                                                    "Infected immediately before arriving" = "inf_upon_arrival"))
                              ),
                              column(6,
                                     sliderInput("rr_asympt",
                                                 "Relative transmission risk of asymptomatic infection",
                                                 min = 0.2,
                                                 max = 1.2,
                                                 value = 0.49),
                                     h4("Probabilities"),
                                     sliderInput("prob_asympt",
                                                 "Asymptomatic infection",
                                                 min = 0,
                                                 max = 1,
                                                 value = .40),
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
                          "Users can toggle between five metrics. The metric \"Days at-risk per 
                          infected traveler\" refers to the average time an infected traveler would be at risk of 
                          infecting members of the community because they are infectious and not in quarantine or 
                          isolation. The metric \"Adjusted days at-risk per infected travelers\" multiplies the 
                          infectious days in the community for those with asymptomatic infections by the relative 
                          transmission risk for asymptomatic infections, which is less than 1 in the base case. 
                          This reflects the belief that symptomatic infections have a higher community transmission risk. 
                          The primary metric \"Percent risk reduced\" is calculated 
                          by comparing the adjusted days at-risk for each policy to that of a 0-day quarantine
                          without testing policy.
                          While these first 3 metrics require only the parameters from the simulation, two others 
                          require additional parameters. The outcome \"Person-days at risk per 10,000 travelers\" 
                          has both infected and uninfected travelers 
                          in the denominator. For this outcome, the user must indicate an estimated or assumed 
                          prevalence of active infection among travelers. The outcome \"Secondary cases\" estimates how
                          many community members you would expect to be infected by travelers. Calculating this outcome 
                          requires an estimate of the rate of secondary infections per person-day that an traveler is infectious and 
                          at-risk in the community.",
                          br(),
                          
                          h1("About the authors"),
                          "This tool was created by Alton Russell, postdoctoral research fellow at the Mass General
                          Hospital Institute for Technology Assessment and Harvard Medical School and visitor
                          at  McGill Clinical & Health Informatics (MCHI), in collaboration
                          with David Buckeridge, Professor of Epidemiology and Biostatistics at McGill University
                          and director of the surveillance lab at MCHI.",
                          br(),
                          br(),
                          "For questions, feel free to email altonr <at> stanford <dot> edu.",
                          br(),
                          br(),
                          "A huge thanks to Aman Verma and Maxime Lavigne for technical assistance with hosting the app.",
                          br(),
                          h2("Links"),
                          tags$a(href="http://mchi.mcgill.ca/", "McGill Clinical Health and Informatics (MCHI)"),
                          br(),
                          tags$a(href="http://altonrussell.com/", "Alton's website"),
                          br(),
                          tags$a(href="https://github.com/altonrus/quarantineTesting/", "Github with code for this project")
                          ) #########

               ) ##Close navbar
) ##Close taglist




server <- function(input, output, session) {

    Values <- reactiveValues(dt_raw = fread("dt_raw.csv"))
    #dt_analysis is conductor between input parameters and output (plot and displayed table)
    dt_analysis <- reactive({
        make_dt_analysis(
                Values$dt_raw,
                input$metric,
                input$prev,
                input$sec_cases_per_day)
    })
    
    #MAKE PLOT
    output$sim_plot <- renderPlot(
        ggplot(dt_analysis())+
            geom_pointrange(aes(x = factor(quarantine_length), 
                                y = q0.5,  ymin = q0.01, ymax = q0.99, 
                                color= testing, shape=testing), 
                            position = position_dodge(width = 0.5))+
            switch(input$metric,
                   "risk"=scale_y_continuous(labels = label_percent(accuracy=1), breaks = seq(0,1,0.1), minor_breaks = seq(0.05, 0.95, 0.1),
                                             limits=c(0,1), name = "% community transmission risk reduced"),
                   "days" = scale_y_continuous(limits=c(0,NA), name = "Days at risk per infected traveler"),
                   "adj-days" = scale_y_continuous(limits=c(0,NA), name = "Adjusted days at risk per infected traveler"),
                   "person-days" = scale_y_continuous(limits=c(0,NA), name = "Person-days at risk per 10,000 travelers"),
                   "sec-cases" = scale_y_continuous(limits=c(0,NA), name = "Secondary cases per 10,000 travelers")
            )+
            xlab("Quarantine length (days)")+
            scale_color_discrete(name = "Testing scenario", labels = test_scenario_labs)+scale_shape(name = "Testing scenario", labels = test_scenario_labs)+
            theme_bw()+
            theme(legend.position = "right",
                  panel.grid.minor = element_line(linetype="dotted"),
                  text = element_text(size=16))#+
            #scale_y_continuous(limits=ifelse(input$metric=="risk", c(0,1), c(0,NA)))
    )
        
        
        # ggplot(dt_analysis())+
        #         geom_pointrange(aes(x = testing, y = q0.5, ymin = q0.01, ymax = q0.99, 
        #                             color = factor(quarantine_length), group = factor(quarantine_length)),
        #                         position = position_dodge(width = 0.2), size=1)+
        #         xlab("")+
        #         ylab( switch(input$metric,
        #                      "days" = "Days at risk per infected traveler",
        #                      "person-days" = "Person-days at risk\nper 10,000 travelers",
        #                      "sec-cases" = "Secondary cases\nper 10,000 travelers",
        #                      "risk"=""))+
        #         scale_color_discrete(name="Quarantine length (days)")+
        #         scale_x_discrete(labels = c("No test", "Test"))+
        #         geom_line(aes(x = testing, y = q0.5, group = factor(quarantine_length), 
        #                       color = factor(quarantine_length)), linetype="dashed",
        #                   position = position_dodge(width = 0.2))+
        #         theme_minimal()+
        #         theme(legend.position = "bottom",
        #               text = element_text(size = 18))+
        #         scale_y_continuous(limits=c(0,NA))
    
    
    #TABLE
    output$sim_data <- renderTable(
        dt_analysis()
        )
    

    
    observe({ ## Hide the prev and secondary cases if outcome doesn't require them
        if (input$metric %in% c("days", "adj-days", "risk")) {
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
            prob_isolate_both = input$prob_isolate_both,
            sn_presympt = input$sn_presympt,
            sn_sympt = input$sn_sympt,
            sn_asympt = input$sn_asympt,
            prob_quarantine_compliance = input$prob_quarantine_compliance,
            dur_presympt_mean_lb = 1.8,
            dur_presympt_mean_ub = 2.8,
            dur_presympt_var_lb = 4.0,
            dur_presympt_var_ub = 6.0,
            dur_presympt_shift = 0.5,
            dur_sympt_mean_lb = 2.6,
            dur_sympt_mean_ub = 3.9,
            dur_sympt_var_lb = 3.0,
            dur_sympt_var_ub = 4.5,
            dur_asympt_mean_lb = 4.0,
            dur_asympt_mean_ub = 6.0,
            dur_asympt_var_lb = 4.0,
            dur_asympt_var_ub = 6.0,
            dur_latent_min = 0.5,
            n_iters = as.numeric(input$n_iters),
            dur_quarantine_alone = as.numeric(input$dur_quarantine_alone),
            dur_quarantine_endtest = as.numeric(input$dur_quarantine_endtest),
            dur_quarantine_arrivetest = as.numeric(input$dur_quarantine_arrivetest),
            dur_quarantine_pretest = as.numeric(input$dur_quarantine_pretest),
            seed = input$RNseed,
            infection_timing = input$infection_timing,
            rr_asympt = input$rr_asympt
        )
        
        #temp
        #print(lapply(sim_params, typeof))
        print(sim_params)
        #Run simulation
        Values$dt_raw <- run_sim(sim_params, dt_incubation_dists_lnorm, progress = TRUE)
        })
    
    #RESET
    observeEvent(input$reset, {
        reset("prev")
        reset("sec_cases_per_day")
        reset("Options")
        reset("dur_quarantine_alone")
        reset("dur_quarantine_endtest")
        reset("dur_quarantine_arrivetest")
        reset("dur_quarantine_pretest")
        reset("RNseed")
        reset("sn_presympt")
        reset("sn_sympt")
        reset("sn_asympt")
        reset("prob_asympt")
        reset("prob_quarantine_compliance")
        reset("prob_isolate_sympt")
        reset("prob_isolate_test")
        reset("prob_isolate_both")
        reset("infection_timing")
        reset("rr_asympt")
        
        Values$dt_raw <- fread("dt_raw.csv")
        
        dt_analysis()
        
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
