#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(forecast)
library(mailtoR)

# Define UI for application
ui <- fluidPage(
        # Give your app a title
        titlePanel("Global Temperature Predictor"),
        
        # Create a sidebar with checkboxes for each independent variable
        sidebarLayout(
                sidebarPanel(
                        checkboxGroupInput(inputId = "independent_variables",
                                           label = "Select independent variables:",
                                           choices = c("CO2", "Solar", "ENSO", "Aerosols"),
                                           selected = c("CO2", "Solar", "ENSO", "Aerosols")),
                        actionButton("goButton", "Go!")
                ),
                
                # Display the temperature plot on the main panel
                mainPanel(
                        h3("Global temperature vs. selected variables"),
                        plotOutput("temperature_plot"),
                        h3("Model results"),
                        verbatimTextOutput("summary_data"),
                        h3("Actual temperature trend"),
                        verbatimTextOutput("actual_trend"),
                        h3("Predicted temperature trend given selected variables"),
                        verbatimTextOutput("predicted_trend")
                ),
        ),
        hr(),
        h5("Created by: Jim Milks"),
        mailtoR(email = "jrmilks@gmail.com",
                "Contact admin"),
        br(),
        "20 April 2023"
)