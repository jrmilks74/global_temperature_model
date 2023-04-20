#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Load the data

source("Data.R")

# Define server logic required to draw a histogram
function(input, output) {
        
        # Define a function to create the temperature model
        create_plot <- function(variables) {
                observeEvent(input$goButton, {
                        variables <- variables
                })
                formula <- paste("Temperature ~", paste(variables, collapse = " + "))
                Climate.fit <- tslm(formula, data = global_temp)
                autoplot(global_temp[, "Temperature"], series = "Actual") + 
                        autolayer(fitted(Climate.fit), series = "Model") +
                        xlab("Year") +
                        ylab("Temperature anomaly (ºC)") +
                        ggtitle("Predicted change in global mean temperature based on selected variables") +
                        guides(colour = guide_legend(title = " "))
        }
        
        output$temperature_plot <- renderPlot({
                # Call the create_plot function with the selected variables
                create_plot(input$independent_variables)
        })
        
        create_summary <- function(variables) {
                observeEvent(input$goButton, {
                        variables <- variables
                })
                formula <- paste("Temperature ~", paste(variables, collapse = " + "))
                Climate.fit <- tslm(formula, data = global_temp)
                summary(Climate.fit)
        }
        
        output$summary_data <- renderPrint({
                create_summary(input$independent_variables)
        })
}
