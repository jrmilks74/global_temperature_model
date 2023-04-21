#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(forecast)

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
                Climate.fit <- lm(formula, data = global_temp)
                trendFitData <- data.frame(time = global_temp$time, temperatureFit = Climate.fit$fitted.values)
                ggplot(data = global_temp, aes(x = time, y = Temperature, colour = "Actual")) + 
                        theme_classic() +
                        geom_line() +
                        geom_smooth(method = "lm", formula = y ~ x, colour = "red") +
                        geom_line(data = trendFitData, aes(x = time, y = temperatureFit, colour = "Predicted")) +
                        geom_smooth(data = trendFitData, aes(x = time, y = temperatureFit), method = "lm", formula = y ~ x, colour = "blue") +
                        labs(x = "Time",
                             y = "Temperature anomaly (ºC)",
                             main = "Actual vs predicted global mean temperature",
                             sub = "Given selected variables")
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
                Climate.fit <- lm(formula, data = global_temp)
                summary(Climate.fit)
        }
        
        output$summary_data <- renderPrint({
                create_summary(input$independent_variables)
        })
        
        output$actual_trend <- renderPrint({
                temp_trend <- lm(Temperature ~ time, data = global_temp)
                summary(temp_trend)
        })
        
        create_predicted_trend <- function(variables) {
                observeEvent(input$goButton, {
                        variables <- variables
                })
                formula <- paste("Temperature ~", paste(variables, collapse = " + "))
                Climate.fit <- lm(formula, data = global_temp)
                trendData <- data.frame(time = global_temp$time, Temperature = Climate.fit$fitted.values)
                trendFit <- lm(Temperature ~ time, data = trendData)
                summary(trendFit)
        }
        
        output$predicted_trend <- renderPrint({
                create_predicted_trend(input$independent_variables)
        })
}
