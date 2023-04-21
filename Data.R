library(tidyverse)
library(tsbox)
library(imputeTS)
library(Hmisc)

# Temperature data
GISS <- read_table("https://data.giss.nasa.gov/gistemp/tabledata_v4/GLB.Ts+dSST.txt",
                   skip = 7) %>%
        slice(1:(n() - 5)) %>%
        filter(!row_number() %in% c(22, 43, 64, 85, 106, 127, 148))
GISS[GISS == "****"] <- NA

GISS <- GISS %>% select("Year", 
                        "Jan", 
                        "Feb", 
                        "Mar", 
                        "Apr", 
                        "May", 
                        "Jun", 
                        "Jul", 
                        "Aug", 
                        "Sep", 
                        "Oct", 
                        "Nov", 
                        "Dec")

GISS$Year <- as.numeric(GISS$Year)
GISS$Jan <- as.numeric(GISS$Jan)
GISS$Feb <- as.numeric(GISS$Feb)
GISS$Mar <- as.numeric(GISS$Mar)
GISS$Apr <- as.numeric(GISS$Apr)
GISS$May <- as.numeric(GISS$May)
GISS$Jun <- as.numeric(GISS$Jun)
GISS$Jul <- as.numeric(GISS$Jul)
GISS$Aug <- as.numeric(GISS$Aug)
GISS$Sep <- as.numeric(GISS$Sep)
GISS$Oct <- as.numeric(GISS$Oct)
GISS$Nov <- as.numeric(GISS$Nov)
GISS$Dec <- as.numeric(GISS$Dec)

GISS <- GISS %>%
        pivot_longer(c(Jan, 
                       Feb, 
                       Mar, 
                       Apr,
                       May,
                       Jun,
                       Jul,
                       Aug,
                       Sep,
                       Oct,
                       Nov,
                       Dec), names_to = "Month", values_to = "anomaly") %>%
        mutate(time = dmy(paste(1, Month, Year))) %>%
        select(time, anomaly)
GISS$anomaly <- GISS$anomaly/100
GISS <- subset(GISS, time >= "1958-03-01") %>%
        rename(Temperature = anomaly)

# CO2 data
co2 <- read_table("https://gml.noaa.gov/webdata/ccgg/trends/co2/co2_mm_mlo.txt",
                  col_names = FALSE, skip = 58)

co2 <- co2 %>%
        select(X1, X2, X4) %>%
        rename(Year = X1, Month = X2, CO2 = X4) %>%
        mutate(time = make_date(year = Year, month = Month, day = 1 )) %>%
        mutate(CO2_rf = 5.35*log(CO2/280)) %>%
        select(time, CO2_rf) %>%
        rename(CO2 = CO2_rf)



# Sunspot data
sunspots <- read_table("https://www.sidc.be/silso/DATA/SN_m_tot_V2.0.txt",
                       col_names = FALSE)
sunspots <- sunspots %>%
        rename(Year = X1, Month = X2, sunspots = X4) %>%
        mutate(time = make_date(year = Year, month = Month, day = 1)) %>%
        select(time, sunspots) %>%
        subset(time >= "1958-03-01") %>%
        rename("Solar" = sunspots)

# ENSO data

ENSO <- read_table("https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt")
ENSO$time <- seq.Date(from = as.Date("1950-01-01"), by = "month", length.out = nrow(ENSO))
ENSO <- subset(ENSO, time >= "1957-09-01") %>%
        rename(ENSO = ANOM) %>%
        select(time, ENSO)
ENSO$ENSO_lagged <- Lag(ENSO$ENSO, 3)
ENSO <- ENSO %>%
        select(time, ENSO_lagged) %>%
        rename(ENSO = ENSO_lagged)

# Aerosols data

aerosols <- read_table("https://data.giss.nasa.gov/modelforce/strataer/tau.line_2012.12.txt",
                       skip = 3) %>%
        rename(time = "year/mon", Aerosols = global) %>%
        mutate_at(c("time", "Aerosols"), as.numeric) %>%
        select(time, Aerosols)

aerosols$time <- format(date_decimal(aerosols$time), "%Y-%m-01") %>%
        as.Date(aerosols$time, format = "%Y-%m-%d")
aerosols <- subset(aerosols, time >= "1955-10-01")

## Extrapolate out from 2012-09-01
length_out <- subset(ENSO, time >= "2012-10-01")
length_out <- length(length_out$time)

extrapolation <- data.frame(
        time = seq(from = as.Date("2012-10-01"), by = "1 month", length = length_out),
        Aerosols = "NaN"
)

aerosols <- rbind(aerosols, extrapolation)
aerosols$Aerosols <- as.numeric(aerosols$Aerosols)
aerosols <- na_interpolation(aerosols, option = "linear")
aerosols$Aerosols_lagged <- Lag(aerosols$Aerosols, 18)
aerosols <- aerosols %>%
        select(time, Aerosols_lagged) %>%
        rename(Aerosols = Aerosols_lagged)
aerosols <- subset(aerosols, time >= "1958-03-01")

# Make data frame

df_list <- list(GISS, co2, sunspots, ENSO, aerosols)
global_temp <- df_list %>%
        reduce(inner_join, by = "time") %>%
        mutate(date = decimal_date(time)) %>%
        select(date, Temperature, CO2, Solar, ENSO, Aerosols) %>%
        rename(time = date)

# Lag function
calculate_lag <- function(timeseries1, timeseries2) {
        cross_corr <- stats::ccf(timeseries1, timeseries2)
        lag <- which.max(cross_corr$acf)
        return(lag)
}
