# app.R — Global Temperature Predictor (clean, JSON-safe, lag-safe)
library(shiny)
library(tidyverse)
library(lubridate)
library(readr)
library(broom)
library(AICcmodavg)
library(lmtest)
library(sandwich)
library(scales)
library(ncdf4)
library(mailtoR)

# Ensure Shiny uses JSON named lists (prevents jsonlite named-vector warnings)
options(shiny.json.use.named.list = TRUE)

# ---------- Small helpers ----------
read_lines_safely <- function(url) tryCatch(readr::read_lines(url, progress = FALSE), error = function(e) character())
read_table_safely <- function(url, ...) tryCatch(readr::read_table(url, ...), error = function(e) tibble())
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---------- Data loaders ----------
# 1) GISTEMP Land+Ocean monthly (robust parser)
load_gistemp <- function(min_date = as.Date("1958-03-01")) {
        url <- "https://data.giss.nasa.gov/gistemp/tabledata_v4/GLB.Ts+dSST.txt"
        lines <- read_lines_safely(url)
        if (!length(lines)) stop("GISTEMP monthly: download failed", call. = FALSE)
        
        rows <- lines[grepl("^\\s*\\d{4}\\s", lines)]
        if (!length(rows)) stop("GISTEMP monthly: no year rows found", call. = FALSE)
        
        split_rows <- strsplit(trimws(rows), "\\s+")
        get_months12 <- function(x) {
                m <- x[-1L]
                if (length(m) < 12L) m <- c(m, rep(NA_character_, 12L - length(m)))
                if (length(m) > 12L) m <- m[1:12]
                as.numeric(ifelse(m == "****", NA, m))
        }
        years <- as.integer(vapply(split_rows, `[[`, "", 1))
        months_mat <- do.call(rbind, lapply(split_rows, get_months12))
        colnames(months_mat) <- month.abb
        
        tibble(Year = years) |>
                bind_cols(as_tibble(months_mat)) |>
                pivot_longer(all_of(month.abb), names_to = "Mon", values_to = "anomaly_x100") |>
                mutate(time = dmy(paste0("01-", Mon, "-", Year)), Temperature = anomaly_x100 / 100) |>
                select(time, Temperature) |>
                arrange(time) |>
                filter(time >= min_date) |>
                tidyr::drop_na() -> out
        
        if (!nrow(out)) stop("GISTEMP monthly: parsed, but no rows after filtering", call. = FALSE)
        out
}

# 2) Solar TSI: SATIRE-T + SATIRE-S → monthly mean
read_satire_t <- function(url = "https://www2.mps.mpg.de/projects/sun-climate/data/SATIRE-T_wu18_tsi.txt") {
        raw <- read_lines_safely(url); raw <- raw[!grepl("^\\s*;", raw)]; raw <- raw[nzchar(trimws(raw))]
        if (!length(raw)) return(tibble(JD = numeric(), TSI = numeric(), date = as.Date(character())))
        readr::read_table(
                paste(raw, collapse = "\n"),
                col_names = c("JD", "TSI"),
                col_types = cols(JD = col_double(), TSI = col_double()),
                na = c("NaN","nan"), show_col_types = FALSE
        ) |>
                mutate(date = as.Date(JD - 2440587.5, origin = "1970-01-01")) |>
                filter(is.finite(TSI), TSI > 1200, TSI < 1500)
}
read_satire_s <- function(url = "https://www2.mps.mpg.de/projects/sun-climate/data/SATIRE/SATIRE-S/SATIRE-S_TSI_latest.txt") {
        raw <- read_lines_safely(url); raw <- raw[!grepl("^\\s*;", raw)]; raw <- raw[nzchar(trimws(raw))]
        if (!length(raw)) return(tibble(JD = numeric(), TSI = numeric(), date = as.Date(character())))
        readr::read_table(
                paste(raw, collapse = "\n"),
                col_names = c("JD", "TSI", "lo", "hi", "source"),
                col_types = cols(JD = col_double(), TSI = col_double(), lo = col_double(), hi = col_double(), source = col_double()),
                na = c("NaN"), show_col_types = FALSE
        ) |>
                mutate(date = as.Date(JD - 2440587.5, origin = "1970-01-01")) |>
                filter(is.finite(JD), is.finite(TSI), TSI > 1200, TSI < 1500)
}
load_satire_tsi <- function(min_date = as.Date("1958-03-01")) {
        sat_t <- read_satire_t(); sat_s <- read_satire_s()
        if (nrow(sat_t) + nrow(sat_s) == 0) stop("SATIRE data not available.", call. = FALSE)
        cutoff <- as.Date("1974-08-22")
        bind_rows(dplyr::filter(sat_t, date <= cutoff), dplyr::filter(sat_s, date > cutoff)) |>
                mutate(time = floor_date(date, "month")) |>
                group_by(time) |>
                summarise(Solar = mean(TSI, na.rm = TRUE), .groups = "drop") |>
                filter(time >= min_date)
}

# 3) ENSO — Niño 3.4 monthly
load_oni_enso <- function() {
        page <- "https://psl.noaa.gov/data/timeseries/month/Nino34/"
        txt  <- read_lines_safely(page)
        rx   <- "^\\s*(\\d{4})\\s+(\\d{1,2})\\s+([-+]?\\d*\\.?\\d+(?:[eE][-+]?\\d+)?)\\s*$"
        if (length(txt) && any(grepl(rx, txt))) {
                mm <- regmatches(txt, regexec(rx, txt, perl = TRUE)); mm <- mm[lengths(mm) > 0]
                yr <- as.integer(vapply(mm, `[[`, "", 2)); mo <- as.integer(vapply(mm, `[[`, "", 3)); val <- as.numeric(vapply(mm, `[[`, "", 4))
                return(tibble(time = make_date(yr, mo, 1), ENSO = val) |> arrange(time) |> filter(time >= as.Date("1957-09-01")))
        }
        url <- "https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt"
        dat <- tryCatch(readr::read_table(url, show_col_types = FALSE), error = function(e) NULL)
        if (!is.null(dat) && "ANOM" %in% names(dat)) {
                n <- nrow(dat)
                tibble(time = seq.Date(as.Date("1950-01-01"), by = "month", length.out = n),
                       ENSO = as.numeric(dat$ANOM)) |> filter(time >= as.Date("1957-09-01"))
        } else stop("ENSO monthly: could not parse from PSL or CPC")
}

# 4) CO2 → radiative forcing
load_co2_rf <- function() {
        url <- "https://gml.noaa.gov/webdata/ccgg/trends/co2/co2_mm_mlo.txt"
        raw <- read_table_safely(url, col_names = FALSE, skip = 58, col_types = cols(.default = col_double()))
        if (!nrow(raw)) stop("CO2 monthly: failed to download", call. = FALSE)
        raw |>
                select(Year = X1, Month = X2, ppm = X4) |>
                mutate(time = make_date(Year, Month, 1), CO2 = 5.35 * log(ppm / 280)) |>
                select(time, CO2) |>
                arrange(time) |>
                tidyr::drop_na()
}

# 5) Aerosols — CREST NetCDF (default) or GISS Sato fallback
load_aerosols_crest <- function(min_date = as.Date("1958-03-01"), use_interpolated = FALSE) {
        url <- "https://fmi.b2share.csc.fi/api/files/19abe20f-7409-40b8-bc85-98759dd0ed7d/ESA-CREST-v1.0-fv001.nc"
        tf  <- tempfile(fileext = ".nc"); utils::download.file(url, tf, mode = "wb", quiet = TRUE)
        nc <- nc_open(tf); on.exit(try(nc_close(nc), silent = TRUE), add = TRUE)
        
        time_v <- ncvar_get(nc, "time")
        units  <- ncatt_get(nc, "time", "units")$value
        origin <- sub(".*since\\s+", "", units); origin <- sub("\\s+0:0:0(\\.0)?$", "", origin)
        dates  <- as.Date(origin) + as.integer(time_v); time_m <- floor_date(dates, "month")
        
        if (!use_interpolated) {
                saod_g <- ncvar_get(nc, "Global_SAOD"); aeros <- tibble(time = time_m, Aerosols = as.numeric(saod_g))
        } else {
                saod_lt <- ncvar_get(nc, "SAOD_interpolated"); lats <- ncvar_get(nc, "latitude_centers")
                w <- cos(lats * pi/180); w <- w / sum(w); glob <- as.numeric(w %*% saod_lt)
                aeros <- tibble(time = time_m, Aerosols = glob)
        }
        
        aeros |> arrange(time) |> mutate(Aerosols = dplyr::lag(Aerosols, 18)) |> filter(time >= min_date) |> tidyr::drop_na()
}
load_aerosols_giss <- function(min_date = as.Date("1958-03-01")) {
        url <- "https://data.giss.nasa.gov/modelforce/strataer/tau.line_2012.12.txt"
        sato <- read_table_safely(url, skip = 3, col_types = cols(.default = col_double()))
        if (!nrow(sato)) stop("Aerosols: failed to download GISS Sato file", call. = FALSE)
        tcol <- names(sato)[1]
        gcol <- if ("global" %in% names(sato)) "global" else names(sato)[grepl("glob", names(sato), ignore.case = TRUE)][1]
        sato |>
                rename(time_dec = all_of(tcol), Aerosols = all_of(gcol)) |>
                mutate(time = date_decimal(time_dec), time = as.Date(format(time, "%Y-%m-01"))) |>
                arrange(time) -> aeros
        last_t <- max(aeros$time, na.rm = TRUE)
        if (last_t < floor_date(Sys.Date(), "month")) {
                ext_times <- seq(from = last_t + months(1), to = floor_date(Sys.Date(), "month"), by = "1 month")
                aeros <- bind_rows(aeros, tibble(time = ext_times, Aerosols = 0.0001))
        }
        aeros |> mutate(Aerosols = dplyr::lag(Aerosols, 18)) |> filter(time >= min_date) |> tidyr::drop_na()
}
load_aerosols <- function(source = c("crest", "giss"),
                          min_date = as.Date("1958-03-01"),
                          use_interpolated = FALSE) {
        source <- match.arg(source)
        if (source == "crest") {
                tryCatch(
                        load_aerosols_crest(min_date = min_date, use_interpolated = use_interpolated),
                        error = function(e) { warning("CREST NetCDF failed (", e$message, "); falling back to GISS Sato.")
                                load_aerosols_giss(min_date = min_date) }
                )
        } else {
                load_aerosols_giss(min_date = min_date)
        }
}

# ---------- Merge ----------
build_dataset <- function(start_year, aerosol_source, use_interpolated) {
        temp <- load_gistemp(); co2 <- load_co2_rf(); sol <- load_satire_tsi(); enso <- load_oni_enso()
        aer  <- load_aerosols(aerosol_source, min_date = as.Date("1958-03-01"), use_interpolated = use_interpolated)
        list(temp, co2, sol, enso, aer) |>
                reduce(inner_join, by = "time") |>
                arrange(time) |>
                mutate(t_num = decimal_date(time)) |>
                filter(year(time) >= start_year) |>
                tidyr::drop_na()
}

# ---------- Modeling ----------
enumerate_models <- function(predictors) {
        if (!length(predictors)) return(character())
        unlist(lapply(1:length(predictors), function(k)
                combn(predictors, k, FUN = function(x) paste(x, collapse = "+"))),
               use.names = FALSE)
}
fit_and_score <- function(df, spec) {
        f <- as.formula(paste("Temperature ~", spec)); m <- lm(f, data = df)
        tibble(formula = deparse(f), adj_r2 = summary(m)$adj.r.squared, AICc = AICc(m), model = list(m))
}
nw_trend <- function(y, tnum, lag = 12) {
        m <- lm(y ~ tnum); co <- coeftest(m, vcov. = NeweyWest(m, lag = lag, prewhite = FALSE)); list(model = m, coefs = co)
}

# ---------- Lags (warning-proof) ----------
apply_lags <- function(df,
                       lag_co2 = 0L, lag_solar = 0L, lag_enso = 0L, lag_aer = 0L,
                       use_vars = c("CO2","Solar","ENSO","Aerosols")) {
        dfl <- df
        if ("CO2" %in% use_vars)      dfl$CO2_L      <- dplyr::lag(dfl$CO2,      n = as.integer(lag_co2))
        if ("Solar" %in% use_vars)    dfl$Solar_L    <- dplyr::lag(dfl$Solar,    n = as.integer(lag_solar))
        if ("ENSO" %in% use_vars)     dfl$ENSO_L     <- dplyr::lag(dfl$ENSO,     n = as.integer(lag_enso))
        if ("Aerosols" %in% use_vars) dfl$Aerosols_L <- dplyr::lag(dfl$Aerosols, n = as.integer(lag_aer))
        
        name_map <- c(CO2 = "CO2_L", Solar = "Solar_L", ENSO = "ENSO_L", Aerosols = "Aerosols_L")
        needed <- unique(c("Temperature", "t_num", unname(name_map[intersect(use_vars, names(name_map))])))
        tidyr::drop_na(dfl, dplyr::all_of(needed))
}
build_lag_grid <- function(preds, minmax, step) {
        seqr <- function(lo, hi, st) seq.int(as.integer(lo), as.integer(hi), by = max(1L, as.integer(st)))
        grids <- list()
        if ("CO2"      %in% preds) grids$lag_co2   <- seqr(minmax$co2_min,   minmax$co2_max,   step)
        if ("Solar"    %in% preds) grids$lag_solar <- seqr(minmax$sol_min,   minmax$sol_max,   step)
        if ("ENSO"     %in% preds) grids$lag_enso  <- seqr(minmax$enso_min,  minmax$enso_max,  step)
        if ("Aerosols" %in% preds) grids$lag_aer   <- seqr(minmax$aer_min,   minmax$aer_max,   step)
        if (length(grids) == 0L) return(tibble(dummy = 0L) |> select(-dummy))
        tidyr::expand_grid(!!!grids)
}
fit_with_lags <- function(df, preds, spec_vars, lags_row) {
        # Safely pull lag values present in this grid row
        get_lag <- function(nm) if (!is.null(lags_row[[nm]])) as.integer(lags_row[[nm]]) else 0L
        
        dfl <- apply_lags(
                df,
                lag_co2   = get_lag("lag_co2"),
                lag_solar = get_lag("lag_solar"),
                lag_enso  = get_lag("lag_enso"),
                lag_aer   = get_lag("lag_aer"),
                use_vars  = spec_vars
        )
        
        name_map <- c(CO2 = "CO2_L", Solar = "Solar_L", ENSO = "ENSO_L", Aerosols = "Aerosols_L")
        use_cols <- unname(name_map[spec_vars])
        
        f <- as.formula(paste("Temperature ~", paste(use_cols, collapse = " + ")))
        m <- lm(f, data = dfl)
        
        tibble(
                formula   = paste("Temperature ~", paste(spec_vars, collapse = " + ")),
                adj_r2    = summary(m)$adj.r.squared,
                AICc      = AICc(m),
                model     = list(m),
                n         = nrow(dfl),
                lag_co2   = get_lag("lag_co2"),
                lag_solar = get_lag("lag_solar"),
                lag_enso  = get_lag("lag_enso"),
                lag_aer   = get_lag("lag_aer")
        )
}

# ---------- UI ----------
ui <- fluidPage(
        titlePanel("Global Temperature Regression"),
        sidebarLayout(
                sidebarPanel(
                        checkboxGroupInput(
                                "iv", "Select independent variables:",
                                choices  = list("CO2"="CO2","Solar"="Solar","ENSO"="ENSO","Aerosols"="Aerosols"),
                                selected = c("CO2","Solar","ENSO","Aerosols")
                        ),
                        selectInput(
                                  "start_year",
                                  label = "Analysis start year:",
                                  choices = 1958:year(Sys.Date()),
                                  selected = 1979
                                ),
                        radioButtons("aer_source", "Aerosol dataset:",
                                     choices = list("CREST (NetCDF, 750 nm)"="crest", "GISS Sato (to 2012, bg extended)"="giss"),
                                     selected = "crest"),
                        checkboxInput("aer_interp", "Use CREST latitude-weighted (interpolated) global mean", FALSE),
                        tags$hr(),
                        strong("Lag search (months)"),
                        checkboxInput("auto_lags", "Auto-optimize lags", TRUE),
                        helpText("Set min/max lags per predictor; larger ranges increase runtime."),
                        fluidRow(
                                column(6, sliderInput("lag_co2_min", "CO₂ min", min = 0, max = 48, value = 0,  step = 1, sep = "")),
                                column(6, sliderInput("lag_co2_max", "CO₂ max", min = 0, max = 48, value = 24, step = 1, sep = ""))
                        ),
                        fluidRow(
                                column(6, sliderInput("lag_solar_min", "Solar min", min = 0, max = 48, value = 1,  step = 1, sep = "")),
                                column(6, sliderInput("lag_solar_max", "Solar max", min = 0, max = 48, value = 24, step = 1, sep = ""))
                        ),
                        fluidRow(
                                column(6, sliderInput("lag_enso_min", "ENSO min", min = 0, max = 24, value = 1,  step = 1, sep = "")),
                                column(6, sliderInput("lag_enso_max", "ENSO max", min = 0, max = 24, value = 12, step = 1, sep = ""))
                        ),
                        fluidRow(
                                column(6, sliderInput("lag_aerosol_min", "Aerosol min", min = 0, max = 48, value = 12, step = 1, sep = "")),
                                column(6, sliderInput("lag_aerosol_max", "Aerosol max", min = 0, max = 48, value = 36, step = 1, sep = ""))
                        ),
                        sliderInput("lag_step", "Lag step", min = 1, max = 6, value = 3, step = 1, sep = "", ticks = TRUE),
                        actionButton("go", "Find Best Model")
                ),
                mainPanel(
                        h4("Observed vs Predicted (best model)"),
                        plotOutput("plot_fit", height = "380px"),
                        br(),
                        tabsetPanel(
                                tabPanel("Best Model", verbatimTextOutput("best_formula"), tableOutput("best_stats"), verbatimTextOutput("best_summary")),
                                tabPanel("Top 10 Models", tableOutput("top_table")),
                                tabPanel("Trends",
                                         h5("Observed temperature trend (°C/year)"), verbatimTextOutput("obs_trend"),
                                         h5("Predicted (fitted) trend from best model (°C/year)"), verbatimTextOutput("fit_trend")),
                                tabPanel("Diagnostics", plotOutput("resid_acf", height = "280px"), plotOutput("resid_hist", height = "280px")),
                                tabPanel("Data coverage", tableOutput("coverage_table")),
                                tabPanel("Selected lags", tableOutput("best_lag_table")),
                                tabPanel("Lag CCF preview",
                                         h6("Cross-correlation (temperature vs predictor) for current date range"),
                                         selectInput("ccf_var", "Predictor:",
                                                     choices = list("CO2"="CO2","Solar"="Solar","ENSO"="ENSO","Aerosols"="Aerosols"),
                                                     selected = "ENSO"),
                                         sliderInput("ccf_maxlag", "Max lag (months)", min = 6, max = 60, value = 24, step = 1, sep = ""),
                                         plotOutput("ccf_plot", height = "280px"),
                                         helpText("Positive lag means predictor leads temperature (temp responds later).")
                                )
                        ),
                        hr(),
                        tags$small(
                                "Sources: NASA/PSL GISTEMP monthly; NOAA GML CO2 monthly; SATIRE-T & SATIRE-S TSI; ",
                                "PSL Niño 3.4 monthly; CREST NetCDF (or GISS Sato) aerosols. Updated ",
                                format(Sys.Date(), "%d %b %Y")
                        )
                )
        ),
        hr(),
        h5("Created by: Jim Milks"),
        mailtoR(email = "jrmilks@gmail.com",
                "Contact admin"),
        br(),
        "9 Sep 2025"
)

# ---------- Server ----------
server <- function(input, output, session) {
        
        base_data <- reactiveVal(NULL)
        observeEvent(input$go, {
                df <- build_dataset(input$start_year, input$aer_source, input$aer_interp)
                validate(need(nrow(df) > 24, "Not enough overlapping monthly records after start year."))
                base_data(df)
        }, ignoreInit = TRUE)
        
        results <- eventReactive(input$go, {
                df <- base_data()
                validate(need(!is.null(df), "Click 'Find Best Model' to run."))
                preds <- intersect(input$iv, c("CO2","Solar","ENSO","Aerosols"))
                validate(need(length(preds) > 0, "Select at least one predictor."))
                
                specs <- enumerate_models(preds)
                
                # Build lag grid respecting selected predictors
                minmax <- list(
                        co2_min  = input$lag_co2_min,   co2_max  = input$lag_co2_max,
                        sol_min  = input$lag_solar_min, sol_max  = input$lag_solar_max,
                        enso_min = input$lag_enso_min,  enso_max = input$lag_enso_max,
                        aer_min  = input$lag_aerosol_min, aer_max = input$lag_aerosol_max
                )
                grid <- build_lag_grid(preds, minmax, step = input$lag_step)
                
                # Guardrail on total fits
                max_models <- 6000L
                combs <- nrow(grid) * max(1L, length(specs))
                if (combs > max_models) {
                        new_step <- ceiling(input$lag_step * sqrt(combs / max_models))
                        grid <- build_lag_grid(preds, minmax, step = new_step)
                }
                
                if (isTRUE(input$auto_lags)) {
                        scored <- purrr::map_dfr(seq_len(nrow(grid)), function(i) {
                                lr <- grid[i, , drop = FALSE]
                                purrr::map_dfr(specs, function(s) {
                                        spec_vars <- trimws(unlist(strsplit(s, "\\+"), use.names = FALSE))
                                        fit_with_lags(df, preds, spec_vars, lr)
                                })
                        }) |> arrange(AICc, desc(adj_r2))
                        
                        best <- scored |> slice(1)
                        best_spec_vars <- trimws(unlist(strsplit(gsub("^Temperature ~\\s*", "", best$formula[[1]]), "\\+"), use.names = FALSE))
                        
                        df_best <- apply_lags(
                                df,
                                lag_co2   = best$lag_co2,
                                lag_solar = best$lag_solar,
                                lag_enso  = best$lag_enso,
                                lag_aer   = best$lag_aer,
                                use_vars  = best_spec_vars
                        )
                        
                        list(df = df_best, scored = scored, best = best)
                } else {
                        scored <- tibble(spec = specs) |>
                                mutate(res = purrr::map(spec, ~fit_and_score(df, .x))) |>
                                tidyr::unnest(res) |>
                                arrange(AICc, desc(adj_r2))
                        best <- scored |> slice(1)
                        list(df = df, scored = scored, best = best)
                }
        })
        
        # Coverage table
        output$coverage_table <- renderTable({
                df <- base_data(); if (is.null(df)) return(NULL)
                data.frame(Series = c("Temperature","CO2","Solar","ENSO","Aerosols"),
                           Start = rep(min(df$time), 5), End = rep(max(df$time), 5))
        })
        
        # Best model outputs
        output$best_formula <- renderText({
                rz <- results(); if (is.null(rz)) return("Click 'Find Best Model' to run.")
                paste("Best model (by AICc):", rz$best$formula[[1]])
        })
        output$best_stats <- renderTable({
                rz <- results(); if (is.null(rz)) return(NULL)
                rz$best |> transmute(AICc = round(AICc, 2), adj_R2 = round(adj_r2, 3))
        })
        output$best_summary <- renderPrint({
                rz <- results(); if (is.null(rz)) return(invisible())
                summary(rz$best$model[[1]])
        })
        
        # Top 10
        output$top_table <- renderTable({
                rz <- results(); if (is.null(rz)) return(NULL)
                rz$scored |>
                        mutate(rank = dplyr::row_number(), AICc = round(AICc, 2), adj_r2 = round(adj_r2, 3)) |>
                        select(rank, formula, AICc, adj_r2) |>
                        slice(1:10)
        })
        
        # Plot observed vs predicted
        output$plot_fit <- renderPlot({
                rz <- results(); if (is.null(rz)) return(invisible())
                m  <- rz$best$model[[1]]
                df <- rz$df |> mutate(Predicted = as.numeric(predict(m, newdata = rz$df)))
                df |>
                        transmute(time, `Actual temperature` = Temperature, `Predicted (model)` = Predicted) |>
                        pivot_longer(-time, names_to = "Series", values_to = "Value") |>
                        ggplot(aes(time, Value, color = Series)) +
                        geom_line(linewidth = 0.8, na.rm = TRUE) +
                        labs(x = NULL, y = "Temperature anomaly (°C)",
                             title = "Observed vs Predicted Global Temperature",
                             subtitle = "Best linear model on selected drivers (AICc)") +
                        theme_classic() +
                        guides(color = guide_legend(title = NULL))
        })
        
        # Diagnostics
        output$resid_acf <- renderPlot({
                rz <- results(); if (is.null(rz)) return(invisible())
                acf(resid(rz$best$model[[1]]), main = "Residual ACF (best model)")
        })
        output$resid_hist <- renderPlot({
                rz <- results(); if (is.null(rz)) return(invisible())
                ggplot(tibble(r = resid(rz$best$model[[1]])), aes(r)) +
                        geom_histogram(bins = 30) +
                        labs(x = "Residuals", y = "Count", title = "Residual histogram") +
                        theme_minimal()
        })
        
        # Trends (Newey–West)
        output$obs_trend <- renderPrint({
                df <- base_data(); if (is.null(df)) return(invisible())
                res <- nw_trend(df$Temperature, df$t_num)
                cat("Slope =", round(coef(res$model)[2], 5), "°C per year\n"); print(res$coefs)
        })
        output$fit_trend <- renderPrint({
                rz <- results(); if (is.null(rz)) return(invisible())
                m  <- rz$best$model[[1]]
                df <- rz$df |> mutate(Pred = as.numeric(predict(m, newdata = rz$df)))
                res <- nw_trend(df$Pred, df$t_num)
                cat("Slope =", round(coef(res$model)[2], 5), "°C per year\n"); print(res$coefs)
        })
        
        # Selected lags table (only variables present in best model)
        output$best_lag_table <- renderTable({
                rz <- results(); if (is.null(rz)) return(NULL)
                spec_vars <- trimws(unlist(strsplit(gsub("^Temperature ~\\s*", "", rz$best$formula[[1]]), "\\+"), use.names = FALSE))
                lag_map <- c(CO2 = "lag_co2", Solar = "lag_solar", ENSO = "lag_enso", Aerosols = "lag_aer")
                tibble::tibble(
                        Variable      = spec_vars,
                        Optimized_lag = vapply(lag_map[spec_vars], function(nm) {
                                if (!is.null(rz$best[[nm]])) as.integer(rz$best[[nm]]) else 0L
                        }, integer(1))
                )
        })
        
        # CCF preview
        output$ccf_plot <- renderPlot({
                df <- base_data(); if (is.null(df)) return(invisible())
                v <- switch(input$ccf_var, "CO2" = df$CO2, "Solar" = df$Solar, "ENSO" = df$ENSO, "Aerosols" = df$Aerosols)
                ok <- is.finite(df$Temperature) & is.finite(v)
                if (sum(ok) < 24) { plot.new(); title("Not enough overlapping data"); return(invisible()) }
                stats::ccf(v[ok], df$Temperature[ok], lag.max = input$ccf_maxlag,
                           main = paste("CCF:", input$ccf_var, "→ Temp"))
        })
}

shinyApp(ui, server)
