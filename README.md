# Global Temperature Predictor — Shiny App

A single-file Shiny app that models global surface temperature anomalies using a set of physical drivers with optional lag optimization:

- **CO₂ radiative forcing** (from Mauna Loa ppm: `RF = 5.35 * ln(C/280)`)
- **Solar irradiance** (SATIRE-T pre-satellite + SATIRE-S satellite era)
- **ENSO** (Niño 3.4 index)
- **Stratospheric aerosols** (CREST NetCDF, with GISS Sato fallback)

The app searches over user-defined lag windows, ranks models by **AICc**, and provides diagnostics, Newey–West trend estimates, and a quick CCF preview.

---

## Features

- 📦 **Self-contained**: one `app.R` file (no external helpers required)
- 🔁 **Auto lag search** per predictor with guardrails on compute size
- 📈 **Observed vs predicted** plot for best model
- 🧪 **Diagnostics**: residual ACF + histogram
- 📊 **Top-10 models** by AICc
- 📐 **Trend estimates** (Observed & Fitted) with Newey–West SEs
- 🔍 **CCF preview** to eyeball lead/lag structure
- 🧰 **Robust loaders** with fallbacks for all data sources
- ✅ **JSON-safe UI** (no `jsonlite` named-vector warnings)
- 🛡️ **Lag safety**: removing predictors will not raise `Unknown ... lag_*` warnings

---

## Data Sources (fetched live at runtime)

- **GISTEMP** Land+Ocean monthly anomalies (NASA)
- **CO₂ (Mauna Loa)** monthly (NOAA GML)
- **SATIRE** TSI: SATIRE-T (historical) + SATIRE-S (satellite era)
- **ENSO Niño 3.4** (NOAA PSL; CPC fallback)
- **Stratospheric aerosols**: **CREST** NetCDF (global SAOD) or **GISS Sato** (fallback)

> The app automatically aligns time coverage and drops missing values after merges.

---

## Requirements

- **R >= 4.1** recommended  
- **Packages**: `shiny`, `tidyverse`, `lubridate`, `readr`, `broom`,  
  `AICcmodavg`, `lmtest`, `sandwich`, `scales`, `ncdf4`

Install once:

```r
install.packages(c(
  "shiny","tidyverse","lubridate","readr","broom",
  "AICcmodavg","lmtest","sandwich","scales","ncdf4"
))
