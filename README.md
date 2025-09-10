# Global Temperature Regression — Shiny App

A single-file Shiny app that models global surface temperature anomalies using a set of physical drivers with optional lag optimization:

- **CO₂ radiative forcing** (from Mauna Loa ppm: `RF = 5.35 * ln(C/280)`)
- **Solar irradiance** (SATIRE-T pre-satellite + SATIRE-S satellite era)
- **ENSO** (Niño 3.4 index)
- **Stratospheric aerosols** (CREST NetCDF, with GISS Sato fallback)
- **Volcanic forcing** (optional analytic pulse regressor built from major eruptions since 1958)

The app searches over user-defined lag windows, ranks models by **AICc**, and provides diagnostics, Newey–West trend estimates, and a quick CCF preview.

---

## Features

- 📦 **Self-contained**: one `app.R` file (no external helpers required)
- 🔁 **Auto lag search** per predictor with guardrails on compute size
- 📈 **Observed vs predicted** plot for best model
- 🌋 **Volcanic signal**:
  - Optional volcanic-pulse regressor with exponential decay (τ months)
  - Plot shading for major eruptions (Agung, El Chichón, Pinatubo, etc.)
- 🧪 **Diagnostics**: residual ACF + histogram
- 📊 **Top-10 models** by AICc
- 📐 **Trend estimates** (Observed & Fitted) with Newey–West SEs
- 🔍 **CCF preview** to eyeball lead/lag structure
- 🧰 **Robust loaders** with fallbacks for all data sources
- ✅ **JSON-safe UI** (no `jsonlite` named-vector warnings)
- 🛡️ **Lag safety**: removing predictors will not raise `Unknown ... lag_*` warnings
- 🎛️ **Improved UX**: dropdown for start year (no confusing sliders)

---

## Data Sources (fetched live at runtime)

- **GISTEMP** Land+Ocean monthly anomalies (NASA)
- **CO₂ (Mauna Loa)** monthly (NOAA GML)
- **SATIRE** TSI: SATIRE-T (historical) + SATIRE-S (satellite era)
- **ENSO Niño 3.4** (NOAA PSL; CPC fallback)
- **Stratospheric aerosols**: **CREST** NetCDF (global SAOD) or **GISS Sato** (fallback)
- **Volcanic events** (hardcoded list of post-1958 eruptions with stratospheric impact):
  - Agung (1963–64), El Chichón (1982), Pinatubo (1991), plus smaller signals (Fuego 1974, Kasatochi 2008, Sarychev 2009, Nabro 2011, Calbuco 2015, Raikoke 2019, Hunga Tonga 2022)

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
```

## How Volcanic Signal Works

The app defines a pulse regressor:

$$
V(t) = \sum_i w_i \,\exp\!\left(-\frac{\max(0,\, t - t_i)}{\tau}\right)
$$

- \(t_i\): eruption onset month  
- \(w_i\): eruption weight (positive = cooling, negative = warming)  
- \(\tau\): decay constant in months (default 18)  

On the plot, eruptions can also be shaded for visual inspection.

> If your markdown viewer doesn’t support LaTeX, use this plain-text form:
> `V(t) = Σ_i w_i * exp(-(max(0, t - t_i))/τ)`

⸻

## Usage
	1.	Select predictors (checkboxes)
	2.	Choose start year (dropdown)
	3.	Pick aerosol dataset (CREST or GISS fallback)
	4.	Adjust lag ranges or enable auto lag search
	5.	Optionally include volcanic regressor and shading
	6.	Click Find Best Model to run analysis

Tabs include Best Model summary, Top 10 models, Trends, Diagnostics, Data coverage, Selected lags, and CCF preview.

⸻

## Attribution
	•	Data: NASA GISTEMP, NOAA GML/PSL, SATIRE (MPS), CREST, GISS Sato
	•	Volcanic eruption list: compiled from peer-reviewed climate literature

Code released under GNU General Public License v3.0.
