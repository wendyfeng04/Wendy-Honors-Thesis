
# =============================================================================\
# US-Fmf FLUXMET HH: Daily and Annual ET, PET, P Summary
# Site: Flagstaff Management Forest, Arizona
# Data: AMF_US-Fmf_FLUXNET_FLUXMET_HH_2005-2010_v1.3_r1.csv
#
# INPUT:
#   Half-hourly FLUXMET data
#
# OUTPUT TABLES:
#   1) daily_summary:
#        daily ET, PET, P
#        cumulative ET, PET, P within each year
#        rain_day: rained = 1, not rained = 0
#        daily missing percentage for ET, P, PET-met
#
#   2) annual_summary:
#        annual ET, PET, P
#        annual average daily ET, PET, P
#        annual rain frequency and mean event depth
#        annual missing-data percentages
#
# UNITS:
#   all water depths in m 
# =============================================================================

library(Evapotranspiration)
library(zoo)
library(dplyr)
library(lubridate)
library(tidyr)
library(plantecophys)

# =============================================================================
# 1. LOAD AND CLEAN
# =============================================================================

file_path <- "~/Desktop/Thesis Fluxnet Data/AMF_US-Fmf_FLUXNET_FLUXMET_HH_2005-2010_v1.3_r1.csv"

raw <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)

raw[raw == -9999] <- NA #missing data

raw$datetime_start <- as.POSIXct(
  as.character(raw$TIMESTAMP_START),
  format = "%Y%m%d%H%M",
  tz = "UTC"
)

raw$datetime_end <- as.POSIXct(
  as.character(raw$TIMESTAMP_END),
  format = "%Y%m%d%H%M",
  tz = "UTC"
)

raw$date  <- as.Date(raw$datetime_start)
raw$year  <- year(raw$date)
raw$month <- month(raw$date)
raw$doy   <- yday(raw$date)
raw$dt_seconds <- 1800

# =============================================================================
# 2. DAILY ET 
#
# Inputs from HH fluxnet data:
# LE_F_MDS: latent heat flux (W m⁻² = J s⁻¹ m⁻²)
# TA_F: air temperature (°C)
#
# ET conversion:
# Latent heat of vaporization: lambda_MJ_kg = 2.501 - 0.002361 * TA_F -> lambda_J_kg  = lambda_MJ_kg * 1e6
# ET_m_step = (LE_F_MDS * dt_seconds / lambda_J_kg) / 1000
# Units check: ((J s-1 m-2) * s / (J kg-1)) / 1000 = kg m-2 / 1000 = mm water / 1000 (water density = 1000 kg m⁻³) = m water
# =============================================================================
daily_ET <- raw %>%
  mutate(
    lambda_MJ_kg = 2.501 - 0.002361 * TA_F,
    lambda_J_kg  = lambda_MJ_kg * 1e6,
    ET_m_step = LE_F_MDS * dt_seconds / lambda_J_kg / 1000,
    ET_m_step = ifelse(!is.na(ET_m_step) & ET_m_step < 0, 0, ET_m_step),
    
    ET_missing = is.na(LE_F_MDS) | is.na(TA_F) | is.na(lambda_J_kg)
  ) %>%
  group_by(date, year, month, doy) %>%
  summarise(
    ET_m_day = sum(ET_m_step, na.rm = TRUE),
    ET_missing_pct_daily = round(mean(ET_missing) * 100, 1),
    mean_lambda_MJ_kg = round(mean(lambda_MJ_kg, na.rm = TRUE), 4),
    n_records_ET = n(),
    .groups = "drop"
  )
  

# =============================================================================
# 3. DAILY PRECIPITATION 
#
# Inputs from HH fluxnet data:
# P_F = precipitation (mm per timesetp)
#
# Daily P = sum half-hourly P_F / 1000 (m per day)
# =============================================================================

rain_threshold_m <- 0.0002 # treat extremely small detected daily rainfall values as noise
daily_P <- raw %>%
  mutate(
    P_m_step = P_F / 1000,
    P_missing = is.na(P_F)
  ) %>%
  group_by(date, year, month, doy) %>%
  summarise(
    P_m_day = sum(P_m_step, na.rm = TRUE),
    P_missing_pct_daily = round(mean(P_missing) * 100, 3),
    n_records_P = n(),
    .groups = "drop"
  ) %>%
  mutate(
    rain_day = ifelse(P_m_day > rain_threshold_m, 1L, 0L)
  )

# =============================================================================
# 4. DAILY PET
#
# R-Package Used: 
# https://cran.r-project.org/web/packages/Evapotranspiration/index.html
#
# Method Used:
# Priestley-Taylor
#
# Package Input Table:
# Date.daily
# Date.monthly
# J = day of year
# i = month
# ndays = days in month
# Tmax = daily maximum air temperature (°C):
# Tmin = daily minimum air temperature (°C):
# max & min half-hourly TA_F of the day
#
# RHmax = daily maximum relative humidity (%)
# RHmin = daily minimum relative humidity (%)
# use package plantecophys (https://cran.r-project.org/web/packages/plantecophys/plantecophys.pdf) to estimate HH Rh from HH temp and VPD (VPD_F,hPa), and extract daily min & max 
#
# Rs = daily solar radiation, MJ m⁻² day⁻¹
# =============================================================================
daily_met <- raw %>%
  mutate(
    SW_IN_MJ_step = SW_IN_F * dt_seconds / 1e6,
    
    # FLUXNET VPD_F is hPa, but plantecophys::VPDtoRH needs kPa
    VPD_kPa_step = VPD_F / 10,
    
    # Estimate half-hourly RH from half-hourly VPD and air temperature
    RH_step = VPDtoRH(VPD = VPD_kPa_step, TdegC = TA_F, Pa = 101),
    RH_step = pmin(100, pmax(0, RH_step))
  ) %>%
  group_by(date, year, month, doy) %>%
  summarise(
    Tmax = ifelse(all(is.na(TA_F)), NA_real_, max(TA_F, na.rm = TRUE)),
    Tmin = ifelse(all(is.na(TA_F)), NA_real_, min(TA_F, na.rm = TRUE)),
    
    RHmax = ifelse(all(is.na(RH_step)), NA_real_, max(RH_step, na.rm = TRUE)),
    RHmin = ifelse(all(is.na(RH_step)), NA_real_, min(RH_step, na.rm = TRUE)),
    
    Rs = ifelse(all(is.na(SW_IN_MJ_step)), NA_real_, sum(SW_IN_MJ_step, na.rm = TRUE)),
    
    TA_missing_pct_daily  = round(sum(is.na(TA_F)) / n() * 100, 1),
    VPD_missing_pct_daily = round(sum(is.na(VPD_F)) / n() * 100, 1),
    SW_missing_pct_daily  = round(sum(is.na(SW_IN_F)) / n() * 100, 1),
    
    .groups = "drop"
  ) %>%
  mutate(
    PET_met_missing = is.na(Tmax) |
      is.na(Tmin) |
      is.na(RHmax) |
      is.na(RHmin) |
      is.na(Rs) |
      Rs < 0
  )
head(daily_met)

daily_met_valid <- daily_met %>%
  filter(!PET_met_missing)

data_PET <- list(
  Date.daily   = zoo(daily_met_valid$date,                  order.by = daily_met_valid$date), # daily date
  Date.monthly = zoo(format(daily_met_valid$date, "%b-%Y"), order.by = daily_met_valid$date), # month-year label
  J            = zoo(daily_met_valid$doy,                   order.by = daily_met_valid$date), # day of year
  i            = zoo(daily_met_valid$month,                 order.by = daily_met_valid$date), # month number
  ndays        = zoo(days_in_month(daily_met_valid$date),   order.by = daily_met_valid$date), # days in that month
  
  Tmax         = zoo(daily_met_valid$Tmax,  order.by = daily_met_valid$date), # daily max air temperature, deg C
  Tmin         = zoo(daily_met_valid$Tmin,  order.by = daily_met_valid$date), # daily min air temperature, deg C
  RHmax        = zoo(daily_met_valid$RHmax, order.by = daily_met_valid$date), # daily max relative humidity, %
  RHmin        = zoo(daily_met_valid$RHmin, order.by = daily_met_valid$date), # daily min relative humidity, %
  Rs           = zoo(daily_met_valid$Rs,    order.by = daily_met_valid$date)  # daily solar radiation, MJ m-2 day-1
)

class(data_PET) <- "PriestleyTaylor"

constants_PT <- list(
  Elev    = 2180,              # site elevation above sea level, m
  lambda  = 2.45,              # latent heat of vaporization used by PET package, MJ kg-1
  lat_rad = 36.9216 * pi / 180,  # site latitude converted from degrees to radians
  Gsc     = 0.0820,            # solar constant, MJ m-2 min-1
  sigma   = 4.903e-9,          # Stefan-Boltzmann constant, MJ K-4 m-2 day-1
  alphaPT = 1.26,              # Priestley-Taylor coefficient
  G       = 0                 # soil heat flux, MJ m-2 day-1; usually 0 for daily timestep
)

pet_results <- ET.PriestleyTaylor(
  data            = data_PET,
  constants       = constants_PT,
  ts              = "daily",
  solar           = "data",    # use measured Rs directly
  alpha           = 0.23,      # albedo of evaporative surface
  message         = "yes",
  AdditionalStats = "yes",
  save.csv        = "no"
)

daily_PET <- data.frame(
  date = as.Date(index(pet_results$ET.Daily)),
  PET_m_day = as.numeric(pet_results$ET.Daily) / 1000
) %>%
  mutate(year = year(date))

daily_PET_missing <- daily_met %>%
  select(
    date, year,
    PET_met_missing,
    TA_missing_pct_daily,
    VPD_missing_pct_daily,
    SW_missing_pct_daily
  ) %>%
  mutate(
    PET_missing_pct_daily = ifelse(PET_met_missing, 100, 0)
  )

# =============================================================================
# 5. DAILY SUMMARY
# =============================================================================

daily_summary <- daily_ET %>%
  full_join(daily_P, by = c("date", "year", "month", "doy")) %>%
  full_join(daily_PET, by = c("date", "year")) %>%
  full_join(daily_PET_missing, by = c("date", "year")) %>%
  arrange(date) %>%
  group_by(year) %>%
  mutate(
    cumulative_ET_m  = cumsum(replace_na(ET_m_day,  0)),
    cumulative_PET_m = cumsum(replace_na(PET_m_day, 0)),
    cumulative_P_m   = cumsum(replace_na(P_m_day,   0))
  ) %>%
  ungroup() %>%
  select(
    date, year, month, doy,
    ET_m_day, PET_m_day, P_m_day,
    cumulative_ET_m, cumulative_PET_m, cumulative_P_m,
    rain_day,
    ET_missing_pct_daily,
    PET_missing_pct_daily,
    P_missing_pct_daily,
    TA_missing_pct_daily,
    VPD_missing_pct_daily,
    SW_missing_pct_daily,
    n_records_ET,
    n_records_P
  )

# =============================================================================
# 6. ANNUAL SUMMARY
# =============================================================================

annual_summary <- daily_summary %>%
  group_by(year) %>%
  summarise(
    ET_total_m  = sum(ET_m_day,  na.rm = TRUE),
    PET_total_m = sum(PET_m_day, na.rm = TRUE),
    P_total_m   = sum(P_m_day,   na.rm = TRUE),
    
    ET_avg_daily_m  = mean(ET_m_day,  na.rm = TRUE),
    PET_avg_daily_m = mean(PET_m_day, na.rm = TRUE),
    P_avg_daily_m   = mean(P_m_day,   na.rm = TRUE),
    
    rain_freq = mean(rain_day == 1, na.rm = TRUE),
    
    mean_event_depth_m = ifelse(
      sum(rain_day == 1, na.rm = TRUE) > 0,
      sum(P_m_day, na.rm = TRUE) / sum(rain_day == 1, na.rm = TRUE),
      NA_real_
    ),
    
    ET_missing_pct  = mean(ET_missing_pct_daily,  na.rm = TRUE),
    PET_missing_pct = mean(PET_missing_pct_daily, na.rm = TRUE),
    P_missing_pct   = mean(P_missing_pct_daily,   na.rm = TRUE),
    
    overall_missing_pct = mean(
      ET_missing_pct_daily > 0 |
        PET_missing_pct_daily > 0 |
        P_missing_pct_daily > 0,
      na.rm = TRUE
    ) * 100,
    
    .groups = "drop"
  ) %>%
  mutate(
    ET_over_PET = ET_total_m / PET_total_m,
    P_over_PET  = P_total_m  / PET_total_m
  ) %>%
  select(
    year,
    ET_total_m, PET_total_m, P_total_m,
    ET_over_PET, P_over_PET,
    ET_avg_daily_m, PET_avg_daily_m, P_avg_daily_m,
    rain_freq, mean_event_depth_m,
    ET_missing_pct, PET_missing_pct, P_missing_pct, overall_missing_pct
  )

# =============================================================================
# 7. EXPORT
# =============================================================================

output_dir <- "~/Desktop/Thesis Fluxnet Data/Processed Data"

write.csv(daily_summary,
          file.path(output_dir, "US_Fmf_daily.csv"),
          row.names = FALSE)

write.csv(annual_summary,
          file.path(output_dir, "US_Fmf_annual.csv"),
          row.names = FALSE)

