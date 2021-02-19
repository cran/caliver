## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, eval = FALSE, fig.width = 7)

## -----------------------------------------------------------------------------
#  install.packages("caliver")

## -----------------------------------------------------------------------------
#  library("caliver")

## -----------------------------------------------------------------------------
#  library("sp")
#  
#  # Get all the GFED4 basis regions
#  BasisRegions <- readRDS(system.file("extdata", "GFED4_BasisRegions.rds",
#                                      package = "caliver"))
#  plot(BasisRegions)

## -----------------------------------------------------------------------------
#  Europe <- BasisRegions[BasisRegions$Region == "EURO",]
#  plot(Europe)

## -----------------------------------------------------------------------------
#  library("raster")
#  # Italy
#  Italy <- raster::getData(name = "GADM", country = "Italy", level = 0)

## -----------------------------------------------------------------------------
#  # Assuming the FWI reanalysis (GEFF reanalysis v3.0, 7.5GB) dataset is downloaded
#  fwi <- raster::brick("~/repos/caliver_extras/data/fwi.nc")

## -----------------------------------------------------------------------------
#  # Global Land Cover 2000
#  # Documentation
#  # https://forobs.jrc.ec.europa.eu/data/products/glc2000/GLC2000_EUR_20849EN.pdf
#  # Get map
#  # system("wget https://forobs.jrc.ec.europa.eu/data/products/glc2000/glc2000_v1_1_Grid.zip")
#  # Unzip the folder in the working directory
#  
#  # Load map
#  fuel <- raster::raster("~/repos/caliver_extras/data/Grid/glc2000_v1_1/w001001.adf")
#  # Save legend in data.frame
#  df <- as.data.frame(levels(fuel))
#  
#  # If fuel and fwi have different extent and/or number of rows/cols, you will need to resample
#  fuel <- raster::resample(fuel, fwi, method = "ngb", progress = "text")
#  
#  # Remove Bare Areas (19), Water Bodies (20), Snow and Ice (21), Artificial surfaces and associated areas (22), No data (23).
#  codes_to_exclude <- 19:23
#  fuel[fuel %in% codes_to_exclude] <- NA
#  levels(fuel) <- df[1:18,]
#  
#  # Finaly, mask the first layer in fwi
#  fwi_masked <- raster::mask(fwi[[1]], fuel)
#  plot(fwi_masked, col = terrain.colors(120))

## -----------------------------------------------------------------------------
#  # Get indices for climatological reference period 1981-2010
#  idx <- which(substr(x = names(fwi), start = 2, stop = 5) %in% 1981:2010)
#  # Mask, crop and subset
#  fwi_Italy <- mask_crop_subset(r = fwi, p = Italy, idx = idx, progress = "text")

## -----------------------------------------------------------------------------
#  # Calculate the daily climatology for the period 10-20 August
#  clima_Italy <- daily_clima(r = fwi_Italy,
#                             dates = seq.Date(from = as.Date("2020-08-10"),
#                                              to = as.Date("2020-08-20"),
#                                              by = "day"))

## -----------------------------------------------------------------------------
#  mapsA <- get_percentile_map(r = fwi_Italy, probs = c(0.50, 0.75, 0.90))

## -----------------------------------------------------------------------------
#  mapsB <- get_percentile_map(r = clima_Italy, probs = c(0.50, 0.75, 0.90))

## ---- fig.width = 7-----------------------------------------------------------
#  # Use the raster plot method
#  raster::plot(mapsA, main = c("FWI 50th perc.", "FWI 75th perc.", "FWI 90th perc."))
#  
#  # Use the caliver plot_percentile_map function
#  plot_percentile_map(maps = mapsA,
#                      main = c("FWI 50th perc.", "FWI 75th perc.", "FWI 90th perc."))

## -----------------------------------------------------------------------------
#  fwi_class <- classify_index(fwi_masked, index = "fwi", thresholds = NULL, labels = NULL)
#  rasterVis::levelplot(fwi_class, col.regions = caliver:::effis_palette, att = "Class")

## -----------------------------------------------------------------------------
#  dataDates <- as.Date(substr(names(fwi), 2, 11), format = "%Y.%m.%d")
#  # Define a function to extract fire seasons in Europe
#  seasons <- get_fire_season(dates = dataDates, zone = "north")
#  
#  # Create an index of fire season dates
#  fireSeasonIndex <- which(seasons == TRUE)

## -----------------------------------------------------------------------------
#  # Mask/Crop/Subset FWI over Europe
#  fwi_euro <- mask_crop_subset(r = fwi, p = Europe, idx = fireSeasonIndex)
#  
#  # Calculate homogenised fire danger levels (or thresholds) for Europe
#  EuropeThr <- get_fire_danger_levels(fire_index = fwi_euro)
#  EuropeThr

## -----------------------------------------------------------------------------
#  # Country level: use a loop to calculate levels for three sample countries.
#  # Please note that in the original paper many more countries were used, here we process
#  # only three to reduce the processing time but the procedure is exactly the same.
#  EUcountries <- c("Spain", "GBR", "Italy")
#  
#  for (country in EUcountries){
#  
#    print(country)
#    country_poly <- raster::getData(name = "GADM",
#                                    country = country,
#                                    level = 0)
#  
#    # Mask/Crop/Subset FWI and generate thresholds for country
#    country_fwi <- mask_crop_subset(r = fwi_euro, p = country_poly)
#    country_thrs <- get_fire_danger_levels(fire_index = country_fwi)
#  
#    # Append values to data.frame
#    if (country == "Spain") {
#      df <- data.frame(matrix(country_thrs, nrow = 1))
#    }else{
#      df <- rbind(df, country_thrs)
#    }
#  
#    print(df)
#  
#  }
#  
#  df_thr <- data.frame(cbind(EUcountries, df, stringsAsFactors=FALSE))
#  names(df_thr) <- c("Country", "Low", "Moderate", "High", "VeryHigh", "Extreme")

## -----------------------------------------------------------------------------
#  countryPDF <- plot_pdf(fire_index = country_fwi,
#                         thresholds = country_thrs,
#                         upper_limit = 60)

## -----------------------------------------------------------------------------
#  library("pROC")
#  
#  # Assuming burned areas was downloaded and combined in a grid file,
#  # this can be loaded as a brick:
#  BurnedAreas <- raster::brick("GFED4_BurnedAreas/BurnedArea.grd")
#  
#  # Mask and crop burned areas over Europe
#  BA <- mask_crop_subset(r = BurnedAreas, p = Europe)
#  
#  # If observations layers have no date, assign it!
#  dataDates <- seq.Date(from = as.Date("2003-01-01"),
#                        to = as.Date("2015-12-31"), by = "day")
#  names(BA) <- dataDates
#  
#  EuroThrHigh <- as.numeric(df_thr[df_thr$Country == "Europe", 4])
#  
#  # The above can be saved and re-loaded as follows:
#  raster::writeRaster(BA, filename="BurnedAreaEurope.grd",
#                      bandorder='BIL', overwrite=TRUE, progress = 'text')
#  BurnedAreaEurope <- raster::brick("BurnedAreaEurope.grd")
#  
#  # For the validation we do not want to subset over the fire season, subset to match days in BurnedAreaEurope
#  FWIEURO <- mask_crop_subset(r = FWI, p = Europe, idx = which(names(FWI) %in% names(BurnedAreaEurope)))
#  # The above can be saved and re-loaded as follows:
#  raster::writeRaster(FWIEURO, filename="FWIEURO.grd", bandorder='BIL', overwrite=TRUE, progress = 'text')
#  FWIEURO <- raster::brick("FWIEURO.grd")
#  
#  # Contingency table for JRC - Europe as a whole
#  x1 <- validate_fire_danger_levels(fire_index = FWIEURO, observation = BurnedAreaEurope,
#                                   fire_threshold = 21.3, obs_threshold = 50)
#  tab_x <- table(pred = x1$pred, obs = x1$obs)
#  hits <- tab_x[2,2]
#  misses <- tab_x[1,2]
#  correct_negatives <- tab_x[1,1]
#  false_alarms <- tab_x[2,1]
#  # POD 47%
#  round(hits/(hits+misses),2)*100
#  roc1 <- pROC::roc(response = x1$obs, predictor = x1$pred)
#  pROC::plot.roc(roc1, print.auc = pROC::auc(roc1), print.auc.x = 0, print.auc.y = 0.9)
#  
#  # Contingency table for caliver - Europe as a whole
#  x2 <- validate_fire_danger_levels(fire_index = FWIEURO, observation = BurnedAreaEurope,
#                                   fire_threshold = EuroThrHigh, obs_threshold = 50)
#  tab_x <- table(pred = x2$pred, obs = x2$obs)
#  hits <- tab_x[2,2]
#  misses <- tab_x[1,2]
#  # POD 65%
#  round(hits/(hits+misses),2)*100
#  roc2 <- pROC::roc(response = x2$obs, predictor = x2$pred)
#  pROC::plot.roc(roc2, col = "red", add = TRUE,
#                 print.auc = pROC::auc(roc2), print.auc.x = 0, print.auc.y = 0.95,
#                 print.auc.col = "red")
#  
#  # Loop throught the countries
#  for (country in df_thr[,"Country"]){
#  
#    print(country)
#  
#    countryPoly <- raster::getData(name = "GADM", country = country, level = 0)
#    countryThr <- as.numeric(df_thr[df_thr$Country == country, 4])
#  
#    # Crop RasterBricks over country of interest
#    BA_country <- mask_crop_subset(r = BurnedAreaEurope, p = countryPoly)
#    FWI_country <- mask_crop_subset(r = FWIEURO, p = countryPoly)
#  
#    JRC <- validate_fire_danger_levels(fire_index = FWI_country,
#                                       observation = BA_country,
#                                       fire_threshold = 21.3,
#                                       obs_threshold = 50)
#    tab_JRC <- data.frame(table(JRC$pred, JRC$obs))
#    caliver1 <- validate_fire_danger_levels(fire_index = FWI_country,
#                                            observation = BA_country,
#                                            fire_threshold = EuroThrHigh,
#                                            obs_threshold = 50)
#    tab_caliver1 <- data.frame(table(caliver1$pred, caliver1$obs))
#    caliver2 <- validate_fire_danger_levels(fire_index = FWI_country,
#                                            observation = BA_country,
#                                            fire_threshold = countryThr,
#                                            obs_threshold = 50)
#    tab_caliver2 <- data.frame(table(caliver2$pred, caliver2$obs))
#  
#    if (country == "Spain") {
#      df_caliver1 <- df_caliver2 <- df_effis <- data.frame("pred" = tab_caliver1$pred, "obs" = tab_caliver1$obs)
#      i <- 3
#    }
#  
#    df_caliver1 <- cbind(df_caliver1, tab_caliver1$Freq)
#    names(df_caliver1)[i] <- country
#    df_caliver2 <- cbind(df_caliver2, tab_caliver2$Freq)
#    names(df_caliver2)[i] <- country
#    df_effis <- cbind(df_effis, tab_JRC$Freq)
#    names(df_effis)[i] <- country
#    i <- i + 1
#  
#    rm(countryPoly, countryThr, BA_country, FWI_country)
#  
#  }
#  
#  # Save contingency tables
#  saveRDS(df_caliver1, "df_caliver1.rds")
#  saveRDS(df_caliver2, "df_caliver2.rds")
#  saveRDS(df_effis, "df_effis.rds")
#  
#  # Europe (EFFIS danger levels)
#  sum(df_effis[4,3:27]) # hits
#  sum(df_effis[3,3:27]) # misses
#  # Europe (averaged danger levels)
#  sum(df_caliver1[4,3:27]) # hits
#  sum(df_caliver1[3,3:27]) # misses
#  # Europe (country-specific danger levels)
#  sum(df_caliver2[4,3:27]) # hits
#  sum(df_caliver2[3,3:27]) # misses
#  
#  # UK (EFFIS danger levels)
#  df_effis[4, which(names(df_caliver2) == "United Kingdom")] # hits
#  df_effis[3, which(names(df_caliver2) == "United Kingdom")] # misses
#  # UK (EU averaged danger levels)
#  df_caliver1[4, which(names(df_caliver2) == "United Kingdom")] # hits
#  df_caliver1[3, which(names(df_caliver2) == "United Kingdom")] # misses
#  # UK (country-specific danger levels)
#  df_caliver2[4, which(names(df_caliver2) == "United Kingdom")] # hits
#  df_caliver2[3, which(names(df_caliver2) == "United Kingdom")] # misses
#  
#  # Spain (EFFIS danger levels)
#  df_effis[4, which(names(df_caliver2) == "Spain")] # hits
#  df_effis[3, which(names(df_caliver2) == "Spain")] # misses
#  # Spain (EU averaged danger levels)
#  df_caliver1[4, which(names(df_caliver2) == "Spain")] # hits
#  df_caliver1[3, which(names(df_caliver2) == "Spain")] # misses
#  # Spain (country-specific danger levels)
#  df_caliver2[4, which(names(df_caliver2) == "Spain")] # hits
#  df_caliver2[3, which(names(df_caliver2) == "Spain")] # misses
#  
#  # Italy (EFFIS danger levels)
#  df_effis[4, which(names(df_caliver2) == "Italy")] # hits
#  df_effis[3, which(names(df_caliver2) == "Italy")] # misses
#  # Italy (EU averaged danger levels)
#  df_caliver1[4, which(names(df_caliver2) == "Italy")] # hits
#  df_caliver1[3, which(names(df_caliver2) == "Italy")] # misses
#  # Italy (country-specific danger levels)
#  df_caliver2[4, which(names(df_caliver2) == "Italy")] # hits
#  df_caliver2[3, which(names(df_caliver2) == "Italy")] # misses

