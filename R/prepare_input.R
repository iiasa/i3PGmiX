#' @title Check and prepare input for running 3-PG model
#' @description Checks and prepares all input tables to be used in \code{\link{run_3PG}}. For detailed descriptions see Forrester (2020).
#'
#' @param site table containing the information about site conditions.
#' \itemize{
#' \item latitude: site latitude in the WGS84 coordinate system.
#' \item altitude: site altitude, m a.s.l.
#' \item soil_class: 1 - Sandy; 2 - Sandy loam; 3 - Clay loam; 4 - Clay; 0 - No effect of asw on production.
#' \item asw_i: initial available soil water (mm).
#' \item asw_min: minimum available soil water (mm).
#' \item asw_max: maximum available soil water (mm).
#' \item from: year and month indicating the start of simulation. Provided in form of year-month. E.g. "2000-01".
#' \item to: year and month indicating the end of simulation. Provided in form of year-month. E.g. "2009-12", will include December 2009 as last simulation month
#' \item soil_carbon: soil carbon stock (tC/ha)
#' \item soil_nitrogen: soil nitrogen stock (tN/ha)
#' \item soil_phosphorus: soil phosphorus stock (tP/ha)
#' \item bb_cbp: spatial scaling for bark beetle infestation
#' \item bb_cst: temporal scaling for bark beetle infestation
#' \item background_infestation: bark beetle background infestation (default = 1)
#' }
#' @param species table containing the information about species level data. Each row corresponds to one species/cohort.
#' \itemize{
#' \item species: species or cohort id/name. It must be consistent with species names in \code{thinning}, \code{parameters} and \code{sizeDist} tables.
#' \item planted: year and month indicating when species was planted. Provided in form of year-month. E.g. "2000-01".
#' \item fertility: soil fertility for a given species. Range from 0 to 1.
#' \item stems_n: number of trees per ha.
#' \item biom_stem: stem biomass for a given species (Mg/ha).
#' \item biom_root: root biomass for a given species (Mg/ha).
#' \item biom_foliage: initial foliage biomass (Mg/ha). If this is a leafless period, provide the spring foliage biomass.
#' \item deadwood_carbon: carbon stock in deadwood (Mg/ha)
#' \item deadwood_nitrogen: nitrogen stock in deadwood (Mg/ha)
#' \item litter_carbon: carbon stock in litter (Mg/ha)
#' \item litter_nitrogen: nitrogen stock in litter (Mg/ha)
#' \item n_sapling: number of saplings in regenerated stand (N/ha)
#' \item sapling_wf: foliage biomass of saplings (Mg/ha)
#' \item sapling_wr: root biomass of saplings (Mg/ha)
#' \item sapling_ws: stem biomass of saplings (Mg/ha)
#' \item cwd_removal: share of removed deadwood
#' \item salvage_wind: share of salvaged stock after wind disturbances
#' \item salvage_biotic: share of salvaged stock after biotic disturbances
#' \item salvage_fire: share of salvaged stock after wildfires
#' }
#' @param climate  table containing the information about monthly values for climatic data. If the climate table has exactly 12 rows it will be replicated for the number of years and months specified by \code{from} - \code{to}. Otherwise, it will be subsetted to the selected time period. More details about preparing climate data are at \code{\link{prepare_climate}}.
#' \itemize{
#' \item year: year of observation (only required for subsetting) (optional).
#' \item month: months of observation (only required for subsetting) (optional).
#' \item tmp_min: monthly mean daily minimum temperature (C).
#' \item tmp_max: monthly mean daily maximum temperature (C).
#' \item tmp_ave: monthly mean daily average temperature (C) (optional).
#' \item prcp: monthly rainfall (mm month-1).
#' \item rh: monthly average relative humidity (%).
#' \item srad: monthly mean daily solar radiation (MJ m-2 d-1).
#' \item frost_days: frost days per month (d month-1).
#' \item vpd_day: frost days per month (mbar) (optional).
#' \item co2: monthly mean atmospheric co2 (ppm), required if calculate_d13c=1 (optional)
#' \item d13catm: monthly mean isotopic composition of air (‰), required if calculate_d13c=1 (optional)
#' \item n_depo: monthly nitrogen deposition (tN/ha) (optional)
#' \item weibull_a: wind speed weibull a parameter, required if wind_dist=1 (optional)
#' \item weibull_k: wind speed weibull k parameter, required if wind_dist=1 (optional)
#' }
#' @param thinning table containing the information about thinnings. If there is no thinning, it must be \code{NULL}.
#' \itemize{
#' \item species: species or cohort id/name. It must be consistent with species names in \code{species}, \code{parameters} and \code{sizeDist} tables.
#' \item age: age when thinning is performed.
#' \item stems_n: number of trees remaining after thinning
#' \item foliage: type of thinning (above/below). Default is 1.
#' \item root: type of thinning (above/below). Default is 1.
#' \item stem: type of thinning (above/below). Default is 1.
#' \item BA: share of basal area removal or target basal area after thinning.
#' \item Vol: share of volume removal or target volume after thinning.
#' \item crop_trees: number of crop trees per hectare.
#' }
#' @param parameters table containing the information about parameters to be modified. Values that are not provided are replaced by defaults.
#' \itemize{
#' \item parameter: name of the parameter, must be consistent in naming with \code{\link{i_parameters}}
#' \item species: each column must correspond to species/cohort id/name, as defined in \code{species} table
#' }
#' @param size_dist table containing the information about size distribution to be modified. Values that are not provided are replaced by defaults.
#' \itemize{
#' \item parameter: name of the parameter, must be consistent in naming with \code{\link{i_sizeDist}}
#' \item species: each column must correspond to species/cohort id/name, as defined in \code{species} table
#' }
#' @param settings a list with settings for the model. Values that are not provided are replaced by defaults.
#' \itemize{
#' \item light_model: `1` - 3-PGpjs (default); `2` - 3-PGmix
#' \item transp_model: `1` - 3-PGpjs (default); `2` - 3-PGmix
#' \item phys_model:  `1` - 3-PGpjs (default); `2` - 3-PGmix
#' \item height_model: `1` - linear (default); `2` - non-linear
#' \item correct_bias: `0` - no (default); `1` - yes
#' \item calculate_d13c: `0` - no (default); `1` - yes
#' \item soil_model: `0` - no soil model (default); `1` - ICBM/2N; `2` - Yasso
#' \item gpp_model: `1` - 3-PGpjs (default); `2` - P-model; `3` - P-hydro model
#' \item method_jmaxlim:  `1` - Wang; `2` - Smith; `3` - none (default)
#' \item nut_lim: `0` - fixed fr (default); `2` - dynamic fr based on N/P demand and soil available N/P
#' \item maint_resp: `1` - fixed NPP/GPP ratio (default); `2` - dynamic maintenance respiration
#' \item canopy_cond: `1` - Jarvis (default); `2` - P-model; `3` - P-hydro model
#' \item mortality_model: `1` - 3-PGpjs (default); `2` - Brandl S(t); `3` - LPJ-GUESS; `4` - iLand
#' \item wind_dist: `0` - no wind disturbance (default); `1` - ForestGALES TMC
#' \item beetle_dist:  `0` no bark beetle disturbance (default); `2` - LandClim bark beetle disturbance
#' \item d_sus_type :  `1` bark beetle vitality effect based on f_sw (default); `2` - bark beetle vitality effect based on growth efficiency
#' \item dist_start :  month when disturbance are initially activated
#' \item spruce_idx :  index of spruce in the species list
#' }
#'
#' @details This function checks and prepares the input data for the \code{\link{run_3PG}}. The output is a list with 7 tables. Each of them corresponds to the one from input.
#'
#' @seealso \code{\link{run_3PG}}, \code{\link{prepare_parameters}}, \code{\link{prepare_sizeDist}}, \code{\link{prepare_thinning}}, \code{\link{prepare_climate}}
#'
#' @return a list with seven tables. Each table corresponds to one of the input tables.
#'
#' @example inst/examples/prepare_input-help.R
#'
#' @references
#' Forrester, D. I., 2020. 3-PG User Manual. Swiss Federal Institute for Forest, Snow and Landscape Research WSL, Birmensdorf, Switzerland. 70 p. Available at the following web site: \url{http://sites.google.com/site/davidforresterssite/home/projects/3PGmix/3pgmixdownload}
#'
#'Sands, P. J., 2010. 3PGpjs user manual. Available at the following web site: \url{http://3pg.sites.olt.ubc.ca/files/2014/04/3PGpjs_UserManual.pdf}
#'
#' @export
#'
prepare_input <- function(
  site,
  species,
  climate,
  thinning = NULL,
  parameters = NULL,
  size_dist = NULL,
  settings = NULL
){

  # Site
  if( !identical( c("latitude","altitude","soil_class","asw_i","asw_min","asw_max","from","to",
                    "soil_carbon","soil_nitrogen","soil_phosphorus"), colnames(site)[1:11]) ){
    stop( 'Columns names of the site table must correspond to: latitude, altitude, soil_class,
          asw_i, asw_min, asw_max, from, to, soil_carbon, soil_nitrogen, soil_phosphorus')
  }

  # Species
  if( !identical(c("species","planted","fertility","stems_n","biom_stem","biom_root","biom_foliage",
                   "deadwood_carbon", "deadwood_nitrogen", "litter_carbon", "litter_nitrogen","n_sapling",
                   "sapling_wf","sapling_wr","sapling_ws","cwd_removal","salvage_wind","salvage_biotic",
                   "salvage_fire"), colnames(species)) ){
    stop( 'Columns names of the species table must correspond to: species, planted, fertility, stems_n,
          biom_stem, biom_root, biom_foliage, deadwood_carbon, deadwood_nitrogen, litter_carbon, litter_nitrogen',
          'n_sapling','sapling_wf','sapling_wr','sapling_ws','cwd_removal','salvage_wind','salvage_biotic','salvage_fire')
  }

  # Settings
  set_def = list(light_model = 1, transp_model = 1, phys_model = 1, height_model = 1, correct_bias = 0,
                 calculate_d13c = 0, soil_model = 0, gpp_model = 1, method_jmaxlim = 3, nut_lim=0, maint_resp = 1,
                 canopy_cond = 1, mortality_model = 1 , wind_dist = 0, beetle_dist = 0, d_sus_type = 1,
                 dist_start = 1, spruce_idx = -999)

  set_def[names(settings)] <- settings


  if (! "bb_cbp" %in% colnames(site) | ! "bb_cst" %in% colnames(site)){
    if (! "bb_cbp" %in% colnames(site)) site$bb_cbp <- 1
    if (! "bb_cst" %in% colnames(site)) site$bb_cst <- 1
    if (set_def$beetle_dist == 1){
      print("!!! missing scaling factors for bark beetle module, using default values (=1)")
    }
  }
  if (! "background_infestation" %in% colnames(site)){
    if (! "background_infestation" %in% colnames(site)) site$background_infestation <- 0.01
    if (set_def$beetle_dist == 1){
      print("!!! missing background infestation for bark beetle module, using default values (=0.01)")
    }
  }
  # Correct ordering of bark beetle data
  if (set_def$beetle_dist == 1 & !identical(c("bb_cbp","bb_cst","background_infestation"),colnames(site)[12:14])){
    site_aux <- site[,1:11]
    site_aux$bb_cbp <- site$bb_cbp
    site_aux$bb_cst <- site$bb_cst
    site_aux$background_infestation <- site$background_infestation
    site <- site_aux
  }

  # Define spruce index in the species list
  if (set_def$beetle_dist == 1){
    s_idx <- grep(".*Picea.*abies.*",species$species,ignore.case=T)
    n_spr <- length(s_idx)
    set_def$spruce_idx <- ifelse(n_spr>0,s_idx,-999)
    if (n_spr == 0) print("!!! bark beetle module is activated but spruce was not found in the species list")
  }


  # Climate
  if( set_def['calculate_d13c'] == 1 ){
    if( !all( c("co2","d13catm") %in% colnames(climate) ) ){
      stop('Please provide forcing data for co2 and d13catm in climate, if calculate_d13c = 1')
    }
  }

  climate = prepare_climate(climate = climate, from = site$from, to = site$to)

  # Thinning
  thinning = prepare_thinning( thinning = thinning, sp_names = species$species)

  # Parameters
  parameters = prepare_parameters( parameters = parameters, sp_names = species$species)

  # Size distribution
  if( set_def['correct_bias'] == 1 & is.null(size_dist) ){
    stop('Please provide sisze_dist table or change the setting to size_dist = 0')
  }
  size_dist = prepare_sizeDist( size_dist = size_dist, sp_names = species$species)


  # return the checked output
  out <- list( site = site, species = species, climate = climate, thinning = thinning, parameters = parameters, size_dist = size_dist, settings = set_def)

  return( out )
}
