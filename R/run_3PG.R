#' @title Runs a model simulation
#'
#' @description Runs the model, with options for the model setup, including the interception approach (3-PG vs 3-PGmix), soil model and disturbances.
#'
#' @param site table as described in \code{\link{prepare_input}} containing the information about site conditions.
#' @param species table as described in \code{\link{prepare_input}} containing the information about species level data. Each row corresponds to one species/cohort.
#' @param climate  table as described in \code{\link{prepare_input}} containing the information about monthly values for climatic data. See also \code{\link{prepare_climate}}
#' @param thinning table as described in \code{\link{prepare_input}} containing the information about thinnings. See also \code{\link{prepare_thinning}}
#' @param parameters table as described in \code{\link{prepare_input}} containing the information about parameters to be modified. See also \code{\link{prepare_parameters}}
#' @param size_dist table as described in \code{\link{prepare_input}} containing the information about size distributions. See also \code{\link{prepare_sizeDist}}
#' @param settings a list as described in \code{\link{prepare_input}} with settings for the model.
#' @param check_input \code{logical} if the input shall be checked for consistency. It will call \code{\link{prepare_input}} function.
#' @param df_out \code{logical} if the output shall be long data.frame (TRUE) the 4-dimensional array (FALSE).
#'
#' @details `i3PGmiX` enhanced 3-PGmix model, based on the r3PG package (https://github.com/trotsiuk/r3PG).i3PGmiX advances the representation of several processes in 3-PGmix, including improved GPP and NPP computation, as well as modules for soil and disturbances.
#'
#' @note The \code{run_3PG} also checks the quality of input data. When names, or structures are not consistent with requirements it will return an error. Turn this off to optimize for speed.
#'
#' @return either a 4-dimentional array or a data.frame, depending on the parameter \code{df_out}. More details on the output is \code{\link{i_output}}
#'
#' @seealso \code{\link{prepare_input}}, \code{\link{prepare_parameters}}, \code{\link{prepare_sizeDist}}, \code{\link{prepare_thinning}}, \code{\link{prepare_climate}}
#'
#'
#' @references
#' Forrester, D. I., 2020. 3-PG User Manual. Swiss Federal Institute for Forest, Snow and Landscape Research WSL, Birmensdorf, Switzerland. 70 p. Available at the following web site: \url{http://sites.google.com/site/davidforresterssite/home/projects/3PGmix/3pgmixdownload}
#'
#' Forrester, D. I., & Tang, X. (2016). Analysing the spatial and temporal dynamics of species interactions in mixed-species forests and the effects of stand density using the 3-PG model. Ecological Modelling, 319, 233–254. \doi{10.1016/j.ecolmodel.2015.07.010}
#'
#' Landsberg, J. J., & Waring, R. H., 1997. A generalised model of forest productivity using simplified concepts of radiation-use efficiency, carbon balance and partitioning. Forest Ecology and Management, 95(3), 209–228. \doi{10.1016/S0378-1127(97)00026-1}
#'
#' Sands, P. J., 2010. 3PGpjs user manual. Available at the following web site: \url{http://3pg.sites.olt.ubc.ca/files/2014/04/3PGpjs_UserManual.pdf}
#'
#' Trotsiuk, V., Hartig, F., & Forrester, D. I. (2020). r3PG–an r package for simulating forest growth using the 3‐PG process‐based model. Methods in Ecology and Evolution, 11(11), 1470-1475.
#'
#' @export
#'
#' @useDynLib i3PGmiX
#'
run_3PG <- function(
  site,
  species,
  climate,
  thinning = NULL,
  parameters = NULL,
  size_dist = NULL,
  settings = NULL,
  check_input = TRUE,
  df_out = TRUE
){

  thinn_null <- is.null(thinning)

  # Check and prepare input if required
  if( check_input ){

    input_checked = prepare_input(site = site, species = species, climate = climate,
      thinning = thinning, parameters = parameters, size_dist = size_dist,
      settings = settings)

    # extract output from the list
    site = input_checked$site
    species = input_checked$species
    climate = input_checked$climate
    thinning = input_checked$thinning
    parameters = input_checked$parameters
    size_dist = input_checked$size_dist
    settings = input_checked$settings

  }

  # Make small adjustments to the input and tranform it to matrix

  # site
  from = as.Date(paste(site$from,"-01",sep=""))
  site$year_i = as.numeric(format(from,'%Y'))
  site$month_i = as.numeric(format(from,'%m'))
  site = site[,c('latitude', 'altitude', 'soil_class', 'asw_i', 'asw_min', 'asw_max', 'year_i', 'month_i','soil_carbon','soil_nitrogen','soil_phosphorus','bb_cbp','bb_cst','background_infestation')]
  site = as.matrix( site, nrow = 1, ncol = 14)

  # species
  n_sp = dim( species )[1]
  sp_names = species$species
  planted = as.Date(paste(species$planted,"-01",sep=""))
  species$year_p = as.numeric(format(planted,'%Y'))
  species$month_p = as.numeric(format(planted,'%m'))
  species = species[,c('year_p', 'month_p', 'fertility', 'stems_n', 'biom_stem', 'biom_root', 'biom_foliage','deadwood_carbon', 'deadwood_nitrogen', 'litter_carbon', 'litter_nitrogen','n_sapling','sapling_wf','sapling_wr','sapling_ws','cwd_removal','salvage_wind','salvage_biotic','salvage_fire')]
  species = as.matrix( species, nrow = n_sp, ncol = 19)

  # climate
  n_m = dim(climate)[1]
  climate = climate[,c('tmp_min', 'tmp_max', 'tmp_ave', 'prcp', 'rh', 'srad', 'frost_days', 'vpd_day', 'co2', 'd13catm','n_depo','weibull_a','weibull_k')]
  climate = as.matrix( climate, nrow = n_m, ncol = 13)

  # thinning
  n_man = dim(thinning)[1]
  if( thinn_null ){
    t_t = 1L
    }else{
      if( dim(thinning)[3] == 1 ){
        t_t = as.integer( length(thinning[,1,]))
      }else{
        t_t = colSums( matrix( !is.na(thinning[,1,]), ncol = dim(thinning)[3]) )
        t_t = as.integer(t_t)
      }
    }

  # Parameters
  parameters = as.matrix( parameters[,-1], nrow = 181, ncol = n_sp)

  # Size distribution
  size_dist = as.matrix( size_dist[,-1], nrow = 30, ncol = n_sp)

  # Settings
  settings <- as.integer( unlist(settings) )

  # Run the simulations
  r3PG_out = .Call('s_3PG_c',
    siteInputs = site,
    speciesInputs = species,
    forcingInputs = climate,
    managementInputs = thinning,
    parameterInputs = parameters,
    sizeDistInputs = size_dist,
    n_sp = n_sp,
    n_m = n_m,
    n_man = n_man,
    t_t = t_t,
    settings = settings)


  # Remove the values for dead cohort or not recruited cohort
  remove_sim_arr <- array(r3PG_out[,,2,2] <= 0 | r3PG_out[,,2,1] < 0, dim=dim(r3PG_out))
  r3PG_out[remove_sim_arr] <- NA_real_


  if( df_out ){
    r3PG_out = transf.out( sim = r3PG_out, sp_names = sp_names, year_i = site[7], month_i = site[8] )
  }

  return( r3PG_out )

}

.onUnload <- function(libpath) {
  library.dynam.unload("i3PGmiX", libpath)
}


transf.out <- function( sim, sp_names, year_i, month_i ){

  # internal variables
  n_ob = dim(sim)[1]
  n_sp = dim(sim)[2]

  sim <- as.data.frame.table( sim, stringsAsFactors = F, responseName = 'value')

  if( month_i == 12){
    year_i = year_i + 1
    month_i = 0
  }

  sim$date <- seq( from = as.Date( paste(year_i, month_i+1, 01, sep = '-') ), by = "month", length.out = n_ob) - 1
  sim$species <- rep(sp_names, each = n_ob)
  sim$group <- rep( unique(i_output$variable_group), each = n_ob * n_sp)
  sim$variable <- rep( i_output$variable_name[order(i_output$variable_id)], each = n_ob * n_sp)

  sim <- sim[!is.na(sim$value),]

  sim <- sim[,c('date', 'species', 'group', 'variable', 'value')]

  return(sim)

}

