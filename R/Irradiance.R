#' Photosynthetic Photon Flux Density PPFD (PPFD) measurements of vine dataset
#'
#' Data were collected during an experiment conducted on a vineyard
#' of the INRAE/Institut Agro campus at Montpellier in 2014 (Syrah vines). The objective
#' of the experiment was to study the influence of the micro-climate (temperature and
#' irradiance) at the grape level on the anthocyanin contents of the berries indicated by
#' the Ferari index. This dataset is related to irradiance measurements in
#' the morning (sunrise to twelve am) between July 24th, 2014 at 09:00 am and August 01,
#' 2014 at 09:00 am. These observations are made at the same time (every 12 minutes) as the
#' temperature observations. The individuals are in columns while the observation times
#' are in rows. The same
#' individuals are also present in the Temperature and FerariIndex_Difference datasets.
#'
#' @format A data frame (of one functionnal variable) with 127 rows (observation times)
#'  and 33 columns: the 1st one is a character vector which corresponds to date-time 
#'  in format "yyyy-mm-dd hh:mm:ss", the others are numeric vectors made of the 
#'  observations of irradiance (PPFD) measured in \eqn{10^{-6}mol.m^{-2}.s^{-1}} on each of 
#'  the 32 statistical individuals Indiv1,...,Indiv32.  Irradiance corresponds to the number of incident 
#'  photons useful for photosynthesis, received per unit of time on a horizontal surface unit. 
#'  
#'
#' @source These data were acquired during the Innovine project, funded by the Seventh
#' Framework Programme of the European Community (FP7/2007-2013), under Grant Agreement
#' No. FP7-311775.

"Irradiance"
