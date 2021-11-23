# National and regional load consumption data at 8 p.m, calendar and meteo variables related to it
#
#
# @format A data frame with 2409 rows and ... variables:
# \describe{
#   \item{Date}{date of observation}
#   \item{Year}{Year of observation}
#   \item{Month}{Month of observation}
#   \item{WeekDays}{Day of the week}
#   \item{WeekEnd}{Indicator of weekd-end}
#   \item{BH}{Indicator of Bank Hollidays}
#   \item{toy}{time of year, value grows linearly from 0 on the 1 of January 00h00 to 1 on the 31 of December 23h30}
#   \item{DLS}{daylight savings (“DLS”) an indicator of winter/summer hour}
#   \item{Summer_break}{summer holidays}
#   \item{Christmas_break}{Christmas holidays}
#   \item{Load}{Electricity consumption (MW) at the national level, lags one day (resp. one week) are denoted the subscripts ".48" (resp. ".336")}
#   \item{Nouvelle_A}{Electricity consumption (MW) of the region Nouvelle Aquitaine}
#   \item{Auvergne_R}{Electricity consumption (MW) of the region Auvergne Rhones-Alpes}
#   \item{Bourgogne}{Electricity consumption (MW) of the region Bourgogne-Franche-Comté}
#   \item{Occitanie}{Electricity consumption (MW) of the region Occitanie}
#   \item{Hauts_de_F}{Electricity consumption (MW) of the region Hauts-de-France}
#   \item{Normandie}{Electricity consumption (MW) of the region Normandie}
#   \item{Bretagne}{Electricity consumption (MW) of the region Bretagne}
#   \item{Centre_Val}{Electricity consumption (MW) of the region Centre-Val de Loire}
#   \item{Ile_de_Fra}{Electricity consumption (MW) of the region Île-de-France}
#   \item{Pays_de_la_Loire}{Electricity consumption (MW) of the region Pays de la Loire}
#   \item{Provence_A}{Electricity consumption (MW) of the region Provence-Alpes-Côte d’Azur}
#   \item{Grand_Est}{Electricity consumption (MW) of the region Grand Est}
#   \item{Occitanie}{Electricity consumption (MW) of the region Occitanie}
#   \item{T_region}{Regional temperature in Celsius, weighted means of meteo stations data with weights proportional to exp(−hd^2)
#   where d is the euclidian distance between the position of the station and the centroid of each regions}
#   \item{T_region_s99}{Smoothed regional temperature in Celsius, exponential smoothing of these temperatures with coefficients α=0.99
#   at a given instant t, T_region_s99_t=αT_region_s99_{t−1}+(1−α)T_region_t}
#   \item{T_region_s95}{Smoothed regional temperature in Celsius, exponential smoothing of these temperatures with coefficients α=0.95}
#   \item{T_region_s99_min}{(_max) minimal (resp. maximal) value of these smooth temperatures for each days}
# }
# @source \url{https://opendata.reseaux-energies.fr/} and \url{https://donneespubliques.meteofrance.fr/}
# "regional_load"