boreholes_number <- 300
area_boreholes <- 150 #meter^2
area_model <- 5400*4000  #meter^2 

kf_boreholes <- 1 #meter/day
kf_natural <- 1E-6 #meter/day


get_realLeakage <- function (area_boreholes = 150, #meter^2
          area_model = 5400*4000,  #meter^2 
          kf_boreholes = 1, #meter/day
          kf_natural = 1E-6 #meter/day
         ) {
(area_boreholes * kf_boreholes + (area_model - area_boreholes) * kf_natural)/area_model 
}


kf_boreholes <- 0.01
kf_natural <- c(1E-5, 1E-6, 1E-7, 1E-8, 1E-9, 1E-10, 1E-15, 1E-100)


get_realLeakage(kf_boreholes = kf_boreholes, kf_natural = kf_natural)

plot(log(kf_natural), get_realLeakage(kf_boreholes = kf_boreholes, kf_natural = kf_natural), type = "l")
