get_realLeakage <- function (area_boreholes = 150, #meter^2
          area_model = 5400*4000,  #meter^2 
          kf_boreholes = 1, #meter/day
          kf_natural = 1E-6 #meter/day
         ) {
(area_boreholes * kf_boreholes + (area_model - area_boreholes) * kf_natural)/area_model 
}

boreholes_number <- 300
area_boreholes <- 150 #meter^2
area_model <- 5400*4000  #meter^2 
kf_boreholes <- 0.01
kf_natural <- c(1E-4, 1E-5, 1E-6, 1E-7, 1E-8, 1E-9, 1E-10, 1E-15)

get_realLeakage(kf_boreholes = kf_boreholes, kf_natural = kf_natural)

plot(log10(kf_natural), 
     log10(get_realLeakage(kf_boreholes = kf_boreholes, 
                           kf_natural = kf_natural)), 
     type = "b", 
     xlab = "Log10 natural leakage", 
     ylab = "Log10 combined leakage (natural + borehole)", 
     main = sprintf("Boreholes (number: %d, total area: %d m2, kf: %1.5f m/s)\nModel area: %d m2", 
                    boreholes_number, 
                    area_boreholes, 
                    kf_boreholes, 
                    area_model), 
     las = 1)
# text(log10(min(kf_natural))+2, 
#      log10(max(get_realLeakage(kf_boreholes = kf_boreholes, kf_natural = kf_natural))), 
#      labels = "area_boreholes * kf_boreholes + (area_model - area_boreholes) * kf_natural)/area_model")