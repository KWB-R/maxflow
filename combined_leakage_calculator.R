get_realLeakage <- function (area_boreholes = 250, #meter^2
          area_model = 5400*4000,  #meter^2 
          kf_boreholes = 1, #meter/day
          kf_natural = 1E-6 #meter/day
         ) {
(area_boreholes * kf_boreholes + (area_model - area_boreholes) * kf_natural)/area_model 
}

boreholes_number <- 500
area_boreholes <- 250 #meter^2
area_model <- 5400*4000  #meter^2 
kf_boreholes <- 0.01
kf_natural <- c(1E-4, 1E-5, 1E-6, 1E-7, 1E-8, 1E-9, 1E-10, 1E-11)

get_realLeakage(kf_boreholes = kf_boreholes, kf_natural = kf_natural)

plot(log10(kf_natural), 
     log10(get_realLeakage(kf_boreholes = kf_boreholes, 
                           kf_natural = kf_natural)), 
     type = "b", 
     xlab = "Log10 geogener kf-Wert (m/s)", 
     ylab = "Log10 mittlerer kf_Wert (m/s)", 
     xlim=c(-11, -4), ylim=c(-10, -4), 
     main = sprintf("Bohrloch-Fläche: %d m2 (Anzahl: %d, mittlere Bohrlochfläche: %1.2f m2)\nGesamtfläche: %d m2", 
                    area_boreholes,
                    boreholes_number, 
                    area_boreholes/boreholes_number, 
                    area_model), 
     col = "blue",
     pch = 1,
     las = 1)
text(-10, -7.2, 
     labels=bquote(paste("kf"[Bohrloch], "= 1E-2 m/s")), 
     col = "blue",
     font=2)
points(log10(kf_natural), 
      log10(get_realLeakage(kf_boreholes = 0.001, 
                            kf_natural = kf_natural)), 
           type="b", 
      col = "darkgreen",
      pch = 2,
      lty = 2)
text(-10, -8.2, 
     labels=bquote(paste("kf"[Bohrloch], "= 1E-3 m/s")), 
     col = "darkgreen",
     font=2)
points(log10(kf_natural), 
       log10(get_realLeakage(kf_boreholes = 0.0001, 
                             kf_natural = kf_natural)), 
       type="b", 
       col = "red",
       pch = 3,
       lty = 3)
text(-10, -9.2, 
     labels=bquote(paste("kf"[Bohrloch], "= 1E-4 m/s")), 
     col = "red",
     font=2)

# text(log10(min(kf_natural))+2, 
#      log10(max(get_realLeakage(kf_boreholes = kf_boreholes, kf_natural = kf_natural))), 
#      labels = "area_boreholes * kf_boreholes + (area_model - area_boreholes) * kf_natural)/area_model")