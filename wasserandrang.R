
wasserandrang <- function(initial_head = 80, #m
                          aquifer_top = 20, #m
                          kf = 1E-4, #m/s
                          well_radius = 1, #m
                          approx = "Sichmann",
                          levels = seq(0,20,5)
                          ) {


if (approx == "Sichmann") const <- 15
if (approx == "Huismann") const <- 30 
  
vf_max <- sqrt(kf)/const 


wad_cbmPerSecond <- 2*well_radius*pi*c(0,aquifer_top)*vf_max


wad_cbmPerDay <- wad_cbmPerSecond*24*3600


label <- sprintf("r_well: %1.1f, kf: %1.5f, vf_max Approx. nach %s", 
                 well_radius,kf,  approx)


plot(c(wad_cbmPerDay, wad_cbmPerDay[2]), c(0,aquifer_top, initial_head), 
     type="l", 
     las=1, 
     xlab = "Wasserandrang (m³/Tag)", 
     ylab = "Grundwasserstand", 
     main = label)
abline(h = aquifer_top, col="red")
abline(h = initial_head, col="blue")
text(x=mean(wad_cbmPerDay), y=aquifer_top+2, labels = "GWL gespannt-ungespannt", col="red")
text(x=mean(wad_cbmPerDay), y=initial_head+2, labels = "Anfangswasserstand", col="blue")

res <- data.frame(levels = levels, wad_cbmPerDay = 2*well_radius*pi*vf_max*levels*24*3600)
return(res )
}

wasserandrang(aquifer_top = 20,kf = 2e-5)
