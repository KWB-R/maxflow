setwd("C:/Users/mrustl/Desktop/WC_Maxflow/trunk/GIS")
library(shapefiles)
library(dplyr)
library(ggplot2)
library(devtools)
if(!require("gganimate")) { 
devtools::install_github("dgrtwo/gganimate", dependencies = TRUE)
}
library(gganimate)

stauer <- shapefiles::read.dbf(dbf.name = "Stauer.dbf")$dbf
entnahme <- shapefiles::read.dbf(dbf.name = "W_Brunnen_bis2015.dbf")$dbf


entnahme_pro_jahr_und_brunnen <- entnahme %>% 
 select(Iz, Qbr, Brgwl, Brkenn, X_WERT, Y_WERT) %>% 
 mutate(Year = Iz + 1970, 
        daysPerYear = lubridate::yday(as.Date(sprintf("%s-12-31", Year, format="%Y-%m-%d"))),
         Q_perYear = Qbr*60*24*daysPerYear, 
        Randbrunnen = stringr::str_detect(string = entnahme$Brkenn,pattern = "T  WR|T  WS")) %>%  
 filter(Brgwl == "xx0987xxxxxx", 
        Randbrunnen == FALSE,
        X_WERT < 2532000,
        X_WERT >= 2528500,
        Y_WERT < 5662600,
        Y_WERT >= 5656900, 
        Year >= 2007) %>% 
 group_by(Brkenn, X_WERT, Y_WERT,Year) %>%  
 summarise(Q_perYear=-sum(Q_perYear)) 


entnahme_pro_jahr <- entnahme_pro_jahr_und_brunnen %>% 
  ungroup() %>% 
  group_by(Year) %>% 
  summarise(Gesamtfoerderung = sum(Q_perYear)) %>% 
  mutate(label = sprintf("Jahr: %d (Gesamtförderung 6B & 6D: %3.1f Millionen m3)", 
                         Year, 
                         round(Gesamtfoerderung/1000000,1)))


entnahme_pro_jahr_und_brunnen <- entnahme_pro_jahr_und_brunnen %>% 
                                  left_join(entnahme_pro_jahr)



p <- ggplot(entnahme_pro_jahr_und_brunnen, aes(x=X_WERT, 
                                   y=Y_WERT, 
                                   size = Q_perYear,  
                                   frame = label)) + 
  geom_point(col = "darkblue") + 
  labs(x = "X Koordinate", 
       y = "Y Koordinate") + 
  theme_bw()


print(p)
  
gganimate(p, "entnahme.html",
          ani.width = 1000, 
          ani.height = 600)

plot(Q_perYear ~ Year, 
     data=entnahme_pro_jahr, 
     type="b", 
     col="blue", 
     las = 1)

# 
#   

