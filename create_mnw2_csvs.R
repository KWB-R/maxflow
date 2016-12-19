library(shapefiles)
library(dplyr)


#setwd("C:/Users/mrustl/Desktop/WC_Maxflow/trunk")

stauer <- shapefiles::read.dbf(dbf.name = "GIS/Stauer.dbf")$dbf
entnahme <- shapefiles::read.dbf(dbf.name = "GIS/W_Brunnen_bis2015.dbf")$dbf


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




wells_nodes <- entnahme_pro_jahr_und_brunnen %>% 
  ungroup() %>% 
  select(Brkenn, X_WERT, Y_WERT) %>% 
  group_by(Brkenn) %>% 
  summarise(X_WERT = min(X_WERT),
            Y_WERT = min(Y_WERT)) %>% 
  mutate(wellid = sprintf("well%d", 1:n()), 
         y = Y_WERT - min(Y_WERT), 
         x = X_WERT - min(X_WERT), 
         ztop = 10,
         zbotm = 0,
         losstype = "thiem",
         pumploc = 0,
         qlimit = 0,
         ppflag = 0,
         pumpcap = 0,
         rw = 0.5,
         hlim = 6,
         qcut = 0)


### L_x und L_y: muessen mit Parametern Ly und Lx in wellfield.py uebereinstimmen 
L_x <- 4000
L_y <- 5400

wells_nodes <- wells_nodes %>% 
              mutate(x = x + (L_x - max(x)) - 200, ### 200 m Abstand vom rechten Rand
                     y = y + (L_y - max(y))/2) ### gleicher  Abstand von oberer/unterer Rand

write.csv(wells_nodes, 
          "wells_nodes.csv",
          row.names = FALSE)

### Benoetigt, da MNW2 erwartet, dass Anzahl der Brunnen je Stressperiode = Gesamtbrunnenanzahl
wells_time_dummy <- function (wells_nodes, pers = 0:9) {
  
  for (myPeriod in pers) {
    tmp <- wells_nodes %>% 
      mutate_(per = ~myPeriod,
              qdes = 0) %>%
      select_(~per,~wellid, ~qdes) 
    if(myPeriod == pers[1]) {
      res <- tmp
    } else {
      res <- rbind(res,tmp)
    }
    
  }
  return(res)  
}

wells_time_dummy(wells_nodes)



wells_times <- entnahme_pro_jahr_und_brunnen %>% 
  left_join(wells_nodes %>% select(Brkenn, wellid)) %>% 
  mutate(per = 1+Year-min(Year), 
         qdes = -Q_perYear/365) %>%  
  ungroup() %>% 
  select(per,wellid, qdes) %>% 
  arrange(per, wellid)



wells_times <- wells_time_dummy(wells_nodes) %>% 
               left_join(wells_times, by = c("per", "wellid")) %>% 
               mutate(qdes = ifelse(is.na(qdes.y), qdes.x, qdes.y)) %>% 
               select(per, wellid, qdes) 
  

write.csv(wells_times, 
          "wells_times.csv",
          row.names = FALSE)
