library(shapefiles)
library(dplyr)
library(readr)
library(tidyr)


#setwd("C:/Users/mrustl/Desktop/WC_Maxflow/trunk")


### Benoetigt, da MNW2 erwartet, dass Anzahl der Brunnen je Stressperiode = Gesamtbrunnenanzahl
wells_time_dummy <- function (wells_nodes, pers = 0:8) {
  
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


create_mnw2_files <- function(###L_x: muessen mit Parametern Ly und Lx in wellfield.py uebereinstimmen 
                              L_x = 4000, 
                              ###L_y: muessen mit Parametern Ly und Lx in wellfield.py uebereinstimmen 
                              L_y = 5400, 
                              ### if FALSE wabis data will be used for MNW2
                              rwe_model = FALSE, 
                              ### if TRUE (one NODES/TIMES csv file will be created for each stress period), 
                              #### if FALSE: export in one csV file only!
                              separate_periods = TRUE,
                              ### should data be exported in for WELL package format for "flopy"?
                              export_wellpackage = FALSE) {


### Using data from RWE model
if (rwe_model == TRUE) {
  stauer <- shapefiles::read.dbf(dbf.name = "GIS/Stauer.dbf")$dbf
  entnahme <- shapefiles::read.dbf(dbf.name = "GIS/W_Brunnen_bis2015.dbf")$dbf
  
  
  entnahme_pro_jahr_und_brunnen <- entnahme %>% 
    select(Iz, Qbr, Brgwl, Brkenn, X_WERT, Y_WERT) %>% 
    mutate(Year = Iz + 1970, 
           daysPerYear = lubridate::yday(as.Date(sprintf("%s-12-31", Year, format="%Y-%m-%d"))),
           Q_perYear = Qbr*60*24*daysPerYear, 
           Randbrunnen = stringr::str_detect(string = entnahme$Brkenn,pattern = "T  WR|T  WS"), 
           Bru_in_6B = stringr::str_detect(string = entnahme$Brgwl,"xx0.*")) %>%  
    filter(Bru_in_6B == TRUE, 
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
    mutate(label = sprintf("Jahr: %d (Gesamtf\u00F6rderung 6B & 6D: %3.1f Millionen m3)", 
                           Year, 
                           round(Gesamtfoerderung/1000000,1)))
  
  
  entnahme_pro_jahr_und_brunnen <- entnahme_pro_jahr_und_brunnen %>% 
    left_join(entnahme_pro_jahr)
  
  entnahme_pro_jahr_und_brunnen %>%  ungroup() %>% group_by(Year) %>%  summarise(n = n())
  
  
  
  wells_nodes <- entnahme_pro_jahr_und_brunnen %>% 
    ungroup() %>% 
    select(Brkenn, X_WERT, Y_WERT) %>% 
    group_by(Brkenn) %>% 
    summarise(X_WERT = min(X_WERT),
              Y_WERT = min(Y_WERT)) %>% 
    mutate(wellid = sprintf("well%d", 1:n()), 
           y = Y_WERT - min(Y_WERT), 
           x = X_WERT - min(X_WERT), 
           k = 2,
           losstype = "thiem",
           pumploc = 0,
           qlimit = 0,
           ppflag = 0,
           pumpcap = 0,
           rw = 0.5,
           hlim = 6,
           qcut = 0)
  
  wells_nodes <- wells_nodes %>% 
    mutate(x = x + (L_x - max(x)) - 200, ### 200 m Abstand vom rechten Rand
           y = max(y)-y + (L_y - max(y))/2) ### gleicher  Abstand von oberer/unterer Rand
  
wells_times <- entnahme_pro_jahr_und_brunnen %>% 
  left_join(wells_nodes %>% select(Brkenn, wellid)) %>% 
  mutate(per = Year-min(entnahme_pro_jahr_und_brunnen$Year), 
         qdes = -Q_perYear/365) %>%  
  ungroup() %>% 
  select(per,wellid, qdes) %>% 
  arrange(per, wellid)




wells_times <- wells_time_dummy(wells_nodes) %>% 
               left_join(wells_times, by = c("per", "wellid")) %>% 
               mutate(qdes = ifelse(is.na(qdes.y), qdes.x, qdes.y)) %>% 
               #mutate(qdes = ifelse(qdes != 0 ,-1500, 0)) %>%  
               select(per, wellid, qdes) 
} else {
#### Using data from RWE WABIS system
wells_times_raw <- read_csv("real_well_times.csv") %>% 
  dplyr::rename(Brkenn = Brunnen) %>%
  tidyr::gather(key = "year", value = "qdes", -Brkenn, -X_WERT, -Y_WERT) %>% 
  dplyr::mutate(per = as.numeric(year) - 2007, 
                qdes =  -qdes*60*24) %>% 
  dplyr::filter(qdes < 0, 
                X_WERT < 2532000,
                X_WERT >= 2528500,
                Y_WERT < 5662600,
                Y_WERT >= 5656900) 


wells_nodes <- wells_times_raw %>% 
  group_by(Brkenn) %>% 
  summarise(X_WERT = min(X_WERT), 
            Y_WERT = min(Y_WERT)) %>% 
  arrange(Brkenn) %>% 
  mutate(wellid = sprintf("well%d", 1:n()), 
         y = Y_WERT - min(Y_WERT), 
         x = X_WERT - min(X_WERT), 
         k = 2,
         losstype = "thiem",
         pumploc = 0,
         qlimit = 0,
         ppflag = 0,
         pumpcap = 0,
         rw = 0.5,
         hlim = 6,
         qcut = 0)

wells_nodes <- wells_nodes %>% 
  mutate(x = x + (L_x - max(x)) - 200, ### 200 m Abstand vom rechten Rand
         y = max(y)-y + (L_y - max(y))/2) ### gleicher  Abstand von oberer/unterer Rand



wells_times  <- wells_times_raw %>% 
  left_join(wells_nodes %>% select(Brkenn, wellid)) %>% 
  select(per,wellid, qdes) %>% 
  arrange(per, wellid)

wells_times <- wells_time_dummy(wells_nodes) %>% 
              left_join(wells_times, by = c("per", "wellid")) %>% 
              mutate(qdes = ifelse(is.na(qdes.y), qdes.x, qdes.y)) %>%  
              select(per, wellid, qdes) 

}

if (separate_periods == TRUE) {
  dir.create("mnw2")
  for (myPer in unique(wells_times_raw$per)) {
    wells_times_per <- wells_times_raw[wells_times_raw$per == myPer, ]
    
    wells_nodes_per <- wells_times_per %>% 
      group_by(Brkenn) %>% 
      summarise(X_WERT = min(X_WERT), 
                Y_WERT = min(Y_WERT)) %>% 
      arrange(Brkenn) %>% 
      mutate(wellid = sprintf("well%d", 1:n()), 
             y = Y_WERT - min(Y_WERT), 
             x = X_WERT - min(X_WERT), 
             k = 2,
             losstype = "thiem",
             pumploc = 0,
             qlimit = 0,
             ppflag = 0,
             pumpcap = 0,
             rw = 0.5,
             hlim = 6,
             qcut = 0)
    write.csv(wells_nodes_per, 
              sprintf("mnw2/wells_nodes_%d.csv", myPer),
              row.names = FALSE)
    
    write.csv(wells_times_per, 
              sprintf("mnw2/wells_times%d.csv", myPer),
              row.names = FALSE)
    
    
  }
} else {

###############################
### Export files for MNW2 
###############################
write.csv(wells_nodes, 
          "wells_nodes.csv",
          row.names = FALSE)

write.csv(wells_times, 
          "wells_times.csv",
          row.names = FALSE)

}
###############################
### Create files for well package 
###############################
if (export_wellpackage == TRUE) {
wells_nodes %>%
  mutate(k = 2, 
         i = y, 
         j = x) %>% 
  select(wellid, k, i, j) %>%
  left_join(wells_times %>% 
              mutate(flux = qdes) %>% 
              select(per, wellid, flux)) %>% 
  filter(flux < 0) %>% 
  arrange(per, wellid) %>% 
  select(per, k, i, j, flux) %>% 
  write.csv("wellpackage.csv",
            row.names = FALSE)
}
}

python <- FALSE

if (python) {
  myArgs <- commandArgs(trailingOnly = TRUE)
  create_mnw2_files(L_x = as.numeric(myArgs[1]),
                    L_y = as.numeric(myArgs[2]),
                    rwe_model = as.logical(myArgs[3]),
                    separate_periods = as.logical(myArgs[4]))
} else {
create_mnw2_files(L_x = 4000,
                  L_y = 5400,
                  rwe_model = TRUE,
                  separate_periods = FALSE)
}

