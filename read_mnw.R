library(data.table)
library(dplyr)
library(ggplot2)
library(lattice)
library(ggplot2)
library(devtools)
library(foreign)
if(!require("gganimate")) { 
  devtools::install_github("dgrtwo/gganimate", dependencies = TRUE)
}
library(gganimate)


### Arbeitsverzeichnis definieren
setwd("C:/Users/RE73858/Desktop/WC_Maxflow/Maxflow-combined_leakage")

### importiere head-file aus FloPy und definiere Teilgebiet
head_master <- read.dbf(file = "test_heads_sp12.dbf", as.is = FALSE) %>%
               dplyr::filter(row >=68,
                             row <=88,
                             column >=40,
                             column <=60)
mean_zonehead <-  head_master %>% # Berechnung des mittleren Wasserstands im Teilgebiet
  summarise(Wasserstand = mean(head002))

### Importiere Brunnenstammdaten
well_master <- read.csv(file = "wells_nodes.csv",header = TRUE) %>% 
               dplyr::filter(k == 2) %>% 
               dplyr::mutate(wellid = toupper(wellid)) %>% 
               dplyr::rename(WELLID = wellid) %>%
               dplyr::select(WELLID, X_WERT, Y_WERT, Brkenn)

### Importiere Brunnenbetriebsdaten
wells <- data.table::fread(input = "wellfield.byn", fill = TRUE)
names(wells) <- c("WELLID", "NODE", "Lay", "Row", "Col", "Totim", "Q_node", "hwell", "hcell", "seepage")
wells <- wells[,c("NODE", "Lay") := NULL]
save(wells, file = "wells.RData")


### Schreibe Brunnenbetriebsdaten
wells_perYear <- wells %>% 
  #filter(WELLID == "WELL26") %>% 
  mutate(Totim = floor((Totim-1)/365)) %>% 
  group_by(WELLID, Totim) %>% 
  summarise(Q_node_cubicmpermin = round(-sum(Q_node)/(365*24*60),2), 
            hwell_median = median(hwell), 
            hcell_median = median(hcell)) %>% 
  right_join(y = well_master)

### Berechnung Gesamtentnahme pro Jahr
entnahme_pro_jahr <-wells_perYear %>% 
  ungroup() %>% 
  group_by(Totim) %>% 
  summarise(Gesamtfoerderung = sum(Q_node_cubicmpermin)) %>% 
  mutate(label = sprintf("Jahr: %d (Gesamtfoerderung 6B: %3.1f Millionen m3)", 
                         Totim+2007, 
                         round((Gesamtfoerderung*24*60*365)/1000000,1)))

### Schreibe Brunnenbetriebsdaten für Teilgebiet
wells_perYear <- wells_perYear %>% 
                 dplyr::left_join(y = entnahme_pro_jahr) %>% 
                  dplyr::filter(X_WERT >= 2529190,
                                X_WERT <=	2530190,
                                Y_WERT <=	5658790,
                                Y_WERT >=	5657790,
                                Totim >=9)

### Gesamtentnahme im Prognosezeitraum für Teilgebiet
entnahme_prognose <- wells_perYear %>%
  ungroup() %>%
  group_by(Totim) %>% 
  summarise(Gesamtfoerderung = sum(Q_node_cubicmpermin))

### momentane Einzelentnahme im Prognosezeitraum für Teilgebiet
Qmom_prognose <- wells_perYear %>%
  ungroup() %>%
  group_by(Totim) %>% 
  summarise(Q_mom = mean(Q_node_cubicmpermin))

### Animation der Einzelentnahme pro Jahr
p <- ggplot(wells_perYear, aes(x=X_WERT, 
                               y=Y_WERT, 
                               size = Q_node_cubicmpermin,  
                               frame = label)) + 
  geom_point(col = "darkblue", ) + 
  labs(x = "X Koordinate", 
       y = "Y Koordinate") + 
  theme_bw()

print(p)

gganimate(p, "entnahme.html",
          ani.width = 1000, 
          ani.height =1080)


######################################
#######Ergebnisdateien ausgeben#######
######################################

write.csv(wells_perYear, 
          "wells_Qperyear.csv",
          row.names = FALSE)

write.csv(head_master, 
          "zonehead_L2_2018.csv",
          row.names = FALSE)

write.csv(mean_zonehead, 
          "mean_head_end.csv",
          row.names = FALSE)

write.csv(entnahme_prognose, 
          "total_Qperyear.csv",
          row.names = FALSE)

write.csv(Qmom_prognose, 
          "Qmomperyear.csv",
          row.names = FALSE)

