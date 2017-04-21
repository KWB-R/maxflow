library(dplyr)
library(reshape2)

setwd("C:/Users/RE73858/Desktop/WC_Maxflow/Maxflow-master")

### Benoetigt, da MNW2 erwartet, dass Anzahl der Brunnen je Stressperiode = Gesamtbrunnenanzahl
wells_time_dummy <- function(wells_nodes, pers = 0:11) {
  
  for (myPeriod in pers) {
    tmp <- wells_nodes %>% 
      mutate_(per = ~myPeriod,
              qdes = 0) %>%
      select_(~per,~wellid, ~Brkenn, ~qdes) 
    if (myPeriod == pers[1]) {
      res <- tmp
    } else {
      res <- rbind(res,tmp)
    }
    
  }
  return(res)  
}

## Auffüllen des Prognosezeitraums mit Förderraten aus dem letzten Jahr der Kalibrierung 
add_last_q <- function(df_old_times, # altbrunnen_times 
                       df_inOperation, #altbrunnen
                       max_per = 11) {
  
  
  
    last_q <- df_old_times %>% 
              dplyr::filter(per == max(per)) 
  
  
  max_per_old <- max(df_old_times$per)
  
  new_pers <- (max_per_old + 1):max_per
  
  for (ind in seq_along(new_pers)) {
    tmp <- last_q %>% 
      dplyr::mutate(per = new_pers[ind],
                    Year = 2007 + per) %>% 
      dplyr::select_("per", 
                     "Year",
                     "wellid",
                     "Brkenn",
                     "qdes")
    
    if (ind == 1) {
      res <- tmp
    } else {
      res <- rbind(res,tmp)
    }
  }
  
  
  in_operation <- df_inOperation %>% dplyr::select_("Brkenn", 
                                                    "endtime")    
  
  for (brkenn in unique(res$Brkenn)) {
    
    endtime <- max(in_operation$endtime[in_operation$Brkenn == brkenn])
    
    rows_to_drop <-  res$per > endtime & res$Brkenn == brkenn   
    res <-  res[!rows_to_drop,]
  } 
  
  
  return(res)
}

## Hilfsfunktion zur Ermittlung des Brunnenabstandes ##
get_well_distances <- function(df_nodes,
                               col_x = "X_WERT",
                               col_y = "Y_WERT",
                               col_labels = "Brkenn") {
  
  distanceMatrix <- round(dist(df_nodes[,c(col_x,col_y)],
                               method = "euclidian", 
                               upper = TRUE,
                               diag = FALSE),
                          1)
  attr(distanceMatrix, "Labels") <- df_nodes[,col_labels] %>% unlist(use.names = FALSE)
  
  
  dist <- melt(as.matrix(distanceMatrix), varnames = c("row", "col"))
  dist$row <- as.character(dist$row)
  dist$col <- as.character(dist$col)
  
  dist <- dist %>% 
    dplyr::filter(value > 0) %>% 
    dplyr::arrange(row, value) 
  
  return(dist)
}

## Vergrößern des Brunnenrasters (nur Brunnen deren Abstand > krit. Distanz ist!)
get_remaining_wells <- function(df_nodes,
                                critDist = 300) {
  
  well_distances <- get_well_distances(df_nodes = df_nodes)
  
  dist_df <- well_distances %>% 
    dplyr::filter(value < critDist) %>% 
    dplyr::arrange(row, value) 
  
  wells_to_remain <- unique(dist_df$row)
  
  if (length(wells_to_remain) == 0) {
    msg <- sprintf("Fehler: keine Neubrunnen erfüllen das Kriterium: Entfernung <= %4.0f m.\n
Bitte einen größeren Wert für Parameter 'critDist' in Funktion get_remaining_wells() wählen!", 
                   critDist)
    stop(msg)
    }
  return(list(names = wells_to_remain,
              distanceMatrix = dist_df))
}


## Einlesen der Brunnen für Kalibrierungszeitraum ##

altbrunnen <- read.csv("GIS/alleBrunnen_BIOS_2007-2015.csv",
                       stringsAsFactors = FALSE) %>% 
              dplyr::mutate(Year = 2007 + per)


altbrunnen_nodes <- altbrunnen %>%  
                    dplyr::select_("Brkenn", 
                                   "X_WERT", 
                                   "Y_WERT", 
                                   "k") %>% 
                    dplyr::group_by_("Brkenn", "k") %>% 
                    dplyr::summarise(X_WERT = min(X_WERT),
                                     Y_WERT = min(Y_WERT)) %>% 
                    dplyr::ungroup() %>% 
                    dplyr::mutate(wellid = sprintf("welll%d",1:n()))

## Einlesen der Förderraten für Kalibrierungszeitraum ##

altbrunnen_times <- altbrunnen %>% 
                    dplyr::select_("Brkenn", 
                                   "k",
                                   "qdes",
                                   "per", 
                                   "Year") %>% 
                    left_join(altbrunnen_nodes %>% select(Brkenn,k,wellid)) %>% 
                    dplyr::select_("per", 
                                   "Year",
                                   "wellid",
                                   "Brkenn",
                                   "qdes"
                                   )  

altbrunnen_times <- rbind(altbrunnen_times,
                          add_last_q(df_old_times = altbrunnen_times,
                                     df_inOperation = altbrunnen,
                                     max_per = 11))

## Einlesen der Brunnen für Prognosezeitraum ##

neubrunnen_nodes <- read.csv("GIS/alleBrunnen_2016-2018.csv", 
                             stringsAsFactors = FALSE) %>% 
  dplyr::select_("Brkenn", 
                 "X_WERT", 
                 "Y_WERT", 
                 "k") %>% 
  dplyr::group_by_("Brkenn", "k") %>% 
  dplyr::summarise(X_WERT = min(X_WERT),
                   Y_WERT = min(Y_WERT)) 


## Einlesen der Förderraten für Prognoseszeitraum ##

neubrunnen_times <- read.csv("GIS/alleBrunnen_Wabis_bis 2020.csv") %>%  
  dplyr::mutate(Year = Iz + 1970) %>% 
  dplyr::filter(Brkenn %in% unique(neubrunnen_nodes$Brkenn),
                Year >= 2016,
                Year <= 2018)


neubrunnen_nodes <- neubrunnen_nodes %>% 
  dplyr::filter(Brkenn %in% unique(neubrunnen_times$Brkenn))

################################################################################
## Vergrößern des Brunnenrasters##
################################################################################

#### Hilfs-Abbildung zur Bestimmung eines optimalen Wertes für 'critDist'  
if (FALSE) {
  
  critDist <- seq(50,300,5)
  
  
  wells <- vector()
  i <- 0
  for (dist in critDist) {
    i <- i + 1 
    wells[i] <- length(get_remaining_wells(df_nodes = neubrunnen_nodes %>% filter(k == 0),
                                           critDist = dist)$names)
  }
  
  plot(x = critDist,y = wells,
       ylim = rev(range(wells)), 
       xlab = "Brunnenabstand [m]", 
       ylab = "Anzahl Brunnen",
       las = 1)
}

# Berücksichtigt nur Brunnen deren Abstände > krit. Entfernung ist 
# (falls Parameter 'critDist' zu klein gewählt -> FEHLER, daher Abschätzung
# zu sinnvollen Wertebereich mittels obigem Plot durchführen!!!!)

critDist <- 150 ### hier ändern um Brunnenraster anzupassen!!!!!!!!!!

neubrunnen_remaining <- get_remaining_wells(df_nodes = neubrunnen_nodes %>% 
                                            filter(k == 0),
                                            critDist = critDist)

## Nur zum Testen -> Entfernungsmatrix
#neubrunnen_remaining$distanceMatrix


### Eliminieren von Neubrunnen die nicht 'critDist' Kriterium erfüllen
neubrunnen_nodes <- neubrunnen_nodes %>%  
                    filter(Brkenn %in% neubrunnen_remaining$names)

neubrunnen_times <- neubrunnen_times %>%  
                    filter(Brkenn %in% neubrunnen_remaining$names)
################################################################################
################################################################################
################################################################################

neubrunnen_nodes <- neubrunnen_nodes %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(wellid = sprintf("welll%d",
                                 (nrow(altbrunnen_nodes) + 1):(nrow(altbrunnen_nodes) + n())))

Q_mom <- 1.0  # Faktor für Förderrate der Neuanlagen aus GWL 6B im Prognosezeitraum

neubrunnen_times <- neubrunnen_times %>% 
                    dplyr::mutate(qdes = Qbr*24*60*Q_mom,  
                                  per = Year - 2007) %>% 
                    dplyr::left_join(neubrunnen_nodes %>% select(Brkenn, wellid)) %>% 
                      dplyr::select_("per", 
                                     "Year",
                                     "wellid",
                                     "Brkenn",
                                     "qdes"
                      )  
## Erstellung der modflow input datei well_nodes ##

wells_nodes <- rbind(altbrunnen_nodes,neubrunnen_nodes) %>% 
               dplyr::mutate(losstype = "thiem",
                            pumploc = 0,
                            qlimit = 0,
                            ppflag = 0,
                            pumpcap = 0,
                            rw = 0.5,
                            hlim = 6,
                            qcut = 0)

## Erstellung der modflow input datei well_times ##

wells_times <- rbind(altbrunnen_times,neubrunnen_times)

max_per <- max(wells_times$per)



wells_times <- wells_time_dummy(wells_nodes,pers = 0:max_per) %>% 
              left_join(wells_times, by = c("per", "wellid", "Brkenn")) %>% 
              mutate(qdes = ifelse(is.na(qdes.y), qdes.x, qdes.y), 
                     Year = per + 2007) %>% 
              select(per, Year, wellid, Brkenn, qdes) 

#wells_times %>% group_by(wellid,per) %>% summarise(n = n()) %>% filter(n > 1) %>% View()




###############################
### Export files for MNW2 
###############################
write.csv(wells_nodes, 
          "wells_nodes.csv",
          row.names = FALSE)

write.csv(wells_times, 
          "wells_times.csv",
          row.names = FALSE)

