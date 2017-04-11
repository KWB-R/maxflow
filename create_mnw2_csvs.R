library(dplyr)
library(reshape2)



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


neubrunnen_nodes <- read.csv("GIS/alleBrunnen_2016-2018.csv", 
                             stringsAsFactors = FALSE) %>% 
  dplyr::select_("Brkenn", 
                 "X_WERT", 
                 "Y_WERT", 
                 "k") %>% 
  dplyr::group_by_("Brkenn", "k") %>% 
  dplyr::summarise(X_WERT = min(X_WERT),
                   Y_WERT = min(Y_WERT)) 


neubrunnen_times <- read.csv("GIS/alleBrunnen_Wabis_bis 2020.csv") %>%  
  dplyr::mutate(Year = Iz + 1970) %>% 
  dplyr::filter(Brkenn %in% unique(neubrunnen_nodes$Brkenn),
                Year >= 2016,
                Year <= 2018)


neubrunnen_nodes <- neubrunnen_nodes %>% 
  dplyr::filter(Brkenn %in% unique(neubrunnen_times$Brkenn)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(wellid = sprintf("welll%d",
                                 (nrow(altbrunnen_nodes) + 1):(nrow(altbrunnen_nodes) + n())))


neubrunnen_times <- neubrunnen_times %>% 
                    dplyr::mutate(qdes = Qbr*24*60,
                                  per = Year - 2007) %>% 
                    dplyr::left_join(neubrunnen_nodes %>% select(Brkenn, wellid)) %>% 
                      dplyr::select_("per", 
                                     "Year",
                                     "wellid",
                                     "Brkenn",
                                     "qdes"
                      )  

wells_nodes <- rbind(altbrunnen_nodes,neubrunnen_nodes) %>% 
               dplyr::mutate(losstype = "thiem",
                            pumploc = 0,
                            qlimit = 0,
                            ppflag = 0,
                            pumpcap = 0,
                            rw = 0.5,
                            hlim = 6,
                            qcut = 0)

wells_times <- rbind(altbrunnen_times,neubrunnen_times)

max_per <- max(wells_times$per)



wells_times <- wells_time_dummy(wells_nodes,pers = 0:max_per) %>% 
              left_join(wells_times, by = c("per", "wellid", "Brkenn")) %>% 
              mutate(qdes = ifelse(is.na(qdes.y), qdes.x, qdes.y), 
                     Year = per + 2007) %>% 
              select(per, Year, wellid, Brkenn, qdes) 

#wells_times %>% group_by(wellid,per) %>% summarise(n = n()) %>% filter(n > 1) %>% View()


if (FALSE) {
  
well_number <- function(critDist = 100,
                        df = neubrunnen_nodes %>% filter(k == 0)) {

  distanceMatrix <- round(dist(df[,c("X_WERT","Y_WERT")],
                             method = "euclidian", 
                             upper = TRUE,
                             diag = FALSE),
                        1)
attr(distanceMatrix, "Labels") <- df$Brkenn
# cat("Well distance matrix table:")
# cat("")
# distanceMatrix



df <- melt(as.matrix(distanceMatrix), varnames = c("row", "col"))
df$row <- as.character(df$row)
df$col <- as.character(df$col)

df <- df %>% dplyr::filter(value > 0) %>% dplyr::arrange(row, value) 

res <- df %>% dplyr::filter(value < critDist) %>% 
  dplyr::arrange(row, value) 

length(unique(unique(res$row), unique(res$col)))

}

critDist <- seq(50,300,5)

wells <- vector()
i <- 0
for(dist in critDist) {
i <- i + 1 
wells[i] <- well_number(dist)
}

plot(x = critDist,y = wells)
}


###############################
### Export files for MNW2 
###############################
write.csv(wells_nodes, 
          "wells_nodes.csv",
          row.names = FALSE)

write.csv(wells_times, 
          "wells_times.csv",
          row.names = FALSE)

