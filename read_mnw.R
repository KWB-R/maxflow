library(data.table)
library(dplyr)
library(ggplot2)
library(lattice)
library(ggplot2)
library(devtools)
if(!require("gganimate")) { 
  devtools::install_github("dgrtwo/gganimate", dependencies = TRUE)
}
library(gganimate)

setwd("C:/Users/cmenz/Desktop/WC_Maxflow/branches/combined_leakage")

well_master <- read.csv(file = "wells_nodes.csv",header = TRUE) %>% 
               dplyr::filter(k == 2) %>% 
               dplyr::mutate(wellid = toupper(wellid)) %>% 
               dplyr::rename(WELLID = wellid) %>%
               dplyr::select(WELLID, X_WERT, Y_WERT, Brkenn)


wells <- data.table::fread(input = "wellfield.byn", fill = TRUE)

names(wells) <- c("WELLID", "NODE", "Lay", "Row", "Col", "Totim", "Q_node", "hwell", "hcell", "seepage")

wells <- wells[,c("NODE", "Lay", "Row", "Col") := NULL]

save(wells, file = "wells.RData")


#load("wells.RData")

wells_perYear <- wells %>% 
  #filter(WELLID == "WELL26") %>% 
  mutate(Totim = floor((Totim-1)/365)) %>% 
  group_by(WELLID, Totim) %>% 
  summarise(Q_node_cubicmpermin = round(-sum(Q_node)/(365*24*60),2), 
            hwell_median = median(hwell), 
            hcell_median = median(hcell)) %>% 
  right_join(y = well_master)


entnahme_pro_jahr <-wells_perYear %>% 
  ungroup() %>% 
  group_by(Totim) %>% 
  summarise(Gesamtfoerderung = sum(Q_node_cubicmpermin)) %>% 
  mutate(label = sprintf("Jahr: %d (Gesamtfoerderung 6B: %3.1f Millionen m3)", 
                         Totim+2007, 
                         round((Gesamtfoerderung*24*60*365)/1000000,1)))

wells_perYear <- wells_perYear %>% 
                 dplyr::left_join(y = entnahme_pro_jahr)

#animation
p <- ggplot(wells_perYear, aes(x=X_WERT, 
                               y=Y_WERT, 
                               size = Q_node_cubicmpermin,  
                               frame = label)) + 
  geom_point(col = "darkblue", ) + 
  labs(x = "X Koordinate", 
       y = "Y Koordinate") + 
  theme_bw() +
  geom_text_repel(aes(label=wellid), size = 3)


print(p)

gganimate(p, "entnahme.html",
          ani.width = 1000, 
          ani.height =1080)

plot(Q_perYear ~ Year, 
     data=entnahme_pro_jahr, 
     type="b", 
     col="blue", 
     las = 1)

# wells_perDay <- wells %>% 
#          #filter(WELLID == "WELL26") %>% 
#          mutate(Totim = round(Totim,0)) %>% 
#          group_by(WELLID, Totim) %>% 
#          summarise(Q_node_median = median(Q_node), 
#                    hwell_median = median(hwell), 
#                    hcell_median = median(hcell))


wells_summary <- wells_perDay %>% 
  ungroup() %>% 
  group_by(WELLID) %>%  
  filter(Q_node_median < 0) %>%
  summarise(Q_node_median = -median(Q_node_median), 
            days_active = n()) %>% 
  arrange(Q_node_median, desc(days_active))

wells_summary %>%  
  ggplot(aes(x = days_active, y = Q_node_median/24/60)) + 
  geom_point(alpha = 0.5, size = 2) + 
  labs(y = "Foerderrate pro Minute (m³/min)") +
  theme_bw()

# Q sum per well for whole simulation
wells_perDay %>% 
  ungroup() %>% 
  mutate(WELLID = as.numeric(gsub("WELL", "", .$WELLID))) %>%  
  group_by(WELLID) %>% 
  summarise(Qsum = -sum(Q_node_median)) %>%  
  plot(Qsum ~ WELLID, data = ., type="b", col="blue")


wells_perDay_tidy <- tidyr::gather(data = wells_perDay,Key, Value, -Totim, -WELLID)

pdf("wells_Q.pdf",height = 7, width = 10)
xyplot(Q_node_median ~ Totim | WELLID,
       data = wells_perDay,
       layout = c(1,1))
dev.off()

pdf("wellfield_Q.pdf",height = 7, width = 10)
xyplot(Q_node_median ~ Totim, group = WELLID,
       data = wells_perDay)
dev.off()

pdf("wells_H.pdf",height = 7, width = 10)
xyplot(Value ~ Totim | WELLID,
       group = Key,
data = wells_perDay_tidy %>% filter(Key !="Q_node_median"),
layout = c(1,1))
dev.off()

# if (FALSE) {
# ggplot(wells_perDay, aes(x = Totim, y = Q_node_median)) +
#   facet_wrap(~ WELLID) +
#   geom_line() +
#   theme_bw()
# 
# 
# plot(well1$Totim, well1$hcell_median,ylim = c(0,80))
# points(well1$Totim, well1$hwell_median)
# 
# plot(well1$Totim, well1$Q_node_median)
# 
# 
# setkey(wells,Totim)
# wells[WELLID == "WELL1",
#       list(Q-node_median=median(Q-node),
#            hwell_median = median(hwell),
#            hcell_median = median(hcell)),
#       by=Totim]
# 
# wells$WELLID <- as.numeric(gsub("WELL", "", wells$WELLID))
# 
# 
# 
# 
# setkey(wells,Totim)
# system.time(wells[,list(Q-node_median=median(Q-node),
#                         hwell_median = median(hwell),
#                         hcell_median = median(hcell)),
#                   by=Totim])
# 
# 
# well1 <- wells[WELLID == "WELL1", ]
# }

