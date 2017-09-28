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

setwd("C:/Users/RE73858/Desktop/WC_Maxflow/Maxflow-combined_leakage")

head_master <- read.dbf(file = "test_heads_sp12.dbf", as.is = FALSE) %>%
               dplyr::filter(row >=68,
                             row <=88,
                             column >=40,
                             column <=60)

mean_zonehead <-  head_master %>%
  summarise(Wasserstand = mean(head002))

well_master <- read.csv(file = "wells_nodes.csv",header = TRUE) %>% 
               dplyr::filter(k == 2) %>% 
               dplyr::mutate(wellid = toupper(wellid)) %>% 
               dplyr::rename(WELLID = wellid) %>%
               dplyr::select(WELLID, X_WERT, Y_WERT, Brkenn)


wells <- data.table::fread(input = "wellfield.byn", fill = TRUE)

names(wells) <- c("WELLID", "NODE", "Lay", "Row", "Col", "Totim", "Q_node", "hwell", "hcell", "seepage")

# wells <- wells[,c("NODE", "Lay", "Row", "Col") := NULL]
wells <- wells[,c("NODE", "Lay") := NULL]

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
                 dplyr::left_join(y = entnahme_pro_jahr) %>% 
                  dplyr::filter(X_WERT >= 2529190,
                                X_WERT <=	2530190,
                                Y_WERT <=	5658790,
                                Y_WERT >=	5657790,
                                Totim >=9)

entnahme_prognose <- wells_perYear %>%
  ungroup() %>%
  group_by(Totim) %>% 
  summarise(Gesamtfoerderung = sum(Q_node_cubicmpermin))

Qmom_prognose <- wells_perYear %>%
  ungroup() %>%
  group_by(Totim) %>% 
  summarise(Q_mom = mean(Q_node_cubicmpermin))

#animation
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

#plot(Q_perYear ~ Year, 
#     data=entnahme_pro_jahr, 
#     type="b", 
#     col="blue", 
#     las = 1)

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

# wells_perDay <- wells %>% 
#          #filter(WELLID == "WELL26") %>% 
#          mutate(Totim = round(Totim,0)) %>% 
#          group_by(WELLID, Totim) %>% 
#          summarise(Q_node_median = median(Q_node), 
#                    hwell_median = median(hwell), 
#                    hcell_median = median(hcell))


# wells_summary <- wells_perDay %>% 
#   ungroup() %>% 
#   group_by(WELLID) %>%  
#   filter(Q_node_median < 0) %>%
#   summarise(Q_node_median = -median(Q_node_median), 
#             days_active = n()) %>% 
#   arrange(Q_node_median, desc(days_active))
# 
# wells_summary %>%  
#   ggplot(aes(x = days_active, y = Q_node_median/24/60)) + 
#   geom_point(alpha = 0.5, size = 2) + 
#   labs(y = "Foerderrate pro Minute (m³/min)") +
#   theme_bw()
# 
# # Q sum per well for whole simulation
# wells_perDay %>% 
#   ungroup() %>% 
#   mutate(WELLID = as.numeric(gsub("WELL", "", .$WELLID))) %>%  
#   group_by(WELLID) %>% 
#   summarise(Qsum = -sum(Q_node_median)) %>%  
#   plot(Qsum ~ WELLID, data = ., type="b", col="blue")
# 
# 
# wells_perDay_tidy <- tidyr::gather(data = wells_perDay,Key, Value, -Totim, -WELLID)
# 
# pdf("wells_Q.pdf",height = 7, width = 10)
# xyplot(Q_node_median ~ Totim | WELLID,
#        data = wells_perDay,
#        layout = c(1,1))
# dev.off()
# 
# pdf("wellfield_Q.pdf",height = 7, width = 10)
# xyplot(Q_node_median ~ Totim, group = WELLID,
#        data = wells_perDay)
# dev.off()
# 
# pdf("wells_H.pdf",height = 7, width = 10)
# xyplot(Value ~ Totim | WELLID,
#        group = Key,
# data = wells_perDay_tidy %>% filter(Key !="Q_node_median"),
# layout = c(1,1))
# dev.off()

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

