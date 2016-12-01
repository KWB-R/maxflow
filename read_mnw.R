library(data.table)
library(dplyr)
library(ggplot2)
library(lattice)
#setwd("C:/Users/mrustl/Desktop/WC_Maxflow/branches/one-layer-south")

wells <- data.table::fread(input = "wellfield.byn",skip = 1)

names(wells) <- c("WELLID", "NODE", "Lay", "Row", "Col", "Totim", "Q_node", "hwell", "hcell")

wells <- wells[,c("NODE", "Lay", "Row", "Col") := NULL]

save(wells, file = "wells.RData")


load("wells.RData")



wells_perDay <- wells %>% 
         #filter(WELLID == "WELL26") %>% 
         mutate(Totim = round(Totim,0)) %>% 
         group_by(WELLID, Totim) %>% 
         summarise(Q_node_median = median(Q_node), 
                   hwell_median = median(hwell), 
                   hcell_median = median(hcell))


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

if (FALSE) {
ggplot(wells_perDay, aes(x = Totim, y = Q_node_median)) +
  facet_wrap(~ WELLID) +
  geom_line() +
  theme_bw()


plot(well1$Totim, well1$hcell_median,ylim = c(0,80))
points(well1$Totim, well1$hwell_median)

plot(well1$Totim, well1$Q_node_median)


setkey(wells,Totim)
wells[WELLID == "WELL1", 
      list(Q-node_median=median(Q-node),
           hwell_median = median(hwell),
           hcell_median = median(hcell)),
      by=Totim]

wells$WELLID <- as.numeric(gsub("WELL", "", wells$WELLID))




setkey(wells,Totim)
system.time(wells[,list(Q-node_median=median(Q-node),
                        hwell_median = median(hwell),
                        hcell_median = median(hcell)),
                  by=Totim])


well1 <- wells[WELLID == "WELL1", ]
}
