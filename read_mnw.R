library(data.table)
setwd("C:/Users/cmenz/Desktop/WC_Maxflow/branches/mnw2_oneLayer")

well <- readLines("wellfield.byn", skip=1, n = 10)
gsub(pattern = "\"", "", well)

wells <- data.table::fread(input = "wellfield.byn", nrows = 500000000,skip = 1)

names(wells) <- c("WELLID", "NODE", "Lay", "Row", "Col", "Totim", "Q-node", "hwell", "hcell")  

wells <- read.table(file = "wellfield.byn", header = FALSE, sep = " ", skip=1,nrow=10)
