source("~/deadenylation.R")
#Simul1 c(0.8, 1, 1.2, 1.5, 1, 0.8, 0.5, 0.5, 0.8, 1, 1.3, 1.5, 1.3, 1, 0.8, 0.8, 0.5, 0.3, 0.3, 0.3)
#Simul2 c(0.3, 0.5, 0.8, 1, 0.8, 0.5, 0.3, 0.5, 0.8, 1, 0.8, 1, 1, 0.8, 0.5, 0.5, 0.3, 0.3, 0.5, 0.8)
parms <- c(0.8, 1, 1.2, 1.5, 1, 0.8, 0.5, 0.5, 0.8, 1, 1.3, 1.5, 1.3, 1, 0.8, 0.8, 0.5, 0.3, 0.3, 0.3)
num.species = 20
heatmap.data <- simulation(parms = parms, l = num.species, t = c(0,8,16,24,32,40,48))
heatmap <- GetHeatMap(heatmap.data, num.species, gtitle = element_blank())
heatmap

