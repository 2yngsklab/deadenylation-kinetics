source("~/deadenylation.R")

filename = './converted.txt'

gtitle     <- "Example"
rawdata    <- read_tsv(filename,skip = 5)
maxpixel   <- 750
num.troughs  <- 24
timepoints <- c(-8,0,2,4,6,8,12,16,24,32,48)


#MetaPlot
predata  <- GetPreData(rawdata, max.pixel = maxpixel)
refdata  <- GetRefData(predata)
metaplot <- GetMetaPlot(refdata, num.troughs)
options(repr.plot.width=8, repr.plot.height=1.5)
metaplot

remove.list <- c()
add.list    <- c(26)


#LanePlot
troughs    <- GetTroughs(refdata, expected.num.troughs = num.troughs, remove.list, add.list)
laneplot <- GetLanePlot(predata, refdata, troughs, lane.no = c(3,4), add.list) #enter which lanes to plot
options(repr.plot.width=8, repr.plot.height=3)
laneplot


#HeatMap
dgtdata <- GetDgtData(predata, refdata, troughs, max.pixel = maxpixel, time.points = timepoints)
heatmap <- GetHeatMap(dgtdata, num.troughs, gtitle)
options(repr.plot.width=3.5, repr.plot.height=4.5)
heatmap


#Deadenylation Rate
parms <- rep(0.8,num.troughs) #initialization
dgtdata.tidy <- dgtdata %>% 
                gather(species,value,-time) %>%
                filter(time >= 0) %>% #remove marker lane
                spread(species,value)
fitval1 <- nls.lm(par=parms, fn=ssq, dgtdata=dgtdata.tidy, control = nls.lm.control(maxiter=35)) #parameter estimation
summary(fitval1)
stepplot <- GetStepPlot(fitval1, num.troughs, gtitle, height = 3.5)
options(repr.plot.width=3.5, repr.plot.height=3)
stepplot

