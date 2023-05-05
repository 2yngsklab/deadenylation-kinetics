library(deSolve) #library for solving differential equations
library(minpack.lm) #library for least squares fit using levenberg-marquart algorithm
library(tidyverse)
library(viridis)

GetPreData <- function(rawdata, invert = FALSE, max.pixel){
        preD <- rawdata %>% 
                rename(pixel = `(pixel)`) %>%
                gather(lane,value,-pixel) %>% 
                filter(!is.na(value)) %>%
                filter(pixel < max.pixel)
        
        if(invert){
          preD <- preD %>% 
                  mutate(value = max(value)-value)
        }
        return (preD)
}

GetRefData <- function(predata){
        refD <- predata %>% 
                group_by(pixel) %>% 
                summarise(value = sum(value))
        return (refD)
}

GetTroughs <- function(refdata, expected.num.troughs, remove.list, add.list){
        auto.troughs  <- which(diff(sign(diff(-refdata$value)))==-2)+1
        troughdat     <- refdata %>% 
                       filter(pixel %in% refdata$pixel[auto.troughs])
        #remove false positives
        if(length(remove.list) > 0){
            troughs.adjusted <- auto.troughs[-remove.list]
        }
        else{
            troughs.adjusted <- auto.troughs
        }
        #add false negatives
        troughs.adjusted <- c(add.list,troughs.adjusted)

        if(length(troughs.adjusted) != expected.num.troughs){
            warning(paste("number of troughs detected is:",length(troughs.adjusted)))
        }
        return (troughs.adjusted)
}

GetMetaPlot <- function(refdata, expected.num.troughs){
        auto.troughs  <- which(diff(sign(diff(-refdata$value)))==-2)+1
        troughdat     <- refdata %>% 
                       filter(pixel %in% refdata$pixel[auto.troughs])
        if(length(auto.troughs) != expected.num.troughs){
            warning(paste("number of troughs detected is:",length(auto.troughs)))
        }
        metaP <- refdata %>%
                 ggplot(aes(pixel,value)) + 
                 geom_line() +
                 theme_gray(base_size=12)+
                 labs(y="Intensity value (a.u.)", x = "Pixel position") +
                 geom_vline(data=troughdat,aes(xintercept = pixel),col='red') 
        return (metaP)
}

GetLanePlot <- function(predata, refdata, troughs, lane.no = c(), add.list){
        for (elem in add.list) {
          troughs <- troughs[troughs != elem]
        }
        troughD <- refdata %>% 
                 filter(pixel %in% refdata$pixel[troughs])
        addTrough <- refdata %>%
                 filter(pixel %in% refdata$pixel[add.list])
        laneP <- predata %>% 
                 group_by(pixel) %>% 
                 filter(lane %in% lane.no) %>%
                 ggplot(aes(pixel,value)) +
                 theme_gray(base_size=12)+
                 labs(y="Intensity value (a.u.)", x="Pixel position") +
                 geom_line() +
                 geom_vline(data=troughD,aes(xintercept = pixel),col='red') +
                 geom_vline(data = addTrough, aes(xintercept = pixel), col = 'red', linetype = "dashed") +
                 facet_wrap(~lane, nrow = length(lane.no))
        return (laneP)
}

GetDgtData <- function(predata, refdata, troughs, max.pixel, time.points){
        bat <- predata %>% 
               ungroup() %>%
               mutate(bin = cut(pixel,c(-1,refdata$pixel[troughs],max.pixel))) %>% 
               group_by(lane,bin) %>%
               summarise(tvalue = max(value)) %>% #discretize
               ungroup() %>%
               mutate(tvalue = (tvalue - min(tvalue))/(max(tvalue) - min(tvalue)) ) %>% #unity-based normalization
               group_by(lane) %>%
               mutate(value = tvalue) %>%
               ungroup()
        
        dgtD <- bat %>% 
                ungroup() %>%
                mutate(time = timepoints[as.numeric(as.character(lane))]) %>%
                mutate(species = as.numeric(bin)) %>% 
                select(time,species,value) %>%
                mutate(value = ifelse(time == 0 & species == 2, 0.008675228, value)) %>% #minor correction
                spread(species,value)
        return (dgtD)
}

GetHeatMap <- function(dgtdata, num.troughs, gtitle){
        heatM <-  as.data.frame(dgtdata) %>%
                  gather(species,value,-time) %>%
                  mutate(value = as.numeric(value)) %>% 
                  mutate(species = factor(species,levels = as.character(rev(seq(1:(num.troughs+1)))))) %>%
                  filter(time >= 0) %>%
                  group_by(time) %>%
                  mutate(value = (value - min(value))/(max(value) - min(value)) ) %>% #unity-based normalization for visualization
                  ggplot(aes(factor(time),species,fill=value)) + geom_tile(col=NA) + 
                  ggtitle(gtitle) + 
                  theme(plot.title = element_text(hjust = 0.5)) +
                  theme(legend.title = element_blank())  +
                  scale_fill_viridis(option = "D") +
                  xlab("Reaction time (min)") +
                  theme(axis.title.y = element_blank())
        return (heatM)
        
}

LaneSpecificCorrection <- function(rawdata, lane.no, correction){
        rawdata <- rawdata %>% 
                   gather(lane,value,-`(pixel)`) %>% 
                   mutate(`(pixel)` = ifelse(lane == lane.no,`(pixel)`+correction, `(pixel)`)) %>% 
                   spread(lane,value,fill = 0)
}

##https://www.r-bloggers.com/learning-r-parameter-fitting-for-models-involving-differential-equations/

rxnrate=function(t,y,parms){
        #parms: vector of L parameters where L is the length of poly(A) tail
        #y is the concentration of RNA species of length L
        #derivatives dy/dt are computed below
        r=rep(0,length(y))
        r[1]=-parms[1]*y[1] #dyA/dt
        for(i in 2:length(y)){
            r[i]=parms[i-1]*y[i-1]-parms[i]*y[i] #dyB/dt
        } 
        r[length(y)]=parms[length(y)-1]*y[length(y)-1] #dyC/dt
        #the computed derivatives are returned as a list
        #order of derivatives needs to be the same as the order of species in y
        return(list(r))
}

simulation <- function(parms, l, t){
        #initialize concentration of each species
        cinit=rep(0,l)
        cinit[1]=1    
        #solve ODE for a given set of parameters
        estimates=ode(y=cinit,times=t,func=rxnrate,parms=parms)
        return(estimates)
}

ssq = function(parms,dgtdata){
        l = ncol(dgtdata) - 1 #number of species
        t = dgtdata$time
        sim.est <- simulation(parms, l, t)    
  
        #calculate residual
        preddf <- as.data.frame(sim.est) %>% #predicted concentration
                  gather(species,value,-time) %>%
                  mutate(species = factor(species, levels=rev(seq(1:l)))) %>% 
                  group_by(time) %>% 
                  mutate(cum.value = cumsum(value)) %>% 
                  ungroup()
  
        expdf <-  as.data.frame(dgtdata) %>% #experimental data concentration 
                  gather(species,value,-time) %>%
                  mutate(species = factor(species, levels=rev(seq(1:l)))) %>% 
                  group_by(time) %>% 
                  mutate(cum.value = cumsum(value)) %>% 
                  ungroup()
  
        combdf <- preddf %>% 
                  inner_join(expdf, by=c("time","species")) %>%
                  mutate(ssqres= value.x - value.y)
        
        ssqres=combdf$ssqres
  
  #return predicted vs experimental residual
  return(ssqres)
}

GetStepPlot <- function(fitval1, num.troughs, gtitle, height){
        data.frame(est=coef(fitval1), se=sqrt(diag(vcov(fitval1))), position=seq(1,num.troughs,1),type=gtitle) %>%
        filter(position < num.troughs+1) %>% 
        ggplot(aes(position,est)) + 
        geom_errorbar(aes(ymin=est-se, ymax=est+se),width=0.5) +
        geom_step(aes(x=position),col="dark red",direction = 'vh') +
        geom_point() + 
        theme_minimal() +
        scale_y_continuous(name = "Deadenylation kinetics (nt / min)", expand = c(0, 0), limits = c(0,height)) +
        scale_x_continuous(name = expression("Single-nucleotide position"), breaks = c(1,7,14,21,28), minor_breaks = NULL) +
        theme_minimal(base_size = 12) +
        theme(strip.text = element_text(color = "black"), axis.text = element_text(color = "black")) +
        theme(axis.ticks = element_line(linetype = "solid",color="black",size = 0.5), axis.line = element_line(linetype = "solid",color="black",size = 0.5))
}




