#JuG on Mon Feb  8 18:26:55 2021
#read RHC dataset

library(tidyverse)
library(utilitR)
library(tableone)

rhc <- read.csv("./data/rhc.csv.gz", header =T, sep=";")
dim(rhc)
table(rhc$swang1)

#select 3000 patients (1500 par groupes)
library(caret)
set.seed(12345)
rhcY<- rhc %>% filter(swang1=="RHC") %>% select(X) %>% slice_sample(.,n=1500)
rhcN<- rhc %>% filter(swang1=="No RHC") %>% select(X) %>% slice_sample(.,n=1500)
indSamp <- sort(c(rhcY$X, rhcN$X))


wRhc <- rhc %>% filter(X %in% indSamp) #working data base (observational)
tabObs <- tableone::CreateTableOne(vars = names( wRhc )[c(-1,-45,-67:-64)], data =  wRhc ,strata = "swang1")
wRhcSMD <- ExtractSmd(tabObs)



rctRhc <-wRhc #working data base (pseudo randomisation groupes)
rctRhc$arm <- sample(x = rep(LETTERS[1:2], each=1500), size = 3000, replace = F)
table(rctRhc$arm)

tabRct <- tableone::CreateTableOne(vars = names( rctRhc )[c(-1,-45,-67:-64, -69)], data =  rctRhc ,strata = "arm")
rctRhcSMD <- ExtractSmd(tabRct)


## Comparaison des SMD obs/rct

dataSMD = data.frame(obs = wRhcSMD[,1],
                     rct = rctRhcSMD[,1]
                    )
rownames(dataSMD) <- rownames(wRhcSMD)

stdDiffDotPlot(dataSMD=dataSMD,col=c("black", "blue"), xlim = c(0,.8),cex=.5 )


