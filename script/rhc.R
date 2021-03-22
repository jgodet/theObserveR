#JuG on Mon Feb  8 18:26:55 2021
#read RHC dataset

library(tidyverse)
library(utilitR)
library(tableone)

rhc <- read.csv("./data/rhc.csv.gz", header =T, sep=";")
dim(rhc)
table(rhc$swang1)

#select 3000 patients (1500 par groupe)
library(caret)
set.seed(12345)
rhcY<- rhc %>% filter(swang1=="RHC") %>% select(X) %>% slice_sample(.,n=1500)
rhcN<- rhc %>% filter(swang1=="No RHC") %>% select(X) %>% slice_sample(.,n=1500)
indSamp <- sort(c(rhcY$X, rhcN$X))


wRhc <- rhc %>% filter(X %in% indSamp) #working data base (observational group)
tabObs <- tableone::CreateTableOne(vars = names( wRhc )[c(-1,-45,-67:-64)], data =  wRhc ,strata = "swang1",smd=TRUE)
wRhcSMD <- ExtractSmd(tabObs)



rctRhc <-wRhc #rct data base (pseudo randomisation group)
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

#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6351359/


############################################################################################################################################
################################################ Comparaison SMD obs: matching/weighting ###################################################
############################################################################################################################################
library(tableone)
library(Matching)
library(survey)
library(reshape2)
library(ggplot2)

attach(rhc)
rhc$swang01<-ifelse(swang1=="RHC",1,0)

vars<-c("cat1", "ca", "sadmdte", "dschdte", 
        "lstctdte", "death", "cardiohx", "chfhx", "dementhx", "psychhx", 
        "chrpulhx", "renalhx", "liverhx", "gibledhx", "malighx", "immunhx", 
        "transhx", "amihx", "age", "sex", "edu", "surv2md1", "das2d3pc", 
        "t3d30", "dth30", "aps1", "scoma1", "meanbp1", "wblc1", "hrt1", 
        "resp1", "temp1", "pafi1", "alb1", "hema1", "bili1", "crea1", 
        "sod1", "pot1", "paco21", "ph1", "dnr1", 
        "ninsclas", "resp", "card", "neuro", "gastr", "renal", "meta", 
        "hema", "seps", "trauma", "ortho", "race", 
        "income")#, "ptid", "pRhc", "pNoRhc", "pAssign", "pMin", "recul")

##### Table 1 (original data)
tabUnmatched <- CreateTableOne(vars = vars, strata = "swang1", data = rhc, test = TRUE)
print(tabUnmatched, nonnormal=c("age","aps"), smd = TRUE)

##### Propensity score (ps) estimation
rhc$ID<-c(1:length(cat1))
rhct<-data.frame(ID,cat1, ca, sadmdte, dschdte, 
                 lstctdte, death, cardiohx, chfhx, dementhx, psychhx, 
                 chrpulhx, renalhx, liverhx, gibledhx, malighx, immunhx, 
                 transhx, amihx, age, sex, edu, surv2md1, das2d3pc, 
                 t3d30, dth30, aps1, scoma1, meanbp1, wblc1, hrt1, 
                 resp1, temp1, pafi1, alb1, hema1, bili1, crea1, 
                 sod1, pot1, paco21, ph1, swang1, swang01, dnr1, 
                 ninsclas, resp, card, neuro, gastr, renal, meta, 
                 hema, seps, trauma, ortho, race, 
                 income)
rhct2<-rhct[complete.cases(rhct),]						

psModel <- glm(formula = swang01~cat1+ ca+ sadmdte+ dschdte+ 
               lstctdte+ death+ cardiohx+ chfhx+ dementhx+ psychhx+ 
               chrpulhx+ renalhx+ liverhx+ gibledhx+ malighx+ immunhx+ 
               transhx+ amihx+ age+ sex+ edu+ surv2md1+ das2d3pc+ 
               t3d30+ dth30+ aps1+ scoma1+ meanbp1+ wblc1+ hrt1+ 
               resp1+ temp1+ pafi1+ alb1+ hema1+ bili1+ crea1+ 
               sod1+ pot1+ paco21+ ph1+ dnr1+ 
               ninsclas+ resp+ card+ neuro+ gastr+ renal+ meta+ 
               hema+ seps+ trauma+ ortho+ race+ 
               income,
               family  = binomial(link = "logit"),
               data    = rhct2)

## Predicted probability of being assigned to treatment (ps)
pRHC <- predict(psModel, type = "response")
## Predicted probability of being assigned to no treatment
pNoRHC <- 1 - pRHC

## Predicted probability of being assigned to the
## treatment actually assigned (either treatment or no treatment)
pAssign <- NA
pAssign[rhct2$swang1 == "RHC"] <- pRHC[rhct2$swang1 == "RHC"]
pAssign[rhct2$swang1 == "No RHC"] <- pNoRHC[rhct2$swang1 == "No RHC"]

## Smaller of pRHC vs pNoRHC for matching weight
pMin <- pmin(pRHC, pNoRHC)

## Create weights (for ATE)
weight <- ifelse(rhct2$swang1=="RHC",1/pRHC,1/(1-pRHC))
PrRHC<-sum(rhct2$swang1=="RHC")/length(rhct2$swang1)
PrNoRHC<-sum(rhct2$swang1=="No RHC")/length(rhct2$swang1)
stweight <- ifelse(rhct2$swang1=="RHC",PrRHC/pRHC,PrNoRHC/(1-pRHC)) # à vérifier


#####################################################################################
##################### Propensity score (ps) matching ################################
#####################################################################################
set.seed(1664)
listMatch <- Match(Tr       = (rhct2$swang01 == 1),      # Need to be in 0,1
                   ## logit of PS,i.e., log(PS/(1-PS)) as matching scale
                   X        = log(pRHC / pNoRHC),
                   ## 1:1 matching
                   M        = 1,
                   ## caliper = 0.1 * SD(logit(PS))
                   caliper  = 0.1,
                   replace  = FALSE,
                   ties     = TRUE,
                   version  = "fast")
## Extract matched data
rhcMatched <- rhct2[unlist(listMatch[c("index.treated","index.control")]), ]
rhct2$ID[unique(c(listMatch$index.treated,listMatch$index.control))]
rhcMatchedb <- rhc[rhct2$ID[unique(c(listMatch$index.treated,listMatch$index.control))], ]

## Construct a table
tabMatched <- CreateTableOne(vars = vars, strata = "swang1", data = rhcMatchedb, test = TRUE)
## Show table with SMD
print(tabMatched, nonnormal=c("age","aps"), smd = TRUE)


#####################################################################################
######################################## IPTW #######################################
#####################################################################################
##  weight
weight <- ifelse(rhct2$swang1=="RHC",1/pRHC,1/(1-pRHC))
## Weighted data
rhct3<-rhc[rhct2$ID,]
rhcSvy <- svydesign(ids = ~ 1, data = rhct3, weights = ~ weight)

## Construct a table (This is a bit slow.)
tabWeighted <- svyCreateTableOne(vars = vars, strata = "swang1", data = rhcSvy, test = TRUE)
## Show table with SMD
print(tabWeighted, smd = TRUE)



#####################################################################################
###################################### IPTWst #######################################
#####################################################################################
##  stabilized weight
stweight <- ifelse(rhct2$swang1=="RHC",PrRHC/pRHC,PrNoRHC/(1-pRHC)) # à vérifier
## Weighted data
rhct3<-rhc[rhct2$ID,]
rhcSvy <- svydesign(ids = ~ 1, data = rhct3, weights = ~ stweight)

## Construct a table (This is a bit slow.)
tabWeightedst <- svyCreateTableOne(vars = vars, strata = "swang1", data = rhcSvy, test = TRUE)
## Show table with SMD
print(tabWeightedst, smd = TRUE)


#####################################################################################
##################### Propensity score matching weight ##############################
#####################################################################################
## matching weight
mw <- pMin / pAssign
## Weighted data
rhct3<-rhc[rhct2$ID,]
rhcSvy <- svydesign(ids = ~ 1, data = rhct3, weights = ~ mw)

## Construct a table (This is a bit slow.)
tabWeightedm <- svyCreateTableOne(vars = vars, strata = "swang1", data = rhcSvy, test = TRUE)
## Show table with SMD
print(tabWeightedm, smd = TRUE)



#####################################################################################
###################### Propensity score overlap weight ##############################
#####################################################################################
## Overlap weight
ow <- (pAssign * (1 - pAssign)) / pAssign
## Weighted data
rhcSvyOw <- svydesign(ids = ~ 1, data = rhct3, weights = ~ ow)

## Construct a table (This is a bit slow.)
tabWeightedOw <- svyCreateTableOne(vars = vars, strata = "swang1", data = rhcSvyOw, test = TRUE)
## Show table with SMD
print(tabWeightedOw, smd = TRUE)


#####################################################################################
############ Assessing balance before and after matching/weighting ##################
#####################################################################################
## Construct a data frame containing variable name and SMD from all methods
dataPlot <- data.frame(variable  = rownames(ExtractSmd(tabUnmatched)),
                       Unmatched = as.numeric(ExtractSmd(tabUnmatched)),
                       Matched   = as.numeric(ExtractSmd(tabMatched)),
                       #ExactMatched = as.numeric(ExtractSmd(tabExMatched)),
                       Weighted = as.numeric(ExtractSmd(tabWeighted)),
                       Weightedst = as.numeric(ExtractSmd(tabWeightedst)),
                       WeightedM  = as.numeric(ExtractSmd(tabWeightedm)),
                       WeightedOw = as.numeric(ExtractSmd(tabWeightedOw)))

## Create long-format data for ggplot2
dataPlotMelt <- melt(data          = dataPlot,
                     id.vars       = c("variable"),
                     variable.name = "Method",
                     value.name    = "SMD")

## Order variable names by magnitude of SMD
varNames <- as.character(dataPlot$variable)[order(dataPlot$Unmatched)]

## Order factor levels in the same order
dataPlotMelt$variable <- factor(dataPlotMelt$variable,
                                levels = varNames)

## Plot using ggplot2
ggplot(data = dataPlotMelt,
       mapping = aes(x = variable, y = SMD, group = Method, color = Method)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.1, color = "black", size = 0.1) +
  coord_flip() +
  theme_bw() + theme(legend.key = element_blank())




