
# Logit Method ------------------------------------------------------------

setwd("C:/R")


# Packages
library(tidyverse)
library(margins)
library(car)
library(pscl)
library(lmtest)
library(lme4)
library(magrittr)
library(corrplot)
library(MatchIt)


# Load the data
IL_ready_data <- read.csv("IL_ready_data.csv", encoding = "latin1")


# Filter for year <= 2011
IL_ready_data %<>% 
  filter(year <= 2011)

# select the variables of interest 
IL_ready_data <- IL_ready_data %>%
  mutate(
    forest = scale(forest),
    def_perc = scale(def_perc),
    def_lag_perc = scale(def_lag_perc),
    buf_def_perc = scale(buf_def_perc),
    buf_def_lag_perc = scale(buf_def_lag_perc),
    towns = scale(towns),
    ports = scale(ports),
    slope = scale(slope),
    mm_Aug_Sep = scale(mm_Aug_Sep),
    temp_Aug_Sep = scale(temp_Aug_Sep)
  )


# Treat year as a factor
IL_ready_data$year <- as.factor(IL_ready_data$year)


# Construct the model wih glmer
model <- glmer(treat ~ 
                 def_perc + def_lag_perc + buf_def_perc + buf_def_lag_perc + ports + soil + slope + temp_Aug_Sep + (1 | year), 
               data = IL_ready_data, family = binomial(link = "logit"))


# View the results
summary(model)

# Verify the marginal contribuition of each variable
summary(margins(model))


# Analyses the VIF
vif(model)








# ###################################################### ------------------


# Results for Propensity Score Weighting and Panel Match for the R --------


# AF With Negatives  ------------------------------------------------------
#database

set.seed(1)

Sys.setlocale(category = "LC_ALL", locale = "Portuguese") #for encoding

rm(list=ls()) 

#change file location
setwd("C:/R/replication/replication/")

require(foreign)
s <- read.dbf("data/GEOFT_TERRA_INDIGENA.dbf")
table(s$fase_ti)

d <- read.csv("data/Indigenous_land_Y1985_2021.csv")[,-c(1,70,71)]
d[0,]
Encoding(d$NAME) <-"UTF-8"
Encoding(d$ORIG_NAME) <-"UTF-8"
colnames(d)[1:37] <- 1985:2021
d <- d[,c(66,1:37,42,52,61,62)]

require(reshape2)
d_long <- melt(d, id=c("WDPAID","GIS_AREA","NAME","STATUS","STATUS_YR"))
colnames(d_long) <- c("ID","GEE_km2","name","status","status_yr","year","forest")
d <- d_long
rm(d_long)
d$GEE_km2 <- round(d$GEE_km2, digits = 2)
d$forest <- round(d$forest, digits = 4) #in km2
d[0,]

buf <- read.csv("data/Indigenous_land_buffer_Y1985_2021.csv")[,-c(1,39,75,77)]
buf[0,]
Encoding(buf$NAME) <-"UTF-8"
Encoding(buf$ORIG_NAME) <-"UTF-8"
colnames(buf)[1:37] <- 1985:2021
buf <- buf[,c(70,1:37,42,52,63,64)]
buf_long <- melt(buf, id=c("WDPAID","GIS_AREA","NAME","STATUS","STATUS_YR"))
colnames(buf_long) <- c("ID","GEE_km2","name","status","status_yr","year","forest")
buf <- buf_long
rm(buf_long)
buf$GEE_km2 <- round(buf$GEE_km2, digits = 2)
buf$forest <- round(buf$forest, digits = 4)
buf[0,]

d$ID2 <- paste(d$ID, d$year)
buf$ID2 <- paste(buf$ID, buf$year)
d$buf_forest <- buf[match(with(d, ID2), with(buf, ID2)),]$forest 
rm(buf)
table(d$buf_forest > d$forest)
d$buf_forest2 <- d$buf_forest - d$forest 
d$buf_forest <- d$buf_forest2
d <- d[,-ncol(d)]

require(dplyr)
require(Hmisc)
list_ID <- unique(d$ID)
d_new <- d[0,]
for (i in 1:length(list_ID)) {
  temp <- subset(d, ID==list_ID[i])
  temp$forest_lag <- Lag(temp$forest, 1)
  temp$buf_forest_lag <- Lag(temp$buf_forest, 1)
  d_new <- rbind(d_new, temp)
}
d <- d_new
rm(d_new, temp)
d$year <- as.integer(as.character(d$year))
d$def <- round((d$forest_lag - d$forest) * 100, 2) #in hectares
d$buf_def <- round((d$buf_forest_lag - d$buf_forest) * 100, 2) #in hectares

d <- subset(d, year > 1985)
#remove lines
#d$def[ d$def < 0] <- 0
#d$buf_def[ d$buf_def < 0] <- 0


list_ID <- unique(d$ID)
d_new <- d[0,]
for (i in 1:length(list_ID)) {
  temp <- subset(d, ID==list_ID[i])
  temp$def_lag <- Lag(temp$def, 1)
  temp$buf_lag <- Lag(temp$buf_def, 1)
  d_new <- rbind(d_new, temp)
}
length(unique(d_new$ID))
d <- d_new
rm(d_new, temp)
d <- subset(d, year > 1985)
d$def_lag[ is.na(d$def_lag) ] <- 0
d$buf_lag[ is.na(d$buf_lag) ] <- 0

require(foreign)
cov <- read.dbf("data/IL_covariates_final.dbf")
d$capitals <- round(cov[match(with(d, ID), with(cov, WDPAID)),]$capitals, digits = 2) 
d$towns <- round(cov[match(with(d, ID), with(cov, WDPAID)),]$towns, digits = 2)  
d$ports <- round(cov[match(with(d, ID), with(cov, WDPAID)),]$ports, digits = 2) 
d$soil <- round(cov[match(with(d, ID), with(cov, WDPAID)),]$soil, digits = 2) 
d$dem <- round(cov[match(with(d, ID), with(cov, WDPAID)),]$dem, digits = 2) 
d$slope <- round(cov[match(with(d, ID), with(cov, WDPAID)),]$slope, digits = 2) 
d$fric <- round(cov[match(with(d, ID), with(cov, WDPAID)),]$fric, digits = 4) 
d$biome <- cov[match(with(d, ID), with(cov, WDPAID)),]$biome
d$region <- cov[match(with(d, ID), with(cov, WDPAID)),]$region
rm(cov)

cov <- read.dbf("data/IL_all.dbf")
d$IL_ha <- cov[match(with(d, ID), with(cov, WDPAID)),]$poly_ha 
cov <- read.dbf("data/IL_all_buffer.dbf")
d$buffer_ha <- cov[match(with(d, ID), with(cov, WDPA_PID)),]$poly_ha 

cov <- read.csv("data/annual_precipit_IL_1981_2022.csv")
cov$ID2 <- paste(paste(cov$WDPA_PID, cov$year))
d$mm_annual <- cov[match(with(d, ID2), with(cov, ID2)),]$mm_annual
cov <- read.csv("data/precipt_Aug_Sep_IL_1981_2022.csv")
cov$ID2 <- paste(paste(cov$WDPA_PID, cov$year))
d$mm_Aug_Sep <- cov[match(with(d, ID2), with(cov, ID2)),]$mm_Aug_Sep
cov <- read.csv("data/max_temp_Aug_Sep_IL_1981_2022.csv")
cov$ID2 <- paste(paste(cov$WDPA_PID, cov$year))
d$temp_Aug_Sep <- cov[match(with(d, ID2), with(cov, ID2)),]$temp_Aug_Sep
rm(cov)

str(d)
d$biome <- as.character(d$biome)
d$biome[ d$biome=="AMAZÔNIA"] <- "Amazon"
d$biome[ d$biome=="CERRADO"] <- "Cerrado"
d$biome[ d$biome=="CAATINGA"] <- "Caatinga"
d$biome[ d$biome=="MATA ATLÂNTICA"] <- "Atlantic Forest"
d$biome[ d$biome=="PAMPA"] <- "Pampa"
d$biome[ d$biome=="PANTANAL"] <- "Pantanal"

d$region <- as.character(d$region)
d$region[ d$region=="Norte"] <- "N"
d$region[ d$region=="Nordeste"] <- "NE"
d$region[ d$region=="Centro-Oeste"] <- "CW"
d$region[ d$region=="Sul"] <- "S"
d$region[ d$region=="Sudeste"] <- "SE"

d$def_perc <- round(d$def/d$IL_ha * 100, digits = 3)
d$def_lag_perc <- round(d$def_lag/d$IL_ha * 100, digits = 3)
d$buf_def_perc <- round(d$buf_def/d$buffer_ha * 100, digits = 3)
d$buf_def_lag_perc <- round(d$buf_lag/d$buffer_ha * 100, digits = 3)

d$Amazon <- ifelse(d$biome=="Amazon", 1, 0)
d$Cerrado <- ifelse(d$biome=="Cerrado", 1, 0)
d$Caatinga <- ifelse(d$biome=="Caatinga", 1, 0)
d$Atlantic <- ifelse(d$biome=="Atlantic Forest", 1, 0)
d$Pampa <- ifelse(d$biome=="Pampa", 1, 0)
d$Pantanal <- ifelse(d$biome=="Pantanal", 1, 0)


require(PanelMatch)
require(ggplot2)

length(unique(d$ID))
length(unique(d$year))
table(d$status)
table(d$status_yr)

d <- subset(d, status_yr != 0 & d$IL_ha >= 100)
length(unique(d$ID))
table(d$biome)
13788/24588
5832/24588
1116/24588
3348/24588
396/24588
108/24588

d$treat <- ifelse(d$status=='Designated' & d$status_yr < d$year, 1, 0)
table(d$treat)

#change file name
write.csv(d, "data/AF_NEG_NEG_IL_ready_data.csv", row.names = F)
d$ID <- as.integer(as.character(d$ID))

#analysis
#filter atlantic biome
#keep only the psw method
#remove biome columns "'Amazon', 'Cerrado', 'Caatinga', 'Atlantic', 'Pampa'"
#from 21 columns to 16 

d <- d %>% filter(biome == "Atlantic Forest")

res <- as.data.frame(matrix(NA,0,7))
colnames(res) <- c('estimate','std.error','2.5%','97.5%','lag','match_size','lead')
res_cov <- as.data.frame(matrix(NA,0,16))
colnames(res_cov) <- c('def_lag', 'buf_def', 'soil','capitals','towns','dem','slope','fric','IL_ha',
                       'mm_annual', 'mm_Aug_Sep', 'temp_Aug_Sep',
                       'lag','match_size','lead','time')
control_list <- as.data.frame(matrix(NA,0,1))
colnames(control_list)[1] <- "weigths"
control_list$ID_match <- as.integer(as.vector(rownames(control_list)))

lags_n <- 5 
size_match_n <- 1
lead_upper_range <- c(1,3,5,7,9)

for (l in 1:length(lags_n)) {
  for (m in 1:length(size_match_n)) {
    for (u in 1:length(lead_upper_range)) {
      
      PM.results <- PanelMatch(lag = lags_n[l],
                               time.id = "year", unit.id = "ID", 
                               treatment = "treat", 
                               refinement.method = "ps.weight",
                               data = d, match.missing = TRUE, 
                               
                               covs.formula = ~def_lag_perc + buf_def_perc + soil + capitals + towns + dem + slope + fric + IL_ha +
                                 mm_annual + mm_Aug_Sep + temp_Aug_Sep,
                               
                               size.match = size_match_n[m],
                               qoi = "att" ,
                               outcome.var = "def_perc",
                               lead = 0:lead_upper_range[u],
                               forbid.treatment.reversal = TRUE, 
                               use.diagonal.variance.matrix = TRUE)
      
      used <- lapply(PM.results$att, attributes) 
      
      for (i in 1:length(used)) {
        temp <- as.data.frame(used[[i]]$weights)
        colnames(temp)[] <- "weigths"
        temp <- subset(temp, weigths > 0)
        temp$ID_match <- as.integer(as.vector(rownames(temp)))
        control_list <- rbind(control_list, temp)
      }
      
      cov <- get_covariate_balance(PM.results$att,
                                   data = d,
                                   covariates = c('def_lag', 'buf_def_perc', 'soil','capitals','towns','dem','slope','fric','IL_ha',    
                                                  'mm_annual', 'mm_Aug_Sep', 'temp_Aug_Sep'),
                                   plot = F, legend = F)
      cov <- as.data.frame(cov)
      cov$lag <- lags_n[l]
      cov$match_size <- size_match_n[m]
      cov$lead <-lead_upper_range[u]
      cov$time <- row.names(cov)
      res_cov <- rbind(res_cov, cov)
      
      
      PE.results <- PanelEstimate(sets = PM.results, data = d,
                                  number.iterations = 1000,
                                  df.adjustment = T,
                                  confidence.level = 0.95)
      
      temp <- summary(PE.results)
      temp <- as.data.frame(temp$summary)
      temp$lag <- lags_n[l]
      temp$match_size <- size_match_n[m]
      temp$lead <-lead_upper_range[u]
      temp$time <- row.names(temp)
      
      res <- rbind(res, temp)
    }
  }
}

#change file name
write.csv(res_cov, "data/AF_NEG_NEG_cov_balance_PSW_indigenous_lands.csv", row.names = F)
write.csv(res, "data/AF_NEG_ATT_PSW_indigenous_lands.csv", row.names = F)

#covariate balance for the indigenous lands in Atlantic Forest 

require(ggplot2)
require(reshape)
require(plyr)

#keep only the psw method
#from all biomes to Atlantic Forest

d <- read.csv('data/AF_NEG_cov_balance_PSW_indigenous_lands.csv')
d <- melt(d, id=c('lag','lead','match_size','time'))
PSW <- ddply(.data=d, .(lag,match_size,lead,variable), .fun=summarise, mean = mean(value))
PSW$Method <- 'B'

#keep only the psw

all <- rbind(PSW)
all$biome <- "Atlantic Forest"

d_plot <- all

d_plot$lead <- d_plot$lead + 1

unique(d_plot$Method)

unique(d_plot$variable)
d_plot$variable <- as.character(d_plot$variable)

#keep only atlantic forest variable

d_plot$variable[ d_plot$variable == 'IL_ha' ] <- 'Size'
d_plot$variable[ d_plot$variable == 'soil' ] <- 'Soil quality'
d_plot$variable[ d_plot$variable == 'capitals' ] <- 'Dist. to capitals'
d_plot$variable[ d_plot$variable == 'towns' ] <- 'Dist. to towns'
d_plot$variable[ d_plot$variable == 'capitals' ] <- 'Dist. capitals'
d_plot$variable[ d_plot$variable == 'dem' ] <- 'Elevation'
d_plot$variable[ d_plot$variable == 'slope' ] <- 'Slope'
d_plot$variable[ d_plot$variable == 'fric' ] <- 'Dist. to cities'
d_plot$variable[ d_plot$variable == 'buf_def_perc' ] <- ' 10-km buffer def.'
d_plot$variable[ d_plot$variable == 'def_lag' ] <- 'Lagged def.'
d_plot$variable[ d_plot$variable == 'mm_annual' ] <- 'Annual precip.'
d_plot$variable[ d_plot$variable == 'mm_Aug_Sep' ] <- 'Aug-Sep. precip.'
d_plot$variable[ d_plot$variable == 'temp_Aug_Sep' ] <- 'Aug-Sep. temp.'
d_plot$variable[ d_plot$variable == 'Atlantic' ] <- '  Atlantic forest biome'

unique(d_plot$lag)
unique(d_plot$match_size)
unique(d_plot$lead)

#from all biomes to Atlantic Forest

g1 <- ggplot(data=subset(d_plot, biome=="Atlantic Forest"),
             aes(x=mean, y=variable, group=as.factor(lead), colour=as.factor(lead), shape=as.factor(lead)))  +
  guides(fill=guide_legend(reverse=TRUE)) +  
  theme_bw() +
  theme(legend.position = 'bottom',
        axis.title.y=element_blank(),
        panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(),  panel.grid.minor.x = element_blank()) +
  facet_wrap(. ~ Method, scales = "free", ncol = 2) +
  geom_vline(xintercept=0, linetype='dashed') + geom_point(size=1.7, stroke = 1.1) +
  scale_shape_manual(values=c(0,1,2,3,4)) +
  labs(colour = "Evaluation period (years)") + labs(shape = "Evaluation period (years)") +
  xlab("Standardized mean difference") + scale_x_continuous(breaks = c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3), limits = c(-0.35,0.35)) +
  annotate("rect", xmin = -0.1, xmax = 0.1,
           ymin = 0, ymax = c(as.numeric(length(unique(d_plot$variable)))+1),alpha = .1, fill = "dodgerblue")

g1
ggsave(g1, units="cm", width=20, height=20, dpi=600, file="figures/AF_NEG_IL_cov_balance_ALL.svg")
#change file name


#plot

require(plyr)
require(ggplot2)

rm(list=ls()) 

#from all biomes to Atlantic Forest
#keep only the p3

p3 <- read.csv("data/AF_NEG_ATT_PSW_indigenous_lands.csv")
p3$method <- "B"
p3$biome <- "Atlantic Forest"

d_plot <- rbind(p3)

d_plot$estimate_minus_SE <- d_plot$estimate + d_plot$std.error
sum(d_plot$estimate_minus_SE < 0)/nrow(d_plot)
ddply(.data=subset(d_plot, estimate_minus_SE < 0), .(time), .fun=summarise, def = mean(estimate))

d_plot$year <- ifelse(d_plot$time=="t+0", 1,
                      ifelse(d_plot$time=="t+1", 2,
                             ifelse(d_plot$time=="t+2", 3,
                                    ifelse(d_plot$time=="t+3", 4,
                                           ifelse(d_plot$time=="t+4", 5,
                                                  ifelse(d_plot$time=="t+5", 6,
                                                         ifelse(d_plot$time=="t+6", 7,
                                                                ifelse(d_plot$time=="t+7", 8,
                                                                       ifelse(d_plot$time=="t+8", 9,
                                                                              ifelse(d_plot$time=="t+9", 10, NA))))))))))

unique(d_plot$lead)
d_plot$lead[ d_plot$lead=='1' ] <- " 2 years"
d_plot$lead[ d_plot$lead=='3' ] <- " 4 years"
d_plot$lead[ d_plot$lead=='5' ] <- " 6 years"
d_plot$lead[ d_plot$lead=='7' ] <- " 8 years"
d_plot$lead[ d_plot$lead=='9' ] <- "10 years"

unique(d_plot$year)
d_plot <- subset(d_plot, year != 1 & year != 3 & year != 5 & year != 7 & year != 9)

unique(d_plot$lag)
d_plot$lag[ d_plot$lag==1 ] <- "1-year lag window"
d_plot$lag[ d_plot$lag==3 ] <- "3-year lag window"
d_plot$lag[ d_plot$lag==5 ] <- "5-year lag window"

d_plot <- subset(d_plot, lag=="5-year lag window")

#keep only Atlantic Forest

d_plot$biome = factor(d_plot$biome, levels=c('Atlantic Forest'))

#remove line
#d_plot <- subset(d_plot, biome != "Caatinga" & biome != "Pampa")

g <- ggplot(d_plot, aes(x=year, y=estimate, colour=as.factor(lead), shape=as.factor(lead))) + 
  theme_bw() + geom_hline(yintercept=0, linetype='dashed') +
  geom_point(position=position_dodge(1)) +
  geom_errorbar(aes(ymin=X2.5., ymax=X97.5.), alpha=0.3, size=2, width=0, position=position_dodge(1)) + #CI figs
  geom_errorbar(aes(ymin=c(estimate-std.error), ymax=c(estimate+std.error)), width=0, position=position_dodge(1)) +
  facet_grid(biome ~ method, scales = "free") +
  labs(colour = "Evaluation period", shape = "Evaluation period") +
  theme(legend.position = "bottom") +
  theme( panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank()) +
  scale_x_continuous(breaks = c(2,4,6,8,10)) +
  labs(x = "Years after designation", y = "Indigenous land designation impact on deforestation (%)")
g

ggsave(g, units="cm", width=21, height=28, file="figures/AF_NEG_ATT_Indigenous_lands.svg", dpi=600)





# AF Without Negativas ----------------------------------------------------

#database

set.seed(1)

Sys.setlocale(category = "LC_ALL", locale = "Portuguese") #for encoding

rm(list=ls()) 

#change file location
setwd("C:/R/replication/replication/")

require(foreign)
s <- read.dbf("data/GEOFT_TERRA_INDIGENA.dbf")
table(s$fase_ti)

d <- read.csv("data/Indigenous_land_Y1985_2021.csv")[,-c(1,70,71)]
d[0,]
Encoding(d$NAME) <-"UTF-8"
Encoding(d$ORIG_NAME) <-"UTF-8"
colnames(d)[1:37] <- 1985:2021
d <- d[,c(66,1:37,42,52,61,62)]

require(reshape2)
d_long <- melt(d, id=c("WDPAID","GIS_AREA","NAME","STATUS","STATUS_YR"))
colnames(d_long) <- c("ID","GEE_km2","name","status","status_yr","year","forest")
d <- d_long
rm(d_long)
d$GEE_km2 <- round(d$GEE_km2, digits = 2)
d$forest <- round(d$forest, digits = 4) #in km2
d[0,]

buf <- read.csv("data/Indigenous_land_buffer_Y1985_2021.csv")[,-c(1,39,75,77)]
buf[0,]
Encoding(buf$NAME) <-"UTF-8"
Encoding(buf$ORIG_NAME) <-"UTF-8"
colnames(buf)[1:37] <- 1985:2021
buf <- buf[,c(70,1:37,42,52,63,64)]
buf_long <- melt(buf, id=c("WDPAID","GIS_AREA","NAME","STATUS","STATUS_YR"))
colnames(buf_long) <- c("ID","GEE_km2","name","status","status_yr","year","forest")
buf <- buf_long
rm(buf_long)
buf$GEE_km2 <- round(buf$GEE_km2, digits = 2)
buf$forest <- round(buf$forest, digits = 4)
buf[0,]

d$ID2 <- paste(d$ID, d$year)
buf$ID2 <- paste(buf$ID, buf$year)
d$buf_forest <- buf[match(with(d, ID2), with(buf, ID2)),]$forest 
rm(buf)
table(d$buf_forest > d$forest)
d$buf_forest2 <- d$buf_forest - d$forest 
d$buf_forest <- d$buf_forest2
d <- d[,-ncol(d)]

require(dplyr)
require(Hmisc)
list_ID <- unique(d$ID)
d_new <- d[0,]
for (i in 1:length(list_ID)) {
  temp <- subset(d, ID==list_ID[i])
  temp$forest_lag <- Lag(temp$forest, 1)
  temp$buf_forest_lag <- Lag(temp$buf_forest, 1)
  d_new <- rbind(d_new, temp)
}
d <- d_new
rm(d_new, temp)
d$year <- as.integer(as.character(d$year))
d$def <- round((d$forest_lag - d$forest) * 100, 2) #in hectares
d$buf_def <- round((d$buf_forest_lag - d$buf_forest) * 100, 2) #in hectares

d <- subset(d, year > 1985)
d$def[ d$def < 0] <- 0
d$buf_def[ d$buf_def < 0] <- 0


list_ID <- unique(d$ID)
d_new <- d[0,]
for (i in 1:length(list_ID)) {
  temp <- subset(d, ID==list_ID[i])
  temp$def_lag <- Lag(temp$def, 1)
  temp$buf_lag <- Lag(temp$buf_def, 1)
  d_new <- rbind(d_new, temp)
}
length(unique(d_new$ID))
d <- d_new
rm(d_new, temp)
d <- subset(d, year > 1985)
d$def_lag[ is.na(d$def_lag) ] <- 0
d$buf_lag[ is.na(d$buf_lag) ] <- 0

require(foreign)
cov <- read.dbf("data/IL_covariates_final.dbf")
d$capitals <- round(cov[match(with(d, ID), with(cov, WDPAID)),]$capitals, digits = 2) 
d$towns <- round(cov[match(with(d, ID), with(cov, WDPAID)),]$towns, digits = 2)  
d$ports <- round(cov[match(with(d, ID), with(cov, WDPAID)),]$ports, digits = 2) 
d$soil <- round(cov[match(with(d, ID), with(cov, WDPAID)),]$soil, digits = 2) 
d$dem <- round(cov[match(with(d, ID), with(cov, WDPAID)),]$dem, digits = 2) 
d$slope <- round(cov[match(with(d, ID), with(cov, WDPAID)),]$slope, digits = 2) 
d$fric <- round(cov[match(with(d, ID), with(cov, WDPAID)),]$fric, digits = 4) 
d$biome <- cov[match(with(d, ID), with(cov, WDPAID)),]$biome
d$region <- cov[match(with(d, ID), with(cov, WDPAID)),]$region
rm(cov)

cov <- read.dbf("data/IL_all.dbf")
d$IL_ha <- cov[match(with(d, ID), with(cov, WDPAID)),]$poly_ha 
cov <- read.dbf("data/IL_all_buffer.dbf")
d$buffer_ha <- cov[match(with(d, ID), with(cov, WDPA_PID)),]$poly_ha 

cov <- read.csv("data/annual_precipit_IL_1981_2022.csv")
cov$ID2 <- paste(paste(cov$WDPA_PID, cov$year))
d$mm_annual <- cov[match(with(d, ID2), with(cov, ID2)),]$mm_annual
cov <- read.csv("data/precipt_Aug_Sep_IL_1981_2022.csv")
cov$ID2 <- paste(paste(cov$WDPA_PID, cov$year))
d$mm_Aug_Sep <- cov[match(with(d, ID2), with(cov, ID2)),]$mm_Aug_Sep
cov <- read.csv("data/max_temp_Aug_Sep_IL_1981_2022.csv")
cov$ID2 <- paste(paste(cov$WDPA_PID, cov$year))
d$temp_Aug_Sep <- cov[match(with(d, ID2), with(cov, ID2)),]$temp_Aug_Sep
rm(cov)

str(d)
d$biome <- as.character(d$biome)
d$biome[ d$biome=="AMAZÔNIA"] <- "Amazon"
d$biome[ d$biome=="CERRADO"] <- "Cerrado"
d$biome[ d$biome=="CAATINGA"] <- "Caatinga"
d$biome[ d$biome=="MATA ATLÂNTICA"] <- "Atlantic Forest"
d$biome[ d$biome=="PAMPA"] <- "Pampa"
d$biome[ d$biome=="PANTANAL"] <- "Pantanal"

d$region <- as.character(d$region)
d$region[ d$region=="Norte"] <- "N"
d$region[ d$region=="Nordeste"] <- "NE"
d$region[ d$region=="Centro-Oeste"] <- "CW"
d$region[ d$region=="Sul"] <- "S"
d$region[ d$region=="Sudeste"] <- "SE"

d$def_perc <- round(d$def/d$IL_ha * 100, digits = 3)
d$def_lag_perc <- round(d$def_lag/d$IL_ha * 100, digits = 3)
d$buf_def_perc <- round(d$buf_def/d$buffer_ha * 100, digits = 3)
d$buf_def_lag_perc <- round(d$buf_lag/d$buffer_ha * 100, digits = 3)

d$Amazon <- ifelse(d$biome=="Amazon", 1, 0)
d$Cerrado <- ifelse(d$biome=="Cerrado", 1, 0)
d$Caatinga <- ifelse(d$biome=="Caatinga", 1, 0)
d$Atlantic <- ifelse(d$biome=="Atlantic Forest", 1, 0)
d$Pampa <- ifelse(d$biome=="Pampa", 1, 0)
d$Pantanal <- ifelse(d$biome=="Pantanal", 1, 0)


require(PanelMatch)
require(ggplot2)

length(unique(d$ID))
length(unique(d$year))
table(d$status)
table(d$status_yr)

d <- subset(d, status_yr != 0 & d$IL_ha >= 100)
length(unique(d$ID))
table(d$biome)
13788/24588
5832/24588
1116/24588
3348/24588
396/24588
108/24588

d$treat <- ifelse(d$status=='Designated' & d$status_yr < d$year, 1, 0)
table(d$treat)

#change file name
write.csv(d, "data/AF_IL_ready_data.csv", row.names = F)
d$ID <- as.integer(as.character(d$ID))

#analysis
#filter atlantic biome
#keep only the psw method
#remove biome columns "'Amazon', 'Cerrado', 'Caatinga', 'Atlantic', 'Pampa'"
#from 21 columns to 16 

d <- d %>% filter(biome == "Atlantic Forest")

res <- as.data.frame(matrix(NA,0,7))
colnames(res) <- c('estimate','std.error','2.5%','97.5%','lag','match_size','lead')
res_cov <- as.data.frame(matrix(NA,0,16))
colnames(res_cov) <- c('def_lag', 'buf_def', 'soil','capitals','towns','dem','slope','fric','IL_ha',
                       'mm_annual', 'mm_Aug_Sep', 'temp_Aug_Sep',
                       'lag','match_size','lead','time')
control_list <- as.data.frame(matrix(NA,0,1))
colnames(control_list)[1] <- "weigths"
control_list$ID_match <- as.integer(as.vector(rownames(control_list)))

lags_n <- 5 
size_match_n <- 1
lead_upper_range <- c(1,3,5,7,9)

for (l in 1:length(lags_n)) {
  for (m in 1:length(size_match_n)) {
    for (u in 1:length(lead_upper_range)) {
      
      PM.results <- PanelMatch(lag = lags_n[l],
                               time.id = "year", unit.id = "ID", 
                               treatment = "treat", 
                               refinement.method = "ps.weight",
                               data = d, match.missing = TRUE, 
                               
                               covs.formula = ~def_lag_perc + buf_def_perc + soil + capitals + towns + dem + slope + fric + IL_ha +
                                 mm_annual + mm_Aug_Sep + temp_Aug_Sep,
                               
                               size.match = size_match_n[m],
                               qoi = "att" ,
                               outcome.var = "def_perc",
                               lead = 0:lead_upper_range[u],
                               forbid.treatment.reversal = TRUE, 
                               use.diagonal.variance.matrix = TRUE)
      
      used <- lapply(PM.results$att, attributes) 
      
      for (i in 1:length(used)) {
        temp <- as.data.frame(used[[i]]$weights)
        colnames(temp)[] <- "weigths"
        temp <- subset(temp, weigths > 0)
        temp$ID_match <- as.integer(as.vector(rownames(temp)))
        control_list <- rbind(control_list, temp)
      }
      
      cov <- get_covariate_balance(PM.results$att,
                                   data = d,
                                   covariates = c('def_lag', 'buf_def_perc', 'soil','capitals','towns','dem','slope','fric','IL_ha',    
                                                  'mm_annual', 'mm_Aug_Sep', 'temp_Aug_Sep'),
                                   plot = F, legend = F)
      cov <- as.data.frame(cov)
      cov$lag <- lags_n[l]
      cov$match_size <- size_match_n[m]
      cov$lead <-lead_upper_range[u]
      cov$time <- row.names(cov)
      res_cov <- rbind(res_cov, cov)
      
      
      PE.results <- PanelEstimate(sets = PM.results, data = d,
                                  number.iterations = 1000,
                                  df.adjustment = T,
                                  confidence.level = 0.95)
      
      temp <- summary(PE.results)
      temp <- as.data.frame(temp$summary)
      temp$lag <- lags_n[l]
      temp$match_size <- size_match_n[m]
      temp$lead <-lead_upper_range[u]
      temp$time <- row.names(temp)
      
      res <- rbind(res, temp)
    }
  }
}

#change file name
write.csv(res_cov, "data/AF_cov_balance_PSW_indigenous_lands.csv", row.names = F)
write.csv(res, "data/AF_ATT_PSW_indigenous_lands.csv", row.names = F)

#covariate balance for the indigenous lands in Atlantic Forest 

require(ggplot2)
require(reshape)
require(plyr)

#keep only the psw method
#from all biomes to Atlantic Forest

d <- read.csv('data/AF_cov_balance_PSW_indigenous_lands.csv')
d <- melt(d, id=c('lag','lead','match_size','time'))
PSW <- ddply(.data=d, .(lag,match_size,lead,variable), .fun=summarise, mean = mean(value))
PSW$Method <- 'A'

#keep only the psw

all <- rbind(PSW)
all$biome <- "Atlantic Forest"

d_plot <- all

d_plot$lead <- d_plot$lead + 1

unique(d_plot$Method)

unique(d_plot$variable)
d_plot$variable <- as.character(d_plot$variable)

#keep only atlantic forest variable

d_plot$variable[ d_plot$variable == 'IL_ha' ] <- 'Size'
d_plot$variable[ d_plot$variable == 'soil' ] <- 'Soil quality'
d_plot$variable[ d_plot$variable == 'capitals' ] <- 'Dist. to capitals'
d_plot$variable[ d_plot$variable == 'towns' ] <- 'Dist. to towns'
d_plot$variable[ d_plot$variable == 'capitals' ] <- 'Dist. capitals'
d_plot$variable[ d_plot$variable == 'dem' ] <- 'Elevation'
d_plot$variable[ d_plot$variable == 'slope' ] <- 'Slope'
d_plot$variable[ d_plot$variable == 'fric' ] <- 'Dist. to cities'
d_plot$variable[ d_plot$variable == 'buf_def_perc' ] <- ' 10-km buffer def.'
d_plot$variable[ d_plot$variable == 'def_lag' ] <- 'Lagged def.'
d_plot$variable[ d_plot$variable == 'mm_annual' ] <- 'Annual precip.'
d_plot$variable[ d_plot$variable == 'mm_Aug_Sep' ] <- 'Aug-Sep. precip.'
d_plot$variable[ d_plot$variable == 'temp_Aug_Sep' ] <- 'Aug-Sep. temp.'
d_plot$variable[ d_plot$variable == 'Atlantic' ] <- '  Atlantic forest biome'

unique(d_plot$lag)
unique(d_plot$match_size)
unique(d_plot$lead)

#from all biomes to Atlantic Forest

g1 <- ggplot(data=subset(d_plot, biome=="Atlantic Forest"),
             aes(x=mean, y=variable, group=as.factor(lead), colour=as.factor(lead), shape=as.factor(lead)))  +
  guides(fill=guide_legend(reverse=TRUE)) +  
  theme_bw() +
  theme(legend.position = 'bottom',
        axis.title.y=element_blank(),
        panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(),  panel.grid.minor.x = element_blank()) +
  facet_wrap(. ~ Method, scales = "free", ncol = 2) +
  geom_vline(xintercept=0, linetype='dashed') + geom_point(size=1.7, stroke = 1.1) +
  scale_shape_manual(values=c(0,1,2,3,4)) +
  labs(colour = "Evaluation period (years)") + labs(shape = "Evaluation period (years)") +
  xlab("Standardized mean difference") + scale_x_continuous(breaks = c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3), limits = c(-0.35,0.35)) +
  annotate("rect", xmin = -0.1, xmax = 0.1,
           ymin = 0, ymax = c(as.numeric(length(unique(d_plot$variable)))+1),alpha = .1, fill = "dodgerblue")

g1
ggsave(g1, units="cm", width=20, height=20, dpi=600, file="figures/AF_IL_cov_balance_ALL.svg")
#change file name


#plot

require(plyr)
require(ggplot2)

rm(list=ls()) 

#from all biomes to Atlantic Forest
#keep only the p3

p3 <- read.csv("data/AF_ATT_PSW_indigenous_lands.csv")
p3$method <- "A"
p3$biome <- "Atlantic Forest"

d_plot <- rbind(p3)

d_plot$estimate_minus_SE <- d_plot$estimate + d_plot$std.error
sum(d_plot$estimate_minus_SE < 0)/nrow(d_plot)
ddply(.data=subset(d_plot, estimate_minus_SE < 0), .(time), .fun=summarise, def = mean(estimate))

d_plot$year <- ifelse(d_plot$time=="t+0", 1,
                      ifelse(d_plot$time=="t+1", 2,
                             ifelse(d_plot$time=="t+2", 3,
                                    ifelse(d_plot$time=="t+3", 4,
                                           ifelse(d_plot$time=="t+4", 5,
                                                  ifelse(d_plot$time=="t+5", 6,
                                                         ifelse(d_plot$time=="t+6", 7,
                                                                ifelse(d_plot$time=="t+7", 8,
                                                                       ifelse(d_plot$time=="t+8", 9,
                                                                              ifelse(d_plot$time=="t+9", 10, NA))))))))))

unique(d_plot$lead)
d_plot$lead[ d_plot$lead=='1' ] <- " 2 years"
d_plot$lead[ d_plot$lead=='3' ] <- " 4 years"
d_plot$lead[ d_plot$lead=='5' ] <- " 6 years"
d_plot$lead[ d_plot$lead=='7' ] <- " 8 years"
d_plot$lead[ d_plot$lead=='9' ] <- "10 years"

unique(d_plot$year)
d_plot <- subset(d_plot, year != 1 & year != 3 & year != 5 & year != 7 & year != 9)

unique(d_plot$lag)
d_plot$lag[ d_plot$lag==1 ] <- "1-year lag window"
d_plot$lag[ d_plot$lag==3 ] <- "3-year lag window"
d_plot$lag[ d_plot$lag==5 ] <- "5-year lag window"

d_plot <- subset(d_plot, lag=="5-year lag window")

#keep only Atlantic Forest

d_plot$biome = factor(d_plot$biome, levels=c('Atlantic Forest'))

#remove line
#d_plot <- subset(d_plot, biome != "Caatinga" & biome != "Pampa")

g <- ggplot(d_plot, aes(x=year, y=estimate, colour=as.factor(lead), shape=as.factor(lead))) + 
  theme_bw() + geom_hline(yintercept=0, linetype='dashed') +
  geom_point(position=position_dodge(1)) +
  geom_errorbar(aes(ymin=X2.5., ymax=X97.5.), alpha=0.3, size=2, width=0, position=position_dodge(1)) + #CI figs
  geom_errorbar(aes(ymin=c(estimate-std.error), ymax=c(estimate+std.error)), width=0, position=position_dodge(1)) +
  facet_grid(biome ~ method, scales = "free") +
  labs(colour = "Evaluation period", shape = "Evaluation period") +
  theme(legend.position = "bottom") +
  theme( panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank()) +
  scale_x_continuous(breaks = c(2,4,6,8,10)) +
  labs(x = "Years after designation", y = "Indigenous land designation impact on deforestation (%)")
g

ggsave(g, units="cm", width=21, height=28, file="figures/AF_ATT_Indigenous_lands.svg", dpi=600)




# ############################################################ ------------







# Synthetic Control -------------------------------------------------------


# Load necessary packages
library(Synth)
library(tidyverse)
library(dplyr)

# Path to the CSV file

# Load the CSV file
data <- read.csv(file_path)

# Check if the file was loaded correctly
if (is.null(data) || nrow(data) == 0) {
  stop("The CSV file could not be loaded or is empty.")
}

# Find IDs that received treatment (treat=1) at any point
treated_ids <- data %>%
  filter(treat == 1) %>%
  distinct(ID)

# Create the 'control' column for all IDs, starting with value 0
data$control <- 0

# Assign value 1 in the 'control' column for IDs that received treatment
data$control[data$ID %in% treated_ids$ID] <- 1

# Invert the values in the 'control' column
data$control <- ifelse(data$control == 1, 0, 1)

# Variables of interest
variables_of_interest <- c("buf_def", "capitals", "towns", "soil", "def")

# Calculate the annual average for each variable where control == 1
annual_averages_control <- data %>%
  filter(control == 1) %>%
  group_by(year) %>%
  summarise(across(all_of(variables_of_interest), mean, na.rm = TRUE))

# Assign a unique numeric ID to the new treatment unit
new_treatment_id <- 555544174

# Add the ID and control columns to the table of annual averages
annual_averages_control <- annual_averages_control %>%
  mutate(ID = new_treatment_id, control = 1)

# Exclude observations where control == 1 from the original dataset
filtered_data <- data %>%
  filter(control == 0)

# Add the new treatment unit to the filtered data frame
data_until_2021 <- bind_rows(filtered_data, annual_averages_control)

# Replace missing data with the mean for the variables of interest
data_until_2021 <- data_until_2021 %>%
  group_by(ID) %>%
  mutate(across(all_of(variables_of_interest), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .))) %>%
  ungroup()

# Ensure that the ID column is numeric
data_until_2021$ID <- as.numeric(data_until_2021$ID)

# Check if all IDs are valid
data_until_2021 <- data_until_2021 %>%
  filter(!is.na(ID))

# Check if the dependent variable ("def") is present and has no missing values for the treatment ID and the periods from 1986 to 2010
missing_def_check <- data_until_2021 %>%
  filter(ID == new_treatment_id & year %in% seq(1986, 2010)) %>%
  summarise(missing_def = sum(is.na(def)))

if (missing_def_check$missing_def > 0) {
  stop("The dependent variable 'def' has missing values for the treatment ID in the periods from 1986 to 2010.")
}

# Check if the ID column is numeric
if (!is.numeric(data_until_2021$ID)) {
  stop("The ID column is not in numeric format.")
}

# Check the structure of the data
str(data_until_2021)

# Prepare the data for the synthetic control method, including lagged years for the dependent variable
dataprep.out <- dataprep(
  foo = data_until_2021,
  predictors = variables_of_interest,
  predictors.op = "mean",
  time.predictors.prior = seq(1986, 2010, 1),
  special.predictors = list(
    list("buf_def", seq(1986, 2010, 1), "mean"),
    list("capitals", seq(1986, 2010, 1), "mean"),
    list("towns", seq(1986, 2010, 1), "mean"),
    list("soil", seq(1986, 2010, 1), "mean"),
    list("def", c(1990, 1997, 2006, 2009), "mean")
  ),
  dependent = "def",
  unit.variable = "ID",
  time.variable = "year",
  treatment.identifier = new_treatment_id,
  controls.identifier = unique(data_until_2021$ID[data_until_2021$control == 0]),
  time.optimize.ssr = seq(1990, 2010, 1),
  time.plot = seq(1986, 2020, 1)
)

# Execute the synthetic control method
synth.out <- synth(dataprep.out)

# Obtain synthetic control results
synth.tables <- synth.tab(dataprep.res = dataprep.out, synth.res = synth.out)

# Visualize the results with a line in 2011
path.plot(dataprep.res = dataprep.out, synth.res = synth.out, Ylab = "Deforestation", Xlab = "Year")
abline(v = 2011, col = "black", lwd = 2, lty = 2)
