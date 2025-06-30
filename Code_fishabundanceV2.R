library(party) 
library(caret) 
library(tidyverse)
library(permimp)
library(pdp)
library(gridExtra)
library(sf)
library(ggspatial)
library(geojsonio)
library(readxl)
library(vegan)
library(BiodiversityR)
library(sp)
library(corrplot)




####import  data ####

#CPUE data by survey
cpue.clean<- read_xlsx("/Users/cindypaquette/Google Drive/Mon disque/Postdoc Cindy/Données poissons/cpue_clean.xlsx") 
cpue.clean$cpue.norm<-log(cpue.clean$CPUE)

#Community Presence data
p.a.survey<- read_xlsx("/Users/cindypaquette/Google Drive/Mon disque/Postdoc Cindy/Données poissons/community_pa.xlsx") %>%
  inner_join(cpue.clean[c(1,3)])

#Predictors data Transformed (log or sqrt) and reduced with Pearson correlation threshold >0.9
env<- read_xlsx("/Users/cindypaquette/Google Drive/Mon disque/Postdoc Cindy/Données environnementales/env_norm_pers_16-05-2024.xlsx") 


#### MAPS ####
#Figure 1

#Eastern Canada
Eastern <- st_read("/Users/cindypaquette/Google Drive/Mon disque/Postdoc Cindy/Données environnementales/src/ref/ne_50m_admin_1_states_provinces_lakes/ne_50m_admin_1_states_provinces_lakes.shp", options = "ENCODING=UTF8")%>%
  subset( name %in% c ("Québec", "New Brunswick", "Prince Edward Island", "Nova Scotia", "Newfoundland and Labrador"))

#Distribution area
sp2 <- geojson_read("/Users/cindypaquette/Google Drive/Mon disque/Postdoc Cindy/Données poissons/Aires_repartition_poisson_eau_douce.geojson", what = "sp")%>%
  subset(NOM_FRANCA %in% c ("doré jaune" ,"touladi", "omble de fontaine" ))%>%
  spTransform(CRS("+init=epsg:4269"))

dore<-subset(sp2, NOM_FRANCA %in% "doré jaune")
touladi<-subset(sp2, NOM_FRANCA %in% "touladi")
omble<-subset(sp2, NOM_FRANCA %in% "omble de fontaine")

map1<-
  ggplot(Eastern)+
  geom_sf(fill="lightgray", colour="darkgray", lwd=0.09)+
  geom_sf(data=st_union(Eastern), fill=NA, lwd=0.5)+
  xlab("Longitude") + ylab("Latitude")+
  annotation_scale(location = "br")+
  coord_sf(ylim=c(44,63))+
  guides( fill = guide_legend())+
  theme_bw() + 
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
         plot.title=element_text(size=20, hjust=0.90, vjust=0.5,margin=margin(t=40,b=-30)),
         legend.position = c(.99, .97),
         legend.justification = c("right", "top"),
         legend.box.just = "right",
         legend.margin = margin(3, 4, 6, 6),
         legend.title = element_text(size=16), 
         legend.text = element_text(size=12))

#SAVI
coords.savi<-env%>%
  filter(survey=="PENDJ")%>%
  left_join(cpue.clean)

map1+
  geom_polygon(data=dore,aes(x = long, y = lat, group = group), colour="darkgoldenrod3", fill="lightgray",inherit.aes = FALSE)+
  geom_point(data = coords.savi, mapping = aes(x = long, y = lat,fill=CPUE, size=CPUE, alpha=CPUE),inherit.aes = FALSE,shape=21)+
  scale_size_continuous(name="Walleye\nCPUE", range=c(1,4)) +
  scale_alpha_continuous(name="Walleye\nCPUE", range=c(0.4, 1))+
  scale_fill_gradient(name="Walleye\nCPUE", low = "darkgoldenrod2", high = "darkgoldenrod4")

#SAFO
coords.safo<-env%>%
  filter(survey=="PENOF")%>%
  left_join(cpue.clean)

map1+
  # geom_polygon(data=omble,aes(x = long, y = lat, group = group), colour="deeppink3", fill=NA,inherit.aes = FALSE)+
  geom_point(data = coords.safo, mapping = aes(x = long, y = lat,fill=CPUE, size=CPUE, alpha=CPUE),inherit.aes = FALSE,shape=21)+
  scale_size_continuous(name="Brook trout\nCPUE", range=c(1,4),) +
  scale_alpha_continuous(name="Brook trout\nCPUE", range=c(0.4, 1))+
  scale_fill_gradient(name="Brook trout\nCPUE",low = "deeppink", high = "deeppink3" )

#SANA
coords.sana<-env%>%
  filter(survey=="PENT")%>%
  left_join(cpue.clean)

map1+
  #geom_polygon(data=touladi,aes(x = long, y = lat, group = group), colour="royalblue3", fill="lightgray",inherit.aes = FALSE)+
  geom_point(data = coords.sana, mapping = aes(x = long, y = lat,fill=CPUE,size=CPUE, alpha=CPUE),inherit.aes = FALSE,shape=21)+
  scale_size_continuous(name="Lake trout\nCPUE", range=c(1,4)) +
  scale_alpha_continuous(name="Lake trout\nCPUE", range=c(0.4, 1))+
  scale_fill_gradient(name="Lake trout\nCPUE",low = "royalblue1", high = "royalblue3" )

#### PCoA ####
#Figure 2

setwd("/Users/cindypaquette/Documents/postdoc/analyses/1.Abundance models/figures/pcoa/cmdscale")

pal_color <- c("darkgoldenrod3", "deeppink3","royalblue3")
names(pal_color) <-c("PENDJ","PENOF", "PENT")

letter<-c("a)","b)","c)")
names(letter) <-c("PENDJ","PENOF", "PENT")

species.axis<-data.frame()

for(i in unique(p.a.survey$survey[!is.na(p.a.survey$survey)])){
  
  species.survey<- p.a.survey %>%
    filter(survey == i)%>%
    column_to_rownames(var="LCE")%>%
    select(c(1:91))
  
  species.survey<-species.survey[,colSums(species.survey)>0]
  species.survey<-species.survey[rowSums(species.survey) > 0,] 
  
  #jaccard dissimilarity 
  species.jac <- vegdist(species.survey,  method = "jac", binary=T)
  
  #pcoa with cmdscale
  species.pcoa <- cmdscale(species.jac , eig = TRUE, k = (nrow(species.survey) - 1)) 
  
  #save relative eigenvalues
  releig<-as.data.frame(species.pcoa$eig / sum(species.pcoa$eig))
  
  #add species scores
  species.scores <- add.spec.scores(species.pcoa, species.survey, 
                                    method='pcoa.scores', Rscale=TRUE, scaling=1, multi=1)
  #extract scores
  species.df6<- as.data.frame(scores(species.scores, display="species", choices=c(1,2),scaling=1))
  species.df7<- as.data.frame(scores(species.scores, display="sites", choices=c(1,2),scaling=1)) 
  
  species.df7<-species.df7 %>%
    rownames_to_column("LCE")%>%
    inner_join(cpue.clean[cpue.clean$survey==i,c(1,3)]) 
  
  #add cpue data
  species.df7<-cpue.clean %>%
    filter(survey == i) %>%
    left_join(species.df7)%>%
    drop_na()
  
  #plot
  ggplot(species.df6, aes(x=Dim1, y=Dim2))+ theme_bw() +
    geom_point(data=species.df7, aes(x=Dim1, y=Dim2, color=survey, size=CPUE, alpha=CPUE))+
    scale_size_continuous(name="CPUE",range=c(1,4)) +
    scale_alpha_continuous(name="CPUE", range=c(0.4, 1))+
    scale_color_manual(name="CPUE",values = pal_color)+
    geom_segment(data=species.df6, aes(x=0, xend=Dim1/35, y=0, yend=Dim2/35),color="darkgray", arrow=arrow(length=unit(0.02,"npc")))+
    geom_text(data=species.df6,aes(x=Dim1/33, y=Dim2/33, label=rownames(species.df6)), color="black", size=6)+
    ggtitle(letter[i])+
    xlab(paste0("Axis 1 (", round(releig[1,]*100,2),"%)"))+
    ylab(paste0("Axis 2 (", round(releig[2,]*100,2),"%)"))+
    theme(plot.title=element_text(hjust=0.05, vjust=0.5,size=20,margin=margin(t=40,b=-30)),
          axis.text=element_text(size=12), axis.title=element_text(size=18))
  
  #ggsave(paste0(i ,"_PCOA",".pdf")) 
  
}

#### Random forest ####

seed <- 2021
mtry.best=ceiling(sqrt(ncol(env)-2)) 

##### SAVI #####

env.savi<-env%>%
  filter(survey=="PENDJ")

dat.savi <- cpue.clean %>%
  filter(species == 'SAVI') %>%
  select(LCE, cpue.norm) %>%
  inner_join(env.savi) %>%
  ungroup() %>%
  select(-survey) %>%
  column_to_rownames("LCE")%>%
  as.data.frame

## Run model on full dataset
rf.savi <- cforest(cpue.norm ~ ., data=dat.savi, controls = cforest_unbiased(ntree = 1000, mtry=mtry.best))

#diagnostic plots
fitted<-predict(rf.savi) 
obs<-dat.savi$cpue.norm
resid<-obs-fitted

par(mfrow = c(1, 2), mar = c(4, 4, 3, 2))

# (a) Residuals vs Fitted
plot(resid ~ fitted,
     xlab = "Fitted values",
     ylab = "Residuals",
     main = "(a) Walleye \n fitted vs residuals",
     cex.main = 0.95,
     pch = 16)
abline(h = 0, lty = 2)

# (b) Observed vs Fitted
plot(obs ~ fitted,
     xlab = "Fitted values",
     ylab = "Observed values",
     main = "(b) Walleye \n fitted vs observed",
     cex.main = 0.95,
     pch = 16)
abline(a = 0, b = 1, lty = 2)

#Full model R2
R2_savi<-round(1 - sum((obs-fitted)^2)/sum((obs-mean(obs))^2),2) 
R2_savi <- bquote(italic(r)^2==.(R2_savi))


## Variable importance
vars.savi <- varimp(rf.savi, conditional = F)
vars.savi <- sort((vars.savi/sum(vars.savi))*100)

#####SAFO #####
env.safo<-env%>%
  filter(survey=="PENOF")

dat.safo <- cpue.clean %>%
  filter(species == 'SAFO') %>%
  select(LCE, cpue.norm) %>%
  inner_join(env.safo) %>%
  ungroup() %>%
  select(-survey) %>%
  column_to_rownames("LCE")%>%
  #   select(-fishing_pressure) #for analysis without the fishing pressure
  as.data.frame

## Run model on full dataset
set.seed(seed)
rf.safo <- cforest(cpue.norm ~ ., data=dat.safo, controls = cforest_unbiased(ntree = 1000, mtry=mtry.best))

#diagnostic plots
fitted<-predict(rf.safo) 
obs<-dat.safo$cpue.norm
resid<-obs-fitted

par(mfrow = c(1, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 1, 0))
# (c) Residuals vs Fitted
plot(resid ~ fitted,
     xlab = "Fitted values",
     ylab = "Residuals",
     main = "(c) Brook trout \n fits vs residuals",
     cex.main = 0.95,
     pch = 16)
abline(h = 0, lty = 2)

# (d) Observed vs Fitted
plot(obs ~ fitted,
     xlab = "Fitted values",
     ylab = "Observed values",
     main = "(d) Brook trout \n fits vs observed",
     cex.main = 0.95,
     pch = 16)
abline(a = 0, b = 1, lty = 2)

#Full model R2
R2_safo<-round( 1 - sum((obs-fitted)^2)/sum((obs-mean(obs))^2),2) 
R2_safo <- bquote(italic(r)^2==.(R2_safo))


## Variable importance
vars.safo <- varimp(rf.safo, conditional = F)
vars.safo <- sort((vars.safo/sum(vars.safo))*100)

##### SANA #####
env.sana<-env%>%
  filter(survey=="PENT")

dat.sana <- cpue.clean %>%
  filter(species == 'SANA') %>%
  select(LCE, cpue.norm) %>%
  inner_join(env.sana) %>%
  ungroup() %>%
  select(-survey) %>%
  column_to_rownames("LCE")%>%
  as.data.frame

## Run model on full dataset
set.seed(seed)
rf.sana <- cforest(cpue.norm ~ ., data=dat.sana, controls = cforest_unbiased(ntree = 1000, mtry=mtry.best))

#diagnostic plots
fitted<-predict(rf.sana) 
obs<-dat.sana$cpue.norm
resid<-obs-fitted

par(mfrow = c(1, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 1, 0))
# (e) Residuals vs Fitted
plot(resid ~ fitted,
     xlab = "Fitted values",
     ylab = "Residuals",
     main = "(e) Lake trout \n fits vs residuals",
     cex.main = 0.95,
     pch = 16)
abline(h = 0, lty = 2)

# (f) Observed vs Fitted
plot(obs ~ fitted,
     xlab = "Fitted values",
     ylab = "Observed values",
     main = "(f) Lake trout \n fits vs observed",
     cex.main = 0.95,
     pch = 16)
abline(a = 0, b = 1, lty = 2)

## Variable importance
vars.sana <- varimp(rf.sana, conditional = F)
vars.sana <- sort((vars.sana/sum(vars.sana))*100)

#Full model R2
R2_sana<-round(1 - sum((obs-fitted)^2)/sum((obs-mean(obs))^2),2) 
R2_sana <- bquote(italic(r)^2==.(R2_sana))


#### RF barplot ####
#Figure 3

#plot only vars >0
par(mfrow=c(1,3))

#SAVI
par(mar=c(4,9,1,1))
barplot(vars.savi[vars.savi>0],col="darkgoldenrod3",border=0,horiz=T,las=1,cex.names=1.2) 
legend(0,1,bty='n',legend=R2_savi, cex=1.2)
title(xlab='relative influence (%)', cex.lab=1.2)

#SAFO
par(mar=c(4,8.5,1,1))
barplot(vars.safo[vars.safo>0],col="deeppink3",border=0,horiz=T,las=1,cex.names=1.2) 
legend(0,1,bty='n',legend=R2_safo, cex=1.2)
title(xlab='relative influence (%)', cex.lab=1.2)

#SANA
par(mar=c(4,8.5,1,1))
barplot(vars.sana[vars.sana>0],col="royalblue3",border=0,horiz=T,las=1,cex.names=1.2) 
legend(0,1,bty='n',legend=R2_sana, cex=1.2)
title(xlab='relative influence (%)', cex.lab=1.2)



#### Supplements ####

##### partial dependency plots #####
#Figures S3-5

#SAVI
pdp_savi<-list()

for (i in names(vars.savi[order(vars.savi, decreasing = T)])[1:10]) {
  message(i)
  pdp_savi[[i]] <- local({
    i <- i
    autoplot(partial(rf.savi, pred.var = i), rug=T, train=dat.savi, size=1.3, col="darkgoldenrod3")
  })}

grid.arrange(grobs = pdp_savi[1:10], ncol = 4)


#SAFO
pdp_safo<-list()

for (i in names(vars.safo[order(vars.safo, decreasing = T)])[1:10]) {
  message(i)
  pdp_safo[[i]] <- local({
    i <- i
    autoplot(partial(rf.safo, pred.var = i), rug=T, train=dat.safo, size=1.3, col="deeppink3")
  })}

grid.arrange(grobs = pdp_safo[1:10], ncol = 4)


#SANA
pdp_sana<-list()

for (i in names(vars.sana[order(vars.sana, decreasing = T)])[1:10]) {
  message(i)
  pdp_sana[[i]] <- local({
    i <- i
    autoplot(partial(rf.sana, pred.var = i), rug=T, train=dat.sana, size=1.3, col="royalblue3")
  })}

grid.arrange(grobs = pdp_sana[1:10], ncol = 4)


##### CPUE time series #####
#Figure S1
cpue<- read_xlsx('/Users/cindypaquette/Google Drive/Mon disque/Postdoc Cindy/Données poissons/PEN_FisHab_Janv 2023.xlsx', sheet=1)%>%
  dplyr::select( 2:4,10, 11)%>%
  filter(!grepl('Réservoir|réservoir', nomlac)) %>%
  subset(!nolac %in% c("01269","04080000","04082000","04300000", "00884","21421")) %>%
  filter(!(cpue == 0.1250000 & sp == "SAVI"))

colnames(cpue) <- c('LCE',"nomlac", 'year','CPUE','species') 

cpue$cpue_start<-str_split(cpue$year,"-",simplify=T)[,1]

mean_cpue2 <- cpue %>%
  mutate(cpue_start = as.integer(cpue_start))%>%
  group_by(species,cpue_start)%>%
  summarize(
    Mean_CPUE = mean(CPUE, na.rm = TRUE),
    SE_CPUE = sd(CPUE, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )


pal_color <- c("darkgoldenrod3", "deeppink3","royalblue3")
names(pal_color) <-c("SAVI","SAFO", "SANA")

ggplot(mean_cpue2, aes(x = cpue_start, y = Mean_CPUE, group=species, color=species)) +
  geom_line(linewidth=2) +
  scale_color_manual(name="CPUE",values = pal_color)+
  geom_point() +
  geom_errorbar(aes(ymin = Mean_CPUE - SE_CPUE, ymax = Mean_CPUE + SE_CPUE), width = 0.2) +
  theme_minimal() +
  labs(
    x = "Year",
    y = "Mean CPUE (±SE)")+
  theme(legend.position="inside", legend.position.inside=c(0.9,0.5),axis.text=element_text(size=12),
        axis.title=element_text(size=14))


##### correlation plot #####
#Figure S2
rcorrDat<- Hmisc::rcorr(as.matrix(env[c(2:36)]), type="pearson") 

M <- rcorrDat$r
par(mfrow=c(1,1))

corrplot(M, 
         method = "square", 
         type = "lower", 
         tl.col = "black", 
         tl.cex=0.5,
         tl.srt = 45, 
         addCoef.col = 'black',number.cex=0.4,number.digits=2)


##### sampling frequency #####
#Table S1
cpue.sampling<-cpue%>%group_by(LCE,species)%>%
  dplyr::summarise(number_samples=n_distinct(year))%>%
  group_by(species, number_samples)%>%
  summarise(count=n())


