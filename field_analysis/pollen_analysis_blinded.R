rm(list=ls())
setwd("field_analysis")
library(tidyverse)
library(furrr)
library(readxl)
library(PROreg)
library(lubridate)
library(glmmTMB)
library(DHARMa)
library(gghalves)
library(egg)
library(MuMIn)
map=purrr::map; select=dplyr::select

#load data. list of focal plants
focal=read_excel('MethodsTable2021_2022-21jan2022.xlsx') %>% 
    rename_all(tolower) %>%
    mutate(genus=sub(" .*","",species),epithet=sub(".* ","",species)) %>% 
    mutate(plant_code=ifelse(genus=='Salix','salix',paste0(tolower(substr(genus,1,3)), tolower(substr(epithet,1,3))))) %>%
    select(-genus,-epithet)
focal[focal$species=="Saxifraga virginiensis",]$plant_code <- 'micvir'
focal[focal$species=="Viburnum opulus var. americanum",]$plant_code <- 'vibtri'

#make table for just the specialists
spec_hosts=focal %>% filter(type == 'Specialist host')
non_hosts=focal %>% filter(type == 'Non-host')
nrow(non_hosts);nrow(spec_hosts)

#load bee species data
genera=read_csv("pollen specialization bee id spreadsheet-19jan2022.csv") %>% #load genus level id of bees
    rename(unique_id=uniqueID)%>%mutate(bee=ifelse(!is.na(species),paste0(genus,"_",species),paste0(genus,"_sp")))

#load pinning log from both years
specs=read_csv( "specialist bee pinning log2020and2021bees-cleaned21jan2022.csv") %>%
    filter(data)%>% 
    left_join(genera %>% select(unique_id,genus,species))


#load prop pollen data
data=read_csv('proportion_pollen_data_19jan2022.csv') %>% 
    left_join(genera %>% select(unique_id,genus)) %>%
    filter(!is.na(genus))


#randomize the pollen data for blind analysis
set.seed(1212)
data_b=data %>% split(.$site_spec) %>%
    map_dfr(function(df) df %>% mutate(host_plant=sample(host_plant))) %>%
    mutate(spec_host=host_plant %in% spec_hosts$plant_code) %>%
    mutate(plant_type=ifelse(spec_host,'host','nonhost'),year_factor=as.factor(year))

# finally, load data of all the plants we collected from,
# organized by site/sampling round 
# (this is so i can include the plants we didn't collect any bees from in 'total pollen' analysis)
sitedate_specplant=read_csv('sitedate_specialisthostplant-21jan2022.csv')
schedule=read_csv("Specialist bee plant schedule 2020_2021-cleaned21jan2022.csv") %>%
    rename(genus_species="plant genus and species") %>%
    distinct(site,genus_species,date) %>% 
    filter(genus_species %in% focal$species)%>% left_join(sitedate_specplant) %>%
    mutate(plant_code=paste0(tolower(substr(genus_species,1,3)),substr(gsub(".* ","",genus_species),1,3)))
schedule[schedule$genus_species=="Saxifraga virginiensis",]$plant_code <- 'micvir'
schedule[schedule$genus_species=="Viburnum opulus var. americanum",]$plant_code <- 'vibtri'
schedule[schedule$genus_species=="Salix sp.",]$plant_code <- 'salix'

##get rid of cucpep sites from schedule (for sites where cucpep was the only specialist host plant)
# for site where we also sampled solidago - just get rid of cucpep from the schedule
head(schedule)
schedule=schedule %>% filter(spec_host!='cucpep' &  plant_code !='cucpep')

#convert roman numeral date to formatted date for calculations
dates_rn=unique(schedule$date)
year=as.numeric(substr(dates_rn,nchar(dates_rn)-3,nchar(dates_rn)))
day_month_rn=substr(dates_rn,1,nchar(dates_rn)-4)
day=ifelse(!is.na(as.numeric(substr(day_month_rn,1,2))),as.numeric(substr(day_month_rn,1,2)),as.numeric(substr(day_month_rn,1,1)))
lgl_ineed=!is.na(as.numeric(substr(day_month_rn,1,2)))
month_rn=ifelse(lgl_ineed,substr(day_month_rn,3,nchar(day_month_rn)),substr(day_month_rn,2,nchar(day_month_rn)))
month_conversion=data.frame(month=c(4,5,6,7,8),month_rn=c("iv","v","vi","vii","viii"))
dates_df=data.frame(month_rn) %>% 
    left_join(month_conversion) %>% 
    mutate(date=dates_rn,year=year,day=day) %>%
    mutate(date_formatted=dmy(paste(day,month,year,sep='-')))

#max distance between dates at a site
schedule %>% left_join(dates_df %>% select(date,date_formatted)) %>%
    split(.$site_spec) %>% map(function(df) max(df$date_formatted)-min(df$date_formatted))


#####model n as an intercept only model with the randomized ata 
# (this is for getting model params for phi, rand bee genus rand effect 
# to use to simulate blinded data)

# mod_bee=BBmm(n~plant_type,random.formula = ~genus,m=200,data=data_b)
# saveRDS(mod_bee,'mod_bee-22dec2021.rds')
mod_bee=readRDS('mod_bee-22dec2021.rds')
summary(mod_bee)


#use these model params (for year and phi) to generate a fake dataset
summary(mod_bee)#intercept=intercept from mod_bee
(beta0=mod_bee$fixed.coef[1])
x1=as.numeric(data_b$plant_type=='host')
beta <- c(beta0,-0.6)
(phi <- mod_bee$phi.coef)
my_genera=as.factor(data_b$genus)
Z=model.matrix(~my_genera-1)

u=rnorm(n_distinct(my_genera),0,.5)
eta <- beta[1]+beta[2]*x1+crossprod(t(Z),u)

#now generate the blinded data
p <- 1/(1+exp(-eta))
y <- rBB(length(my_genera),200,p,phi)
df_faked=data.frame(n=y,plant_type=data_b$plant_type,year_factor=data_b$year_factor,
                    genus=data_b$genus,site=data_b$site,site_spec=data_b$site_spec,plant_code=data_b$host_plant) %>%
    mutate(prop=n/200)

df_faked2=df_faked %>% mutate(n_nonhost=200-n)


# ###########
# ###uncomment me to run analysis without: Echium vulgare, Convolvulus arvensis,
# #Sororia sorbifolia, Prunus susquehanae and Vicia cracca
# #note I think viccra was only sampled at a cucpep so shouldn't be in any of the dfs
# schedule=schedule %>% filter(!plant_code %in% c("prusus","sorsor","conarv","echvul","viccra"))
# data=data %>% filter(!host_plant %in% c("prusus","sorsor","conarv","echvul","viccra"))
# data_b=data_b %>% filter(!host_plant %in% c("prusus","sorsor","conarv","echvul","viccra"))
# df_faked2=df_faked2 %>% filter(!plant_code %in% c("prusus","sorsor","conarv","echvul","viccra"))
# 

#conduct blind analysis on the simulated data
#start with model with plant type as fixed effect and genus and site as random effects
mod_start=glmmTMB(cbind(n, n_nonhost)~plant_type + (1|genus) +(1|plant_code),family=betabinomial,data=df_faked2)
mod_updated1=update(mod_start,.~.+(1|site))
mod_updated2=update(mod_start,.~.+year_factor)
mod_updated3=update(mod_start,.~.+year_factor+(1|site))

(aic_tab=(AIC(mod_start,mod_updated1,mod_updated2,mod_updated3)))

aic_tab$mod=row.names(aic_tab)
aic_tab=as_tibble(aic_tab) %>% select(mod,AIC,df)
aic_tab_ms=aic_tab %>% 
    mutate("predictor variable added"=c('starting model','+ (1|site)','+ year','+ year + (1|site)'),
           "response variable"="proportion visited plant's pollen") %>%
    select('response variable','predictor variable added', AIC,df) %>%
    rename('AIC/AICc'=AIC)

##i tried running these to fix the convergence issues of mod_updated1 -didn't work
## these lines of code are commented out because they take a while to run
# mod_updated1_optim=update(mod_updated1,control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS")))
# mod_updated1_optim2=update(mod_updated1,control=glmmTMBControl(optimizer=optim,optArgs=list(method="CG")))


#pick the final model using AIC vals
mod_final=mod_start

aic_modstart=aic_tab[aic_tab$mod=='mod_start',]$AIC
aic_tab$delta_startmod=aic_tab$AIC-aic_modstart
better_mods=aic_tab %>% filter(delta_startmod<=-2)

#go with more complicated mod if it improves the AIC by more than 2
if(nrow(better_mods) !=0){
    mod_final_name=as.name(better_mods[better_mods$AIC==min(better_mods$AIC),]$mod)
    mod_final=eval(mod_final_name)
    }

#check final model meets assumptions:
simulationOutput_finalmod <- simulateResiduals(fittedModel = mod_final)
plot(simulationOutput_finalmod)#

#see if assumptions are met on randomized data as well
data_b$obs=1:nrow(data_b)
a_mod=glmmTMB(cbind(n, 200-n)~plant_type + (1|genus) +(1|host_plant) ,family=betabinomial,data=data_b)
simout=simulateResiduals(fittedModel = a_mod) 
plot(simout) #nope...

#now get effect sizes:
fixed_effs=as.data.frame(coef(summary(mod_final))$cond)[,1]

# how many pollen grains does a generalist collect from host v nonhost plants? 
# (note that this code may need to change a bit, 
#if year is included in the final model when the data gets unblinded)
log_odds_nonhost=sum(fixed_effs)
(prob_nonhost=exp(log_odds_nonhost)/(1+exp(log_odds_nonhost)))
(numb_nonhosts=prob_nonhost*200) # number of pollen grains from a nonhost plant

#odds of bee collected from non-host plant in 2020
(log_odds_host=fixed_effs[1])
(prob_host=exp(log_odds_host)/(1+exp(log_odds_host)))
(numb_hosts=prob_host*200)

percent_diff=(numb_nonhosts-numb_hosts)/numb_hosts*100

paste0("Generalist bees visiting plants used by specialist bees carried an 
       estimated ", round(numb_hosts)," grains of pollen out of the 200 maximum possible that we counted.
       By contrast generalist bees visiting other plants carried an estimated ",round(numb_nonhosts),
       ", which amounted to a percent increase of ", round(percent_diff,2), "%",")")
       
coef_table=data.frame(coef(summary(mod_final))$cond) %>%
    mutate(coef=row.names(.)) %>% rename(estimate=Estimate, std_error=Std..Error,
                                         zval=z.value,pval=Pr...z..)
plant_type_stats=coef_table %>% filter(coef=="plant_typenonhost")
pval_char=as.character(round(plant_type_stats$pval,2))
if(plant_type_stats$pval<0.01){
    pval_char=as.character(plant_type_stats$pval)
    
    pval_char=paste0(as.character(round(as.numeric(gsub("*e.*","",pval_char)),2)),
    'e',
    gsub(".*e","",pval_char))
    
}
paste0("(plant type estimate +- SE = ", round(plant_type_stats$estimate,2),"+-",round(plant_type_stats$std_error,2),
       ", z-value = ",round(plant_type_stats$zval,2),', p-value = ',pval_char,')')

### next analysis: look at the total pollen removed
# calculate the total pollen removed 
# blind the data by shuffling plants within sites

#code below assigns  zero values to plants that we didn't collect any bees off of
i=20
prop_per_plantsite=1:(n_distinct(data$site_spec)) %>% map(function(i){ 
    
    #get the site round
    the_site_spec=unique(data$site_spec)[i]
    
    #plants_atsite includes all the plants we sampled from at a site
    #even if we didn't collect any bees from them
    plants_atsite=schedule %>% filter(site_spec==the_site_spec) %>% 
        distinct(plant_code)
    
    #filter data to be just that site round
    df=data %>% filter(site_spec==the_site_spec)
    
    total_pollen_df=df %>% 
        group_by(host_plant) %>% 
        summarize(mean_prop=mean(prop),sum_prop=sum(prop),sum_n=sum(n)) 
    
    #if 
    if(sum(!plants_atsite$plant_code %in% total_pollen_df$host_plant)>0){
        add_me_plants=plants_atsite$plant_code[!plants_atsite$plant_code %in% total_pollen_df$host_plant]
        total_pollen_df=total_pollen_df %>%
            bind_rows(data.frame(host_plant=add_me_plants,mean_prop=0,sum_prop=0,sum_n=0))
        }
    total_pollen_blinded=total_pollen_df %>%
        mutate(host_plant=sample(host_plant))%>% #blind the data
        mutate(spec_host=host_plant %in% spec_hosts$plant_code,
               site_spec=the_site_spec,year=as.character(df$year[1]),
               site=df$site[1]) %>%
        mutate(sum_prop2=sum_prop+.001,plant_type=ifelse(spec_host,'host','nonhost'))
    
    })

#bind everything together to make a df for the analysis
total=prop_per_plantsite %>% bind_rows %>% 
    mutate(id=1:n(),year_factor=as.factor(year)) 

#plot some graphs
hist(total$sum_prop2)
hist(total$sum_n)

median(total$sum_prop2)
mean(total$sum_n)

# write_csv(total,"field_data/total_randomized_20dec2021.csv")

#run the models 'vismod' for visitation model
vismod_start=glmmTMB(sum_prop2~plant_type+(1|host_plant) ,family=Gamma(link='log'),data=total)
vismod_updated1=update(vismod_start,.~.+(1|site))
vismod_updated2=update(vismod_start,.~.+year_factor)
vismod_updated3=update(vismod_start,.~.+year_factor+(1|site))
aic_tab_vis=data.frame(AICc(vismod_start,vismod_updated1,vismod_updated2,vismod_updated3)) %>%
    mutate(mod=row.names(.))

#pick the final model using AIC vals
vismod_final=vismod_start

aic_vismodstart=aic_tab_vis[aic_tab_vis$mod=='vismod_start',]$AIC
aic_tab_vis$delta_startmod=aic_tab_vis$AIC-aic_vismodstart
better_vismods=aic_tab_vis %>% filter(delta_startmod<=-2)

#make aic table for the supplement
aic_tab_ms
aic_tab_vis

aic_tab_vis_ms=aic_tab_vis %>% 
    mutate("predictor variable added"=c('starting model','+ (1|site)','+ year','+ year + (1|site)'),
           "response variable"="total pollen removed") %>%
    select('response variable','predictor variable added', AICc,df) %>%
    rename('AIC/AICc'=AICc)
aic_tab_ms_both=aic_tab_ms %>%bind_rows(aic_tab_vis_ms)
# write_csv(aic_tab_ms_both,'figures/aic_table-blinded.csv')
#go with more complicated mod if it improves the AIC by more than 2
if(nrow(better_vismods) !=0){
    vismod_final_name=as.name(better_vismods[better_vismods$AIC==min(better_vismods$AIC),]$mod)
    vismod_final=eval(vismod_final_name)
}

#check model meets assumptions:
simulationOutput <- simulateResiduals(fittedModel = vismod_final)
plot(simulationOutput)
summary(vismod_final)

#report summary stats and effect sizes
paste0("Generalist bees visiting host plants used by specialist bees removed on average,
x% more pollen grains than ones visiting nonhost plants not used by any specialists")


coef_vistable=data.frame(coef(summary(vismod_final))$cond) %>%
    mutate(coef=row.names(.)) %>% rename(estimate=Estimate, std_error=Std..Error,
                                         zval=z.value,pval=Pr...z..)
total_poll_host=exp(coef_vistable$estimate[1])
total_poll_nonhost=exp(sum(coef_vistable$estimate))

if(coef_vistable$estimate[2]>0){
    moreorless='fewer'
    percent_diffvis=(total_poll_nonhost-total_poll_host)/total_poll_nonhost
}
if(coef_vistable$estimate[2]<0){
    moreorless='more'
    percent_diffvis=(total_poll_host-total_poll_nonhost)/total_poll_host
}

#report summary stats and effect sizes
paste0("Generalist bees visiting plants used by specialist bees removed an estimated ",
round(percent_diffvis,2),"% ",moreorless," pollen grains than ones visiting plants not used by any specialists")


plant_type_visstats=coef_vistable %>% filter(coef=="plant_typenonhost")
vis_pval_char=as.character(round(plant_type_visstats$pval,2))
if(plant_type_visstats$pval<0.01){
    vis_pval_char=as.character(plant_type_visstats$pval)
    
    vis_pval_char=paste0(as.character(round(as.numeric(gsub("*e.*","",vis_pval_char)),2)),
                     'e',
                     gsub(".*e","",vis_pval_char))
    
}
paste0("(plant type estimate +- SE = ", round(plant_type_visstats$estimate,2),"+-",round(plant_type_visstats$std_error,2),
       ", z-value = ",round(plant_type_visstats$zval,2),', p-value = ',vis_pval_char,')')



##now plot the results
text_size=20
with(total,boxplot(sum_prop2~plant_type))
with(df_faked,boxplot(prop~plant_type))
my_cols=RColorBrewer::brewer.pal(8,"Accent")
prop_pollen_plot=ggplot(df_faked,aes(x=plant_type,y=prop,color=plant_type))+
    geom_half_boxplot(aes(fill = plant_type, alpha = 0.8),  color= "black", nudge = 0.01, outlier.color = NA) +
    geom_half_violin(aes(fill = plant_type, alpha = 0.8), color= "black",
                     side = "r", nudge = 0.01)+
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text = element_text(face = "italic"),
          strip.background = element_blank(),
          text = element_text(size = text_size),
          legend.position = "none") +
    xlab("Plant type") +
    ylab("Proportion of visited plant's pollen in pollen load")+
    scale_fill_manual(values = c(my_cols[1], my_cols[2])) 
(total_pollen_plot=ggplot(total,aes(x=plant_type,y=sum_prop2,color=plant_type))+
    geom_half_boxplot(aes(fill = plant_type, alpha = 0.8),  color= "black", nudge = 0.01, outlier.color = NA) +
    geom_half_violin(aes(fill = plant_type, alpha = 0.8), color= "black",
                     side = "r", nudge = 0.01)+
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text = element_text(face = "italic"),
          strip.background = element_blank(),
          text = element_text(size = text_size),
          legend.position = "none") +
    xlab("Plant type") +
    ylab("Total pollen removed")+
    scale_fill_manual(values = c(my_cols[1], my_cols[2])) )


    # theme_bw(base_size=12)+theme(
    #     # Hide panel borders and remove grid lines
    #     panel.border = element_blank(),
    #     panel.grid.major = element_blank(),
    #     panel.grid.minor = element_blank(),
    #     # Change axis line
    #     axis.line = element_line(colour = "black"),
    #     text = element_text(size=20)
    # )+

# pdf('blinded_results_pollendata.pdf',width=12)
ggarrange(prop_pollen_plot,  total_pollen_plot,
          labels = c("A", "B"),
          ncol = 2, nrow = 1)
# dev.off()


#old
library(ggpol)
a=ggplot(df_faked,aes(x=plant_type,y=prop,fill=plant_type))+
    geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA, 
                    errorbar.draw = TRUE) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text = element_text(face = "italic"),
          strip.background = element_blank(),
          text = element_text(size = 15),
          legend.position = "none") +
    xlab("Plant type") +
    ylab("Proportion of visited plant's pollen in pollen load")+
    scale_fill_manual(values = c(my_cols[1], my_cols[2])) #+


#plot the results:
ggplot(nshift, aes(x = survey, y = cwm, color = survey)) +
    facet_wrap("habitat ~ .", scales = "free") +
    geom_half_boxplot(aes(fill = survey, alpha = 0.8),  color= "grey50", nudge = 0.05, outlier.color = NA) +
    geom_half_violin(aes(fill = survey, alpha = 0.8), color= "grey50",
                     side = "r", nudge = 0.01) +
    stat_summary(fun = mean, geom = "point", shape = 24, aes(fill = survey, color = survey), size = 4) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text = element_text(face = "italic"),
          strip.background = element_blank(),
          text = element_text(size = 15),
          legend.position = "none") +
    xlab("Survey period") +
    ylab("CWM N-number")+
    scale_fill_manual(values = c("lightcyan", "thistle1")) +
    scale_color_manual(values = c("#29e786", "#E7298A")) -> cwmn_supp

#######OLD:
(m <- glmer(sum_prop2 ~ plant_type +year_factor+ (1 | site) + (1|id) , family = Gamma,data=total))
(m2 <- glmer(sum_n ~ plant_type +year+ (1 | site) + (1|host_plant)+(1|id) , family = poisson,data=total))
(m3 <- glmer(sum_n ~ plant_type +year+ (1 | site) + (1|host_plant), family = negative.binomial(2),data=total))
(m4 <- glmer(sum_n ~ plant_type +year+ (1 | site) + (1|host_plant)+(1|id), family = negative.binomial(2),data=total))

(m <- glmer(sum_n ~ plant_type +year+ (1 | site)  , family = poisson,data=total))
(m <- glmer(sum_n ~ plant_type +year + (1|host_plant) , family = poisson,data=total))



library("blmeco") 
dispersion_glmer(m2) # over 1.4 means there's overdispersion
dispersion_glmer(m3) 
plot(m2)
plot(m3)

par(mfrow=c(1,3))
plot(fitted(m2) ~ total$sum_n, col="darkgrey", 
     xlab="Y (response)", ylab="Fitted Values",main='poisson with orle')
abline(a=0, b=1, col="red")
plot(fitted(m3) ~ total$sum_n, col="darkgrey",
     xlab="Y (response)", ylab="Fitted Values",main='negative binomial (no orle)')
abline(a=0, b=1, col="red")
plot(fitted(m4) ~ total$sum_n, col="darkgrey",
     xlab="Y (response)", ylab="Fitted Values", main='negative binomial with orle')
abline(a=0, b=1, col="red")

#From: https://rpubs.com/INBOstats/OLRE
#"Models with an OLRE should be used carefully. Because OLRE can have a very 
# strong influence on the model. One should always check the standard deviation 
# of the ORLE."
summary(m2) #std deviation of orle =2.7 (vs other random effects 0.7 and 0.9)
summary(m3) #std deviation of orle =2.5 
summary(m4) #std deviation of orle =2.5 

fittedModel <- m2
simulationOutput <- simulateResiduals(fittedModel = fittedModel)
plot(simulationOutput)

#----
fittedModel <- m3
simulationOutput <- simulateResiduals(fittedModel = fittedModel)
plot(simulationOutput)

#----
fittedModel <- m4
simulationOutput <- simulateResiduals(fittedModel = fittedModel)
plot(simulationOutput)


#issues with singularity so we'll go with a glm
total_pollen_mod<- glm(sum_prop2 ~ plant_type +year, family = Gamma,data=total)
plot(total_pollen_mod)
summary(total_pollen_mod)
with(total,boxplot(sum_prop2~plant_type))

fittedModel <- total_pollen_mod
simulationOutput <- simulateResiduals(fittedModel = fittedModel)
plot(simulationOutput)

range(total$sum_n)  # Actual response values Y range from 0 to 247

set.seed(1234567)
glmer_sim1 <- simulate(m2, nsim = 1000); 
glmer_sim2 <- simulate(m3, nsim = 1000)
hist(total$sum_n)
par(mfrow=c(2,3))
for(i in 1:6) hist(glmer_sim1[,i])
for(i in 1:6) hist(glmer_sim2[,i])

glm_sim1=simulate(total_pollen_mod,nsim=6)
hist(total$sum_prop2)
par(mfrow=c(2,3))
for(i in 1:6) hist(glm_sim1[,i])



set.seed(11)
n_perm=999




