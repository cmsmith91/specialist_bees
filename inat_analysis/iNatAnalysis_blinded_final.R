rm(list=ls())
setwd("inat_analysis")
library(tidyverse)
library(furrr)
library(phylolm)
library(lme4)
select=dplyr::select; map =purrr::map

#upload data formatted for phyloglm analysis
data=read.csv("df_inat_02nov2021.csv")
rownames(data)=data$X
data_df = data 

#get rid of wind pollinated plants
new_df = data_df[data_df$pollination=='biotic',] 
new_df$abund_scaled=as.vector(scale(new_df$effort_corrected_abund))


#upload and format phylogenetic tree of the data
angio_tree=readRDS("tree_inat_angiosperms_30april2021.rds")
scenario3 = angio_tree$scenario.3
scenario3$node.label<-NULL
new_tree=keep.tip(scenario3,gsub(" ","_",new_df$species_underscore))

#simulate data for running the blinded analysis
my_X=cbind(rep(1,nrow(new_df)),new_df$abund_scaled)
row.names(my_X)=sub(" ","_",new_df$species_underscore)
use_me=read_csv("rbintrait_model_coefficients.csv") %>% 
    filter(focal_alpha)
beta1=0.3
(phylo_sig=unique(use_me$alpha)[2])
beta0=use_me[use_me$alpha==phylo_sig & use_me$beta1==beta1,]$beta0
set.seed(56)#Set seed
ys=rbinTrait(n=1,phy=new_tree,beta=c(beta0,beta1),alpha=phylo_sig, X=my_X)
y_val=as.vector(ys)
new_df$blind_y=y_val
with(new_df,boxplot(abund_scaled~blind_y))
table(new_df$category)

#upload data for the paired analysis
sisters=read_csv("sisters_02nov2021.csv") %>% 
    filter(keep & !focal_a_crop) #get rid of plants that don't meet criteria

#set seed for blinding the data: randomizing the predictor variables within pairs 
set.seed(100) 

#change data into a long format (separate rows for each plant in a pair)
sisters_long=sisters %>% mutate(pair_id=1:n()) %>% 
    select(focal_genus,sister_genus,dist,pair_id) %>%
    pivot_longer(!dist & !pair_id,names_to='type', values_to='plant_name') %>%
    select(-type) %>% select(plant_name,pair_id,dist) %>% 
    split(.$pair_id) %>% map_dfr(function(df){ #split by pair id, and randomize which plant is the host plant and which is a non-host
        df %>% 
            mutate(category_rand=sample(c('specialist_host','control'),2))})
        
#make a vector of all focal genera in the paired analysis
focal_genera=c(sisters$focal_genus,sisters$sister_genus)

#format the data for the paired analysis
paired_genera=new_df %>% 
    filter(genus %in% focal_genera) %>%  
    dplyr::select(genus,family,category,pollination,effort_corrected_abund,raw_abund) %>%
    mutate(plant_name=genus,rank='genus')
paired_df=paired_genera  %>% 
    left_join(sisters_long) %>%mutate(pair_id=as.factor(pair_id))

#run paired analysis
#add a column for ln abundance
paired_df=paired_df %>% 
  mutate(log_abund=log(effort_corrected_abund)) #natural log abundance

# code the mixed model below results in a singular fit: Probably because there are
# only 2 observations per level of random effect?
paired_mod=lmer(effort_corrected_abund~category_rand + (1| pair_id),data=paired_df)

plot(paired_mod)
boxplot(paired_df$effort_corrected_abund) #there's are some outliers
err=as.vector(residuals(paired_mod))
boxplot(err~paired_df$category_rand)
summary(paired_mod)

#other options:
#paired t-test

#reformat the data to visualize if differences are normally distributed
test_mapfun=paired_df %>% split(.$pair_id)
df=test_mapfun[[1]]
df_differences=paired_df %>% split(.$pair_id) %>% map_dfr(function(df) {
    
    the_diff=df[df$category_rand=='specialist_host',]$effort_corrected_abund-df[df$category_rand=='control',]$effort_corrected_abund
    the_diff_log=log10(df[df$category_rand=='specialist_host',]$effort_corrected_abund)-log10(df[df$category_rand=='control',]$effort_corrected_abund)
    
    data.frame(pair_id=df$pair_id[1],diff=the_diff,diff_log_abund=the_diff_log)
    
    })
hist(df_differences$diff)
hist(df_differences$diff_log_abund) #looks better - though this may change upon unblinding data though

library(ggpubr)
ggqqplot(df_differences$diff_log_abund) #looks ok
shapiro.test(df_differences$diff_log_abund) #not sig different from a normal distribution


#reformat the data so they're in a wide format for the t-test
pairs_wide=paired_df %>% mutate(abund_log10=log10(effort_corrected_abund)) %>%
    select(category_rand,pair_id,abund_log10) %>% 
    pivot_wider(names_from=category_rand,values_from=abund_log10)
t.test(Pair(specialist_host, control) ~ 1, ,alternative='greater',data = pairs_wide)
less_or_more=ifelse(mean(df_differences$diff)<0,'less','more')
print(paste0('Plants hosting specialist bees were, on average, ', abs(round(mean(df_differences$diff),3)*100), '% ',less_or_more,' abundant than close relatives not hosting specialist bees (median percent difference = ',round(median(df_differences$diff),3)*100,"%)"))

# if the data are not normally distributed after unblinding, 
# we will use a wilcoxin signed rank test
wilcox.test(Pair(specialist_host, control) ~ 1, ,alternative='greater',data = pairs_wide)

#make the figure
paired_df_blinded=paired_df %>% mutate(abund_log10=log10(effort_corrected_abund))

#make transparent color for figure
transparent_col=rgb(225,225,225,max=255,alpha=0)
col2rgb(4)
percent_col=40
my_alpha=( percent_col) * 255 / 100
blue_transparent=rgb(34,151,230,max=255,alpha=my_alpha)
loc_nonhost=1.2; loc_host=1.8

# pdf('paired_analysis_boxplot_02nov2021.pdf')
par(mar=c(4.5,5.5,3,2),cex.lab=1.5,cex.axis=1.4,cex=1.5)
with(paired_df_blinded,
     stripchart(abund_log10~category_rand,
                ylab='log10(effort-corrected abundance)',
                xlab='Plant type',
                group.names=c('non-host','host'),
                vertical=T,pch=16,col=4,cex=.8,at=c(loc_nonhost,loc_host)))
for(i in 1:nrow(pairs_wide)){segments(loc_nonhost, pairs_wide$control[i], loc_host, pairs_wide$specialist_host[i],lty=2,col=blue_transparent)}
with(paired_df_blinded,boxplot(abund_log10~category_rand,boxwex=c(.35,.35),
                               xaxt = "n" ,xlab='Plant type',pch=1,col=transparent_col,alpha=.1,at=c(loc_nonhost,loc_host),add=T))
# dev.off()


#run blind analysis on all angiosperm genera 
obs_mod=phyloglm(blind_y~abund_scaled,phy=new_tree,data=new_df,btol=40,method = c("logistic_MPLE"))
summary(obs_mod)
obs_beta=coef(obs_mod)[2]

# run permutations to assess the significance of model coefficients
n_perm=999
plan(multisession,workers=6)
permut_output=1:n_perm %>% future_map(possibly(function(i){
    new_y=sample(y_val)
    mod=phyloglm(new_y~abund_scaled,phy=new_tree,data=new_df,btol=30,method = c("logistic_MPLE"))
    null_beta=coef(mod)[2]
    converged=mod$convergence==0
    # 
    tibble(null_beta=null_beta,converged=converged)
},otherwise=NULL),.options=furrr_options(seed=88)) %>% bind_rows

# check whether number of permutations = n_perm
# (will likely not be becuase the model doesn't always converge)
# if not, I'll do more permutations in units of 10
# until I get the correct number of converged permutations
while(sum(permut_output$converged)<n_perm){
    add_me=1:10 %>% future_map(possibly(function(i){
        new_y=sample(y_val)
        mod=phyloglm(new_y~abund_scaled,phy=new_tree,data=new_df,btol=30,method = c("logistic_MPLE"))
        null_beta=coef(mod)[2]
        converged=mod$convergence==0
        # 
        tibble(null_beta=null_beta,converged=converged)
    },otherwise=NULL),.options=furrr_options(seed=765)) %>% bind_rows
    permut_output=permut_output %>% bind_rows(add_me)
}

#calculate pval from permutations
null_betas=as.vector(permut_output$null_beta[1:n_perm])
r=sum(null_betas>obs_beta)
(pval=(r+1)/(n_perm+1))

#now plot the model prediction and results
# get probabilities from the model coefficients
obs_intercept=coef(obs_mod)[1]
prob_host=exp(obs_intercept+obs_beta*new_df$abund_scaled)/(1+exp(obs_intercept+obs_beta*new_df$abund_scaled))
prob_host=prob_host[order(prob_host)]
abund_ordered=new_df$abund_scaled[order(new_df$abund_scaled)]

#get effect sizes
percent_increase=round((max(prob_host)-min(prob_host))/min(prob_host)*100,1)
percent_decrease=round((max(prob_host)-min(prob_host))/max(prob_host)*100,1)

if(obs_beta>0){
    print(paste0('The least abundant angiosperm in our data had an ',round(min(prob_host)*100),'% chance of hosting pollen specialist bee species, whereas the most abundant angiosperm in our data had an ', round(max(prob_host)*100),'% chance.'))
    }
if(obs_beta<0){
    print(paste0('The least abundant angiosperm in our data had an ',round(max(prob_host)*100),'% chance of hosting pollen specialist bee species, whereas the most abundant angiosperm in our data had an ', round(min(prob_host)*100),'% chance.'))
    }

# pdf('phyloglm_results_blinded.pdf',width=8)
par(cex.lab=1.8,cex.axis=1.5,mar=c(4.8,5.1,3.1,2.1))
plot(abund_ordered,prob_host,type='n',ylim=c(0,1),ylab='Probability of hosting a specialist bee',xlab='Abundance (corrected for effort and scaled)')
lines(abund_ordered,prob_host,type = "l",lwd=3)
points(new_df$abund_scaled,y_val,pch="|")
# dev.off()
