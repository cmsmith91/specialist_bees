rm(list=ls())
library(tidyverse)
library(furrr)
library(phylolm)
setwd("/Users/colleen/Dropbox/github/specialist_bees")
select=dplyr::select; map =purrr::map

#upload data formatted for phyloglm analysis
data=read.csv("inat_analysis/df_inat_02nov2021.csv")
rownames(data)=data$X
data_df = data 

#get rid of wind pollinated plants
new_df = data_df[data_df$pollination=='biotic',] 
new_df$abund_scaled=as.vector(scale(new_df$effort_corrected_abund))

#upload and format phylogenetic tree of the data
angio_tree=readRDS("inat_analysis/tree_inat_angiosperms_30april2021.rds")
scenario3 = angio_tree$scenario.3
scenario3$node.label<-NULL
new_tree=keep.tip(scenario3,gsub(" ","_",new_df$species_underscore))


#upload data for the paired analysis
sisters=read_csv("inat_analysis/sisters_02nov2021.csv") %>% 
    filter(keep & !focal_a_crop) #get rid of plants that don't meet criteria

set.seed(100) #set seed for blinding the data: randomizing the predictor variables within pairs 

#change data into a long format (separate rows for each plant in a pair)
sisters_long=sisters %>% mutate(pair_id=1:n()) %>% 
    select(focal_genus,sister_genus,dist,pair_id) %>%
    pivot_longer(!dist & !pair_id,names_to='type', values_to='plant_name') %>%
    select(-type) %>% select(plant_name,pair_id,dist) 

#make a vector of all focal genera in the paired analysis, then separate by family pair and genera pairs
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
paired_df=paired_df %>% mutate(abund_log10=log10(effort_corrected_abund))
#run paired t-test

#reformat the data to visualize if differences are normally distributed
test_mapfun=paired_df %>% split(.$pair_id)
df=test_mapfun[[30]]
df_differences=paired_df %>% split(.$pair_id) %>% map_dfr(function(df) {
    spec=df[df$category=='specialist_host',]$effort_corrected_abund
    not_spec=df[df$category=='nonhost',]$effort_corrected_abund
    
    the_diff=df[df$category=='specialist_host',]$effort_corrected_abund-df[df$category=='nonhost',]$effort_corrected_abund
    the_diff_log=log10(df[df$category=='specialist_host',]$effort_corrected_abund)-log10(df[df$category=='nonhost',]$effort_corrected_abund)
    percent_diff=(spec-not_spec)/spec*100
    data.frame(pair_id=df$pair_id[1],diff=the_diff,diff_log_abund=the_diff_log,percent_diff=percent_diff)
    
    })
hist(df_differences$percent_diff)
median(df_differences$percent_diff)
paired_df %>% split(.$pair_id) %>% map_dfr(function(df) {
    
    
    data.frame(specialist_host=df[df$category=='specialist_host',]$genus,nonhost=df[df$category=='nonhost',]$genus)
    
})

# pdf('figures/histograms_differences.pdf',width=11)
par(mfrow=c(1,2))
hist(df_differences$diff,xlab='difference between pairs in abundance \n(raw data)',main='')
hist(df_differences$diff_log_abund,xlab='difference between pairs in abundance \n(log-transformed data)',main='') #looks better - though this may change upon unblinding data though
# dev.off()

library(ggpubr)
ggqqplot(df_differences$diff)
shapiro.test(df_differences$diff)
ggqqplot(df_differences$diff_log_abund) #looks ok
shapiro.test(df_differences$diff_log_abund) #not sig different from a normal distribution
hist(df_differences$diff_log_abund)

#reformat the data so they're in a wide foramt for the t-test
pairs_wide=paired_df  %>%
    select(category,pair_id,abund_log10) %>% 
    pivot_wider(names_from=category,values_from=abund_log10)
t.test(Pair(specialist_host, nonhost) ~ 1, ,alternative='greater',data = pairs_wide)
less_or_more=ifelse(mean(df_differences$diff)<0,'less','more')


print(paste0('Plants hosting specialist bees were, on average, ', 
             abs(round(mean(df_differences$diff),3)*100), '% ',
             less_or_more,' abundant than close relatives not hosting specialist bees (median percent difference = ',round(median(df_differences$diff),3)*100,"%)"))
df_differences

#need to report diffs on a raw scale
host_abund=mean(paired_df[paired_df$category=='specialist_host',]$effort_corrected_abund)
nhost_abund=mean(paired_df[paired_df$category=='nonhost',]$effort_corrected_abund)

(host_abund-nhost_abund)/(nhost_abund)*100 #percent increase
(host_abund)/(nhost_abund)#times abundant



#what is the percent difference for the pair with the median difference?
diff=pairs_wide$specialist_host-pairs_wide$nonhost
med_diff=median(diff)
pairs_wide[which(diff==med_diff),]
middle_pair=paired_df %>% filter(pair_id==22)
host_abund2=middle_pair[middle_pair$category=='specialist_host',]$effort_corrected_abund
nonhost_abund2=middle_pair[middle_pair$category=='nonhost',]$effort_corrected_abund

(host_abund2-nonhost_abund2)/(nonhost_abund2)*100 #percent increase
(host_abund2)/(nonhost_abund2) #percent increase

#make transparent color for figure
transparent_col=rgb(225,225,225,max=255,alpha=0)
col2rgb(4)
percent_col=40
my_alpha=( percent_col) * 255 / 100
my_alpha2=60*225/100
blue_transparent=rgb(34,151,230,max=255,alpha=my_alpha)
blue_transparent2=rgb(34,151,230,max=255,alpha=my_alpha2)

loc_nonhost=1.2; loc_host=1.8
par(mfrow=c(1,1))

# pdf('figures/paired_analysis_boxplot_25feb2022.pdf')
par(mar=c(4.5,5.5,3,2),cex.lab=1.5,cex.axis=1.4,cex=1.5)
with(paired_df,
     stripchart(abund_log10~(category),
                ylab=expression('log'[10]*'(effort-corrected abundance)'),
                xlab='plant type',
                group.names=c('nonhost','host'),
                vertical=T,pch=16,col=blue_transparent2,cex=.5,at=c(loc_nonhost,loc_host)))
for(i in 1:nrow(pairs_wide)){segments(loc_nonhost, pairs_wide$nonhost[i], loc_host, pairs_wide$specialist_host[i],lty=2,col=blue_transparent)}
with(paired_df,boxplot(abund_log10~category,boxwex=c(.35,.35),
                               xaxt = "n" ,xlab='Plant type',pch=1,col=transparent_col,alpha=.1,at=c(loc_nonhost,loc_host),add=T))
# dev.off()


 # pdf('figures/paired_analysis_boxplot_raw.pdf')
par(mar=c(4.5,5.5,3,2),cex.lab=1.5,cex.axis=1.4,cex=1.5)
with(paired_df,
     stripchart(effort_corrected_abund~category,
                ylab='log10(effort-corrected abundance)',
                xlab='Plant type',
                group.names=c('non-host','host'),
                vertical=T,pch=16,col=4,cex=.8,at=c(loc_nonhost,loc_host)))
for(i in 1:nrow(pairs_wide)){segments(loc_nonhost, 10^pairs_wide$nonhost[i], loc_host, 10^pairs_wide$specialist_host[i],lty=2,col=blue_transparent)}
with(paired_df,boxplot(effort_corrected_abund~category,boxwex=c(.35,.35),
                       xaxt = "n" ,xlab='Plant type',pch=1,col=transparent_col,alpha=.1,at=c(loc_nonhost,loc_host),add=T))
 # dev.off()

 paired_df %>% group_by(category) %>% summarize(med= median(effort_corrected_abund))
 
 paired_df %>% group_by(category) %>% summarize(med= median(abund_scaled))
data.frame(pairs_wide)

#summary of new_df
paste0(sum(new_df$raw_abund),' observations')
paste0(nrow(new_df),' plant genera') #number of genera


#run all angiosperm analysis on blinded data
obs_mod=phyloglm(category_binary~abund_scaled,phy=new_tree,data=new_df,btol=40,method = c("logistic_MPLE"))
summary(obs_mod)
obs_beta=coef(obs_mod)[2]

# run permutations to assess the significance of model coefficients
n_perm=999
plan(multisession,workers=6)
permut_output=1:n_perm %>% future_map(possibly(function(i){
    new_y=sample(new_df$category_binary)
    mod=phyloglm(new_y~abund_scaled,phy=new_tree,data=new_df,btol=30,method = c("logistic_MPLE"))
    null_beta=coef(mod)[2]
    converged=mod$convergence==0
    # 
    tibble(null_beta=null_beta,converged=converged)
},otherwise=NULL),.options=furrr_options(seed=88)) %>% bind_rows

# check whether number of permutations =n_perm
# if not, do more permtuations in units of 10
# until you get the correct number of converged permutations
while(sum(permut_output$converged)<n_perm){
    add_me=1:10 %>% future_map(possibly(function(i){
        new_y=sample(new_df$category_binary)
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

#bootstrap confidence intervals for the slope
#resample with replacement from the data 1000 times
boot_n=2000
plan(multisession(workers=6))
boot_output=1:boot_n %>% future_map(function(i){
    indices=sample(1:nrow(new_df),size=nrow(new_df),replace=T)
    boot_data=new_df[indices,]
    boot_glm=phyloglm(category_binary~abund_scaled,phy=new_tree,data=boot_data,btol=40,method = c("logistic_MPLE"))
    coefs=boot_glm$coefficients 
    names(coefs)<-c('intercept','abund_scaled')
    return(coefs)
}, .options = furrr_options(seed = 175456))

#order the bootoutput
boot_arranged=boot_output %>% bind_rows %>% arrange(abund_scaled)
lower_coefs=boot_arranged[50,]
lower_intercept=as.numeric(lower_coefs[1]);lower_slope=as.numeric(lower_coefs[2])

upper_coefs=boot_arranged[1950,]
upper_intercept=as.numeric(upper_coefs[1]);upper_slope=as.numeric(upper_coefs[2])


#now plot the model prediction and results
# get probabilities from the model coefficients
obs_intercept=coef(obs_mod)[1]
prob_host=exp(obs_intercept+obs_beta*new_df$abund_scaled)/(1+exp(obs_intercept+obs_beta*new_df$abund_scaled))
prob_host=prob_host[order(prob_host)]
abund_ordered=new_df$abund_scaled[order(new_df$abund_scaled)]

#get y-vals for upper and lower ci lines
lower_y=exp(lower_intercept+lower_slope*new_df$abund_scaled)/(1+exp(lower_intercept+lower_slope*new_df$abund_scaled))
upper_y=exp(upper_intercept+upper_slope*new_df$abund_scaled)/(1+exp(upper_intercept+upper_slope*new_df$abund_scaled))
ci_df=data.frame(upper_y=upper_y,lower_y=lower_y,abund_scaled=new_df$abund_scaled) %>% arrange(abund_scaled)


#get effect sizes
percent_increase=round((max(prob_host)-min(prob_host))/min(prob_host)*100,1)
percent_decrease=round((max(prob_host)-min(prob_host))/max(prob_host)*100,1)

if(obs_beta>0){
    print(paste0('The least abundant angiosperm in our data had an ',round(min(prob_host)*100),'% chance of hosting pollen specialist bee species, whereas the most abundant angiosperm in our data had an ', round(max(prob_host)*100),'% chance.'))
    }
if(obs_beta<0){
    print(paste0('The least abundant angiosperm in our data had an ',round(max(prob_host)*100),'% chance of hosting pollen specialist bee species, whereas the most abundant angiosperm in our data had an ', round(min(prob_host)*100),'% chance.'))
}

#basic stats about sample size
sum(new_df$raw_abund)
nrow(paired_df)
hosts=read_csv("~/Dropbox/Fall_2014/postdoc/inat project/data_current/bee_hosts_allUS.csv") %>% 
    filter(plant_rank=="genus" & region=='east')
data.frame(hosts %>% filter(host_plant %in% paired_df$genus) %>% group_by(host_plant) %>% 
    summarize(n_distinct(bee)))
paired_df%>% filter(category=="specialist_host" & !genus %in% hosts$host_plant)

# pdf('figures/phyloglm_results_22nov2021.pdf',width=8)
par(cex.lab=1.8,cex.axis=1.5,mar=c(4.8,5.1,3.1,2.1))
plot(abund_ordered,prob_host,type='n',ylim=c(0,1),ylab='Probability of hosting a specialist bee',xlab='Abundance (corrected for effort and scaled)')
for(vec in boot_output) {
    y_val=exp(vec[1]+vec[2]*new_df$abund_scaled)/(1+exp(vec[1]+vec[2]*new_df$abund_scaled))
    y_val=y_val[order(y_val)]
    lines(x=ci_df$abund_scaled,y_val,col=adjustcolor('gray',alpha=.5))
}
lines(abund_ordered,prob_host,type = "l",lwd=3)
with(new_df,points(abund_scaled,category_binary,pch="|",cex=.6))

#add boot-strapped cis to visualize uncertainty
boot_output
abund_seq=seq(min(new_df$abund_scaled),max(new_df$abund_scaled),length.out=1000)
boot_pred=boot_output %>% map(function(coef_vec){
    intercept=coef_vec[1];slope=coef_vec[2]
    data.frame(y_pred=exp(intercept+slope*abund_seq)/(1+exp(intercept+slope*abund_seq)))
})
boot_pred_df=boot_pred %>% bind_cols
boot_pred_cis=apply(boot_pred_df,1,function(x) x[order(x)][c(50,1950)])

#add predictions to graph as well
prob_host_seq=exp(obs_intercept+obs_beta*abund_seq)/(1+exp(obs_intercept+obs_beta*abund_seq))

new_ci_df=data.frame(t(boot_pred_cis)) %>% rename(lower_ci=X1,upper_ci=X2) %>%
    mutate(abund=abund_seq,prob_host=prob_host_seq)
ci_df

# dev.off()
#which are the 10 most common plants in the data
head(new_df)
new_df_arranged=new_df %>% arrange(desc(abund_scaled))
new_df_arranged[1:20,]

break_n=40
abund_min=min(new_df$abund_scaled); abund_max=max(new_df$abund_scaled)
my_breaks=seq(abund_min,abund_max,by=(abund_max-abund_min)/break_n) 

densities_df=new_df %>% split(.$category_binary) %>% map_dfr(function(df){
    dens=hist(df$abund_scaled,plot=F,breaks=my_breaks)
    percent_dens=  dens$density/sum(dens$density)
    
    
    data.frame(abund_scaled=dens$mid,percent_dens=percent_dens) %>% 
        mutate(category_binary=df$category_binary[1])
    
}) 

hist_df=densities_df %>% mutate(pct=ifelse(category_binary,1-percent_dens,percent_dens))
my_cols=RColorBrewer::brewer.pal(8,"Accent")

field_plants=c('Salix','Claytonia','Vaccinium','Cucurbita','Cornus','Lysimachia','Pontederia','Solidago')
new_df %>% filter(genus %in% field_plants)  %>% arrange(abund_scaled)

#get conf intervals around line
summary(obs_mod)

 # pdf('figures/new_abund_plot_08dec2021.pdf')
ggplot() +
    geom_segment(data=hist_df[hist_df$category_binary==1,], size=4, show.legend=FALSE,colour=my_cols[1],
                 aes(x=abund_scaled, xend=abund_scaled, y=category_binary, yend=pct)) +
    geom_segment(data=hist_df[hist_df$category_binary==0,], size=4, show.legend=FALSE,colour=my_cols[2],
                 aes(x=abund_scaled, xend=abund_scaled, y=category_binary, yend=pct))+
    geom_segment(dat=new_df[new_df$category_binary==0,], aes(x=abund_scaled, xend=abund_scaled, y=0, yend=-0.02), size=0.2, colour="grey30") +
    geom_segment(dat=new_df[new_df$category_binary==1,] %>% filter(!genus %in% field_plants), aes(x=abund_scaled, xend=abund_scaled, y=1, yend=1.02), size=0.2, colour="grey30") +
    geom_segment(dat=new_df[new_df$category_binary==1,] %>% filter(genus %in% field_plants), aes(x=abund_scaled, xend=abund_scaled, y=1, yend=1.02), size=0.2, colour=my_cols[6]) +
    # geom_line(data=data.frame(x=abund_ordered, 
    #                           y=prob_host), 
    #           aes(x,y), colour="black", lwd=1) +
    geom_ribbon(data=new_ci_df,aes(x=abund_seq,ymin = lower_ci, ymax = upper_ci), fill = adjustcolor("grey90",.5)) + 
    geom_line(data=new_ci_df,aes(x=abund_seq,y=prob_host),color='black',lwd=1)+
    #geom_line(data=new_ci_df,aes(abund_seq,upper_ci),color='gray',lwd=1)+
    #geom_line(data=new_ci_df,aes(abund_seq,lower_ci),color='gray',lwd=1)+
    theme_bw(base_size=12)+theme(
        # Hide panel borders and remove grid lines
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Change axis line
        axis.line = element_line(colour = "black"),
        text = element_text(size=20)
    )+
    labs(x='abundance (corrected for effort and scaled)',y='probability of hosting a specialist bee')

 # dev.off()


