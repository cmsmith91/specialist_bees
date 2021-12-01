rm(list=ls())
library(tidyverse)
library(furrr)
library(phylolm)
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
df=test_mapfun[[1]]
df_differences=paired_df %>% split(.$pair_id) %>% map_dfr(function(df) {
    
    the_diff=df[df$category=='specialist_host',]$effort_corrected_abund-df[df$category=='nonhost',]$effort_corrected_abund
    the_diff_log=log10(df[df$category=='specialist_host',]$effort_corrected_abund)-log10(df[df$category=='nonhost',]$effort_corrected_abund)
    
    data.frame(pair_id=df$pair_id[1],diff=the_diff,diff_log_abund=the_diff_log)
    
    })
hist(df_differences$diff)
hist(df_differences$diff_log_abund) #looks better - though this may change upon unblinding data though
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


print(paste0('Plants hosting specialist bees were, on average, ', abs(round(mean(df_differences$diff),3)*100), '% ',less_or_more,' abundant than close relatives not hosting specialist bees (median percent difference = ',round(median(df_differences$diff),3)*100,"%)"))
df_differences
#need to report diffs on a raw scale
mean(pairs_wide$specialist_host)
mean(pairs_wide$nonhost)


#make transparent color for figure
transparent_col=rgb(225,225,225,max=255,alpha=0)
col2rgb(4)
percent_col=40
my_alpha=( percent_col) * 255 / 100
blue_transparent=rgb(34,151,230,max=255,alpha=my_alpha)
loc_nonhost=1.2; loc_host=1.8

# pdf('figures/paired_analysis_boxplot_16nov2021.pdf')
par(mar=c(4.5,5.5,3,2),cex.lab=1.5,cex.axis=1.4,cex=1.5)
with(paired_df,
     stripchart(abund_log10~category,
                ylab='log10(effort-corrected abundance)',
                xlab='Plant type',
                group.names=c('non-host','host'),
                vertical=T,pch=16,col=4,cex=.8,at=c(loc_nonhost,loc_host)))
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

# pdf('figures/phyloglm_results_22nov2021.pdf',width=8)
par(cex.lab=1.8,cex.axis=1.5,mar=c(4.8,5.1,3.1,2.1))
plot(abund_ordered,prob_host,type='n',ylim=c(0,1),ylab='Probability of hosting a specialist bee',xlab='Abundance (corrected for effort and scaled)')
lines(abund_ordered,prob_host,type = "l",lwd=3)
with(new_df,points(abund_scaled,category_binary,pch="|",cex=.6))
# dev.off()

densities_df=new_df %>% split(.$category_binary) %>% map_dfr(function(df){
    dens=hist(df$abund_scaled,plot=F,breaks=my_breaks)
    percent_dens=  dens$density/sum(dens$density)
    
    
    data.frame(abund_scaled=dens$mid,percent_dens=percent_dens) %>% 
        mutate(category_binary=df$category_binary[1])
    
}) %>% mutate(my_alphas=255*(percent_dens+(1-max(percent_dens)))) %>% #make the max value 1 
    mutate(col=rgb(34,151,230,max=255,alpha=my_alphas)) #have colors be  transparencnt basded on density


hist_df=densities_df %>% mutate(pct=ifelse(category_binary,1-percent_dens,percent_dens))
my_cols=RColorBrewer::brewer.pal(3,"Accent")
pdf('figures/new_abund_plot.pdf')
ggplot() +
    geom_segment(data=h[h$category_binary==1,], size=4, show.legend=FALSE,colour=my_cols[1],
                 aes(x=abund_scaled, xend=abund_scaled, y=category_binary, yend=pct)) +
    geom_segment(data=h[h$category_binary==0,], size=4, show.legend=FALSE,colour=my_cols[2],
                 aes(x=abund_scaled, xend=abund_scaled, y=category_binary, yend=pct))+
    geom_segment(dat=new_df[new_df$category_binary==0,], aes(x=abund_scaled, xend=abund_scaled, y=0, yend=-0.02), size=0.2, colour="grey30") +
    geom_segment(dat=new_df[new_df$category_binary==1,], aes(x=abund_scaled, xend=abund_scaled, y=1, yend=1.02), size=0.2, colour="grey30") +
    geom_line(data=data.frame(x=abund_ordered, 
                              y=prob_host), 
              aes(x,y), colour="grey50", lwd=1) +
    theme_bw(base_size=12)+theme(
        # Hide panel borders and remove grid lines
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Change axis line
        axis.line = element_line(colour = "black")
    )+
    labs(x='abundance (corrected for effort and scaled)',y='probability of hosting a specialist bee')

dev.off()


#try doing the color gradient 
my_alphas=dens_vec*255
transparent_cols=rgb(34,151,230,max=255,alpha=my_alphas)
plot(col=transparent_cols)
b=new_df %>% split(.$category_binary)
df= b[[1]]
break_n=40
abund_min=min(new_df$abund_scaled); abund_max=max(new_df$abund_scaled)
my_breaks=seq(abund_min,abund_max,by=(abund_max-abund_min)/break_n) 

densities_df=new_df %>% split(.$category_binary) %>% map_dfr(function(df){
    dens=hist(df$abund_scaled,plot=F,breaks=my_breaks)
    percent_dens=  dens$density/sum(dens$density)

    
    data.frame(abund_scaled=dens$mid,percent_dens=percent_dens) %>% 
        mutate(category_binary=df$category_binary[1])
    
}) %>% mutate(my_alphas=255*(percent_dens+(1-max(percent_dens)))) %>% #make the max value 1 
    mutate(col=rgb(34,151,230,max=255,alpha=my_alphas)) #have colors be  transparencnt basded on density

densities_df[densities_df$percent_dens==0,]$col='white'
hist(log(densities_df$percent_dens))
min_dens=min(densities_df[densities_df$percent_dens!=0,]$percent_dens)
max_dens=max(densities_df$percent_dens)
new_col_seq=seq(min_dens,max_dens,by=(max_dens-min_dens)/8)
col_gradient=RColorBrewer::brewer.pal(9,'YlOrRd')
col_gradient1=RColorBrewer::brewer.pal(9,'YlOrRd')

col_gradient2=col_gradient[c(3,5,7,9)]

densities_df$new_col=col_gradient[1]
for(i in 2:length(col_gradient1)){
    densities_df[densities_df$percent_dens>=new_col_seq[i],]$new_col=col_gradient[i]
    
}
densities_df[densities_df$percent_dens==0,]$new_col='white'
densities_df$line_col=ifelse(densities_df$new_col=='white',"white",'black')

new_col_seq2=seq(min_dens,max_dens,by=(max_dens-min_dens)/length(col_gradient2))

mid_pts_density=c(mean(c(new_col_seq2[1],new_col_seq2[2])))
labs_vec=c(paste0(round(new_col_seq2[1],3)," - ",round(new_col_seq2[2],2)))

densities_df$new_col2=col_gradient2[1]
for(i in 2:length(col_gradient2)){
    densities_df[densities_df$percent_dens>=new_col_seq2[i],]$new_col2=col_gradient2[i]
    mid_pts_density[i]=mean(c(new_col_seq2[i],new_col_seq2[i+1]))
    labs_vec[i]=paste0(round(new_col_seq2[i],2)," - ",round(new_col_seq2[i+1],2))
    
    
}
densities_df[densities_df$percent_dens==0,]$new_col2='white'

# pdf('figures/phyloglm_results_22nov2021.pdf')
par(cex.lab=1.8,cex.axis=1.5,xpd=T,mar=c(4,5,5,9),mfrow=c(1,1))
plot(abund_ordered,prob_host,type='n',ylim=c(0,1),ylab='Probability of hosting a specialist bee',xlab='Abundance (corrected for effort and scaled)')
lines(abund_ordered,prob_host,type = "l",lwd=3)
with(densities_df,points(abund_scaled,category_binary,pch=22,cex=1.2,bg=new_col2,col=new_col2))
with(densities_df,legend('topright',title='Proportion \nof observations',legend=labs_vec,pch=22,cex=1.2,bty='n',
                         pt.bg=col_gradient2,col=line_col,inset=c(-0.4,0)))
# dev.off()


#old
col_df=data.frame(col=c(unique(densities_df[densities_df$percent_dens!=0,]$col),'white'),
           new_col=c(RColorBrewer::brewer.pal(8,'YlOrRd'),"white") )

densities_df_col=densities_df %>% left_join(col_df)

par(cex.lab=1.8,cex.axis=1.5,mar=c(4.8,5.1,3.1,2.1))
plot(abund_ordered,prob_host,type='n',ylim=c(0,1),ylab='Probability of hosting a specialist bee',xlab='Abundance (corrected for effort and scaled)')
lines(abund_ordered,prob_host,type = "l",lwd=3)
with(densities_df_col,points(abund_scaled,category_binary,pch=15,cex=.7,col=new_col))

RColorBrewer::brewer.pal(9,'YlOrRd')
unique(densities_df$col)
    
    
