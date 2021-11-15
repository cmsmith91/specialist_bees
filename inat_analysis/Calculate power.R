rm(list=ls())
library(tidyverse)
library(furrr)
library(phyclust) #for calculating tree hieght
source("error_analysis_functions.R")
setwd("~/Dropbox/Fall_2014/postdoc/inat project/data_current")

#upload data formatted for phyloglm analysis
data=read.csv("/Users/colleen/Dropbox/github/inat_project/processed_data/df_phyloglm_01oct2021.csv")
data$abund_scaled=as.vector(scale(data$effort_corrected_abund))
rownames(data)=data$X
data_df = data 
new_df = data_df[data_df$pollination=='biotic',] 
new_df$abund_scaled=as.vector(scale(new_df$effort_corrected_abund))

#what percentage of plant genera in the data host pollen specialist bees?
paste0(round(mean(new_df$category_binary)*100,2),'% of the plants in our data host pollen specialist bee species')

#upload phylogenetic tree of the data
angio_tree = readRDS("/Users/colleen/Dropbox/github/inat_project/processed_data/tree_inat_angiosperms_30april2021.rds")
scenario3 = angio_tree$scenario.3
scenario3$node.label<-NULL
new_tree=keep.tip(scenario3,gsub(" ","_",new_df$species_underscore))


#  bet the model intercepts to use: vary beta0 between zero and one
b0s=seq(-2.85,-1.65,.025)
b1s=c(0,0.3)
tree_height=get.rooted.tree.height(new_tree)
focal_alphas=c(log(2)/(.5*tree_height),0.1) #alpha where half life is .5 the height of the tree, and alpha =0.1
vary_alphas=seq(log(2)/(tree_height),0.393,by=(0.393-log(2)/(.5*tree_height))/10)
alphas=c(focal_alphas,vary_alphas)
xs=new_df$abund_scaled

n_paramcombos=5000

my_X=cbind(rep(1,nrow(new_df)),new_df$abund_scaled)
row.names(my_X)=sub(" ","_",new_df$species_underscore)

# below is the code for the parameter exploration to pick the different values of the 
# model intercept for each parameter combination

# plan(multisession,workers=6)
# check_props=b0s %>% map(function(beta0){
#   b1s %>% map(function(beta1){
#     alphas %>% future_map(function(my_alpha){
#         a=rbinTrait(n=n_paramcombos,phy=new_tree,beta=c(beta0,beta1),alpha=my_alpha, X=my_X)
# 
#         prop_hosts=colSums(a)/nrow(a)
# 
#         data.frame("prop_hosts"=prop_hosts,"beta0"=beta0,"beta1"=beta1,alpha=my_alpha)
# 
#     },.options=furrr_options(seed=T)) %>% bind_rows
#   }) %>% bind_rows
# }) %>% bind_rows
# future:::ClusterRegistry("stop") #stop any existing cores
# 
# # for(df in check_props){
# #     hist(df$prop_hosts)
# # }
# 
# # RColorBrewer::display.brewer.all()
# some_cols=c(RColorBrewer::brewer.pal(10,"Set3"),RColorBrewer::brewer.pal(10,"Paired"),
#             RColorBrewer::brewer.pal(9,"Pastel1"),RColorBrewer::brewer.pal(8,"Accent"),
#             RColorBrewer::brewer.pal(8,"Dark2"),RColorBrewer::brewer.pal(9,"Set1"),
#             RColorBrewer::brewer.pal(9,"Spectral"),RColorBrewer::brewer.pal(9,"YlOrRd"),
#             RColorBrewer::brewer.pal(9,"Blues"),RColorBrewer::brewer.pal(9,"Greens"),
#             RColorBrewer::brewer.pal(9,"Greys"),RColorBrewer::brewer.pal(9,"BrBG"),
#             RColorBrewer::brewer.pal(9,"PiYG"),RColorBrewer::brewer.pal(9,"BuPu"))
# 
# 
# 
# 
# prop_df=check_props #%>% bind_rows
# obs_value=mean(new_df$category_binary)
# prop_summ=prop_df %>% group_by(beta0,beta1,alpha) %>%
#     summarize(m_prop_hosts=median(prop_hosts),sd=sd(prop_hosts))  %>%
#     rename(prop_hosts=m_prop_hosts) %>%
#     mutate(distance=abs((obs_value-prop_hosts)))
# 
# prop_summ$col=some_cols[as.factor(prop_summ$beta0)]
# par(mfrow=c(1,2))
# for(df in prop_summ %>% split(.$beta1)){
#     with(df,plot(alpha,prop_hosts,col=col,pch=16,main=paste0('beta1 = ',beta1[1])))
#     abline(h=mean(new_df$category_binary))
#     
# }
# 
# #add legend
# col_df=prop_summ %>% distinct(beta0,col)
# #with(col_df,plot(rep(1,nrow(col_df))~beta0,col=col,pch=16,cex=2))
# par(mfrow=c(1,1))
# plot.new()
# legend('center',legend=unique(col_df$beta0)[1:floor(n_distinct(b0s)/2)],col=unique(col_df$col)[1:floor(n_distinct(b0s)/2)],pch=16 )
# legend('right',legend=unique(col_df$beta0)[(floor(n_distinct(b0s)/2)+1):n_distinct(b0s)],col=unique(col_df$col)[(floor(n_distinct(b0s)/2)+1):n_distinct(b0s)],pch=16 )
# 
# 
# min_dists=prop_summ %>% group_by(beta1,alpha) %>% summarize(min_score=min(distance)) %>%
#     mutate(combo=paste0(beta1,alpha,min_score))
# use_me=prop_summ %>%
#     filter(paste0(beta1,alpha,distance) %in% min_dists$combo) %>%
#     arrange(beta1,alpha) %>%
#     mutate(beta1_alpha=paste0(beta1,"_",alpha))
# set.seed(454) # if multiple intercept values are of the same minimum distance randomly pick 1:
# use_me=use_me %>% split(.$beta1_alpha) %>% map_dfr(function(df){
#     if(nrow(df) !=1) {df=df[sample(nrow(df),1),]}
#     return(df)
# } ) %>%
#     mutate(focal_alpha=alpha %in% focal_alphas)
# use_me$b1_index=1:nrow(use_me)

# write_csv(use_me,"/Users/colleen/Dropbox/github/inat_project/processed_data/coefficient_model_siml_rbintrait_variedalpha.csv")
use_me=read_csv("/Users/colleen/Dropbox/github/inat_project/processed_data/coefficient_model_siml_rbintrait_variedalpha.csv")

#calculate the power when b1 = 0.3
#generate the data
N=130
start_alpha_val='null'

focal_use_me=use_me %>% filter(focal_alpha)
mod_ls=focal_use_me%>% split(1:sum(use_me$focal_alpha)) %>%
    map(function(df){
        beta0= df$beta0
        beta1=df$beta1
        phylo_sig=df$alpha
        
        #simulate response using binaryPGLMM.sim function
        y_mat=rbinTrait(n=N,phy=new_tree,beta=c(beta0,beta1),alpha=phylo_sig, X=my_X)
        
        #analyze the simulated datasets with phloglm
        phyloglm_wald(y_mat,new_df$abund_scaled,new_tree,n_reps=100,new_df,'logistic_MPLE',start_alpha=start_alpha_val)
        
        
    })

#get the medians of alpha:
median_alphas=mod_ls %>% map_dbl(function(phyloglm_output){
    median(phyloglm_output$mod_info[phyloglm_output$mod_info$converged,]$alpha[1:100])
    
})
range(median_alphas)
median(phyloglm_output_centered$mod_info[phyloglm_output_centered$mod_info$converged,]$alpha)

#plot histograms of the slopes
# pdf('slope_distributions.pdf')
par(mfrow=c(1,1),cex.lab=1.3,cex.axis=1.3)
for(i in 1:length(mod_ls)){
    true_slope=focal_use_me[i,]$beta1
    true_alpha=focal_use_me[i,]$alpha
    
    phyloglm_output=mod_ls[[i]]
    
    hist(phyloglm_output$mod_info[phyloglm_output$mod_info$converged,]$b1_est[1:100],
         main=paste0('b1 = ',true_slope,' \nalpha = ',round(true_alpha,3)),xlab='fitted slope estimate')
    abline(v=true_slope,col='red')
    
}
# dev.off()

mod_ls %>% map(function(ls){
    ls$analysis_summary
})

#now run error analysis
(my_alpha=focal_alphas[2]) #pick the parameter combination to test
(my_beta=b1s[2])

i=which(focal_use_me$beta1==my_beta & focal_use_me$alpha==my_alpha)
phyloglm_output=mod_ls[[i]]

# i=4
beta0= focal_use_me[i,]$beta0
beta1=focal_use_me[i,]$beta1
phylo_sig=focal_use_me[i,]$alpha
b1_index=focal_use_me[i,]$b1_index

# start error analysis 
folder_name='error_analysis_bigN_perms_26oct2021' # folder for saving the results
if(!dir.exists(folder_name)) dir.create(folder_name)

path=paste0(folder_name,'/phyloglm_typeiscaled_b1_',beta1,"_alpha_",round(my_alpha,3),"_startalpha_",ifelse(is.null(start_alpha_val),"NULL",start_alpha_val),'.rds')
start_time=proc.time()
permute_output=phyloglm_permute(phyloglm_output,999)
(time_elapsed_minutes=(proc.time()-start_time)[3]/60)
saveRDS(permute_output,path)

permute_output$analysis_summary
# permute_output=readRDS("error_analysis_bigN_perms_26oct2021/phyloglm_typeiscaled_b1index_4_alpha_0.007_startalpha.rds") #b1 of 0.4

hist(permute_output$permut_list %>% map_dbl(function(df) df$obs_beta[1]))

with(phyloglm_output$mod_info %>% filter(converged),hist(b1_est[1:100]))

#########
#how do the estimates for alpha look across the range of alpha tested?
other_use_me=use_me %>% filter(!focal_alpha)
big_mod_ls=other_use_me %>% split(1:sum(!use_me$focal_alpha))%>%
    map(function(df){
        beta0= df$beta0
        beta1=df$beta1
        phylo_sig=df$alpha
        
        #simulate response using binaryPGLMM.sim function
        y_mat=rbinTrait(n=N,phy=new_tree,beta=c(beta0,beta1),alpha=phylo_sig, X=my_X)
        
        #analyze the simulated datasets with phloglm
        phyloglm_wald(y_mat,new_df$abund_scaled,new_tree,n_reps=100,new_df,'logistic_MPLE',start_alpha=start_alpha_val)
        
        
    })
median_alphas=1:length(big_mod_ls) %>% map_dfr(function(i){
    phyloglm_output=big_mod_ls[[i]]
    fitted_alphas=phyloglm_output$mod_info[phyloglm_output$mod_info$converged,]$alpha[1:100]
    
    lower=quantile(fitted_alphas,.025)
    upper=quantile(fitted_alphas,.975)
    
    the_median=median(fitted_alphas)
    other_use_me[i,] %>% mutate(median_alpha=the_median,lower=lower,upper=upper)
})

range(median_alphas)
# pdf('plot_of_alphas.pdf')
par(mfrow=c(1,2),pch=16)
with(median_alphas %>% filter(beta1==0),
     plot(alpha,median_alpha,ylim=c(0,.45),ylab='alpha',main=paste0('beta1 = ',beta1[1])))
with(median_alphas %>% filter(beta1==0),lines(alpha,upper,col='red'))
with(median_alphas %>% filter(beta1==0),lines(alpha,lower,col='red'))

with(median_alphas %>% filter(beta1==.3),
     plot(alpha,median_alpha,ylim=c(0,.45),ylab='alpha',main=paste0('beta1 = ',beta1[1])))
with(median_alphas %>% filter(beta1==.3),lines(alpha,upper,col='red'))
with(median_alphas %>% filter(beta1==.3),lines(alpha,lower,col='red'))
# dev.off()
unique(median_alphas$alpha)
range(median_alphas$median_alpha)




#code looking at parameter combos one by one:
#     #get model intercept and slope coefficients to use
(my_alpha=focal_alphas[1])
(my_beta=b1s[2])

i=which(use_me$beta1==my_beta & use_me$alpha==my_alpha)
# i=4
beta0= use_me[i,]$beta0
beta1=use_me[i,]$beta1
phylo_sig=use_me[i,]$alpha

#get a measure of effect size
#if b1=0.4 what is the difference in probability between the most abundant plant and the least.
minmax_x=c(min(new_df$abund_scaled),max(new_df$abund_scaled))
(p=exp(beta0+beta1*minmax_x)/(1+exp(beta0+beta1*minmax_x)))
(p[2]-p[1])/(p[1])*100

#prob at model intercept
(exp(beta0)/(1+exp(beta0)))

#simulate response using binaryPGLMM.sim function
y_mat=rbinTrait(n=N,phy=new_tree,beta=c(beta0,beta1),alpha=phylo_sig, X=my_X)


#run phyloglm 
start_alpha_val='null'
# start_alpha_val=0.000132
# start_alpha_val=0.1

plan(multisession,workers=6)
ys=y_mat[,i]
mod=phyloglm(ys~xs,phy=new_tree,data=new_df,btol=40,method = c("logistic_MPLE"),start.alpha=0.00001)



phyloglm_output_centered = phyloglm_wald(y_mat,new_df$abund_scaled,new_tree,n_reps=100,new_df,'logistic_MPLE',start_alpha=start_alpha_val)

median(phyloglm_output_centered$mod_info[phyloglm_output_centered$mod_info$converged,]$alpha)
hist(phyloglm_output_centered$mod_info[phyloglm_output_centered$mod_info$converged,]$alpha)

# pdf('slope_hist-whenphylosignalis0.007.pdf')
median(phyloglm_output_centered$mod_info[phyloglm_output_centered$mod_info$converged,]$b1_est)
hist(phyloglm_output_centered$mod_info[phyloglm_output_centered$mod_info$converged,]$b1_est,
     main=paste0('distribution of slope estimate \nalpha = ',round(moderate_alpha,3)),xlab='b1')
abline(v=beta1,col='red')
# dev.off()

phyloglm_output_centered$analysis_summary


#check to make sure all had at least 100 that converged and that were used in error analysis
sum(phyloglm_output_centered$mod_info$in_error_analysis)

#start error analysis with 99 permutations; start with beta1=0 and sigma=0.5 first
folder_name='error_analysis_bigN_perms_26oct2021'
if(!dir.exists(folder_name)) dir.create(folder_name)

path=paste0(folder_name,'/phyloglm_typeiscaled_b1index_',i,"_alpha_",round(use_me[i,]$alpha,3),"_startalpha_",ifelse(is.null(start_alpha_val),"NULL",start_alpha_val),'.rds')
start_time=proc.time()
permute_output=phyloglm_permute(phyloglm_output_centered,999)
(time_elapsed_minutes=(proc.time()-start_time)[3]/60)
saveRDS(permute_output,path)

# permute_output=readRDS("error_analysis_bigN_perms_26oct2021/phyloglm_typeiscaled_b1index_4_alpha_0.007_startalpha.rds") #b1 of 0.4

"error_analysis_bigN_perms_26oct2021/typei_b1index_1_sigma_0.5.rds"
permute_output[[1]]
permute_output$analysis_summary

mean(permute_output$mod_info[1:100,]$pval<0.05)
phyloglm_output$analysis_summary

mod_info$analysis_summary
permute_output$mod_info

