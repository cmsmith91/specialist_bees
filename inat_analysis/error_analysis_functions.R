#load packages
library(tidyverse)
library(phylolm)
library(furrr) 
map=purrr::map; select=dplyr::select

phyloglm_wald=function(y_mat,xs,tree,n_reps,my_data,the_method, start_alpha='null'){
    plan(multisession,workers=6)
    mod_info=1:ncol(y_mat) %>% future_map(function(j){
        ys=y_mat[,j]
        if(start_alpha != 'null') mod=phyloglm(ys~xs,phy=tree,data=my_data,btol=40,method = c(the_method),start.alpha = start_alpha)
        if(start_alpha == 'null') mod=phyloglm(ys~xs,phy=tree,data=my_data,btol=40,method = c(the_method))
        

        tibble(
            b0_est=mod$coefficients[1],
            b1_est=mod$coefficients[2],
            pval=coef(summary(mod))[2,4],
            alpha=mod$alpha,
            alpha_warning=mod$alphaWarn,
            converged=mod$convergence==0,
            mod_index=j
            
        )
    }) %>% bind_rows
    future:::ClusterRegistry("stop")
    
    #only estimate error rates with first of n_reps models that converged
    test_me=mod_info$mod_index[mod_info$converged][1:n_reps]
    mod_info$in_error_analysis=mod_info$mod_index %in% test_me
    
    if(sum(mod_info$converged)<n_reps) warning('there were less than n_reps models that converged')
    
    mod_summary=with(mod_info %>% filter(in_error_analysis),data.frame(
        mean_b0=mean(b0_est),
        sd_b0=sd(b0_est),
        mean_b1=mean(b1_est),
        sd_b1=sd(b1_est),
        prop_sig=mean(pval<.05),
        method='phyloglm_wald',
        n=length(test_me) ,
        start_alpha_val=start_alpha,
        alg=the_method
    ))
    output=list(mod_info=mod_info,analysis_summary=mod_summary,ys=y_mat,xs=xs,
                tree=tree,n_reps=n_reps,data=my_data)
    return(output)
}


phyloglm_permute=function(wald_object,n_perm){
    
    #only do permutations for models used in wald error analysis 
    mod_info=wald_object$mod_info %>% filter(in_error_analysis)
    y_mat=wald_object$ys
    xs=wald_object$xs
    n_reps=wald_object$n_reps
    my_tree=wald_object$tree
    my_data=wald_object$data
    start_alpha=wald_object$analysis_summary$start_alpha_val
    the_method=wald_object$analysis_summary$alg
    
    

    plan(multisession,workers=6)
    #do more permutation tests than you need bc some models won't converge
    permut_detailed=mod_info$mod_index %>% purrr::map(function(my_mod_index){
        
        #get ys and
        ys=y_mat[,my_mod_index]
        
        #get the slope estimate from the phylglm
        obs_beta=mod_info[mod_info$mod_index==my_mod_index,]$b1_est
        
        #for each of n_try replicates, permute data and refit model
        # save null beta, and whether the model converged
        
        permut_output=1:n_perm %>% future_map(possibly(function(i){
            new_y=sample(ys)
            if(start_alpha != 'null') mod=phyloglm(new_y~xs,phy=my_tree,data=my_data,btol=40,method = c(the_method),start.alpha = start_alpha)
            if(start_alpha == 'null') mod=phyloglm(new_y~xs,phy=my_tree,data=my_data,btol=40,method = c(the_method))
            
            null_beta=coef(mod)[2]
            converged=mod$convergence==0
            # 
            tibble(null_beta=null_beta,converged=converged,mod_index=my_mod_index,obs_beta=obs_beta)
        },otherwise=NULL),.options=furrr_options(seed=T)) %>% bind_rows
        
        # check whether number of permutations =n_perm
        # if not, do more permtuations in units of 10
        # until you get the correct number of converged permutations
        while(sum(permut_output$converged)<n_perm){
            add_me=1:10 %>% future_map(possibly(function(i){
                new_y=sample(ys)
                if(start_alpha != 'null') mod=phyloglm(new_y~xs,phy=my_tree,data=my_data,btol=40,method = c(the_method),start.alpha = start_alpha)
                if(start_alpha == 'null') mod=phyloglm(new_y~xs,phy=my_tree,data=my_data,btol=40,method = c(the_method))
                
                null_beta=coef(mod)[2]
                converged=mod$convergence==0
                #
                tibble(null_beta=null_beta,converged=converged,mod_index=my_mod_index,obs_beta=obs_beta)
            },otherwise=NULL),.options=furrr_options(seed=T)) %>% bind_rows
            permut_output=permut_output %>% bind_rows(add_me)

        }
        
        
        # ##to get p-value compare observed beta to null beta
        # used the first n_perm models that converged
        permut_output$unique_perm=1:nrow(permut_output) # give each permutation a unique id
        test_me=permut_output$unique_perm[permut_output$converged][1:n_perm]
        permut_output$used_in_permanalysis=permut_output$unique_perm %in% test_me
        
        
        return(permut_output)
        
    })
    permut_pvals=permut_detailed %>% map_dfr(function(permut_df){
        null_betas=permut_df[permut_df$used_in_permanalysis,]$null_beta
        obs_beta=permut_df$obs_beta[1]
        r=sum(null_betas>obs_beta)
        pval=(r+1)/(n_perm+1)
        
        permut_result=tibble(
            pval=pval,n_used_for_p=length(null_betas),
            mod_index=permut_df$mod_index[1]
        )
    })
    
    permut_summary=tibble(
        prop_sig=mean(permut_pvals[1:n_reps,]$pval<.05),
        method='phyloglm_permute',
        n=nrow(permut_pvals)
        
    )
    future:::ClusterRegistry("stop")
    output=list(mod_info=permut_pvals,analysis_summary=permut_summary,permut_list=permut_detailed)
    return(output)

    }
