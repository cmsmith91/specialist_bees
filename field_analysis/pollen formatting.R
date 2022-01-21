#pollen formatting
rm(list=ls())
setwd("/Users/colleen/Dropbox/github/inat_project/")
library(tidyverse)
library(readxl)
library(furrr)
library(readxl)

select=dplyr::select;map=purrr::map


sitedate_specplant=read_csv('field_data/sitedate_specialisthostplant-17nov2021.csv')

# name of old 2020 file: 'field_data/pollen_id_SL.csv'
pollen=read_csv('field_data/pollen_id_SL_2020bees-cleaned.csv') %>%
    bind_rows(read_csv('field_data/pollen_id_SL_2021bees-cleaned.csv')) %>%
    rename(unique_id="Special ID",pollen=`Pollen Type`,n=`Number of Pollen Grains`) %>%
    select(unique_id,pollen,n,notes) %>%
    mutate(pollen=tolower(pollen))

#load list of specialist bees in the US
specialist_bees=read_csv("/Users/colleen/Dropbox/Fall_2014/postdoc/inat\ project/data_current/bee_hosts_allUS.csv") %>%
    distinct(bee,host_plant)

#load and format pollen keys for both years
pollen_key=read_excel("field_data/Pollen Reference Guide Specialists 2020.xlsx") %>%
    rename(pollen=`Pollen Name`,species=Genus,family=Family) %>% select(family,pollen,species) %>%
    mutate(pollen=tolower(pollen),epithet=gsub(".* ","",species)) %>%
    filter(!is.na(pollen)) %>%
    mutate(plant_code=tolower(paste0(substr(species,1,3),substr(epithet,1,3))))
pollen_key[pollen_key$species=="Viburnum opulus var. americanum" ,]$plant_code<-'vibtri'
pollen_key[pollen_key$plant_code=='saxvir',]$plant_code <- 'micvir'

pollen_key2021 <- read_excel("field_data/Pollen Reference Guide Specialists Spring 2021.xlsx")%>%
    rename(pollen=`Pollen Name`,species=Genus,family=Family) %>% select(family,pollen,species) %>%
    mutate(pollen=tolower(pollen),epithet=gsub(".* ","",species)) %>%
    mutate(plant_code=ifelse(grepl('Salix',species),'salix',tolower(paste0(substr(species,1,3),substr(epithet,1,3)))))
pollen_key2021[pollen_key2021$plant_code=='saxvir',]$plant_code <- 'micvir'


#since they were called 'rham,'glec and micr in the pollen doc
#remove the other codes from the pollen_key
pollen_key_both=pollen_key %>% bind_rows(pollen_key2021) %>%
    distinct()
pollen_key_both %>% filter(pollen=='saxi')
pollen_key %>% filter(pollen=='saxi')
pollen_key2021 %>% filter(pollen=='saxi')

pollen_key_both %>% filter(pollen=='micr')

pollen_key_both %>% filter(grepl('Micranthes',species))
pollen_key_both=pollen_key_both[!pollen_key_both$pollen %in% c('unstained','ajug','saxi'),]



#load list of focal plants
focal=read_excel('field_data/MethodsTable2021_2022.xlsx') %>% 
    rename_all(tolower) %>%
    mutate(genus=sub(" .*","",species),epithet=sub(".* ","",species)) %>% 
    mutate(plant_code=ifelse(genus=='Salix','salix',paste0(tolower(substr(genus,1,3)), tolower(substr(epithet,1,3))))) %>%
    select(-genus,-epithet)
focal[focal$species=="Saxifraga virginiensis",]$plant_code <- 'micvir'
focal[focal$species=="Viburnum opulus var. americanum",]$plant_code <- 'vibtri'


##
genera=read_csv("field_data/pollen specialization bee id spreadsheet-19jan2022.csv") %>% #load genus level id of bees
    rename(unique_id=uniqueID)%>%mutate(bee=ifelse(!is.na(species),paste0(genus,"_",species),paste0(genus,"_sp")))

#load pinning log from both years
specs=read_csv( "field_data/specialist bee pinning log2020and2021bees-cleaned.csv") %>%
    filter(data)%>% 
    left_join(genera %>% select(unique_id,genus,species))

#get the pollen codes of the focal plants
pollen_key_focal=focal %>% select(species,plant_code) %>% 
    left_join(pollen_key_both %>% select(plant_code,pollen))  


#
specialist_ids=genera  %>% 
    filter(bee %in% specialist_bees$bee )  #get rid of specialist bees from the dataset
unsure_ids=genera %>% filter(bee=='Melissodes_sp.')
data.frame(specialist_ids %>% distinct(bee) %>% left_join(specialist_bees))
male_ids=genera %>% filter(sex=='M')
parasite_ids=genera %>% filter(genus=='Sphecodes')

rm_me=unique(c(specialist_ids$unique_id,male_ids$unique_id,parasite_ids$unique_id,unsure_ids$unique_id))

#make data frame of the plant visited by each bee and the associated pollen code
plants_visited=specs %>% filter(!unique_id %in% rm_me) %>%
    select(unique_id,plant_code) %>%
    left_join(pollen_key_focal %>% select(-species))%>%
    rename(host_plant=plant_code,pollen_of_hostplant=pollen)# %>%
nrow(plants_visited); nrow(specs %>% filter(!unique_id %in% rm_me) )


#now filter pollen data to just be data bees, and join with the plants visited df
specs_pollen=pollen %>% filter(unique_id %in% plants_visited$unique_id)  %>%
    left_join(plants_visited )

#get the proportion of each pollen type in each pollen load
pollen_summ=specs_pollen %>%group_by(unique_id,pollen,host_plant,pollen_of_hostplant) %>%
    summarize(n=sum(n)) %>%
    left_join(specs_pollen %>% group_by(unique_id) %>% summarize(total_n=sum(n))) %>%
    mutate(prop=n/total_n,prop200=n/200) %>%
    mutate(host_pollen=pollen==pollen_of_hostplant)

# get the proportion of the pollen of the plant visited pollen in each pollen load
prop_host=1:length(unique(pollen_summ$unique_id)) %>% map_dfr(function(i){
    the_bee=unique(pollen_summ$unique_id)[i]
    the_df=pollen_summ %>% filter(unique_id==the_bee)
    
    if(sum(the_df$host_pollen) > 0) {
        
        new_df = data.frame(the_df %>% filter(host_pollen)) %>%
            select(unique_id,host_plant,prop200,n=n) %>% 
            rename(prop=prop200)
        
    }else{
        new_df1=data.frame(the_df) %>% select(unique_id,host_plant) %>%
            mutate(prop=0,n=0)
        new_df=new_df1[1,]
    }
    return(new_df)
    
})  

# next at site and specialist host plant data to prop_host
data=prop_host %>%
    left_join(specs %>% select(unique_id,site,date)) %>% 
    left_join(sitedate_specplant)
head(data)
##get summary stats (first with Cucurbita pepo sites, then without)
paste0("We collected ", nrow(data), " individual bees total from generalist species")

#get rid of cucurbita pepo sites bc the specialist wasn't
# even collecting pollen (see below)
#at one of the cucpep sites, we were also sampling solcan - don't get rid of that site,
#just get rid of bees sampled from cucurbita at that site
head(data)
data=data %>% filter(spec_host!='cucpep' & host_plant !='cucpep')
#double check that mil is still there, but that I removed 
data %>% filter(grepl('cucpep',site_spec)) %>% distinct(host_plant) #l

#get summary stats without the cucpep sites
paste0("The sample size of generalist bees after this filtering step was " ,nrow(data),' individual bees')

# write_csv(data,'field_data/proportion_pollen_data_19jan2022.csv')

##################Specialist bee data formatting below
#finally get the proportion of host plant pollen found in the pollen loads of specialist bees
#make data frame of the plant visited by each bee and the associated pollen code
plants_visited=specs %>% filter(unique_id %in% specialist_ids$unique_id) %>%
    select(unique_id,plant_code) %>%
    left_join(pollen_key_focal %>% select(-species))%>%
    rename(host_plant=plant_code,pollen_of_hostplant=pollen)# %>%
nrow(plants_visited); nrow(specs %>% filter(unique_id %in% specialist_ids$unique_id) )

specialist_specs=specs %>% filter(unique_id %in% specialist_ids$unique_id) %>%
    mutate(genus_species=paste0(genus,"_",species))
n_distinct(specialist_specs$genus_species)
nrow(specialist_specs)

#now filter pollen data to just be data bees, and join with the plants visited df
specs_pollen= pollen %>% 
    filter(unique_id %in% plants_visited$unique_id)  %>%
    left_join(plants_visited )

#get the proportion of each pollen type in each pollen load
pollen_summ=specs_pollen %>%group_by(unique_id,pollen,host_plant,pollen_of_hostplant) %>%
    summarize(n=sum(n)) %>%
    left_join(specs_pollen %>% group_by(unique_id) %>% summarize(total_n=sum(n))) %>%
    mutate(prop=n/total_n,prop200=n/200) %>%
    mutate(host_pollen=pollen==pollen_of_hostplant)

# get the proportion of the pollen of the plant visited pollen in each pollen load
prop_host=1:length(unique(pollen_summ$unique_id)) %>% map_dfr(function(i){
    the_bee=unique(pollen_summ$unique_id)[i]
    the_df=pollen_summ %>% filter(unique_id==the_bee)
    
    if(sum(the_df$host_pollen) > 0) {
        
        new_df = data.frame(the_df %>% filter(host_pollen)) %>%
            select(unique_id,host_plant,prop200,n=n) %>% 
            rename(prop=prop200)
        
    }else{
        new_df1=data.frame(the_df) %>% select(unique_id,host_plant) %>%
            mutate(prop=0,n=0)
        new_df=new_df1[1,]
    }
    return(new_df)
    
})  

# next add site and specialist host plant data to prop_host
data_specs=prop_host %>%
    left_join(specs %>% select(unique_id,site,date)) %>% 
    left_join(sitedate_specplant) %>%
    left_join(genera %>% select(unique_id,genus,species)) %>%
    mutate(genus_species=paste0(genus,'_',species))
hist(data_specs$prop)
spec_hosts=c('clacar','salix','vacmyr','cucpep','corser','solcan','lyscil','poncor')
data_specs %>% filter(!host_plant %in% spec_hosts) %>%
    select(unique_id,genus_species,host_plant,spec_host,prop)
data_specs %>% filter(host_plant =="salix") %>%
    select(genus_species,host_plant,spec_host,prop,n)

par(mfrow=c(2,5))
for(the_plant in unique(data_specs$host_plant)){
    with(data_specs %>% filter(host_plant==the_plant),hist(prop,main=the_plant))

    } #how is the prop for anecan negative?? 
#also, need to change so that it's only specialists on their host plant
#use specialist_bees
data_specs2=data_specs %>% left_join(specialist_bees %>% rename(genus_species=bee,spec_host_plant=host_plant))
unique(data_specs2$spec_host_plant)

#modify specialist_bees df so there's one host per bee
#go with the lowest tax resolution (so Solidago if the bee specializes on solidago vs also asteraceae)
our_specialist_bees=specialist_bees %>% filter(bee %in% data_specs$genus_species) %>% split(.$bee) %>%
    map_dfr(function(df){
        if("Solidago" %in% df$host_plant){
            concise_host="Solidago"
        }else{
            if("Asteraceae" %in% df$host_plant){
                concise_host='Asteraceae'
            }else{
                concise_host=df$host_plant[1]
            }
            
        }
        data.frame(bee=df$bee[1],host=concise_host)
    })
#associate the host plant with it's pollen /  plant code
hostplant_toplantcode=data.frame(host=unique(our_specialist_bees$host),
           host_plant_code=c("vacmyr","clacar","salix","corser","solcan","poncor","lyscil","solcan","cucpep"))
code_to_species=data.frame(ho)
the_plant=unique(hostplant_toplantcode$host_plant_code)[5]
pdf('figures/specialist_foraging.pdf',width=11)
par(mfrow=c(2,4),cex.lab=1.5,mar=c(5,4.5,4,2),cex.main=1.5,cex.axis=1.3)
for(the_plant in unique(hostplant_toplantcode$host_plant_code)){
    plant_genus=hostplant_toplantcode[hostplant_toplantcode$host_plant_code==the_plant,]$host[1]
    if(plant_genus=='Ericaceae') plant_genus='Vaccinium'
    host_df=hostplant_toplantcode %>% filter(host_plant_code==the_plant)
    spec_bees=our_specialist_bees %>% filter(host==host_df$host[1])
    plot_df=data_specs %>% filter(host_plant==the_plant & genus_species %in% spec_bees$bee)
    
    n=nrow(plot_df)
    the_main=paste0(plant_genus,"\n",'n = ',n)
    with(plot_df,hist(n,main=the_main,xlab='number of pollen grains'))
    
} 
dev.off()

#for cucurbita check if there are any differences by site
the_plant='cucpep'
host_df=hostplant_toplantcode %>% filter(host_plant_code==the_plant)
spec_bees=our_specialist_bees %>% filter(host==host_df$host[1])
plot_df=data_specs %>% filter(host_plant==the_plant & genus_species %in% spec_bees$bee)
plot_df_ls=plot_df %>% split(.$site_spec)

par(mfrow=c(1,3))
for(new_plot_df in plot_df_ls){
    n=nrow(new_plot_df)
    
    the_main=paste0(new_plot_df$site[1],"\n",'n = ',n)
    with(new_plot_df,hist(n,main=the_main,xlab='number of pollen grains'))
    
    
}

#get the median number carried at each site
plot_df_ls %>% map(function(df) median(df$n))

#specialist_specs and specialist_ids are different sizes. why?
specialist_ids %>% filter(!unique_id %in% specialist_specs$unique_id)
all_specs=read_csv( "field_data/specialist bee pinning log2020and2021bees-cleaned.csv")
all_specs %>% filter(unique_id==1164) #ahh collected off of an achmil, a non-data plant
pollen %>% filter(unique_id==1164)
