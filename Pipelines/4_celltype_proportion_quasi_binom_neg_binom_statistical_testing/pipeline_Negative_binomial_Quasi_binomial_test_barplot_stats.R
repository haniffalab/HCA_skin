#Author: Issac Goh
#Date: 010620
#Modified: 010620

###############################################################################
###############################################################################
#Tips

# x is a dataframe with the following columns in this specific order:
# x$variable # celltype
# x$value # cell count
# x$batch # condition_donor
# x$timepoint # condition (e.g. adult, fetal)"

#If you are importing a metafata dataframe from another object, 
#use below to import and subset the data to columns of interest

###############################################################################
###############################################################################
#Data Preperation module#

# - Give adresses of metadata to concat
adrr<-("/Users/issac/Documents/Projects/Fetal Bone Marrow/BM_12052020_stats/meta_data_files/")

# - Concat data of interest
concat_adrr_book<-function(){
adress_book<- grep(list.files(path=adrr), pattern='concatenated_metadata.csv', inv=T, value=T)
if ("met_concat"%in%ls()){
  rm(met_concat)}
for (met_adr_temp in adress_book){
  met_temp <- read.csv(paste0(adrr,met_adr_temp))
  if (!"met_concat"%in%ls()){
    met_concat <- met_temp
  }else{
    met_concat <- plyr::rbind.fill(met_concat,met_temp)
  }
}
return(met_concat)
}
met_concat<-concat_adrr_book()
write.csv(met_concat,paste0(adrr,"concatenated_metadata.csv"))


# - Read cell_list for each fig and define which is current figure
fig_cell_list_key <- read.csv("/Users/issac/Documents/Projects/Fetal Bone Marrow/BM_12052020_stats/cell.labels by figures.csv",stringsAsFactors = F)
fig_of_interest<-"Figure.S6"
cells_list<- as.character(unlist(fig_cell_list_key[fig_of_interest][!fig_cell_list_key[fig_of_interest]==""]))

###############################################################################
###############################################################################
#Variables input module#

# - choice of Negative binom of quasi binomial distribution ('quasi' or 'neg')
dist <- "quasi" 

# - metadata adress
meta_add<-"/Users/issac/Documents/Projects/Fetal Bone Marrow/BM_12052020_stats/fBM_overall_fig1/fig1d_barplot_meta_20200513.csv"

# - output directory to create
output_dr<-"/Users/issac/Documents/Projects/Fetal Bone Marrow/BM_12052020_stats/pcw_stroma_figs6"

# - input dependent cateogircal variable containing cell.label information
dep_var<-"cell.labels"

# - input independent cateogircal variable containing the test conditions (e.g stage/pcw/sex)
inde_var<-"stage"

# - input categorical variable containing batch information (a column to compute dispersion coefficients)
batch<-"orig.ident"

# - run name
run_id <- "pcw_stroma_fbm"

# - give list of cell labels in verbatim to metadata to subset the data by (note that the cell.labels of analougous populations must match across data)
cells_keep<-cells_list #(give a list of cell.labels)

# - if input data needs to be subset give list of indpendent categroricals to be kept (e.g stage/pcw/sex or "all")
inde_keep<-c("all")

###############################################################################
###############################################################################
####Env setup and functions module####
{
#Function to install/load libraries required
libs<-c("base","dplyr","data.table","plyr","MASS","stringr","gdata","reshape2")
pkg_check <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg,repos = "https://cloud.r-project.org", dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
pkg_check(libs)

#create and Set working directory
dir.create(output_dr)
setwd(output_dr)

#Read the metadata and subset to independent and dependent categoricals of interest
if (!"meta"%in%ls()){
meta<-read.csv(meta_add,stringsAsFactors = F)
}
meta<- meta[,c(batch,dep_var,inde_var)]
cells_keep <- toupper(cells_keep)
meta[,dep_var]<-toupper(meta[,dep_var])

if(inde_keep=="all"){
  print("option to use all independent variables chosen, automatically assigning comparison for everything")
}else{
  print(paste0("option to subset of independent variables chosen, automatically assigning comparison for:",inde_keep))
  meta<-meta[(meta[,inde_var]%in%inde_keep),]
}
meta<-meta[(meta[,dep_var]%in%cells_keep),]

#Check if unique labels count is equal to cell_keep
if(!length(unique(meta[,dep_var])) == length(cells_keep)){
  warning("length of cell_list and unique labels is not equal, check spelling")
}

#Assign a freq count value per cell
meta$freq = 1

#Aggregate the metadata in by grouping variable ()
meta[,batch][is.na(meta[,batch])] <- "not_specified"
x <- aggregate(meta[,4], meta[,1:3], FUN = sum)
#convert all labels to uppercase (minimises errors)
x <- as.data.frame(x%>%apply(2,toupper))
#relabel and reorder x
colnames(x)<- c("batch","variable","timepoint","value")
col_order <- c("value","timepoint","variable","batch")
x <- x[, col_order]
x$value<-as.integer(x$value)

#Convert data to proportional data if quasibinom option chosen
if(dist=="quasi"){
  x <- x %>% dplyr::group_by(timepoint) %>% dplyr::mutate(percent = value/sum(value))
}

# Define all your conditions that group your replicates as the "batch variable" (e.g donors), this goes into df of conditional "variates"
len = length(unique(x$batch))
condition = as.data.frame(1: len)
colnames(condition)<-"len"
condition$len <-gsub("^","cond_",condition$len)
condition$cond <-unique(x$batch)

# Define the variates in "condition$cond" across dataset
conditions<-as.character(unique(x$timepoint))

# Define as many empty variables as many conditions you want to test:
n<-length(conditions)
test_vecs<-as.data.frame(gtools::combinations(n,2,conditions))

#Create a unpopulated df
score_results<- data.frame(row.names = paste0(test_vecs$V1,"_vs_",test_vecs$V2))
}
###############################################################################
###############################################################################
####compute module####
{
# Compute the p-values for increasing labeled fraction (using negative binomial regression)
if("score_results_concat"%in%ls()){
  remove(score_results_concat)
}
pv_tot <-vector()
for(i in 1:length(unique(x$variable))) {
  #pv_temp<-vector()
  pv_temp<-rep("NA",NROW(score_results))
  temp_test_list<-vector()
  cell = levels(x$variable)[i]
  totals = aggregate(x$value, by = list(x$batch), FUN = sum)
  colnames(totals) <- c("batch", "total")
  d = data.frame(subset(x, variable == cell))
  d <- merge(d, totals, by = "batch")
#QC block
  if(any(d$total == 0)) {
    warning("Removing sample that didn't detect this celltype")
    d = subset(d, total > 0)
  }
  if(length(unique(d$timepoint))<=1){
    warning("One or fewer unique independent variables detected in loop, proceeding with next loop")
    }else{
    for (ref in 1:NROW(test_vecs)){
      temp_vec<-test_vecs[ref,]
      reference <- as.character(unlist(temp_vec[1]))
      test <-as.character(unlist(temp_vec[2]))
      test_name<-paste0(reference,"_v_",test)
      temp_test_list[ref] = test_name
      #test block
      if(reference%in%d$timepoint&&test%in%d$timepoint){
        x1 = d[d$timepoint %in% c(reference, test), ]
        if(dist=="quasi"){
          nb = glm(formula = percent ~ timepoint, data = x1, family=quasibinomial)
        }else{
          nb = MASS::glm.nb(formula = value ~ timepoint + offset(log(as.numeric(x1$total))), data = x1, maxit = 1000)#, control=glm.control(trace = 3))
        }
        tryCatch(anova(nb, test = "LRT")$`Pr(>Chi)`[2],error=function(e){pv_temp[ref]<<-"NA"})
        pv_temp[ref] = anova(nb, test = "LRT")$`Pr(>Chi)`[2]
       score_results[,cell] <- pv_temp
      }
    }
      #Compute overall score
      if(dist=="quasi"){
        nb_tot = glm(formula = percent ~ timepoint, data = d, family=quasibinomial)
      }else{
        nb_tot = MASS::glm.nb(formula = value ~ timepoint + offset(log(as.numeric(d$total))), data = d, maxit = 1000)#, control=glm.control(trace = 3))
      }
      pav_tot_val <- anova(nb_tot, test = "LRT")$`Pr(>Chi)`[2]
      if(is.na(pav_tot_val)==FALSE){
      pv_tot[i] = pav_tot_val
      }else{
      pv_tot[i] <- "NA"
      }
      
    }
  #if(!"score_results_concat"%in%ls()){
  #  score_results_concat<-score_results
  #}else{
  #  score_results_concat<-cbind(score_results_concat,score_results)
  #}
}
score_results_final<-rbind(score_results,pv_tot)
base::row.names(score_results_final[NROW(score_results_final),])<-"overall"
}

###############################################################################
###############################################################################
####results module####
date<-gsub( " .*$", "", Sys.time())
write.csv(score_results_final,paste0(output_dr,"/",run_id,"_",date,".csv"))


###############################################################################
###############################################################################
####module to add:: polarity####

#Decide if prop goes up or down (use non-parametric trend test)

direction = NULL
cell_list = NULL
rho = NULL
for(i in 1:length(unique(x$variable))) {
  cell = levels(x$variable)[i]
  d = data.frame(subset(x, variable == cell))
  totals = aggregate(d$value, by = list(d$timepoint), FUN = sum)
  colnames(totals) <- c("timepoint", "total")
  d <- merge(d, totals, by = "timepoint")
  
  #Create order for factors
  order(unique(d$timepoint))
  d$factor_levels <- ordered(d$timepoint, levels = c("STAGE 1", "STAGE 2", "STAGE 3", "STAGE 4"))
  d1<-d[!duplicated(d$timepoint),]
  d1$time_factor <- order(d1$factor_levels)
  
  if(length(unique(d$timepoint))>=2){
    print("parameters met to use Bartels trend test statistic")
    #d1$time_factor<-order(d1$timepoint)
    ts_d<-as.ts(d1)
    #sig<-bartels.test(ts_d)
    sig <- cor.test(d1$value,d1$time_factor,method = "spearman")
    cell_list <- append(cell_list,cell)
    rho<- append(rho,sig$estimate)
    if(as.numeric(sig$estimate) == '1' ){
      direction <- append(direction,"No_direction")
    }else if(sig$estimate > '0' ){
      direction<-append(direction,"up")
    }else{
      direction<-append(direction,"down")
    }
    
  }else{
    #If no statistic can be applied, use simple arithmic and compare largest to smallest factor
    reference <-as.character(d1$timepoint[d1$time_factor== min(d1$time_factor)]) # !
    test <- as.character(d1$timepoint[d1$time_factor== max(d1$time_factor)]) # !
    x1 = d[d$time %in% c(reference, test), ]
    x_sums_test <- sum(x1$value[x1$timepoint==test])
    x_sums_ref<- sum(x1$value[x1$timepoint==reference])
    
    cell_list <- append(cell_list,cell)
    rho<- append(rho,"NA")
    if(x_sums_test > x_sums_ref){
      direction <- append(direction,"up")
    }else{
      direction<-append(direction,"down")
    }
  }
  

}
polarity<-data.frame(cell_list,direction,rho)
write.csv(polarity,"/Users/issac/Documents/Projects/Fetal Bone Marrow/BM_12052020_stats/pcw_stroma_figs6/fig_S6_polarity.csv")
#pv$polarity_of_change <- direction
