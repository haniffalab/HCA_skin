#Automatic selection of stats test based on data given parametric/non-para

{ ##Variables module
  #1- Provide an output directory and run name
  filpath ="/Users/issac/Documents/Projects/adult-skin/"
  run_name<-"F13A1"
  
  #2- Provide the input data file as a csv
  my_data <- read.csv("/Users/issac/Documents/Projects/adult-skin/beth_count_stats.csv",stringsAsFactors = F)

  #3- Provide the column ids to identifiers
  # These are columns which you should place equivilent columns in order below
  #"count","group","patient_id","counter_id","image_id"
  set_equiv_column_ids <- c("count","group","patient_id","counter_id","image_id") #Give the col names in your data which match above

  #4- Would you like to remove any patients from your data
  removal_id <- c("Patient 2")
}

{#Run from here

##env setup module##
###################################################################################################################################################
{ # start of setup and functions module
  gc()
  #Define Libraries and communication axis
  libs<-c("dplyr","ggplot2","plyr","BiocManager","multcomp","data.table","ggpubr","gridExtra","rstatix")
  
  #Function to install/load libraries required
  pkg_check <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
      install.packages(new.pkg,repos = "https://cloud.r-project.org", dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
  }
  pkg_check(libs)
  #create and Set working directory
  save_location <- filpath
  dir.create(save_location)
  #setwd(save_location)
  
  #Set columns names appropriately according to lookup
  setnames(my_data, set_equiv_column_ids, col_identifiers)
  
  #Remove any patients decided to remove
  my_data<-my_data[(!my_data$patient_id%in%(c(removal_id))),]
  
}

##general descriptive stats module##
###################################################################################################################################################
{
  #Generate and write descriptive stats
  #Summary stats
  summary_stats<-group_by(my_data, patient_id) %>%
    summarise(
      mean = mean(count, na.rm = TRUE),
      sd = sd(count, na.rm = TRUE),
      median = median(count, na.rm = TRUE),
      IQR = IQR(count, na.rm = TRUE)
    )
  write.csv(summary_stats,paste0(filpath,run_name,"_summary_stats.csv"))
}


##Normaility testing##
###################################################################################################################################################
{
  #Normality tests for each variable
  p <- list()
  normality_rec<-NULL
  for (i in unique(my_data$group)){
    temp <- my_data[(my_data$group==i),]
    shapiro<-shapiro.test(temp$count)
    normality_rec<-append(normality_rec,shapiro$p.value)
    print(i)
    if(shapiro$p.value>0.05){
      test_stat <- paste("distribution is normal p=",shapiro$p.value)
    }else{
      test_stat <- paste("distribution is not normal p=",shapiro$p.value)
    }
    p[[i]] <-ggqqplot(temp$count,xlab = "Theoretical",ylab = "count",title = paste("Oberservational distribution for", i),subtitle=paste(test_stat))
  }
  normnality_test<-do.call(grid.arrange,p)
  jpeg(paste0(filpath,run_name,"_observation_normality_tests",".jpg"),width = 1500, height = 800, units = "px")
  plot(normnality_test)
  dev.off()
}

##Parametric or non_para decision forking##
###################################################################################################################################################
{
  #Check if any values in normality test < 0.05 if so, no normailty can be assumed and non-para should be used
  rec <- normality_rec > 0.05
  if("FALSE"%in%rec){
    print("Non-normal discribution detected in one or more vairables, proceeding to apply non-parametric testing")
    use_para <- FALSE 
  }else{
    use_para <- TRUE
  }
}

##Parametric##
###################################################################################################################################################

if(use_para==TRUE){print("parametric tests..")
  # Compute the analysis of variance
  my_data$group <- as.factor(my_data$group)
  res.aov <- aov(count ~ group, data = my_data)
  # Summary of the analysis
  summary(res.aov)
  #TukeyHSD(res.aov)
  #summary(glht(res.aov, linfct = mcp(patient_id = "Tukey")))
  counter_summary<-summary(glht(res.aov, linfct = mcp(group= "Tukey")))
  multi_comp_dat <- as.data.frame(counter_summary$test$coefficients)
  multi_comp_dat$pval <- counter_summary$test$pvalues

  anova<-head((summary(res.aov)[[1]][4:5]),1)
  rownames(anova)<-"overall"
  colnames(anova) <- colnames(multi_comp_dat)
  #anova<-append("overall",as.character(head(anova,1)))
  multi_comp_dat<-rbind(multi_comp_dat,anova)
  pval_sig <- add_significance(
    multi_comp_dat ,
    p.col = 'pval',
    output.col = NULL,
    cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1),
    symbols = c("****", "***", "**", "*", "ns")
  )
  multi_comp_dat$sig<-pval_sig$pval.signif
  #pairwise.t.test(my_data$count, my_data$group, p.adjust.method = "BH")
  write.csv(multi_comp_dat,paste0(filpath,run_name,"_score_anova.csv"))
  multi_comp_dat
  ###################################
}


##non-Parametric##
###################################################################################################################################################
if(use_para==FALSE){print("non-parametric tests..")
  
  #Multiple wilcox non-parametric pairwise-comparison between groups:
  pairwise <- pairwise.wilcox.test(my_data$count, my_data$group,
                                 p.adjust.method = "BH")
  #overall kruskal wallis
  kruskal<-kruskal.test(my_data$count, my_data$group,
               p.adjust.method = "BH")
  
  
  #melt results
  pairwaise_p<-as.data.frame(pairwise$p.value)
  rownames<-paste0(rownames(pairwaise_p)[row(pairwaise_p)]," - ", colnames(pairwaise_p)[col(pairwaise_p)])
  pairwaise_p<-(melt(pairwaise_p))
  rownames(pairwaise_p)<-rownames
  pairwaise_p<-pairwaise_p[2]
  colnames(pairwaise_p)<-"pval"
  #add overall by kruskal wallis
  overall<-as.data.frame(kruskal$p.value,row.names=c("overall"),col.names=(colnames(pairwaise_p)))
  #append
  pairwaise_p <- rbind(pairwaise_p,anova)
  
  pval_sig <- add_significance(
    pairwaise_p ,
    p.col = 'pval',
    output.col = NULL,
    cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1),
    symbols = c("****", "***", "**", "*", "ns")
  )
  
  pairwaise_p$sig<-pval_sig$pval.signif
  write.csv(pairwaise_p,paste0(filpath,run_name,"_score_anova.csv"))
  pairwaise_p
}

##Plotting module#
###################################################################################################################################################
{
  my_data$group <- as.factor(my_data$group)
  MinMeanSEMMax <- function(x) {
    v <- c(min(x), mean(x) - sd(x)/sqrt(length(x)), mean(x), mean(x) + sd(x)/sqrt(length(x)), max(x))
    names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
    v
  }
  plot <- ggplot(my_data, aes( x = group, y = count, color=group,shape=group)) + 
    scale_shape_manual(values=c(16, 17, 15,18,12)) +
    scale_color_manual(values=c("#00AFBB", "#E7B800", "#FC4E07", "#DC90B3","#2A7946"))+
    #geom_boxplot(inherit.aes = TRUE) +
    stat_summary(fun.data=MinMeanSEMMax, geom="boxplot",width = 0.5) +
    geom_jitter(position=position_jitter(0.2),size = 2) +
    stat_summary(fun=mean, geom="line", aes(group=1),color="black", size=0.5)  + 
    #geom_smooth(aes(group = 1), method="lm", size = 2, se = F) +
    #stat_summary(aes(x=group, y=count), geom="errorbar", fun.data="mean_se",fun.args=list(mult=1), size=0.5, color="black", width=.1) +
    ggtitle('EASI score by visit')+  xlab('Measurement conditional') +  ylab("EASI score")
  plot <- plot + theme_classic(base_size = 20) 
  plot
  jpeg(paste0(filpath,run_name,"_score_by_visit_jitter",".jpg"),width = 800, height = 600, units = "px",pointsize = 20)
  plot(plot)
  dev.off()
}

}#end of loop
