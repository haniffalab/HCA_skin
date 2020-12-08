#Author: Issac Goh
#Written on: 050620
#last modified: 150620

###################################################################################################################################################
{ ##Variables module
#1-  Input ranked gene list ids
gene_list_ids <- c("downs_fbm_conserved",
                   "DE_downs",
                   "DE_fbm")

#2-  Input ranked gene list adresses
gene_list_adresses <- c("/Users/issac/Documents/Projects/Fetal Bone Marrow/network_gsea/HSC_downs_fbm/conserved_hsc_markers_fbm_downs.csv",
                        "/Users/issac/Documents/Projects/Fetal Bone Marrow/network_gsea/HSC_downs_fbm/DEGS_downs_hsc.csv",
                        "/Users/issac/Documents/Projects/Fetal Bone Marrow/network_gsea/HSC_downs_fbm/DEGS_fbm_hsc.csv") 

#3- Input ranked gene column id
gene_id <- c("X",
             "genes",
             "genes")

#- Provide run name
run_id <- "downs_Fbm_comp"

#- Provide number of ranked genes for each set
ngene <- 1000

#- Provide minimum intersect size for Geneset
intersect = 5

# - Provide a location you wish to save your files to
save_location = "/Users/issac/Documents/Projects/Fetal Bone Marrow/network_gsea/HSC_downs_fbm/"
}

##env setup module##
###################################################################################################################################################
{ # start of setup and functions module
  gc()
  #Define Libraries and communication axis
  libs<-c("RCy3","zoo","dplyr","data.table","plyr","BiocManager","gProfileR","gprofiler2","textnets","MCL","RColorBrewer")
  cyto_app_toinstall <- c("clustermaker2", "enrichmentmap", "autoannotate", "wordcloud", "stringapp", "aMatReader")
  
#Function to install/load libraries required
  pkg_check <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
   if (length(new.pkg)) 
     install.packages(new.pkg,repos = "https://cloud.r-project.org", dependencies = TRUE)
   sapply(pkg, require, character.only = TRUE)
}
  pkg_check(libs)
#create and Set working directory
  dir.create(save_location)
  setwd(save_location)
#Connect to Cytoscape instance
  cytoscapePing ()
  cytoscapeVersionInfo ()

#Get list of cytoscape dependencies and install if required
  installation_responses <- cyto_app_toinstall
  cytoscape_version <- unlist(strsplit( cytoscapeVersionInfo()['cytoscapeVersion'],split = "\\."))
  if(length(cytoscape_version) == 3 && as.numeric(cytoscape_version[1]>=3) 
     && as.numeric(cytoscape_version[2]>=7)){
   for(i in 1:length(cyto_app_toinstall)){
      #check to see if the app is installed.  Only install it if it hasn't been installed
     if(!grep(commandsGET(paste("apps status app=\"", cyto_app_toinstall[i],"\"", sep="")), 
              pattern = "status: Installed")){
        installation_response <-commandsGET(paste("apps install app=\"", 
                                                 cyto_app_toinstall[i],"\"", sep=""))
       installation_responses <- c(installation_responses,installation_response)
      } else{
       installation_responses <- c(installation_responses,"already installed")
      }
   }
    installation_summary <- data.frame(name = cyto_app_toinstall, 
                                       status = installation_responses)
  
   knitr::kable(list(installation_summary),
                booktabs = TRUE, caption = 'A Summary of automated app installation'
   )
  }

# Troubleshooting cytoscape commands
#help(package=RCy3)
#RCy3::commandsHelp("help string")
#RCy3::commandsHelp("string protein query")
##DEGs reading and modification module##
runGprofiler2 <- function(genes,current_organism = "hsapiens",filter_gs_size_min = intersect)
{
  gprofiler_results <- gost(genes,organism = current_organism, ordered_query = TRUE,significant = TRUE,exclude_iea = F,user_threshold = 0.05,correction_method ="false_discovery_rate",sources = c("GO:BP"),evcodes = TRUE)
  #filter gprofiler results
  gprofiler_meta = gprofiler_results$meta
  gprofiler_results = gprofiler_results$result
  gprofiler_results <- gprofiler_results[which(gprofiler_results[,'term_size'] >= 3 & gprofiler_results[,'term_size'] <= 200 & gprofiler_results[,'intersection_size'] >= filter_gs_size_min ),]
  # gProfiler returns corrected p-values only as such, use pvalue as FDR
  if(dim(gprofiler_results)[1] > 0){
    gprofiler_results[, c("term_id","term_name")]
    em_results <- cbind(gprofiler_results[, c("term_id","term_name","term_size","p_value","intersection_size","parents","evidence_codes","intersection")])
    colnames(em_results) <- c("Name","Description","intersection_size", "pvalue","term_size","parents","evidence_codes","intersection")
    return(em_results)
  } else {
    return("no gprofiler results for supplied query")
  }
}

# Quick Check if files exist
check_files<-function(){
  for (exist in (ad_dic$gene_list_adresses)){
    if (file.exists(exist)==FALSE){
      warning(paste("File corresponding to",i,"does not exist!"))
      break
    }
  }
  print(paste("All files located, proceeding with network generation"))
}

# function to sumit enrichment query
submit_enrich_query<-function(network_name=i){
#remove non-significant results from gprofiler results, this should already be removed if "option of significant = tue was already selected"
  em_results_filtered<-em_results[(em_results$pvalue<0.05),]
#rank em results by pvalue of ranking
  em_results_filtered<-em_results[order(em_results$pvalue,decreasing=FALSE),]
## write out the g:Profiler results
  em_results_filtered$parents<-as.character(em_results_filtered$parents)
  em_results_filtered$genes<-as.character(em_results_filtered$intersection)
  em_results_filename <-file.path(getwd(),"enr_results.txt")
  write.table(em_results_filtered,em_results_filename,col.name=TRUE,sep="\t",row.names=FALSE,quote=FALSE)
#output for enr map
  em_results_filtered$qvalue<-em_results_filtered$pvalue
  em_results_filtered$phenotype<-1
  em_results_enrichment<- cbind(em_results_filtered[, c("Name","Description", "pvalue","qvalue","phenotype","genes","list_id","run_id")])
  em_results_enrichment_filename <-file.path(getwd(),"enr_results_for_enrichment_map.txt")
  write.table(em_results_enrichment,em_results_enrichment_filename,col.name=TRUE,sep="\t",row.names=FALSE,quote=FALSE)
# Submit EnrichmentMap query
  em_command = paste('enrichmentmap build analysisType="generic" ', 'pvalue=',"0.05", 'qvalue=',"0.05",'similaritycutoff=',"0.25",'coeffecients=',"JACCARD",'enrichmentsDataset1=',em_results_enrichment_filename ,sep=" ")
#enrichment map command will return the suid of newly created network.
  em_network_suid <- commandsRun(em_command)
#get the column from the nodetable and node table
  nodetable_colnames <- getTableColumnNames(table="node",  network =  as.numeric(em_network_suid))
  descr_attrib <- nodetable_colnames[grep(nodetable_colnames, pattern = "GS_DESCR")]
#Rename the network according to current loop
  current_network_id = getNetworkName()
  renameNetwork(title = network_name, network=current_network_id)
}

# Run Auto_annotate with clustering for individual networks
auto_annotate<-function(cluster_boosted,cluster_id = "__mclCluster"){
  current_network_id = getNetworkName()
  nodetable_colnames<-getTableColumnNames(table="node",  network =  current_network_id)
  descr_attrib <- nodetable_colnames[grep(nodetable_colnames, pattern = "GS_DESCR")]
  clus_id <- cluster_id
  if(cluster_boosted == TRUE){
    autoannotate_url <- paste("autoannotate annotate-clusterBoosted", " labelColumn=", descr_attrib," maxWords=4 ", sep="")
    current_name <-commandsGET(autoannotate_url)
#Create a unique network to cluster identifier column
    nodetable <- as.data.frame(getTableColumns(table= "node", network = current_network_id))
#load new columns into cytoscape
    deleteTableColumn(column = cluster_id )
    nodetable[,cluster_id] <- as.integer(nodetable[,cluster_id]+1)
    nodetable$run_id <- run_id
    nodetable$net_id <- current_network_id
    nodetable$clust_id <- paste0(as.character(nodetable[,cluster_id]),"_",current_network_id)
    loadTableData(nodetable,data.key.column="name")
  }else{
    autoannotate_url <- paste("autoannotate annotate-clusterBoosted clusterIdColumn=",clus_id , " labelColumn=", descr_attrib," maxWords=4 ","useClusterMaker=false ", sep="")
    current_name <-commandsGET(autoannotate_url)
  }
# Re-layout the network
  getLayoutNames()
  getLayoutPropertyNames(layout.name='attributes-layout')
  getTableColumnNames(table="node")
  layout.url  = paste0("layout attributes-layout nodeAttribute=",clus_id)
  commandsGET(layout.url)
# Set style name
  style.name = "new_fig_style"
  setVisualStyle(style.name=style.name)
}

# Create summary network and output both full nodetable and summary table
create_summary_network <-function(network_id=i){
  autoannotate_url <- paste("autoannotate summary includeUnclustered='false network=",network_id)
  current_name <-commandsGET(autoannotate_url)
  summary_nodes <- getTableColumns(table="node",columns=c("name"))
  #clearNodePropertyBypass(node.names = summary_nodes$name,visual.property = "NODE_SIZE")
  style.name = "new_fig_style_summary"
  setVisualStyle(style.name=style.name)
  renameNetwork(title = paste0(network_id,"_summary"), getNetworkName())
  
  #Get the pre-summary node table
  network_id<-getNetworkName()
  nodetable <- as.data.frame(getTableColumns(table= "node" , network = network_id))
  nodetable$`EnrichmentMap::Genes`<-as.character(nodetable$`EnrichmentMap::Genes`)
  #nodetable$EnrichmentMap_Genes<-as.character(nodetable$EnrichmentMap_Genes)
  
  # output the node and edge table of summary network publishing
  current_network_id = getNetworkName()
  nodetable_summary <- as.data.frame(getTableColumns(table= "node" , network = current_network_id))
  nodetable_summary$selected = NULL
  nodetable_summary$`EnrichmentMap::Genes`<-as.character(nodetable_summary$`EnrichmentMap::Genes`)
  #nodetable_summary$EnrichmentMap_Genes<-as.character(nodetable_summary$EnrichmentMap_Genes)
  
  specific_save_local<-paste0(save_location,"/",network_id,"/")
  dir.create(specific_save_local)
  node_file_name = paste0(specific_save_local,network_id,"all_nodes.csv")
  write.csv(nodetable,node_file_name)
  node_file_name = paste0(specific_save_local,network_id,"_summary_node.csv")
  write.csv(nodetable_summary,node_file_name)
}

#Function to merge all networks in disctionary
merge_all_in_ad_dic<-function(){
  print("commencing union of networks selected")
  union_list <-ad_dic$gene_list_ids[1:2]
  if(length(ad_dic$gene_list_ids)>2){
    mergeNetworks(union_list,title = "temp_merged",nodesOnly = TRUE)
    for (z in (ad_dic$gene_list_ids[!ad_dic$gene_list_ids%in%union_list])){
      mergeNetworks(c(getNetworkName(),z),title = "temp_merged")
    }
    current_network_id = getNetworkName()
    renameNetwork(title = paste0(run_id,"_merged_networks"), network=current_network_id)
    rm_list<-(grep(pattern = "temp",value = TRUE,getNetworkList()))
    for (o in rm_list){
      deleteNetwork(network = o)
    }
  }else{
    mergeNetworks(union_list,title = paste0(run_id,"_merged_networks"))
  }
  getLayoutNames()
  getLayoutPropertyNames(layout.name='attributes-layout')
  getTableColumnNames(table="node")
  layout.url  = paste0("layout attributes-layout nodeAttribute=clust_id")
  commandsGET(layout.url)
  # Set style name
  style.name = "new_fig_style"
  setVisualStyle(style.name=style.name)
}

}##End of functions

##Start of compute modules##
###################################################################################################################################################
{  ##Step1:: Plot individual networks for each genelist
  
#List adress dictionary
  ad_dic<-data.frame(gene_list_ids,gene_list_adresses,gene_id,stringsAsFactors = FALSE)
  check_files()
#Remove any prior output
  if ("em_results_concat"%in%ls()){
    rm(em_results_concat)}
  for (i in ad_dic$gene_list_ids){   ##Start of individual netowkr plotting loop
    print(paste("Preparing network for", i))
    DEG_adrr<-ad_dic$gene_list_adresses[ad_dic$gene_list_ids==i]
#Read CSV file and assign variables
    DE<-read.csv(DEG_adrr)
    DEGS = head(DE,n = ngene)
    topics_use = ad_dic$gene_id[ad_dic$gene_list_ids==i]
# Pull ranked list and submit for ORA
    DEGS_topic = DEGS %>% pull(topics_use)
    genes = as.character(DEGS_topic)
    em_results_temp <- runGprofiler2(genes)
    #Create a concatenated copy of all em results for output
    em_results_temp$run_id <- run_id
    em_results_temp$list_id <- i
    if (!"em_results_concat"%in%ls()){
     em_results_concat <- em_results_temp
    }else{
       em_results_concat <- rbind(em_results_concat,em_results_temp)
    }
    em_results<-em_results_temp
#Start of plotting module
    submit_enrich_query()
# Run Auto_annotate with clustering for individual networks
    auto_annotate(cluster_boosted = TRUE)

  }#End of the individual network plotting loop
  #specific_save_local<-paste0(save_location,"/",i,"/")
  #em_results_concat$`EnrichmentMap::Genes`<-as.character(em_results_concat$`EnrichmentMap::Genes`)
  #write.csv(em_results_concat,paste0(specific_save_local,"concat_nodes_all.csv"))
}
#########################################################################################################
# Prune each individual network manually
##Go to cytoscape to curate results now##

# Nested loop to create a summary network of current networks and output the node table as csv files
for (z in ad_dic$gene_list_ids){
create_summary_network(network_id = z)
}

#########################################################################################################
# Loop to union networks of interest, and compute union summary this will re-run autoannotate
merge_all_in_ad_dic()

## Curate and manually arrange in cytoscape from here ##

#after union autoannotate
auto_annotate(cluster_boosted = FALSE,cluster_id = "clust_id")

#Style network for ease of arrnagement
column <- "net_id"
cols <- colorRampPalette(brewer.pal(12,"Paired"))(length(ad_dic$gene_list_ids))
setNodeColorMapping(column,mapping.type = "d",colors = cols, table.column.values = ad_dic$gene_list_ids)
setVisualStyle(style.name="default")


#Arrange the clusters in order of interest then set visual style back to normal
# Set style name
style.name = "new_fig_style"
setVisualStyle(style.name=style.name)


#Create summary network for output
create_summary_network(network_id = getNetworkName())


# Create summary network and output both full nodetable and summary table
create_summary_network <-function(network_id=i){
  autoannotate_url <- paste("autoannotate summary includeUnclustered='false network=",network_id)
  current_name <-commandsGET(autoannotate_url)
  summary_nodes <- getTableColumns(table="node",columns=c("name"))
  #clearNodePropertyBypass(node.names = summary_nodes$name,visual.property = "NODE_SIZE")
  style.name = "new_fig_style_summary"
  setVisualStyle(style.name=style.name)
  renameNetwork(title = paste0(network_id,"_summary"), getNetworkName())
  
  #Get the pre-summary node table
  network_id<-getNetworkName()
  nodetable <- as.data.frame(getTableColumns(table= "node" , network = network_id))
  #nodetable$`EnrichmentMap::Genes`<-as.character(nodetable$`EnrichmentMap::Genes`)
  nodetable$EnrichmentMap_Genes<-as.character(nodetable$EnrichmentMap_Genes)
  
  # output the node and edge table of summary network publishing
  current_network_id = getNetworkName()
  nodetable_summary <- as.data.frame(getTableColumns(table= "node" , network = current_network_id))
  nodetable_summary$selected = NULL
  #nodetable_summary$`EnrichmentMap::Genes`<-as.character(nodetable_summary$`EnrichmentMap::Genes`)
  nodetable_summary$EnrichmentMap_Genes<-as.character(nodetable_summary$EnrichmentMap_Genes)
  
  specific_save_local<-paste0(save_location,"/",network_id,"/")
  dir.create(specific_save_local)
  node_file_name = paste0(specific_save_local,network_id,"all_nodes.csv")
  write.csv(nodetable,node_file_name)
  node_file_name = paste0(specific_save_local,network_id,"_summary_node.csv")
  write.csv(nodetable_summary,node_file_name)
}
