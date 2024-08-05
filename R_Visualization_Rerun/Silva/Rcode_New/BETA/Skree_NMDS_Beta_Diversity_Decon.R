rm(list = ls())
library(devtools)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)
library(lubridate)
library(broom)
library(lmtest)
library(car)
library(xlsx)
library(tidyverse)
library(qiime2R)
library(reshape2)
library(readxl)
library(dplyr)
library(openxlsx)
library(DataCombine)
library(ggplot2)
library(tidyr)
library(vegan)
library(tidyverse)
library(ape)
library(dplyr)
library(ggplot2)
library(gplots)
library(lme4)
library(phangorn)
library(plotly)
library(tidyr)
library(vegan)
library(VennDiagram)
library(phyloseq)
library(plyr)
library(hrbrthemes)
library(viridis)
library(gridExtra)
library(shadowtext)
library(ggpubr)
library(gghighlight)
setwd("C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures")
getwd()

#Theme 
#Themes
theme_Publication <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1.2)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            #axis.text.x = element_text(angle = 45),
            axis.text = element_text(size=rel(1.2)), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.box="vertical", 
            legend.margin=margin(),
            legend.direction = "horizontal",
            legend.key.size= unit(2, "cm"),
            legend.text = element_text(size=rel(1.2)),
            #legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic",size=rel(1.2)),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}


scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#7fc97f","#fdb462","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#7fc97f","#fdb462","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

####
####None rarified inputs
#from final table that was created... note this is not the rarefied data, if rarification is need can use vegan package - https://www.rdocumentation.org/packages/vegan/versions/2.4-2/topics/rarefy
table_all<-read_qza("C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/Qiime_Outputs_Rerun/Silva/Removed_classsilva_clust99_Comparison__table-no-mitochondria-chloroplast.qza")

#taxonomy, classification
taxonomy_all<-read_qza("C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/Qiime_Outputs_Rerun/Silva/Classify_silva_clust99_Comparison_classifications.qza")

#alpha diversity evenness
evenness_RL<-read_qza("C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/Qiime_Outputs_Rerun/Silva/Analysis/betadiversity_3500-results/evenness_vector.qza")

#alpha diversity faithpd
faith_RL<-read_qza("C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/Qiime_Outputs_Rerun/Silva/Analysis/phylo-betadiversity_3500-results/faith_pd_vector.qza")

#alpha diversity observed otus
observedotus_RL<-read_qza("C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/Qiime_Outputs_Rerun/Silva/Analysis/phylo-betadiversity_3500-results/observed_features_vector.qza")

#alpha diversity shannon
shannon_RL<-read_qza("C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/Qiime_Outputs_Rerun/Silva/Analysis/betadiversity_3500-results/shannon_vector.qza")

#beta diversity weighted and unweighted unifrac
NMDS_unweighted_unifrac_RL<-read_qza("C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/Qiime_Outputs_Rerun/Silva/Analysis/phylo-betadiversity_3500-results/unweighted_unifrac_NMDS_results.qza")

NMDS_weighted_unifrac_RL<-read_qza("C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/Qiime_Outputs_Rerun/Silva/Analysis/phylo-betadiversity_3500-results/weighted_unifrac_NMDS_results.qza")

#rooted tree
comp_rooted_tree<-read_qza("C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/Qiime_Outputs_Rerun/Silva/Analysis/Comparison_Rerun_silva-rooted-tree.qza")

#Alpha diversity
alpha_diversity<-read_qza("C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/Qiime_Outputs_Rerun/Silva/Analysis/alpha-diversity/alpha_diversity.qza")


#for all comparisons 
rownames(metadata_raw) <- metadata_raw$SampleID
metadata_raw$SampleID <- NULL

#extract data version, one for RA and one for PA 
table_all_data<-as.data.frame(table_all$data)
taxonomy_all_data<-taxonomy_all$data

#Filter OTUs from Deontam
#Decontam
Decontam_OTU<-read.xlsx("C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/QAQC/Decontam_all.xlsx", sheet= 'Sheet1')
Decontam_OTU_filter<-Decontam_OTU[Decontam_OTU$contaminant == 'TRUE',]#can add Blank

table_all_data <- tibble::rownames_to_column(table_all_data, "OTU")
table_all_data <- table_all_data[!(table_all_data$OTU %in% Decontam_OTU_filter$OTU),]
rownames(table_all_data) <- table_all_data[,1]
table_all_data[,1] <- NULL

table_all_data_PA<-table_all_data
table_all_data_PA[]<-as.integer((table_all_data_PA[]!=0))

#parse out taxonomy data so that all columns are individual 
taxonomy_all_data<-parse_taxonomy(taxonomy_all_data)

#merge taxonomy into data table 
Combined_Tax_Data<-transform(merge(taxonomy_all_data,table_all_data,by=0), row.names=Row.names, Row.names=NULL)
Combined_Tax_Data_PA<-transform(merge(taxonomy_all_data,table_all_data_PA,by=0), row.names=Row.names, Row.names=NULL)

#Removed Outliers, note could maybe provide a normalization especially for PR1WEEF
Test.ConnorPrimer58<-as.data.frame(table_all_data$ConnerPrimer58)  ###Removed Connorprimer58 due to low number of reads (3) 

Test.PR1WEEF<-as.data.frame(table_all_data$PR1WEEF)  ###Removed PR1WEEF due to High number of reads (256687) 
sum(Test.PR1WEEF$`table_all_data$PR1WEEF`)


#Read Meta Data 
Metadata_POCPOU<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_DS_Final', rowNames = TRUE)
Metadata_POC<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_Final', rowNames = TRUE)
Metadata_POU<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'DS_Final', rowNames = TRUE)

#Order Meta Data 
Metadata_POCPOU$Classification_1 = factor(Metadata_POCPOU$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))
Metadata_POC$Classification_1 = factor(Metadata_POC$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))
Metadata_POU$Classification_1 = factor(Metadata_POU$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))

######All Samples, POCPOU
Meta.All.POCPOU<-Metadata_POCPOU
Meta.All.POCPOU<-Meta.All.POCPOU[Meta.All.POCPOU$Description != 'ConnerPrimer58',] #removed because of low read count
Meta.All.POCPOU<-Meta.All.POCPOU[Meta.All.POCPOU$Description != 'PR1WEEF',] #removed because of high read count 

Meta.All.POCPOU.Filter<-Meta.All.POCPOU[Meta.All.POCPOU$Classification_1 == 'Potable'| Meta.All.POCPOU$Classification_1 =='Non-potable Reuse'|Meta.All.POCPOU$Classification_1 =='Potable Reuse'|Meta.All.POCPOU$Classification_1 =='',]#can add Blank
#Meta.All.POCPOU.Filter<-Meta.All.POCPOU.Filter[Meta.All.POCPOU.Filter$Matrix == 'water'| Meta.All.POCPOU.Filter$Matrix == 'EXTRACTIONBLANK'| Meta.All.POCPOU.Filter$Matrix == 'FIELDBLANK'| Meta.All.POCPOU.Filter$Matrix == 'PCRBLANK'| Meta.All.POCPOU.Filter$Matrix == 'NA',]
Meta.All.POCPOU.RowNames<-as.data.frame(rownames(Meta.All.POCPOU.Filter))
Meta.All.POCPOU.List<-dplyr::pull(Meta.All.POCPOU.RowNames,1)

OTU.Table.Master<-table_all_data
OTU.Table.Master.Trans<- as.data.frame(t(as.matrix(OTU.Table.Master)))
OTU.All.POCPOU.Filter<- OTU.Table.Master.Trans %>% filter(row.names(OTU.Table.Master.Trans) %in% Meta.All.POCPOU.List)
Meta.All.POCPOU.Filter<-Meta.All.POCPOU.Filter %>% filter(row.names(Meta.All.POCPOU.Filter) %in% row.names(OTU.All.POCPOU.Filter))

Tax.Master<-taxonomy_all$data
Tax.Master<-parse_taxonomy(Tax.Master)

OTU.phy.All.POCPOU=otu_table(as.matrix(OTU.All.POCPOU.Filter), taxa_are_rows=FALSE)
Tax.phy.All.POCPOU = tax_table(as.matrix(Tax.Master))
Meta.phy.All.POCPOU = sample_data(Meta.All.POCPOU.Filter)

physeq.All.POCPOU = phyloseq(OTU.phy.All.POCPOU, Tax.phy.All.POCPOU, Meta.phy.All.POCPOU)

Tree.Master<-comp_rooted_tree$data
edges=phy_tree(Tree.Master)$edge
mycounts = table(edges[,1]) # Source nodes; 1st column of edge matrix
length(mycounts[mycounts ==2]) # Number of nodes with exactly 2 children
length(mycounts[mycounts !=2]) # Number of nodes with more or fewer children
mycounts[mycounts !=2] # How many nodes each of the above has

Tree.Master.V2 <- ape::multi2di(Tree.Master)
edges=phy_tree(Tree.Master.V2)$edge
mycounts = table(edges[,1]) # Source nodes; 1st column of edge matrix
length(mycounts[mycounts ==2]) # Number of nodes with exactly 2 children
length(mycounts[mycounts !=2]) # Number of nodes with more or fewer children
mycounts[mycounts !=2] # How many nodes each of the above has


physeq.combined.All.POCPOU = merge_phyloseq(physeq.All.POCPOU,Tree.Master.V2)
physeq.combined.All.POCPOU

#remove if needed
phy_tree(physeq.combined.All.POCPOU)<-phangorn::midpoint(phy_tree(physeq.combined.All.POCPOU))

physeq.combined.All.POCPOU.rarified3500<-rarefy_even_depth(physeq.combined.All.POCPOU, sample.size = 3500,
                                                        rngseed = 711, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)


######All Samples, POC
Meta.All.POC<-Metadata_POC
Meta.All.POC<-Meta.All.POC[Meta.All.POC$Description != 'ConnerPrimer58',] #removed because of low read count
Meta.All.POC<-Meta.All.POC[Meta.All.POC$Description != 'PR1WEEF',] #removed because of high read count 

Meta.All.POC.Filter<-Meta.All.POC[Meta.All.POC$Classification_1 == 'Potable'| Meta.All.POC$Classification_1 =='Non-potable Reuse'|Meta.All.POC$Classification_1 =='Potable Reuse'|Meta.All.POC$Classification_1 =='',]#can add Blank
#Meta.All.POC.Filter<-Meta.All.POC.Filter[Meta.All.POC.Filter$Matrix == 'water'| Meta.All.POC.Filter$Matrix == 'EXTRACTIONBLANK'| Meta.All.POC.Filter$Matrix == 'FIELDBLANK'| Meta.All.POC.Filter$Matrix == 'PCRBLANK'| Meta.All.POC.Filter$Matrix == 'NA',]
Meta.All.POC.RowNames<-as.data.frame(rownames(Meta.All.POC.Filter))
Meta.All.POC.List<-dplyr::pull(Meta.All.POC.RowNames,1)

OTU.Table.Master<-table_all_data
OTU.Table.Master.Trans<- as.data.frame(t(as.matrix(OTU.Table.Master)))
OTU.All.POC.Filter<- OTU.Table.Master.Trans %>% filter(row.names(OTU.Table.Master.Trans) %in% Meta.All.POC.List)
Meta.All.POC.Filter<-Meta.All.POC.Filter %>% filter(row.names(Meta.All.POC.Filter) %in% row.names(OTU.All.POC.Filter))

Tax.Master<-taxonomy_all$data
Tax.Master<-parse_taxonomy(Tax.Master)

OTU.phy.All.POC=otu_table(as.matrix(OTU.All.POC.Filter), taxa_are_rows=FALSE)
Tax.phy.All.POC = tax_table(as.matrix(Tax.Master))
Meta.phy.All.POC = sample_data(Meta.All.POC.Filter)

physeq.All.POC = phyloseq(OTU.phy.All.POC, Tax.phy.All.POC, Meta.phy.All.POC)

Tree.Master<-comp_rooted_tree$data
edges=phy_tree(Tree.Master)$edge
mycounts = table(edges[,1]) # Source nodes; 1st column of edge matrix
length(mycounts[mycounts ==2]) # Number of nodes with exactly 2 children
length(mycounts[mycounts !=2]) # Number of nodes with more or fewer children
mycounts[mycounts !=2] # How many nodes each of the above has

Tree.Master.V2 <- ape::multi2di(Tree.Master)
edges=phy_tree(Tree.Master.V2)$edge
mycounts = table(edges[,1]) # Source nodes; 1st column of edge matrix
length(mycounts[mycounts ==2]) # Number of nodes with exactly 2 children
length(mycounts[mycounts !=2]) # Number of nodes with more or fewer children
mycounts[mycounts !=2] # How many nodes each of the above has


physeq.combined.All.POC = merge_phyloseq(physeq.All.POC,Tree.Master.V2)
physeq.combined.All.POC

#remove if needed
phy_tree(physeq.combined.All.POC)<-phangorn::midpoint(phy_tree(physeq.combined.All.POC))

physeq.combined.All.POC.rarified3500<-rarefy_even_depth(physeq.combined.All.POC, sample.size = 3500,
                                                           rngseed = 711, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

######All Samples, POU
Meta.All.POU<-Metadata_POU
Meta.All.POU<-Meta.All.POU[Meta.All.POU$Description != 'ConnerPrimer58',] #removed because of low read count
Meta.All.POU<-Meta.All.POU[Meta.All.POU$Description != 'PR1WEEF',] #removed because of high read count 

Meta.All.POU.Filter<-Meta.All.POU[Meta.All.POU$Classification_1 == 'Potable'| Meta.All.POU$Classification_1 =='Non-potable Reuse'|Meta.All.POU$Classification_1 =='Potable Reuse'|Meta.All.POU$Classification_1 =='',]#can add Blank
#Meta.All.POU.Filter<-Meta.All.POU.Filter[Meta.All.POU.Filter$Matrix == 'water'| Meta.All.POU.Filter$Matrix == 'EXTRACTIONBLANK'| Meta.All.POU.Filter$Matrix == 'FIELDBLANK'| Meta.All.POU.Filter$Matrix == 'PCRBLANK'| Meta.All.POU.Filter$Matrix == 'NA',]
Meta.All.POU.RowNames<-as.data.frame(rownames(Meta.All.POU.Filter))
Meta.All.POU.List<-dplyr::pull(Meta.All.POU.RowNames,1)

OTU.Table.Master<-table_all_data
OTU.Table.Master.Trans<- as.data.frame(t(as.matrix(OTU.Table.Master)))
OTU.All.POU.Filter<- OTU.Table.Master.Trans %>% filter(row.names(OTU.Table.Master.Trans) %in% Meta.All.POU.List)
Meta.All.POU.Filter<-Meta.All.POU.Filter %>% filter(row.names(Meta.All.POU.Filter) %in% row.names(OTU.All.POU.Filter))

Tax.Master<-taxonomy_all$data
Tax.Master<-parse_taxonomy(Tax.Master)

OTU.phy.All.POU=otu_table(as.matrix(OTU.All.POU.Filter), taxa_are_rows=FALSE)
Tax.phy.All.POU = tax_table(as.matrix(Tax.Master))
Meta.phy.All.POU = sample_data(Meta.All.POU.Filter)

physeq.All.POU = phyloseq(OTU.phy.All.POU, Tax.phy.All.POU, Meta.phy.All.POU)

Tree.Master<-comp_rooted_tree$data
edges=phy_tree(Tree.Master)$edge
mycounts = table(edges[,1]) # Source nodes; 1st column of edge matrix
length(mycounts[mycounts ==2]) # Number of nodes with exactly 2 children
length(mycounts[mycounts !=2]) # Number of nodes with more or fewer children
mycounts[mycounts !=2] # How many nodes each of the above has

Tree.Master.V2 <- ape::multi2di(Tree.Master)
edges=phy_tree(Tree.Master.V2)$edge
mycounts = table(edges[,1]) # Source nodes; 1st column of edge matrix
length(mycounts[mycounts ==2]) # Number of nodes with exactly 2 children
length(mycounts[mycounts !=2]) # Number of nodes with more or fewer children
mycounts[mycounts !=2] # How many nodes each of the above has


physeq.combined.All.POU = merge_phyloseq(physeq.All.POU,Tree.Master.V2)
physeq.combined.All.POU

#remove if needed
phy_tree(physeq.combined.All.POU)<-phangorn::midpoint(phy_tree(physeq.combined.All.POU))

physeq.combined.All.POU.rarified3500<-rarefy_even_depth(physeq.combined.All.POU, sample.size = 3500,
                                                           rngseed = 711, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)


#Read Meta Data 
Metadata_PotableReusePOCPOU<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_DS_Special', rowNames = TRUE)
Metadata_POC<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_Final', rowNames = TRUE)
Metadata_POU<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'DS_Final', rowNames = TRUE)

#Order Meta Data 
Metadata_PotableReusePOCPOU$Classification_1 = factor(Metadata_PotableReusePOCPOU$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))
Metadata_POC$Classification_1 = factor(Metadata_POC$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))
Metadata_POU$Classification_1 = factor(Metadata_POU$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))

######All Samples, PotableReusePOCPOU
Meta.All.PotableReusePOCPOU<-Metadata_PotableReusePOCPOU
Meta.All.PotableReusePOCPOU<-Meta.All.PotableReusePOCPOU[Meta.All.PotableReusePOCPOU$Description != 'ConnerPrimer58',] #removed because of low read count
Meta.All.PotableReusePOCPOU<-Meta.All.PotableReusePOCPOU[Meta.All.PotableReusePOCPOU$Description != 'PR1WEEF',] #removed because of high read count 

Meta.All.PotableReusePOCPOU.Filter<-Meta.All.PotableReusePOCPOU[Meta.All.PotableReusePOCPOU$Classification_1 == 'Potable Reuse',]#can add Blank
#Meta.All.PotableReusePOCPOU.Filter<-Meta.All.PotableReusePOCPOU.Filter[Meta.All.PotableReusePOCPOU.Filter$Matrix == 'water'| Meta.All.PotableReusePOCPOU.Filter$Matrix == 'EXTRACTIONBLANK'| Meta.All.PotableReusePOCPOU.Filter$Matrix == 'FIELDBLANK'| Meta.All.PotableReusePOCPOU.Filter$Matrix == 'PCRBLANK'| Meta.All.PotableReusePOCPOU.Filter$Matrix == 'NA',]
Meta.All.PotableReusePOCPOU.RowNames<-as.data.frame(rownames(Meta.All.PotableReusePOCPOU.Filter))
Meta.All.PotableReusePOCPOU.List<-dplyr::pull(Meta.All.PotableReusePOCPOU.RowNames,1)

OTU.Table.Master<-table_all_data
OTU.Table.Master.Trans<- as.data.frame(t(as.matrix(OTU.Table.Master)))
OTU.All.PotableReusePOCPOU.Filter<- OTU.Table.Master.Trans %>% filter(row.names(OTU.Table.Master.Trans) %in% Meta.All.PotableReusePOCPOU.List)
Meta.All.PotableReusePOCPOU.Filter<-Meta.All.PotableReusePOCPOU.Filter %>% filter(row.names(Meta.All.PotableReusePOCPOU.Filter) %in% row.names(OTU.All.PotableReusePOCPOU.Filter))

Tax.Master<-taxonomy_all$data
Tax.Master<-parse_taxonomy(Tax.Master)

OTU.phy.All.PotableReusePOCPOU=otu_table(as.matrix(OTU.All.PotableReusePOCPOU.Filter), taxa_are_rows=FALSE)
Tax.phy.All.PotableReusePOCPOU = tax_table(as.matrix(Tax.Master))
Meta.phy.All.PotableReusePOCPOU = sample_data(Meta.All.PotableReusePOCPOU.Filter)

physeq.All.PotableReusePOCPOU = phyloseq(OTU.phy.All.PotableReusePOCPOU, Tax.phy.All.PotableReusePOCPOU, Meta.phy.All.PotableReusePOCPOU)

Tree.Master<-comp_rooted_tree$data
edges=phy_tree(Tree.Master)$edge
mycounts = table(edges[,1]) # Source nodes; 1st column of edge matrix
length(mycounts[mycounts ==2]) # Number of nodes with exactly 2 children
length(mycounts[mycounts !=2]) # Number of nodes with more or fewer children
mycounts[mycounts !=2] # How many nodes each of the above has

Tree.Master.V2 <- ape::multi2di(Tree.Master)
edges=phy_tree(Tree.Master.V2)$edge
mycounts = table(edges[,1]) # Source nodes; 1st column of edge matrix
length(mycounts[mycounts ==2]) # Number of nodes with exactly 2 children
length(mycounts[mycounts !=2]) # Number of nodes with more or fewer children
mycounts[mycounts !=2] # How many nodes each of the above has


physeq.combined.All.PotableReusePOCPOU = merge_phyloseq(physeq.All.PotableReusePOCPOU,Tree.Master.V2)
physeq.combined.All.PotableReusePOCPOU

#remove if needed
phy_tree(physeq.combined.All.PotableReusePOCPOU)<-phangorn::midpoint(phy_tree(physeq.combined.All.PotableReusePOCPOU))

physeq.combined.All.PotableReusePOCPOU.rarified3500<-rarefy_even_depth(physeq.combined.All.PotableReusePOCPOU, sample.size = 3500,
                                                                       rngseed = 711, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)



#Read Meta Data 
Metadata_NonPotableReusePOCPOU<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_DS_Special', rowNames = TRUE)
Metadata_POC<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_Final', rowNames = TRUE)
Metadata_POU<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'DS_Final', rowNames = TRUE)

#Order Meta Data 
Metadata_NonPotableReusePOCPOU$Classification_1 = factor(Metadata_NonPotableReusePOCPOU$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))
Metadata_POC$Classification_1 = factor(Metadata_POC$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))
Metadata_POU$Classification_1 = factor(Metadata_POU$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))

######All Samples, NonPotableReusePOCPOU
Meta.All.NonPotableReusePOCPOU<-Metadata_NonPotableReusePOCPOU
Meta.All.NonPotableReusePOCPOU<-Meta.All.NonPotableReusePOCPOU[Meta.All.NonPotableReusePOCPOU$Description != 'ConnerPrimer58',] #removed because of low read count
Meta.All.NonPotableReusePOCPOU<-Meta.All.NonPotableReusePOCPOU[Meta.All.NonPotableReusePOCPOU$Description != 'PR1WEEF',] #removed because of high read count 

Meta.All.NonPotableReusePOCPOU.Filter<-Meta.All.NonPotableReusePOCPOU[Meta.All.NonPotableReusePOCPOU$Classification_1 == 'Non-potable Reuse',]#can add Blank
#Meta.All.NonPotableReusePOCPOU.Filter<-Meta.All.NonPotableReusePOCPOU.Filter[Meta.All.NonPotableReusePOCPOU.Filter$Matrix == 'water'| Meta.All.NonPotableReusePOCPOU.Filter$Matrix == 'EXTRACTIONBLANK'| Meta.All.NonPotableReusePOCPOU.Filter$Matrix == 'FIELDBLANK'| Meta.All.NonPotableReusePOCPOU.Filter$Matrix == 'PCRBLANK'| Meta.All.NonPotableReusePOCPOU.Filter$Matrix == 'NA',]
Meta.All.NonPotableReusePOCPOU.RowNames<-as.data.frame(rownames(Meta.All.NonPotableReusePOCPOU.Filter))
Meta.All.NonPotableReusePOCPOU.List<-dplyr::pull(Meta.All.NonPotableReusePOCPOU.RowNames,1)

OTU.Table.Master<-table_all_data
OTU.Table.Master.Trans<- as.data.frame(t(as.matrix(OTU.Table.Master)))
OTU.All.NonPotableReusePOCPOU.Filter<- OTU.Table.Master.Trans %>% filter(row.names(OTU.Table.Master.Trans) %in% Meta.All.NonPotableReusePOCPOU.List)
Meta.All.NonPotableReusePOCPOU.Filter<-Meta.All.NonPotableReusePOCPOU.Filter %>% filter(row.names(Meta.All.NonPotableReusePOCPOU.Filter) %in% row.names(OTU.All.NonPotableReusePOCPOU.Filter))

Tax.Master<-taxonomy_all$data
Tax.Master<-parse_taxonomy(Tax.Master)

OTU.phy.All.NonPotableReusePOCPOU=otu_table(as.matrix(OTU.All.NonPotableReusePOCPOU.Filter), taxa_are_rows=FALSE)
Tax.phy.All.NonPotableReusePOCPOU = tax_table(as.matrix(Tax.Master))
Meta.phy.All.NonPotableReusePOCPOU = sample_data(Meta.All.NonPotableReusePOCPOU.Filter)

physeq.All.NonPotableReusePOCPOU = phyloseq(OTU.phy.All.NonPotableReusePOCPOU, Tax.phy.All.NonPotableReusePOCPOU, Meta.phy.All.NonPotableReusePOCPOU)

Tree.Master<-comp_rooted_tree$data
edges=phy_tree(Tree.Master)$edge
mycounts = table(edges[,1]) # Source nodes; 1st column of edge matrix
length(mycounts[mycounts ==2]) # Number of nodes with exactly 2 children
length(mycounts[mycounts !=2]) # Number of nodes with more or fewer children
mycounts[mycounts !=2] # How many nodes each of the above has

Tree.Master.V2 <- ape::multi2di(Tree.Master)
edges=phy_tree(Tree.Master.V2)$edge
mycounts = table(edges[,1]) # Source nodes; 1st column of edge matrix
length(mycounts[mycounts ==2]) # Number of nodes with exactly 2 children
length(mycounts[mycounts !=2]) # Number of nodes with more or fewer children
mycounts[mycounts !=2] # How many nodes each of the above has


physeq.combined.All.NonPotableReusePOCPOU = merge_phyloseq(physeq.All.NonPotableReusePOCPOU,Tree.Master.V2)
physeq.combined.All.NonPotableReusePOCPOU

#remove if needed
phy_tree(physeq.combined.All.NonPotableReusePOCPOU)<-phangorn::midpoint(phy_tree(physeq.combined.All.NonPotableReusePOCPOU))

physeq.combined.All.NonPotableReusePOCPOU.rarified3500<-rarefy_even_depth(physeq.combined.All.NonPotableReusePOCPOU, sample.size = 3500,
                                                                          rngseed = 711, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)


#Read Meta Data 
Metadata_PotablePOCPOU<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_DS_Special', rowNames = TRUE)
Metadata_POC<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_Final', rowNames = TRUE)
Metadata_POU<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'DS_Final', rowNames = TRUE)

#Order Meta Data 
Metadata_PotablePOCPOU$Classification_1 = factor(Metadata_PotablePOCPOU$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))
Metadata_POC$Classification_1 = factor(Metadata_POC$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))
Metadata_POU$Classification_1 = factor(Metadata_POU$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))

######All Samples, PotablePOCPOU
Meta.All.PotablePOCPOU<-Metadata_PotablePOCPOU
Meta.All.PotablePOCPOU<-Meta.All.PotablePOCPOU[Meta.All.PotablePOCPOU$Description != 'ConnerPrimer58',] #removed because of low read count
Meta.All.PotablePOCPOU<-Meta.All.PotablePOCPOU[Meta.All.PotablePOCPOU$Description != 'PR1WEEF',] #removed because of high read count 

Meta.All.PotablePOCPOU.Filter<-Meta.All.PotablePOCPOU[Meta.All.PotablePOCPOU$Classification_1 == 'Potable',]#can add Blank
#Meta.All.PotablePOCPOU.Filter<-Meta.All.PotablePOCPOU.Filter[Meta.All.PotablePOCPOU.Filter$Matrix == 'water'| Meta.All.PotablePOCPOU.Filter$Matrix == 'EXTRACTIONBLANK'| Meta.All.PotablePOCPOU.Filter$Matrix == 'FIELDBLANK'| Meta.All.PotablePOCPOU.Filter$Matrix == 'PCRBLANK'| Meta.All.PotablePOCPOU.Filter$Matrix == 'NA',]
Meta.All.PotablePOCPOU.RowNames<-as.data.frame(rownames(Meta.All.PotablePOCPOU.Filter))
Meta.All.PotablePOCPOU.List<-dplyr::pull(Meta.All.PotablePOCPOU.RowNames,1)

OTU.Table.Master<-table_all_data
OTU.Table.Master.Trans<- as.data.frame(t(as.matrix(OTU.Table.Master)))
OTU.All.PotablePOCPOU.Filter<- OTU.Table.Master.Trans %>% filter(row.names(OTU.Table.Master.Trans) %in% Meta.All.PotablePOCPOU.List)
Meta.All.PotablePOCPOU.Filter<-Meta.All.PotablePOCPOU.Filter %>% filter(row.names(Meta.All.PotablePOCPOU.Filter) %in% row.names(OTU.All.PotablePOCPOU.Filter))

Tax.Master<-taxonomy_all$data
Tax.Master<-parse_taxonomy(Tax.Master)

OTU.phy.All.PotablePOCPOU=otu_table(as.matrix(OTU.All.PotablePOCPOU.Filter), taxa_are_rows=FALSE)
Tax.phy.All.PotablePOCPOU = tax_table(as.matrix(Tax.Master))
Meta.phy.All.PotablePOCPOU = sample_data(Meta.All.PotablePOCPOU.Filter)

physeq.All.PotablePOCPOU = phyloseq(OTU.phy.All.PotablePOCPOU, Tax.phy.All.PotablePOCPOU, Meta.phy.All.PotablePOCPOU)

Tree.Master<-comp_rooted_tree$data
edges=phy_tree(Tree.Master)$edge
mycounts = table(edges[,1]) # Source nodes; 1st column of edge matrix
length(mycounts[mycounts ==2]) # Number of nodes with exactly 2 children
length(mycounts[mycounts !=2]) # Number of nodes with more or fewer children
mycounts[mycounts !=2] # How many nodes each of the above has

Tree.Master.V2 <- ape::multi2di(Tree.Master)
edges=phy_tree(Tree.Master.V2)$edge
mycounts = table(edges[,1]) # Source nodes; 1st column of edge matrix
length(mycounts[mycounts ==2]) # Number of nodes with exactly 2 children
length(mycounts[mycounts !=2]) # Number of nodes with more or fewer children
mycounts[mycounts !=2] # How many nodes each of the above has


physeq.combined.All.PotablePOCPOU = merge_phyloseq(physeq.All.PotablePOCPOU,Tree.Master.V2)
physeq.combined.All.PotablePOCPOU

#remove if needed
phy_tree(physeq.combined.All.PotablePOCPOU)<-phangorn::midpoint(phy_tree(physeq.combined.All.PotablePOCPOU))

physeq.combined.All.PotablePOCPOU.rarified3500<-rarefy_even_depth(physeq.combined.All.PotablePOCPOU, sample.size = 3500,
                                                                  rngseed = 711, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)


#Read Meta Data 
Metadata_PotableBothPOCPOU<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_DS_Special', rowNames = TRUE)
Metadata_POC<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_Final', rowNames = TRUE)
Metadata_POU<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'DS_Final', rowNames = TRUE)

#Order Meta Data 
Metadata_PotableBothPOCPOU$Classification_1 = factor(Metadata_PotableBothPOCPOU$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))
Metadata_POC$Classification_1 = factor(Metadata_POC$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))
Metadata_POU$Classification_1 = factor(Metadata_POU$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))

######All Samples, PotableBothPOCPOU
Meta.All.PotableBothPOCPOU<-Metadata_PotableBothPOCPOU
Meta.All.PotableBothPOCPOU<-Meta.All.PotableBothPOCPOU[Meta.All.PotableBothPOCPOU$Description != 'ConnerPrimer58',] #removed because of low read count
Meta.All.PotableBothPOCPOU<-Meta.All.PotableBothPOCPOU[Meta.All.PotableBothPOCPOU$Description != 'PR1WEEF',] #removed because of high read count 

Meta.All.PotableBothPOCPOU.Filter<-Meta.All.PotableBothPOCPOU[Meta.All.PotableBothPOCPOU$Classification_1 == 'Potable'| Meta.All.PotableBothPOCPOU$Classification_1 == 'Potable Reuse',]#can add Blank
#Meta.All.PotableBothPOCPOU.Filter<-Meta.All.PotableBothPOCPOU.Filter[Meta.All.PotableBothPOCPOU.Filter$Matrix == 'water'| Meta.All.PotableBothPOCPOU.Filter$Matrix == 'EXTRACTIONBLANK'| Meta.All.PotableBothPOCPOU.Filter$Matrix == 'FIELDBLANK'| Meta.All.PotableBothPOCPOU.Filter$Matrix == 'PCRBLANK'| Meta.All.PotableBothPOCPOU.Filter$Matrix == 'NA',]
Meta.All.PotableBothPOCPOU.RowNames<-as.data.frame(rownames(Meta.All.PotableBothPOCPOU.Filter))
Meta.All.PotableBothPOCPOU.List<-dplyr::pull(Meta.All.PotableBothPOCPOU.RowNames,1)

OTU.Table.Master<-table_all_data
OTU.Table.Master.Trans<- as.data.frame(t(as.matrix(OTU.Table.Master)))
OTU.All.PotableBothPOCPOU.Filter<- OTU.Table.Master.Trans %>% filter(row.names(OTU.Table.Master.Trans) %in% Meta.All.PotableBothPOCPOU.List)
Meta.All.PotableBothPOCPOU.Filter<-Meta.All.PotableBothPOCPOU.Filter %>% filter(row.names(Meta.All.PotableBothPOCPOU.Filter) %in% row.names(OTU.All.PotableBothPOCPOU.Filter))

Tax.Master<-taxonomy_all$data
Tax.Master<-parse_taxonomy(Tax.Master)

OTU.phy.All.PotableBothPOCPOU=otu_table(as.matrix(OTU.All.PotableBothPOCPOU.Filter), taxa_are_rows=FALSE)
Tax.phy.All.PotableBothPOCPOU = tax_table(as.matrix(Tax.Master))
Meta.phy.All.PotableBothPOCPOU = sample_data(Meta.All.PotableBothPOCPOU.Filter)

physeq.All.PotableBothPOCPOU = phyloseq(OTU.phy.All.PotableBothPOCPOU, Tax.phy.All.PotableBothPOCPOU, Meta.phy.All.PotableBothPOCPOU)

Tree.Master<-comp_rooted_tree$data
edges=phy_tree(Tree.Master)$edge
mycounts = table(edges[,1]) # Source nodes; 1st column of edge matrix
length(mycounts[mycounts ==2]) # Number of nodes with exactly 2 children
length(mycounts[mycounts !=2]) # Number of nodes with more or fewer children
mycounts[mycounts !=2] # How many nodes each of the above has

Tree.Master.V2 <- ape::multi2di(Tree.Master)
edges=phy_tree(Tree.Master.V2)$edge
mycounts = table(edges[,1]) # Source nodes; 1st column of edge matrix
length(mycounts[mycounts ==2]) # Number of nodes with exactly 2 children
length(mycounts[mycounts !=2]) # Number of nodes with more or fewer children
mycounts[mycounts !=2] # How many nodes each of the above has


physeq.combined.All.PotableBothPOCPOU = merge_phyloseq(physeq.All.PotableBothPOCPOU,Tree.Master.V2)
physeq.combined.All.PotableBothPOCPOU

#remove if needed
phy_tree(physeq.combined.All.PotableBothPOCPOU)<-phangorn::midpoint(phy_tree(physeq.combined.All.PotableBothPOCPOU))

physeq.combined.All.PotableBothPOCPOU.rarified3500<-rarefy_even_depth(physeq.combined.All.PotableBothPOCPOU, sample.size = 3500,
                                                                      rngseed = 711, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)


###Run Skree Plots 
#devtools::install_github("fvlampe/goeveg")
library(goeveg)

skree_all_POCPOU<-screeplot_NMDS(
  OTU.phy.All.POCPOU,
  distance = "bray",
  k = 6,
  trymax = 20,
  autotransform = TRUE
)


skree_all_POC<-screeplot_NMDS(
  OTU.phy.All.POC,
  distance = "bray",
  k = 6,
  trymax = 20,
  autotransform = TRUE
)

skree_all_POU<-screeplot_NMDS(
  OTU.phy.All.POU,
  distance = "bray",
  k = 6,
  trymax = 20,
  autotransform = TRUE
)

skree_NonpotableReuse_POCPOU<-screeplot_NMDS(
  OTU.phy.All.NonPotableReusePOCPOU,
  distance = "bray",
  k = 6,
  trymax = 20,
  autotransform = TRUE
)

skree_Potable_POCPOU<-screeplot_NMDS(
  OTU.phy.All.PotablePOCPOU,
  distance = "bray",
  k = 6,
  trymax = 20,
  autotransform = TRUE
)

skree_PotableBoth_POCPOU<-screeplot_NMDS(
  OTU.phy.All.PotableBothPOCPOU,
  distance = "bray",
  k = 6,
  trymax = 20,
  autotransform = TRUE
)

skree_PotableReuse_POCPOU<-screeplot_NMDS(
  OTU.phy.All.PotableReusePOCPOU,
  distance = "bray",
  k = 6,
  trymax = 20,
  autotransform = TRUE
)



##Export
skree_all_POCPOU
skree_all_POC
skree_all_POU
skree_NonpotableReuse_POCPOU
skree_Potable_POCPOU
skree_PotableBoth_POCPOU
skree_PotableReuse_POCPOU


df.skree<-data.frame(skree_all_POCPOU,
                    skree_all_POC,
                    skree_all_POU,
                    skree_NonpotableReuse_POCPOU,
                    skree_Potable_POCPOU,
                    skree_PotableBoth_POCPOU,
                    skree_PotableReuse_POCPOU)



##Export
write.xlsx(df.skree, file ="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/Skree/skree.xlsx", sheetName = "Sheet1", 
           colNames = TRUE, rowNames = TRUE, append = FALSE)


