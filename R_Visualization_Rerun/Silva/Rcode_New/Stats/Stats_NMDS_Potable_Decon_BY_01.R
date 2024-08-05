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
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#7fc97f","#fdb462","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33","#6d4b4b","#333333","#599e94")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#7fc97f","#fdb462","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33","#6d4b4b","#333333","#599e94")), ...)
  
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
Metadata_PotablePOCPOU<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_DS_Special_V2', rowNames = TRUE)
Metadata_POC<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_Final', rowNames = TRUE)
Metadata_POU<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'DS_Final', rowNames = TRUE)

#Order Meta Data 
Metadata_PotablePOCPOU$Classification_1 = factor(Metadata_PotablePOCPOU$Classification_1, levels=c("Potable Conventional","Potable Reuse","Non-potable Reuse","Blank"))
Metadata_POC$Classification_1 = factor(Metadata_POC$Classification_1, levels=c("Potable Conventional","Potable Reuse","Non-potable Reuse","Blank"))
Metadata_POU$Classification_1 = factor(Metadata_POU$Classification_1, levels=c("Potable Conventional","Potable Reuse","Non-potable Reuse","Blank"))

######All Samples, PotablePOCPOU
Meta.All.PotablePOCPOU<-Metadata_PotablePOCPOU
Meta.All.PotablePOCPOU<-Meta.All.PotablePOCPOU[Meta.All.PotablePOCPOU$Description != 'ConnerPrimer58',] #removed because of low read count
Meta.All.PotablePOCPOU<-Meta.All.PotablePOCPOU[Meta.All.PotablePOCPOU$Description != 'PR1WEEF',] #removed because of high read count 

Meta.All.PotablePOCPOU.Filter<-Meta.All.PotablePOCPOU[Meta.All.PotablePOCPOU$Classification_1 == 'Potable Conventional',]#can add Blank
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

#################
#################
#################
##################################
#################
#################
##################################
#################
#################
##################################
#################
#################
##################################
#################
#################
#################
#Need to update Metadata file 

#Anosim/Adonis/Betadispersity 
#Can adjust permutations for higher accuracy. 
permutations=1000
#Setup Input Files 

All.PotablePOCPOU.dist.bray<- phyloseq::distance(physeq.combined.All.PotablePOCPOU, method = "bray")
All.PotablePOCPOU.meta<- data.frame(sample_data(physeq.combined.All.PotablePOCPOU))

ANOSIM.GroupedP.All.PotablePOCPOU.bray<-numeric() 
ADONIS2.GroupedP.All.PotablePOCPOU.bray<-numeric()
Betadis.GroupedP.All.PotablePOCPOU.bray<-numeric()
ANOSIM.GroupedR.All.PotablePOCPOU.bray<-numeric() 
ADONIS2.GroupedR.All.PotablePOCPOU.bray<-numeric()


# #################
# #################
# #################
# #################
# #Class1_Pot_NonPot 
# ANOSIM.All.PotablePOCPOU.bray.Class1_Pot_NonPot<-anosim(All.PotablePOCPOU.dist.bray,
#                                All.PotablePOCPOU.meta$Class1_Pot_NonPot, 
#                                permutations = permutations)
# 
# ADONIS2.All.PotablePOCPOU.bray.Class1_Pot_NonPot<-adonis2(All.PotablePOCPOU.dist.bray ~ Class1_Pot_NonPot, 
#                                  data = All.PotablePOCPOU.meta, 
#                                  permutations = permutations)
# 
# Betadis.All.PotablePOCPOU.bray.Class1_Pot_NonPot<-anova(betadisper(All.PotablePOCPOU.dist.bray,
#                                                       All.PotablePOCPOU.meta$Class1_Pot_NonPot))
# 
# 
# 
# ANOSIM.Pval.All.PotablePOCPOU.bray.Class1_Pot_NonPot<-ANOSIM.All.PotablePOCPOU.bray.Class1_Pot_NonPot$signif
# ANOSIM.Pval.All.PotablePOCPOU.bray.Class1_Pot_NonPot
# 
# ADONIS2.Pval.All.PotablePOCPOU.bray.Class1_Pot_NonPot<-ADONIS2.All.PotablePOCPOU.bray.Class1_Pot_NonPot$`Pr(>F)`[1]
# ADONIS2.Pval.All.PotablePOCPOU.bray.Class1_Pot_NonPot  
# 
# Betadis.Pval.All.PotablePOCPOU.bray.Class1_Pot_NonPot<-Betadis.All.PotablePOCPOU.bray.Class1_Pot_NonPot$`Pr(>F)`[1]
# Betadis.Pval.All.PotablePOCPOU.bray.Class1_Pot_NonPot
# 
# ANOSIM.R.All.PotablePOCPOU.bray.Class1_Pot_NonPot<-ANOSIM.All.PotablePOCPOU.bray.Class1_Pot_NonPot$statistic
# ANOSIM.R.All.PotablePOCPOU.bray.Class1_Pot_NonPot
# 
# ADONIS2.R2.All.PotablePOCPOU.bray.Class1_Pot_NonPot<-ADONIS2.All.PotablePOCPOU.bray.Class1_Pot_NonPot$R2[1]
# ADONIS2.R2.All.PotablePOCPOU.bray.Class1_Pot_NonPot  
# 
# #Pull Together Pvalues to adjust 
# ANOSIM.GroupedP.All.PotablePOCPOU.bray["Class1_Pot_NonPot"]<-ANOSIM.Pval.All.PotablePOCPOU.bray.Class1_Pot_NonPot
# ADONIS2.GroupedP.All.PotablePOCPOU.bray["Class1_Pot_NonPot"]<-ADONIS2.Pval.All.PotablePOCPOU.bray.Class1_Pot_NonPot
# 
# Betadis.GroupedP.All.PotablePOCPOU.bray["Class1_Pot_NonPot"]<-Betadis.Pval.All.PotablePOCPOU.bray.Class1_Pot_NonPot
# 
# ANOSIM.GroupedR.All.PotablePOCPOU.bray["Class1_Pot_NonPot"]<-ANOSIM.R.All.PotablePOCPOU.bray.Class1_Pot_NonPot
# ADONIS2.GroupedR.All.PotablePOCPOU.bray["Class1_Pot_NonPot"]<-ADONIS2.R2.All.PotablePOCPOU.bray.Class1_Pot_NonPot
# 
# #Pairwise if needed, check if interesting or not. Remove if not. 
# count(All.PotablePOCPOU.meta$Class1_Pot_NonPot)
# 
# ADONIS2.Pairwise.All.PotablePOCPOU.bray.Class1_Pot_NonPot<-pairwise.adonis2(All.PotablePOCPOU.dist.bray ~ Class1_Pot_NonPot, 
#                                                                     data = All.PotablePOCPOU.meta, 
#                                                                     permutations = permutations)
# ADONIS2.Pairwise.All.PotablePOCPOU.bray.Class1_Pot_NonPot
# 
# 
# 
# 
# #################
# #################
# #################
# #################
# #Classification_1 
# ANOSIM.All.PotablePOCPOU.bray.Classification_1<-anosim(All.PotablePOCPOU.dist.bray,
#                                                  All.PotablePOCPOU.meta$Classification_1, 
#                                                  permutations = permutations)
# 
# ADONIS2.All.PotablePOCPOU.bray.Classification_1<-adonis2(All.PotablePOCPOU.dist.bray ~ Classification_1, 
#                                                    data = All.PotablePOCPOU.meta, 
#                                                    permutations = permutations)
# 
# Betadis.All.PotablePOCPOU.bray.Classification_1<-anova(betadisper(All.PotablePOCPOU.dist.bray,
#                                                             All.PotablePOCPOU.meta$Classification_1))
# 
# 
# 
# ANOSIM.Pval.All.PotablePOCPOU.bray.Classification_1<-ANOSIM.All.PotablePOCPOU.bray.Classification_1$signif
# ANOSIM.Pval.All.PotablePOCPOU.bray.Classification_1
# 
# ADONIS2.Pval.All.PotablePOCPOU.bray.Classification_1<-ADONIS2.All.PotablePOCPOU.bray.Classification_1$`Pr(>F)`[1]
# ADONIS2.Pval.All.PotablePOCPOU.bray.Classification_1  
# 
# Betadis.Pval.All.PotablePOCPOU.bray.Classification_1<-Betadis.All.PotablePOCPOU.bray.Classification_1$`Pr(>F)`[1]
# Betadis.Pval.All.PotablePOCPOU.bray.Classification_1
# 
# ANOSIM.R.All.PotablePOCPOU.bray.Classification_1<-ANOSIM.All.PotablePOCPOU.bray.Classification_1$statistic
# ANOSIM.R.All.PotablePOCPOU.bray.Classification_1
# 
# ADONIS2.R2.All.PotablePOCPOU.bray.Classification_1<-ADONIS2.All.PotablePOCPOU.bray.Classification_1$R2[1]
# ADONIS2.R2.All.PotablePOCPOU.bray.Classification_1  
# 
# #Pull Together Pvalues to adjust 
# ANOSIM.GroupedP.All.PotablePOCPOU.bray["Classification_1"]<-ANOSIM.Pval.All.PotablePOCPOU.bray.Classification_1
# ADONIS2.GroupedP.All.PotablePOCPOU.bray["Classification_1"]<-ADONIS2.Pval.All.PotablePOCPOU.bray.Classification_1
# 
# Betadis.GroupedP.All.PotablePOCPOU.bray["Classification_1"]<-Betadis.Pval.All.PotablePOCPOU.bray.Classification_1
# 
# ANOSIM.GroupedR.All.PotablePOCPOU.bray["Classification_1"]<-ANOSIM.R.All.PotablePOCPOU.bray.Classification_1
# ADONIS2.GroupedR.All.PotablePOCPOU.bray["Classification_1"]<-ADONIS2.R2.All.PotablePOCPOU.bray.Classification_1
# 
# #Pairwise if needed, check if interesting or not. Remove if not. 
# count(All.PotablePOCPOU.meta$Classification_1)
# 
# ADONIS2.Pairwise.All.PotablePOCPOU.bray.Classification_1<-pairwise.adonis2(All.PotablePOCPOU.dist.bray ~ Classification_1, 
#                                                                      data = All.PotablePOCPOU.meta, 
#                                                                      permutations = permutations)
# ADONIS2.Pairwise.All.PotablePOCPOU.bray.Classification_1
# 
# 
# 
# 
# 
# #################
# #################
# #################
# #################
# #Classification_1_mod 
# ANOSIM.All.PotablePOCPOU.bray.Classification_1_mod<-anosim(All.PotablePOCPOU.dist.bray,
#                                                  All.PotablePOCPOU.meta$Classification_1_mod, 
#                                                  permutations = permutations)
# 
# ADONIS2.All.PotablePOCPOU.bray.Classification_1_mod<-adonis2(All.PotablePOCPOU.dist.bray ~ Classification_1_mod, 
#                                                    data = All.PotablePOCPOU.meta, 
#                                                    permutations = permutations)
# 
# Betadis.All.PotablePOCPOU.bray.Classification_1_mod<-anova(betadisper(All.PotablePOCPOU.dist.bray,
#                                                             All.PotablePOCPOU.meta$Classification_1_mod))
# 
# 
# 
# ANOSIM.Pval.All.PotablePOCPOU.bray.Classification_1_mod<-ANOSIM.All.PotablePOCPOU.bray.Classification_1_mod$signif
# ANOSIM.Pval.All.PotablePOCPOU.bray.Classification_1_mod
# 
# ADONIS2.Pval.All.PotablePOCPOU.bray.Classification_1_mod<-ADONIS2.All.PotablePOCPOU.bray.Classification_1_mod$`Pr(>F)`[1]
# ADONIS2.Pval.All.PotablePOCPOU.bray.Classification_1_mod  
# 
# Betadis.Pval.All.PotablePOCPOU.bray.Classification_1_mod<-Betadis.All.PotablePOCPOU.bray.Classification_1_mod$`Pr(>F)`[1]
# Betadis.Pval.All.PotablePOCPOU.bray.Classification_1_mod
# 
# ANOSIM.R.All.PotablePOCPOU.bray.Classification_1_mod<-ANOSIM.All.PotablePOCPOU.bray.Classification_1_mod$statistic
# ANOSIM.R.All.PotablePOCPOU.bray.Classification_1_mod
# 
# ADONIS2.R2.All.PotablePOCPOU.bray.Classification_1_mod<-ADONIS2.All.PotablePOCPOU.bray.Classification_1_mod$R2[1]
# ADONIS2.R2.All.PotablePOCPOU.bray.Classification_1_mod  
# 
# #Pull Together Pvalues to adjust 
# ANOSIM.GroupedP.All.PotablePOCPOU.bray["Classification_1_mod"]<-ANOSIM.Pval.All.PotablePOCPOU.bray.Classification_1_mod
# ADONIS2.GroupedP.All.PotablePOCPOU.bray["Classification_1_mod"]<-ADONIS2.Pval.All.PotablePOCPOU.bray.Classification_1_mod
# 
# Betadis.GroupedP.All.PotablePOCPOU.bray["Classification_1_mod"]<-Betadis.Pval.All.PotablePOCPOU.bray.Classification_1_mod
# 
# ANOSIM.GroupedR.All.PotablePOCPOU.bray["Classification_1_mod"]<-ANOSIM.R.All.PotablePOCPOU.bray.Classification_1_mod
# ADONIS2.GroupedR.All.PotablePOCPOU.bray["Classification_1_mod"]<-ADONIS2.R2.All.PotablePOCPOU.bray.Classification_1_mod
# 
# #Pairwise if needed, check if interesting or not. Remove if not. 
# count(All.PotablePOCPOU.meta$Classification_1_mod)
# 
# ADONIS2.Pairwise.All.PotablePOCPOU.bray.Classification_1_mod<-pairwise.adonis2(All.PotablePOCPOU.dist.bray ~ Classification_1_mod, 
#                                                                      data = All.PotablePOCPOU.meta, 
#                                                                      permutations = permutations)
# ADONIS2.Pairwise.All.PotablePOCPOU.bray.Classification_1_mod
# 
# 
# 
# 
# 
# #################
# #################
# #################
# #################
# #Classification_2 
# ANOSIM.All.PotablePOCPOU.bray.Classification_2<-anosim(All.PotablePOCPOU.dist.bray,
#                                                  All.PotablePOCPOU.meta$Classification_2, 
#                                                  permutations = permutations)
# 
# ADONIS2.All.PotablePOCPOU.bray.Classification_2<-adonis2(All.PotablePOCPOU.dist.bray ~ Classification_2, 
#                                                    data = All.PotablePOCPOU.meta, 
#                                                    permutations = permutations)
# 
# Betadis.All.PotablePOCPOU.bray.Classification_2<-anova(betadisper(All.PotablePOCPOU.dist.bray,
#                                                             All.PotablePOCPOU.meta$Classification_2))
# 
# 
# 
# ANOSIM.Pval.All.PotablePOCPOU.bray.Classification_2<-ANOSIM.All.PotablePOCPOU.bray.Classification_2$signif
# ANOSIM.Pval.All.PotablePOCPOU.bray.Classification_2
# 
# ADONIS2.Pval.All.PotablePOCPOU.bray.Classification_2<-ADONIS2.All.PotablePOCPOU.bray.Classification_2$`Pr(>F)`[1]
# ADONIS2.Pval.All.PotablePOCPOU.bray.Classification_2  
# 
# Betadis.Pval.All.PotablePOCPOU.bray.Classification_2<-Betadis.All.PotablePOCPOU.bray.Classification_2$`Pr(>F)`[1]
# Betadis.Pval.All.PotablePOCPOU.bray.Classification_2
# 
# ANOSIM.R.All.PotablePOCPOU.bray.Classification_2<-ANOSIM.All.PotablePOCPOU.bray.Classification_2$statistic
# ANOSIM.R.All.PotablePOCPOU.bray.Classification_2
# 
# ADONIS2.R2.All.PotablePOCPOU.bray.Classification_2<-ADONIS2.All.PotablePOCPOU.bray.Classification_2$R2[1]
# ADONIS2.R2.All.PotablePOCPOU.bray.Classification_2  
# 
# #Pull Together Pvalues to adjust 
# ANOSIM.GroupedP.All.PotablePOCPOU.bray["Classification_2"]<-ANOSIM.Pval.All.PotablePOCPOU.bray.Classification_2
# ADONIS2.GroupedP.All.PotablePOCPOU.bray["Classification_2"]<-ADONIS2.Pval.All.PotablePOCPOU.bray.Classification_2
# 
# Betadis.GroupedP.All.PotablePOCPOU.bray["Classification_2"]<-Betadis.Pval.All.PotablePOCPOU.bray.Classification_2
# 
# ANOSIM.GroupedR.All.PotablePOCPOU.bray["Classification_2"]<-ANOSIM.R.All.PotablePOCPOU.bray.Classification_2
# ADONIS2.GroupedR.All.PotablePOCPOU.bray["Classification_2"]<-ADONIS2.R2.All.PotablePOCPOU.bray.Classification_2
# 
# #Pairwise if needed, check if interesting or not. Remove if not. 
# count(All.PotablePOCPOU.meta$Classification_2)
# 
# ADONIS2.Pairwise.All.PotablePOCPOU.bray.Classification_2<-pairwise.adonis2(All.PotablePOCPOU.dist.bray ~ Classification_2, 
#                                                                      data = All.PotablePOCPOU.meta, 
#                                                                      permutations = permutations)
# ADONIS2.Pairwise.All.PotablePOCPOU.bray.Classification_2
# #Classification_2





# #################
# #################
# #################
# #################
# #Climate_Region 
# ANOSIM.All.PotablePOCPOU.bray.Climate_Region<-anosim(All.PotablePOCPOU.dist.bray,
#                                                  All.PotablePOCPOU.meta$Climate_Region, 
#                                                  permutations = permutations)
# 
# ADONIS2.All.PotablePOCPOU.bray.Climate_Region<-adonis2(All.PotablePOCPOU.dist.bray ~ Climate_Region, 
#                                                    data = All.PotablePOCPOU.meta, 
#                                                    permutations = permutations)
# 
# Betadis.All.PotablePOCPOU.bray.Climate_Region<-anova(betadisper(All.PotablePOCPOU.dist.bray,
#                                                             All.PotablePOCPOU.meta$Climate_Region))
# 
# 
# 
# ANOSIM.Pval.All.PotablePOCPOU.bray.Climate_Region<-ANOSIM.All.PotablePOCPOU.bray.Climate_Region$signif
# ANOSIM.Pval.All.PotablePOCPOU.bray.Climate_Region
# 
# ADONIS2.Pval.All.PotablePOCPOU.bray.Climate_Region<-ADONIS2.All.PotablePOCPOU.bray.Climate_Region$`Pr(>F)`[1]
# ADONIS2.Pval.All.PotablePOCPOU.bray.Climate_Region  
# 
# Betadis.Pval.All.PotablePOCPOU.bray.Climate_Region<-Betadis.All.PotablePOCPOU.bray.Climate_Region$`Pr(>F)`[1]
# Betadis.Pval.All.PotablePOCPOU.bray.Climate_Region
# 
# ANOSIM.R.All.PotablePOCPOU.bray.Climate_Region<-ANOSIM.All.PotablePOCPOU.bray.Climate_Region$statistic
# ANOSIM.R.All.PotablePOCPOU.bray.Climate_Region
# 
# ADONIS2.R2.All.PotablePOCPOU.bray.Climate_Region<-ADONIS2.All.PotablePOCPOU.bray.Climate_Region$R2[1]
# ADONIS2.R2.All.PotablePOCPOU.bray.Climate_Region  
# 
# #Pull Together Pvalues to adjust 
# ANOSIM.GroupedP.All.PotablePOCPOU.bray["Climate_Region"]<-ANOSIM.Pval.All.PotablePOCPOU.bray.Climate_Region
# ADONIS2.GroupedP.All.PotablePOCPOU.bray["Climate_Region"]<-ADONIS2.Pval.All.PotablePOCPOU.bray.Climate_Region
# 
# Betadis.GroupedP.All.PotablePOCPOU.bray["Climate_Region"]<-Betadis.Pval.All.PotablePOCPOU.bray.Climate_Region
# 
# ANOSIM.GroupedR.All.PotablePOCPOU.bray["Climate_Region"]<-ANOSIM.R.All.PotablePOCPOU.bray.Climate_Region
# ADONIS2.GroupedR.All.PotablePOCPOU.bray["Climate_Region"]<-ADONIS2.R2.All.PotablePOCPOU.bray.Climate_Region
# 
# #Pairwise if needed, check if interesting or not. Remove if not. 
# count(All.PotablePOCPOU.meta$Climate_Region)
# 
# ADONIS2.Pairwise.All.PotablePOCPOU.bray.Climate_Region<-pairwise.adonis2(All.PotablePOCPOU.dist.bray ~ Climate_Region, 
#                                                                      data = All.PotablePOCPOU.meta, 
#                                                                      permutations = permutations)
# ADONIS2.Pairwise.All.PotablePOCPOU.bray.Climate_Region
# #Climate_Region





#################
#################
#################
#################
#Climate 
ANOSIM.All.PotablePOCPOU.bray.Climate<-anosim(All.PotablePOCPOU.dist.bray,
                                                 All.PotablePOCPOU.meta$Climate, 
                                                 permutations = permutations)

ADONIS2.All.PotablePOCPOU.bray.Climate<-adonis2(All.PotablePOCPOU.dist.bray ~ Climate, 
                                                   data = All.PotablePOCPOU.meta, 
                                                   permutations = permutations)

Betadis.All.PotablePOCPOU.bray.Climate<-anova(betadisper(All.PotablePOCPOU.dist.bray,
                                                            All.PotablePOCPOU.meta$Climate))



ANOSIM.Pval.All.PotablePOCPOU.bray.Climate<-ANOSIM.All.PotablePOCPOU.bray.Climate$signif
ANOSIM.Pval.All.PotablePOCPOU.bray.Climate

ADONIS2.Pval.All.PotablePOCPOU.bray.Climate<-ADONIS2.All.PotablePOCPOU.bray.Climate$`Pr(>F)`[1]
ADONIS2.Pval.All.PotablePOCPOU.bray.Climate  

Betadis.Pval.All.PotablePOCPOU.bray.Climate<-Betadis.All.PotablePOCPOU.bray.Climate$`Pr(>F)`[1]
Betadis.Pval.All.PotablePOCPOU.bray.Climate

ANOSIM.R.All.PotablePOCPOU.bray.Climate<-ANOSIM.All.PotablePOCPOU.bray.Climate$statistic
ANOSIM.R.All.PotablePOCPOU.bray.Climate

ADONIS2.R2.All.PotablePOCPOU.bray.Climate<-ADONIS2.All.PotablePOCPOU.bray.Climate$R2[1]
ADONIS2.R2.All.PotablePOCPOU.bray.Climate  

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.PotablePOCPOU.bray["Climate"]<-ANOSIM.Pval.All.PotablePOCPOU.bray.Climate
ADONIS2.GroupedP.All.PotablePOCPOU.bray["Climate"]<-ADONIS2.Pval.All.PotablePOCPOU.bray.Climate

Betadis.GroupedP.All.PotablePOCPOU.bray["Climate"]<-Betadis.Pval.All.PotablePOCPOU.bray.Climate

ANOSIM.GroupedR.All.PotablePOCPOU.bray["Climate"]<-ANOSIM.R.All.PotablePOCPOU.bray.Climate
ADONIS2.GroupedR.All.PotablePOCPOU.bray["Climate"]<-ADONIS2.R2.All.PotablePOCPOU.bray.Climate

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.PotablePOCPOU.meta$Climate)

ADONIS2.Pairwise.All.PotablePOCPOU.bray.Climate<-pairwise.adonis2(All.PotablePOCPOU.dist.bray ~ Climate, 
                                                                     data = All.PotablePOCPOU.meta, 
                                                                     permutations = permutations)
ADONIS2.Pairwise.All.PotablePOCPOU.bray.Climate
#Climate





#################
#################
#################
#################
#Region 
ANOSIM.All.PotablePOCPOU.bray.Region<-anosim(All.PotablePOCPOU.dist.bray,
                                                 All.PotablePOCPOU.meta$Region, 
                                                 permutations = permutations)

ADONIS2.All.PotablePOCPOU.bray.Region<-adonis2(All.PotablePOCPOU.dist.bray ~ Region, 
                                                   data = All.PotablePOCPOU.meta, 
                                                   permutations = permutations)

Betadis.All.PotablePOCPOU.bray.Region<-anova(betadisper(All.PotablePOCPOU.dist.bray,
                                                            All.PotablePOCPOU.meta$Region))



ANOSIM.Pval.All.PotablePOCPOU.bray.Region<-ANOSIM.All.PotablePOCPOU.bray.Region$signif
ANOSIM.Pval.All.PotablePOCPOU.bray.Region

ADONIS2.Pval.All.PotablePOCPOU.bray.Region<-ADONIS2.All.PotablePOCPOU.bray.Region$`Pr(>F)`[1]
ADONIS2.Pval.All.PotablePOCPOU.bray.Region  

Betadis.Pval.All.PotablePOCPOU.bray.Region<-Betadis.All.PotablePOCPOU.bray.Region$`Pr(>F)`[1]
Betadis.Pval.All.PotablePOCPOU.bray.Region

ANOSIM.R.All.PotablePOCPOU.bray.Region<-ANOSIM.All.PotablePOCPOU.bray.Region$statistic
ANOSIM.R.All.PotablePOCPOU.bray.Region

ADONIS2.R2.All.PotablePOCPOU.bray.Region<-ADONIS2.All.PotablePOCPOU.bray.Region$R2[1]
ADONIS2.R2.All.PotablePOCPOU.bray.Region  

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.PotablePOCPOU.bray["Region"]<-ANOSIM.Pval.All.PotablePOCPOU.bray.Region
ADONIS2.GroupedP.All.PotablePOCPOU.bray["Region"]<-ADONIS2.Pval.All.PotablePOCPOU.bray.Region

Betadis.GroupedP.All.PotablePOCPOU.bray["Region"]<-Betadis.Pval.All.PotablePOCPOU.bray.Region

ANOSIM.GroupedR.All.PotablePOCPOU.bray["Region"]<-ANOSIM.R.All.PotablePOCPOU.bray.Region
ADONIS2.GroupedR.All.PotablePOCPOU.bray["Region"]<-ADONIS2.R2.All.PotablePOCPOU.bray.Region

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.PotablePOCPOU.meta$Region)

ADONIS2.Pairwise.All.PotablePOCPOU.bray.Region<-pairwise.adonis2(All.PotablePOCPOU.dist.bray ~ Region, 
                                                                     data = All.PotablePOCPOU.meta, 
                                                                     permutations = permutations)
ADONIS2.Pairwise.All.PotablePOCPOU.bray.Region
#Region





#################
#################
#################
#################
#Sample_Local 
ANOSIM.All.PotablePOCPOU.bray.Sample_Local<-anosim(All.PotablePOCPOU.dist.bray,
                                                 All.PotablePOCPOU.meta$Sample_Local, 
                                                 permutations = permutations)

ADONIS2.All.PotablePOCPOU.bray.Sample_Local<-adonis2(All.PotablePOCPOU.dist.bray ~ Sample_Local, 
                                                   data = All.PotablePOCPOU.meta, 
                                                   permutations = permutations)

Betadis.All.PotablePOCPOU.bray.Sample_Local<-anova(betadisper(All.PotablePOCPOU.dist.bray,
                                                            All.PotablePOCPOU.meta$Sample_Local))



ANOSIM.Pval.All.PotablePOCPOU.bray.Sample_Local<-ANOSIM.All.PotablePOCPOU.bray.Sample_Local$signif
ANOSIM.Pval.All.PotablePOCPOU.bray.Sample_Local

ADONIS2.Pval.All.PotablePOCPOU.bray.Sample_Local<-ADONIS2.All.PotablePOCPOU.bray.Sample_Local$`Pr(>F)`[1]
ADONIS2.Pval.All.PotablePOCPOU.bray.Sample_Local  

Betadis.Pval.All.PotablePOCPOU.bray.Sample_Local<-Betadis.All.PotablePOCPOU.bray.Sample_Local$`Pr(>F)`[1]
Betadis.Pval.All.PotablePOCPOU.bray.Sample_Local

ANOSIM.R.All.PotablePOCPOU.bray.Sample_Local<-ANOSIM.All.PotablePOCPOU.bray.Sample_Local$statistic
ANOSIM.R.All.PotablePOCPOU.bray.Sample_Local

ADONIS2.R2.All.PotablePOCPOU.bray.Sample_Local<-ADONIS2.All.PotablePOCPOU.bray.Sample_Local$R2[1]
ADONIS2.R2.All.PotablePOCPOU.bray.Sample_Local  

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.PotablePOCPOU.bray["Sample_Local"]<-ANOSIM.Pval.All.PotablePOCPOU.bray.Sample_Local
ADONIS2.GroupedP.All.PotablePOCPOU.bray["Sample_Local"]<-ADONIS2.Pval.All.PotablePOCPOU.bray.Sample_Local

Betadis.GroupedP.All.PotablePOCPOU.bray["Sample_Local"]<-Betadis.Pval.All.PotablePOCPOU.bray.Sample_Local

ANOSIM.GroupedR.All.PotablePOCPOU.bray["Sample_Local"]<-ANOSIM.R.All.PotablePOCPOU.bray.Sample_Local
ADONIS2.GroupedR.All.PotablePOCPOU.bray["Sample_Local"]<-ADONIS2.R2.All.PotablePOCPOU.bray.Sample_Local

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.PotablePOCPOU.meta$Sample_Local)

ADONIS2.Pairwise.All.PotablePOCPOU.bray.Sample_Local<-pairwise.adonis2(All.PotablePOCPOU.dist.bray ~ Sample_Local, 
                                                                     data = All.PotablePOCPOU.meta, 
                                                                     permutations = permutations)
ADONIS2.Pairwise.All.PotablePOCPOU.bray.Sample_Local
#Sample_Local





#################
#################
#################
#################
#Final_Disinfection_Residual 
ANOSIM.All.PotablePOCPOU.bray.Final_Disinfection_Residual<-anosim(All.PotablePOCPOU.dist.bray,
                                                 All.PotablePOCPOU.meta$Final_Disinfection_Residual, 
                                                 permutations = permutations)

ADONIS2.All.PotablePOCPOU.bray.Final_Disinfection_Residual<-adonis2(All.PotablePOCPOU.dist.bray ~ Final_Disinfection_Residual, 
                                                   data = All.PotablePOCPOU.meta, 
                                                   permutations = permutations)

Betadis.All.PotablePOCPOU.bray.Final_Disinfection_Residual<-anova(betadisper(All.PotablePOCPOU.dist.bray,
                                                            All.PotablePOCPOU.meta$Final_Disinfection_Residual))



ANOSIM.Pval.All.PotablePOCPOU.bray.Final_Disinfection_Residual<-ANOSIM.All.PotablePOCPOU.bray.Final_Disinfection_Residual$signif
ANOSIM.Pval.All.PotablePOCPOU.bray.Final_Disinfection_Residual

ADONIS2.Pval.All.PotablePOCPOU.bray.Final_Disinfection_Residual<-ADONIS2.All.PotablePOCPOU.bray.Final_Disinfection_Residual$`Pr(>F)`[1]
ADONIS2.Pval.All.PotablePOCPOU.bray.Final_Disinfection_Residual  

Betadis.Pval.All.PotablePOCPOU.bray.Final_Disinfection_Residual<-Betadis.All.PotablePOCPOU.bray.Final_Disinfection_Residual$`Pr(>F)`[1]
Betadis.Pval.All.PotablePOCPOU.bray.Final_Disinfection_Residual

ANOSIM.R.All.PotablePOCPOU.bray.Final_Disinfection_Residual<-ANOSIM.All.PotablePOCPOU.bray.Final_Disinfection_Residual$statistic
ANOSIM.R.All.PotablePOCPOU.bray.Final_Disinfection_Residual

ADONIS2.R2.All.PotablePOCPOU.bray.Final_Disinfection_Residual<-ADONIS2.All.PotablePOCPOU.bray.Final_Disinfection_Residual$R2[1]
ADONIS2.R2.All.PotablePOCPOU.bray.Final_Disinfection_Residual  

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.PotablePOCPOU.bray["Final_Disinfection_Residual"]<-ANOSIM.Pval.All.PotablePOCPOU.bray.Final_Disinfection_Residual
ADONIS2.GroupedP.All.PotablePOCPOU.bray["Final_Disinfection_Residual"]<-ADONIS2.Pval.All.PotablePOCPOU.bray.Final_Disinfection_Residual

Betadis.GroupedP.All.PotablePOCPOU.bray["Final_Disinfection_Residual"]<-Betadis.Pval.All.PotablePOCPOU.bray.Final_Disinfection_Residual

ANOSIM.GroupedR.All.PotablePOCPOU.bray["Final_Disinfection_Residual"]<-ANOSIM.R.All.PotablePOCPOU.bray.Final_Disinfection_Residual
ADONIS2.GroupedR.All.PotablePOCPOU.bray["Final_Disinfection_Residual"]<-ADONIS2.R2.All.PotablePOCPOU.bray.Final_Disinfection_Residual

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.PotablePOCPOU.meta$Final_Disinfection_Residual)

ADONIS2.Pairwise.All.PotablePOCPOU.bray.Final_Disinfection_Residual<-pairwise.adonis2(All.PotablePOCPOU.dist.bray ~ Final_Disinfection_Residual, 
                                                                     data = All.PotablePOCPOU.meta, 
                                                                     permutations = permutations)
ADONIS2.Pairwise.All.PotablePOCPOU.bray.Final_Disinfection_Residual
#Final_Disinfection_Residual





# #################
# #################
# #################
# #################
# #Primers 
# ANOSIM.All.PotablePOCPOU.bray.Primers<-anosim(All.PotablePOCPOU.dist.bray,
#                                                  All.PotablePOCPOU.meta$Primers, 
#                                                  permutations = permutations)
# 
# ADONIS2.All.PotablePOCPOU.bray.Primers<-adonis2(All.PotablePOCPOU.dist.bray ~ Primers, 
#                                                    data = All.PotablePOCPOU.meta, 
#                                                    permutations = permutations)
# 
# Betadis.All.PotablePOCPOU.bray.Primers<-anova(betadisper(All.PotablePOCPOU.dist.bray,
#                                                             All.PotablePOCPOU.meta$Primers))
# 
# 
# 
# ANOSIM.Pval.All.PotablePOCPOU.bray.Primers<-ANOSIM.All.PotablePOCPOU.bray.Primers$signif
# ANOSIM.Pval.All.PotablePOCPOU.bray.Primers
# 
# ADONIS2.Pval.All.PotablePOCPOU.bray.Primers<-ADONIS2.All.PotablePOCPOU.bray.Primers$`Pr(>F)`[1]
# ADONIS2.Pval.All.PotablePOCPOU.bray.Primers  
# 
# Betadis.Pval.All.PotablePOCPOU.bray.Primers<-Betadis.All.PotablePOCPOU.bray.Primers$`Pr(>F)`[1]
# Betadis.Pval.All.PotablePOCPOU.bray.Primers
# 
# ANOSIM.R.All.PotablePOCPOU.bray.Primers<-ANOSIM.All.PotablePOCPOU.bray.Primers$statistic
# ANOSIM.R.All.PotablePOCPOU.bray.Primers
# 
# ADONIS2.R2.All.PotablePOCPOU.bray.Primers<-ADONIS2.All.PotablePOCPOU.bray.Primers$R2[1]
# ADONIS2.R2.All.PotablePOCPOU.bray.Primers  
# 
# #Pull Together Pvalues to adjust 
# ANOSIM.GroupedP.All.PotablePOCPOU.bray["Primers"]<-ANOSIM.Pval.All.PotablePOCPOU.bray.Primers
# ADONIS2.GroupedP.All.PotablePOCPOU.bray["Primers"]<-ADONIS2.Pval.All.PotablePOCPOU.bray.Primers
# 
# Betadis.GroupedP.All.PotablePOCPOU.bray["Primers"]<-Betadis.Pval.All.PotablePOCPOU.bray.Primers
# 
# ANOSIM.GroupedR.All.PotablePOCPOU.bray["Primers"]<-ANOSIM.R.All.PotablePOCPOU.bray.Primers
# ADONIS2.GroupedR.All.PotablePOCPOU.bray["Primers"]<-ADONIS2.R2.All.PotablePOCPOU.bray.Primers
# 
# #Pairwise if needed, check if interesting or not. Remove if not. 
# count(All.PotablePOCPOU.meta$Primers)
# 
# ADONIS2.Pairwise.All.PotablePOCPOU.bray.Primers<-pairwise.adonis2(All.PotablePOCPOU.dist.bray ~ Primers, 
#                                                                      data = All.PotablePOCPOU.meta, 
#                                                                      permutations = permutations)
# ADONIS2.Pairwise.All.PotablePOCPOU.bray.Primers
# #Primers
# 
# 
# 
# 
# 
# #################
# #################
# #################
# #################
# #Trimmed 
# ANOSIM.All.PotablePOCPOU.bray.Trimmed<-anosim(All.PotablePOCPOU.dist.bray,
#                                                  All.PotablePOCPOU.meta$Trimmed, 
#                                                  permutations = permutations)
# 
# ADONIS2.All.PotablePOCPOU.bray.Trimmed<-adonis2(All.PotablePOCPOU.dist.bray ~ Trimmed, 
#                                                    data = All.PotablePOCPOU.meta, 
#                                                    permutations = permutations)
# 
# Betadis.All.PotablePOCPOU.bray.Trimmed<-anova(betadisper(All.PotablePOCPOU.dist.bray,
#                                                             All.PotablePOCPOU.meta$Trimmed))
# 
# 
# 
# ANOSIM.Pval.All.PotablePOCPOU.bray.Trimmed<-ANOSIM.All.PotablePOCPOU.bray.Trimmed$signif
# ANOSIM.Pval.All.PotablePOCPOU.bray.Trimmed
# 
# ADONIS2.Pval.All.PotablePOCPOU.bray.Trimmed<-ADONIS2.All.PotablePOCPOU.bray.Trimmed$`Pr(>F)`[1]
# ADONIS2.Pval.All.PotablePOCPOU.bray.Trimmed  
# 
# Betadis.Pval.All.PotablePOCPOU.bray.Trimmed<-Betadis.All.PotablePOCPOU.bray.Trimmed$`Pr(>F)`[1]
# Betadis.Pval.All.PotablePOCPOU.bray.Trimmed
# 
# ANOSIM.R.All.PotablePOCPOU.bray.Trimmed<-ANOSIM.All.PotablePOCPOU.bray.Trimmed$statistic
# ANOSIM.R.All.PotablePOCPOU.bray.Trimmed
# 
# ADONIS2.R2.All.PotablePOCPOU.bray.Trimmed<-ADONIS2.All.PotablePOCPOU.bray.Trimmed$R2[1]
# ADONIS2.R2.All.PotablePOCPOU.bray.Trimmed  
# 
# #Pull Together Pvalues to adjust 
# ANOSIM.GroupedP.All.PotablePOCPOU.bray["Trimmed"]<-ANOSIM.Pval.All.PotablePOCPOU.bray.Trimmed
# ADONIS2.GroupedP.All.PotablePOCPOU.bray["Trimmed"]<-ADONIS2.Pval.All.PotablePOCPOU.bray.Trimmed
# 
# Betadis.GroupedP.All.PotablePOCPOU.bray["Trimmed"]<-Betadis.Pval.All.PotablePOCPOU.bray.Trimmed
# 
# ANOSIM.GroupedR.All.PotablePOCPOU.bray["Trimmed"]<-ANOSIM.R.All.PotablePOCPOU.bray.Trimmed
# ADONIS2.GroupedR.All.PotablePOCPOU.bray["Trimmed"]<-ADONIS2.R2.All.PotablePOCPOU.bray.Trimmed
# 
# #Pairwise if needed, check if interesting or not. Remove if not. 
# count(All.PotablePOCPOU.meta$Trimmed)
# 
# ADONIS2.Pairwise.All.PotablePOCPOU.bray.Trimmed<-pairwise.adonis2(All.PotablePOCPOU.dist.bray ~ Trimmed, 
#                                                                      data = All.PotablePOCPOU.meta, 
#                                                                      permutations = permutations)
# ADONIS2.Pairwise.All.PotablePOCPOU.bray.Trimmed
# #Trimmed
# 
# 
# 
# 
# 
# #################
# #################
# #################
# #################
# #Extraction_Type 
# ANOSIM.All.PotablePOCPOU.bray.Extraction_Type<-anosim(All.PotablePOCPOU.dist.bray,
#                                                  All.PotablePOCPOU.meta$Extraction_Type, 
#                                                  permutations = permutations)
# 
# ADONIS2.All.PotablePOCPOU.bray.Extraction_Type<-adonis2(All.PotablePOCPOU.dist.bray ~ Extraction_Type, 
#                                                    data = All.PotablePOCPOU.meta, 
#                                                    permutations = permutations)
# 
# Betadis.All.PotablePOCPOU.bray.Extraction_Type<-anova(betadisper(All.PotablePOCPOU.dist.bray,
#                                                             All.PotablePOCPOU.meta$Extraction_Type))
# 
# 
# 
# ANOSIM.Pval.All.PotablePOCPOU.bray.Extraction_Type<-ANOSIM.All.PotablePOCPOU.bray.Extraction_Type$signif
# ANOSIM.Pval.All.PotablePOCPOU.bray.Extraction_Type
# 
# ADONIS2.Pval.All.PotablePOCPOU.bray.Extraction_Type<-ADONIS2.All.PotablePOCPOU.bray.Extraction_Type$`Pr(>F)`[1]
# ADONIS2.Pval.All.PotablePOCPOU.bray.Extraction_Type  
# 
# Betadis.Pval.All.PotablePOCPOU.bray.Extraction_Type<-Betadis.All.PotablePOCPOU.bray.Extraction_Type$`Pr(>F)`[1]
# Betadis.Pval.All.PotablePOCPOU.bray.Extraction_Type
# 
# ANOSIM.R.All.PotablePOCPOU.bray.Extraction_Type<-ANOSIM.All.PotablePOCPOU.bray.Extraction_Type$statistic
# ANOSIM.R.All.PotablePOCPOU.bray.Extraction_Type
# 
# ADONIS2.R2.All.PotablePOCPOU.bray.Extraction_Type<-ADONIS2.All.PotablePOCPOU.bray.Extraction_Type$R2[1]
# ADONIS2.R2.All.PotablePOCPOU.bray.Extraction_Type  
# 
# #Pull Together Pvalues to adjust 
# ANOSIM.GroupedP.All.PotablePOCPOU.bray["Extraction_Type"]<-ANOSIM.Pval.All.PotablePOCPOU.bray.Extraction_Type
# ADONIS2.GroupedP.All.PotablePOCPOU.bray["Extraction_Type"]<-ADONIS2.Pval.All.PotablePOCPOU.bray.Extraction_Type
# 
# Betadis.GroupedP.All.PotablePOCPOU.bray["Extraction_Type"]<-Betadis.Pval.All.PotablePOCPOU.bray.Extraction_Type
# 
# ANOSIM.GroupedR.All.PotablePOCPOU.bray["Extraction_Type"]<-ANOSIM.R.All.PotablePOCPOU.bray.Extraction_Type
# ADONIS2.GroupedR.All.PotablePOCPOU.bray["Extraction_Type"]<-ADONIS2.R2.All.PotablePOCPOU.bray.Extraction_Type
# 
# #Pairwise if needed, check if interesting or not. Remove if not. 
# count(All.PotablePOCPOU.meta$Extraction_Type)
# 
# ADONIS2.Pairwise.All.PotablePOCPOU.bray.Extraction_Type<-pairwise.adonis2(All.PotablePOCPOU.dist.bray ~ Extraction_Type, 
#                                                                      data = All.PotablePOCPOU.meta, 
#                                                                      permutations = permutations)
# ADONIS2.Pairwise.All.PotablePOCPOU.bray.Extraction_Type
# #Extraction_Type





#################
#################
#################
#################
#Potable_Source_Water 
ANOSIM.All.PotablePOCPOU.bray.Potable_Source_Water<-anosim(All.PotablePOCPOU.dist.bray,
                                                      All.PotablePOCPOU.meta$Potable_Source_Water, 
                                                      permutations = permutations)

ADONIS2.All.PotablePOCPOU.bray.Potable_Source_Water<-adonis2(All.PotablePOCPOU.dist.bray ~ Potable_Source_Water, 
                                                        data = All.PotablePOCPOU.meta, 
                                                        permutations = permutations)

Betadis.All.PotablePOCPOU.bray.Potable_Source_Water<-anova(betadisper(All.PotablePOCPOU.dist.bray,
                                                                 All.PotablePOCPOU.meta$Potable_Source_Water))



ANOSIM.Pval.All.PotablePOCPOU.bray.Potable_Source_Water<-ANOSIM.All.PotablePOCPOU.bray.Potable_Source_Water$signif
ANOSIM.Pval.All.PotablePOCPOU.bray.Potable_Source_Water

ADONIS2.Pval.All.PotablePOCPOU.bray.Potable_Source_Water<-ADONIS2.All.PotablePOCPOU.bray.Potable_Source_Water$`Pr(>F)`[1]
ADONIS2.Pval.All.PotablePOCPOU.bray.Potable_Source_Water  

Betadis.Pval.All.PotablePOCPOU.bray.Potable_Source_Water<-Betadis.All.PotablePOCPOU.bray.Potable_Source_Water$`Pr(>F)`[1]
Betadis.Pval.All.PotablePOCPOU.bray.Potable_Source_Water

ANOSIM.R.All.PotablePOCPOU.bray.Potable_Source_Water<-ANOSIM.All.PotablePOCPOU.bray.Potable_Source_Water$statistic
ANOSIM.R.All.PotablePOCPOU.bray.Potable_Source_Water

ADONIS2.R2.All.PotablePOCPOU.bray.Potable_Source_Water<-ADONIS2.All.PotablePOCPOU.bray.Potable_Source_Water$R2[1]
ADONIS2.R2.All.PotablePOCPOU.bray.Potable_Source_Water  

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.PotablePOCPOU.bray["Potable_Source_Water"]<-ANOSIM.Pval.All.PotablePOCPOU.bray.Potable_Source_Water
ADONIS2.GroupedP.All.PotablePOCPOU.bray["Potable_Source_Water"]<-ADONIS2.Pval.All.PotablePOCPOU.bray.Potable_Source_Water

Betadis.GroupedP.All.PotablePOCPOU.bray["Potable_Source_Water"]<-Betadis.Pval.All.PotablePOCPOU.bray.Potable_Source_Water

ANOSIM.GroupedR.All.PotablePOCPOU.bray["Potable_Source_Water"]<-ANOSIM.R.All.PotablePOCPOU.bray.Potable_Source_Water
ADONIS2.GroupedR.All.PotablePOCPOU.bray["Potable_Source_Water"]<-ADONIS2.R2.All.PotablePOCPOU.bray.Potable_Source_Water

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.PotablePOCPOU.meta$Potable_Source_Water)

ADONIS2.Pairwise.All.PotablePOCPOU.bray.Potable_Source_Water<-pairwise.adonis2(All.PotablePOCPOU.dist.bray ~ Potable_Source_Water, 
                                                                          data = All.PotablePOCPOU.meta, 
                                                                          permutations = permutations)
ADONIS2.Pairwise.All.PotablePOCPOU.bray.Potable_Source_Water
#Potable_Source_Water
##Potable_Source_Water




#################
#################
#################
#################
#Potable_Disinfection_Primary  
ANOSIM.All.PotablePOCPOU.bray.Potable_Disinfection_Primary <-anosim(All.PotablePOCPOU.dist.bray,
                                                      All.PotablePOCPOU.meta$Potable_Disinfection_Primary , 
                                                      permutations = permutations)

ADONIS2.All.PotablePOCPOU.bray.Potable_Disinfection_Primary <-adonis2(All.PotablePOCPOU.dist.bray ~ Potable_Disinfection_Primary , 
                                                        data = All.PotablePOCPOU.meta, 
                                                        permutations = permutations)

Betadis.All.PotablePOCPOU.bray.Potable_Disinfection_Primary <-anova(betadisper(All.PotablePOCPOU.dist.bray,
                                                                 All.PotablePOCPOU.meta$Potable_Disinfection_Primary ))



ANOSIM.Pval.All.PotablePOCPOU.bray.Potable_Disinfection_Primary <-ANOSIM.All.PotablePOCPOU.bray.Potable_Disinfection_Primary $signif
ANOSIM.Pval.All.PotablePOCPOU.bray.Potable_Disinfection_Primary 

ADONIS2.Pval.All.PotablePOCPOU.bray.Potable_Disinfection_Primary <-ADONIS2.All.PotablePOCPOU.bray.Potable_Disinfection_Primary $`Pr(>F)`[1]
ADONIS2.Pval.All.PotablePOCPOU.bray.Potable_Disinfection_Primary   

Betadis.Pval.All.PotablePOCPOU.bray.Potable_Disinfection_Primary <-Betadis.All.PotablePOCPOU.bray.Potable_Disinfection_Primary $`Pr(>F)`[1]
Betadis.Pval.All.PotablePOCPOU.bray.Potable_Disinfection_Primary 

ANOSIM.R.All.PotablePOCPOU.bray.Potable_Disinfection_Primary <-ANOSIM.All.PotablePOCPOU.bray.Potable_Disinfection_Primary $statistic
ANOSIM.R.All.PotablePOCPOU.bray.Potable_Disinfection_Primary 

ADONIS2.R2.All.PotablePOCPOU.bray.Potable_Disinfection_Primary <-ADONIS2.All.PotablePOCPOU.bray.Potable_Disinfection_Primary $R2[1]
ADONIS2.R2.All.PotablePOCPOU.bray.Potable_Disinfection_Primary   

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.PotablePOCPOU.bray["Potable_Disinfection_Primary "]<-ANOSIM.Pval.All.PotablePOCPOU.bray.Potable_Disinfection_Primary 
ADONIS2.GroupedP.All.PotablePOCPOU.bray["Potable_Disinfection_Primary "]<-ADONIS2.Pval.All.PotablePOCPOU.bray.Potable_Disinfection_Primary 

Betadis.GroupedP.All.PotablePOCPOU.bray["Potable_Disinfection_Primary "]<-Betadis.Pval.All.PotablePOCPOU.bray.Potable_Disinfection_Primary 

ANOSIM.GroupedR.All.PotablePOCPOU.bray["Potable_Disinfection_Primary "]<-ANOSIM.R.All.PotablePOCPOU.bray.Potable_Disinfection_Primary 
ADONIS2.GroupedR.All.PotablePOCPOU.bray["Potable_Disinfection_Primary "]<-ADONIS2.R2.All.PotablePOCPOU.bray.Potable_Disinfection_Primary 

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.PotablePOCPOU.meta$Potable_Disinfection_Primary )

ADONIS2.Pairwise.All.PotablePOCPOU.bray.Potable_Disinfection_Primary <-pairwise.adonis2(All.PotablePOCPOU.dist.bray ~ Potable_Disinfection_Primary , 
                                                                          data = All.PotablePOCPOU.meta, 
                                                                          permutations = permutations)
ADONIS2.Pairwise.All.PotablePOCPOU.bray.Potable_Disinfection_Primary 
#Potable_Disinfection_Primary 
##Potable_Disinfection_Primary 




#################
#################
#################
#################
#Potable_Treatment_Type 
ANOSIM.All.PotablePOCPOU.bray.Potable_Treatment_Type<-anosim(All.PotablePOCPOU.dist.bray,
                                                      All.PotablePOCPOU.meta$Potable_Treatment_Type, 
                                                      permutations = permutations)

ADONIS2.All.PotablePOCPOU.bray.Potable_Treatment_Type<-adonis2(All.PotablePOCPOU.dist.bray ~ Potable_Treatment_Type, 
                                                        data = All.PotablePOCPOU.meta, 
                                                        permutations = permutations)

Betadis.All.PotablePOCPOU.bray.Potable_Treatment_Type<-anova(betadisper(All.PotablePOCPOU.dist.bray,
                                                                 All.PotablePOCPOU.meta$Potable_Treatment_Type))



ANOSIM.Pval.All.PotablePOCPOU.bray.Potable_Treatment_Type<-ANOSIM.All.PotablePOCPOU.bray.Potable_Treatment_Type$signif
ANOSIM.Pval.All.PotablePOCPOU.bray.Potable_Treatment_Type

ADONIS2.Pval.All.PotablePOCPOU.bray.Potable_Treatment_Type<-ADONIS2.All.PotablePOCPOU.bray.Potable_Treatment_Type$`Pr(>F)`[1]
ADONIS2.Pval.All.PotablePOCPOU.bray.Potable_Treatment_Type  

Betadis.Pval.All.PotablePOCPOU.bray.Potable_Treatment_Type<-Betadis.All.PotablePOCPOU.bray.Potable_Treatment_Type$`Pr(>F)`[1]
Betadis.Pval.All.PotablePOCPOU.bray.Potable_Treatment_Type

ANOSIM.R.All.PotablePOCPOU.bray.Potable_Treatment_Type<-ANOSIM.All.PotablePOCPOU.bray.Potable_Treatment_Type$statistic
ANOSIM.R.All.PotablePOCPOU.bray.Potable_Treatment_Type

ADONIS2.R2.All.PotablePOCPOU.bray.Potable_Treatment_Type<-ADONIS2.All.PotablePOCPOU.bray.Potable_Treatment_Type$R2[1]
ADONIS2.R2.All.PotablePOCPOU.bray.Potable_Treatment_Type  

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.PotablePOCPOU.bray["Potable_Treatment_Type"]<-ANOSIM.Pval.All.PotablePOCPOU.bray.Potable_Treatment_Type
ADONIS2.GroupedP.All.PotablePOCPOU.bray["Potable_Treatment_Type"]<-ADONIS2.Pval.All.PotablePOCPOU.bray.Potable_Treatment_Type

Betadis.GroupedP.All.PotablePOCPOU.bray["Potable_Treatment_Type"]<-Betadis.Pval.All.PotablePOCPOU.bray.Potable_Treatment_Type

ANOSIM.GroupedR.All.PotablePOCPOU.bray["Potable_Treatment_Type"]<-ANOSIM.R.All.PotablePOCPOU.bray.Potable_Treatment_Type
ADONIS2.GroupedR.All.PotablePOCPOU.bray["Potable_Treatment_Type"]<-ADONIS2.R2.All.PotablePOCPOU.bray.Potable_Treatment_Type

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.PotablePOCPOU.meta$Potable_Treatment_Type)

ADONIS2.Pairwise.All.PotablePOCPOU.bray.Potable_Treatment_Type<-pairwise.adonis2(All.PotablePOCPOU.dist.bray ~ Potable_Treatment_Type, 
                                                                          data = All.PotablePOCPOU.meta, 
                                                                          permutations = permutations)
ADONIS2.Pairwise.All.PotablePOCPOU.bray.Potable_Treatment_Type
#Potable_Treatment_Type
##Potable_Treatment_Type 




#################
#################
#################
#################
#Disinfection_Total 
ANOSIM.All.PotablePOCPOU.bray.Disinfection_Total<-anosim(All.PotablePOCPOU.dist.bray,
                                                      All.PotablePOCPOU.meta$Disinfection_Total, 
                                                      permutations = permutations)

ADONIS2.All.PotablePOCPOU.bray.Disinfection_Total<-adonis2(All.PotablePOCPOU.dist.bray ~ Disinfection_Total, 
                                                        data = All.PotablePOCPOU.meta, 
                                                        permutations = permutations)

Betadis.All.PotablePOCPOU.bray.Disinfection_Total<-anova(betadisper(All.PotablePOCPOU.dist.bray,
                                                                 All.PotablePOCPOU.meta$Disinfection_Total))



ANOSIM.Pval.All.PotablePOCPOU.bray.Disinfection_Total<-ANOSIM.All.PotablePOCPOU.bray.Disinfection_Total$signif
ANOSIM.Pval.All.PotablePOCPOU.bray.Disinfection_Total

ADONIS2.Pval.All.PotablePOCPOU.bray.Disinfection_Total<-ADONIS2.All.PotablePOCPOU.bray.Disinfection_Total$`Pr(>F)`[1]
ADONIS2.Pval.All.PotablePOCPOU.bray.Disinfection_Total  

Betadis.Pval.All.PotablePOCPOU.bray.Disinfection_Total<-Betadis.All.PotablePOCPOU.bray.Disinfection_Total$`Pr(>F)`[1]
Betadis.Pval.All.PotablePOCPOU.bray.Disinfection_Total

ANOSIM.R.All.PotablePOCPOU.bray.Disinfection_Total<-ANOSIM.All.PotablePOCPOU.bray.Disinfection_Total$statistic
ANOSIM.R.All.PotablePOCPOU.bray.Disinfection_Total

ADONIS2.R2.All.PotablePOCPOU.bray.Disinfection_Total<-ADONIS2.All.PotablePOCPOU.bray.Disinfection_Total$R2[1]
ADONIS2.R2.All.PotablePOCPOU.bray.Disinfection_Total  

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.PotablePOCPOU.bray["Disinfection_Total"]<-ANOSIM.Pval.All.PotablePOCPOU.bray.Disinfection_Total
ADONIS2.GroupedP.All.PotablePOCPOU.bray["Disinfection_Total"]<-ADONIS2.Pval.All.PotablePOCPOU.bray.Disinfection_Total

Betadis.GroupedP.All.PotablePOCPOU.bray["Disinfection_Total"]<-Betadis.Pval.All.PotablePOCPOU.bray.Disinfection_Total

ANOSIM.GroupedR.All.PotablePOCPOU.bray["Disinfection_Total"]<-ANOSIM.R.All.PotablePOCPOU.bray.Disinfection_Total
ADONIS2.GroupedR.All.PotablePOCPOU.bray["Disinfection_Total"]<-ADONIS2.R2.All.PotablePOCPOU.bray.Disinfection_Total

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.PotablePOCPOU.meta$Disinfection_Total)

ADONIS2.Pairwise.All.PotablePOCPOU.bray.Disinfection_Total<-pairwise.adonis2(All.PotablePOCPOU.dist.bray ~ Disinfection_Total, 
                                                                          data = All.PotablePOCPOU.meta, 
                                                                          permutations = permutations)
ADONIS2.Pairwise.All.PotablePOCPOU.bray.Disinfection_Total
#Disinfection_Total
##Disinfection_Total 




#################
#################
#################
#################
#Treatment_Total  
ANOSIM.All.PotablePOCPOU.bray.Treatment_Total <-anosim(All.PotablePOCPOU.dist.bray,
                                                      All.PotablePOCPOU.meta$Treatment_Total , 
                                                      permutations = permutations)

ADONIS2.All.PotablePOCPOU.bray.Treatment_Total <-adonis2(All.PotablePOCPOU.dist.bray ~ Treatment_Total , 
                                                        data = All.PotablePOCPOU.meta, 
                                                        permutations = permutations)

Betadis.All.PotablePOCPOU.bray.Treatment_Total <-anova(betadisper(All.PotablePOCPOU.dist.bray,
                                                                 All.PotablePOCPOU.meta$Treatment_Total ))



ANOSIM.Pval.All.PotablePOCPOU.bray.Treatment_Total <-ANOSIM.All.PotablePOCPOU.bray.Treatment_Total $signif
ANOSIM.Pval.All.PotablePOCPOU.bray.Treatment_Total 

ADONIS2.Pval.All.PotablePOCPOU.bray.Treatment_Total <-ADONIS2.All.PotablePOCPOU.bray.Treatment_Total $`Pr(>F)`[1]
ADONIS2.Pval.All.PotablePOCPOU.bray.Treatment_Total   

Betadis.Pval.All.PotablePOCPOU.bray.Treatment_Total <-Betadis.All.PotablePOCPOU.bray.Treatment_Total $`Pr(>F)`[1]
Betadis.Pval.All.PotablePOCPOU.bray.Treatment_Total 

ANOSIM.R.All.PotablePOCPOU.bray.Treatment_Total <-ANOSIM.All.PotablePOCPOU.bray.Treatment_Total $statistic
ANOSIM.R.All.PotablePOCPOU.bray.Treatment_Total 

ADONIS2.R2.All.PotablePOCPOU.bray.Treatment_Total <-ADONIS2.All.PotablePOCPOU.bray.Treatment_Total $R2[1]
ADONIS2.R2.All.PotablePOCPOU.bray.Treatment_Total   

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.PotablePOCPOU.bray["Treatment_Total "]<-ANOSIM.Pval.All.PotablePOCPOU.bray.Treatment_Total 
ADONIS2.GroupedP.All.PotablePOCPOU.bray["Treatment_Total "]<-ADONIS2.Pval.All.PotablePOCPOU.bray.Treatment_Total 

Betadis.GroupedP.All.PotablePOCPOU.bray["Treatment_Total "]<-Betadis.Pval.All.PotablePOCPOU.bray.Treatment_Total 

ANOSIM.GroupedR.All.PotablePOCPOU.bray["Treatment_Total "]<-ANOSIM.R.All.PotablePOCPOU.bray.Treatment_Total 
ADONIS2.GroupedR.All.PotablePOCPOU.bray["Treatment_Total "]<-ADONIS2.R2.All.PotablePOCPOU.bray.Treatment_Total 

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.PotablePOCPOU.meta$Treatment_Total )

ADONIS2.Pairwise.All.PotablePOCPOU.bray.Treatment_Total <-pairwise.adonis2(All.PotablePOCPOU.dist.bray ~ Treatment_Total , 
                                                                          data = All.PotablePOCPOU.meta, 
                                                                          permutations = permutations)
ADONIS2.Pairwise.All.PotablePOCPOU.bray.Treatment_Total 
#Treatment_Total 
##Treatment_Total 




#################
#################
#################
#################
#Treatment_Classification_Number 
ANOSIM.All.PotablePOCPOU.bray.Treatment_Classification_Number<-anosim(All.PotablePOCPOU.dist.bray,
                                                      All.PotablePOCPOU.meta$Treatment_Classification_Number, 
                                                      permutations = permutations)

ADONIS2.All.PotablePOCPOU.bray.Treatment_Classification_Number<-adonis2(All.PotablePOCPOU.dist.bray ~ Treatment_Classification_Number, 
                                                        data = All.PotablePOCPOU.meta, 
                                                        permutations = permutations)

Betadis.All.PotablePOCPOU.bray.Treatment_Classification_Number<-anova(betadisper(All.PotablePOCPOU.dist.bray,
                                                                 All.PotablePOCPOU.meta$Treatment_Classification_Number))



ANOSIM.Pval.All.PotablePOCPOU.bray.Treatment_Classification_Number<-ANOSIM.All.PotablePOCPOU.bray.Treatment_Classification_Number$signif
ANOSIM.Pval.All.PotablePOCPOU.bray.Treatment_Classification_Number

ADONIS2.Pval.All.PotablePOCPOU.bray.Treatment_Classification_Number<-ADONIS2.All.PotablePOCPOU.bray.Treatment_Classification_Number$`Pr(>F)`[1]
ADONIS2.Pval.All.PotablePOCPOU.bray.Treatment_Classification_Number  

Betadis.Pval.All.PotablePOCPOU.bray.Treatment_Classification_Number<-Betadis.All.PotablePOCPOU.bray.Treatment_Classification_Number$`Pr(>F)`[1]
Betadis.Pval.All.PotablePOCPOU.bray.Treatment_Classification_Number

ANOSIM.R.All.PotablePOCPOU.bray.Treatment_Classification_Number<-ANOSIM.All.PotablePOCPOU.bray.Treatment_Classification_Number$statistic
ANOSIM.R.All.PotablePOCPOU.bray.Treatment_Classification_Number

ADONIS2.R2.All.PotablePOCPOU.bray.Treatment_Classification_Number<-ADONIS2.All.PotablePOCPOU.bray.Treatment_Classification_Number$R2[1]
ADONIS2.R2.All.PotablePOCPOU.bray.Treatment_Classification_Number  

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.PotablePOCPOU.bray["Treatment_Classification_Number"]<-ANOSIM.Pval.All.PotablePOCPOU.bray.Treatment_Classification_Number
ADONIS2.GroupedP.All.PotablePOCPOU.bray["Treatment_Classification_Number"]<-ADONIS2.Pval.All.PotablePOCPOU.bray.Treatment_Classification_Number

Betadis.GroupedP.All.PotablePOCPOU.bray["Treatment_Classification_Number"]<-Betadis.Pval.All.PotablePOCPOU.bray.Treatment_Classification_Number

ANOSIM.GroupedR.All.PotablePOCPOU.bray["Treatment_Classification_Number"]<-ANOSIM.R.All.PotablePOCPOU.bray.Treatment_Classification_Number
ADONIS2.GroupedR.All.PotablePOCPOU.bray["Treatment_Classification_Number"]<-ADONIS2.R2.All.PotablePOCPOU.bray.Treatment_Classification_Number

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.PotablePOCPOU.meta$Treatment_Classification_Number)

ADONIS2.Pairwise.All.PotablePOCPOU.bray.Treatment_Classification_Number<-pairwise.adonis2(All.PotablePOCPOU.dist.bray ~ Treatment_Classification_Number, 
                                                                          data = All.PotablePOCPOU.meta, 
                                                                          permutations = permutations)
ADONIS2.Pairwise.All.PotablePOCPOU.bray.Treatment_Classification_Number
#Treatment_Classification_Number
##Treatment_Classification_Number




#################
#################
#################
#################
#Treatment_Classification_Cat 
ANOSIM.All.PotablePOCPOU.bray.Treatment_Classification_Cat<-anosim(All.PotablePOCPOU.dist.bray,
                                                      All.PotablePOCPOU.meta$Treatment_Classification_Cat, 
                                                      permutations = permutations)

ADONIS2.All.PotablePOCPOU.bray.Treatment_Classification_Cat<-adonis2(All.PotablePOCPOU.dist.bray ~ Treatment_Classification_Cat, 
                                                        data = All.PotablePOCPOU.meta, 
                                                        permutations = permutations)

Betadis.All.PotablePOCPOU.bray.Treatment_Classification_Cat<-anova(betadisper(All.PotablePOCPOU.dist.bray,
                                                                 All.PotablePOCPOU.meta$Treatment_Classification_Cat))



ANOSIM.Pval.All.PotablePOCPOU.bray.Treatment_Classification_Cat<-ANOSIM.All.PotablePOCPOU.bray.Treatment_Classification_Cat$signif
ANOSIM.Pval.All.PotablePOCPOU.bray.Treatment_Classification_Cat

ADONIS2.Pval.All.PotablePOCPOU.bray.Treatment_Classification_Cat<-ADONIS2.All.PotablePOCPOU.bray.Treatment_Classification_Cat$`Pr(>F)`[1]
ADONIS2.Pval.All.PotablePOCPOU.bray.Treatment_Classification_Cat  

Betadis.Pval.All.PotablePOCPOU.bray.Treatment_Classification_Cat<-Betadis.All.PotablePOCPOU.bray.Treatment_Classification_Cat$`Pr(>F)`[1]
Betadis.Pval.All.PotablePOCPOU.bray.Treatment_Classification_Cat

ANOSIM.R.All.PotablePOCPOU.bray.Treatment_Classification_Cat<-ANOSIM.All.PotablePOCPOU.bray.Treatment_Classification_Cat$statistic
ANOSIM.R.All.PotablePOCPOU.bray.Treatment_Classification_Cat

ADONIS2.R2.All.PotablePOCPOU.bray.Treatment_Classification_Cat<-ADONIS2.All.PotablePOCPOU.bray.Treatment_Classification_Cat$R2[1]
ADONIS2.R2.All.PotablePOCPOU.bray.Treatment_Classification_Cat  

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.PotablePOCPOU.bray["Treatment_Classification_Cat"]<-ANOSIM.Pval.All.PotablePOCPOU.bray.Treatment_Classification_Cat
ADONIS2.GroupedP.All.PotablePOCPOU.bray["Treatment_Classification_Cat"]<-ADONIS2.Pval.All.PotablePOCPOU.bray.Treatment_Classification_Cat

Betadis.GroupedP.All.PotablePOCPOU.bray["Treatment_Classification_Cat"]<-Betadis.Pval.All.PotablePOCPOU.bray.Treatment_Classification_Cat

ANOSIM.GroupedR.All.PotablePOCPOU.bray["Treatment_Classification_Cat"]<-ANOSIM.R.All.PotablePOCPOU.bray.Treatment_Classification_Cat
ADONIS2.GroupedR.All.PotablePOCPOU.bray["Treatment_Classification_Cat"]<-ADONIS2.R2.All.PotablePOCPOU.bray.Treatment_Classification_Cat

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.PotablePOCPOU.meta$Treatment_Classification_Cat)

ADONIS2.Pairwise.All.PotablePOCPOU.bray.Treatment_Classification_Cat<-pairwise.adonis2(All.PotablePOCPOU.dist.bray ~ Treatment_Classification_Cat, 
                                                                          data = All.PotablePOCPOU.meta, 
                                                                          permutations = permutations)
ADONIS2.Pairwise.All.PotablePOCPOU.bray.Treatment_Classification_Cat
#Treatment_Classification_Cat
##Treatment_Classification_Cat




#################
#################
#################
#################
#Mixed_Source_Water 
ANOSIM.All.PotablePOCPOU.bray.Mixed_Source_Water<-anosim(All.PotablePOCPOU.dist.bray,
                                                      All.PotablePOCPOU.meta$Mixed_Source_Water, 
                                                      permutations = permutations)

ADONIS2.All.PotablePOCPOU.bray.Mixed_Source_Water<-adonis2(All.PotablePOCPOU.dist.bray ~ Mixed_Source_Water, 
                                                        data = All.PotablePOCPOU.meta, 
                                                        permutations = permutations)

Betadis.All.PotablePOCPOU.bray.Mixed_Source_Water<-anova(betadisper(All.PotablePOCPOU.dist.bray,
                                                                 All.PotablePOCPOU.meta$Mixed_Source_Water))



ANOSIM.Pval.All.PotablePOCPOU.bray.Mixed_Source_Water<-ANOSIM.All.PotablePOCPOU.bray.Mixed_Source_Water$signif
ANOSIM.Pval.All.PotablePOCPOU.bray.Mixed_Source_Water

ADONIS2.Pval.All.PotablePOCPOU.bray.Mixed_Source_Water<-ADONIS2.All.PotablePOCPOU.bray.Mixed_Source_Water$`Pr(>F)`[1]
ADONIS2.Pval.All.PotablePOCPOU.bray.Mixed_Source_Water  

Betadis.Pval.All.PotablePOCPOU.bray.Mixed_Source_Water<-Betadis.All.PotablePOCPOU.bray.Mixed_Source_Water$`Pr(>F)`[1]
Betadis.Pval.All.PotablePOCPOU.bray.Mixed_Source_Water

ANOSIM.R.All.PotablePOCPOU.bray.Mixed_Source_Water<-ANOSIM.All.PotablePOCPOU.bray.Mixed_Source_Water$statistic
ANOSIM.R.All.PotablePOCPOU.bray.Mixed_Source_Water

ADONIS2.R2.All.PotablePOCPOU.bray.Mixed_Source_Water<-ADONIS2.All.PotablePOCPOU.bray.Mixed_Source_Water$R2[1]
ADONIS2.R2.All.PotablePOCPOU.bray.Mixed_Source_Water  

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.PotablePOCPOU.bray["Mixed_Source_Water"]<-ANOSIM.Pval.All.PotablePOCPOU.bray.Mixed_Source_Water
ADONIS2.GroupedP.All.PotablePOCPOU.bray["Mixed_Source_Water"]<-ADONIS2.Pval.All.PotablePOCPOU.bray.Mixed_Source_Water

Betadis.GroupedP.All.PotablePOCPOU.bray["Mixed_Source_Water"]<-Betadis.Pval.All.PotablePOCPOU.bray.Mixed_Source_Water

ANOSIM.GroupedR.All.PotablePOCPOU.bray["Mixed_Source_Water"]<-ANOSIM.R.All.PotablePOCPOU.bray.Mixed_Source_Water
ADONIS2.GroupedR.All.PotablePOCPOU.bray["Mixed_Source_Water"]<-ADONIS2.R2.All.PotablePOCPOU.bray.Mixed_Source_Water

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.PotablePOCPOU.meta$Mixed_Source_Water)

ADONIS2.Pairwise.All.PotablePOCPOU.bray.Mixed_Source_Water<-pairwise.adonis2(All.PotablePOCPOU.dist.bray ~ Mixed_Source_Water, 
                                                                          data = All.PotablePOCPOU.meta, 
                                                                          permutations = permutations)
ADONIS2.Pairwise.All.PotablePOCPOU.bray.Mixed_Source_Water
#Mixed_Source_Water
##Mixed_Source_Water




#################
#################
#################
#################
#Potable_Overlap 
ANOSIM.All.PotablePOCPOU.bray.Potable_Overlap<-anosim(All.PotablePOCPOU.dist.bray,
                                                      All.PotablePOCPOU.meta$Potable_Overlap, 
                                                      permutations = permutations)

ADONIS2.All.PotablePOCPOU.bray.Potable_Overlap<-adonis2(All.PotablePOCPOU.dist.bray ~ Potable_Overlap, 
                                                        data = All.PotablePOCPOU.meta, 
                                                        permutations = permutations)

Betadis.All.PotablePOCPOU.bray.Potable_Overlap<-anova(betadisper(All.PotablePOCPOU.dist.bray,
                                                                 All.PotablePOCPOU.meta$Potable_Overlap))



ANOSIM.Pval.All.PotablePOCPOU.bray.Potable_Overlap<-ANOSIM.All.PotablePOCPOU.bray.Potable_Overlap$signif
ANOSIM.Pval.All.PotablePOCPOU.bray.Potable_Overlap

ADONIS2.Pval.All.PotablePOCPOU.bray.Potable_Overlap<-ADONIS2.All.PotablePOCPOU.bray.Potable_Overlap$`Pr(>F)`[1]
ADONIS2.Pval.All.PotablePOCPOU.bray.Potable_Overlap  

Betadis.Pval.All.PotablePOCPOU.bray.Potable_Overlap<-Betadis.All.PotablePOCPOU.bray.Potable_Overlap$`Pr(>F)`[1]
Betadis.Pval.All.PotablePOCPOU.bray.Potable_Overlap

ANOSIM.R.All.PotablePOCPOU.bray.Potable_Overlap<-ANOSIM.All.PotablePOCPOU.bray.Potable_Overlap$statistic
ANOSIM.R.All.PotablePOCPOU.bray.Potable_Overlap

ADONIS2.R2.All.PotablePOCPOU.bray.Potable_Overlap<-ADONIS2.All.PotablePOCPOU.bray.Potable_Overlap$R2[1]
ADONIS2.R2.All.PotablePOCPOU.bray.Potable_Overlap  

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.PotablePOCPOU.bray["Potable_Overlap"]<-ANOSIM.Pval.All.PotablePOCPOU.bray.Potable_Overlap
ADONIS2.GroupedP.All.PotablePOCPOU.bray["Potable_Overlap"]<-ADONIS2.Pval.All.PotablePOCPOU.bray.Potable_Overlap

Betadis.GroupedP.All.PotablePOCPOU.bray["Potable_Overlap"]<-Betadis.Pval.All.PotablePOCPOU.bray.Potable_Overlap

ANOSIM.GroupedR.All.PotablePOCPOU.bray["Potable_Overlap"]<-ANOSIM.R.All.PotablePOCPOU.bray.Potable_Overlap
ADONIS2.GroupedR.All.PotablePOCPOU.bray["Potable_Overlap"]<-ADONIS2.R2.All.PotablePOCPOU.bray.Potable_Overlap

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.PotablePOCPOU.meta$Potable_Overlap)

ADONIS2.Pairwise.All.PotablePOCPOU.bray.Potable_Overlap<-pairwise.adonis2(All.PotablePOCPOU.dist.bray ~ Potable_Overlap, 
                                                                          data = All.PotablePOCPOU.meta, 
                                                                          permutations = permutations)
ADONIS2.Pairwise.All.PotablePOCPOU.bray.Potable_Overlap
#Potable_Overlap
##Potable_Overlap




#################
#################
#################
#################
#Includes_Ozone 
ANOSIM.All.PotablePOCPOU.bray.Includes_Ozone<-anosim(All.PotablePOCPOU.dist.bray,
                                                      All.PotablePOCPOU.meta$Includes_Ozone, 
                                                      permutations = permutations)

ADONIS2.All.PotablePOCPOU.bray.Includes_Ozone<-adonis2(All.PotablePOCPOU.dist.bray ~ Includes_Ozone, 
                                                        data = All.PotablePOCPOU.meta, 
                                                        permutations = permutations)

Betadis.All.PotablePOCPOU.bray.Includes_Ozone<-anova(betadisper(All.PotablePOCPOU.dist.bray,
                                                                 All.PotablePOCPOU.meta$Includes_Ozone))



ANOSIM.Pval.All.PotablePOCPOU.bray.Includes_Ozone<-ANOSIM.All.PotablePOCPOU.bray.Includes_Ozone$signif
ANOSIM.Pval.All.PotablePOCPOU.bray.Includes_Ozone

ADONIS2.Pval.All.PotablePOCPOU.bray.Includes_Ozone<-ADONIS2.All.PotablePOCPOU.bray.Includes_Ozone$`Pr(>F)`[1]
ADONIS2.Pval.All.PotablePOCPOU.bray.Includes_Ozone  

Betadis.Pval.All.PotablePOCPOU.bray.Includes_Ozone<-Betadis.All.PotablePOCPOU.bray.Includes_Ozone$`Pr(>F)`[1]
Betadis.Pval.All.PotablePOCPOU.bray.Includes_Ozone

ANOSIM.R.All.PotablePOCPOU.bray.Includes_Ozone<-ANOSIM.All.PotablePOCPOU.bray.Includes_Ozone$statistic
ANOSIM.R.All.PotablePOCPOU.bray.Includes_Ozone

ADONIS2.R2.All.PotablePOCPOU.bray.Includes_Ozone<-ADONIS2.All.PotablePOCPOU.bray.Includes_Ozone$R2[1]
ADONIS2.R2.All.PotablePOCPOU.bray.Includes_Ozone  

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.PotablePOCPOU.bray["Includes_Ozone"]<-ANOSIM.Pval.All.PotablePOCPOU.bray.Includes_Ozone
ADONIS2.GroupedP.All.PotablePOCPOU.bray["Includes_Ozone"]<-ADONIS2.Pval.All.PotablePOCPOU.bray.Includes_Ozone

Betadis.GroupedP.All.PotablePOCPOU.bray["Includes_Ozone"]<-Betadis.Pval.All.PotablePOCPOU.bray.Includes_Ozone

ANOSIM.GroupedR.All.PotablePOCPOU.bray["Includes_Ozone"]<-ANOSIM.R.All.PotablePOCPOU.bray.Includes_Ozone
ADONIS2.GroupedR.All.PotablePOCPOU.bray["Includes_Ozone"]<-ADONIS2.R2.All.PotablePOCPOU.bray.Includes_Ozone

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.PotablePOCPOU.meta$Includes_Ozone)

ADONIS2.Pairwise.All.PotablePOCPOU.bray.Includes_Ozone<-pairwise.adonis2(All.PotablePOCPOU.dist.bray ~ Includes_Ozone, 
                                                                          data = All.PotablePOCPOU.meta, 
                                                                          permutations = permutations)
ADONIS2.Pairwise.All.PotablePOCPOU.bray.Includes_Ozone
#Includes_Ozone
##Includes_Ozone 




#################
#################
#################
#################
#Includes_RO_UF
ANOSIM.All.PotablePOCPOU.bray.Includes_RO_UF<-anosim(All.PotablePOCPOU.dist.bray,
                                                      All.PotablePOCPOU.meta$Includes_RO_UF,
                                                      permutations = permutations)

ADONIS2.All.PotablePOCPOU.bray.Includes_RO_UF<-adonis2(All.PotablePOCPOU.dist.bray ~ Includes_RO_UF,
                                                        data = All.PotablePOCPOU.meta,
                                                        permutations = permutations)

Betadis.All.PotablePOCPOU.bray.Includes_RO_UF<-anova(betadisper(All.PotablePOCPOU.dist.bray,
                                                                 All.PotablePOCPOU.meta$Includes_RO_UF))



ANOSIM.Pval.All.PotablePOCPOU.bray.Includes_RO_UF<-ANOSIM.All.PotablePOCPOU.bray.Includes_RO_UF$signif
ANOSIM.Pval.All.PotablePOCPOU.bray.Includes_RO_UF

ADONIS2.Pval.All.PotablePOCPOU.bray.Includes_RO_UF<-ADONIS2.All.PotablePOCPOU.bray.Includes_RO_UF$`Pr(>F)`[1]
ADONIS2.Pval.All.PotablePOCPOU.bray.Includes_RO_UF

Betadis.Pval.All.PotablePOCPOU.bray.Includes_RO_UF<-Betadis.All.PotablePOCPOU.bray.Includes_RO_UF$`Pr(>F)`[1]
Betadis.Pval.All.PotablePOCPOU.bray.Includes_RO_UF

ANOSIM.R.All.PotablePOCPOU.bray.Includes_RO_UF<-ANOSIM.All.PotablePOCPOU.bray.Includes_RO_UF$statistic
ANOSIM.R.All.PotablePOCPOU.bray.Includes_RO_UF

ADONIS2.R2.All.PotablePOCPOU.bray.Includes_RO_UF<-ADONIS2.All.PotablePOCPOU.bray.Includes_RO_UF$R2[1]
ADONIS2.R2.All.PotablePOCPOU.bray.Includes_RO_UF

#Pull Together Pvalues to adjust
ANOSIM.GroupedP.All.PotablePOCPOU.bray["Includes_RO_UF"]<-ANOSIM.Pval.All.PotablePOCPOU.bray.Includes_RO_UF
ADONIS2.GroupedP.All.PotablePOCPOU.bray["Includes_RO_UF"]<-ADONIS2.Pval.All.PotablePOCPOU.bray.Includes_RO_UF

Betadis.GroupedP.All.PotablePOCPOU.bray["Includes_RO_UF"]<-Betadis.Pval.All.PotablePOCPOU.bray.Includes_RO_UF

ANOSIM.GroupedR.All.PotablePOCPOU.bray["Includes_RO_UF"]<-ANOSIM.R.All.PotablePOCPOU.bray.Includes_RO_UF
ADONIS2.GroupedR.All.PotablePOCPOU.bray["Includes_RO_UF"]<-ADONIS2.R2.All.PotablePOCPOU.bray.Includes_RO_UF

#Pairwise if needed, check if interesting or not. Remove if not.
count(All.PotablePOCPOU.meta$Includes_RO_UF)

ADONIS2.Pairwise.All.PotablePOCPOU.bray.Includes_RO_UF<-pairwise.adonis2(All.PotablePOCPOU.dist.bray ~ Includes_RO_UF,
                                                                          data = All.PotablePOCPOU.meta,
                                                                          permutations = permutations)
ADONIS2.Pairwise.All.PotablePOCPOU.bray.Includes_RO_UF
#Includes_RO_UF
##Includes_RO_UF




#################
#################
#################
#################
#Includes_UV 
ANOSIM.All.PotablePOCPOU.bray.Includes_UV<-anosim(All.PotablePOCPOU.dist.bray,
                                                      All.PotablePOCPOU.meta$Includes_UV, 
                                                      permutations = permutations)

ADONIS2.All.PotablePOCPOU.bray.Includes_UV<-adonis2(All.PotablePOCPOU.dist.bray ~ Includes_UV, 
                                                        data = All.PotablePOCPOU.meta, 
                                                        permutations = permutations)

Betadis.All.PotablePOCPOU.bray.Includes_UV<-anova(betadisper(All.PotablePOCPOU.dist.bray,
                                                                 All.PotablePOCPOU.meta$Includes_UV))



ANOSIM.Pval.All.PotablePOCPOU.bray.Includes_UV<-ANOSIM.All.PotablePOCPOU.bray.Includes_UV$signif
ANOSIM.Pval.All.PotablePOCPOU.bray.Includes_UV

ADONIS2.Pval.All.PotablePOCPOU.bray.Includes_UV<-ADONIS2.All.PotablePOCPOU.bray.Includes_UV$`Pr(>F)`[1]
ADONIS2.Pval.All.PotablePOCPOU.bray.Includes_UV  

Betadis.Pval.All.PotablePOCPOU.bray.Includes_UV<-Betadis.All.PotablePOCPOU.bray.Includes_UV$`Pr(>F)`[1]
Betadis.Pval.All.PotablePOCPOU.bray.Includes_UV

ANOSIM.R.All.PotablePOCPOU.bray.Includes_UV<-ANOSIM.All.PotablePOCPOU.bray.Includes_UV$statistic
ANOSIM.R.All.PotablePOCPOU.bray.Includes_UV

ADONIS2.R2.All.PotablePOCPOU.bray.Includes_UV<-ADONIS2.All.PotablePOCPOU.bray.Includes_UV$R2[1]
ADONIS2.R2.All.PotablePOCPOU.bray.Includes_UV  

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.PotablePOCPOU.bray["Includes_UV"]<-ANOSIM.Pval.All.PotablePOCPOU.bray.Includes_UV
ADONIS2.GroupedP.All.PotablePOCPOU.bray["Includes_UV"]<-ADONIS2.Pval.All.PotablePOCPOU.bray.Includes_UV

Betadis.GroupedP.All.PotablePOCPOU.bray["Includes_UV"]<-Betadis.Pval.All.PotablePOCPOU.bray.Includes_UV

ANOSIM.GroupedR.All.PotablePOCPOU.bray["Includes_UV"]<-ANOSIM.R.All.PotablePOCPOU.bray.Includes_UV
ADONIS2.GroupedR.All.PotablePOCPOU.bray["Includes_UV"]<-ADONIS2.R2.All.PotablePOCPOU.bray.Includes_UV

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.PotablePOCPOU.meta$Includes_UV)

ADONIS2.Pairwise.All.PotablePOCPOU.bray.Includes_UV<-pairwise.adonis2(All.PotablePOCPOU.dist.bray ~ Includes_UV, 
                                                                          data = All.PotablePOCPOU.meta, 
                                                                          permutations = permutations)
ADONIS2.Pairwise.All.PotablePOCPOU.bray.Includes_UV
#Includes_UV
##Includes_UV




# #################
# #################
# #################
# #################
# #Type_DS  
# ANOSIM.All.PotablePOCPOU.bray.Type_DS <-anosim(All.PotablePOCPOU.dist.bray,
#                                                       All.PotablePOCPOU.meta$Type_DS , 
#                                                       permutations = permutations)
# 
# ADONIS2.All.PotablePOCPOU.bray.Type_DS <-adonis2(All.PotablePOCPOU.dist.bray ~ Type_DS , 
#                                                         data = All.PotablePOCPOU.meta, 
#                                                         permutations = permutations)
# 
# Betadis.All.PotablePOCPOU.bray.Type_DS <-anova(betadisper(All.PotablePOCPOU.dist.bray,
#                                                                  All.PotablePOCPOU.meta$Type_DS ))
# 
# 
# 
# ANOSIM.Pval.All.PotablePOCPOU.bray.Type_DS <-ANOSIM.All.PotablePOCPOU.bray.Type_DS $signif
# ANOSIM.Pval.All.PotablePOCPOU.bray.Type_DS 
# 
# ADONIS2.Pval.All.PotablePOCPOU.bray.Type_DS <-ADONIS2.All.PotablePOCPOU.bray.Type_DS $`Pr(>F)`[1]
# ADONIS2.Pval.All.PotablePOCPOU.bray.Type_DS   
# 
# Betadis.Pval.All.PotablePOCPOU.bray.Type_DS <-Betadis.All.PotablePOCPOU.bray.Type_DS $`Pr(>F)`[1]
# Betadis.Pval.All.PotablePOCPOU.bray.Type_DS 
# 
# ANOSIM.R.All.PotablePOCPOU.bray.Type_DS <-ANOSIM.All.PotablePOCPOU.bray.Type_DS $statistic
# ANOSIM.R.All.PotablePOCPOU.bray.Type_DS 
# 
# ADONIS2.R2.All.PotablePOCPOU.bray.Type_DS <-ADONIS2.All.PotablePOCPOU.bray.Type_DS $R2[1]
# ADONIS2.R2.All.PotablePOCPOU.bray.Type_DS   
# 
# #Pull Together Pvalues to adjust 
# ANOSIM.GroupedP.All.PotablePOCPOU.bray["Type_DS "]<-ANOSIM.Pval.All.PotablePOCPOU.bray.Type_DS 
# ADONIS2.GroupedP.All.PotablePOCPOU.bray["Type_DS "]<-ADONIS2.Pval.All.PotablePOCPOU.bray.Type_DS 
# 
# Betadis.GroupedP.All.PotablePOCPOU.bray["Type_DS "]<-Betadis.Pval.All.PotablePOCPOU.bray.Type_DS 
# 
# ANOSIM.GroupedR.All.PotablePOCPOU.bray["Type_DS "]<-ANOSIM.R.All.PotablePOCPOU.bray.Type_DS 
# ADONIS2.GroupedR.All.PotablePOCPOU.bray["Type_DS "]<-ADONIS2.R2.All.PotablePOCPOU.bray.Type_DS 
# 
# #Pairwise if needed, check if interesting or not. Remove if not. 
# count(All.PotablePOCPOU.meta$Type_DS )
# 
# ADONIS2.Pairwise.All.PotablePOCPOU.bray.Type_DS <-pairwise.adonis2(All.PotablePOCPOU.dist.bray ~ Type_DS , 
#                                                                           data = All.PotablePOCPOU.meta, 
#                                                                           permutations = permutations)
# ADONIS2.Pairwise.All.PotablePOCPOU.bray.Type_DS 
# #Type_DS 
# ##Type_DS 










##Potable_Source_Water
##Potable_Disinfection_Primary 
##Potable_Treatment_Type 
##Disinfection_Total 
##Treatment_Total 
##Treatment_Classification_Number
##Treatment_Classification_Cat
##Mixed_Source_Water
##Potable_Overlap
##Includes_Ozone 
##Includes_RO_UF 
##Includes_UV
##Type_DS 

#Adjust PValue, Check Significance
ANOSIM.GroupedP.Adjusted.All.PotablePOCPOU.bray<-p.adjust(ANOSIM.GroupedP.All.PotablePOCPOU.bray,method="BY")
p.adjust(ANOSIM.GroupedP.Adjusted.All.PotablePOCPOU.bray,method="BY")<.05

ADONIS2.GroupedP.Adjusted.All.PotablePOCPOU.bray<-p.adjust(ADONIS2.GroupedP.All.PotablePOCPOU.bray,method="BY")
p.adjust(ADONIS2.GroupedP.Adjusted.All.PotablePOCPOU.bray,method="BY")<.05

Betadis.GroupedP.All.PotablePOCPOU.bray<0.05
Betadis.GroupedP.Adjusted.All.PotablePOCPOU.bray<-p.adjust(Betadis.GroupedP.All.PotablePOCPOU.bray,method="BY")
p.adjust(Betadis.GroupedP.Adjusted.All.PotablePOCPOU.bray,method="BY")<.05


#Check Value of Pairwise Sub comparisons 
ADONIS2.Pairwise.All.PotablePOCPOU.bray.Class1_Pot_NonPot
ADONIS2.Pairwise.All.PotablePOCPOU.bray.Classification_1
ADONIS2.Pairwise.All.PotablePOCPOU.bray.Classification_1_mod
ADONIS2.Pairwise.All.PotablePOCPOU.bray.Classification_2
ADONIS2.Pairwise.All.PotablePOCPOU.bray.Climate
ADONIS2.Pairwise.All.PotablePOCPOU.bray.Climate_Region
ADONIS2.Pairwise.All.PotablePOCPOU.bray.Final_Disinfection_Residual
ADONIS2.Pairwise.All.PotablePOCPOU.bray.Primers
ADONIS2.Pairwise.All.PotablePOCPOU.bray.Region
ADONIS2.Pairwise.All.PotablePOCPOU.bray.Sample_Local
ADONIS2.Pairwise.All.PotablePOCPOU.bray.Trimmed
ADONIS2.Pairwise.All.PotablePOCPOU.bray.Extraction_Type
ADONIS2.Pairwise.All.PotablePOCPOU.bray.Potable_Source_Water
ADONIS2.Pairwise.All.PotablePOCPOU.bray.Potable_Disinfection_Primary
ADONIS2.Pairwise.All.PotablePOCPOU.bray.Potable_Treatment_Type
ADONIS2.Pairwise.All.PotablePOCPOU.bray.Disinfection_Total
ADONIS2.Pairwise.All.PotablePOCPOU.bray.Treatment_Total
ADONIS2.Pairwise.All.PotablePOCPOU.bray.Treatment_Classification_Number
ADONIS2.Pairwise.All.PotablePOCPOU.bray.Treatment_Classification_Cat
ADONIS2.Pairwise.All.PotablePOCPOU.bray.Mixed_Source_Water
ADONIS2.Pairwise.All.PotablePOCPOU.bray.Potable_Overlap
ADONIS2.Pairwise.All.PotablePOCPOU.bray.Includes_Ozone
ADONIS2.Pairwise.All.PotablePOCPOU.bray.Includes_RO_UF
ADONIS2.Pairwise.All.PotablePOCPOU.bray.Includes_UV
ADONIS2.Pairwise.All.PotablePOCPOU.bray.Type_DS

#Permanova on Multiple Variables that are significant and not heterogeneous (not sig from Beta, and sig from anosim)
ADONIS2.All.PotablePOCPOU.bray.multiple_comparisons<-adonis2(All.PotablePOCPOU.dist.bray ~ 
                                                        Classification_1+Climate, 
                                                 data = All.PotablePOCPOU.meta, 
                                                 permutations = permutations)
ADONIS2.All.PotablePOCPOU.bray.multiple_comparisons


#################
#################
##################################
#################
##################################
#################
##################################
#################
##################################
#################
##################################
#################
##################################
#################
#################
#EXPORT
#ANOSIM
All.PotablePOCPOU.Stats<-data.frame(
  "ANOSIM.P.All.PotablePOCPOU.bray"= ANOSIM.GroupedP.All.PotablePOCPOU.bray,
  "ANOSIM.P.Adj.All.PotablePOCPOU.bray"= ANOSIM.GroupedP.Adjusted.All.PotablePOCPOU.bray,
  "ADONIS2.P.All.PotablePOCPOU.bray"= ADONIS2.GroupedP.All.PotablePOCPOU.bray,
  "ADONIS2.P.Adj.All.PotablePOCPOU.bray"= ADONIS2.GroupedP.Adjusted.All.PotablePOCPOU.bray,
  "Betadis.P.All.PotablePOCPOU.bray"= Betadis.GroupedP.All.PotablePOCPOU.bray,
  "Betadis.P.Adj.All.PotablePOCPOU.bray"= Betadis.GroupedP.Adjusted.All.PotablePOCPOU.bray,
  "ANOSIM.R.All.PotablePOCPOU.bray"=  ANOSIM.GroupedR.All.PotablePOCPOU.bray,
  "ADONIS2.R.All.PotablePOCPOU.bray"= ADONIS2.GroupedR.All.PotablePOCPOU.bray)


write.xlsx(All.PotablePOCPOU.Stats, file = "C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Stats/Beta/All.PotablePOCPOU.List.xlsx", rowNames = TRUE)

