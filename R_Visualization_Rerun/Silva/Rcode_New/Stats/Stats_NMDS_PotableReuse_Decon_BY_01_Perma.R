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
Metadata_PotableReusePOCPOU<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_DS_Special_V2', rowNames = TRUE)
Metadata_POC<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_Final', rowNames = TRUE)
Metadata_POU<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'DS_Final', rowNames = TRUE)

#Order Meta Data 
Metadata_PotableReusePOCPOU$Classification_1 = factor(Metadata_PotableReusePOCPOU$Classification_1, levels=c("Potable Conventional","Potable Reuse","Non-potable Reuse","Blank"))
Metadata_POC$Classification_1 = factor(Metadata_POC$Classification_1, levels=c("Potable Conventional","Potable Reuse","Non-potable Reuse","Blank"))
Metadata_POU$Classification_1 = factor(Metadata_POU$Classification_1, levels=c("Potable Conventional","Potable Reuse","Non-potable Reuse","Blank"))

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


#Anosim/Adonis/Betadispersity 
#Can adjust permutations for higher accuracy. 
permutations=1000


All.PotableReusePOCPOU.dist.bray<- phyloseq::distance(physeq.combined.All.PotableReusePOCPOU, method = "bray")
All.PotableReusePOCPOU.meta<- data.frame(sample_data(physeq.combined.All.PotableReusePOCPOU))

ANOSIM.GroupedP.All.PotableReusePOCPOU.bray<-numeric() 
ADONIS2.GroupedP.All.PotableReusePOCPOU.bray<-numeric()
Betadis.GroupedP.All.PotableReusePOCPOU.bray<-numeric()
ANOSIM.GroupedR.All.PotableReusePOCPOU.bray<-numeric() 
ADONIS2.GroupedR.All.PotableReusePOCPOU.bray<-numeric()


# #################
# #################
# #################
# #################
# #Class1_Pot_NonPot
# ANOSIM.All.PotableReusePOCPOU.bray.Class1_Pot_NonPot<-anosim(All.PotableReusePOCPOU.dist.bray,
#                                All.PotableReusePOCPOU.meta$Class1_Pot_NonPot,
#                                permutations = permutations)
#
# ADONIS2.All.PotableReusePOCPOU.bray.Class1_Pot_NonPot<-adonis2(All.PotableReusePOCPOU.dist.bray ~ Class1_Pot_NonPot,
#                                  data = All.PotableReusePOCPOU.meta,
#                                  permutations = permutations)
#
# Betadis.All.PotableReusePOCPOU.bray.Class1_Pot_NonPot<-anova(betadisper(All.PotableReusePOCPOU.dist.bray,
#                                                       All.PotableReusePOCPOU.meta$Class1_Pot_NonPot))
#
#
#
# ANOSIM.Pval.All.PotableReusePOCPOU.bray.Class1_Pot_NonPot<-ANOSIM.All.PotableReusePOCPOU.bray.Class1_Pot_NonPot$signif
# ANOSIM.Pval.All.PotableReusePOCPOU.bray.Class1_Pot_NonPot
#
# ADONIS2.Pval.All.PotableReusePOCPOU.bray.Class1_Pot_NonPot<-ADONIS2.All.PotableReusePOCPOU.bray.Class1_Pot_NonPot$`Pr(>F)`[1]
# ADONIS2.Pval.All.PotableReusePOCPOU.bray.Class1_Pot_NonPot
#
# Betadis.Pval.All.PotableReusePOCPOU.bray.Class1_Pot_NonPot<-Betadis.All.PotableReusePOCPOU.bray.Class1_Pot_NonPot$`Pr(>F)`[1]
# Betadis.Pval.All.PotableReusePOCPOU.bray.Class1_Pot_NonPot
#
# ANOSIM.R.All.PotableReusePOCPOU.bray.Class1_Pot_NonPot<-ANOSIM.All.PotableReusePOCPOU.bray.Class1_Pot_NonPot$statistic
# ANOSIM.R.All.PotableReusePOCPOU.bray.Class1_Pot_NonPot
#
# ADONIS2.R2.All.PotableReusePOCPOU.bray.Class1_Pot_NonPot<-ADONIS2.All.PotableReusePOCPOU.bray.Class1_Pot_NonPot$R2[1]
# ADONIS2.R2.All.PotableReusePOCPOU.bray.Class1_Pot_NonPot
#
# #Pull Together Pvalues to adjust
# ANOSIM.GroupedP.All.PotableReusePOCPOU.bray["Class1_Pot_NonPot"]<-ANOSIM.Pval.All.PotableReusePOCPOU.bray.Class1_Pot_NonPot
# ADONIS2.GroupedP.All.PotableReusePOCPOU.bray["Class1_Pot_NonPot"]<-ADONIS2.Pval.All.PotableReusePOCPOU.bray.Class1_Pot_NonPot
#
# Betadis.GroupedP.All.PotableReusePOCPOU.bray["Class1_Pot_NonPot"]<-Betadis.Pval.All.PotableReusePOCPOU.bray.Class1_Pot_NonPot
#
# ANOSIM.GroupedR.All.PotableReusePOCPOU.bray["Class1_Pot_NonPot"]<-ANOSIM.R.All.PotableReusePOCPOU.bray.Class1_Pot_NonPot
# ADONIS2.GroupedR.All.PotableReusePOCPOU.bray["Class1_Pot_NonPot"]<-ADONIS2.R2.All.PotableReusePOCPOU.bray.Class1_Pot_NonPot
#
# #Pairwise if needed, check if interesting or not. Remove if not.
# count(All.PotableReusePOCPOU.meta$Class1_Pot_NonPot)
#
# ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Class1_Pot_NonPot<-pairwise.adonis2(All.PotableReusePOCPOU.dist.bray ~ Class1_Pot_NonPot,
#                                                                     data = All.PotableReusePOCPOU.meta,
#                                                                     permutations = permutations)
# ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Class1_Pot_NonPot
#
#
#
#
# #################
# #################
# #################
# #################
# #Classification_1
# ANOSIM.All.PotableReusePOCPOU.bray.Classification_1<-anosim(All.PotableReusePOCPOU.dist.bray,
#                                                  All.PotableReusePOCPOU.meta$Classification_1,
#                                                  permutations = permutations)
#
# ADONIS2.All.PotableReusePOCPOU.bray.Classification_1<-adonis2(All.PotableReusePOCPOU.dist.bray ~ Classification_1,
#                                                    data = All.PotableReusePOCPOU.meta,
#                                                    permutations = permutations)
#
# Betadis.All.PotableReusePOCPOU.bray.Classification_1<-anova(betadisper(All.PotableReusePOCPOU.dist.bray,
#                                                             All.PotableReusePOCPOU.meta$Classification_1))
#
#
#
# ANOSIM.Pval.All.PotableReusePOCPOU.bray.Classification_1<-ANOSIM.All.PotableReusePOCPOU.bray.Classification_1$signif
# ANOSIM.Pval.All.PotableReusePOCPOU.bray.Classification_1
#
# ADONIS2.Pval.All.PotableReusePOCPOU.bray.Classification_1<-ADONIS2.All.PotableReusePOCPOU.bray.Classification_1$`Pr(>F)`[1]
# ADONIS2.Pval.All.PotableReusePOCPOU.bray.Classification_1
#
# Betadis.Pval.All.PotableReusePOCPOU.bray.Classification_1<-Betadis.All.PotableReusePOCPOU.bray.Classification_1$`Pr(>F)`[1]
# Betadis.Pval.All.PotableReusePOCPOU.bray.Classification_1
#
# ANOSIM.R.All.PotableReusePOCPOU.bray.Classification_1<-ANOSIM.All.PotableReusePOCPOU.bray.Classification_1$statistic
# ANOSIM.R.All.PotableReusePOCPOU.bray.Classification_1
#
# ADONIS2.R2.All.PotableReusePOCPOU.bray.Classification_1<-ADONIS2.All.PotableReusePOCPOU.bray.Classification_1$R2[1]
# ADONIS2.R2.All.PotableReusePOCPOU.bray.Classification_1
#
# #Pull Together Pvalues to adjust
# ANOSIM.GroupedP.All.PotableReusePOCPOU.bray["Classification_1"]<-ANOSIM.Pval.All.PotableReusePOCPOU.bray.Classification_1
# ADONIS2.GroupedP.All.PotableReusePOCPOU.bray["Classification_1"]<-ADONIS2.Pval.All.PotableReusePOCPOU.bray.Classification_1
#
# Betadis.GroupedP.All.PotableReusePOCPOU.bray["Classification_1"]<-Betadis.Pval.All.PotableReusePOCPOU.bray.Classification_1
#
# ANOSIM.GroupedR.All.PotableReusePOCPOU.bray["Classification_1"]<-ANOSIM.R.All.PotableReusePOCPOU.bray.Classification_1
# ADONIS2.GroupedR.All.PotableReusePOCPOU.bray["Classification_1"]<-ADONIS2.R2.All.PotableReusePOCPOU.bray.Classification_1
#
# #Pairwise if needed, check if interesting or not. Remove if not.
# count(All.PotableReusePOCPOU.meta$Classification_1)
#
# ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Classification_1<-pairwise.adonis2(All.PotableReusePOCPOU.dist.bray ~ Classification_1,
#                                                                      data = All.PotableReusePOCPOU.meta,
#                                                                      permutations = permutations)
# ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Classification_1
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
# ANOSIM.All.PotableReusePOCPOU.bray.Classification_1_mod<-anosim(All.PotableReusePOCPOU.dist.bray,
#                                                  All.PotableReusePOCPOU.meta$Classification_1_mod,
#                                                  permutations = permutations)
#
# ADONIS2.All.PotableReusePOCPOU.bray.Classification_1_mod<-adonis2(All.PotableReusePOCPOU.dist.bray ~ Classification_1_mod,
#                                                    data = All.PotableReusePOCPOU.meta,
#                                                    permutations = permutations)
#
# Betadis.All.PotableReusePOCPOU.bray.Classification_1_mod<-anova(betadisper(All.PotableReusePOCPOU.dist.bray,
#                                                             All.PotableReusePOCPOU.meta$Classification_1_mod))
#
#
#
# ANOSIM.Pval.All.PotableReusePOCPOU.bray.Classification_1_mod<-ANOSIM.All.PotableReusePOCPOU.bray.Classification_1_mod$signif
# ANOSIM.Pval.All.PotableReusePOCPOU.bray.Classification_1_mod
#
# ADONIS2.Pval.All.PotableReusePOCPOU.bray.Classification_1_mod<-ADONIS2.All.PotableReusePOCPOU.bray.Classification_1_mod$`Pr(>F)`[1]
# ADONIS2.Pval.All.PotableReusePOCPOU.bray.Classification_1_mod
#
# Betadis.Pval.All.PotableReusePOCPOU.bray.Classification_1_mod<-Betadis.All.PotableReusePOCPOU.bray.Classification_1_mod$`Pr(>F)`[1]
# Betadis.Pval.All.PotableReusePOCPOU.bray.Classification_1_mod
#
# ANOSIM.R.All.PotableReusePOCPOU.bray.Classification_1_mod<-ANOSIM.All.PotableReusePOCPOU.bray.Classification_1_mod$statistic
# ANOSIM.R.All.PotableReusePOCPOU.bray.Classification_1_mod
#
# ADONIS2.R2.All.PotableReusePOCPOU.bray.Classification_1_mod<-ADONIS2.All.PotableReusePOCPOU.bray.Classification_1_mod$R2[1]
# ADONIS2.R2.All.PotableReusePOCPOU.bray.Classification_1_mod
#
# #Pull Together Pvalues to adjust
# ANOSIM.GroupedP.All.PotableReusePOCPOU.bray["Classification_1_mod"]<-ANOSIM.Pval.All.PotableReusePOCPOU.bray.Classification_1_mod
# ADONIS2.GroupedP.All.PotableReusePOCPOU.bray["Classification_1_mod"]<-ADONIS2.Pval.All.PotableReusePOCPOU.bray.Classification_1_mod
#
# Betadis.GroupedP.All.PotableReusePOCPOU.bray["Classification_1_mod"]<-Betadis.Pval.All.PotableReusePOCPOU.bray.Classification_1_mod
#
# ANOSIM.GroupedR.All.PotableReusePOCPOU.bray["Classification_1_mod"]<-ANOSIM.R.All.PotableReusePOCPOU.bray.Classification_1_mod
# ADONIS2.GroupedR.All.PotableReusePOCPOU.bray["Classification_1_mod"]<-ADONIS2.R2.All.PotableReusePOCPOU.bray.Classification_1_mod
#
# #Pairwise if needed, check if interesting or not. Remove if not.
# count(All.PotableReusePOCPOU.meta$Classification_1_mod)
#
# ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Classification_1_mod<-pairwise.adonis2(All.PotableReusePOCPOU.dist.bray ~ Classification_1_mod,
#                                                                      data = All.PotableReusePOCPOU.meta,
#                                                                      permutations = permutations)
# ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Classification_1_mod





#################
#################
#################
#################
#Classification_2v2
ANOSIM.All.PotableReusePOCPOU.bray.Classification_2<-anosim(All.PotableReusePOCPOU.dist.bray,
                                                 All.PotableReusePOCPOU.meta$Classification_2v2,
                                                 permutations = permutations)

ADONIS2.All.PotableReusePOCPOU.bray.Classification_2<-adonis2(All.PotableReusePOCPOU.dist.bray ~ Classification_2,
                                                   data = All.PotableReusePOCPOU.meta,
                                                   permutations = permutations)

Betadis.All.PotableReusePOCPOU.bray.Classification_2<-anova(betadisper(All.PotableReusePOCPOU.dist.bray,
                                                            All.PotableReusePOCPOU.meta$Classification_2v2))



ANOSIM.Pval.All.PotableReusePOCPOU.bray.Classification_2<-ANOSIM.All.PotableReusePOCPOU.bray.Classification_2$signif
ANOSIM.Pval.All.PotableReusePOCPOU.bray.Classification_2

ADONIS2.Pval.All.PotableReusePOCPOU.bray.Classification_2<-ADONIS2.All.PotableReusePOCPOU.bray.Classification_2$`Pr(>F)`[1]
ADONIS2.Pval.All.PotableReusePOCPOU.bray.Classification_2

Betadis.Pval.All.PotableReusePOCPOU.bray.Classification_2<-Betadis.All.PotableReusePOCPOU.bray.Classification_2$`Pr(>F)`[1]
Betadis.Pval.All.PotableReusePOCPOU.bray.Classification_2

ANOSIM.R.All.PotableReusePOCPOU.bray.Classification_2<-ANOSIM.All.PotableReusePOCPOU.bray.Classification_2$statistic
ANOSIM.R.All.PotableReusePOCPOU.bray.Classification_2

ADONIS2.R2.All.PotableReusePOCPOU.bray.Classification_2<-ADONIS2.All.PotableReusePOCPOU.bray.Classification_2$R2[1]
ADONIS2.R2.All.PotableReusePOCPOU.bray.Classification_2

#Pull Together Pvalues to adjust
ANOSIM.GroupedP.All.PotableReusePOCPOU.bray["Classification_2"]<-ANOSIM.Pval.All.PotableReusePOCPOU.bray.Classification_2
ADONIS2.GroupedP.All.PotableReusePOCPOU.bray["Classification_2"]<-ADONIS2.Pval.All.PotableReusePOCPOU.bray.Classification_2

Betadis.GroupedP.All.PotableReusePOCPOU.bray["Classification_2"]<-Betadis.Pval.All.PotableReusePOCPOU.bray.Classification_2

ANOSIM.GroupedR.All.PotableReusePOCPOU.bray["Classification_2"]<-ANOSIM.R.All.PotableReusePOCPOU.bray.Classification_2
ADONIS2.GroupedR.All.PotableReusePOCPOU.bray["Classification_2"]<-ADONIS2.R2.All.PotableReusePOCPOU.bray.Classification_2

#Pairwise if needed, check if interesting or not. Remove if not.
count(All.PotableReusePOCPOU.meta$Classification_2v2)

ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Classification_2<-pairwise.adonis2(All.PotableReusePOCPOU.dist.bray ~ Classification_2,
                                                                     data = All.PotableReusePOCPOU.meta,
                                                                     permutations = permutations)
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Classification_2
#Classification_2





# #################
# #################
# #################
# #################
# #Climate_Region
# ANOSIM.All.PotableReusePOCPOU.bray.Climate_Region<-anosim(All.PotableReusePOCPOU.dist.bray,
#                                                  All.PotableReusePOCPOU.meta$Climate_Region,
#                                                  permutations = permutations)
#
# ADONIS2.All.PotableReusePOCPOU.bray.Climate_Region<-adonis2(All.PotableReusePOCPOU.dist.bray ~ Climate_Region,
#                                                    data = All.PotableReusePOCPOU.meta,
#                                                    permutations = permutations)
#
# Betadis.All.PotableReusePOCPOU.bray.Climate_Region<-anova(betadisper(All.PotableReusePOCPOU.dist.bray,
#                                                             All.PotableReusePOCPOU.meta$Climate_Region))
#
#
#
# ANOSIM.Pval.All.PotableReusePOCPOU.bray.Climate_Region<-ANOSIM.All.PotableReusePOCPOU.bray.Climate_Region$signif
# ANOSIM.Pval.All.PotableReusePOCPOU.bray.Climate_Region
#
# ADONIS2.Pval.All.PotableReusePOCPOU.bray.Climate_Region<-ADONIS2.All.PotableReusePOCPOU.bray.Climate_Region$`Pr(>F)`[1]
# ADONIS2.Pval.All.PotableReusePOCPOU.bray.Climate_Region
#
# Betadis.Pval.All.PotableReusePOCPOU.bray.Climate_Region<-Betadis.All.PotableReusePOCPOU.bray.Climate_Region$`Pr(>F)`[1]
# Betadis.Pval.All.PotableReusePOCPOU.bray.Climate_Region
#
# ANOSIM.R.All.PotableReusePOCPOU.bray.Climate_Region<-ANOSIM.All.PotableReusePOCPOU.bray.Climate_Region$statistic
# ANOSIM.R.All.PotableReusePOCPOU.bray.Climate_Region
#
# ADONIS2.R2.All.PotableReusePOCPOU.bray.Climate_Region<-ADONIS2.All.PotableReusePOCPOU.bray.Climate_Region$R2[1]
# ADONIS2.R2.All.PotableReusePOCPOU.bray.Climate_Region
#
# #Pull Together Pvalues to adjust
# ANOSIM.GroupedP.All.PotableReusePOCPOU.bray["Climate_Region"]<-ANOSIM.Pval.All.PotableReusePOCPOU.bray.Climate_Region
# ADONIS2.GroupedP.All.PotableReusePOCPOU.bray["Climate_Region"]<-ADONIS2.Pval.All.PotableReusePOCPOU.bray.Climate_Region
#
# Betadis.GroupedP.All.PotableReusePOCPOU.bray["Climate_Region"]<-Betadis.Pval.All.PotableReusePOCPOU.bray.Climate_Region
#
# ANOSIM.GroupedR.All.PotableReusePOCPOU.bray["Climate_Region"]<-ANOSIM.R.All.PotableReusePOCPOU.bray.Climate_Region
# ADONIS2.GroupedR.All.PotableReusePOCPOU.bray["Climate_Region"]<-ADONIS2.R2.All.PotableReusePOCPOU.bray.Climate_Region
#
# #Pairwise if needed, check if interesting or not. Remove if not.
# count(All.PotableReusePOCPOU.meta$Climate_Region)
#
# ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Climate_Region<-pairwise.adonis2(All.PotableReusePOCPOU.dist.bray ~ Climate_Region,
#                                                                      data = All.PotableReusePOCPOU.meta,
#                                                                      permutations = permutations)
# ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Climate_Region
# #Climate_Region





#################
#################
#################
#################
#Climate
ANOSIM.All.PotableReusePOCPOU.bray.Climate<-anosim(All.PotableReusePOCPOU.dist.bray,
                                                 All.PotableReusePOCPOU.meta$Climate,
                                                 permutations = permutations)

ADONIS2.All.PotableReusePOCPOU.bray.Climate<-adonis2(All.PotableReusePOCPOU.dist.bray ~ Climate,
                                                   data = All.PotableReusePOCPOU.meta,
                                                   permutations = permutations)

Betadis.All.PotableReusePOCPOU.bray.Climate<-anova(betadisper(All.PotableReusePOCPOU.dist.bray,
                                                            All.PotableReusePOCPOU.meta$Climate))



ANOSIM.Pval.All.PotableReusePOCPOU.bray.Climate<-ANOSIM.All.PotableReusePOCPOU.bray.Climate$signif
ANOSIM.Pval.All.PotableReusePOCPOU.bray.Climate

ADONIS2.Pval.All.PotableReusePOCPOU.bray.Climate<-ADONIS2.All.PotableReusePOCPOU.bray.Climate$`Pr(>F)`[1]
ADONIS2.Pval.All.PotableReusePOCPOU.bray.Climate

Betadis.Pval.All.PotableReusePOCPOU.bray.Climate<-Betadis.All.PotableReusePOCPOU.bray.Climate$`Pr(>F)`[1]
Betadis.Pval.All.PotableReusePOCPOU.bray.Climate

ANOSIM.R.All.PotableReusePOCPOU.bray.Climate<-ANOSIM.All.PotableReusePOCPOU.bray.Climate$statistic
ANOSIM.R.All.PotableReusePOCPOU.bray.Climate

ADONIS2.R2.All.PotableReusePOCPOU.bray.Climate<-ADONIS2.All.PotableReusePOCPOU.bray.Climate$R2[1]
ADONIS2.R2.All.PotableReusePOCPOU.bray.Climate

#Pull Together Pvalues to adjust
ANOSIM.GroupedP.All.PotableReusePOCPOU.bray["Climate"]<-ANOSIM.Pval.All.PotableReusePOCPOU.bray.Climate
ADONIS2.GroupedP.All.PotableReusePOCPOU.bray["Climate"]<-ADONIS2.Pval.All.PotableReusePOCPOU.bray.Climate

Betadis.GroupedP.All.PotableReusePOCPOU.bray["Climate"]<-Betadis.Pval.All.PotableReusePOCPOU.bray.Climate

ANOSIM.GroupedR.All.PotableReusePOCPOU.bray["Climate"]<-ANOSIM.R.All.PotableReusePOCPOU.bray.Climate
ADONIS2.GroupedR.All.PotableReusePOCPOU.bray["Climate"]<-ADONIS2.R2.All.PotableReusePOCPOU.bray.Climate

#Pairwise if needed, check if interesting or not. Remove if not.
count(All.PotableReusePOCPOU.meta$Climate)

ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Climate<-pairwise.adonis2(All.PotableReusePOCPOU.dist.bray ~ Climate,
                                                                     data = All.PotableReusePOCPOU.meta,
                                                                     permutations = permutations)
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Climate
#Climate





#################
#################
#################
#################
#Region
ANOSIM.All.PotableReusePOCPOU.bray.Region<-anosim(All.PotableReusePOCPOU.dist.bray,
                                                 All.PotableReusePOCPOU.meta$Region,
                                                 permutations = permutations)

ADONIS2.All.PotableReusePOCPOU.bray.Region<-adonis2(All.PotableReusePOCPOU.dist.bray ~ Region,
                                                   data = All.PotableReusePOCPOU.meta,
                                                   permutations = permutations)

Betadis.All.PotableReusePOCPOU.bray.Region<-anova(betadisper(All.PotableReusePOCPOU.dist.bray,
                                                            All.PotableReusePOCPOU.meta$Region))



ANOSIM.Pval.All.PotableReusePOCPOU.bray.Region<-ANOSIM.All.PotableReusePOCPOU.bray.Region$signif
ANOSIM.Pval.All.PotableReusePOCPOU.bray.Region

ADONIS2.Pval.All.PotableReusePOCPOU.bray.Region<-ADONIS2.All.PotableReusePOCPOU.bray.Region$`Pr(>F)`[1]
ADONIS2.Pval.All.PotableReusePOCPOU.bray.Region

Betadis.Pval.All.PotableReusePOCPOU.bray.Region<-Betadis.All.PotableReusePOCPOU.bray.Region$`Pr(>F)`[1]
Betadis.Pval.All.PotableReusePOCPOU.bray.Region

ANOSIM.R.All.PotableReusePOCPOU.bray.Region<-ANOSIM.All.PotableReusePOCPOU.bray.Region$statistic
ANOSIM.R.All.PotableReusePOCPOU.bray.Region

ADONIS2.R2.All.PotableReusePOCPOU.bray.Region<-ADONIS2.All.PotableReusePOCPOU.bray.Region$R2[1]
ADONIS2.R2.All.PotableReusePOCPOU.bray.Region

#Pull Together Pvalues to adjust
ANOSIM.GroupedP.All.PotableReusePOCPOU.bray["Region"]<-ANOSIM.Pval.All.PotableReusePOCPOU.bray.Region
ADONIS2.GroupedP.All.PotableReusePOCPOU.bray["Region"]<-ADONIS2.Pval.All.PotableReusePOCPOU.bray.Region

Betadis.GroupedP.All.PotableReusePOCPOU.bray["Region"]<-Betadis.Pval.All.PotableReusePOCPOU.bray.Region

ANOSIM.GroupedR.All.PotableReusePOCPOU.bray["Region"]<-ANOSIM.R.All.PotableReusePOCPOU.bray.Region
ADONIS2.GroupedR.All.PotableReusePOCPOU.bray["Region"]<-ADONIS2.R2.All.PotableReusePOCPOU.bray.Region

#Pairwise if needed, check if interesting or not. Remove if not.
count(All.PotableReusePOCPOU.meta$Region)

ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Region<-pairwise.adonis2(All.PotableReusePOCPOU.dist.bray ~ Region,
                                                                     data = All.PotableReusePOCPOU.meta,
                                                                     permutations = permutations)
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Region
#Region





#################
#################
#################
#################
#Sample_Local
ANOSIM.All.PotableReusePOCPOU.bray.Sample_Local<-anosim(All.PotableReusePOCPOU.dist.bray,
                                                 All.PotableReusePOCPOU.meta$Sample_Local,
                                                 permutations = permutations)

ADONIS2.All.PotableReusePOCPOU.bray.Sample_Local<-adonis2(All.PotableReusePOCPOU.dist.bray ~ Sample_Local,
                                                   data = All.PotableReusePOCPOU.meta,
                                                   permutations = permutations)

Betadis.All.PotableReusePOCPOU.bray.Sample_Local<-anova(betadisper(All.PotableReusePOCPOU.dist.bray,
                                                            All.PotableReusePOCPOU.meta$Sample_Local))



ANOSIM.Pval.All.PotableReusePOCPOU.bray.Sample_Local<-ANOSIM.All.PotableReusePOCPOU.bray.Sample_Local$signif
ANOSIM.Pval.All.PotableReusePOCPOU.bray.Sample_Local

ADONIS2.Pval.All.PotableReusePOCPOU.bray.Sample_Local<-ADONIS2.All.PotableReusePOCPOU.bray.Sample_Local$`Pr(>F)`[1]
ADONIS2.Pval.All.PotableReusePOCPOU.bray.Sample_Local

Betadis.Pval.All.PotableReusePOCPOU.bray.Sample_Local<-Betadis.All.PotableReusePOCPOU.bray.Sample_Local$`Pr(>F)`[1]
Betadis.Pval.All.PotableReusePOCPOU.bray.Sample_Local

ANOSIM.R.All.PotableReusePOCPOU.bray.Sample_Local<-ANOSIM.All.PotableReusePOCPOU.bray.Sample_Local$statistic
ANOSIM.R.All.PotableReusePOCPOU.bray.Sample_Local

ADONIS2.R2.All.PotableReusePOCPOU.bray.Sample_Local<-ADONIS2.All.PotableReusePOCPOU.bray.Sample_Local$R2[1]
ADONIS2.R2.All.PotableReusePOCPOU.bray.Sample_Local

#Pull Together Pvalues to adjust
ANOSIM.GroupedP.All.PotableReusePOCPOU.bray["Sample_Local"]<-ANOSIM.Pval.All.PotableReusePOCPOU.bray.Sample_Local
ADONIS2.GroupedP.All.PotableReusePOCPOU.bray["Sample_Local"]<-ADONIS2.Pval.All.PotableReusePOCPOU.bray.Sample_Local

Betadis.GroupedP.All.PotableReusePOCPOU.bray["Sample_Local"]<-Betadis.Pval.All.PotableReusePOCPOU.bray.Sample_Local

ANOSIM.GroupedR.All.PotableReusePOCPOU.bray["Sample_Local"]<-ANOSIM.R.All.PotableReusePOCPOU.bray.Sample_Local
ADONIS2.GroupedR.All.PotableReusePOCPOU.bray["Sample_Local"]<-ADONIS2.R2.All.PotableReusePOCPOU.bray.Sample_Local

#Pairwise if needed, check if interesting or not. Remove if not.
count(All.PotableReusePOCPOU.meta$Sample_Local)

ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Sample_Local<-pairwise.adonis2(All.PotableReusePOCPOU.dist.bray ~ Sample_Local,
                                                                     data = All.PotableReusePOCPOU.meta,
                                                                     permutations = permutations)
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Sample_Local
#Sample_Local





#################
#################
#################
#################
#Final_Disinfection_Residual
ANOSIM.All.PotableReusePOCPOU.bray.Final_Disinfection_Residual<-anosim(All.PotableReusePOCPOU.dist.bray,
                                                 All.PotableReusePOCPOU.meta$Final_Disinfection_Residual,
                                                 permutations = permutations)

ADONIS2.All.PotableReusePOCPOU.bray.Final_Disinfection_Residual<-adonis2(All.PotableReusePOCPOU.dist.bray ~ Final_Disinfection_Residual,
                                                   data = All.PotableReusePOCPOU.meta,
                                                   permutations = permutations)

Betadis.All.PotableReusePOCPOU.bray.Final_Disinfection_Residual<-anova(betadisper(All.PotableReusePOCPOU.dist.bray,
                                                            All.PotableReusePOCPOU.meta$Final_Disinfection_Residual))



ANOSIM.Pval.All.PotableReusePOCPOU.bray.Final_Disinfection_Residual<-ANOSIM.All.PotableReusePOCPOU.bray.Final_Disinfection_Residual$signif
ANOSIM.Pval.All.PotableReusePOCPOU.bray.Final_Disinfection_Residual

ADONIS2.Pval.All.PotableReusePOCPOU.bray.Final_Disinfection_Residual<-ADONIS2.All.PotableReusePOCPOU.bray.Final_Disinfection_Residual$`Pr(>F)`[1]
ADONIS2.Pval.All.PotableReusePOCPOU.bray.Final_Disinfection_Residual

Betadis.Pval.All.PotableReusePOCPOU.bray.Final_Disinfection_Residual<-Betadis.All.PotableReusePOCPOU.bray.Final_Disinfection_Residual$`Pr(>F)`[1]
Betadis.Pval.All.PotableReusePOCPOU.bray.Final_Disinfection_Residual

ANOSIM.R.All.PotableReusePOCPOU.bray.Final_Disinfection_Residual<-ANOSIM.All.PotableReusePOCPOU.bray.Final_Disinfection_Residual$statistic
ANOSIM.R.All.PotableReusePOCPOU.bray.Final_Disinfection_Residual

ADONIS2.R2.All.PotableReusePOCPOU.bray.Final_Disinfection_Residual<-ADONIS2.All.PotableReusePOCPOU.bray.Final_Disinfection_Residual$R2[1]
ADONIS2.R2.All.PotableReusePOCPOU.bray.Final_Disinfection_Residual

#Pull Together Pvalues to adjust
ANOSIM.GroupedP.All.PotableReusePOCPOU.bray["Final_Disinfection_Residual"]<-ANOSIM.Pval.All.PotableReusePOCPOU.bray.Final_Disinfection_Residual
ADONIS2.GroupedP.All.PotableReusePOCPOU.bray["Final_Disinfection_Residual"]<-ADONIS2.Pval.All.PotableReusePOCPOU.bray.Final_Disinfection_Residual

Betadis.GroupedP.All.PotableReusePOCPOU.bray["Final_Disinfection_Residual"]<-Betadis.Pval.All.PotableReusePOCPOU.bray.Final_Disinfection_Residual

ANOSIM.GroupedR.All.PotableReusePOCPOU.bray["Final_Disinfection_Residual"]<-ANOSIM.R.All.PotableReusePOCPOU.bray.Final_Disinfection_Residual
ADONIS2.GroupedR.All.PotableReusePOCPOU.bray["Final_Disinfection_Residual"]<-ADONIS2.R2.All.PotableReusePOCPOU.bray.Final_Disinfection_Residual

#Pairwise if needed, check if interesting or not. Remove if not.
count(All.PotableReusePOCPOU.meta$Final_Disinfection_Residual)

ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Final_Disinfection_Residual<-pairwise.adonis2(All.PotableReusePOCPOU.dist.bray ~ Final_Disinfection_Residual,
                                                                     data = All.PotableReusePOCPOU.meta,
                                                                     permutations = permutations)
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Final_Disinfection_Residual
#Final_Disinfection_Residual





# #################
# #################
# #################
# #################
# #Primers
# ANOSIM.All.PotableReusePOCPOU.bray.Primers<-anosim(All.PotableReusePOCPOU.dist.bray,
#                                                  All.PotableReusePOCPOU.meta$Primers,
#                                                  permutations = permutations)
#
# ADONIS2.All.PotableReusePOCPOU.bray.Primers<-adonis2(All.PotableReusePOCPOU.dist.bray ~ Primers,
#                                                    data = All.PotableReusePOCPOU.meta,
#                                                    permutations = permutations)
#
# Betadis.All.PotableReusePOCPOU.bray.Primers<-anova(betadisper(All.PotableReusePOCPOU.dist.bray,
#                                                             All.PotableReusePOCPOU.meta$Primers))
#
#
#
# ANOSIM.Pval.All.PotableReusePOCPOU.bray.Primers<-ANOSIM.All.PotableReusePOCPOU.bray.Primers$signif
# ANOSIM.Pval.All.PotableReusePOCPOU.bray.Primers
#
# ADONIS2.Pval.All.PotableReusePOCPOU.bray.Primers<-ADONIS2.All.PotableReusePOCPOU.bray.Primers$`Pr(>F)`[1]
# ADONIS2.Pval.All.PotableReusePOCPOU.bray.Primers
#
# Betadis.Pval.All.PotableReusePOCPOU.bray.Primers<-Betadis.All.PotableReusePOCPOU.bray.Primers$`Pr(>F)`[1]
# Betadis.Pval.All.PotableReusePOCPOU.bray.Primers
#
# ANOSIM.R.All.PotableReusePOCPOU.bray.Primers<-ANOSIM.All.PotableReusePOCPOU.bray.Primers$statistic
# ANOSIM.R.All.PotableReusePOCPOU.bray.Primers
#
# ADONIS2.R2.All.PotableReusePOCPOU.bray.Primers<-ADONIS2.All.PotableReusePOCPOU.bray.Primers$R2[1]
# ADONIS2.R2.All.PotableReusePOCPOU.bray.Primers
#
# #Pull Together Pvalues to adjust
# ANOSIM.GroupedP.All.PotableReusePOCPOU.bray["Primers"]<-ANOSIM.Pval.All.PotableReusePOCPOU.bray.Primers
# ADONIS2.GroupedP.All.PotableReusePOCPOU.bray["Primers"]<-ADONIS2.Pval.All.PotableReusePOCPOU.bray.Primers
#
# Betadis.GroupedP.All.PotableReusePOCPOU.bray["Primers"]<-Betadis.Pval.All.PotableReusePOCPOU.bray.Primers
#
# ANOSIM.GroupedR.All.PotableReusePOCPOU.bray["Primers"]<-ANOSIM.R.All.PotableReusePOCPOU.bray.Primers
# ADONIS2.GroupedR.All.PotableReusePOCPOU.bray["Primers"]<-ADONIS2.R2.All.PotableReusePOCPOU.bray.Primers
#
# #Pairwise if needed, check if interesting or not. Remove if not.
# count(All.PotableReusePOCPOU.meta$Primers)
#
# ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Primers<-pairwise.adonis2(All.PotableReusePOCPOU.dist.bray ~ Primers,
#                                                                      data = All.PotableReusePOCPOU.meta,
#                                                                      permutations = permutations)
# ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Primers
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
# ANOSIM.All.PotableReusePOCPOU.bray.Trimmed<-anosim(All.PotableReusePOCPOU.dist.bray,
#                                                  All.PotableReusePOCPOU.meta$Trimmed,
#                                                  permutations = permutations)
#
# ADONIS2.All.PotableReusePOCPOU.bray.Trimmed<-adonis2(All.PotableReusePOCPOU.dist.bray ~ Trimmed,
#                                                    data = All.PotableReusePOCPOU.meta,
#                                                    permutations = permutations)
#
# Betadis.All.PotableReusePOCPOU.bray.Trimmed<-anova(betadisper(All.PotableReusePOCPOU.dist.bray,
#                                                             All.PotableReusePOCPOU.meta$Trimmed))
#
#
#
# ANOSIM.Pval.All.PotableReusePOCPOU.bray.Trimmed<-ANOSIM.All.PotableReusePOCPOU.bray.Trimmed$signif
# ANOSIM.Pval.All.PotableReusePOCPOU.bray.Trimmed
#
# ADONIS2.Pval.All.PotableReusePOCPOU.bray.Trimmed<-ADONIS2.All.PotableReusePOCPOU.bray.Trimmed$`Pr(>F)`[1]
# ADONIS2.Pval.All.PotableReusePOCPOU.bray.Trimmed
#
# Betadis.Pval.All.PotableReusePOCPOU.bray.Trimmed<-Betadis.All.PotableReusePOCPOU.bray.Trimmed$`Pr(>F)`[1]
# Betadis.Pval.All.PotableReusePOCPOU.bray.Trimmed
#
# ANOSIM.R.All.PotableReusePOCPOU.bray.Trimmed<-ANOSIM.All.PotableReusePOCPOU.bray.Trimmed$statistic
# ANOSIM.R.All.PotableReusePOCPOU.bray.Trimmed
#
# ADONIS2.R2.All.PotableReusePOCPOU.bray.Trimmed<-ADONIS2.All.PotableReusePOCPOU.bray.Trimmed$R2[1]
# ADONIS2.R2.All.PotableReusePOCPOU.bray.Trimmed
#
# #Pull Together Pvalues to adjust
# ANOSIM.GroupedP.All.PotableReusePOCPOU.bray["Trimmed"]<-ANOSIM.Pval.All.PotableReusePOCPOU.bray.Trimmed
# ADONIS2.GroupedP.All.PotableReusePOCPOU.bray["Trimmed"]<-ADONIS2.Pval.All.PotableReusePOCPOU.bray.Trimmed
#
# Betadis.GroupedP.All.PotableReusePOCPOU.bray["Trimmed"]<-Betadis.Pval.All.PotableReusePOCPOU.bray.Trimmed
#
# ANOSIM.GroupedR.All.PotableReusePOCPOU.bray["Trimmed"]<-ANOSIM.R.All.PotableReusePOCPOU.bray.Trimmed
# ADONIS2.GroupedR.All.PotableReusePOCPOU.bray["Trimmed"]<-ADONIS2.R2.All.PotableReusePOCPOU.bray.Trimmed
#
# #Pairwise if needed, check if interesting or not. Remove if not.
# count(All.PotableReusePOCPOU.meta$Trimmed)
#
# ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Trimmed<-pairwise.adonis2(All.PotableReusePOCPOU.dist.bray ~ Trimmed,
#                                                                      data = All.PotableReusePOCPOU.meta,
#                                                                      permutations = permutations)
# ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Trimmed
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
# ANOSIM.All.PotableReusePOCPOU.bray.Extraction_Type<-anosim(All.PotableReusePOCPOU.dist.bray,
#                                                  All.PotableReusePOCPOU.meta$Extraction_Type,
#                                                  permutations = permutations)
#
# ADONIS2.All.PotableReusePOCPOU.bray.Extraction_Type<-adonis2(All.PotableReusePOCPOU.dist.bray ~ Extraction_Type,
#                                                    data = All.PotableReusePOCPOU.meta,
#                                                    permutations = permutations)
#
# Betadis.All.PotableReusePOCPOU.bray.Extraction_Type<-anova(betadisper(All.PotableReusePOCPOU.dist.bray,
#                                                             All.PotableReusePOCPOU.meta$Extraction_Type))
#
#
#
# ANOSIM.Pval.All.PotableReusePOCPOU.bray.Extraction_Type<-ANOSIM.All.PotableReusePOCPOU.bray.Extraction_Type$signif
# ANOSIM.Pval.All.PotableReusePOCPOU.bray.Extraction_Type
#
# ADONIS2.Pval.All.PotableReusePOCPOU.bray.Extraction_Type<-ADONIS2.All.PotableReusePOCPOU.bray.Extraction_Type$`Pr(>F)`[1]
# ADONIS2.Pval.All.PotableReusePOCPOU.bray.Extraction_Type
#
# Betadis.Pval.All.PotableReusePOCPOU.bray.Extraction_Type<-Betadis.All.PotableReusePOCPOU.bray.Extraction_Type$`Pr(>F)`[1]
# Betadis.Pval.All.PotableReusePOCPOU.bray.Extraction_Type
#
# ANOSIM.R.All.PotableReusePOCPOU.bray.Extraction_Type<-ANOSIM.All.PotableReusePOCPOU.bray.Extraction_Type$statistic
# ANOSIM.R.All.PotableReusePOCPOU.bray.Extraction_Type
#
# ADONIS2.R2.All.PotableReusePOCPOU.bray.Extraction_Type<-ADONIS2.All.PotableReusePOCPOU.bray.Extraction_Type$R2[1]
# ADONIS2.R2.All.PotableReusePOCPOU.bray.Extraction_Type
#
# #Pull Together Pvalues to adjust
# ANOSIM.GroupedP.All.PotableReusePOCPOU.bray["Extraction_Type"]<-ANOSIM.Pval.All.PotableReusePOCPOU.bray.Extraction_Type
# ADONIS2.GroupedP.All.PotableReusePOCPOU.bray["Extraction_Type"]<-ADONIS2.Pval.All.PotableReusePOCPOU.bray.Extraction_Type
#
# Betadis.GroupedP.All.PotableReusePOCPOU.bray["Extraction_Type"]<-Betadis.Pval.All.PotableReusePOCPOU.bray.Extraction_Type
#
# ANOSIM.GroupedR.All.PotableReusePOCPOU.bray["Extraction_Type"]<-ANOSIM.R.All.PotableReusePOCPOU.bray.Extraction_Type
# ADONIS2.GroupedR.All.PotableReusePOCPOU.bray["Extraction_Type"]<-ADONIS2.R2.All.PotableReusePOCPOU.bray.Extraction_Type
#
# #Pairwise if needed, check if interesting or not. Remove if not.
# count(All.PotableReusePOCPOU.meta$Extraction_Type)
#
# ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Extraction_Type<-pairwise.adonis2(All.PotableReusePOCPOU.dist.bray ~ Extraction_Type,
#                                                                      data = All.PotableReusePOCPOU.meta,
#                                                                      permutations = permutations)
# ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Extraction_Type
# #Extraction_Type





#################
#################
#################
#################
#Potable_Source_Water
ANOSIM.All.PotableReusePOCPOU.bray.Potable_Source_Water<-anosim(All.PotableReusePOCPOU.dist.bray,
                                                      All.PotableReusePOCPOU.meta$Potable_Source_Water,
                                                      permutations = permutations)

ADONIS2.All.PotableReusePOCPOU.bray.Potable_Source_Water<-adonis2(All.PotableReusePOCPOU.dist.bray ~ Potable_Source_Water,
                                                        data = All.PotableReusePOCPOU.meta,
                                                        permutations = permutations)

Betadis.All.PotableReusePOCPOU.bray.Potable_Source_Water<-anova(betadisper(All.PotableReusePOCPOU.dist.bray,
                                                                 All.PotableReusePOCPOU.meta$Potable_Source_Water))



ANOSIM.Pval.All.PotableReusePOCPOU.bray.Potable_Source_Water<-ANOSIM.All.PotableReusePOCPOU.bray.Potable_Source_Water$signif
ANOSIM.Pval.All.PotableReusePOCPOU.bray.Potable_Source_Water

ADONIS2.Pval.All.PotableReusePOCPOU.bray.Potable_Source_Water<-ADONIS2.All.PotableReusePOCPOU.bray.Potable_Source_Water$`Pr(>F)`[1]
ADONIS2.Pval.All.PotableReusePOCPOU.bray.Potable_Source_Water

Betadis.Pval.All.PotableReusePOCPOU.bray.Potable_Source_Water<-Betadis.All.PotableReusePOCPOU.bray.Potable_Source_Water$`Pr(>F)`[1]
Betadis.Pval.All.PotableReusePOCPOU.bray.Potable_Source_Water

ANOSIM.R.All.PotableReusePOCPOU.bray.Potable_Source_Water<-ANOSIM.All.PotableReusePOCPOU.bray.Potable_Source_Water$statistic
ANOSIM.R.All.PotableReusePOCPOU.bray.Potable_Source_Water

ADONIS2.R2.All.PotableReusePOCPOU.bray.Potable_Source_Water<-ADONIS2.All.PotableReusePOCPOU.bray.Potable_Source_Water$R2[1]
ADONIS2.R2.All.PotableReusePOCPOU.bray.Potable_Source_Water

#Pull Together Pvalues to adjust
ANOSIM.GroupedP.All.PotableReusePOCPOU.bray["Potable_Source_Water"]<-ANOSIM.Pval.All.PotableReusePOCPOU.bray.Potable_Source_Water
ADONIS2.GroupedP.All.PotableReusePOCPOU.bray["Potable_Source_Water"]<-ADONIS2.Pval.All.PotableReusePOCPOU.bray.Potable_Source_Water

Betadis.GroupedP.All.PotableReusePOCPOU.bray["Potable_Source_Water"]<-Betadis.Pval.All.PotableReusePOCPOU.bray.Potable_Source_Water

ANOSIM.GroupedR.All.PotableReusePOCPOU.bray["Potable_Source_Water"]<-ANOSIM.R.All.PotableReusePOCPOU.bray.Potable_Source_Water
ADONIS2.GroupedR.All.PotableReusePOCPOU.bray["Potable_Source_Water"]<-ADONIS2.R2.All.PotableReusePOCPOU.bray.Potable_Source_Water

#Pairwise if needed, check if interesting or not. Remove if not.
count(All.PotableReusePOCPOU.meta$Potable_Source_Water)

ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Potable_Source_Water<-pairwise.adonis2(All.PotableReusePOCPOU.dist.bray ~ Potable_Source_Water,
                                                                          data = All.PotableReusePOCPOU.meta,
                                                                          permutations = permutations)
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Potable_Source_Water
#Potable_Source_Water
##Potable_Source_Water




#################
#################
#################
#################
#Potable_Disinfection_Primary
ANOSIM.All.PotableReusePOCPOU.bray.Potable_Disinfection_Primary <-anosim(All.PotableReusePOCPOU.dist.bray,
                                                      All.PotableReusePOCPOU.meta$Potable_Disinfection_Primary ,
                                                      permutations = permutations)

ADONIS2.All.PotableReusePOCPOU.bray.Potable_Disinfection_Primary <-adonis2(All.PotableReusePOCPOU.dist.bray ~ Potable_Disinfection_Primary ,
                                                        data = All.PotableReusePOCPOU.meta,
                                                        permutations = permutations)

Betadis.All.PotableReusePOCPOU.bray.Potable_Disinfection_Primary <-anova(betadisper(All.PotableReusePOCPOU.dist.bray,
                                                                 All.PotableReusePOCPOU.meta$Potable_Disinfection_Primary ))



ANOSIM.Pval.All.PotableReusePOCPOU.bray.Potable_Disinfection_Primary <-ANOSIM.All.PotableReusePOCPOU.bray.Potable_Disinfection_Primary $signif
ANOSIM.Pval.All.PotableReusePOCPOU.bray.Potable_Disinfection_Primary

ADONIS2.Pval.All.PotableReusePOCPOU.bray.Potable_Disinfection_Primary <-ADONIS2.All.PotableReusePOCPOU.bray.Potable_Disinfection_Primary $`Pr(>F)`[1]
ADONIS2.Pval.All.PotableReusePOCPOU.bray.Potable_Disinfection_Primary

Betadis.Pval.All.PotableReusePOCPOU.bray.Potable_Disinfection_Primary <-Betadis.All.PotableReusePOCPOU.bray.Potable_Disinfection_Primary $`Pr(>F)`[1]
Betadis.Pval.All.PotableReusePOCPOU.bray.Potable_Disinfection_Primary

ANOSIM.R.All.PotableReusePOCPOU.bray.Potable_Disinfection_Primary <-ANOSIM.All.PotableReusePOCPOU.bray.Potable_Disinfection_Primary $statistic
ANOSIM.R.All.PotableReusePOCPOU.bray.Potable_Disinfection_Primary

ADONIS2.R2.All.PotableReusePOCPOU.bray.Potable_Disinfection_Primary <-ADONIS2.All.PotableReusePOCPOU.bray.Potable_Disinfection_Primary $R2[1]
ADONIS2.R2.All.PotableReusePOCPOU.bray.Potable_Disinfection_Primary

#Pull Together Pvalues to adjust
ANOSIM.GroupedP.All.PotableReusePOCPOU.bray["Potable_Disinfection_Primary "]<-ANOSIM.Pval.All.PotableReusePOCPOU.bray.Potable_Disinfection_Primary
ADONIS2.GroupedP.All.PotableReusePOCPOU.bray["Potable_Disinfection_Primary "]<-ADONIS2.Pval.All.PotableReusePOCPOU.bray.Potable_Disinfection_Primary

Betadis.GroupedP.All.PotableReusePOCPOU.bray["Potable_Disinfection_Primary "]<-Betadis.Pval.All.PotableReusePOCPOU.bray.Potable_Disinfection_Primary

ANOSIM.GroupedR.All.PotableReusePOCPOU.bray["Potable_Disinfection_Primary "]<-ANOSIM.R.All.PotableReusePOCPOU.bray.Potable_Disinfection_Primary
ADONIS2.GroupedR.All.PotableReusePOCPOU.bray["Potable_Disinfection_Primary "]<-ADONIS2.R2.All.PotableReusePOCPOU.bray.Potable_Disinfection_Primary

#Pairwise if needed, check if interesting or not. Remove if not.
count(All.PotableReusePOCPOU.meta$Potable_Disinfection_Primary )

ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Potable_Disinfection_Primary <-pairwise.adonis2(All.PotableReusePOCPOU.dist.bray ~ Potable_Disinfection_Primary ,
                                                                          data = All.PotableReusePOCPOU.meta,
                                                                          permutations = permutations)
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Potable_Disinfection_Primary
#Potable_Disinfection_Primary
##Potable_Disinfection_Primary




#################
#################
#################
#################
#Potable_Treatment_Type
ANOSIM.All.PotableReusePOCPOU.bray.Potable_Treatment_Type<-anosim(All.PotableReusePOCPOU.dist.bray,
                                                      All.PotableReusePOCPOU.meta$Potable_Treatment_Type,
                                                      permutations = permutations)

ADONIS2.All.PotableReusePOCPOU.bray.Potable_Treatment_Type<-adonis2(All.PotableReusePOCPOU.dist.bray ~ Potable_Treatment_Type,
                                                        data = All.PotableReusePOCPOU.meta,
                                                        permutations = permutations)

Betadis.All.PotableReusePOCPOU.bray.Potable_Treatment_Type<-anova(betadisper(All.PotableReusePOCPOU.dist.bray,
                                                                 All.PotableReusePOCPOU.meta$Potable_Treatment_Type))



ANOSIM.Pval.All.PotableReusePOCPOU.bray.Potable_Treatment_Type<-ANOSIM.All.PotableReusePOCPOU.bray.Potable_Treatment_Type$signif
ANOSIM.Pval.All.PotableReusePOCPOU.bray.Potable_Treatment_Type

ADONIS2.Pval.All.PotableReusePOCPOU.bray.Potable_Treatment_Type<-ADONIS2.All.PotableReusePOCPOU.bray.Potable_Treatment_Type$`Pr(>F)`[1]
ADONIS2.Pval.All.PotableReusePOCPOU.bray.Potable_Treatment_Type

Betadis.Pval.All.PotableReusePOCPOU.bray.Potable_Treatment_Type<-Betadis.All.PotableReusePOCPOU.bray.Potable_Treatment_Type$`Pr(>F)`[1]
Betadis.Pval.All.PotableReusePOCPOU.bray.Potable_Treatment_Type

ANOSIM.R.All.PotableReusePOCPOU.bray.Potable_Treatment_Type<-ANOSIM.All.PotableReusePOCPOU.bray.Potable_Treatment_Type$statistic
ANOSIM.R.All.PotableReusePOCPOU.bray.Potable_Treatment_Type

ADONIS2.R2.All.PotableReusePOCPOU.bray.Potable_Treatment_Type<-ADONIS2.All.PotableReusePOCPOU.bray.Potable_Treatment_Type$R2[1]
ADONIS2.R2.All.PotableReusePOCPOU.bray.Potable_Treatment_Type

#Pull Together Pvalues to adjust
ANOSIM.GroupedP.All.PotableReusePOCPOU.bray["Potable_Treatment_Type"]<-ANOSIM.Pval.All.PotableReusePOCPOU.bray.Potable_Treatment_Type
ADONIS2.GroupedP.All.PotableReusePOCPOU.bray["Potable_Treatment_Type"]<-ADONIS2.Pval.All.PotableReusePOCPOU.bray.Potable_Treatment_Type

Betadis.GroupedP.All.PotableReusePOCPOU.bray["Potable_Treatment_Type"]<-Betadis.Pval.All.PotableReusePOCPOU.bray.Potable_Treatment_Type

ANOSIM.GroupedR.All.PotableReusePOCPOU.bray["Potable_Treatment_Type"]<-ANOSIM.R.All.PotableReusePOCPOU.bray.Potable_Treatment_Type
ADONIS2.GroupedR.All.PotableReusePOCPOU.bray["Potable_Treatment_Type"]<-ADONIS2.R2.All.PotableReusePOCPOU.bray.Potable_Treatment_Type

#Pairwise if needed, check if interesting or not. Remove if not.
count(All.PotableReusePOCPOU.meta$Potable_Treatment_Type)

ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Potable_Treatment_Type<-pairwise.adonis2(All.PotableReusePOCPOU.dist.bray ~ Potable_Treatment_Type,
                                                                          data = All.PotableReusePOCPOU.meta,
                                                                          permutations = permutations)
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Potable_Treatment_Type
#Potable_Treatment_Type
##Potable_Treatment_Type




#################
#################
#################
#################
#Disinfection_Total
ANOSIM.All.PotableReusePOCPOU.bray.Disinfection_Total<-anosim(All.PotableReusePOCPOU.dist.bray,
                                                      All.PotableReusePOCPOU.meta$Disinfection_Total,
                                                      permutations = permutations)

ADONIS2.All.PotableReusePOCPOU.bray.Disinfection_Total<-adonis2(All.PotableReusePOCPOU.dist.bray ~ Disinfection_Total,
                                                        data = All.PotableReusePOCPOU.meta,
                                                        permutations = permutations)

Betadis.All.PotableReusePOCPOU.bray.Disinfection_Total<-anova(betadisper(All.PotableReusePOCPOU.dist.bray,
                                                                 All.PotableReusePOCPOU.meta$Disinfection_Total))



ANOSIM.Pval.All.PotableReusePOCPOU.bray.Disinfection_Total<-ANOSIM.All.PotableReusePOCPOU.bray.Disinfection_Total$signif
ANOSIM.Pval.All.PotableReusePOCPOU.bray.Disinfection_Total

ADONIS2.Pval.All.PotableReusePOCPOU.bray.Disinfection_Total<-ADONIS2.All.PotableReusePOCPOU.bray.Disinfection_Total$`Pr(>F)`[1]
ADONIS2.Pval.All.PotableReusePOCPOU.bray.Disinfection_Total

Betadis.Pval.All.PotableReusePOCPOU.bray.Disinfection_Total<-Betadis.All.PotableReusePOCPOU.bray.Disinfection_Total$`Pr(>F)`[1]
Betadis.Pval.All.PotableReusePOCPOU.bray.Disinfection_Total

ANOSIM.R.All.PotableReusePOCPOU.bray.Disinfection_Total<-ANOSIM.All.PotableReusePOCPOU.bray.Disinfection_Total$statistic
ANOSIM.R.All.PotableReusePOCPOU.bray.Disinfection_Total

ADONIS2.R2.All.PotableReusePOCPOU.bray.Disinfection_Total<-ADONIS2.All.PotableReusePOCPOU.bray.Disinfection_Total$R2[1]
ADONIS2.R2.All.PotableReusePOCPOU.bray.Disinfection_Total

#Pull Together Pvalues to adjust
ANOSIM.GroupedP.All.PotableReusePOCPOU.bray["Disinfection_Total"]<-ANOSIM.Pval.All.PotableReusePOCPOU.bray.Disinfection_Total
ADONIS2.GroupedP.All.PotableReusePOCPOU.bray["Disinfection_Total"]<-ADONIS2.Pval.All.PotableReusePOCPOU.bray.Disinfection_Total

Betadis.GroupedP.All.PotableReusePOCPOU.bray["Disinfection_Total"]<-Betadis.Pval.All.PotableReusePOCPOU.bray.Disinfection_Total

ANOSIM.GroupedR.All.PotableReusePOCPOU.bray["Disinfection_Total"]<-ANOSIM.R.All.PotableReusePOCPOU.bray.Disinfection_Total
ADONIS2.GroupedR.All.PotableReusePOCPOU.bray["Disinfection_Total"]<-ADONIS2.R2.All.PotableReusePOCPOU.bray.Disinfection_Total

#Pairwise if needed, check if interesting or not. Remove if not.
count(All.PotableReusePOCPOU.meta$Disinfection_Total)

ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Disinfection_Total<-pairwise.adonis2(All.PotableReusePOCPOU.dist.bray ~ Disinfection_Total,
                                                                          data = All.PotableReusePOCPOU.meta,
                                                                          permutations = permutations)
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Disinfection_Total
#Disinfection_Total
##Disinfection_Total




#################
#################
#################
#################
#Treatment_Total
ANOSIM.All.PotableReusePOCPOU.bray.Treatment_Total <-anosim(All.PotableReusePOCPOU.dist.bray,
                                                      All.PotableReusePOCPOU.meta$Treatment_Total ,
                                                      permutations = permutations)

ADONIS2.All.PotableReusePOCPOU.bray.Treatment_Total <-adonis2(All.PotableReusePOCPOU.dist.bray ~ Treatment_Total ,
                                                        data = All.PotableReusePOCPOU.meta,
                                                        permutations = permutations)

Betadis.All.PotableReusePOCPOU.bray.Treatment_Total <-anova(betadisper(All.PotableReusePOCPOU.dist.bray,
                                                                 All.PotableReusePOCPOU.meta$Treatment_Total ))



ANOSIM.Pval.All.PotableReusePOCPOU.bray.Treatment_Total <-ANOSIM.All.PotableReusePOCPOU.bray.Treatment_Total $signif
ANOSIM.Pval.All.PotableReusePOCPOU.bray.Treatment_Total

ADONIS2.Pval.All.PotableReusePOCPOU.bray.Treatment_Total <-ADONIS2.All.PotableReusePOCPOU.bray.Treatment_Total $`Pr(>F)`[1]
ADONIS2.Pval.All.PotableReusePOCPOU.bray.Treatment_Total

Betadis.Pval.All.PotableReusePOCPOU.bray.Treatment_Total <-Betadis.All.PotableReusePOCPOU.bray.Treatment_Total $`Pr(>F)`[1]
Betadis.Pval.All.PotableReusePOCPOU.bray.Treatment_Total

ANOSIM.R.All.PotableReusePOCPOU.bray.Treatment_Total <-ANOSIM.All.PotableReusePOCPOU.bray.Treatment_Total $statistic
ANOSIM.R.All.PotableReusePOCPOU.bray.Treatment_Total

ADONIS2.R2.All.PotableReusePOCPOU.bray.Treatment_Total <-ADONIS2.All.PotableReusePOCPOU.bray.Treatment_Total $R2[1]
ADONIS2.R2.All.PotableReusePOCPOU.bray.Treatment_Total

#Pull Together Pvalues to adjust
ANOSIM.GroupedP.All.PotableReusePOCPOU.bray["Treatment_Total "]<-ANOSIM.Pval.All.PotableReusePOCPOU.bray.Treatment_Total
ADONIS2.GroupedP.All.PotableReusePOCPOU.bray["Treatment_Total "]<-ADONIS2.Pval.All.PotableReusePOCPOU.bray.Treatment_Total

Betadis.GroupedP.All.PotableReusePOCPOU.bray["Treatment_Total "]<-Betadis.Pval.All.PotableReusePOCPOU.bray.Treatment_Total

ANOSIM.GroupedR.All.PotableReusePOCPOU.bray["Treatment_Total "]<-ANOSIM.R.All.PotableReusePOCPOU.bray.Treatment_Total
ADONIS2.GroupedR.All.PotableReusePOCPOU.bray["Treatment_Total "]<-ADONIS2.R2.All.PotableReusePOCPOU.bray.Treatment_Total

#Pairwise if needed, check if interesting or not. Remove if not.
count(All.PotableReusePOCPOU.meta$Treatment_Total )

ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Treatment_Total <-pairwise.adonis2(All.PotableReusePOCPOU.dist.bray ~ Treatment_Total ,
                                                                          data = All.PotableReusePOCPOU.meta,
                                                                          permutations = permutations)
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Treatment_Total
#Treatment_Total
##Treatment_Total




#################
#################
#################
#################
#Treatment_Classification_Number
ANOSIM.All.PotableReusePOCPOU.bray.Treatment_Classification_Number<-anosim(All.PotableReusePOCPOU.dist.bray,
                                                      All.PotableReusePOCPOU.meta$Treatment_Classification_Number,
                                                      permutations = permutations)

ADONIS2.All.PotableReusePOCPOU.bray.Treatment_Classification_Number<-adonis2(All.PotableReusePOCPOU.dist.bray ~ Treatment_Classification_Number,
                                                        data = All.PotableReusePOCPOU.meta,
                                                        permutations = permutations)

Betadis.All.PotableReusePOCPOU.bray.Treatment_Classification_Number<-anova(betadisper(All.PotableReusePOCPOU.dist.bray,
                                                                 All.PotableReusePOCPOU.meta$Treatment_Classification_Number))



ANOSIM.Pval.All.PotableReusePOCPOU.bray.Treatment_Classification_Number<-ANOSIM.All.PotableReusePOCPOU.bray.Treatment_Classification_Number$signif
ANOSIM.Pval.All.PotableReusePOCPOU.bray.Treatment_Classification_Number

ADONIS2.Pval.All.PotableReusePOCPOU.bray.Treatment_Classification_Number<-ADONIS2.All.PotableReusePOCPOU.bray.Treatment_Classification_Number$`Pr(>F)`[1]
ADONIS2.Pval.All.PotableReusePOCPOU.bray.Treatment_Classification_Number

Betadis.Pval.All.PotableReusePOCPOU.bray.Treatment_Classification_Number<-Betadis.All.PotableReusePOCPOU.bray.Treatment_Classification_Number$`Pr(>F)`[1]
Betadis.Pval.All.PotableReusePOCPOU.bray.Treatment_Classification_Number

ANOSIM.R.All.PotableReusePOCPOU.bray.Treatment_Classification_Number<-ANOSIM.All.PotableReusePOCPOU.bray.Treatment_Classification_Number$statistic
ANOSIM.R.All.PotableReusePOCPOU.bray.Treatment_Classification_Number

ADONIS2.R2.All.PotableReusePOCPOU.bray.Treatment_Classification_Number<-ADONIS2.All.PotableReusePOCPOU.bray.Treatment_Classification_Number$R2[1]
ADONIS2.R2.All.PotableReusePOCPOU.bray.Treatment_Classification_Number

#Pull Together Pvalues to adjust
ANOSIM.GroupedP.All.PotableReusePOCPOU.bray["Treatment_Classification_Number"]<-ANOSIM.Pval.All.PotableReusePOCPOU.bray.Treatment_Classification_Number
ADONIS2.GroupedP.All.PotableReusePOCPOU.bray["Treatment_Classification_Number"]<-ADONIS2.Pval.All.PotableReusePOCPOU.bray.Treatment_Classification_Number

Betadis.GroupedP.All.PotableReusePOCPOU.bray["Treatment_Classification_Number"]<-Betadis.Pval.All.PotableReusePOCPOU.bray.Treatment_Classification_Number

ANOSIM.GroupedR.All.PotableReusePOCPOU.bray["Treatment_Classification_Number"]<-ANOSIM.R.All.PotableReusePOCPOU.bray.Treatment_Classification_Number
ADONIS2.GroupedR.All.PotableReusePOCPOU.bray["Treatment_Classification_Number"]<-ADONIS2.R2.All.PotableReusePOCPOU.bray.Treatment_Classification_Number

#Pairwise if needed, check if interesting or not. Remove if not.
count(All.PotableReusePOCPOU.meta$Treatment_Classification_Number)

ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Treatment_Classification_Number<-pairwise.adonis2(All.PotableReusePOCPOU.dist.bray ~ Treatment_Classification_Number,
                                                                          data = All.PotableReusePOCPOU.meta,
                                                                          permutations = permutations)
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Treatment_Classification_Number
#Treatment_Classification_Number
##Treatment_Classification_Number




#################
#################
#################
#################
#Treatment_Classification_Cat
ANOSIM.All.PotableReusePOCPOU.bray.Treatment_Classification_Cat<-anosim(All.PotableReusePOCPOU.dist.bray,
                                                      All.PotableReusePOCPOU.meta$Treatment_Classification_Cat,
                                                      permutations = permutations)

ADONIS2.All.PotableReusePOCPOU.bray.Treatment_Classification_Cat<-adonis2(All.PotableReusePOCPOU.dist.bray ~ Treatment_Classification_Cat,
                                                        data = All.PotableReusePOCPOU.meta,
                                                        permutations = permutations)

Betadis.All.PotableReusePOCPOU.bray.Treatment_Classification_Cat<-anova(betadisper(All.PotableReusePOCPOU.dist.bray,
                                                                 All.PotableReusePOCPOU.meta$Treatment_Classification_Cat))



ANOSIM.Pval.All.PotableReusePOCPOU.bray.Treatment_Classification_Cat<-ANOSIM.All.PotableReusePOCPOU.bray.Treatment_Classification_Cat$signif
ANOSIM.Pval.All.PotableReusePOCPOU.bray.Treatment_Classification_Cat

ADONIS2.Pval.All.PotableReusePOCPOU.bray.Treatment_Classification_Cat<-ADONIS2.All.PotableReusePOCPOU.bray.Treatment_Classification_Cat$`Pr(>F)`[1]
ADONIS2.Pval.All.PotableReusePOCPOU.bray.Treatment_Classification_Cat

Betadis.Pval.All.PotableReusePOCPOU.bray.Treatment_Classification_Cat<-Betadis.All.PotableReusePOCPOU.bray.Treatment_Classification_Cat$`Pr(>F)`[1]
Betadis.Pval.All.PotableReusePOCPOU.bray.Treatment_Classification_Cat

ANOSIM.R.All.PotableReusePOCPOU.bray.Treatment_Classification_Cat<-ANOSIM.All.PotableReusePOCPOU.bray.Treatment_Classification_Cat$statistic
ANOSIM.R.All.PotableReusePOCPOU.bray.Treatment_Classification_Cat

ADONIS2.R2.All.PotableReusePOCPOU.bray.Treatment_Classification_Cat<-ADONIS2.All.PotableReusePOCPOU.bray.Treatment_Classification_Cat$R2[1]
ADONIS2.R2.All.PotableReusePOCPOU.bray.Treatment_Classification_Cat

#Pull Together Pvalues to adjust
ANOSIM.GroupedP.All.PotableReusePOCPOU.bray["Treatment_Classification_Cat"]<-ANOSIM.Pval.All.PotableReusePOCPOU.bray.Treatment_Classification_Cat
ADONIS2.GroupedP.All.PotableReusePOCPOU.bray["Treatment_Classification_Cat"]<-ADONIS2.Pval.All.PotableReusePOCPOU.bray.Treatment_Classification_Cat

Betadis.GroupedP.All.PotableReusePOCPOU.bray["Treatment_Classification_Cat"]<-Betadis.Pval.All.PotableReusePOCPOU.bray.Treatment_Classification_Cat

ANOSIM.GroupedR.All.PotableReusePOCPOU.bray["Treatment_Classification_Cat"]<-ANOSIM.R.All.PotableReusePOCPOU.bray.Treatment_Classification_Cat
ADONIS2.GroupedR.All.PotableReusePOCPOU.bray["Treatment_Classification_Cat"]<-ADONIS2.R2.All.PotableReusePOCPOU.bray.Treatment_Classification_Cat

#Pairwise if needed, check if interesting or not. Remove if not.
count(All.PotableReusePOCPOU.meta$Treatment_Classification_Cat)

ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Treatment_Classification_Cat<-pairwise.adonis2(All.PotableReusePOCPOU.dist.bray ~ Treatment_Classification_Cat,
                                                                          data = All.PotableReusePOCPOU.meta,
                                                                          permutations = permutations)
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Treatment_Classification_Cat
#Treatment_Classification_Cat
##Treatment_Classification_Cat




#################
#################
#################
#################
#Mixed_Source_Water
ANOSIM.All.PotableReusePOCPOU.bray.Mixed_Source_Water<-anosim(All.PotableReusePOCPOU.dist.bray,
                                                      All.PotableReusePOCPOU.meta$Mixed_Source_Water,
                                                      permutations = permutations)

ADONIS2.All.PotableReusePOCPOU.bray.Mixed_Source_Water<-adonis2(All.PotableReusePOCPOU.dist.bray ~ Mixed_Source_Water,
                                                        data = All.PotableReusePOCPOU.meta,
                                                        permutations = permutations)

Betadis.All.PotableReusePOCPOU.bray.Mixed_Source_Water<-anova(betadisper(All.PotableReusePOCPOU.dist.bray,
                                                                 All.PotableReusePOCPOU.meta$Mixed_Source_Water))



ANOSIM.Pval.All.PotableReusePOCPOU.bray.Mixed_Source_Water<-ANOSIM.All.PotableReusePOCPOU.bray.Mixed_Source_Water$signif
ANOSIM.Pval.All.PotableReusePOCPOU.bray.Mixed_Source_Water

ADONIS2.Pval.All.PotableReusePOCPOU.bray.Mixed_Source_Water<-ADONIS2.All.PotableReusePOCPOU.bray.Mixed_Source_Water$`Pr(>F)`[1]
ADONIS2.Pval.All.PotableReusePOCPOU.bray.Mixed_Source_Water

Betadis.Pval.All.PotableReusePOCPOU.bray.Mixed_Source_Water<-Betadis.All.PotableReusePOCPOU.bray.Mixed_Source_Water$`Pr(>F)`[1]
Betadis.Pval.All.PotableReusePOCPOU.bray.Mixed_Source_Water

ANOSIM.R.All.PotableReusePOCPOU.bray.Mixed_Source_Water<-ANOSIM.All.PotableReusePOCPOU.bray.Mixed_Source_Water$statistic
ANOSIM.R.All.PotableReusePOCPOU.bray.Mixed_Source_Water

ADONIS2.R2.All.PotableReusePOCPOU.bray.Mixed_Source_Water<-ADONIS2.All.PotableReusePOCPOU.bray.Mixed_Source_Water$R2[1]
ADONIS2.R2.All.PotableReusePOCPOU.bray.Mixed_Source_Water

#Pull Together Pvalues to adjust
ANOSIM.GroupedP.All.PotableReusePOCPOU.bray["Mixed_Source_Water"]<-ANOSIM.Pval.All.PotableReusePOCPOU.bray.Mixed_Source_Water
ADONIS2.GroupedP.All.PotableReusePOCPOU.bray["Mixed_Source_Water"]<-ADONIS2.Pval.All.PotableReusePOCPOU.bray.Mixed_Source_Water

Betadis.GroupedP.All.PotableReusePOCPOU.bray["Mixed_Source_Water"]<-Betadis.Pval.All.PotableReusePOCPOU.bray.Mixed_Source_Water

ANOSIM.GroupedR.All.PotableReusePOCPOU.bray["Mixed_Source_Water"]<-ANOSIM.R.All.PotableReusePOCPOU.bray.Mixed_Source_Water
ADONIS2.GroupedR.All.PotableReusePOCPOU.bray["Mixed_Source_Water"]<-ADONIS2.R2.All.PotableReusePOCPOU.bray.Mixed_Source_Water

#Pairwise if needed, check if interesting or not. Remove if not.
count(All.PotableReusePOCPOU.meta$Mixed_Source_Water)

ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Mixed_Source_Water<-pairwise.adonis2(All.PotableReusePOCPOU.dist.bray ~ Mixed_Source_Water,
                                                                          data = All.PotableReusePOCPOU.meta,
                                                                          permutations = permutations)
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Mixed_Source_Water
#Mixed_Source_Water
##Mixed_Source_Water




#################
#################
#################
#################
#Potable_Overlap
ANOSIM.All.PotableReusePOCPOU.bray.Potable_Overlap<-anosim(All.PotableReusePOCPOU.dist.bray,
                                                      All.PotableReusePOCPOU.meta$Potable_Overlap,
                                                      permutations = permutations)

ADONIS2.All.PotableReusePOCPOU.bray.Potable_Overlap<-adonis2(All.PotableReusePOCPOU.dist.bray ~ Potable_Overlap,
                                                        data = All.PotableReusePOCPOU.meta,
                                                        permutations = permutations)

Betadis.All.PotableReusePOCPOU.bray.Potable_Overlap<-anova(betadisper(All.PotableReusePOCPOU.dist.bray,
                                                                 All.PotableReusePOCPOU.meta$Potable_Overlap))



ANOSIM.Pval.All.PotableReusePOCPOU.bray.Potable_Overlap<-ANOSIM.All.PotableReusePOCPOU.bray.Potable_Overlap$signif
ANOSIM.Pval.All.PotableReusePOCPOU.bray.Potable_Overlap

ADONIS2.Pval.All.PotableReusePOCPOU.bray.Potable_Overlap<-ADONIS2.All.PotableReusePOCPOU.bray.Potable_Overlap$`Pr(>F)`[1]
ADONIS2.Pval.All.PotableReusePOCPOU.bray.Potable_Overlap

Betadis.Pval.All.PotableReusePOCPOU.bray.Potable_Overlap<-Betadis.All.PotableReusePOCPOU.bray.Potable_Overlap$`Pr(>F)`[1]
Betadis.Pval.All.PotableReusePOCPOU.bray.Potable_Overlap

ANOSIM.R.All.PotableReusePOCPOU.bray.Potable_Overlap<-ANOSIM.All.PotableReusePOCPOU.bray.Potable_Overlap$statistic
ANOSIM.R.All.PotableReusePOCPOU.bray.Potable_Overlap

ADONIS2.R2.All.PotableReusePOCPOU.bray.Potable_Overlap<-ADONIS2.All.PotableReusePOCPOU.bray.Potable_Overlap$R2[1]
ADONIS2.R2.All.PotableReusePOCPOU.bray.Potable_Overlap

#Pull Together Pvalues to adjust
ANOSIM.GroupedP.All.PotableReusePOCPOU.bray["Potable_Overlap"]<-ANOSIM.Pval.All.PotableReusePOCPOU.bray.Potable_Overlap
ADONIS2.GroupedP.All.PotableReusePOCPOU.bray["Potable_Overlap"]<-ADONIS2.Pval.All.PotableReusePOCPOU.bray.Potable_Overlap

Betadis.GroupedP.All.PotableReusePOCPOU.bray["Potable_Overlap"]<-Betadis.Pval.All.PotableReusePOCPOU.bray.Potable_Overlap

ANOSIM.GroupedR.All.PotableReusePOCPOU.bray["Potable_Overlap"]<-ANOSIM.R.All.PotableReusePOCPOU.bray.Potable_Overlap
ADONIS2.GroupedR.All.PotableReusePOCPOU.bray["Potable_Overlap"]<-ADONIS2.R2.All.PotableReusePOCPOU.bray.Potable_Overlap

#Pairwise if needed, check if interesting or not. Remove if not.
count(All.PotableReusePOCPOU.meta$Potable_Overlap)

ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Potable_Overlap<-pairwise.adonis2(All.PotableReusePOCPOU.dist.bray ~ Potable_Overlap,
                                                                          data = All.PotableReusePOCPOU.meta,
                                                                          permutations = permutations)
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Potable_Overlap
#Potable_Overlap
##Potable_Overlap




#################
#################
#################
#################
#Includes_Ozone
ANOSIM.All.PotableReusePOCPOU.bray.Includes_Ozone<-anosim(All.PotableReusePOCPOU.dist.bray,
                                                      All.PotableReusePOCPOU.meta$Includes_Ozone,
                                                      permutations = permutations)

ADONIS2.All.PotableReusePOCPOU.bray.Includes_Ozone<-adonis2(All.PotableReusePOCPOU.dist.bray ~ Includes_Ozone,
                                                        data = All.PotableReusePOCPOU.meta,
                                                        permutations = permutations)

Betadis.All.PotableReusePOCPOU.bray.Includes_Ozone<-anova(betadisper(All.PotableReusePOCPOU.dist.bray,
                                                                 All.PotableReusePOCPOU.meta$Includes_Ozone))



ANOSIM.Pval.All.PotableReusePOCPOU.bray.Includes_Ozone<-ANOSIM.All.PotableReusePOCPOU.bray.Includes_Ozone$signif
ANOSIM.Pval.All.PotableReusePOCPOU.bray.Includes_Ozone

ADONIS2.Pval.All.PotableReusePOCPOU.bray.Includes_Ozone<-ADONIS2.All.PotableReusePOCPOU.bray.Includes_Ozone$`Pr(>F)`[1]
ADONIS2.Pval.All.PotableReusePOCPOU.bray.Includes_Ozone

Betadis.Pval.All.PotableReusePOCPOU.bray.Includes_Ozone<-Betadis.All.PotableReusePOCPOU.bray.Includes_Ozone$`Pr(>F)`[1]
Betadis.Pval.All.PotableReusePOCPOU.bray.Includes_Ozone

ANOSIM.R.All.PotableReusePOCPOU.bray.Includes_Ozone<-ANOSIM.All.PotableReusePOCPOU.bray.Includes_Ozone$statistic
ANOSIM.R.All.PotableReusePOCPOU.bray.Includes_Ozone

ADONIS2.R2.All.PotableReusePOCPOU.bray.Includes_Ozone<-ADONIS2.All.PotableReusePOCPOU.bray.Includes_Ozone$R2[1]
ADONIS2.R2.All.PotableReusePOCPOU.bray.Includes_Ozone

#Pull Together Pvalues to adjust
ANOSIM.GroupedP.All.PotableReusePOCPOU.bray["Includes_Ozone"]<-ANOSIM.Pval.All.PotableReusePOCPOU.bray.Includes_Ozone
ADONIS2.GroupedP.All.PotableReusePOCPOU.bray["Includes_Ozone"]<-ADONIS2.Pval.All.PotableReusePOCPOU.bray.Includes_Ozone

Betadis.GroupedP.All.PotableReusePOCPOU.bray["Includes_Ozone"]<-Betadis.Pval.All.PotableReusePOCPOU.bray.Includes_Ozone

ANOSIM.GroupedR.All.PotableReusePOCPOU.bray["Includes_Ozone"]<-ANOSIM.R.All.PotableReusePOCPOU.bray.Includes_Ozone
ADONIS2.GroupedR.All.PotableReusePOCPOU.bray["Includes_Ozone"]<-ADONIS2.R2.All.PotableReusePOCPOU.bray.Includes_Ozone

#Pairwise if needed, check if interesting or not. Remove if not.
count(All.PotableReusePOCPOU.meta$Includes_Ozone)

ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Includes_Ozone<-pairwise.adonis2(All.PotableReusePOCPOU.dist.bray ~ Includes_Ozone,
                                                                          data = All.PotableReusePOCPOU.meta,
                                                                          permutations = permutations)
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Includes_Ozone
#Includes_Ozone
##Includes_Ozone




#################
#################
#################
#################
#Includes_RO_UF
ANOSIM.All.PotableReusePOCPOU.bray.Includes_RO_UF<-anosim(All.PotableReusePOCPOU.dist.bray,
                                                      All.PotableReusePOCPOU.meta$Includes_RO_UF,
                                                      permutations = permutations)

ADONIS2.All.PotableReusePOCPOU.bray.Includes_RO_UF<-adonis2(All.PotableReusePOCPOU.dist.bray ~ Includes_RO_UF,
                                                        data = All.PotableReusePOCPOU.meta,
                                                        permutations = permutations)

Betadis.All.PotableReusePOCPOU.bray.Includes_RO_UF<-anova(betadisper(All.PotableReusePOCPOU.dist.bray,
                                                                 All.PotableReusePOCPOU.meta$Includes_RO_UF))



ANOSIM.Pval.All.PotableReusePOCPOU.bray.Includes_RO_UF<-ANOSIM.All.PotableReusePOCPOU.bray.Includes_RO_UF$signif
ANOSIM.Pval.All.PotableReusePOCPOU.bray.Includes_RO_UF

ADONIS2.Pval.All.PotableReusePOCPOU.bray.Includes_RO_UF<-ADONIS2.All.PotableReusePOCPOU.bray.Includes_RO_UF$`Pr(>F)`[1]
ADONIS2.Pval.All.PotableReusePOCPOU.bray.Includes_RO_UF

Betadis.Pval.All.PotableReusePOCPOU.bray.Includes_RO_UF<-Betadis.All.PotableReusePOCPOU.bray.Includes_RO_UF$`Pr(>F)`[1]
Betadis.Pval.All.PotableReusePOCPOU.bray.Includes_RO_UF

ANOSIM.R.All.PotableReusePOCPOU.bray.Includes_RO_UF<-ANOSIM.All.PotableReusePOCPOU.bray.Includes_RO_UF$statistic
ANOSIM.R.All.PotableReusePOCPOU.bray.Includes_RO_UF

ADONIS2.R2.All.PotableReusePOCPOU.bray.Includes_RO_UF<-ADONIS2.All.PotableReusePOCPOU.bray.Includes_RO_UF$R2[1]
ADONIS2.R2.All.PotableReusePOCPOU.bray.Includes_RO_UF

#Pull Together Pvalues to adjust
ANOSIM.GroupedP.All.PotableReusePOCPOU.bray["Includes_RO_UF"]<-ANOSIM.Pval.All.PotableReusePOCPOU.bray.Includes_RO_UF
ADONIS2.GroupedP.All.PotableReusePOCPOU.bray["Includes_RO_UF"]<-ADONIS2.Pval.All.PotableReusePOCPOU.bray.Includes_RO_UF

Betadis.GroupedP.All.PotableReusePOCPOU.bray["Includes_RO_UF"]<-Betadis.Pval.All.PotableReusePOCPOU.bray.Includes_RO_UF

ANOSIM.GroupedR.All.PotableReusePOCPOU.bray["Includes_RO_UF"]<-ANOSIM.R.All.PotableReusePOCPOU.bray.Includes_RO_UF
ADONIS2.GroupedR.All.PotableReusePOCPOU.bray["Includes_RO_UF"]<-ADONIS2.R2.All.PotableReusePOCPOU.bray.Includes_RO_UF

#Pairwise if needed, check if interesting or not. Remove if not.
count(All.PotableReusePOCPOU.meta$Includes_RO_UF)

ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Includes_RO_UF<-pairwise.adonis2(All.PotableReusePOCPOU.dist.bray ~ Includes_RO_UF,
                                                                          data = All.PotableReusePOCPOU.meta,
                                                                          permutations = permutations)
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Includes_RO_UF
#Includes_RO_UF
##Includes_RO_UF




#################
#################
#################
#################
#Includes_UV
ANOSIM.All.PotableReusePOCPOU.bray.Includes_UV<-anosim(All.PotableReusePOCPOU.dist.bray,
                                                      All.PotableReusePOCPOU.meta$Includes_UV,
                                                      permutations = permutations)

ADONIS2.All.PotableReusePOCPOU.bray.Includes_UV<-adonis2(All.PotableReusePOCPOU.dist.bray ~ Includes_UV,
                                                        data = All.PotableReusePOCPOU.meta,
                                                        permutations = permutations)

Betadis.All.PotableReusePOCPOU.bray.Includes_UV<-anova(betadisper(All.PotableReusePOCPOU.dist.bray,
                                                                 All.PotableReusePOCPOU.meta$Includes_UV))



ANOSIM.Pval.All.PotableReusePOCPOU.bray.Includes_UV<-ANOSIM.All.PotableReusePOCPOU.bray.Includes_UV$signif
ANOSIM.Pval.All.PotableReusePOCPOU.bray.Includes_UV

ADONIS2.Pval.All.PotableReusePOCPOU.bray.Includes_UV<-ADONIS2.All.PotableReusePOCPOU.bray.Includes_UV$`Pr(>F)`[1]
ADONIS2.Pval.All.PotableReusePOCPOU.bray.Includes_UV

Betadis.Pval.All.PotableReusePOCPOU.bray.Includes_UV<-Betadis.All.PotableReusePOCPOU.bray.Includes_UV$`Pr(>F)`[1]
Betadis.Pval.All.PotableReusePOCPOU.bray.Includes_UV

ANOSIM.R.All.PotableReusePOCPOU.bray.Includes_UV<-ANOSIM.All.PotableReusePOCPOU.bray.Includes_UV$statistic
ANOSIM.R.All.PotableReusePOCPOU.bray.Includes_UV

ADONIS2.R2.All.PotableReusePOCPOU.bray.Includes_UV<-ADONIS2.All.PotableReusePOCPOU.bray.Includes_UV$R2[1]
ADONIS2.R2.All.PotableReusePOCPOU.bray.Includes_UV

#Pull Together Pvalues to adjust
ANOSIM.GroupedP.All.PotableReusePOCPOU.bray["Includes_UV"]<-ANOSIM.Pval.All.PotableReusePOCPOU.bray.Includes_UV
ADONIS2.GroupedP.All.PotableReusePOCPOU.bray["Includes_UV"]<-ADONIS2.Pval.All.PotableReusePOCPOU.bray.Includes_UV

Betadis.GroupedP.All.PotableReusePOCPOU.bray["Includes_UV"]<-Betadis.Pval.All.PotableReusePOCPOU.bray.Includes_UV

ANOSIM.GroupedR.All.PotableReusePOCPOU.bray["Includes_UV"]<-ANOSIM.R.All.PotableReusePOCPOU.bray.Includes_UV
ADONIS2.GroupedR.All.PotableReusePOCPOU.bray["Includes_UV"]<-ADONIS2.R2.All.PotableReusePOCPOU.bray.Includes_UV

#Pairwise if needed, check if interesting or not. Remove if not.
count(All.PotableReusePOCPOU.meta$Includes_UV)

ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Includes_UV<-pairwise.adonis2(All.PotableReusePOCPOU.dist.bray ~ Includes_UV,
                                                                          data = All.PotableReusePOCPOU.meta,
                                                                          permutations = permutations)
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Includes_UV
#Includes_UV
##Includes_UV




# #################
# #################
# #################
# #################
# #Type_DS
# ANOSIM.All.PotableReusePOCPOU.bray.Type_DS <-anosim(All.PotableReusePOCPOU.dist.bray,
#                                                       All.PotableReusePOCPOU.meta$Type_DS ,
#                                                       permutations = permutations)
#
# ADONIS2.All.PotableReusePOCPOU.bray.Type_DS <-adonis2(All.PotableReusePOCPOU.dist.bray ~ Type_DS ,
#                                                         data = All.PotableReusePOCPOU.meta,
#                                                         permutations = permutations)
#
# Betadis.All.PotableReusePOCPOU.bray.Type_DS <-anova(betadisper(All.PotableReusePOCPOU.dist.bray,
#                                                                  All.PotableReusePOCPOU.meta$Type_DS ))
#
#
#
# ANOSIM.Pval.All.PotableReusePOCPOU.bray.Type_DS <-ANOSIM.All.PotableReusePOCPOU.bray.Type_DS $signif
# ANOSIM.Pval.All.PotableReusePOCPOU.bray.Type_DS
#
# ADONIS2.Pval.All.PotableReusePOCPOU.bray.Type_DS <-ADONIS2.All.PotableReusePOCPOU.bray.Type_DS $`Pr(>F)`[1]
# ADONIS2.Pval.All.PotableReusePOCPOU.bray.Type_DS
#
# Betadis.Pval.All.PotableReusePOCPOU.bray.Type_DS <-Betadis.All.PotableReusePOCPOU.bray.Type_DS $`Pr(>F)`[1]
# Betadis.Pval.All.PotableReusePOCPOU.bray.Type_DS
#
# ANOSIM.R.All.PotableReusePOCPOU.bray.Type_DS <-ANOSIM.All.PotableReusePOCPOU.bray.Type_DS $statistic
# ANOSIM.R.All.PotableReusePOCPOU.bray.Type_DS
#
# ADONIS2.R2.All.PotableReusePOCPOU.bray.Type_DS <-ADONIS2.All.PotableReusePOCPOU.bray.Type_DS $R2[1]
# ADONIS2.R2.All.PotableReusePOCPOU.bray.Type_DS
#
# #Pull Together Pvalues to adjust
# ANOSIM.GroupedP.All.PotableReusePOCPOU.bray["Type_DS "]<-ANOSIM.Pval.All.PotableReusePOCPOU.bray.Type_DS
# ADONIS2.GroupedP.All.PotableReusePOCPOU.bray["Type_DS "]<-ADONIS2.Pval.All.PotableReusePOCPOU.bray.Type_DS
#
# Betadis.GroupedP.All.PotableReusePOCPOU.bray["Type_DS "]<-Betadis.Pval.All.PotableReusePOCPOU.bray.Type_DS
#
# ANOSIM.GroupedR.All.PotableReusePOCPOU.bray["Type_DS "]<-ANOSIM.R.All.PotableReusePOCPOU.bray.Type_DS
# ADONIS2.GroupedR.All.PotableReusePOCPOU.bray["Type_DS "]<-ADONIS2.R2.All.PotableReusePOCPOU.bray.Type_DS
#
# #Pairwise if needed, check if interesting or not. Remove if not.
# count(All.PotableReusePOCPOU.meta$Type_DS )
#
# ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Type_DS <-pairwise.adonis2(All.PotableReusePOCPOU.dist.bray ~ Type_DS ,
#                                                                           data = All.PotableReusePOCPOU.meta,
#                                                                           permutations = permutations)
# ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Type_DS
# #Type_DS
# ##Type_DS


#################
#################
#################
#################
#Reuse.Reclaimed_Treatment
ANOSIM.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment <-anosim(All.PotableReusePOCPOU.dist.bray,
                                                    All.PotableReusePOCPOU.meta$Reuse.Reclaimed_Treatment ,
                                                    permutations = permutations)

ADONIS2.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment <-adonis2(All.PotableReusePOCPOU.dist.bray ~ Reuse.Reclaimed_Treatment ,
                                                      data = All.PotableReusePOCPOU.meta,
                                                      permutations = permutations)

Betadis.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment <-anova(betadisper(All.PotableReusePOCPOU.dist.bray,
                                                               All.PotableReusePOCPOU.meta$Reuse.Reclaimed_Treatment ))



ANOSIM.Pval.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment <-ANOSIM.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment $signif
ANOSIM.Pval.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment

ADONIS2.Pval.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment <-ADONIS2.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment $`Pr(>F)`[1]
ADONIS2.Pval.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment

Betadis.Pval.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment <-Betadis.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment $`Pr(>F)`[1]
Betadis.Pval.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment

ANOSIM.R.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment <-ANOSIM.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment $statistic
ANOSIM.R.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment

ADONIS2.R2.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment <-ADONIS2.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment $R2[1]
ADONIS2.R2.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment

#Pull Together Pvalues to adjust
ANOSIM.GroupedP.All.PotableReusePOCPOU.bray["Reuse.Reclaimed_Treatment "]<-ANOSIM.Pval.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment
ADONIS2.GroupedP.All.PotableReusePOCPOU.bray["Reuse.Reclaimed_Treatment "]<-ADONIS2.Pval.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment

Betadis.GroupedP.All.PotableReusePOCPOU.bray["Reuse.Reclaimed_Treatment "]<-Betadis.Pval.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment

ANOSIM.GroupedR.All.PotableReusePOCPOU.bray["Reuse.Reclaimed_Treatment "]<-ANOSIM.R.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment
ADONIS2.GroupedR.All.PotableReusePOCPOU.bray["Reuse.Reclaimed_Treatment "]<-ADONIS2.R2.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment

#Pairwise if needed, check if interesting or not. Remove if not.
count(All.PotableReusePOCPOU.meta$Reuse.Reclaimed_Treatment )

ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment <-pairwise.adonis2(All.PotableReusePOCPOU.dist.bray ~ Reuse.Reclaimed_Treatment ,
                                                                        data = All.PotableReusePOCPOU.meta,
                                                                        permutations = permutations)
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment
#Reuse.Reclaimed_Treatment
##Reuse.Reclaimed_Disinfection




#################
#################
#################
#################
#Reuse.Reclaimed_Treatment
ANOSIM.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment <-anosim(All.PotableReusePOCPOU.dist.bray,
                                                    All.PotableReusePOCPOU.meta$Reuse.Reclaimed_Treatment ,
                                                    permutations = permutations)

ADONIS2.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment <-adonis2(All.PotableReusePOCPOU.dist.bray ~ Reuse.Reclaimed_Treatment ,
                                                      data = All.PotableReusePOCPOU.meta,
                                                      permutations = permutations)

Betadis.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment <-anova(betadisper(All.PotableReusePOCPOU.dist.bray,
                                                               All.PotableReusePOCPOU.meta$Reuse.Reclaimed_Treatment ))



ANOSIM.Pval.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment <-ANOSIM.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment $signif
ANOSIM.Pval.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment

ADONIS2.Pval.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment <-ADONIS2.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment $`Pr(>F)`[1]
ADONIS2.Pval.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment

Betadis.Pval.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment <-Betadis.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment $`Pr(>F)`[1]
Betadis.Pval.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment

ANOSIM.R.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment <-ANOSIM.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment $statistic
ANOSIM.R.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment

ADONIS2.R2.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment <-ADONIS2.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment $R2[1]
ADONIS2.R2.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment

#Pull Together Pvalues to adjust
ANOSIM.GroupedP.All.PotableReusePOCPOU.bray["Reuse.Reclaimed_Treatment "]<-ANOSIM.Pval.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment
ADONIS2.GroupedP.All.PotableReusePOCPOU.bray["Reuse.Reclaimed_Treatment "]<-ADONIS2.Pval.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment

Betadis.GroupedP.All.PotableReusePOCPOU.bray["Reuse.Reclaimed_Treatment "]<-Betadis.Pval.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment

ANOSIM.GroupedR.All.PotableReusePOCPOU.bray["Reuse.Reclaimed_Treatment "]<-ANOSIM.R.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment
ADONIS2.GroupedR.All.PotableReusePOCPOU.bray["Reuse.Reclaimed_Treatment "]<-ADONIS2.R2.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment

#Pairwise if needed, check if interesting or not. Remove if not.
count(All.PotableReusePOCPOU.meta$Reuse.Reclaimed_Treatment )

ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment <-pairwise.adonis2(All.PotableReusePOCPOU.dist.bray ~ Reuse.Reclaimed_Treatment ,
                                                                        data = All.PotableReusePOCPOU.meta,
                                                                        permutations = permutations)
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment
#Reuse.Reclaimed_Treatment
##Reuse.Reclaimed_Treatment




#################
#################
#################
#################
#Reuse.Reclaimed_Blending
ANOSIM.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Blending<-anosim(All.PotableReusePOCPOU.dist.bray,
                                                    All.PotableReusePOCPOU.meta$Reuse.Reclaimed_Blending,
                                                    permutations = permutations)

ADONIS2.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Blending<-adonis2(All.PotableReusePOCPOU.dist.bray ~ Reuse.Reclaimed_Blending,
                                                      data = All.PotableReusePOCPOU.meta,
                                                      permutations = permutations)

Betadis.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Blending<-anova(betadisper(All.PotableReusePOCPOU.dist.bray,
                                                               All.PotableReusePOCPOU.meta$Reuse.Reclaimed_Blending))



ANOSIM.Pval.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Blending<-ANOSIM.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Blending$signif
ANOSIM.Pval.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Blending

ADONIS2.Pval.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Blending<-ADONIS2.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Blending$`Pr(>F)`[1]
ADONIS2.Pval.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Blending

Betadis.Pval.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Blending<-Betadis.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Blending$`Pr(>F)`[1]
Betadis.Pval.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Blending

ANOSIM.R.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Blending<-ANOSIM.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Blending$statistic
ANOSIM.R.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Blending

ADONIS2.R2.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Blending<-ADONIS2.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Blending$R2[1]
ADONIS2.R2.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Blending

#Pull Together Pvalues to adjust
ANOSIM.GroupedP.All.PotableReusePOCPOU.bray["Reuse.Reclaimed_Blending"]<-ANOSIM.Pval.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Blending
ADONIS2.GroupedP.All.PotableReusePOCPOU.bray["Reuse.Reclaimed_Blending"]<-ADONIS2.Pval.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Blending

Betadis.GroupedP.All.PotableReusePOCPOU.bray["Reuse.Reclaimed_Blending"]<-Betadis.Pval.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Blending

ANOSIM.GroupedR.All.PotableReusePOCPOU.bray["Reuse.Reclaimed_Blending"]<-ANOSIM.R.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Blending
ADONIS2.GroupedR.All.PotableReusePOCPOU.bray["Reuse.Reclaimed_Blending"]<-ADONIS2.R2.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Blending

#Pairwise if needed, check if interesting or not. Remove if not.
count(All.PotableReusePOCPOU.meta$Reuse.Reclaimed_Blending)

ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Blending<-pairwise.adonis2(All.PotableReusePOCPOU.dist.bray ~ Reuse.Reclaimed_Blending,
                                                                        data = All.PotableReusePOCPOU.meta,
                                                                        permutations = permutations)
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Blending
#Reuse.Reclaimed_Blending
##Reuse.Reclaimed_Blending




#################
#################
#################
#################
#PostBlend_Disinfection
ANOSIM.All.PotableReusePOCPOU.bray.PostBlend_Disinfection<-anosim(All.PotableReusePOCPOU.dist.bray,
                                                    All.PotableReusePOCPOU.meta$PostBlend_Disinfection,
                                                    permutations = permutations)

ADONIS2.All.PotableReusePOCPOU.bray.PostBlend_Disinfection<-adonis2(All.PotableReusePOCPOU.dist.bray ~ PostBlend_Disinfection,
                                                      data = All.PotableReusePOCPOU.meta,
                                                      permutations = permutations)

Betadis.All.PotableReusePOCPOU.bray.PostBlend_Disinfection<-anova(betadisper(All.PotableReusePOCPOU.dist.bray,
                                                               All.PotableReusePOCPOU.meta$PostBlend_Disinfection))



ANOSIM.Pval.All.PotableReusePOCPOU.bray.PostBlend_Disinfection<-ANOSIM.All.PotableReusePOCPOU.bray.PostBlend_Disinfection$signif
ANOSIM.Pval.All.PotableReusePOCPOU.bray.PostBlend_Disinfection

ADONIS2.Pval.All.PotableReusePOCPOU.bray.PostBlend_Disinfection<-ADONIS2.All.PotableReusePOCPOU.bray.PostBlend_Disinfection$`Pr(>F)`[1]
ADONIS2.Pval.All.PotableReusePOCPOU.bray.PostBlend_Disinfection

Betadis.Pval.All.PotableReusePOCPOU.bray.PostBlend_Disinfection<-Betadis.All.PotableReusePOCPOU.bray.PostBlend_Disinfection$`Pr(>F)`[1]
Betadis.Pval.All.PotableReusePOCPOU.bray.PostBlend_Disinfection

ANOSIM.R.All.PotableReusePOCPOU.bray.PostBlend_Disinfection<-ANOSIM.All.PotableReusePOCPOU.bray.PostBlend_Disinfection$statistic
ANOSIM.R.All.PotableReusePOCPOU.bray.PostBlend_Disinfection

ADONIS2.R2.All.PotableReusePOCPOU.bray.PostBlend_Disinfection<-ADONIS2.All.PotableReusePOCPOU.bray.PostBlend_Disinfection$R2[1]
ADONIS2.R2.All.PotableReusePOCPOU.bray.PostBlend_Disinfection

#Pull Together Pvalues to adjust
ANOSIM.GroupedP.All.PotableReusePOCPOU.bray["PostBlend_Disinfection"]<-ANOSIM.Pval.All.PotableReusePOCPOU.bray.PostBlend_Disinfection
ADONIS2.GroupedP.All.PotableReusePOCPOU.bray["PostBlend_Disinfection"]<-ADONIS2.Pval.All.PotableReusePOCPOU.bray.PostBlend_Disinfection

Betadis.GroupedP.All.PotableReusePOCPOU.bray["PostBlend_Disinfection"]<-Betadis.Pval.All.PotableReusePOCPOU.bray.PostBlend_Disinfection

ANOSIM.GroupedR.All.PotableReusePOCPOU.bray["PostBlend_Disinfection"]<-ANOSIM.R.All.PotableReusePOCPOU.bray.PostBlend_Disinfection
ADONIS2.GroupedR.All.PotableReusePOCPOU.bray["PostBlend_Disinfection"]<-ADONIS2.R2.All.PotableReusePOCPOU.bray.PostBlend_Disinfection

#Pairwise if needed, check if interesting or not. Remove if not.
count(All.PotableReusePOCPOU.meta$PostBlend_Disinfection)

ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.PostBlend_Disinfection<-pairwise.adonis2(All.PotableReusePOCPOU.dist.bray ~ PostBlend_Disinfection,
                                                                        data = All.PotableReusePOCPOU.meta,
                                                                        permutations = permutations)
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.PostBlend_Disinfection
#PostBlend_Disinfection
##PostBlend_Disinfection




#################
#################
#################
#################
#PostBlend_Treatment
ANOSIM.All.PotableReusePOCPOU.bray.PostBlend_Treatment<-anosim(All.PotableReusePOCPOU.dist.bray,
                                                    All.PotableReusePOCPOU.meta$PostBlend_Treatment,
                                                    permutations = permutations)

ADONIS2.All.PotableReusePOCPOU.bray.PostBlend_Treatment<-adonis2(All.PotableReusePOCPOU.dist.bray ~ PostBlend_Treatment,
                                                      data = All.PotableReusePOCPOU.meta,
                                                      permutations = permutations)

Betadis.All.PotableReusePOCPOU.bray.PostBlend_Treatment<-anova(betadisper(All.PotableReusePOCPOU.dist.bray,
                                                               All.PotableReusePOCPOU.meta$PostBlend_Treatment))



ANOSIM.Pval.All.PotableReusePOCPOU.bray.PostBlend_Treatment<-ANOSIM.All.PotableReusePOCPOU.bray.PostBlend_Treatment$signif
ANOSIM.Pval.All.PotableReusePOCPOU.bray.PostBlend_Treatment

ADONIS2.Pval.All.PotableReusePOCPOU.bray.PostBlend_Treatment<-ADONIS2.All.PotableReusePOCPOU.bray.PostBlend_Treatment$`Pr(>F)`[1]
ADONIS2.Pval.All.PotableReusePOCPOU.bray.PostBlend_Treatment

Betadis.Pval.All.PotableReusePOCPOU.bray.PostBlend_Treatment<-Betadis.All.PotableReusePOCPOU.bray.PostBlend_Treatment$`Pr(>F)`[1]
Betadis.Pval.All.PotableReusePOCPOU.bray.PostBlend_Treatment

ANOSIM.R.All.PotableReusePOCPOU.bray.PostBlend_Treatment<-ANOSIM.All.PotableReusePOCPOU.bray.PostBlend_Treatment$statistic
ANOSIM.R.All.PotableReusePOCPOU.bray.PostBlend_Treatment

ADONIS2.R2.All.PotableReusePOCPOU.bray.PostBlend_Treatment<-ADONIS2.All.PotableReusePOCPOU.bray.PostBlend_Treatment$R2[1]
ADONIS2.R2.All.PotableReusePOCPOU.bray.PostBlend_Treatment

#Pull Together Pvalues to adjust
ANOSIM.GroupedP.All.PotableReusePOCPOU.bray["PostBlend_Treatment"]<-ANOSIM.Pval.All.PotableReusePOCPOU.bray.PostBlend_Treatment
ADONIS2.GroupedP.All.PotableReusePOCPOU.bray["PostBlend_Treatment"]<-ADONIS2.Pval.All.PotableReusePOCPOU.bray.PostBlend_Treatment

Betadis.GroupedP.All.PotableReusePOCPOU.bray["PostBlend_Treatment"]<-Betadis.Pval.All.PotableReusePOCPOU.bray.PostBlend_Treatment

ANOSIM.GroupedR.All.PotableReusePOCPOU.bray["PostBlend_Treatment"]<-ANOSIM.R.All.PotableReusePOCPOU.bray.PostBlend_Treatment
ADONIS2.GroupedR.All.PotableReusePOCPOU.bray["PostBlend_Treatment"]<-ADONIS2.R2.All.PotableReusePOCPOU.bray.PostBlend_Treatment

#Pairwise if needed, check if interesting or not. Remove if not.
count(All.PotableReusePOCPOU.meta$PostBlend_Treatment)

ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.PostBlend_Treatment<-pairwise.adonis2(All.PotableReusePOCPOU.dist.bray ~ PostBlend_Treatment,
                                                                        data = All.PotableReusePOCPOU.meta,
                                                                        permutations = permutations)
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.PostBlend_Treatment
#PostBlend_Treatment
##PostBlend_Treatment





#Adjust PValue, Check Significance
ANOSIM.GroupedP.Adjusted.All.PotableReusePOCPOU.bray<-p.adjust(ANOSIM.GroupedP.All.PotableReusePOCPOU.bray,method="BY")
p.adjust(ANOSIM.GroupedP.Adjusted.All.PotableReusePOCPOU.bray,method="BY")<.05

ADONIS2.GroupedP.Adjusted.All.PotableReusePOCPOU.bray<-p.adjust(ADONIS2.GroupedP.All.PotableReusePOCPOU.bray,method="BY")
p.adjust(ADONIS2.GroupedP.Adjusted.All.PotableReusePOCPOU.bray,method="BY")<.05

Betadis.GroupedP.All.PotableReusePOCPOU.bray<0.05
Betadis.GroupedP.Adjusted.All.PotableReusePOCPOU.bray<-p.adjust(Betadis.GroupedP.All.PotableReusePOCPOU.bray,method="BY")
p.adjust(Betadis.GroupedP.Adjusted.All.PotableReusePOCPOU.bray,method="BY")<.05


#Check Value of Pairwise Sub comparisons
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Class1_Pot_NonPot
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Classification_1
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Classification_1_mod
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Classification_2
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Climate
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Climate_Region
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Final_Disinfection_Residual
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Primers
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Region
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Sample_Local
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Trimmed
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Extraction_Type
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Potable_Source_Water
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Potable_Disinfection_Primary
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Potable_Treatment_Type
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Disinfection_Total
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Treatment_Total
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Treatment_Classification_Number
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Treatment_Classification_Cat
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Mixed_Source_Water
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Potable_Overlap
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Includes_Ozone
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Includes_RO_UF
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Includes_UV
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Type_DS
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Disinfection
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Treatment
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.Reuse.Reclaimed_Blending
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.PostBlend_Disinfection
ADONIS2.Pairwise.All.PotableReusePOCPOU.bray.PostBlend_Treatment

#Permanova on Multiple Variables that are significant and not heterogeneous (not sig from Beta, and sig from anosim)
ADONIS2.All.PotableReusePOCPOU.bray.multiple_comparisons<-adonis2(All.PotableReusePOCPOU.dist.bray ~ 
                                                                  Treatment_Total+ 
                                                                  Disinfection_Total+
                                                                  Final_Disinfection_Residual+
                                                                  Sample_Local+
                                                                  Climate,  
                                                 data = All.PotableReusePOCPOU.meta, 
                                                 permutations = permutations,
                                                 by="terms")
ADONIS2.All.PotableReusePOCPOU.bray.multiple_comparisons


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
All.PotableReusePOCPOU.Stats<-data.frame(
  "ANOSIM.P.All.PotableReusePOCPOU.bray"= ANOSIM.GroupedP.All.PotableReusePOCPOU.bray,
  "ANOSIM.P.Adj.All.PotableReusePOCPOU.bray"= ANOSIM.GroupedP.Adjusted.All.PotableReusePOCPOU.bray,
  "ADONIS2.P.All.PotableReusePOCPOU.bray"= ADONIS2.GroupedP.All.PotableReusePOCPOU.bray,
  "ADONIS2.P.Adj.All.PotableReusePOCPOU.bray"= ADONIS2.GroupedP.Adjusted.All.PotableReusePOCPOU.bray,
  "Betadis.P.All.PotableReusePOCPOU.bray"= Betadis.GroupedP.All.PotableReusePOCPOU.bray,
  "Betadis.P.Adj.All.PotableReusePOCPOU.bray"= Betadis.GroupedP.Adjusted.All.PotableReusePOCPOU.bray,
  "ANOSIM.R.All.PotableReusePOCPOU.bray"=  ANOSIM.GroupedR.All.PotableReusePOCPOU.bray,
  "ADONIS2.R.All.PotableReusePOCPOU.bray"= ADONIS2.GroupedR.All.PotableReusePOCPOU.bray)


write.xlsx(All.PotableReusePOCPOU.Stats, file = "C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Stats/Beta/All.PotableReusePOCPOU.List.xlsx", rowNames = TRUE)

