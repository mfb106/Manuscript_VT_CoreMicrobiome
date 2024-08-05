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
Metadata_PotableBothPOCPOU<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_DS_Special_V2', rowNames = TRUE)
Metadata_POC<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_Final', rowNames = TRUE)
Metadata_POU<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'DS_Final', rowNames = TRUE)

#Order Meta Data 
Metadata_PotableBothPOCPOU$Classification_1 = factor(Metadata_PotableBothPOCPOU$Classification_1, levels=c("Potable Conventional","Potable Reuse","Non-potable Reuse","Blank"))
Metadata_POC$Classification_1 = factor(Metadata_POC$Classification_1, levels=c("Potable Conventional","Potable Reuse","Non-potable Reuse","Blank"))
Metadata_POU$Classification_1 = factor(Metadata_POU$Classification_1, levels=c("Potable Conventional","Potable Reuse","Non-potable Reuse","Blank"))

######All Samples, PotableBothPOCPOU
Meta.All.PotableBothPOCPOU<-Metadata_PotableBothPOCPOU
Meta.All.PotableBothPOCPOU<-Meta.All.PotableBothPOCPOU[Meta.All.PotableBothPOCPOU$Description != 'ConnerPrimer58',] #removed because of low read count
Meta.All.PotableBothPOCPOU<-Meta.All.PotableBothPOCPOU[Meta.All.PotableBothPOCPOU$Description != 'PR1WEEF',] #removed because of high read count 

Meta.All.PotableBothPOCPOU.Filter<-Meta.All.PotableBothPOCPOU[Meta.All.PotableBothPOCPOU$Classification_1 == 'Potable Conventional'| Meta.All.PotableBothPOCPOU$Classification_1 == 'Potable Reuse',]#can add Blank
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

All.PotableBothPOCPOU.dist.bray<- phyloseq::distance(physeq.combined.All.PotableBothPOCPOU, method = "bray")
All.PotableBothPOCPOU.meta<- data.frame(sample_data(physeq.combined.All.PotableBothPOCPOU))

ANOSIM.GroupedP.All.PotableBothPOCPOU.bray<-numeric() 
ADONIS2.GroupedP.All.PotableBothPOCPOU.bray<-numeric()
Betadis.GroupedP.All.PotableBothPOCPOU.bray<-numeric()
ANOSIM.GroupedR.All.PotableBothPOCPOU.bray<-numeric() 
ADONIS2.GroupedR.All.PotableBothPOCPOU.bray<-numeric()


# #################
# #################
# #################
# #################
# #Class1_Pot_NonPot 
# ANOSIM.All.PotableBothPOCPOU.bray.Class1_Pot_NonPot<-anosim(All.PotableBothPOCPOU.dist.bray,
#                                                              All.PotableBothPOCPOU.meta$Class1_Pot_NonPot, 
#                                                              permutations = permutations)
# 
# ADONIS2.All.PotableBothPOCPOU.bray.Class1_Pot_NonPot<-adonis2(All.PotableBothPOCPOU.dist.bray ~ Class1_Pot_NonPot, 
#                                                                data = All.PotableBothPOCPOU.meta, 
#                                                                permutations = permutations)
# 
# Betadis.All.PotableBothPOCPOU.bray.Class1_Pot_NonPot<-anova(betadisper(All.PotableBothPOCPOU.dist.bray,
#                                                                         All.PotableBothPOCPOU.meta$Class1_Pot_NonPot))
# 
# 
# 
# ANOSIM.Pval.All.PotableBothPOCPOU.bray.Class1_Pot_NonPot<-ANOSIM.All.PotableBothPOCPOU.bray.Class1_Pot_NonPot$signif
# ANOSIM.Pval.All.PotableBothPOCPOU.bray.Class1_Pot_NonPot
# 
# ADONIS2.Pval.All.PotableBothPOCPOU.bray.Class1_Pot_NonPot<-ADONIS2.All.PotableBothPOCPOU.bray.Class1_Pot_NonPot$`Pr(>F)`[1]
# ADONIS2.Pval.All.PotableBothPOCPOU.bray.Class1_Pot_NonPot  
# 
# Betadis.Pval.All.PotableBothPOCPOU.bray.Class1_Pot_NonPot<-Betadis.All.PotableBothPOCPOU.bray.Class1_Pot_NonPot$`Pr(>F)`[1]
# Betadis.Pval.All.PotableBothPOCPOU.bray.Class1_Pot_NonPot
# 
# ANOSIM.R.All.PotableBothPOCPOU.bray.Class1_Pot_NonPot<-ANOSIM.All.PotableBothPOCPOU.bray.Class1_Pot_NonPot$statistic
# ANOSIM.R.All.PotableBothPOCPOU.bray.Class1_Pot_NonPot
# 
# ADONIS2.R2.All.PotableBothPOCPOU.bray.Class1_Pot_NonPot<-ADONIS2.All.PotableBothPOCPOU.bray.Class1_Pot_NonPot$R2[1]
# ADONIS2.R2.All.PotableBothPOCPOU.bray.Class1_Pot_NonPot  
# 
# #Pull Together Pvalues to adjust 
# ANOSIM.GroupedP.All.PotableBothPOCPOU.bray["Class1_Pot_NonPot"]<-ANOSIM.Pval.All.PotableBothPOCPOU.bray.Class1_Pot_NonPot
# ADONIS2.GroupedP.All.PotableBothPOCPOU.bray["Class1_Pot_NonPot"]<-ADONIS2.Pval.All.PotableBothPOCPOU.bray.Class1_Pot_NonPot
# 
# Betadis.GroupedP.All.PotableBothPOCPOU.bray["Class1_Pot_NonPot"]<-Betadis.Pval.All.PotableBothPOCPOU.bray.Class1_Pot_NonPot
# 
# ANOSIM.GroupedR.All.PotableBothPOCPOU.bray["Class1_Pot_NonPot"]<-ANOSIM.R.All.PotableBothPOCPOU.bray.Class1_Pot_NonPot
# ADONIS2.GroupedR.All.PotableBothPOCPOU.bray["Class1_Pot_NonPot"]<-ADONIS2.R2.All.PotableBothPOCPOU.bray.Class1_Pot_NonPot
# 
# #Pairwise if needed, check if interesting or not. Remove if not. 
# count(All.PotableBothPOCPOU.meta$Class1_Pot_NonPot)
# 
# ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Class1_Pot_NonPot<-pairwise.adonis2(All.PotableBothPOCPOU.dist.bray ~ Class1_Pot_NonPot, 
#                                                                                  data = All.PotableBothPOCPOU.meta, 
#                                                                                  permutations = permutations)
# ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Class1_Pot_NonPot




#################
#################
#################
#################
#Classification_1 
ANOSIM.All.PotableBothPOCPOU.bray.Classification_1<-anosim(All.PotableBothPOCPOU.dist.bray,
                                                            All.PotableBothPOCPOU.meta$Classification_1, 
                                                            permutations = permutations)

ADONIS2.All.PotableBothPOCPOU.bray.Classification_1<-adonis2(All.PotableBothPOCPOU.dist.bray ~ Classification_1, 
                                                              data = All.PotableBothPOCPOU.meta, 
                                                              permutations = permutations)

Betadis.All.PotableBothPOCPOU.bray.Classification_1<-anova(betadisper(All.PotableBothPOCPOU.dist.bray,
                                                                       All.PotableBothPOCPOU.meta$Classification_1))



ANOSIM.Pval.All.PotableBothPOCPOU.bray.Classification_1<-ANOSIM.All.PotableBothPOCPOU.bray.Classification_1$signif
ANOSIM.Pval.All.PotableBothPOCPOU.bray.Classification_1

ADONIS2.Pval.All.PotableBothPOCPOU.bray.Classification_1<-ADONIS2.All.PotableBothPOCPOU.bray.Classification_1$`Pr(>F)`[1]
ADONIS2.Pval.All.PotableBothPOCPOU.bray.Classification_1  

Betadis.Pval.All.PotableBothPOCPOU.bray.Classification_1<-Betadis.All.PotableBothPOCPOU.bray.Classification_1$`Pr(>F)`[1]
Betadis.Pval.All.PotableBothPOCPOU.bray.Classification_1

ANOSIM.R.All.PotableBothPOCPOU.bray.Classification_1<-ANOSIM.All.PotableBothPOCPOU.bray.Classification_1$statistic
ANOSIM.R.All.PotableBothPOCPOU.bray.Classification_1

ADONIS2.R2.All.PotableBothPOCPOU.bray.Classification_1<-ADONIS2.All.PotableBothPOCPOU.bray.Classification_1$R2[1]
ADONIS2.R2.All.PotableBothPOCPOU.bray.Classification_1  

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.PotableBothPOCPOU.bray["Classification_1"]<-ANOSIM.Pval.All.PotableBothPOCPOU.bray.Classification_1
ADONIS2.GroupedP.All.PotableBothPOCPOU.bray["Classification_1"]<-ADONIS2.Pval.All.PotableBothPOCPOU.bray.Classification_1

Betadis.GroupedP.All.PotableBothPOCPOU.bray["Classification_1"]<-Betadis.Pval.All.PotableBothPOCPOU.bray.Classification_1

ANOSIM.GroupedR.All.PotableBothPOCPOU.bray["Classification_1"]<-ANOSIM.R.All.PotableBothPOCPOU.bray.Classification_1
ADONIS2.GroupedR.All.PotableBothPOCPOU.bray["Classification_1"]<-ADONIS2.R2.All.PotableBothPOCPOU.bray.Classification_1

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.PotableBothPOCPOU.meta$Classification_1)

ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Classification_1<-pairwise.adonis2(All.PotableBothPOCPOU.dist.bray ~ Classification_1, 
                                                                                data = All.PotableBothPOCPOU.meta, 
                                                                                permutations = permutations)
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Classification_1





# #################
# #################
# #################
# #################
# #Classification_1_mod 
# ANOSIM.All.PotableBothPOCPOU.bray.Classification_1_mod<-anosim(All.PotableBothPOCPOU.dist.bray,
#                                                                 All.PotableBothPOCPOU.meta$Classification_1_mod, 
#                                                                 permutations = permutations)
# 
# ADONIS2.All.PotableBothPOCPOU.bray.Classification_1_mod<-adonis2(All.PotableBothPOCPOU.dist.bray ~ Classification_1_mod, 
#                                                                   data = All.PotableBothPOCPOU.meta, 
#                                                                   permutations = permutations)
# 
# Betadis.All.PotableBothPOCPOU.bray.Classification_1_mod<-anova(betadisper(All.PotableBothPOCPOU.dist.bray,
#                                                                            All.PotableBothPOCPOU.meta$Classification_1_mod))
# 
# 
# 
# ANOSIM.Pval.All.PotableBothPOCPOU.bray.Classification_1_mod<-ANOSIM.All.PotableBothPOCPOU.bray.Classification_1_mod$signif
# ANOSIM.Pval.All.PotableBothPOCPOU.bray.Classification_1_mod
# 
# ADONIS2.Pval.All.PotableBothPOCPOU.bray.Classification_1_mod<-ADONIS2.All.PotableBothPOCPOU.bray.Classification_1_mod$`Pr(>F)`[1]
# ADONIS2.Pval.All.PotableBothPOCPOU.bray.Classification_1_mod  
# 
# Betadis.Pval.All.PotableBothPOCPOU.bray.Classification_1_mod<-Betadis.All.PotableBothPOCPOU.bray.Classification_1_mod$`Pr(>F)`[1]
# Betadis.Pval.All.PotableBothPOCPOU.bray.Classification_1_mod
# 
# ANOSIM.R.All.PotableBothPOCPOU.bray.Classification_1_mod<-ANOSIM.All.PotableBothPOCPOU.bray.Classification_1_mod$statistic
# ANOSIM.R.All.PotableBothPOCPOU.bray.Classification_1_mod
# 
# ADONIS2.R2.All.PotableBothPOCPOU.bray.Classification_1_mod<-ADONIS2.All.PotableBothPOCPOU.bray.Classification_1_mod$R2[1]
# ADONIS2.R2.All.PotableBothPOCPOU.bray.Classification_1_mod  
# 
# #Pull Together Pvalues to adjust 
# ANOSIM.GroupedP.All.PotableBothPOCPOU.bray["Classification_1_mod"]<-ANOSIM.Pval.All.PotableBothPOCPOU.bray.Classification_1_mod
# ADONIS2.GroupedP.All.PotableBothPOCPOU.bray["Classification_1_mod"]<-ADONIS2.Pval.All.PotableBothPOCPOU.bray.Classification_1_mod
# 
# Betadis.GroupedP.All.PotableBothPOCPOU.bray["Classification_1_mod"]<-Betadis.Pval.All.PotableBothPOCPOU.bray.Classification_1_mod
# 
# ANOSIM.GroupedR.All.PotableBothPOCPOU.bray["Classification_1_mod"]<-ANOSIM.R.All.PotableBothPOCPOU.bray.Classification_1_mod
# ADONIS2.GroupedR.All.PotableBothPOCPOU.bray["Classification_1_mod"]<-ADONIS2.R2.All.PotableBothPOCPOU.bray.Classification_1_mod
# 
# #Pairwise if needed, check if interesting or not. Remove if not. 
# count(All.PotableBothPOCPOU.meta$Classification_1_mod)
# 
# ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Classification_1_mod<-pairwise.adonis2(All.PotableBothPOCPOU.dist.bray ~ Classification_1_mod, 
#                                                                                     data = All.PotableBothPOCPOU.meta, 
#                                                                                     permutations = permutations)
# ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Classification_1_mod
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
# ANOSIM.All.PotableBothPOCPOU.bray.Classification_2<-anosim(All.PotableBothPOCPOU.dist.bray,
#                                                             All.PotableBothPOCPOU.meta$Classification_2, 
#                                                             permutations = permutations)
# 
# ADONIS2.All.PotableBothPOCPOU.bray.Classification_2<-adonis2(All.PotableBothPOCPOU.dist.bray ~ Classification_2, 
#                                                               data = All.PotableBothPOCPOU.meta, 
#                                                               permutations = permutations)
# 
# Betadis.All.PotableBothPOCPOU.bray.Classification_2<-anova(betadisper(All.PotableBothPOCPOU.dist.bray,
#                                                                        All.PotableBothPOCPOU.meta$Classification_2))
# 
# 
# 
# ANOSIM.Pval.All.PotableBothPOCPOU.bray.Classification_2<-ANOSIM.All.PotableBothPOCPOU.bray.Classification_2$signif
# ANOSIM.Pval.All.PotableBothPOCPOU.bray.Classification_2
# 
# ADONIS2.Pval.All.PotableBothPOCPOU.bray.Classification_2<-ADONIS2.All.PotableBothPOCPOU.bray.Classification_2$`Pr(>F)`[1]
# ADONIS2.Pval.All.PotableBothPOCPOU.bray.Classification_2  
# 
# Betadis.Pval.All.PotableBothPOCPOU.bray.Classification_2<-Betadis.All.PotableBothPOCPOU.bray.Classification_2$`Pr(>F)`[1]
# Betadis.Pval.All.PotableBothPOCPOU.bray.Classification_2
# 
# ANOSIM.R.All.PotableBothPOCPOU.bray.Classification_2<-ANOSIM.All.PotableBothPOCPOU.bray.Classification_2$statistic
# ANOSIM.R.All.PotableBothPOCPOU.bray.Classification_2
# 
# ADONIS2.R2.All.PotableBothPOCPOU.bray.Classification_2<-ADONIS2.All.PotableBothPOCPOU.bray.Classification_2$R2[1]
# ADONIS2.R2.All.PotableBothPOCPOU.bray.Classification_2  
# 
# #Pull Together Pvalues to adjust 
# ANOSIM.GroupedP.All.PotableBothPOCPOU.bray["Classification_2"]<-ANOSIM.Pval.All.PotableBothPOCPOU.bray.Classification_2
# ADONIS2.GroupedP.All.PotableBothPOCPOU.bray["Classification_2"]<-ADONIS2.Pval.All.PotableBothPOCPOU.bray.Classification_2
# 
# Betadis.GroupedP.All.PotableBothPOCPOU.bray["Classification_2"]<-Betadis.Pval.All.PotableBothPOCPOU.bray.Classification_2
# 
# ANOSIM.GroupedR.All.PotableBothPOCPOU.bray["Classification_2"]<-ANOSIM.R.All.PotableBothPOCPOU.bray.Classification_2
# ADONIS2.GroupedR.All.PotableBothPOCPOU.bray["Classification_2"]<-ADONIS2.R2.All.PotableBothPOCPOU.bray.Classification_2
# 
# #Pairwise if needed, check if interesting or not. Remove if not. 
# count(All.PotableBothPOCPOU.meta$Classification_2)
# 
# ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Classification_2<-pairwise.adonis2(All.PotableBothPOCPOU.dist.bray ~ Classification_2, 
#                                                                                 data = All.PotableBothPOCPOU.meta, 
#                                                                                 permutations = permutations)
# ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Classification_2
# #Classification_2
# 
# 
# 
# 
# 
# #################
# #################
# #################
# #################
# #Climate_Region 
# ANOSIM.All.PotableBothPOCPOU.bray.Climate_Region<-anosim(All.PotableBothPOCPOU.dist.bray,
#                                                           All.PotableBothPOCPOU.meta$Climate_Region, 
#                                                           permutations = permutations)
# 
# ADONIS2.All.PotableBothPOCPOU.bray.Climate_Region<-adonis2(All.PotableBothPOCPOU.dist.bray ~ Climate_Region, 
#                                                             data = All.PotableBothPOCPOU.meta, 
#                                                             permutations = permutations)
# 
# Betadis.All.PotableBothPOCPOU.bray.Climate_Region<-anova(betadisper(All.PotableBothPOCPOU.dist.bray,
#                                                                      All.PotableBothPOCPOU.meta$Climate_Region))
# 
# 
# 
# ANOSIM.Pval.All.PotableBothPOCPOU.bray.Climate_Region<-ANOSIM.All.PotableBothPOCPOU.bray.Climate_Region$signif
# ANOSIM.Pval.All.PotableBothPOCPOU.bray.Climate_Region
# 
# ADONIS2.Pval.All.PotableBothPOCPOU.bray.Climate_Region<-ADONIS2.All.PotableBothPOCPOU.bray.Climate_Region$`Pr(>F)`[1]
# ADONIS2.Pval.All.PotableBothPOCPOU.bray.Climate_Region  
# 
# Betadis.Pval.All.PotableBothPOCPOU.bray.Climate_Region<-Betadis.All.PotableBothPOCPOU.bray.Climate_Region$`Pr(>F)`[1]
# Betadis.Pval.All.PotableBothPOCPOU.bray.Climate_Region
# 
# ANOSIM.R.All.PotableBothPOCPOU.bray.Climate_Region<-ANOSIM.All.PotableBothPOCPOU.bray.Climate_Region$statistic
# ANOSIM.R.All.PotableBothPOCPOU.bray.Climate_Region
# 
# ADONIS2.R2.All.PotableBothPOCPOU.bray.Climate_Region<-ADONIS2.All.PotableBothPOCPOU.bray.Climate_Region$R2[1]
# ADONIS2.R2.All.PotableBothPOCPOU.bray.Climate_Region  
# 
# #Pull Together Pvalues to adjust 
# ANOSIM.GroupedP.All.PotableBothPOCPOU.bray["Climate_Region"]<-ANOSIM.Pval.All.PotableBothPOCPOU.bray.Climate_Region
# ADONIS2.GroupedP.All.PotableBothPOCPOU.bray["Climate_Region"]<-ADONIS2.Pval.All.PotableBothPOCPOU.bray.Climate_Region
# 
# Betadis.GroupedP.All.PotableBothPOCPOU.bray["Climate_Region"]<-Betadis.Pval.All.PotableBothPOCPOU.bray.Climate_Region
# 
# ANOSIM.GroupedR.All.PotableBothPOCPOU.bray["Climate_Region"]<-ANOSIM.R.All.PotableBothPOCPOU.bray.Climate_Region
# ADONIS2.GroupedR.All.PotableBothPOCPOU.bray["Climate_Region"]<-ADONIS2.R2.All.PotableBothPOCPOU.bray.Climate_Region
# 
# #Pairwise if needed, check if interesting or not. Remove if not. 
# count(All.PotableBothPOCPOU.meta$Climate_Region)
# 
# ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Climate_Region<-pairwise.adonis2(All.PotableBothPOCPOU.dist.bray ~ Climate_Region, 
#                                                                               data = All.PotableBothPOCPOU.meta, 
#                                                                               permutations = permutations)
# ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Climate_Region
# #Climate_Region





#################
#################
#################
#################
#Climate 
ANOSIM.All.PotableBothPOCPOU.bray.Climate<-anosim(All.PotableBothPOCPOU.dist.bray,
                                                   All.PotableBothPOCPOU.meta$Climate, 
                                                   permutations = permutations)

ADONIS2.All.PotableBothPOCPOU.bray.Climate<-adonis2(All.PotableBothPOCPOU.dist.bray ~ Climate, 
                                                     data = All.PotableBothPOCPOU.meta, 
                                                     permutations = permutations)

Betadis.All.PotableBothPOCPOU.bray.Climate<-anova(betadisper(All.PotableBothPOCPOU.dist.bray,
                                                              All.PotableBothPOCPOU.meta$Climate))



ANOSIM.Pval.All.PotableBothPOCPOU.bray.Climate<-ANOSIM.All.PotableBothPOCPOU.bray.Climate$signif
ANOSIM.Pval.All.PotableBothPOCPOU.bray.Climate

ADONIS2.Pval.All.PotableBothPOCPOU.bray.Climate<-ADONIS2.All.PotableBothPOCPOU.bray.Climate$`Pr(>F)`[1]
ADONIS2.Pval.All.PotableBothPOCPOU.bray.Climate  

Betadis.Pval.All.PotableBothPOCPOU.bray.Climate<-Betadis.All.PotableBothPOCPOU.bray.Climate$`Pr(>F)`[1]
Betadis.Pval.All.PotableBothPOCPOU.bray.Climate

ANOSIM.R.All.PotableBothPOCPOU.bray.Climate<-ANOSIM.All.PotableBothPOCPOU.bray.Climate$statistic
ANOSIM.R.All.PotableBothPOCPOU.bray.Climate

ADONIS2.R2.All.PotableBothPOCPOU.bray.Climate<-ADONIS2.All.PotableBothPOCPOU.bray.Climate$R2[1]
ADONIS2.R2.All.PotableBothPOCPOU.bray.Climate  

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.PotableBothPOCPOU.bray["Climate"]<-ANOSIM.Pval.All.PotableBothPOCPOU.bray.Climate
ADONIS2.GroupedP.All.PotableBothPOCPOU.bray["Climate"]<-ADONIS2.Pval.All.PotableBothPOCPOU.bray.Climate

Betadis.GroupedP.All.PotableBothPOCPOU.bray["Climate"]<-Betadis.Pval.All.PotableBothPOCPOU.bray.Climate

ANOSIM.GroupedR.All.PotableBothPOCPOU.bray["Climate"]<-ANOSIM.R.All.PotableBothPOCPOU.bray.Climate
ADONIS2.GroupedR.All.PotableBothPOCPOU.bray["Climate"]<-ADONIS2.R2.All.PotableBothPOCPOU.bray.Climate

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.PotableBothPOCPOU.meta$Climate)

ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Climate<-pairwise.adonis2(All.PotableBothPOCPOU.dist.bray ~ Climate, 
                                                                       data = All.PotableBothPOCPOU.meta, 
                                                                       permutations = permutations)
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Climate
#Climate





#################
#################
#################
#################
#Region 
ANOSIM.All.PotableBothPOCPOU.bray.Region<-anosim(All.PotableBothPOCPOU.dist.bray,
                                                  All.PotableBothPOCPOU.meta$Region, 
                                                  permutations = permutations)

ADONIS2.All.PotableBothPOCPOU.bray.Region<-adonis2(All.PotableBothPOCPOU.dist.bray ~ Region, 
                                                    data = All.PotableBothPOCPOU.meta, 
                                                    permutations = permutations)

Betadis.All.PotableBothPOCPOU.bray.Region<-anova(betadisper(All.PotableBothPOCPOU.dist.bray,
                                                             All.PotableBothPOCPOU.meta$Region))



ANOSIM.Pval.All.PotableBothPOCPOU.bray.Region<-ANOSIM.All.PotableBothPOCPOU.bray.Region$signif
ANOSIM.Pval.All.PotableBothPOCPOU.bray.Region

ADONIS2.Pval.All.PotableBothPOCPOU.bray.Region<-ADONIS2.All.PotableBothPOCPOU.bray.Region$`Pr(>F)`[1]
ADONIS2.Pval.All.PotableBothPOCPOU.bray.Region  

Betadis.Pval.All.PotableBothPOCPOU.bray.Region<-Betadis.All.PotableBothPOCPOU.bray.Region$`Pr(>F)`[1]
Betadis.Pval.All.PotableBothPOCPOU.bray.Region

ANOSIM.R.All.PotableBothPOCPOU.bray.Region<-ANOSIM.All.PotableBothPOCPOU.bray.Region$statistic
ANOSIM.R.All.PotableBothPOCPOU.bray.Region

ADONIS2.R2.All.PotableBothPOCPOU.bray.Region<-ADONIS2.All.PotableBothPOCPOU.bray.Region$R2[1]
ADONIS2.R2.All.PotableBothPOCPOU.bray.Region  

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.PotableBothPOCPOU.bray["Region"]<-ANOSIM.Pval.All.PotableBothPOCPOU.bray.Region
ADONIS2.GroupedP.All.PotableBothPOCPOU.bray["Region"]<-ADONIS2.Pval.All.PotableBothPOCPOU.bray.Region

Betadis.GroupedP.All.PotableBothPOCPOU.bray["Region"]<-Betadis.Pval.All.PotableBothPOCPOU.bray.Region

ANOSIM.GroupedR.All.PotableBothPOCPOU.bray["Region"]<-ANOSIM.R.All.PotableBothPOCPOU.bray.Region
ADONIS2.GroupedR.All.PotableBothPOCPOU.bray["Region"]<-ADONIS2.R2.All.PotableBothPOCPOU.bray.Region

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.PotableBothPOCPOU.meta$Region)

ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Region<-pairwise.adonis2(All.PotableBothPOCPOU.dist.bray ~ Region, 
                                                                      data = All.PotableBothPOCPOU.meta, 
                                                                      permutations = permutations)
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Region
#Region





#################
#################
#################
#################
#Sample_Local 
ANOSIM.All.PotableBothPOCPOU.bray.Sample_Local<-anosim(All.PotableBothPOCPOU.dist.bray,
                                                        All.PotableBothPOCPOU.meta$Sample_Local, 
                                                        permutations = permutations)

ADONIS2.All.PotableBothPOCPOU.bray.Sample_Local<-adonis2(All.PotableBothPOCPOU.dist.bray ~ Sample_Local, 
                                                          data = All.PotableBothPOCPOU.meta, 
                                                          permutations = permutations)

Betadis.All.PotableBothPOCPOU.bray.Sample_Local<-anova(betadisper(All.PotableBothPOCPOU.dist.bray,
                                                                   All.PotableBothPOCPOU.meta$Sample_Local))



ANOSIM.Pval.All.PotableBothPOCPOU.bray.Sample_Local<-ANOSIM.All.PotableBothPOCPOU.bray.Sample_Local$signif
ANOSIM.Pval.All.PotableBothPOCPOU.bray.Sample_Local

ADONIS2.Pval.All.PotableBothPOCPOU.bray.Sample_Local<-ADONIS2.All.PotableBothPOCPOU.bray.Sample_Local$`Pr(>F)`[1]
ADONIS2.Pval.All.PotableBothPOCPOU.bray.Sample_Local  

Betadis.Pval.All.PotableBothPOCPOU.bray.Sample_Local<-Betadis.All.PotableBothPOCPOU.bray.Sample_Local$`Pr(>F)`[1]
Betadis.Pval.All.PotableBothPOCPOU.bray.Sample_Local

ANOSIM.R.All.PotableBothPOCPOU.bray.Sample_Local<-ANOSIM.All.PotableBothPOCPOU.bray.Sample_Local$statistic
ANOSIM.R.All.PotableBothPOCPOU.bray.Sample_Local

ADONIS2.R2.All.PotableBothPOCPOU.bray.Sample_Local<-ADONIS2.All.PotableBothPOCPOU.bray.Sample_Local$R2[1]
ADONIS2.R2.All.PotableBothPOCPOU.bray.Sample_Local  

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.PotableBothPOCPOU.bray["Sample_Local"]<-ANOSIM.Pval.All.PotableBothPOCPOU.bray.Sample_Local
ADONIS2.GroupedP.All.PotableBothPOCPOU.bray["Sample_Local"]<-ADONIS2.Pval.All.PotableBothPOCPOU.bray.Sample_Local

Betadis.GroupedP.All.PotableBothPOCPOU.bray["Sample_Local"]<-Betadis.Pval.All.PotableBothPOCPOU.bray.Sample_Local

ANOSIM.GroupedR.All.PotableBothPOCPOU.bray["Sample_Local"]<-ANOSIM.R.All.PotableBothPOCPOU.bray.Sample_Local
ADONIS2.GroupedR.All.PotableBothPOCPOU.bray["Sample_Local"]<-ADONIS2.R2.All.PotableBothPOCPOU.bray.Sample_Local

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.PotableBothPOCPOU.meta$Sample_Local)

ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Sample_Local<-pairwise.adonis2(All.PotableBothPOCPOU.dist.bray ~ Sample_Local, 
                                                                            data = All.PotableBothPOCPOU.meta, 
                                                                            permutations = permutations)
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Sample_Local
#Sample_Local





#################
#################
#################
#################
#Final_Disinfection_Residual 
ANOSIM.All.PotableBothPOCPOU.bray.Final_Disinfection_Residual<-anosim(All.PotableBothPOCPOU.dist.bray,
                                                                       All.PotableBothPOCPOU.meta$Final_Disinfection_Residual, 
                                                                       permutations = permutations)

ADONIS2.All.PotableBothPOCPOU.bray.Final_Disinfection_Residual<-adonis2(All.PotableBothPOCPOU.dist.bray ~ Final_Disinfection_Residual, 
                                                                         data = All.PotableBothPOCPOU.meta, 
                                                                         permutations = permutations)

Betadis.All.PotableBothPOCPOU.bray.Final_Disinfection_Residual<-anova(betadisper(All.PotableBothPOCPOU.dist.bray,
                                                                                  All.PotableBothPOCPOU.meta$Final_Disinfection_Residual))



ANOSIM.Pval.All.PotableBothPOCPOU.bray.Final_Disinfection_Residual<-ANOSIM.All.PotableBothPOCPOU.bray.Final_Disinfection_Residual$signif
ANOSIM.Pval.All.PotableBothPOCPOU.bray.Final_Disinfection_Residual

ADONIS2.Pval.All.PotableBothPOCPOU.bray.Final_Disinfection_Residual<-ADONIS2.All.PotableBothPOCPOU.bray.Final_Disinfection_Residual$`Pr(>F)`[1]
ADONIS2.Pval.All.PotableBothPOCPOU.bray.Final_Disinfection_Residual  

Betadis.Pval.All.PotableBothPOCPOU.bray.Final_Disinfection_Residual<-Betadis.All.PotableBothPOCPOU.bray.Final_Disinfection_Residual$`Pr(>F)`[1]
Betadis.Pval.All.PotableBothPOCPOU.bray.Final_Disinfection_Residual

ANOSIM.R.All.PotableBothPOCPOU.bray.Final_Disinfection_Residual<-ANOSIM.All.PotableBothPOCPOU.bray.Final_Disinfection_Residual$statistic
ANOSIM.R.All.PotableBothPOCPOU.bray.Final_Disinfection_Residual

ADONIS2.R2.All.PotableBothPOCPOU.bray.Final_Disinfection_Residual<-ADONIS2.All.PotableBothPOCPOU.bray.Final_Disinfection_Residual$R2[1]
ADONIS2.R2.All.PotableBothPOCPOU.bray.Final_Disinfection_Residual  

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.PotableBothPOCPOU.bray["Final_Disinfection_Residual"]<-ANOSIM.Pval.All.PotableBothPOCPOU.bray.Final_Disinfection_Residual
ADONIS2.GroupedP.All.PotableBothPOCPOU.bray["Final_Disinfection_Residual"]<-ADONIS2.Pval.All.PotableBothPOCPOU.bray.Final_Disinfection_Residual

Betadis.GroupedP.All.PotableBothPOCPOU.bray["Final_Disinfection_Residual"]<-Betadis.Pval.All.PotableBothPOCPOU.bray.Final_Disinfection_Residual

ANOSIM.GroupedR.All.PotableBothPOCPOU.bray["Final_Disinfection_Residual"]<-ANOSIM.R.All.PotableBothPOCPOU.bray.Final_Disinfection_Residual
ADONIS2.GroupedR.All.PotableBothPOCPOU.bray["Final_Disinfection_Residual"]<-ADONIS2.R2.All.PotableBothPOCPOU.bray.Final_Disinfection_Residual

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.PotableBothPOCPOU.meta$Final_Disinfection_Residual)

ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Final_Disinfection_Residual<-pairwise.adonis2(All.PotableBothPOCPOU.dist.bray ~ Final_Disinfection_Residual, 
                                                                                           data = All.PotableBothPOCPOU.meta, 
                                                                                           permutations = permutations)
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Final_Disinfection_Residual
#Final_Disinfection_Residual





# #################
# #################
# #################
# #################
# #Primers 
# ANOSIM.All.PotableBothPOCPOU.bray.Primers<-anosim(All.PotableBothPOCPOU.dist.bray,
#                                                    All.PotableBothPOCPOU.meta$Primers, 
#                                                    permutations = permutations)
# 
# ADONIS2.All.PotableBothPOCPOU.bray.Primers<-adonis2(All.PotableBothPOCPOU.dist.bray ~ Primers, 
#                                                      data = All.PotableBothPOCPOU.meta, 
#                                                      permutations = permutations)
# 
# Betadis.All.PotableBothPOCPOU.bray.Primers<-anova(betadisper(All.PotableBothPOCPOU.dist.bray,
#                                                               All.PotableBothPOCPOU.meta$Primers))
# 
# 
# 
# ANOSIM.Pval.All.PotableBothPOCPOU.bray.Primers<-ANOSIM.All.PotableBothPOCPOU.bray.Primers$signif
# ANOSIM.Pval.All.PotableBothPOCPOU.bray.Primers
# 
# ADONIS2.Pval.All.PotableBothPOCPOU.bray.Primers<-ADONIS2.All.PotableBothPOCPOU.bray.Primers$`Pr(>F)`[1]
# ADONIS2.Pval.All.PotableBothPOCPOU.bray.Primers  
# 
# Betadis.Pval.All.PotableBothPOCPOU.bray.Primers<-Betadis.All.PotableBothPOCPOU.bray.Primers$`Pr(>F)`[1]
# Betadis.Pval.All.PotableBothPOCPOU.bray.Primers
# 
# ANOSIM.R.All.PotableBothPOCPOU.bray.Primers<-ANOSIM.All.PotableBothPOCPOU.bray.Primers$statistic
# ANOSIM.R.All.PotableBothPOCPOU.bray.Primers
# 
# ADONIS2.R2.All.PotableBothPOCPOU.bray.Primers<-ADONIS2.All.PotableBothPOCPOU.bray.Primers$R2[1]
# ADONIS2.R2.All.PotableBothPOCPOU.bray.Primers  
# 
# #Pull Together Pvalues to adjust 
# ANOSIM.GroupedP.All.PotableBothPOCPOU.bray["Primers"]<-ANOSIM.Pval.All.PotableBothPOCPOU.bray.Primers
# ADONIS2.GroupedP.All.PotableBothPOCPOU.bray["Primers"]<-ADONIS2.Pval.All.PotableBothPOCPOU.bray.Primers
# 
# Betadis.GroupedP.All.PotableBothPOCPOU.bray["Primers"]<-Betadis.Pval.All.PotableBothPOCPOU.bray.Primers
# 
# ANOSIM.GroupedR.All.PotableBothPOCPOU.bray["Primers"]<-ANOSIM.R.All.PotableBothPOCPOU.bray.Primers
# ADONIS2.GroupedR.All.PotableBothPOCPOU.bray["Primers"]<-ADONIS2.R2.All.PotableBothPOCPOU.bray.Primers
# 
# #Pairwise if needed, check if interesting or not. Remove if not. 
# count(All.PotableBothPOCPOU.meta$Primers)
# 
# ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Primers<-pairwise.adonis2(All.PotableBothPOCPOU.dist.bray ~ Primers, 
#                                                                        data = All.PotableBothPOCPOU.meta, 
#                                                                        permutations = permutations)
# ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Primers
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
# ANOSIM.All.PotableBothPOCPOU.bray.Trimmed<-anosim(All.PotableBothPOCPOU.dist.bray,
#                                                    All.PotableBothPOCPOU.meta$Trimmed, 
#                                                    permutations = permutations)
# 
# ADONIS2.All.PotableBothPOCPOU.bray.Trimmed<-adonis2(All.PotableBothPOCPOU.dist.bray ~ Trimmed, 
#                                                      data = All.PotableBothPOCPOU.meta, 
#                                                      permutations = permutations)
# 
# Betadis.All.PotableBothPOCPOU.bray.Trimmed<-anova(betadisper(All.PotableBothPOCPOU.dist.bray,
#                                                               All.PotableBothPOCPOU.meta$Trimmed))
# 
# 
# 
# ANOSIM.Pval.All.PotableBothPOCPOU.bray.Trimmed<-ANOSIM.All.PotableBothPOCPOU.bray.Trimmed$signif
# ANOSIM.Pval.All.PotableBothPOCPOU.bray.Trimmed
# 
# ADONIS2.Pval.All.PotableBothPOCPOU.bray.Trimmed<-ADONIS2.All.PotableBothPOCPOU.bray.Trimmed$`Pr(>F)`[1]
# ADONIS2.Pval.All.PotableBothPOCPOU.bray.Trimmed  
# 
# Betadis.Pval.All.PotableBothPOCPOU.bray.Trimmed<-Betadis.All.PotableBothPOCPOU.bray.Trimmed$`Pr(>F)`[1]
# Betadis.Pval.All.PotableBothPOCPOU.bray.Trimmed
# 
# ANOSIM.R.All.PotableBothPOCPOU.bray.Trimmed<-ANOSIM.All.PotableBothPOCPOU.bray.Trimmed$statistic
# ANOSIM.R.All.PotableBothPOCPOU.bray.Trimmed
# 
# ADONIS2.R2.All.PotableBothPOCPOU.bray.Trimmed<-ADONIS2.All.PotableBothPOCPOU.bray.Trimmed$R2[1]
# ADONIS2.R2.All.PotableBothPOCPOU.bray.Trimmed  
# 
# #Pull Together Pvalues to adjust 
# ANOSIM.GroupedP.All.PotableBothPOCPOU.bray["Trimmed"]<-ANOSIM.Pval.All.PotableBothPOCPOU.bray.Trimmed
# ADONIS2.GroupedP.All.PotableBothPOCPOU.bray["Trimmed"]<-ADONIS2.Pval.All.PotableBothPOCPOU.bray.Trimmed
# 
# Betadis.GroupedP.All.PotableBothPOCPOU.bray["Trimmed"]<-Betadis.Pval.All.PotableBothPOCPOU.bray.Trimmed
# 
# ANOSIM.GroupedR.All.PotableBothPOCPOU.bray["Trimmed"]<-ANOSIM.R.All.PotableBothPOCPOU.bray.Trimmed
# ADONIS2.GroupedR.All.PotableBothPOCPOU.bray["Trimmed"]<-ADONIS2.R2.All.PotableBothPOCPOU.bray.Trimmed
# 
# #Pairwise if needed, check if interesting or not. Remove if not. 
# count(All.PotableBothPOCPOU.meta$Trimmed)
# 
# ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Trimmed<-pairwise.adonis2(All.PotableBothPOCPOU.dist.bray ~ Trimmed, 
#                                                                        data = All.PotableBothPOCPOU.meta, 
#                                                                        permutations = permutations)
# ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Trimmed
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
# ANOSIM.All.PotableBothPOCPOU.bray.Extraction_Type<-anosim(All.PotableBothPOCPOU.dist.bray,
#                                                            All.PotableBothPOCPOU.meta$Extraction_Type, 
#                                                            permutations = permutations)
# 
# ADONIS2.All.PotableBothPOCPOU.bray.Extraction_Type<-adonis2(All.PotableBothPOCPOU.dist.bray ~ Extraction_Type, 
#                                                              data = All.PotableBothPOCPOU.meta, 
#                                                              permutations = permutations)
# 
# Betadis.All.PotableBothPOCPOU.bray.Extraction_Type<-anova(betadisper(All.PotableBothPOCPOU.dist.bray,
#                                                                       All.PotableBothPOCPOU.meta$Extraction_Type))
# 
# 
# 
# ANOSIM.Pval.All.PotableBothPOCPOU.bray.Extraction_Type<-ANOSIM.All.PotableBothPOCPOU.bray.Extraction_Type$signif
# ANOSIM.Pval.All.PotableBothPOCPOU.bray.Extraction_Type
# 
# ADONIS2.Pval.All.PotableBothPOCPOU.bray.Extraction_Type<-ADONIS2.All.PotableBothPOCPOU.bray.Extraction_Type$`Pr(>F)`[1]
# ADONIS2.Pval.All.PotableBothPOCPOU.bray.Extraction_Type  
# 
# Betadis.Pval.All.PotableBothPOCPOU.bray.Extraction_Type<-Betadis.All.PotableBothPOCPOU.bray.Extraction_Type$`Pr(>F)`[1]
# Betadis.Pval.All.PotableBothPOCPOU.bray.Extraction_Type
# 
# ANOSIM.R.All.PotableBothPOCPOU.bray.Extraction_Type<-ANOSIM.All.PotableBothPOCPOU.bray.Extraction_Type$statistic
# ANOSIM.R.All.PotableBothPOCPOU.bray.Extraction_Type
# 
# ADONIS2.R2.All.PotableBothPOCPOU.bray.Extraction_Type<-ADONIS2.All.PotableBothPOCPOU.bray.Extraction_Type$R2[1]
# ADONIS2.R2.All.PotableBothPOCPOU.bray.Extraction_Type  
# 
# #Pull Together Pvalues to adjust 
# ANOSIM.GroupedP.All.PotableBothPOCPOU.bray["Extraction_Type"]<-ANOSIM.Pval.All.PotableBothPOCPOU.bray.Extraction_Type
# ADONIS2.GroupedP.All.PotableBothPOCPOU.bray["Extraction_Type"]<-ADONIS2.Pval.All.PotableBothPOCPOU.bray.Extraction_Type
# 
# Betadis.GroupedP.All.PotableBothPOCPOU.bray["Extraction_Type"]<-Betadis.Pval.All.PotableBothPOCPOU.bray.Extraction_Type
# 
# ANOSIM.GroupedR.All.PotableBothPOCPOU.bray["Extraction_Type"]<-ANOSIM.R.All.PotableBothPOCPOU.bray.Extraction_Type
# ADONIS2.GroupedR.All.PotableBothPOCPOU.bray["Extraction_Type"]<-ADONIS2.R2.All.PotableBothPOCPOU.bray.Extraction_Type
# 
# #Pairwise if needed, check if interesting or not. Remove if not. 
# count(All.PotableBothPOCPOU.meta$Extraction_Type)
# 
# ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Extraction_Type<-pairwise.adonis2(All.PotableBothPOCPOU.dist.bray ~ Extraction_Type, 
#                                                                                data = All.PotableBothPOCPOU.meta, 
#                                                                                permutations = permutations)
# ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Extraction_Type
# #Extraction_Type





#################
#################
#################
#################
#Potable_Source_Water 
ANOSIM.All.PotableBothPOCPOU.bray.Potable_Source_Water<-anosim(All.PotableBothPOCPOU.dist.bray,
                                                                All.PotableBothPOCPOU.meta$Potable_Source_Water, 
                                                                permutations = permutations)

ADONIS2.All.PotableBothPOCPOU.bray.Potable_Source_Water<-adonis2(All.PotableBothPOCPOU.dist.bray ~ Potable_Source_Water, 
                                                                  data = All.PotableBothPOCPOU.meta, 
                                                                  permutations = permutations)

Betadis.All.PotableBothPOCPOU.bray.Potable_Source_Water<-anova(betadisper(All.PotableBothPOCPOU.dist.bray,
                                                                           All.PotableBothPOCPOU.meta$Potable_Source_Water))



ANOSIM.Pval.All.PotableBothPOCPOU.bray.Potable_Source_Water<-ANOSIM.All.PotableBothPOCPOU.bray.Potable_Source_Water$signif
ANOSIM.Pval.All.PotableBothPOCPOU.bray.Potable_Source_Water

ADONIS2.Pval.All.PotableBothPOCPOU.bray.Potable_Source_Water<-ADONIS2.All.PotableBothPOCPOU.bray.Potable_Source_Water$`Pr(>F)`[1]
ADONIS2.Pval.All.PotableBothPOCPOU.bray.Potable_Source_Water  

Betadis.Pval.All.PotableBothPOCPOU.bray.Potable_Source_Water<-Betadis.All.PotableBothPOCPOU.bray.Potable_Source_Water$`Pr(>F)`[1]
Betadis.Pval.All.PotableBothPOCPOU.bray.Potable_Source_Water

ANOSIM.R.All.PotableBothPOCPOU.bray.Potable_Source_Water<-ANOSIM.All.PotableBothPOCPOU.bray.Potable_Source_Water$statistic
ANOSIM.R.All.PotableBothPOCPOU.bray.Potable_Source_Water

ADONIS2.R2.All.PotableBothPOCPOU.bray.Potable_Source_Water<-ADONIS2.All.PotableBothPOCPOU.bray.Potable_Source_Water$R2[1]
ADONIS2.R2.All.PotableBothPOCPOU.bray.Potable_Source_Water  

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.PotableBothPOCPOU.bray["Potable_Source_Water"]<-ANOSIM.Pval.All.PotableBothPOCPOU.bray.Potable_Source_Water
ADONIS2.GroupedP.All.PotableBothPOCPOU.bray["Potable_Source_Water"]<-ADONIS2.Pval.All.PotableBothPOCPOU.bray.Potable_Source_Water

Betadis.GroupedP.All.PotableBothPOCPOU.bray["Potable_Source_Water"]<-Betadis.Pval.All.PotableBothPOCPOU.bray.Potable_Source_Water

ANOSIM.GroupedR.All.PotableBothPOCPOU.bray["Potable_Source_Water"]<-ANOSIM.R.All.PotableBothPOCPOU.bray.Potable_Source_Water
ADONIS2.GroupedR.All.PotableBothPOCPOU.bray["Potable_Source_Water"]<-ADONIS2.R2.All.PotableBothPOCPOU.bray.Potable_Source_Water

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.PotableBothPOCPOU.meta$Potable_Source_Water)

ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Potable_Source_Water<-pairwise.adonis2(All.PotableBothPOCPOU.dist.bray ~ Potable_Source_Water, 
                                                                                    data = All.PotableBothPOCPOU.meta, 
                                                                                    permutations = permutations)
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Potable_Source_Water
#Potable_Source_Water
##Potable_Source_Water




#################
#################
#################
#################
#Potable_Disinfection_Primary  
ANOSIM.All.PotableBothPOCPOU.bray.Potable_Disinfection_Primary <-anosim(All.PotableBothPOCPOU.dist.bray,
                                                                         All.PotableBothPOCPOU.meta$Potable_Disinfection_Primary , 
                                                                         permutations = permutations)

ADONIS2.All.PotableBothPOCPOU.bray.Potable_Disinfection_Primary <-adonis2(All.PotableBothPOCPOU.dist.bray ~ Potable_Disinfection_Primary , 
                                                                           data = All.PotableBothPOCPOU.meta, 
                                                                           permutations = permutations)

Betadis.All.PotableBothPOCPOU.bray.Potable_Disinfection_Primary <-anova(betadisper(All.PotableBothPOCPOU.dist.bray,
                                                                                    All.PotableBothPOCPOU.meta$Potable_Disinfection_Primary ))



ANOSIM.Pval.All.PotableBothPOCPOU.bray.Potable_Disinfection_Primary <-ANOSIM.All.PotableBothPOCPOU.bray.Potable_Disinfection_Primary $signif
ANOSIM.Pval.All.PotableBothPOCPOU.bray.Potable_Disinfection_Primary 

ADONIS2.Pval.All.PotableBothPOCPOU.bray.Potable_Disinfection_Primary <-ADONIS2.All.PotableBothPOCPOU.bray.Potable_Disinfection_Primary $`Pr(>F)`[1]
ADONIS2.Pval.All.PotableBothPOCPOU.bray.Potable_Disinfection_Primary   

Betadis.Pval.All.PotableBothPOCPOU.bray.Potable_Disinfection_Primary <-Betadis.All.PotableBothPOCPOU.bray.Potable_Disinfection_Primary $`Pr(>F)`[1]
Betadis.Pval.All.PotableBothPOCPOU.bray.Potable_Disinfection_Primary 

ANOSIM.R.All.PotableBothPOCPOU.bray.Potable_Disinfection_Primary <-ANOSIM.All.PotableBothPOCPOU.bray.Potable_Disinfection_Primary $statistic
ANOSIM.R.All.PotableBothPOCPOU.bray.Potable_Disinfection_Primary 

ADONIS2.R2.All.PotableBothPOCPOU.bray.Potable_Disinfection_Primary <-ADONIS2.All.PotableBothPOCPOU.bray.Potable_Disinfection_Primary $R2[1]
ADONIS2.R2.All.PotableBothPOCPOU.bray.Potable_Disinfection_Primary   

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.PotableBothPOCPOU.bray["Potable_Disinfection_Primary "]<-ANOSIM.Pval.All.PotableBothPOCPOU.bray.Potable_Disinfection_Primary 
ADONIS2.GroupedP.All.PotableBothPOCPOU.bray["Potable_Disinfection_Primary "]<-ADONIS2.Pval.All.PotableBothPOCPOU.bray.Potable_Disinfection_Primary 

Betadis.GroupedP.All.PotableBothPOCPOU.bray["Potable_Disinfection_Primary "]<-Betadis.Pval.All.PotableBothPOCPOU.bray.Potable_Disinfection_Primary 

ANOSIM.GroupedR.All.PotableBothPOCPOU.bray["Potable_Disinfection_Primary "]<-ANOSIM.R.All.PotableBothPOCPOU.bray.Potable_Disinfection_Primary 
ADONIS2.GroupedR.All.PotableBothPOCPOU.bray["Potable_Disinfection_Primary "]<-ADONIS2.R2.All.PotableBothPOCPOU.bray.Potable_Disinfection_Primary 

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.PotableBothPOCPOU.meta$Potable_Disinfection_Primary )

ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Potable_Disinfection_Primary <-pairwise.adonis2(All.PotableBothPOCPOU.dist.bray ~ Potable_Disinfection_Primary , 
                                                                                             data = All.PotableBothPOCPOU.meta, 
                                                                                             permutations = permutations)
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Potable_Disinfection_Primary 
#Potable_Disinfection_Primary 
##Potable_Disinfection_Primary 




#################
#################
#################
#################
#Potable_Treatment_Type 
ANOSIM.All.PotableBothPOCPOU.bray.Potable_Treatment_Type<-anosim(All.PotableBothPOCPOU.dist.bray,
                                                                  All.PotableBothPOCPOU.meta$Potable_Treatment_Type, 
                                                                  permutations = permutations)

ADONIS2.All.PotableBothPOCPOU.bray.Potable_Treatment_Type<-adonis2(All.PotableBothPOCPOU.dist.bray ~ Potable_Treatment_Type, 
                                                                    data = All.PotableBothPOCPOU.meta, 
                                                                    permutations = permutations)

Betadis.All.PotableBothPOCPOU.bray.Potable_Treatment_Type<-anova(betadisper(All.PotableBothPOCPOU.dist.bray,
                                                                             All.PotableBothPOCPOU.meta$Potable_Treatment_Type))



ANOSIM.Pval.All.PotableBothPOCPOU.bray.Potable_Treatment_Type<-ANOSIM.All.PotableBothPOCPOU.bray.Potable_Treatment_Type$signif
ANOSIM.Pval.All.PotableBothPOCPOU.bray.Potable_Treatment_Type

ADONIS2.Pval.All.PotableBothPOCPOU.bray.Potable_Treatment_Type<-ADONIS2.All.PotableBothPOCPOU.bray.Potable_Treatment_Type$`Pr(>F)`[1]
ADONIS2.Pval.All.PotableBothPOCPOU.bray.Potable_Treatment_Type  

Betadis.Pval.All.PotableBothPOCPOU.bray.Potable_Treatment_Type<-Betadis.All.PotableBothPOCPOU.bray.Potable_Treatment_Type$`Pr(>F)`[1]
Betadis.Pval.All.PotableBothPOCPOU.bray.Potable_Treatment_Type

ANOSIM.R.All.PotableBothPOCPOU.bray.Potable_Treatment_Type<-ANOSIM.All.PotableBothPOCPOU.bray.Potable_Treatment_Type$statistic
ANOSIM.R.All.PotableBothPOCPOU.bray.Potable_Treatment_Type

ADONIS2.R2.All.PotableBothPOCPOU.bray.Potable_Treatment_Type<-ADONIS2.All.PotableBothPOCPOU.bray.Potable_Treatment_Type$R2[1]
ADONIS2.R2.All.PotableBothPOCPOU.bray.Potable_Treatment_Type  

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.PotableBothPOCPOU.bray["Potable_Treatment_Type"]<-ANOSIM.Pval.All.PotableBothPOCPOU.bray.Potable_Treatment_Type
ADONIS2.GroupedP.All.PotableBothPOCPOU.bray["Potable_Treatment_Type"]<-ADONIS2.Pval.All.PotableBothPOCPOU.bray.Potable_Treatment_Type

Betadis.GroupedP.All.PotableBothPOCPOU.bray["Potable_Treatment_Type"]<-Betadis.Pval.All.PotableBothPOCPOU.bray.Potable_Treatment_Type

ANOSIM.GroupedR.All.PotableBothPOCPOU.bray["Potable_Treatment_Type"]<-ANOSIM.R.All.PotableBothPOCPOU.bray.Potable_Treatment_Type
ADONIS2.GroupedR.All.PotableBothPOCPOU.bray["Potable_Treatment_Type"]<-ADONIS2.R2.All.PotableBothPOCPOU.bray.Potable_Treatment_Type

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.PotableBothPOCPOU.meta$Potable_Treatment_Type)

ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Potable_Treatment_Type<-pairwise.adonis2(All.PotableBothPOCPOU.dist.bray ~ Potable_Treatment_Type, 
                                                                                      data = All.PotableBothPOCPOU.meta, 
                                                                                      permutations = permutations)
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Potable_Treatment_Type
#Potable_Treatment_Type
##Potable_Treatment_Type 




#################
#################
#################
#################
#Disinfection_Total 
ANOSIM.All.PotableBothPOCPOU.bray.Disinfection_Total<-anosim(All.PotableBothPOCPOU.dist.bray,
                                                              All.PotableBothPOCPOU.meta$Disinfection_Total, 
                                                              permutations = permutations)

ADONIS2.All.PotableBothPOCPOU.bray.Disinfection_Total<-adonis2(All.PotableBothPOCPOU.dist.bray ~ Disinfection_Total, 
                                                                data = All.PotableBothPOCPOU.meta, 
                                                                permutations = permutations)

Betadis.All.PotableBothPOCPOU.bray.Disinfection_Total<-anova(betadisper(All.PotableBothPOCPOU.dist.bray,
                                                                         All.PotableBothPOCPOU.meta$Disinfection_Total))



ANOSIM.Pval.All.PotableBothPOCPOU.bray.Disinfection_Total<-ANOSIM.All.PotableBothPOCPOU.bray.Disinfection_Total$signif
ANOSIM.Pval.All.PotableBothPOCPOU.bray.Disinfection_Total

ADONIS2.Pval.All.PotableBothPOCPOU.bray.Disinfection_Total<-ADONIS2.All.PotableBothPOCPOU.bray.Disinfection_Total$`Pr(>F)`[1]
ADONIS2.Pval.All.PotableBothPOCPOU.bray.Disinfection_Total  

Betadis.Pval.All.PotableBothPOCPOU.bray.Disinfection_Total<-Betadis.All.PotableBothPOCPOU.bray.Disinfection_Total$`Pr(>F)`[1]
Betadis.Pval.All.PotableBothPOCPOU.bray.Disinfection_Total

ANOSIM.R.All.PotableBothPOCPOU.bray.Disinfection_Total<-ANOSIM.All.PotableBothPOCPOU.bray.Disinfection_Total$statistic
ANOSIM.R.All.PotableBothPOCPOU.bray.Disinfection_Total

ADONIS2.R2.All.PotableBothPOCPOU.bray.Disinfection_Total<-ADONIS2.All.PotableBothPOCPOU.bray.Disinfection_Total$R2[1]
ADONIS2.R2.All.PotableBothPOCPOU.bray.Disinfection_Total  

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.PotableBothPOCPOU.bray["Disinfection_Total"]<-ANOSIM.Pval.All.PotableBothPOCPOU.bray.Disinfection_Total
ADONIS2.GroupedP.All.PotableBothPOCPOU.bray["Disinfection_Total"]<-ADONIS2.Pval.All.PotableBothPOCPOU.bray.Disinfection_Total

Betadis.GroupedP.All.PotableBothPOCPOU.bray["Disinfection_Total"]<-Betadis.Pval.All.PotableBothPOCPOU.bray.Disinfection_Total

ANOSIM.GroupedR.All.PotableBothPOCPOU.bray["Disinfection_Total"]<-ANOSIM.R.All.PotableBothPOCPOU.bray.Disinfection_Total
ADONIS2.GroupedR.All.PotableBothPOCPOU.bray["Disinfection_Total"]<-ADONIS2.R2.All.PotableBothPOCPOU.bray.Disinfection_Total

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.PotableBothPOCPOU.meta$Disinfection_Total)

ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Disinfection_Total<-pairwise.adonis2(All.PotableBothPOCPOU.dist.bray ~ Disinfection_Total, 
                                                                                  data = All.PotableBothPOCPOU.meta, 
                                                                                  permutations = permutations)
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Disinfection_Total
#Disinfection_Total
##Disinfection_Total 




#################
#################
#################
#################
#Treatment_Total  
ANOSIM.All.PotableBothPOCPOU.bray.Treatment_Total <-anosim(All.PotableBothPOCPOU.dist.bray,
                                                            All.PotableBothPOCPOU.meta$Treatment_Total , 
                                                            permutations = permutations)

ADONIS2.All.PotableBothPOCPOU.bray.Treatment_Total <-adonis2(All.PotableBothPOCPOU.dist.bray ~ Treatment_Total , 
                                                              data = All.PotableBothPOCPOU.meta, 
                                                              permutations = permutations)

Betadis.All.PotableBothPOCPOU.bray.Treatment_Total <-anova(betadisper(All.PotableBothPOCPOU.dist.bray,
                                                                       All.PotableBothPOCPOU.meta$Treatment_Total ))



ANOSIM.Pval.All.PotableBothPOCPOU.bray.Treatment_Total <-ANOSIM.All.PotableBothPOCPOU.bray.Treatment_Total $signif
ANOSIM.Pval.All.PotableBothPOCPOU.bray.Treatment_Total 

ADONIS2.Pval.All.PotableBothPOCPOU.bray.Treatment_Total <-ADONIS2.All.PotableBothPOCPOU.bray.Treatment_Total $`Pr(>F)`[1]
ADONIS2.Pval.All.PotableBothPOCPOU.bray.Treatment_Total   

Betadis.Pval.All.PotableBothPOCPOU.bray.Treatment_Total <-Betadis.All.PotableBothPOCPOU.bray.Treatment_Total $`Pr(>F)`[1]
Betadis.Pval.All.PotableBothPOCPOU.bray.Treatment_Total 

ANOSIM.R.All.PotableBothPOCPOU.bray.Treatment_Total <-ANOSIM.All.PotableBothPOCPOU.bray.Treatment_Total $statistic
ANOSIM.R.All.PotableBothPOCPOU.bray.Treatment_Total 

ADONIS2.R2.All.PotableBothPOCPOU.bray.Treatment_Total <-ADONIS2.All.PotableBothPOCPOU.bray.Treatment_Total $R2[1]
ADONIS2.R2.All.PotableBothPOCPOU.bray.Treatment_Total   

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.PotableBothPOCPOU.bray["Treatment_Total "]<-ANOSIM.Pval.All.PotableBothPOCPOU.bray.Treatment_Total 
ADONIS2.GroupedP.All.PotableBothPOCPOU.bray["Treatment_Total "]<-ADONIS2.Pval.All.PotableBothPOCPOU.bray.Treatment_Total 

Betadis.GroupedP.All.PotableBothPOCPOU.bray["Treatment_Total "]<-Betadis.Pval.All.PotableBothPOCPOU.bray.Treatment_Total 

ANOSIM.GroupedR.All.PotableBothPOCPOU.bray["Treatment_Total "]<-ANOSIM.R.All.PotableBothPOCPOU.bray.Treatment_Total 
ADONIS2.GroupedR.All.PotableBothPOCPOU.bray["Treatment_Total "]<-ADONIS2.R2.All.PotableBothPOCPOU.bray.Treatment_Total 

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.PotableBothPOCPOU.meta$Treatment_Total )

ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Treatment_Total <-pairwise.adonis2(All.PotableBothPOCPOU.dist.bray ~ Treatment_Total , 
                                                                                data = All.PotableBothPOCPOU.meta, 
                                                                                permutations = permutations)
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Treatment_Total 
#Treatment_Total 
##Treatment_Total 




#################
#################
#################
#################
#Treatment_Classification_Number 
ANOSIM.All.PotableBothPOCPOU.bray.Treatment_Classification_Number<-anosim(All.PotableBothPOCPOU.dist.bray,
                                                                           All.PotableBothPOCPOU.meta$Treatment_Classification_Number, 
                                                                           permutations = permutations)

ADONIS2.All.PotableBothPOCPOU.bray.Treatment_Classification_Number<-adonis2(All.PotableBothPOCPOU.dist.bray ~ Treatment_Classification_Number, 
                                                                             data = All.PotableBothPOCPOU.meta, 
                                                                             permutations = permutations)

Betadis.All.PotableBothPOCPOU.bray.Treatment_Classification_Number<-anova(betadisper(All.PotableBothPOCPOU.dist.bray,
                                                                                      All.PotableBothPOCPOU.meta$Treatment_Classification_Number))



ANOSIM.Pval.All.PotableBothPOCPOU.bray.Treatment_Classification_Number<-ANOSIM.All.PotableBothPOCPOU.bray.Treatment_Classification_Number$signif
ANOSIM.Pval.All.PotableBothPOCPOU.bray.Treatment_Classification_Number

ADONIS2.Pval.All.PotableBothPOCPOU.bray.Treatment_Classification_Number<-ADONIS2.All.PotableBothPOCPOU.bray.Treatment_Classification_Number$`Pr(>F)`[1]
ADONIS2.Pval.All.PotableBothPOCPOU.bray.Treatment_Classification_Number  

Betadis.Pval.All.PotableBothPOCPOU.bray.Treatment_Classification_Number<-Betadis.All.PotableBothPOCPOU.bray.Treatment_Classification_Number$`Pr(>F)`[1]
Betadis.Pval.All.PotableBothPOCPOU.bray.Treatment_Classification_Number

ANOSIM.R.All.PotableBothPOCPOU.bray.Treatment_Classification_Number<-ANOSIM.All.PotableBothPOCPOU.bray.Treatment_Classification_Number$statistic
ANOSIM.R.All.PotableBothPOCPOU.bray.Treatment_Classification_Number

ADONIS2.R2.All.PotableBothPOCPOU.bray.Treatment_Classification_Number<-ADONIS2.All.PotableBothPOCPOU.bray.Treatment_Classification_Number$R2[1]
ADONIS2.R2.All.PotableBothPOCPOU.bray.Treatment_Classification_Number  

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.PotableBothPOCPOU.bray["Treatment_Classification_Number"]<-ANOSIM.Pval.All.PotableBothPOCPOU.bray.Treatment_Classification_Number
ADONIS2.GroupedP.All.PotableBothPOCPOU.bray["Treatment_Classification_Number"]<-ADONIS2.Pval.All.PotableBothPOCPOU.bray.Treatment_Classification_Number

Betadis.GroupedP.All.PotableBothPOCPOU.bray["Treatment_Classification_Number"]<-Betadis.Pval.All.PotableBothPOCPOU.bray.Treatment_Classification_Number

ANOSIM.GroupedR.All.PotableBothPOCPOU.bray["Treatment_Classification_Number"]<-ANOSIM.R.All.PotableBothPOCPOU.bray.Treatment_Classification_Number
ADONIS2.GroupedR.All.PotableBothPOCPOU.bray["Treatment_Classification_Number"]<-ADONIS2.R2.All.PotableBothPOCPOU.bray.Treatment_Classification_Number

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.PotableBothPOCPOU.meta$Treatment_Classification_Number)

ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Treatment_Classification_Number<-pairwise.adonis2(All.PotableBothPOCPOU.dist.bray ~ Treatment_Classification_Number, 
                                                                                               data = All.PotableBothPOCPOU.meta, 
                                                                                               permutations = permutations)
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Treatment_Classification_Number
#Treatment_Classification_Number
##Treatment_Classification_Number




#################
#################
#################
#################
#Treatment_Classification_Cat 
ANOSIM.All.PotableBothPOCPOU.bray.Treatment_Classification_Cat<-anosim(All.PotableBothPOCPOU.dist.bray,
                                                                        All.PotableBothPOCPOU.meta$Treatment_Classification_Cat, 
                                                                        permutations = permutations)

ADONIS2.All.PotableBothPOCPOU.bray.Treatment_Classification_Cat<-adonis2(All.PotableBothPOCPOU.dist.bray ~ Treatment_Classification_Cat, 
                                                                          data = All.PotableBothPOCPOU.meta, 
                                                                          permutations = permutations)

Betadis.All.PotableBothPOCPOU.bray.Treatment_Classification_Cat<-anova(betadisper(All.PotableBothPOCPOU.dist.bray,
                                                                                   All.PotableBothPOCPOU.meta$Treatment_Classification_Cat))



ANOSIM.Pval.All.PotableBothPOCPOU.bray.Treatment_Classification_Cat<-ANOSIM.All.PotableBothPOCPOU.bray.Treatment_Classification_Cat$signif
ANOSIM.Pval.All.PotableBothPOCPOU.bray.Treatment_Classification_Cat

ADONIS2.Pval.All.PotableBothPOCPOU.bray.Treatment_Classification_Cat<-ADONIS2.All.PotableBothPOCPOU.bray.Treatment_Classification_Cat$`Pr(>F)`[1]
ADONIS2.Pval.All.PotableBothPOCPOU.bray.Treatment_Classification_Cat  

Betadis.Pval.All.PotableBothPOCPOU.bray.Treatment_Classification_Cat<-Betadis.All.PotableBothPOCPOU.bray.Treatment_Classification_Cat$`Pr(>F)`[1]
Betadis.Pval.All.PotableBothPOCPOU.bray.Treatment_Classification_Cat

ANOSIM.R.All.PotableBothPOCPOU.bray.Treatment_Classification_Cat<-ANOSIM.All.PotableBothPOCPOU.bray.Treatment_Classification_Cat$statistic
ANOSIM.R.All.PotableBothPOCPOU.bray.Treatment_Classification_Cat

ADONIS2.R2.All.PotableBothPOCPOU.bray.Treatment_Classification_Cat<-ADONIS2.All.PotableBothPOCPOU.bray.Treatment_Classification_Cat$R2[1]
ADONIS2.R2.All.PotableBothPOCPOU.bray.Treatment_Classification_Cat  

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.PotableBothPOCPOU.bray["Treatment_Classification_Cat"]<-ANOSIM.Pval.All.PotableBothPOCPOU.bray.Treatment_Classification_Cat
ADONIS2.GroupedP.All.PotableBothPOCPOU.bray["Treatment_Classification_Cat"]<-ADONIS2.Pval.All.PotableBothPOCPOU.bray.Treatment_Classification_Cat

Betadis.GroupedP.All.PotableBothPOCPOU.bray["Treatment_Classification_Cat"]<-Betadis.Pval.All.PotableBothPOCPOU.bray.Treatment_Classification_Cat

ANOSIM.GroupedR.All.PotableBothPOCPOU.bray["Treatment_Classification_Cat"]<-ANOSIM.R.All.PotableBothPOCPOU.bray.Treatment_Classification_Cat
ADONIS2.GroupedR.All.PotableBothPOCPOU.bray["Treatment_Classification_Cat"]<-ADONIS2.R2.All.PotableBothPOCPOU.bray.Treatment_Classification_Cat

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.PotableBothPOCPOU.meta$Treatment_Classification_Cat)

ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Treatment_Classification_Cat<-pairwise.adonis2(All.PotableBothPOCPOU.dist.bray ~ Treatment_Classification_Cat, 
                                                                                            data = All.PotableBothPOCPOU.meta, 
                                                                                            permutations = permutations)
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Treatment_Classification_Cat
#Treatment_Classification_Cat
##Treatment_Classification_Cat




#################
#################
#################
#################
#Mixed_Source_Water 
ANOSIM.All.PotableBothPOCPOU.bray.Mixed_Source_Water<-anosim(All.PotableBothPOCPOU.dist.bray,
                                                              All.PotableBothPOCPOU.meta$Mixed_Source_Water, 
                                                              permutations = permutations)

ADONIS2.All.PotableBothPOCPOU.bray.Mixed_Source_Water<-adonis2(All.PotableBothPOCPOU.dist.bray ~ Mixed_Source_Water, 
                                                                data = All.PotableBothPOCPOU.meta, 
                                                                permutations = permutations)

Betadis.All.PotableBothPOCPOU.bray.Mixed_Source_Water<-anova(betadisper(All.PotableBothPOCPOU.dist.bray,
                                                                         All.PotableBothPOCPOU.meta$Mixed_Source_Water))



ANOSIM.Pval.All.PotableBothPOCPOU.bray.Mixed_Source_Water<-ANOSIM.All.PotableBothPOCPOU.bray.Mixed_Source_Water$signif
ANOSIM.Pval.All.PotableBothPOCPOU.bray.Mixed_Source_Water

ADONIS2.Pval.All.PotableBothPOCPOU.bray.Mixed_Source_Water<-ADONIS2.All.PotableBothPOCPOU.bray.Mixed_Source_Water$`Pr(>F)`[1]
ADONIS2.Pval.All.PotableBothPOCPOU.bray.Mixed_Source_Water  

Betadis.Pval.All.PotableBothPOCPOU.bray.Mixed_Source_Water<-Betadis.All.PotableBothPOCPOU.bray.Mixed_Source_Water$`Pr(>F)`[1]
Betadis.Pval.All.PotableBothPOCPOU.bray.Mixed_Source_Water

ANOSIM.R.All.PotableBothPOCPOU.bray.Mixed_Source_Water<-ANOSIM.All.PotableBothPOCPOU.bray.Mixed_Source_Water$statistic
ANOSIM.R.All.PotableBothPOCPOU.bray.Mixed_Source_Water

ADONIS2.R2.All.PotableBothPOCPOU.bray.Mixed_Source_Water<-ADONIS2.All.PotableBothPOCPOU.bray.Mixed_Source_Water$R2[1]
ADONIS2.R2.All.PotableBothPOCPOU.bray.Mixed_Source_Water  

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.PotableBothPOCPOU.bray["Mixed_Source_Water"]<-ANOSIM.Pval.All.PotableBothPOCPOU.bray.Mixed_Source_Water
ADONIS2.GroupedP.All.PotableBothPOCPOU.bray["Mixed_Source_Water"]<-ADONIS2.Pval.All.PotableBothPOCPOU.bray.Mixed_Source_Water

Betadis.GroupedP.All.PotableBothPOCPOU.bray["Mixed_Source_Water"]<-Betadis.Pval.All.PotableBothPOCPOU.bray.Mixed_Source_Water

ANOSIM.GroupedR.All.PotableBothPOCPOU.bray["Mixed_Source_Water"]<-ANOSIM.R.All.PotableBothPOCPOU.bray.Mixed_Source_Water
ADONIS2.GroupedR.All.PotableBothPOCPOU.bray["Mixed_Source_Water"]<-ADONIS2.R2.All.PotableBothPOCPOU.bray.Mixed_Source_Water

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.PotableBothPOCPOU.meta$Mixed_Source_Water)

ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Mixed_Source_Water<-pairwise.adonis2(All.PotableBothPOCPOU.dist.bray ~ Mixed_Source_Water, 
                                                                                  data = All.PotableBothPOCPOU.meta, 
                                                                                  permutations = permutations)
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Mixed_Source_Water
#Mixed_Source_Water
##Mixed_Source_Water




#################
#################
#################
#################
#Potable_Overlap 
ANOSIM.All.PotableBothPOCPOU.bray.Potable_Overlap<-anosim(All.PotableBothPOCPOU.dist.bray,
                                                           All.PotableBothPOCPOU.meta$Potable_Overlap, 
                                                           permutations = permutations)

ADONIS2.All.PotableBothPOCPOU.bray.Potable_Overlap<-adonis2(All.PotableBothPOCPOU.dist.bray ~ Potable_Overlap, 
                                                             data = All.PotableBothPOCPOU.meta, 
                                                             permutations = permutations)

Betadis.All.PotableBothPOCPOU.bray.Potable_Overlap<-anova(betadisper(All.PotableBothPOCPOU.dist.bray,
                                                                      All.PotableBothPOCPOU.meta$Potable_Overlap))



ANOSIM.Pval.All.PotableBothPOCPOU.bray.Potable_Overlap<-ANOSIM.All.PotableBothPOCPOU.bray.Potable_Overlap$signif
ANOSIM.Pval.All.PotableBothPOCPOU.bray.Potable_Overlap

ADONIS2.Pval.All.PotableBothPOCPOU.bray.Potable_Overlap<-ADONIS2.All.PotableBothPOCPOU.bray.Potable_Overlap$`Pr(>F)`[1]
ADONIS2.Pval.All.PotableBothPOCPOU.bray.Potable_Overlap  

Betadis.Pval.All.PotableBothPOCPOU.bray.Potable_Overlap<-Betadis.All.PotableBothPOCPOU.bray.Potable_Overlap$`Pr(>F)`[1]
Betadis.Pval.All.PotableBothPOCPOU.bray.Potable_Overlap

ANOSIM.R.All.PotableBothPOCPOU.bray.Potable_Overlap<-ANOSIM.All.PotableBothPOCPOU.bray.Potable_Overlap$statistic
ANOSIM.R.All.PotableBothPOCPOU.bray.Potable_Overlap

ADONIS2.R2.All.PotableBothPOCPOU.bray.Potable_Overlap<-ADONIS2.All.PotableBothPOCPOU.bray.Potable_Overlap$R2[1]
ADONIS2.R2.All.PotableBothPOCPOU.bray.Potable_Overlap  

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.PotableBothPOCPOU.bray["Potable_Overlap"]<-ANOSIM.Pval.All.PotableBothPOCPOU.bray.Potable_Overlap
ADONIS2.GroupedP.All.PotableBothPOCPOU.bray["Potable_Overlap"]<-ADONIS2.Pval.All.PotableBothPOCPOU.bray.Potable_Overlap

Betadis.GroupedP.All.PotableBothPOCPOU.bray["Potable_Overlap"]<-Betadis.Pval.All.PotableBothPOCPOU.bray.Potable_Overlap

ANOSIM.GroupedR.All.PotableBothPOCPOU.bray["Potable_Overlap"]<-ANOSIM.R.All.PotableBothPOCPOU.bray.Potable_Overlap
ADONIS2.GroupedR.All.PotableBothPOCPOU.bray["Potable_Overlap"]<-ADONIS2.R2.All.PotableBothPOCPOU.bray.Potable_Overlap

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.PotableBothPOCPOU.meta$Potable_Overlap)

ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Potable_Overlap<-pairwise.adonis2(All.PotableBothPOCPOU.dist.bray ~ Potable_Overlap, 
                                                                               data = All.PotableBothPOCPOU.meta, 
                                                                               permutations = permutations)
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Potable_Overlap
#Potable_Overlap
##Potable_Overlap




#################
#################
#################
#################
#Includes_Ozone 
ANOSIM.All.PotableBothPOCPOU.bray.Includes_Ozone<-anosim(All.PotableBothPOCPOU.dist.bray,
                                                          All.PotableBothPOCPOU.meta$Includes_Ozone, 
                                                          permutations = permutations)

ADONIS2.All.PotableBothPOCPOU.bray.Includes_Ozone<-adonis2(All.PotableBothPOCPOU.dist.bray ~ Includes_Ozone, 
                                                            data = All.PotableBothPOCPOU.meta, 
                                                            permutations = permutations)

Betadis.All.PotableBothPOCPOU.bray.Includes_Ozone<-anova(betadisper(All.PotableBothPOCPOU.dist.bray,
                                                                     All.PotableBothPOCPOU.meta$Includes_Ozone))



ANOSIM.Pval.All.PotableBothPOCPOU.bray.Includes_Ozone<-ANOSIM.All.PotableBothPOCPOU.bray.Includes_Ozone$signif
ANOSIM.Pval.All.PotableBothPOCPOU.bray.Includes_Ozone

ADONIS2.Pval.All.PotableBothPOCPOU.bray.Includes_Ozone<-ADONIS2.All.PotableBothPOCPOU.bray.Includes_Ozone$`Pr(>F)`[1]
ADONIS2.Pval.All.PotableBothPOCPOU.bray.Includes_Ozone  

Betadis.Pval.All.PotableBothPOCPOU.bray.Includes_Ozone<-Betadis.All.PotableBothPOCPOU.bray.Includes_Ozone$`Pr(>F)`[1]
Betadis.Pval.All.PotableBothPOCPOU.bray.Includes_Ozone

ANOSIM.R.All.PotableBothPOCPOU.bray.Includes_Ozone<-ANOSIM.All.PotableBothPOCPOU.bray.Includes_Ozone$statistic
ANOSIM.R.All.PotableBothPOCPOU.bray.Includes_Ozone

ADONIS2.R2.All.PotableBothPOCPOU.bray.Includes_Ozone<-ADONIS2.All.PotableBothPOCPOU.bray.Includes_Ozone$R2[1]
ADONIS2.R2.All.PotableBothPOCPOU.bray.Includes_Ozone  

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.PotableBothPOCPOU.bray["Includes_Ozone"]<-ANOSIM.Pval.All.PotableBothPOCPOU.bray.Includes_Ozone
ADONIS2.GroupedP.All.PotableBothPOCPOU.bray["Includes_Ozone"]<-ADONIS2.Pval.All.PotableBothPOCPOU.bray.Includes_Ozone

Betadis.GroupedP.All.PotableBothPOCPOU.bray["Includes_Ozone"]<-Betadis.Pval.All.PotableBothPOCPOU.bray.Includes_Ozone

ANOSIM.GroupedR.All.PotableBothPOCPOU.bray["Includes_Ozone"]<-ANOSIM.R.All.PotableBothPOCPOU.bray.Includes_Ozone
ADONIS2.GroupedR.All.PotableBothPOCPOU.bray["Includes_Ozone"]<-ADONIS2.R2.All.PotableBothPOCPOU.bray.Includes_Ozone

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.PotableBothPOCPOU.meta$Includes_Ozone)

ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Includes_Ozone<-pairwise.adonis2(All.PotableBothPOCPOU.dist.bray ~ Includes_Ozone, 
                                                                              data = All.PotableBothPOCPOU.meta, 
                                                                              permutations = permutations)
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Includes_Ozone
#Includes_Ozone
##Includes_Ozone 




#################
#################
#################
#################
#Includes_RO_UF 
ANOSIM.All.PotableBothPOCPOU.bray.Includes_RO_UF<-anosim(All.PotableBothPOCPOU.dist.bray,
                                                          All.PotableBothPOCPOU.meta$Includes_RO_UF, 
                                                          permutations = permutations)

ADONIS2.All.PotableBothPOCPOU.bray.Includes_RO_UF<-adonis2(All.PotableBothPOCPOU.dist.bray ~ Includes_RO_UF, 
                                                            data = All.PotableBothPOCPOU.meta, 
                                                            permutations = permutations)

Betadis.All.PotableBothPOCPOU.bray.Includes_RO_UF<-anova(betadisper(All.PotableBothPOCPOU.dist.bray,
                                                                     All.PotableBothPOCPOU.meta$Includes_RO_UF))



ANOSIM.Pval.All.PotableBothPOCPOU.bray.Includes_RO_UF<-ANOSIM.All.PotableBothPOCPOU.bray.Includes_RO_UF$signif
ANOSIM.Pval.All.PotableBothPOCPOU.bray.Includes_RO_UF

ADONIS2.Pval.All.PotableBothPOCPOU.bray.Includes_RO_UF<-ADONIS2.All.PotableBothPOCPOU.bray.Includes_RO_UF$`Pr(>F)`[1]
ADONIS2.Pval.All.PotableBothPOCPOU.bray.Includes_RO_UF  

Betadis.Pval.All.PotableBothPOCPOU.bray.Includes_RO_UF<-Betadis.All.PotableBothPOCPOU.bray.Includes_RO_UF$`Pr(>F)`[1]
Betadis.Pval.All.PotableBothPOCPOU.bray.Includes_RO_UF

ANOSIM.R.All.PotableBothPOCPOU.bray.Includes_RO_UF<-ANOSIM.All.PotableBothPOCPOU.bray.Includes_RO_UF$statistic
ANOSIM.R.All.PotableBothPOCPOU.bray.Includes_RO_UF

ADONIS2.R2.All.PotableBothPOCPOU.bray.Includes_RO_UF<-ADONIS2.All.PotableBothPOCPOU.bray.Includes_RO_UF$R2[1]
ADONIS2.R2.All.PotableBothPOCPOU.bray.Includes_RO_UF  

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.PotableBothPOCPOU.bray["Includes_RO_UF"]<-ANOSIM.Pval.All.PotableBothPOCPOU.bray.Includes_RO_UF
ADONIS2.GroupedP.All.PotableBothPOCPOU.bray["Includes_RO_UF"]<-ADONIS2.Pval.All.PotableBothPOCPOU.bray.Includes_RO_UF

Betadis.GroupedP.All.PotableBothPOCPOU.bray["Includes_RO_UF"]<-Betadis.Pval.All.PotableBothPOCPOU.bray.Includes_RO_UF

ANOSIM.GroupedR.All.PotableBothPOCPOU.bray["Includes_RO_UF"]<-ANOSIM.R.All.PotableBothPOCPOU.bray.Includes_RO_UF
ADONIS2.GroupedR.All.PotableBothPOCPOU.bray["Includes_RO_UF"]<-ADONIS2.R2.All.PotableBothPOCPOU.bray.Includes_RO_UF

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.PotableBothPOCPOU.meta$Includes_RO_UF)

ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Includes_RO_UF<-pairwise.adonis2(All.PotableBothPOCPOU.dist.bray ~ Includes_RO_UF, 
                                                                              data = All.PotableBothPOCPOU.meta, 
                                                                              permutations = permutations)
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Includes_RO_UF
#Includes_RO_UF
##Includes_RO_UF 




#################
#################
#################
#################
#Includes_UV 
ANOSIM.All.PotableBothPOCPOU.bray.Includes_UV<-anosim(All.PotableBothPOCPOU.dist.bray,
                                                       All.PotableBothPOCPOU.meta$Includes_UV, 
                                                       permutations = permutations)

ADONIS2.All.PotableBothPOCPOU.bray.Includes_UV<-adonis2(All.PotableBothPOCPOU.dist.bray ~ Includes_UV, 
                                                         data = All.PotableBothPOCPOU.meta, 
                                                         permutations = permutations)

Betadis.All.PotableBothPOCPOU.bray.Includes_UV<-anova(betadisper(All.PotableBothPOCPOU.dist.bray,
                                                                  All.PotableBothPOCPOU.meta$Includes_UV))



ANOSIM.Pval.All.PotableBothPOCPOU.bray.Includes_UV<-ANOSIM.All.PotableBothPOCPOU.bray.Includes_UV$signif
ANOSIM.Pval.All.PotableBothPOCPOU.bray.Includes_UV

ADONIS2.Pval.All.PotableBothPOCPOU.bray.Includes_UV<-ADONIS2.All.PotableBothPOCPOU.bray.Includes_UV$`Pr(>F)`[1]
ADONIS2.Pval.All.PotableBothPOCPOU.bray.Includes_UV  

Betadis.Pval.All.PotableBothPOCPOU.bray.Includes_UV<-Betadis.All.PotableBothPOCPOU.bray.Includes_UV$`Pr(>F)`[1]
Betadis.Pval.All.PotableBothPOCPOU.bray.Includes_UV

ANOSIM.R.All.PotableBothPOCPOU.bray.Includes_UV<-ANOSIM.All.PotableBothPOCPOU.bray.Includes_UV$statistic
ANOSIM.R.All.PotableBothPOCPOU.bray.Includes_UV

ADONIS2.R2.All.PotableBothPOCPOU.bray.Includes_UV<-ADONIS2.All.PotableBothPOCPOU.bray.Includes_UV$R2[1]
ADONIS2.R2.All.PotableBothPOCPOU.bray.Includes_UV  

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.PotableBothPOCPOU.bray["Includes_UV"]<-ANOSIM.Pval.All.PotableBothPOCPOU.bray.Includes_UV
ADONIS2.GroupedP.All.PotableBothPOCPOU.bray["Includes_UV"]<-ADONIS2.Pval.All.PotableBothPOCPOU.bray.Includes_UV

Betadis.GroupedP.All.PotableBothPOCPOU.bray["Includes_UV"]<-Betadis.Pval.All.PotableBothPOCPOU.bray.Includes_UV

ANOSIM.GroupedR.All.PotableBothPOCPOU.bray["Includes_UV"]<-ANOSIM.R.All.PotableBothPOCPOU.bray.Includes_UV
ADONIS2.GroupedR.All.PotableBothPOCPOU.bray["Includes_UV"]<-ADONIS2.R2.All.PotableBothPOCPOU.bray.Includes_UV

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.PotableBothPOCPOU.meta$Includes_UV)

ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Includes_UV<-pairwise.adonis2(All.PotableBothPOCPOU.dist.bray ~ Includes_UV, 
                                                                           data = All.PotableBothPOCPOU.meta, 
                                                                           permutations = permutations)
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Includes_UV
#Includes_UV
##Includes_UV




# #################
# #################
# #################
# #################
# #Type_DS  
# ANOSIM.All.PotableBothPOCPOU.bray.Type_DS <-anosim(All.PotableBothPOCPOU.dist.bray,
#                                                     All.PotableBothPOCPOU.meta$Type_DS , 
#                                                     permutations = permutations)
# 
# ADONIS2.All.PotableBothPOCPOU.bray.Type_DS <-adonis2(All.PotableBothPOCPOU.dist.bray ~ Type_DS , 
#                                                       data = All.PotableBothPOCPOU.meta, 
#                                                       permutations = permutations)
# 
# Betadis.All.PotableBothPOCPOU.bray.Type_DS <-anova(betadisper(All.PotableBothPOCPOU.dist.bray,
#                                                                All.PotableBothPOCPOU.meta$Type_DS ))
# 
# 
# 
# ANOSIM.Pval.All.PotableBothPOCPOU.bray.Type_DS <-ANOSIM.All.PotableBothPOCPOU.bray.Type_DS $signif
# ANOSIM.Pval.All.PotableBothPOCPOU.bray.Type_DS 
# 
# ADONIS2.Pval.All.PotableBothPOCPOU.bray.Type_DS <-ADONIS2.All.PotableBothPOCPOU.bray.Type_DS $`Pr(>F)`[1]
# ADONIS2.Pval.All.PotableBothPOCPOU.bray.Type_DS   
# 
# Betadis.Pval.All.PotableBothPOCPOU.bray.Type_DS <-Betadis.All.PotableBothPOCPOU.bray.Type_DS $`Pr(>F)`[1]
# Betadis.Pval.All.PotableBothPOCPOU.bray.Type_DS 
# 
# ANOSIM.R.All.PotableBothPOCPOU.bray.Type_DS <-ANOSIM.All.PotableBothPOCPOU.bray.Type_DS $statistic
# ANOSIM.R.All.PotableBothPOCPOU.bray.Type_DS 
# 
# ADONIS2.R2.All.PotableBothPOCPOU.bray.Type_DS <-ADONIS2.All.PotableBothPOCPOU.bray.Type_DS $R2[1]
# ADONIS2.R2.All.PotableBothPOCPOU.bray.Type_DS   
# 
# #Pull Together Pvalues to adjust 
# ANOSIM.GroupedP.All.PotableBothPOCPOU.bray["Type_DS "]<-ANOSIM.Pval.All.PotableBothPOCPOU.bray.Type_DS 
# ADONIS2.GroupedP.All.PotableBothPOCPOU.bray["Type_DS "]<-ADONIS2.Pval.All.PotableBothPOCPOU.bray.Type_DS 
# 
# Betadis.GroupedP.All.PotableBothPOCPOU.bray["Type_DS "]<-Betadis.Pval.All.PotableBothPOCPOU.bray.Type_DS 
# 
# ANOSIM.GroupedR.All.PotableBothPOCPOU.bray["Type_DS "]<-ANOSIM.R.All.PotableBothPOCPOU.bray.Type_DS 
# ADONIS2.GroupedR.All.PotableBothPOCPOU.bray["Type_DS "]<-ADONIS2.R2.All.PotableBothPOCPOU.bray.Type_DS 
# 
# #Pairwise if needed, check if interesting or not. Remove if not. 
# count(All.PotableBothPOCPOU.meta$Type_DS )
# 
# ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Type_DS <-pairwise.adonis2(All.PotableBothPOCPOU.dist.bray ~ Type_DS , 
#                                                                         data = All.PotableBothPOCPOU.meta, 
#                                                                         permutations = permutations)
# ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Type_DS 
# #Type_DS 
# ##Type_DS 


#################
#################
#################
#################
#Reuse.Reclaimed_Treatment  
ANOSIM.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment <-anosim(All.PotableBothPOCPOU.dist.bray,
                                                                      All.PotableBothPOCPOU.meta$Reuse.Reclaimed_Treatment , 
                                                                      permutations = permutations)

ADONIS2.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment <-adonis2(All.PotableBothPOCPOU.dist.bray ~ Reuse.Reclaimed_Treatment , 
                                                                        data = All.PotableBothPOCPOU.meta, 
                                                                        permutations = permutations)

Betadis.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment <-anova(betadisper(All.PotableBothPOCPOU.dist.bray,
                                                                                 All.PotableBothPOCPOU.meta$Reuse.Reclaimed_Treatment ))



ANOSIM.Pval.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment <-ANOSIM.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment $signif
ANOSIM.Pval.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment 

ADONIS2.Pval.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment <-ADONIS2.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment $`Pr(>F)`[1]
ADONIS2.Pval.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment   

Betadis.Pval.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment <-Betadis.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment $`Pr(>F)`[1]
Betadis.Pval.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment 

ANOSIM.R.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment <-ANOSIM.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment $statistic
ANOSIM.R.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment 

ADONIS2.R2.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment <-ADONIS2.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment $R2[1]
ADONIS2.R2.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment   

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.PotableBothPOCPOU.bray["Reuse.Reclaimed_Treatment "]<-ANOSIM.Pval.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment 
ADONIS2.GroupedP.All.PotableBothPOCPOU.bray["Reuse.Reclaimed_Treatment "]<-ADONIS2.Pval.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment 

Betadis.GroupedP.All.PotableBothPOCPOU.bray["Reuse.Reclaimed_Treatment "]<-Betadis.Pval.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment 

ANOSIM.GroupedR.All.PotableBothPOCPOU.bray["Reuse.Reclaimed_Treatment "]<-ANOSIM.R.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment 
ADONIS2.GroupedR.All.PotableBothPOCPOU.bray["Reuse.Reclaimed_Treatment "]<-ADONIS2.R2.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment 

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.PotableBothPOCPOU.meta$Reuse.Reclaimed_Treatment )

ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment <-pairwise.adonis2(All.PotableBothPOCPOU.dist.bray ~ Reuse.Reclaimed_Treatment , 
                                                                                          data = All.PotableBothPOCPOU.meta, 
                                                                                          permutations = permutations)
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment 
#Reuse.Reclaimed_Treatment 
##Reuse.Reclaimed_Disinfection




#################
#################
#################
#################
#Reuse.Reclaimed_Treatment  
ANOSIM.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment <-anosim(All.PotableBothPOCPOU.dist.bray,
                                                                      All.PotableBothPOCPOU.meta$Reuse.Reclaimed_Treatment , 
                                                                      permutations = permutations)

ADONIS2.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment <-adonis2(All.PotableBothPOCPOU.dist.bray ~ Reuse.Reclaimed_Treatment , 
                                                                        data = All.PotableBothPOCPOU.meta, 
                                                                        permutations = permutations)

Betadis.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment <-anova(betadisper(All.PotableBothPOCPOU.dist.bray,
                                                                                 All.PotableBothPOCPOU.meta$Reuse.Reclaimed_Treatment ))



ANOSIM.Pval.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment <-ANOSIM.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment $signif
ANOSIM.Pval.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment 

ADONIS2.Pval.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment <-ADONIS2.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment $`Pr(>F)`[1]
ADONIS2.Pval.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment   

Betadis.Pval.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment <-Betadis.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment $`Pr(>F)`[1]
Betadis.Pval.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment 

ANOSIM.R.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment <-ANOSIM.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment $statistic
ANOSIM.R.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment 

ADONIS2.R2.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment <-ADONIS2.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment $R2[1]
ADONIS2.R2.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment   

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.PotableBothPOCPOU.bray["Reuse.Reclaimed_Treatment "]<-ANOSIM.Pval.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment 
ADONIS2.GroupedP.All.PotableBothPOCPOU.bray["Reuse.Reclaimed_Treatment "]<-ADONIS2.Pval.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment 

Betadis.GroupedP.All.PotableBothPOCPOU.bray["Reuse.Reclaimed_Treatment "]<-Betadis.Pval.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment 

ANOSIM.GroupedR.All.PotableBothPOCPOU.bray["Reuse.Reclaimed_Treatment "]<-ANOSIM.R.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment 
ADONIS2.GroupedR.All.PotableBothPOCPOU.bray["Reuse.Reclaimed_Treatment "]<-ADONIS2.R2.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment 

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.PotableBothPOCPOU.meta$Reuse.Reclaimed_Treatment )

ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment <-pairwise.adonis2(All.PotableBothPOCPOU.dist.bray ~ Reuse.Reclaimed_Treatment , 
                                                                                          data = All.PotableBothPOCPOU.meta, 
                                                                                          permutations = permutations)
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment 
#Reuse.Reclaimed_Treatment 
##Reuse.Reclaimed_Treatment 




#################
#################
#################
#################
#Reuse.Reclaimed_Blending 
ANOSIM.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Blending<-anosim(All.PotableBothPOCPOU.dist.bray,
                                                                    All.PotableBothPOCPOU.meta$Reuse.Reclaimed_Blending, 
                                                                    permutations = permutations)

ADONIS2.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Blending<-adonis2(All.PotableBothPOCPOU.dist.bray ~ Reuse.Reclaimed_Blending, 
                                                                      data = All.PotableBothPOCPOU.meta, 
                                                                      permutations = permutations)

Betadis.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Blending<-anova(betadisper(All.PotableBothPOCPOU.dist.bray,
                                                                               All.PotableBothPOCPOU.meta$Reuse.Reclaimed_Blending))



ANOSIM.Pval.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Blending<-ANOSIM.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Blending$signif
ANOSIM.Pval.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Blending

ADONIS2.Pval.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Blending<-ADONIS2.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Blending$`Pr(>F)`[1]
ADONIS2.Pval.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Blending  

Betadis.Pval.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Blending<-Betadis.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Blending$`Pr(>F)`[1]
Betadis.Pval.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Blending

ANOSIM.R.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Blending<-ANOSIM.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Blending$statistic
ANOSIM.R.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Blending

ADONIS2.R2.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Blending<-ADONIS2.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Blending$R2[1]
ADONIS2.R2.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Blending  

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.PotableBothPOCPOU.bray["Reuse.Reclaimed_Blending"]<-ANOSIM.Pval.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Blending
ADONIS2.GroupedP.All.PotableBothPOCPOU.bray["Reuse.Reclaimed_Blending"]<-ADONIS2.Pval.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Blending

Betadis.GroupedP.All.PotableBothPOCPOU.bray["Reuse.Reclaimed_Blending"]<-Betadis.Pval.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Blending

ANOSIM.GroupedR.All.PotableBothPOCPOU.bray["Reuse.Reclaimed_Blending"]<-ANOSIM.R.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Blending
ADONIS2.GroupedR.All.PotableBothPOCPOU.bray["Reuse.Reclaimed_Blending"]<-ADONIS2.R2.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Blending

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.PotableBothPOCPOU.meta$Reuse.Reclaimed_Blending)

ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Blending<-pairwise.adonis2(All.PotableBothPOCPOU.dist.bray ~ Reuse.Reclaimed_Blending, 
                                                                                        data = All.PotableBothPOCPOU.meta, 
                                                                                        permutations = permutations)
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Blending
#Reuse.Reclaimed_Blending
##Reuse.Reclaimed_Blending




#################
#################
#################
#################
#PostBlend_Disinfection 
ANOSIM.All.PotableBothPOCPOU.bray.PostBlend_Disinfection<-anosim(All.PotableBothPOCPOU.dist.bray,
                                                                  All.PotableBothPOCPOU.meta$PostBlend_Disinfection, 
                                                                  permutations = permutations)

ADONIS2.All.PotableBothPOCPOU.bray.PostBlend_Disinfection<-adonis2(All.PotableBothPOCPOU.dist.bray ~ PostBlend_Disinfection, 
                                                                    data = All.PotableBothPOCPOU.meta, 
                                                                    permutations = permutations)

Betadis.All.PotableBothPOCPOU.bray.PostBlend_Disinfection<-anova(betadisper(All.PotableBothPOCPOU.dist.bray,
                                                                             All.PotableBothPOCPOU.meta$PostBlend_Disinfection))



ANOSIM.Pval.All.PotableBothPOCPOU.bray.PostBlend_Disinfection<-ANOSIM.All.PotableBothPOCPOU.bray.PostBlend_Disinfection$signif
ANOSIM.Pval.All.PotableBothPOCPOU.bray.PostBlend_Disinfection

ADONIS2.Pval.All.PotableBothPOCPOU.bray.PostBlend_Disinfection<-ADONIS2.All.PotableBothPOCPOU.bray.PostBlend_Disinfection$`Pr(>F)`[1]
ADONIS2.Pval.All.PotableBothPOCPOU.bray.PostBlend_Disinfection  

Betadis.Pval.All.PotableBothPOCPOU.bray.PostBlend_Disinfection<-Betadis.All.PotableBothPOCPOU.bray.PostBlend_Disinfection$`Pr(>F)`[1]
Betadis.Pval.All.PotableBothPOCPOU.bray.PostBlend_Disinfection

ANOSIM.R.All.PotableBothPOCPOU.bray.PostBlend_Disinfection<-ANOSIM.All.PotableBothPOCPOU.bray.PostBlend_Disinfection$statistic
ANOSIM.R.All.PotableBothPOCPOU.bray.PostBlend_Disinfection

ADONIS2.R2.All.PotableBothPOCPOU.bray.PostBlend_Disinfection<-ADONIS2.All.PotableBothPOCPOU.bray.PostBlend_Disinfection$R2[1]
ADONIS2.R2.All.PotableBothPOCPOU.bray.PostBlend_Disinfection  

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.PotableBothPOCPOU.bray["PostBlend_Disinfection"]<-ANOSIM.Pval.All.PotableBothPOCPOU.bray.PostBlend_Disinfection
ADONIS2.GroupedP.All.PotableBothPOCPOU.bray["PostBlend_Disinfection"]<-ADONIS2.Pval.All.PotableBothPOCPOU.bray.PostBlend_Disinfection

Betadis.GroupedP.All.PotableBothPOCPOU.bray["PostBlend_Disinfection"]<-Betadis.Pval.All.PotableBothPOCPOU.bray.PostBlend_Disinfection

ANOSIM.GroupedR.All.PotableBothPOCPOU.bray["PostBlend_Disinfection"]<-ANOSIM.R.All.PotableBothPOCPOU.bray.PostBlend_Disinfection
ADONIS2.GroupedR.All.PotableBothPOCPOU.bray["PostBlend_Disinfection"]<-ADONIS2.R2.All.PotableBothPOCPOU.bray.PostBlend_Disinfection

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.PotableBothPOCPOU.meta$PostBlend_Disinfection)

ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.PostBlend_Disinfection<-pairwise.adonis2(All.PotableBothPOCPOU.dist.bray ~ PostBlend_Disinfection, 
                                                                                      data = All.PotableBothPOCPOU.meta, 
                                                                                      permutations = permutations)
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.PostBlend_Disinfection
#PostBlend_Disinfection
##PostBlend_Disinfection




#################
#################
#################
#################
#PostBlend_Treatment 
ANOSIM.All.PotableBothPOCPOU.bray.PostBlend_Treatment<-anosim(All.PotableBothPOCPOU.dist.bray,
                                                               All.PotableBothPOCPOU.meta$PostBlend_Treatment, 
                                                               permutations = permutations)

ADONIS2.All.PotableBothPOCPOU.bray.PostBlend_Treatment<-adonis2(All.PotableBothPOCPOU.dist.bray ~ PostBlend_Treatment, 
                                                                 data = All.PotableBothPOCPOU.meta, 
                                                                 permutations = permutations)

Betadis.All.PotableBothPOCPOU.bray.PostBlend_Treatment<-anova(betadisper(All.PotableBothPOCPOU.dist.bray,
                                                                          All.PotableBothPOCPOU.meta$PostBlend_Treatment))



ANOSIM.Pval.All.PotableBothPOCPOU.bray.PostBlend_Treatment<-ANOSIM.All.PotableBothPOCPOU.bray.PostBlend_Treatment$signif
ANOSIM.Pval.All.PotableBothPOCPOU.bray.PostBlend_Treatment

ADONIS2.Pval.All.PotableBothPOCPOU.bray.PostBlend_Treatment<-ADONIS2.All.PotableBothPOCPOU.bray.PostBlend_Treatment$`Pr(>F)`[1]
ADONIS2.Pval.All.PotableBothPOCPOU.bray.PostBlend_Treatment  

Betadis.Pval.All.PotableBothPOCPOU.bray.PostBlend_Treatment<-Betadis.All.PotableBothPOCPOU.bray.PostBlend_Treatment$`Pr(>F)`[1]
Betadis.Pval.All.PotableBothPOCPOU.bray.PostBlend_Treatment

ANOSIM.R.All.PotableBothPOCPOU.bray.PostBlend_Treatment<-ANOSIM.All.PotableBothPOCPOU.bray.PostBlend_Treatment$statistic
ANOSIM.R.All.PotableBothPOCPOU.bray.PostBlend_Treatment

ADONIS2.R2.All.PotableBothPOCPOU.bray.PostBlend_Treatment<-ADONIS2.All.PotableBothPOCPOU.bray.PostBlend_Treatment$R2[1]
ADONIS2.R2.All.PotableBothPOCPOU.bray.PostBlend_Treatment  

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.PotableBothPOCPOU.bray["PostBlend_Treatment"]<-ANOSIM.Pval.All.PotableBothPOCPOU.bray.PostBlend_Treatment
ADONIS2.GroupedP.All.PotableBothPOCPOU.bray["PostBlend_Treatment"]<-ADONIS2.Pval.All.PotableBothPOCPOU.bray.PostBlend_Treatment

Betadis.GroupedP.All.PotableBothPOCPOU.bray["PostBlend_Treatment"]<-Betadis.Pval.All.PotableBothPOCPOU.bray.PostBlend_Treatment

ANOSIM.GroupedR.All.PotableBothPOCPOU.bray["PostBlend_Treatment"]<-ANOSIM.R.All.PotableBothPOCPOU.bray.PostBlend_Treatment
ADONIS2.GroupedR.All.PotableBothPOCPOU.bray["PostBlend_Treatment"]<-ADONIS2.R2.All.PotableBothPOCPOU.bray.PostBlend_Treatment

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.PotableBothPOCPOU.meta$PostBlend_Treatment)

ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.PostBlend_Treatment<-pairwise.adonis2(All.PotableBothPOCPOU.dist.bray ~ PostBlend_Treatment, 
                                                                                   data = All.PotableBothPOCPOU.meta, 
                                                                                   permutations = permutations)
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.PostBlend_Treatment
#PostBlend_Treatment
##PostBlend_Treatment





#Adjust PValue, Check Significance
ANOSIM.GroupedP.Adjusted.All.PotableBothPOCPOU.bray<-p.adjust(ANOSIM.GroupedP.All.PotableBothPOCPOU.bray,method="BY")
p.adjust(ANOSIM.GroupedP.Adjusted.All.PotableBothPOCPOU.bray,method="BY")<.05

ADONIS2.GroupedP.Adjusted.All.PotableBothPOCPOU.bray<-p.adjust(ADONIS2.GroupedP.All.PotableBothPOCPOU.bray,method="BY")
p.adjust(ADONIS2.GroupedP.Adjusted.All.PotableBothPOCPOU.bray,method="BY")<.05

Betadis.GroupedP.All.PotableBothPOCPOU.bray<0.05
Betadis.GroupedP.Adjusted.All.PotableBothPOCPOU.bray<-p.adjust(Betadis.GroupedP.All.PotableBothPOCPOU.bray,method="BY")
p.adjust(Betadis.GroupedP.Adjusted.All.PotableBothPOCPOU.bray,method="BY")<.05


#Check Value of Pairwise Sub comparisons 
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Class1_Pot_NonPot
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Classification_1
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Classification_1_mod
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Classification_2
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Climate
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Climate_Region
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Final_Disinfection_Residual
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Primers
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Region
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Sample_Local
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Trimmed
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Extraction_Type
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Potable_Source_Water
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Potable_Disinfection_Primary
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Potable_Treatment_Type
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Disinfection_Total
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Treatment_Total
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Treatment_Classification_Number
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Treatment_Classification_Cat
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Mixed_Source_Water
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Potable_Overlap
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Includes_Ozone
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Includes_RO_UF
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Includes_UV
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Type_DS
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Disinfection
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Treatment
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.Reuse.Reclaimed_Blending
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.PostBlend_Disinfection
ADONIS2.Pairwise.All.PotableBothPOCPOU.bray.PostBlend_Treatment

#Permanova on Multiple Variables that are significant and not heterogeneous (not sig from Beta, and sig from anosim)
ADONIS2.All.PotableBothPOCPOU.bray.multiple_comparisons<-adonis2(All.PotableBothPOCPOU.dist.bray ~ 
                                                                    Classification_1+Climate, 
                                                                  data = All.PotableBothPOCPOU.meta, 
                                                                  permutations = permutations)
ADONIS2.All.PotableBothPOCPOU.bray.multiple_comparisons


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
All.PotableBothPOCPOU.Stats<-data.frame(
  "ANOSIM.P.All.PotableBothPOCPOU.bray"= ANOSIM.GroupedP.All.PotableBothPOCPOU.bray,
  "ANOSIM.P.Adj.All.PotableBothPOCPOU.bray"= ANOSIM.GroupedP.Adjusted.All.PotableBothPOCPOU.bray,
  "ADONIS2.P.All.PotableBothPOCPOU.bray"= ADONIS2.GroupedP.All.PotableBothPOCPOU.bray,
  "ADONIS2.P.Adj.All.PotableBothPOCPOU.bray"= ADONIS2.GroupedP.Adjusted.All.PotableBothPOCPOU.bray,
  "Betadis.P.All.PotableBothPOCPOU.bray"= Betadis.GroupedP.All.PotableBothPOCPOU.bray,
  "Betadis.P.Adj.All.PotableBothPOCPOU.bray"= Betadis.GroupedP.Adjusted.All.PotableBothPOCPOU.bray,
  "ANOSIM.R.All.PotableBothPOCPOU.bray"=  ANOSIM.GroupedR.All.PotableBothPOCPOU.bray,
  "ADONIS2.R.All.PotableBothPOCPOU.bray"= ADONIS2.GroupedR.All.PotableBothPOCPOU.bray)


write.xlsx(All.PotableBothPOCPOU.Stats, file = "C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Stats/Beta/All.PotableBothPOCPOU.List.xlsx", rowNames = TRUE)

