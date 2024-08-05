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
library(pairwiseAdonis)
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

########Sample_Local, BW removed 

#Read Meta Data 
Metadata_POCPOU<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_DS_Special', rowNames = TRUE)
#Metadata_POC<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_Final', rowNames = TRUE)
Metadata_POC<-Metadata_POCPOU[Metadata_POCPOU$Sample_Local == "POC",]
#Metadata_POU<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'DS_Final', rowNames = TRUE)
Metadata_POU<-Metadata_POCPOU[Metadata_POCPOU$Sample_Local == "POU",]

#Order Meta Data 
Metadata_POCPOU$Classification_1 = factor(Metadata_POCPOU$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))
Metadata_POC$Classification_1 = factor(Metadata_POC$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))
Metadata_POU$Classification_1 = factor(Metadata_POU$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))

######All Samples, POCPOU
Meta.All.POCPOU<-Metadata_POCPOU
Meta.All.POCPOU<-Meta.All.POCPOU[Meta.All.POCPOU$Description != 'ConnerPrimer58',] #removed because of low read count
Meta.All.POCPOU<-Meta.All.POCPOU[Meta.All.POCPOU$Description != 'PR1WEEF',] #removed because of high read count 

Meta.All.POCPOU.Filter<-Meta.All.POCPOU[Meta.All.POCPOU$Classification_1 == 'Potable'| Meta.All.POCPOU$Classification_1 =='Non-potable Reuse'|Meta.All.POCPOU$Classification_1 =='Potable Reuse'|Meta.All.POCPOU$Classification_1 =='',]#can add Blank
Meta.All.POCPOU.Filter<-Meta.All.POCPOU.Filter[Meta.All.POCPOU.Filter$Sample_Local != 'Bottled Water',]
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
All.POCPOU.dist.bray<- phyloseq::distance(physeq.combined.All.POCPOU, method = "bray")
All.POCPOU.meta<- data.frame(sample_data(physeq.combined.All.POCPOU))

ANOSIM.GroupedP.All.POCPOU.bray<-numeric() 
ADONIS2.GroupedP.All.POCPOU.bray<-numeric()
Betadis.GroupedP.All.POCPOU.bray<-numeric()
ANOSIM.GroupedR.All.POCPOU.bray<-numeric() 
ADONIS2.GroupedR.All.POCPOU.bray<-numeric()

#################
#################
#################
#################
#Sample_Local 
ANOSIM.All.POCPOU.bray.Sample_Local<-anosim(All.POCPOU.dist.bray,
                                                 All.POCPOU.meta$Sample_Local, 
                                                 permutations = permutations)

ADONIS2.All.POCPOU.bray.Sample_Local<-adonis2(All.POCPOU.dist.bray ~ Sample_Local, 
                                                   data = All.POCPOU.meta, 
                                                   permutations = permutations)

Betadis.All.POCPOU.bray.Sample_Local<-anova(betadisper(All.POCPOU.dist.bray,
                                                            All.POCPOU.meta$Sample_Local))



ANOSIM.Pval.All.POCPOU.bray.Sample_Local<-ANOSIM.All.POCPOU.bray.Sample_Local$signif
ANOSIM.Pval.All.POCPOU.bray.Sample_Local

ADONIS2.Pval.All.POCPOU.bray.Sample_Local<-ADONIS2.All.POCPOU.bray.Sample_Local$`Pr(>F)`[1]
ADONIS2.Pval.All.POCPOU.bray.Sample_Local  

Betadis.Pval.All.POCPOU.bray.Sample_Local<-Betadis.All.POCPOU.bray.Sample_Local$`Pr(>F)`[1]
Betadis.Pval.All.POCPOU.bray.Sample_Local

ANOSIM.R.All.POCPOU.bray.Sample_Local<-ANOSIM.All.POCPOU.bray.Sample_Local$statistic
ANOSIM.R.All.POCPOU.bray.Sample_Local

ADONIS2.R2.All.POCPOU.bray.Sample_Local<-ADONIS2.All.POCPOU.bray.Sample_Local$R2[1]
ADONIS2.R2.All.POCPOU.bray.Sample_Local  

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.POCPOU.bray["Sample_Local"]<-ANOSIM.Pval.All.POCPOU.bray.Sample_Local
ADONIS2.GroupedP.All.POCPOU.bray["Sample_Local"]<-ADONIS2.Pval.All.POCPOU.bray.Sample_Local

Betadis.GroupedP.All.POCPOU.bray["Sample_Local"]<-Betadis.Pval.All.POCPOU.bray.Sample_Local

ANOSIM.GroupedR.All.POCPOU.bray["Sample_Local"]<-ANOSIM.R.All.POCPOU.bray.Sample_Local
ADONIS2.GroupedR.All.POCPOU.bray["Sample_Local"]<-ADONIS2.R2.All.POCPOU.bray.Sample_Local

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.POCPOU.meta$Sample_Local)

ADONIS2.Pairwise.All.POCPOU.bray.Sample_Local<-pairwise.adonis2(All.POCPOU.dist.bray ~ Sample_Local, 
                                                                     data = All.POCPOU.meta, 
                                                                     permutations = permutations)
ADONIS2.Pairwise.All.POCPOU.bray.Sample_Local
#Sample_Local














########Region 

#Read Meta Data 
Metadata_POCPOU<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_DS_Special', rowNames = TRUE)
#Metadata_POC<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_Final', rowNames = TRUE)
Metadata_POC<-Metadata_POCPOU[Metadata_POCPOU$Region == "POC",]
#Metadata_POU<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'DS_Final', rowNames = TRUE)
Metadata_POU<-Metadata_POCPOU[Metadata_POCPOU$Region == "POU",]

#Order Meta Data 
Metadata_POCPOU$Classification_1 = factor(Metadata_POCPOU$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))
Metadata_POC$Classification_1 = factor(Metadata_POC$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))
Metadata_POU$Classification_1 = factor(Metadata_POU$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))

######All Samples, POCPOU
Meta.All.POCPOU<-Metadata_POCPOU
Meta.All.POCPOU<-Meta.All.POCPOU[Meta.All.POCPOU$Description != 'ConnerPrimer58',] #removed because of low read count
Meta.All.POCPOU<-Meta.All.POCPOU[Meta.All.POCPOU$Description != 'PR1WEEF',] #removed because of high read count 

Meta.All.POCPOU.Filter<-Meta.All.POCPOU[Meta.All.POCPOU$Classification_1 == 'Potable'| Meta.All.POCPOU$Classification_1 =='Non-potable Reuse'|Meta.All.POCPOU$Classification_1 =='Potable Reuse'|Meta.All.POCPOU$Classification_1 =='',]#can add Blank
Meta.All.POCPOU.Filter<-Meta.All.POCPOU.Filter[Meta.All.POCPOU.Filter$Region == 'Southeast'|Meta.All.POCPOU.Filter$Region == 'West'|Meta.All.POCPOU.Filter$Region == 'Southwest',]
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
All.POCPOU.dist.bray<- phyloseq::distance(physeq.combined.All.POCPOU, method = "bray")
All.POCPOU.meta<- data.frame(sample_data(physeq.combined.All.POCPOU))

# ANOSIM.GroupedP.All.POCPOU.bray<-numeric() 
# ADONIS2.GroupedP.All.POCPOU.bray<-numeric()
# Betadis.GroupedP.All.POCPOU.bray<-numeric()
# ANOSIM.GroupedR.All.POCPOU.bray<-numeric() 
# ADONIS2.GroupedR.All.POCPOU.bray<-numeric()

#################
#################
#################
#################
#Region 
ANOSIM.All.POCPOU.bray.Region<-anosim(All.POCPOU.dist.bray,
                                            All.POCPOU.meta$Region, 
                                            permutations = permutations)

ADONIS2.All.POCPOU.bray.Region<-adonis2(All.POCPOU.dist.bray ~ Region, 
                                              data = All.POCPOU.meta, 
                                              permutations = permutations)

Betadis.All.POCPOU.bray.Region<-anova(betadisper(All.POCPOU.dist.bray,
                                                       All.POCPOU.meta$Region))



ANOSIM.Pval.All.POCPOU.bray.Region<-ANOSIM.All.POCPOU.bray.Region$signif
ANOSIM.Pval.All.POCPOU.bray.Region

ADONIS2.Pval.All.POCPOU.bray.Region<-ADONIS2.All.POCPOU.bray.Region$`Pr(>F)`[1]
ADONIS2.Pval.All.POCPOU.bray.Region  

Betadis.Pval.All.POCPOU.bray.Region<-Betadis.All.POCPOU.bray.Region$`Pr(>F)`[1]
Betadis.Pval.All.POCPOU.bray.Region

ANOSIM.R.All.POCPOU.bray.Region<-ANOSIM.All.POCPOU.bray.Region$statistic
ANOSIM.R.All.POCPOU.bray.Region

ADONIS2.R2.All.POCPOU.bray.Region<-ADONIS2.All.POCPOU.bray.Region$R2[1]
ADONIS2.R2.All.POCPOU.bray.Region  

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.POCPOU.bray["Region"]<-ANOSIM.Pval.All.POCPOU.bray.Region
ADONIS2.GroupedP.All.POCPOU.bray["Region"]<-ADONIS2.Pval.All.POCPOU.bray.Region

Betadis.GroupedP.All.POCPOU.bray["Region"]<-Betadis.Pval.All.POCPOU.bray.Region

ANOSIM.GroupedR.All.POCPOU.bray["Region"]<-ANOSIM.R.All.POCPOU.bray.Region
ADONIS2.GroupedR.All.POCPOU.bray["Region"]<-ADONIS2.R2.All.POCPOU.bray.Region

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.POCPOU.meta$Region)

ADONIS2.Pairwise.All.POCPOU.bray.Region<-pairwise.adonis2(All.POCPOU.dist.bray ~ Region, 
                                                                data = All.POCPOU.meta, 
                                                                permutations = permutations)
ADONIS2.Pairwise.All.POCPOU.bray.Region
#Region










########Climate 
#Read Meta Data 
Metadata_POCPOU<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_DS_Special', rowNames = TRUE)
#Metadata_POC<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_Final', rowNames = TRUE)
Metadata_POC<-Metadata_POCPOU[Metadata_POCPOU$Climate == "POC",]
#Metadata_POU<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'DS_Final', rowNames = TRUE)
Metadata_POU<-Metadata_POCPOU[Metadata_POCPOU$Climate == "POU",]

#Order Meta Data 
Metadata_POCPOU$Classification_1 = factor(Metadata_POCPOU$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))
Metadata_POC$Classification_1 = factor(Metadata_POC$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))
Metadata_POU$Classification_1 = factor(Metadata_POU$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))

######All Samples, POCPOU
Meta.All.POCPOU<-Metadata_POCPOU
Meta.All.POCPOU<-Meta.All.POCPOU[Meta.All.POCPOU$Description != 'ConnerPrimer58',] #removed because of low read count
Meta.All.POCPOU<-Meta.All.POCPOU[Meta.All.POCPOU$Description != 'PR1WEEF',] #removed because of high read count 

Meta.All.POCPOU.Filter<-Meta.All.POCPOU[Meta.All.POCPOU$Classification_1 == 'Potable'| Meta.All.POCPOU$Classification_1 =='Non-potable Reuse'|Meta.All.POCPOU$Classification_1 =='Potable Reuse'|Meta.All.POCPOU$Classification_1 =='',]#can add Blank
Meta.All.POCPOU.Filter<-Meta.All.POCPOU.Filter[Meta.All.POCPOU.Filter$Climate == 'Cfa'|Meta.All.POCPOU.Filter$Climate == 'Bsh'|Meta.All.POCPOU.Filter$Climate == 'Csb',]
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
All.POCPOU.dist.bray<- phyloseq::distance(physeq.combined.All.POCPOU, method = "bray")
All.POCPOU.meta<- data.frame(sample_data(physeq.combined.All.POCPOU))

# ANOSIM.GroupedP.All.POCPOU.bray<-numeric() 
# ADONIS2.GroupedP.All.POCPOU.bray<-numeric()
# Betadis.GroupedP.All.POCPOU.bray<-numeric()
# ANOSIM.GroupedR.All.POCPOU.bray<-numeric() 
# ADONIS2.GroupedR.All.POCPOU.bray<-numeric()

#################
#################
#################
#################
#Climate 
ANOSIM.All.POCPOU.bray.Climate<-anosim(All.POCPOU.dist.bray,
                                      All.POCPOU.meta$Climate, 
                                      permutations = permutations)

ADONIS2.All.POCPOU.bray.Climate<-adonis2(All.POCPOU.dist.bray ~ Climate, 
                                        data = All.POCPOU.meta, 
                                        permutations = permutations)

Betadis.All.POCPOU.bray.Climate<-anova(betadisper(All.POCPOU.dist.bray,
                                                 All.POCPOU.meta$Climate))



ANOSIM.Pval.All.POCPOU.bray.Climate<-ANOSIM.All.POCPOU.bray.Climate$signif
ANOSIM.Pval.All.POCPOU.bray.Climate

ADONIS2.Pval.All.POCPOU.bray.Climate<-ADONIS2.All.POCPOU.bray.Climate$`Pr(>F)`[1]
ADONIS2.Pval.All.POCPOU.bray.Climate  

Betadis.Pval.All.POCPOU.bray.Climate<-Betadis.All.POCPOU.bray.Climate$`Pr(>F)`[1]
Betadis.Pval.All.POCPOU.bray.Climate

ANOSIM.R.All.POCPOU.bray.Climate<-ANOSIM.All.POCPOU.bray.Climate$statistic
ANOSIM.R.All.POCPOU.bray.Climate

ADONIS2.R2.All.POCPOU.bray.Climate<-ADONIS2.All.POCPOU.bray.Climate$R2[1]
ADONIS2.R2.All.POCPOU.bray.Climate  

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.POCPOU.bray["Climate"]<-ANOSIM.Pval.All.POCPOU.bray.Climate
ADONIS2.GroupedP.All.POCPOU.bray["Climate"]<-ADONIS2.Pval.All.POCPOU.bray.Climate

Betadis.GroupedP.All.POCPOU.bray["Climate"]<-Betadis.Pval.All.POCPOU.bray.Climate

ANOSIM.GroupedR.All.POCPOU.bray["Climate"]<-ANOSIM.R.All.POCPOU.bray.Climate
ADONIS2.GroupedR.All.POCPOU.bray["Climate"]<-ADONIS2.R2.All.POCPOU.bray.Climate

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.POCPOU.meta$Climate)

ADONIS2.Pairwise.All.POCPOU.bray.Climate<-pairwise.adonis2(All.POCPOU.dist.bray ~ Climate, 
                                                          data = All.POCPOU.meta, 
                                                          permutations = permutations)
ADONIS2.Pairwise.All.POCPOU.bray.Climate
#Climate





########Climate_Region 
#Read Meta Data 
Metadata_POCPOU<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_DS_Special', rowNames = TRUE)
#Metadata_POC<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_Final', rowNames = TRUE)
Metadata_POC<-Metadata_POCPOU[Metadata_POCPOU$Climate_Region == "POC",]
#Metadata_POU<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'DS_Final', rowNames = TRUE)
Metadata_POU<-Metadata_POCPOU[Metadata_POCPOU$Climate_Region == "POU",]

#Order Meta Data 
Metadata_POCPOU$Classification_1 = factor(Metadata_POCPOU$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))
Metadata_POC$Classification_1 = factor(Metadata_POC$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))
Metadata_POU$Classification_1 = factor(Metadata_POU$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))

######All Samples, POCPOU
Meta.All.POCPOU<-Metadata_POCPOU
Meta.All.POCPOU<-Meta.All.POCPOU[Meta.All.POCPOU$Description != 'ConnerPrimer58',] #removed because of low read count
Meta.All.POCPOU<-Meta.All.POCPOU[Meta.All.POCPOU$Description != 'PR1WEEF',] #removed because of high read count 

Meta.All.POCPOU.Filter<-Meta.All.POCPOU[Meta.All.POCPOU$Classification_1 == 'Potable'| Meta.All.POCPOU$Classification_1 =='Non-potable Reuse'|Meta.All.POCPOU$Classification_1 =='Potable Reuse'|Meta.All.POCPOU$Classification_1 =='',]#can add Blank
Meta.All.POCPOU.Filter<-Meta.All.POCPOU.Filter[Meta.All.POCPOU.Filter$Climate_Region == 'Southeast_Cfa'|Meta.All.POCPOU.Filter$Climate_Region == 'Southwest__Csb'|Meta.All.POCPOU.Filter$Climate_Region == 'Southwest_Bsh',]
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
All.POCPOU.dist.bray<- phyloseq::distance(physeq.combined.All.POCPOU, method = "bray")
All.POCPOU.meta<- data.frame(sample_data(physeq.combined.All.POCPOU))

# ANOSIM.GroupedP.All.POCPOU.bray<-numeric() 
# ADONIS2.GroupedP.All.POCPOU.bray<-numeric()
# Betadis.GroupedP.All.POCPOU.bray<-numeric()
# ANOSIM.GroupedR.All.POCPOU.bray<-numeric() 
# ADONIS2.GroupedR.All.POCPOU.bray<-numeric()

#################
#################
#################
#################
#Climate_Region 
ANOSIM.All.POCPOU.bray.Climate_Region<-anosim(All.POCPOU.dist.bray,
                                       All.POCPOU.meta$Climate_Region, 
                                       permutations = permutations)

ADONIS2.All.POCPOU.bray.Climate_Region<-adonis2(All.POCPOU.dist.bray ~ Climate_Region, 
                                         data = All.POCPOU.meta, 
                                         permutations = permutations)

Betadis.All.POCPOU.bray.Climate_Region<-anova(betadisper(All.POCPOU.dist.bray,
                                                  All.POCPOU.meta$Climate_Region))



ANOSIM.Pval.All.POCPOU.bray.Climate_Region<-ANOSIM.All.POCPOU.bray.Climate_Region$signif
ANOSIM.Pval.All.POCPOU.bray.Climate_Region

ADONIS2.Pval.All.POCPOU.bray.Climate_Region<-ADONIS2.All.POCPOU.bray.Climate_Region$`Pr(>F)`[1]
ADONIS2.Pval.All.POCPOU.bray.Climate_Region  

Betadis.Pval.All.POCPOU.bray.Climate_Region<-Betadis.All.POCPOU.bray.Climate_Region$`Pr(>F)`[1]
Betadis.Pval.All.POCPOU.bray.Climate_Region

ANOSIM.R.All.POCPOU.bray.Climate_Region<-ANOSIM.All.POCPOU.bray.Climate_Region$statistic
ANOSIM.R.All.POCPOU.bray.Climate_Region

ADONIS2.R2.All.POCPOU.bray.Climate_Region<-ADONIS2.All.POCPOU.bray.Climate_Region$R2[1]
ADONIS2.R2.All.POCPOU.bray.Climate_Region  

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.POCPOU.bray["Climate_Region"]<-ANOSIM.Pval.All.POCPOU.bray.Climate_Region
ADONIS2.GroupedP.All.POCPOU.bray["Climate_Region"]<-ADONIS2.Pval.All.POCPOU.bray.Climate_Region

Betadis.GroupedP.All.POCPOU.bray["Climate_Region"]<-Betadis.Pval.All.POCPOU.bray.Climate_Region

ANOSIM.GroupedR.All.POCPOU.bray["Climate_Region"]<-ANOSIM.R.All.POCPOU.bray.Climate_Region
ADONIS2.GroupedR.All.POCPOU.bray["Climate_Region"]<-ADONIS2.R2.All.POCPOU.bray.Climate_Region

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.POCPOU.meta$Climate_Region)

ADONIS2.Pairwise.All.POCPOU.bray.Climate_Region<-pairwise.adonis2(All.POCPOU.dist.bray ~ Climate_Region, 
                                                           data = All.POCPOU.meta, 
                                                           permutations = permutations)
ADONIS2.Pairwise.All.POCPOU.bray.Climate_Region
#Climate_Region


#Adjust PValue, Check Significance
ANOSIM.GroupedP.Adjusted.All.POCPOU.bray<-p.adjust(ANOSIM.GroupedP.All.POCPOU.bray,method="bonferroni")
p.adjust(ANOSIM.GroupedP.Adjusted.All.POCPOU.bray,method="bonferroni")<.05

ADONIS2.GroupedP.Adjusted.All.POCPOU.bray<-p.adjust(ADONIS2.GroupedP.All.POCPOU.bray,method="bonferroni")
p.adjust(ADONIS2.GroupedP.Adjusted.All.POCPOU.bray,method="bonferroni")<.05

Betadis.GroupedP.All.POCPOU.bray<0.05
Betadis.GroupedP.Adjusted.All.POCPOU.bray<-p.adjust(Betadis.GroupedP.All.POCPOU.bray,method="bonferroni")
p.adjust(Betadis.GroupedP.Adjusted.All.POCPOU.bray,method="bonferroni")<.05


#Check Value of Pairwise Sub comparisons 
ADONIS2.Pairwise.All.POCPOU.bray.Class1_Pot_NonPot
ADONIS2.Pairwise.All.POCPOU.bray.Classification_1
ADONIS2.Pairwise.All.POCPOU.bray.Classification_1_mod
ADONIS2.Pairwise.All.POCPOU.bray.Classification_2
ADONIS2.Pairwise.All.POCPOU.bray.Climate
ADONIS2.Pairwise.All.POCPOU.bray.Climate_Region
ADONIS2.Pairwise.All.POCPOU.bray.Final_Disinfection_Residual
ADONIS2.Pairwise.All.POCPOU.bray.Primers
ADONIS2.Pairwise.All.POCPOU.bray.Region
ADONIS2.Pairwise.All.POCPOU.bray.Sample_Local
ADONIS2.Pairwise.All.POCPOU.bray.Trimmed
ADONIS2.Pairwise.All.POCPOU.bray.Extraction_Type

#Permanova on Multiple Variables that are significant and not heterogeneous (not sig from Beta, and sig from anosim)
ADONIS2.All.POCPOU.bray.multiple_comparisons<-adonis2(All.POCPOU.dist.bray ~ 
                                                        Classification_1+Climate, 
                                                 data = All.POCPOU.meta, 
                                                 permutations = permutations)
ADONIS2.All.POCPOU.bray.multiple_comparisons
#Note, in this particular case, I think the permanova approach is not great for the All.POC.POU condition. Focus on it for future tests. 


########Sample_Local, BW removed 

#Read Meta Data 
Metadata_POC<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_DS_Special', rowNames = TRUE)
#Metadata_POC<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_Final', rowNames = TRUE)
Metadata_POC<-Metadata_POC[Metadata_POC$Sample_Local == "POC",]
#Metadata_POU<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'DS_Final', rowNames = TRUE)
Metadata_POU<-Metadata_POC[Metadata_POC$Sample_Local == "POU",]

#Order Meta Data 
Metadata_POC$Classification_1 = factor(Metadata_POC$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))
Metadata_POC$Classification_1 = factor(Metadata_POC$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))
Metadata_POU$Classification_1 = factor(Metadata_POU$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))

######All Samples, POC
Meta.All.POC<-Metadata_POC
Meta.All.POC<-Meta.All.POC[Meta.All.POC$Description != 'ConnerPrimer58',] #removed because of low read count
Meta.All.POC<-Meta.All.POC[Meta.All.POC$Description != 'PR1WEEF',] #removed because of high read count 

Meta.All.POC.Filter<-Meta.All.POC[Meta.All.POC$Classification_1 == 'Potable'| Meta.All.POC$Classification_1 =='Non-potable Reuse'|Meta.All.POC$Classification_1 =='Potable Reuse'|Meta.All.POC$Classification_1 =='',]#can add Blank
Meta.All.POC.Filter<-Meta.All.POC.Filter[Meta.All.POC.Filter$Sample_Local != 'Bottled Water',]
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
All.POC.dist.bray<- phyloseq::distance(physeq.combined.All.POC, method = "bray")
All.POC.meta<- data.frame(sample_data(physeq.combined.All.POC))

ANOSIM.GroupedP.All.POC.bray<-numeric()
ADONIS2.GroupedP.All.POC.bray<-numeric()
Betadis.GroupedP.All.POC.bray<-numeric()
ANOSIM.GroupedR.All.POC.bray<-numeric()
ADONIS2.GroupedR.All.POC.bray<-numeric()

#################
#################
#################
#################
#Sample_Local 
ANOSIM.All.POC.bray.Sample_Local<-anosim(All.POC.dist.bray,
                                            All.POC.meta$Sample_Local, 
                                            permutations = permutations)

ADONIS2.All.POC.bray.Sample_Local<-adonis2(All.POC.dist.bray ~ Sample_Local, 
                                              data = All.POC.meta, 
                                              permutations = permutations)

Betadis.All.POC.bray.Sample_Local<-anova(betadisper(All.POC.dist.bray,
                                                       All.POC.meta$Sample_Local))



ANOSIM.Pval.All.POC.bray.Sample_Local<-ANOSIM.All.POC.bray.Sample_Local$signif
ANOSIM.Pval.All.POC.bray.Sample_Local

ADONIS2.Pval.All.POC.bray.Sample_Local<-ADONIS2.All.POC.bray.Sample_Local$`Pr(>F)`[1]
ADONIS2.Pval.All.POC.bray.Sample_Local  

Betadis.Pval.All.POC.bray.Sample_Local<-Betadis.All.POC.bray.Sample_Local$`Pr(>F)`[1]
Betadis.Pval.All.POC.bray.Sample_Local

ANOSIM.R.All.POC.bray.Sample_Local<-ANOSIM.All.POC.bray.Sample_Local$statistic
ANOSIM.R.All.POC.bray.Sample_Local

ADONIS2.R2.All.POC.bray.Sample_Local<-ADONIS2.All.POC.bray.Sample_Local$R2[1]
ADONIS2.R2.All.POC.bray.Sample_Local  

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.POC.bray["Sample_Local"]<-ANOSIM.Pval.All.POC.bray.Sample_Local
ADONIS2.GroupedP.All.POC.bray["Sample_Local"]<-ADONIS2.Pval.All.POC.bray.Sample_Local

Betadis.GroupedP.All.POC.bray["Sample_Local"]<-Betadis.Pval.All.POC.bray.Sample_Local

ANOSIM.GroupedR.All.POC.bray["Sample_Local"]<-ANOSIM.R.All.POC.bray.Sample_Local
ADONIS2.GroupedR.All.POC.bray["Sample_Local"]<-ADONIS2.R2.All.POC.bray.Sample_Local

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.POC.meta$Sample_Local)

ADONIS2.Pairwise.All.POC.bray.Sample_Local<-pairwise.adonis2(All.POC.dist.bray ~ Sample_Local, 
                                                                data = All.POC.meta, 
                                                                permutations = permutations)
ADONIS2.Pairwise.All.POC.bray.Sample_Local
#Sample_Local














########Region 

#Read Meta Data 
Metadata_POC<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_DS_Special', rowNames = TRUE)
#Metadata_POC<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_Final', rowNames = TRUE)
Metadata_POC<-Metadata_POC[Metadata_POC$Region == "POC",]
#Metadata_POU<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'DS_Final', rowNames = TRUE)
Metadata_POU<-Metadata_POC[Metadata_POC$Region == "POU",]

#Order Meta Data 
Metadata_POC$Classification_1 = factor(Metadata_POC$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))
Metadata_POC$Classification_1 = factor(Metadata_POC$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))
Metadata_POU$Classification_1 = factor(Metadata_POU$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))

######All Samples, POC
Meta.All.POC<-Metadata_POC
Meta.All.POC<-Meta.All.POC[Meta.All.POC$Description != 'ConnerPrimer58',] #removed because of low read count
Meta.All.POC<-Meta.All.POC[Meta.All.POC$Description != 'PR1WEEF',] #removed because of high read count 

Meta.All.POC.Filter<-Meta.All.POC[Meta.All.POC$Classification_1 == 'Potable'| Meta.All.POC$Classification_1 =='Non-potable Reuse'|Meta.All.POC$Classification_1 =='Potable Reuse'|Meta.All.POC$Classification_1 =='',]#can add Blank
Meta.All.POC.Filter<-Meta.All.POC.Filter[Meta.All.POC.Filter$Region == 'Southeast'|Meta.All.POC.Filter$Region == 'West'|Meta.All.POC.Filter$Region == 'Southwest',]
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
All.POC.dist.bray<- phyloseq::distance(physeq.combined.All.POC, method = "bray")
All.POC.meta<- data.frame(sample_data(physeq.combined.All.POC))

# ANOSIM.GroupedP.All.POC.bray<-numeric() 
# ADONIS2.GroupedP.All.POC.bray<-numeric()
# Betadis.GroupedP.All.POC.bray<-numeric()
# ANOSIM.GroupedR.All.POC.bray<-numeric() 
# ADONIS2.GroupedR.All.POC.bray<-numeric()

#################
#################
#################
#################
#Region 
ANOSIM.All.POC.bray.Region<-anosim(All.POC.dist.bray,
                                      All.POC.meta$Region, 
                                      permutations = permutations)

ADONIS2.All.POC.bray.Region<-adonis2(All.POC.dist.bray ~ Region, 
                                        data = All.POC.meta, 
                                        permutations = permutations)

Betadis.All.POC.bray.Region<-anova(betadisper(All.POC.dist.bray,
                                                 All.POC.meta$Region))



ANOSIM.Pval.All.POC.bray.Region<-ANOSIM.All.POC.bray.Region$signif
ANOSIM.Pval.All.POC.bray.Region

ADONIS2.Pval.All.POC.bray.Region<-ADONIS2.All.POC.bray.Region$`Pr(>F)`[1]
ADONIS2.Pval.All.POC.bray.Region  

Betadis.Pval.All.POC.bray.Region<-Betadis.All.POC.bray.Region$`Pr(>F)`[1]
Betadis.Pval.All.POC.bray.Region

ANOSIM.R.All.POC.bray.Region<-ANOSIM.All.POC.bray.Region$statistic
ANOSIM.R.All.POC.bray.Region

ADONIS2.R2.All.POC.bray.Region<-ADONIS2.All.POC.bray.Region$R2[1]
ADONIS2.R2.All.POC.bray.Region  

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.POC.bray["Region"]<-ANOSIM.Pval.All.POC.bray.Region
ADONIS2.GroupedP.All.POC.bray["Region"]<-ADONIS2.Pval.All.POC.bray.Region

Betadis.GroupedP.All.POC.bray["Region"]<-Betadis.Pval.All.POC.bray.Region

ANOSIM.GroupedR.All.POC.bray["Region"]<-ANOSIM.R.All.POC.bray.Region
ADONIS2.GroupedR.All.POC.bray["Region"]<-ADONIS2.R2.All.POC.bray.Region

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.POC.meta$Region)

ADONIS2.Pairwise.All.POC.bray.Region<-pairwise.adonis2(All.POC.dist.bray ~ Region, 
                                                          data = All.POC.meta, 
                                                          permutations = permutations)
ADONIS2.Pairwise.All.POC.bray.Region
#Region










########Climate 
#Read Meta Data 
Metadata_POC<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_DS_Special', rowNames = TRUE)
#Metadata_POC<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_Final', rowNames = TRUE)
Metadata_POC<-Metadata_POC[Metadata_POC$Climate == "POC",]
#Metadata_POU<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'DS_Final', rowNames = TRUE)
Metadata_POU<-Metadata_POC[Metadata_POC$Climate == "POU",]

#Order Meta Data 
Metadata_POC$Classification_1 = factor(Metadata_POC$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))
Metadata_POC$Classification_1 = factor(Metadata_POC$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))
Metadata_POU$Classification_1 = factor(Metadata_POU$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))

######All Samples, POC
Meta.All.POC<-Metadata_POC
Meta.All.POC<-Meta.All.POC[Meta.All.POC$Description != 'ConnerPrimer58',] #removed because of low read count
Meta.All.POC<-Meta.All.POC[Meta.All.POC$Description != 'PR1WEEF',] #removed because of high read count 

Meta.All.POC.Filter<-Meta.All.POC[Meta.All.POC$Classification_1 == 'Potable'| Meta.All.POC$Classification_1 =='Non-potable Reuse'|Meta.All.POC$Classification_1 =='Potable Reuse'|Meta.All.POC$Classification_1 =='',]#can add Blank
Meta.All.POC.Filter<-Meta.All.POC.Filter[Meta.All.POC.Filter$Climate == 'Cfa'|Meta.All.POC.Filter$Climate == 'Bsh'|Meta.All.POC.Filter$Climate == 'Csb',]
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
All.POC.dist.bray<- phyloseq::distance(physeq.combined.All.POC, method = "bray")
All.POC.meta<- data.frame(sample_data(physeq.combined.All.POC))

# ANOSIM.GroupedP.All.POC.bray<-numeric() 
# ADONIS2.GroupedP.All.POC.bray<-numeric()
# Betadis.GroupedP.All.POC.bray<-numeric()
# ANOSIM.GroupedR.All.POC.bray<-numeric() 
# ADONIS2.GroupedR.All.POC.bray<-numeric()

#################
#################
#################
#################
#Climate 
ANOSIM.All.POC.bray.Climate<-anosim(All.POC.dist.bray,
                                       All.POC.meta$Climate, 
                                       permutations = permutations)

ADONIS2.All.POC.bray.Climate<-adonis2(All.POC.dist.bray ~ Climate, 
                                         data = All.POC.meta, 
                                         permutations = permutations)

Betadis.All.POC.bray.Climate<-anova(betadisper(All.POC.dist.bray,
                                                  All.POC.meta$Climate))



ANOSIM.Pval.All.POC.bray.Climate<-ANOSIM.All.POC.bray.Climate$signif
ANOSIM.Pval.All.POC.bray.Climate

ADONIS2.Pval.All.POC.bray.Climate<-ADONIS2.All.POC.bray.Climate$`Pr(>F)`[1]
ADONIS2.Pval.All.POC.bray.Climate  

Betadis.Pval.All.POC.bray.Climate<-Betadis.All.POC.bray.Climate$`Pr(>F)`[1]
Betadis.Pval.All.POC.bray.Climate

ANOSIM.R.All.POC.bray.Climate<-ANOSIM.All.POC.bray.Climate$statistic
ANOSIM.R.All.POC.bray.Climate

ADONIS2.R2.All.POC.bray.Climate<-ADONIS2.All.POC.bray.Climate$R2[1]
ADONIS2.R2.All.POC.bray.Climate  

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.POC.bray["Climate"]<-ANOSIM.Pval.All.POC.bray.Climate
ADONIS2.GroupedP.All.POC.bray["Climate"]<-ADONIS2.Pval.All.POC.bray.Climate

Betadis.GroupedP.All.POC.bray["Climate"]<-Betadis.Pval.All.POC.bray.Climate

ANOSIM.GroupedR.All.POC.bray["Climate"]<-ANOSIM.R.All.POC.bray.Climate
ADONIS2.GroupedR.All.POC.bray["Climate"]<-ADONIS2.R2.All.POC.bray.Climate

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.POC.meta$Climate)

ADONIS2.Pairwise.All.POC.bray.Climate<-pairwise.adonis2(All.POC.dist.bray ~ Climate, 
                                                           data = All.POC.meta, 
                                                           permutations = permutations)
ADONIS2.Pairwise.All.POC.bray.Climate
#Climate





########Climate_Region 
#Read Meta Data 
Metadata_POC<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_DS_Special', rowNames = TRUE)
#Metadata_POC<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_Final', rowNames = TRUE)
Metadata_POC<-Metadata_POC[Metadata_POC$Climate_Region == "POC",]
#Metadata_POU<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'DS_Final', rowNames = TRUE)
Metadata_POU<-Metadata_POC[Metadata_POC$Climate_Region == "POU",]

#Order Meta Data 
Metadata_POC$Classification_1 = factor(Metadata_POC$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))
Metadata_POC$Classification_1 = factor(Metadata_POC$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))
Metadata_POU$Classification_1 = factor(Metadata_POU$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))

######All Samples, POC
Meta.All.POC<-Metadata_POC
Meta.All.POC<-Meta.All.POC[Meta.All.POC$Description != 'ConnerPrimer58',] #removed because of low read count
Meta.All.POC<-Meta.All.POC[Meta.All.POC$Description != 'PR1WEEF',] #removed because of high read count 

Meta.All.POC.Filter<-Meta.All.POC[Meta.All.POC$Classification_1 == 'Potable'| Meta.All.POC$Classification_1 =='Non-potable Reuse'|Meta.All.POC$Classification_1 =='Potable Reuse'|Meta.All.POC$Classification_1 =='',]#can add Blank
Meta.All.POC.Filter<-Meta.All.POC.Filter[Meta.All.POC.Filter$Climate_Region == 'Southeast_Cfa'|Meta.All.POC.Filter$Climate_Region == 'Southwest__Csb'|Meta.All.POC.Filter$Climate_Region == 'Southwest_Bsh',]
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
All.POC.dist.bray<- phyloseq::distance(physeq.combined.All.POC, method = "bray")
All.POC.meta<- data.frame(sample_data(physeq.combined.All.POC))

# ANOSIM.GroupedP.All.POC.bray<-numeric() 
# ADONIS2.GroupedP.All.POC.bray<-numeric()
# Betadis.GroupedP.All.POC.bray<-numeric()
# ANOSIM.GroupedR.All.POC.bray<-numeric() 
# ADONIS2.GroupedR.All.POC.bray<-numeric()

#################
#################
#################
#################
#Climate_Region 
ANOSIM.All.POC.bray.Climate_Region<-anosim(All.POC.dist.bray,
                                              All.POC.meta$Climate_Region, 
                                              permutations = permutations)

ADONIS2.All.POC.bray.Climate_Region<-adonis2(All.POC.dist.bray ~ Climate_Region, 
                                                data = All.POC.meta, 
                                                permutations = permutations)

Betadis.All.POC.bray.Climate_Region<-anova(betadisper(All.POC.dist.bray,
                                                         All.POC.meta$Climate_Region))



ANOSIM.Pval.All.POC.bray.Climate_Region<-ANOSIM.All.POC.bray.Climate_Region$signif
ANOSIM.Pval.All.POC.bray.Climate_Region

ADONIS2.Pval.All.POC.bray.Climate_Region<-ADONIS2.All.POC.bray.Climate_Region$`Pr(>F)`[1]
ADONIS2.Pval.All.POC.bray.Climate_Region  

Betadis.Pval.All.POC.bray.Climate_Region<-Betadis.All.POC.bray.Climate_Region$`Pr(>F)`[1]
Betadis.Pval.All.POC.bray.Climate_Region

ANOSIM.R.All.POC.bray.Climate_Region<-ANOSIM.All.POC.bray.Climate_Region$statistic
ANOSIM.R.All.POC.bray.Climate_Region

ADONIS2.R2.All.POC.bray.Climate_Region<-ADONIS2.All.POC.bray.Climate_Region$R2[1]
ADONIS2.R2.All.POC.bray.Climate_Region  

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.POC.bray["Climate_Region"]<-ANOSIM.Pval.All.POC.bray.Climate_Region
ADONIS2.GroupedP.All.POC.bray["Climate_Region"]<-ADONIS2.Pval.All.POC.bray.Climate_Region

Betadis.GroupedP.All.POC.bray["Climate_Region"]<-Betadis.Pval.All.POC.bray.Climate_Region

ANOSIM.GroupedR.All.POC.bray["Climate_Region"]<-ANOSIM.R.All.POC.bray.Climate_Region
ADONIS2.GroupedR.All.POC.bray["Climate_Region"]<-ADONIS2.R2.All.POC.bray.Climate_Region

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.POC.meta$Climate_Region)

ADONIS2.Pairwise.All.POC.bray.Climate_Region<-pairwise.adonis2(All.POC.dist.bray ~ Climate_Region, 
                                                                  data = All.POC.meta, 
                                                                  permutations = permutations)
ADONIS2.Pairwise.All.POC.bray.Climate_Region
#Climate_Region


#Adjust PValue, Check Significance
ANOSIM.GroupedP.Adjusted.All.POC.bray<-p.adjust(ANOSIM.GroupedP.All.POC.bray,method="bonferroni")
p.adjust(ANOSIM.GroupedP.Adjusted.All.POC.bray,method="bonferroni")<.05

ADONIS2.GroupedP.Adjusted.All.POC.bray<-p.adjust(ADONIS2.GroupedP.All.POC.bray,method="bonferroni")
p.adjust(ADONIS2.GroupedP.Adjusted.All.POC.bray,method="bonferroni")<.05

Betadis.GroupedP.All.POC.bray<0.05
Betadis.GroupedP.Adjusted.All.POC.bray<-p.adjust(Betadis.GroupedP.All.POC.bray,method="bonferroni")
p.adjust(Betadis.GroupedP.Adjusted.All.POC.bray,method="bonferroni")<.05


#Check Value of Pairwise Sub comparisons 
ADONIS2.Pairwise.All.POC.bray.Class1_Pot_NonPot
ADONIS2.Pairwise.All.POC.bray.Classification_1
ADONIS2.Pairwise.All.POC.bray.Classification_1_mod
ADONIS2.Pairwise.All.POC.bray.Classification_2
ADONIS2.Pairwise.All.POC.bray.Climate
ADONIS2.Pairwise.All.POC.bray.Climate_Region
ADONIS2.Pairwise.All.POC.bray.Final_Disinfection_Residual
ADONIS2.Pairwise.All.POC.bray.Primers
ADONIS2.Pairwise.All.POC.bray.Region
ADONIS2.Pairwise.All.POC.bray.Sample_Local
ADONIS2.Pairwise.All.POC.bray.Trimmed
ADONIS2.Pairwise.All.POC.bray.Extraction_Type

#Permanova on Multiple Variables that are significant and not heterogeneous (not sig from Beta, and sig from anosim)
ADONIS2.All.POC.bray.multiple_comparisons<-adonis2(All.POC.dist.bray ~ 
                                                        Classification_1+Climate, 
                                                      data = All.POC.meta, 
                                                      permutations = permutations)
ADONIS2.All.POC.bray.multiple_comparisons
#Note, in this particular case, I think the permanova approach is not great for the All.POC.POU condition. Focus on it for future tests. 




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
######All Samples, POU
########Sample_Local, BW removed 

#Read Meta Data 
Metadata_POU<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_DS_Special', rowNames = TRUE)
#Metadata_POC<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_Final', rowNames = TRUE)
Metadata_POC<-Metadata_POU[Metadata_POU$Sample_Local == "POC",]
#Metadata_POU<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'DS_Final', rowNames = TRUE)
Metadata_POU<-Metadata_POU[Metadata_POU$Sample_Local == "POU",]

#Order Meta Data 
Metadata_POU$Classification_1 = factor(Metadata_POU$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))
Metadata_POC$Classification_1 = factor(Metadata_POC$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))
Metadata_POU$Classification_1 = factor(Metadata_POU$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))

######All Samples, POU
Meta.All.POU<-Metadata_POU
Meta.All.POU<-Meta.All.POU[Meta.All.POU$Description != 'ConnerPrimer58',] #removed because of low read count
Meta.All.POU<-Meta.All.POU[Meta.All.POU$Description != 'PR1WEEF',] #removed because of high read count 

Meta.All.POU.Filter<-Meta.All.POU[Meta.All.POU$Classification_1 == 'Potable'| Meta.All.POU$Classification_1 =='Non-potable Reuse'|Meta.All.POU$Classification_1 =='Potable Reuse'|Meta.All.POU$Classification_1 =='',]#can add Blank
Meta.All.POU.Filter<-Meta.All.POU.Filter[Meta.All.POU.Filter$Sample_Local != 'Bottled Water',]
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
All.POU.dist.bray<- phyloseq::distance(physeq.combined.All.POU, method = "bray")
All.POU.meta<- data.frame(sample_data(physeq.combined.All.POU))

ANOSIM.GroupedP.All.POU.bray<-numeric()
ADONIS2.GroupedP.All.POU.bray<-numeric()
Betadis.GroupedP.All.POU.bray<-numeric()
ANOSIM.GroupedR.All.POU.bray<-numeric()
ADONIS2.GroupedR.All.POU.bray<-numeric()

#################
#################
#################
#################
#Sample_Local 
ANOSIM.All.POU.bray.Sample_Local<-anosim(All.POU.dist.bray,
                                            All.POU.meta$Sample_Local, 
                                            permutations = permutations)

ADONIS2.All.POU.bray.Sample_Local<-adonis2(All.POU.dist.bray ~ Sample_Local, 
                                              data = All.POU.meta, 
                                              permutations = permutations)

Betadis.All.POU.bray.Sample_Local<-anova(betadisper(All.POU.dist.bray,
                                                       All.POU.meta$Sample_Local))



ANOSIM.Pval.All.POU.bray.Sample_Local<-ANOSIM.All.POU.bray.Sample_Local$signif
ANOSIM.Pval.All.POU.bray.Sample_Local

ADONIS2.Pval.All.POU.bray.Sample_Local<-ADONIS2.All.POU.bray.Sample_Local$`Pr(>F)`[1]
ADONIS2.Pval.All.POU.bray.Sample_Local  

Betadis.Pval.All.POU.bray.Sample_Local<-Betadis.All.POU.bray.Sample_Local$`Pr(>F)`[1]
Betadis.Pval.All.POU.bray.Sample_Local

ANOSIM.R.All.POU.bray.Sample_Local<-ANOSIM.All.POU.bray.Sample_Local$statistic
ANOSIM.R.All.POU.bray.Sample_Local

ADONIS2.R2.All.POU.bray.Sample_Local<-ADONIS2.All.POU.bray.Sample_Local$R2[1]
ADONIS2.R2.All.POU.bray.Sample_Local  

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.POU.bray["Sample_Local"]<-ANOSIM.Pval.All.POU.bray.Sample_Local
ADONIS2.GroupedP.All.POU.bray["Sample_Local"]<-ADONIS2.Pval.All.POU.bray.Sample_Local

Betadis.GroupedP.All.POU.bray["Sample_Local"]<-Betadis.Pval.All.POU.bray.Sample_Local

ANOSIM.GroupedR.All.POU.bray["Sample_Local"]<-ANOSIM.R.All.POU.bray.Sample_Local
ADONIS2.GroupedR.All.POU.bray["Sample_Local"]<-ADONIS2.R2.All.POU.bray.Sample_Local

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.POU.meta$Sample_Local)

ADONIS2.Pairwise.All.POU.bray.Sample_Local<-pairwise.adonis2(All.POU.dist.bray ~ Sample_Local, 
                                                                data = All.POU.meta, 
                                                                permutations = permutations)
ADONIS2.Pairwise.All.POU.bray.Sample_Local
#Sample_Local














########Region 

#Read Meta Data 
Metadata_POU<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_DS_Special', rowNames = TRUE)
#Metadata_POC<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_Final', rowNames = TRUE)
Metadata_POC<-Metadata_POU[Metadata_POU$Region == "POC",]
#Metadata_POU<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'DS_Final', rowNames = TRUE)
Metadata_POU<-Metadata_POU[Metadata_POU$Region == "POU",]

#Order Meta Data 
Metadata_POU$Classification_1 = factor(Metadata_POU$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))
Metadata_POC$Classification_1 = factor(Metadata_POC$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))
Metadata_POU$Classification_1 = factor(Metadata_POU$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))

######All Samples, POU
Meta.All.POU<-Metadata_POU
Meta.All.POU<-Meta.All.POU[Meta.All.POU$Description != 'ConnerPrimer58',] #removed because of low read count
Meta.All.POU<-Meta.All.POU[Meta.All.POU$Description != 'PR1WEEF',] #removed because of high read count 

Meta.All.POU.Filter<-Meta.All.POU[Meta.All.POU$Classification_1 == 'Potable'| Meta.All.POU$Classification_1 =='Non-potable Reuse'|Meta.All.POU$Classification_1 =='Potable Reuse'|Meta.All.POU$Classification_1 =='',]#can add Blank
Meta.All.POU.Filter<-Meta.All.POU.Filter[Meta.All.POU.Filter$Region == 'Southeast'|Meta.All.POU.Filter$Region == 'West'|Meta.All.POU.Filter$Region == 'Southwest',]
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
All.POU.dist.bray<- phyloseq::distance(physeq.combined.All.POU, method = "bray")
All.POU.meta<- data.frame(sample_data(physeq.combined.All.POU))
# 
# ANOSIM.GroupedP.All.POU.bray<-numeric() 
# ADONIS2.GroupedP.All.POU.bray<-numeric()
# Betadis.GroupedP.All.POU.bray<-numeric()
# ANOSIM.GroupedR.All.POU.bray<-numeric() 
# ADONIS2.GroupedR.All.POU.bray<-numeric()

#################
#################
#################
#################
#Region 
ANOSIM.All.POU.bray.Region<-anosim(All.POU.dist.bray,
                                      All.POU.meta$Region, 
                                      permutations = permutations)

ADONIS2.All.POU.bray.Region<-adonis2(All.POU.dist.bray ~ Region, 
                                        data = All.POU.meta, 
                                        permutations = permutations)

Betadis.All.POU.bray.Region<-anova(betadisper(All.POU.dist.bray,
                                                 All.POU.meta$Region))



ANOSIM.Pval.All.POU.bray.Region<-ANOSIM.All.POU.bray.Region$signif
ANOSIM.Pval.All.POU.bray.Region

ADONIS2.Pval.All.POU.bray.Region<-ADONIS2.All.POU.bray.Region$`Pr(>F)`[1]
ADONIS2.Pval.All.POU.bray.Region  

Betadis.Pval.All.POU.bray.Region<-Betadis.All.POU.bray.Region$`Pr(>F)`[1]
Betadis.Pval.All.POU.bray.Region

ANOSIM.R.All.POU.bray.Region<-ANOSIM.All.POU.bray.Region$statistic
ANOSIM.R.All.POU.bray.Region

ADONIS2.R2.All.POU.bray.Region<-ADONIS2.All.POU.bray.Region$R2[1]
ADONIS2.R2.All.POU.bray.Region  

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.POU.bray["Region"]<-ANOSIM.Pval.All.POU.bray.Region
ADONIS2.GroupedP.All.POU.bray["Region"]<-ADONIS2.Pval.All.POU.bray.Region

Betadis.GroupedP.All.POU.bray["Region"]<-Betadis.Pval.All.POU.bray.Region

ANOSIM.GroupedR.All.POU.bray["Region"]<-ANOSIM.R.All.POU.bray.Region
ADONIS2.GroupedR.All.POU.bray["Region"]<-ADONIS2.R2.All.POU.bray.Region

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.POU.meta$Region)

ADONIS2.Pairwise.All.POU.bray.Region<-pairwise.adonis2(All.POU.dist.bray ~ Region, 
                                                          data = All.POU.meta, 
                                                          permutations = permutations)
ADONIS2.Pairwise.All.POU.bray.Region
#Region










########Climate 
#Read Meta Data 
Metadata_POU<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_DS_Special', rowNames = TRUE)
#Metadata_POC<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_Final', rowNames = TRUE)
Metadata_POC<-Metadata_POU[Metadata_POU$Climate == "POC",]
#Metadata_POU<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'DS_Final', rowNames = TRUE)
Metadata_POU<-Metadata_POU[Metadata_POU$Climate == "POU",]

#Order Meta Data 
Metadata_POU$Classification_1 = factor(Metadata_POU$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))
Metadata_POC$Classification_1 = factor(Metadata_POC$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))
Metadata_POU$Classification_1 = factor(Metadata_POU$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))

######All Samples, POU
Meta.All.POU<-Metadata_POU
Meta.All.POU<-Meta.All.POU[Meta.All.POU$Description != 'ConnerPrimer58',] #removed because of low read count
Meta.All.POU<-Meta.All.POU[Meta.All.POU$Description != 'PR1WEEF',] #removed because of high read count 

Meta.All.POU.Filter<-Meta.All.POU[Meta.All.POU$Classification_1 == 'Potable'| Meta.All.POU$Classification_1 =='Non-potable Reuse'|Meta.All.POU$Classification_1 =='Potable Reuse'|Meta.All.POU$Classification_1 =='',]#can add Blank
Meta.All.POU.Filter<-Meta.All.POU.Filter[Meta.All.POU.Filter$Climate == 'Cfa'|Meta.All.POU.Filter$Climate == 'Bsh'|Meta.All.POU.Filter$Climate == 'Csb',]
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
All.POU.dist.bray<- phyloseq::distance(physeq.combined.All.POU, method = "bray")
All.POU.meta<- data.frame(sample_data(physeq.combined.All.POU))

# ANOSIM.GroupedP.All.POU.bray<-numeric() 
# ADONIS2.GroupedP.All.POU.bray<-numeric()
# Betadis.GroupedP.All.POU.bray<-numeric()
# ANOSIM.GroupedR.All.POU.bray<-numeric() 
# ADONIS2.GroupedR.All.POU.bray<-numeric()

#################
#################
#################
#################
#Climate 
ANOSIM.All.POU.bray.Climate<-anosim(All.POU.dist.bray,
                                       All.POU.meta$Climate, 
                                       permutations = permutations)

ADONIS2.All.POU.bray.Climate<-adonis2(All.POU.dist.bray ~ Climate, 
                                         data = All.POU.meta, 
                                         permutations = permutations)

Betadis.All.POU.bray.Climate<-anova(betadisper(All.POU.dist.bray,
                                                  All.POU.meta$Climate))



ANOSIM.Pval.All.POU.bray.Climate<-ANOSIM.All.POU.bray.Climate$signif
ANOSIM.Pval.All.POU.bray.Climate

ADONIS2.Pval.All.POU.bray.Climate<-ADONIS2.All.POU.bray.Climate$`Pr(>F)`[1]
ADONIS2.Pval.All.POU.bray.Climate  

Betadis.Pval.All.POU.bray.Climate<-Betadis.All.POU.bray.Climate$`Pr(>F)`[1]
Betadis.Pval.All.POU.bray.Climate

ANOSIM.R.All.POU.bray.Climate<-ANOSIM.All.POU.bray.Climate$statistic
ANOSIM.R.All.POU.bray.Climate

ADONIS2.R2.All.POU.bray.Climate<-ADONIS2.All.POU.bray.Climate$R2[1]
ADONIS2.R2.All.POU.bray.Climate  

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.POU.bray["Climate"]<-ANOSIM.Pval.All.POU.bray.Climate
ADONIS2.GroupedP.All.POU.bray["Climate"]<-ADONIS2.Pval.All.POU.bray.Climate

Betadis.GroupedP.All.POU.bray["Climate"]<-Betadis.Pval.All.POU.bray.Climate

ANOSIM.GroupedR.All.POU.bray["Climate"]<-ANOSIM.R.All.POU.bray.Climate
ADONIS2.GroupedR.All.POU.bray["Climate"]<-ADONIS2.R2.All.POU.bray.Climate

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.POU.meta$Climate)

ADONIS2.Pairwise.All.POU.bray.Climate<-pairwise.adonis2(All.POU.dist.bray ~ Climate, 
                                                           data = All.POU.meta, 
                                                           permutations = permutations)
ADONIS2.Pairwise.All.POU.bray.Climate
#Climate





########Climate_Region 
#Read Meta Data 
Metadata_POU<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_DS_Special', rowNames = TRUE)
#Metadata_POC<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_Final', rowNames = TRUE)
Metadata_POC<-Metadata_POU[Metadata_POU$Climate_Region == "POC",]
#Metadata_POU<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'DS_Final', rowNames = TRUE)
Metadata_POU<-Metadata_POU[Metadata_POU$Climate_Region == "POU",]

#Order Meta Data 
Metadata_POU$Classification_1 = factor(Metadata_POU$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))
Metadata_POC$Classification_1 = factor(Metadata_POC$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))
Metadata_POU$Classification_1 = factor(Metadata_POU$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))

######All Samples, POU
Meta.All.POU<-Metadata_POU
Meta.All.POU<-Meta.All.POU[Meta.All.POU$Description != 'ConnerPrimer58',] #removed because of low read count
Meta.All.POU<-Meta.All.POU[Meta.All.POU$Description != 'PR1WEEF',] #removed because of high read count 

Meta.All.POU.Filter<-Meta.All.POU[Meta.All.POU$Classification_1 == 'Potable'| Meta.All.POU$Classification_1 =='Non-potable Reuse'|Meta.All.POU$Classification_1 =='Potable Reuse'|Meta.All.POU$Classification_1 =='',]#can add Blank
Meta.All.POU.Filter<-Meta.All.POU.Filter[Meta.All.POU.Filter$Climate_Region == 'Southeast_Cfa'|Meta.All.POU.Filter$Climate_Region == 'Southwest__Csb'|Meta.All.POU.Filter$Climate_Region == 'Southwest_Bsh',]
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
All.POU.dist.bray<- phyloseq::distance(physeq.combined.All.POU, method = "bray")
All.POU.meta<- data.frame(sample_data(physeq.combined.All.POU))

# ANOSIM.GroupedP.All.POU.bray<-numeric() 
# ADONIS2.GroupedP.All.POU.bray<-numeric()
# Betadis.GroupedP.All.POU.bray<-numeric()
# ANOSIM.GroupedR.All.POU.bray<-numeric() 
# ADONIS2.GroupedR.All.POU.bray<-numeric()

#################
#################
#################
#################
#Climate_Region 
ANOSIM.All.POU.bray.Climate_Region<-anosim(All.POU.dist.bray,
                                              All.POU.meta$Climate_Region, 
                                              permutations = permutations)

ADONIS2.All.POU.bray.Climate_Region<-adonis2(All.POU.dist.bray ~ Climate_Region, 
                                                data = All.POU.meta, 
                                                permutations = permutations)

Betadis.All.POU.bray.Climate_Region<-anova(betadisper(All.POU.dist.bray,
                                                         All.POU.meta$Climate_Region))



ANOSIM.Pval.All.POU.bray.Climate_Region<-ANOSIM.All.POU.bray.Climate_Region$signif
ANOSIM.Pval.All.POU.bray.Climate_Region

ADONIS2.Pval.All.POU.bray.Climate_Region<-ADONIS2.All.POU.bray.Climate_Region$`Pr(>F)`[1]
ADONIS2.Pval.All.POU.bray.Climate_Region  

Betadis.Pval.All.POU.bray.Climate_Region<-Betadis.All.POU.bray.Climate_Region$`Pr(>F)`[1]
Betadis.Pval.All.POU.bray.Climate_Region

ANOSIM.R.All.POU.bray.Climate_Region<-ANOSIM.All.POU.bray.Climate_Region$statistic
ANOSIM.R.All.POU.bray.Climate_Region

ADONIS2.R2.All.POU.bray.Climate_Region<-ADONIS2.All.POU.bray.Climate_Region$R2[1]
ADONIS2.R2.All.POU.bray.Climate_Region  

#Pull Together Pvalues to adjust 
ANOSIM.GroupedP.All.POU.bray["Climate_Region"]<-ANOSIM.Pval.All.POU.bray.Climate_Region
ADONIS2.GroupedP.All.POU.bray["Climate_Region"]<-ADONIS2.Pval.All.POU.bray.Climate_Region

Betadis.GroupedP.All.POU.bray["Climate_Region"]<-Betadis.Pval.All.POU.bray.Climate_Region

ANOSIM.GroupedR.All.POU.bray["Climate_Region"]<-ANOSIM.R.All.POU.bray.Climate_Region
ADONIS2.GroupedR.All.POU.bray["Climate_Region"]<-ADONIS2.R2.All.POU.bray.Climate_Region

#Pairwise if needed, check if interesting or not. Remove if not. 
count(All.POU.meta$Climate_Region)

ADONIS2.Pairwise.All.POU.bray.Climate_Region<-pairwise.adonis2(All.POU.dist.bray ~ Climate_Region, 
                                                                  data = All.POU.meta, 
                                                                  permutations = permutations)
ADONIS2.Pairwise.All.POU.bray.Climate_Region
#Climate_Region


#Adjust PValue, Check Significance
ANOSIM.GroupedP.Adjusted.All.POU.bray<-p.adjust(ANOSIM.GroupedP.All.POU.bray,method="bonferroni")
p.adjust(ANOSIM.GroupedP.Adjusted.All.POU.bray,method="bonferroni")<.05

ADONIS2.GroupedP.Adjusted.All.POU.bray<-p.adjust(ADONIS2.GroupedP.All.POU.bray,method="bonferroni")
p.adjust(ADONIS2.GroupedP.Adjusted.All.POU.bray,method="bonferroni")<.05

Betadis.GroupedP.All.POU.bray<0.05
Betadis.GroupedP.Adjusted.All.POU.bray<-p.adjust(Betadis.GroupedP.All.POU.bray,method="bonferroni")
p.adjust(Betadis.GroupedP.Adjusted.All.POU.bray,method="bonferroni")<.05


#Check Value of Pairwise Sub comparisons 
ADONIS2.Pairwise.All.POU.bray.Class1_Pot_NonPot
ADONIS2.Pairwise.All.POU.bray.Classification_1
ADONIS2.Pairwise.All.POU.bray.Classification_1_mod
ADONIS2.Pairwise.All.POU.bray.Classification_2
ADONIS2.Pairwise.All.POU.bray.Climate
ADONIS2.Pairwise.All.POU.bray.Climate_Region
ADONIS2.Pairwise.All.POU.bray.Final_Disinfection_Residual
ADONIS2.Pairwise.All.POU.bray.Primers
ADONIS2.Pairwise.All.POU.bray.Region
ADONIS2.Pairwise.All.POU.bray.Sample_Local
ADONIS2.Pairwise.All.POU.bray.Trimmed
ADONIS2.Pairwise.All.POU.bray.Extraction_Type

#Permanova on Multiple Variables that are significant and not heterogeneous (not sig from Beta, and sig from anosim)
ADONIS2.All.POU.bray.multiple_comparisons<-adonis2(All.POU.dist.bray ~ 
                                                        Classification_1+Climate, 
                                                      data = All.POU.meta, 
                                                      permutations = permutations)
ADONIS2.All.POU.bray.multiple_comparisons






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
All.POCPOU.Stats<-data.frame(
  "ANOSIM.P.All.POCPOU.bray"= ANOSIM.GroupedP.All.POCPOU.bray,
  "ANOSIM.P.Adj.All.POCPOU.bray"= ANOSIM.GroupedP.Adjusted.All.POCPOU.bray,
  "ADONIS2.P.All.POCPOU.bray"= ADONIS2.GroupedP.All.POCPOU.bray,
  "ADONIS2.P.Adj.All.POCPOU.bray"= ADONIS2.GroupedP.Adjusted.All.POCPOU.bray,
  "Betadis.P.All.POCPOU.bray"= Betadis.GroupedP.All.POCPOU.bray,
  "Betadis.P.Adj.All.POCPOU.bray"= Betadis.GroupedP.Adjusted.All.POCPOU.bray,
  "ANOSIM.R.All.POCPOU.bray"=  ANOSIM.GroupedR.All.POCPOU.bray,
  "ADONIS2.R.All.POCPOU.bray"= ADONIS2.GroupedR.All.POCPOU.bray)

# All.POCPOU.List<-list(
# "ANOSIM.P.All.POCPOU.bray"= ANOSIM.GroupedP.All.POCPOU.bray,
# "ANOSIM.P.Adj.All.POCPOU.bray"= ANOSIM.GroupedP.Adjusted.All.POCPOU.bray,
# "ADONIS2.P.All.POCPOU.bray"= ADONIS2.GroupedP.All.POCPOU.bray,
# "ADONIS2.P.Adj.All.POCPOU.bray"= ADONIS2.GroupedP.Adjusted.All.POCPOU.bray,
# "Betadis.P.All.POCPOU.bray"= Betadis.GroupedP.All.POCPOU.bray,
# "ANOSIM.R.All.POCPOU.bray"=  ANOSIM.GroupedR.All.POCPOU.bray,
# "ADONIS2.R.All.POCPOU.bray"= ADONIS2.GroupedR.All.POCPOU.bray)


All.POC.Stats<-data.frame(
  "ANOSIM.P.All.POC.bray"= ANOSIM.GroupedP.All.POC.bray,
  "ANOSIM.P.Adjusted.All.POC.bray"= ANOSIM.GroupedP.Adjusted.All.POC.bray,
  "ADONIS2.P.All.POC.bray"= ADONIS2.GroupedP.All.POC.bray,
  "ADONIS2.P.Adj.All.POC.bray"= ADONIS2.GroupedP.Adjusted.All.POC.bray,
  "Betadis.P.All.POC.bray"= Betadis.GroupedP.All.POC.bray,
  "Betadis.P.Adj.All.POC.bray"= Betadis.GroupedP.Adjusted.All.POC.bray,
  "ANOSIM.R.All.POC.bray"=  ANOSIM.GroupedR.All.POC.bray,
  "ADONIS2.R.All.POC.bray"= ADONIS2.GroupedR.All.POC.bray)

All.POU.Stats<-data.frame(
  "ANOSIM.P.All.POU.bray"= ANOSIM.GroupedP.All.POU.bray,
  "ANOSIM.P.Adj.All.POU.bray"= ANOSIM.GroupedP.Adjusted.All.POU.bray,
  "ADONIS2.P.All.POU.bray"= ADONIS2.GroupedP.All.POU.bray,
  "ADONIS2.P.Adj.All.POU.bray"= ADONIS2.GroupedP.Adjusted.All.POU.bray,
  "Betadis.P.All.POU.bray"= Betadis.GroupedP.All.POU.bray,
  "Betadis.P.Adj.All.POU.bray"= Betadis.GroupedP.Adjusted.All.POU.bray,
  "ANOSIM.R.All.POU.bray"=  ANOSIM.GroupedR.All.POU.bray,
  "ADONIS2.R.All.POU.bray"= ADONIS2.GroupedR.All.POU.bray)

write.xlsx(All.POCPOU.Stats, file = "C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Stats/Beta/All.POCPOU.List_special.xlsx", rowNames = TRUE)
write.xlsx(All.POC.Stats, file = "C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Stats/Beta/All.POC.List_special.xlsx", rowNames = TRUE)
write.xlsx(All.POU.Stats, file = "C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Stats/Beta/All.POU.List_special.xlsx", rowNames = TRUE)

