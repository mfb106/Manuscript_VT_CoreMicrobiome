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
            axis.title = element_text(face = "bold",size = rel(1)),
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
            legend.key.size= unit(.1, "cm"),
            legend.text = element_text(size=rel(1)),
            #legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic",size=rel(1)),
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

#resources to solve problems
# https://github.com/joey711/phyloseq/issues/936
# https://rdrr.io/github/vallenderlab/MicrobiomeR/man/root_phyloseq_tree.html
# https://forum.qiime2.org/t/qiime2r-trouble-with-phyloseq-and-rooted-tree/12052
#https://github.com/joey711/phyloseq/issues/936 


bray.weight.unrare.NMDS.PotableReusePOCPOU.ordu <- ordinate(physeq.combined.All.PotableReusePOCPOU, method="NMDS", distance="bray", weighted=TRUE, k=2, maxit = 5000, trymax = 2000 ,sfgrmin=1e-9, sratmax=0.999999999) #previous.best
bray.weight.unrare.NMDS.PotableReusePOCPOU.ordu

##Climate 
PotableReusePOCPOU.climate.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.PotableReusePOCPOU, 
                                                   bray.weight.unrare.NMDS.PotableReusePOCPOU.ordu, 
                                                   type="sites", 
                                                   color="Climate",
                                                   shape="Class1_Pot_NonPot") + 
  geom_point(size=3) +
  theme(text = element_text(size=25, face = "bold"), 
        axis.text.x = element_text( hjust = 1)) +
  #ggtitle("Unweighted UniFrac Potable Source Water")+
  theme_ipsum()+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  labs(shape="Potability",col="Climate")
PotableReusePOCPOU.climate.bray.unrare.NMDS

##Region 
PotableReusePOCPOU.region.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.PotableReusePOCPOU, 
                                                 bray.weight.unrare.NMDS.PotableReusePOCPOU.ordu, 
                                                 type="sites", 
                                                 color="Region",
                                                 shape="Class1_Pot_NonPot") + 
  geom_point(size=3) +
  theme(text = element_text(size=25, face = "bold"), 
        axis.text.x = element_text( hjust = 1)) +
  #ggtitle("Unweighted UniFrac Potable Source Water")+
  theme_ipsum()+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  labs(shape="Potability",col="Region")
PotableReusePOCPOU.region.bray.unrare.NMDS

##Final Disinfection Residual 
PotableReusePOCPOU.disinfectionresidual.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.PotableReusePOCPOU, 
                                                bray.weight.unrare.NMDS.PotableReusePOCPOU.ordu, 
                                                type="sites", 
                                                color="Final_Disinfection_Residual",
                                                shape="Class1_Pot_NonPot") + 
  geom_point(size=3) +
  theme(text = element_text(size=25, face = "bold"), 
        axis.text.x = element_text( hjust = 1)) +
  #ggtitle("Unweighted UniFrac Potable Source Water")+
  theme_ipsum()+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  labs(shape="Potability",col="Disinfection Residual")
PotableReusePOCPOU.disinfectionresidual.bray.unrare.NMDS

##Sample Local 
PotableReusePOCPOU.samplelocal.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.PotableReusePOCPOU, 
                                                     bray.weight.unrare.NMDS.PotableReusePOCPOU.ordu, 
                                                     type="sites", 
                                                     color="Sample_Local",
                                                     shape="Class1_Pot_NonPot") + 
  geom_point(size=3) +
  theme(text = element_text(size=25, face = "bold"), 
        axis.text.x = element_text( hjust = 1)) +
  #ggtitle("Unweighted UniFrac Potable Source Water")+
  theme_ipsum()+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  labs(shape="Potability",col="Sample Type")
PotableReusePOCPOU.samplelocal.bray.unrare.NMDS

##Reuse.Reclaimed_Disinfection
PotableReusePOCPOU.Reuse.Reclaimed_Disinfection.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.PotableReusePOCPOU, 
                                                     bray.weight.unrare.NMDS.PotableReusePOCPOU.ordu, 
                                                     type="sites", 
                                                     color="Reuse.Reclaimed_Disinfection",
                                                     shape="Sample_Local") + 
  geom_point(size=3) +
  theme(text = element_text(size=25, face = "bold"), 
        axis.text.x = element_text( hjust = 1)) +
  #ggtitle("Unweighted UniFrac Potable Source Water")+
  theme_ipsum()+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  labs(shape="Sample Type",col="Disinfectant")
PotableReusePOCPOU.Reuse.Reclaimed_Disinfection.bray.unrare.NMDS


##Reuse.Reclaimed_Treatment 
PotableReusePOCPOU.Reuse.Reclaimed_Treatment.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.PotableReusePOCPOU, 
                                                     bray.weight.unrare.NMDS.PotableReusePOCPOU.ordu, 
                                                     type="sites", 
                                                     color="Reuse.Reclaimed_Treatment",
                                                     shape="Sample_Local") + 
  geom_point(size=3) +
  theme(text = element_text(size=25, face = "bold"), 
        axis.text.x = element_text( hjust = 1)) +
  #ggtitle("Unweighted UniFrac Potable Source Water")+
  theme_ipsum()+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  labs(shape="Sample Type",col="Treatment Type")
PotableReusePOCPOU.Reuse.Reclaimed_Treatment.bray.unrare.NMDS


##Reuse.Reclaimed_Blending 
PotableReusePOCPOU.Reuse.Reclaimed_Blending.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.PotableReusePOCPOU, 
                                                              bray.weight.unrare.NMDS.PotableReusePOCPOU.ordu, 
                                                              type="sites", 
                                                              color="Reuse.Reclaimed_Blending",
                                                              shape="Sample_Local") + 
  geom_point(size=3) +
  theme(text = element_text(size=25, face = "bold"), 
        axis.text.x = element_text( hjust = 1)) +
  #ggtitle("Unweighted UniFrac Potable Source Water")+
  theme_ipsum()+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  labs(shape="Sample Type",col="Blending")
PotableReusePOCPOU.Reuse.Reclaimed_Blending.bray.unrare.NMDS

##PostBlend_Disinfection
PotableReusePOCPOU.PostBlend_Disinfection.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.PotableReusePOCPOU, 
                                                                              bray.weight.unrare.NMDS.PotableReusePOCPOU.ordu, 
                                                                              type="sites", 
                                                                              color="PostBlend_Disinfection",
                                                                              shape="Sample_Local") + 
  geom_point(size=3) +
  theme(text = element_text(size=25, face = "bold"), 
        axis.text.x = element_text( hjust = 1)) +
  #ggtitle("Unweighted UniFrac Potable Source Water")+
  theme_ipsum()+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  labs(shape="Sample Type",col="Post Blend Disinfection")
PotableReusePOCPOU.PostBlend_Disinfection.bray.unrare.NMDS

##PostBlend_Treatment
PotableReusePOCPOU.PostBlend_Treatment.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.PotableReusePOCPOU, 
                                                                            bray.weight.unrare.NMDS.PotableReusePOCPOU.ordu, 
                                                                            type="sites", 
                                                                            color="PostBlend_Treatment",
                                                                            shape="Sample_Local") + 
  geom_point(size=3) +
  theme(text = element_text(size=25, face = "bold"), 
        axis.text.x = element_text( hjust = 1)) +
  #ggtitle("Unweighted UniFrac Potable Source Water")+
  theme_ipsum()+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  labs(shape="Sample Type",col="Post Blend Treatment")
PotableReusePOCPOU.PostBlend_Treatment.bray.unrare.NMDS

##Disinfection_Total 
PotableReusePOCPOU.Disinfection_Total.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.PotableReusePOCPOU, 
                                                                bray.weight.unrare.NMDS.PotableReusePOCPOU.ordu, 
                                                                type="sites", 
                                                                color="Disinfection_Total",
                                                                shape="Sample_Local") + 
  geom_point(size=3) +
  theme(text = element_text(size=25, face = "bold"), 
        axis.text.x = element_text( hjust = 1)) +
  #ggtitle("Unweighted UniFrac Potable Source Water")+
  theme_ipsum()+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  labs(shape="Sample Type",col="Total Disinfection")+
  theme(legend.position='bottom')+
  guides(color = guide_legend(nrow = 11))
PotableReusePOCPOU.Disinfection_Total.bray.unrare.NMDS


##Treatment_Total 
PotableReusePOCPOU.Treatment_Total.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.PotableReusePOCPOU, 
                                                            bray.weight.unrare.NMDS.PotableReusePOCPOU.ordu, 
                                                            type="sites", 
                                                            color="Treatment_Total",
                                                            shape="Sample_Local") + 
  geom_point(size=3) +
  theme(text = element_text(size=25, face = "bold"), 
        axis.text.x = element_text( hjust = 1)) +
  #ggtitle("Unweighted UniFrac Potable Source Water")+
  theme_ipsum()+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  labs(shape="Sample Type",col="Total Treatment")+
  theme(legend.position='bottom')+
  guides(color = guide_legend(nrow = 11))
PotableReusePOCPOU.Treatment_Total.bray.unrare.NMDS

##Treatment_Classification_Number 
PotableReusePOCPOU.Treatment_Classification_Number.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.PotableReusePOCPOU, 
                                                                                bray.weight.unrare.NMDS.PotableReusePOCPOU.ordu, 
                                                                                type="sites", 
                                                                                color="Treatment_Classification_Number",
                                                                                shape="Sample_Local") + 
  geom_point(size=3) +
  theme(text = element_text(size=25, face = "bold"), 
        axis.text.x = element_text( hjust = 1)) +
  #ggtitle("Unweighted UniFrac Potable Source Water")+
  theme_ipsum()+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  labs(shape="Sample Type",col="Number of Treatment Proccessess")
PotableReusePOCPOU.Treatment_Classification_Number.bray.unrare.NMDS

##Treatment_Classification_Cat 
PotableReusePOCPOU.Treatment_Classification_Cat.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.PotableReusePOCPOU, 
                                                                             bray.weight.unrare.NMDS.PotableReusePOCPOU.ordu, 
                                                                             type="sites", 
                                                                             color="Treatment_Classification_Cat",
                                                                             shape="Sample_Local") + 
  geom_point(size=3) +
  theme(text = element_text(size=25, face = "bold"), 
        axis.text.x = element_text( hjust = 1)) +
  #ggtitle("Unweighted UniFrac Potable Source Water")+
  theme_ipsum()+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  labs(shape="Sample Type",col="Classified Treatment")
PotableReusePOCPOU.Treatment_Classification_Cat.bray.unrare.NMDS

##Mixed_Source_Water 
PotableReusePOCPOU.Mixed_Source_Water.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.PotableReusePOCPOU, 
                                                                   bray.weight.unrare.NMDS.PotableReusePOCPOU.ordu, 
                                                                   type="sites", 
                                                                   color="Mixed_Source_Water",
                                                                   shape="Sample_Local") + 
  geom_point(size=3) +
  theme(text = element_text(size=25, face = "bold"), 
        axis.text.x = element_text( hjust = 1)) +
  #ggtitle("Unweighted UniFrac Potable Source Water")+
  theme_ipsum()+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  labs(shape="Sample Type",col="Mixed Source Water")
PotableReusePOCPOU.Mixed_Source_Water.bray.unrare.NMDS

##Potable_Overlap 
PotableReusePOCPOU.Potable_Overlap.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.PotableReusePOCPOU, 
                                                                bray.weight.unrare.NMDS.PotableReusePOCPOU.ordu, 
                                                                type="sites", 
                                                                color="Potable_Overlap",
                                                                shape="Sample_Local") + 
  geom_point(size=3) +
  theme(text = element_text(size=25, face = "bold"), 
        axis.text.x = element_text( hjust = 1)) +
  #ggtitle("Unweighted UniFrac Potable Source Water")+
  theme_ipsum()+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  labs(shape="Sample Type",col="Impacted Potable Water")+ 
  guides(color = guide_legend(nrow = 2))
PotableReusePOCPOU.Potable_Overlap.bray.unrare.NMDS

##Includes_Ozone 
PotableReusePOCPOU.Includes_Ozone.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.PotableReusePOCPOU, 
                                                               bray.weight.unrare.NMDS.PotableReusePOCPOU.ordu, 
                                                               type="sites", 
                                                               color="Includes_Ozone",
                                                               shape="Sample_Local") + 
  geom_point(size=3) +
  theme(text = element_text(size=25, face = "bold"), 
        axis.text.x = element_text( hjust = 1)) +
  #ggtitle("Unweighted UniFrac Potable Source Water")+
  theme_ipsum()+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  labs(shape="Sample Type",col="Includes Ozonation")
PotableReusePOCPOU.Includes_Ozone.bray.unrare.NMDS

##Includes_RO_UF 
PotableReusePOCPOU.Includes_RO_UF.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.PotableReusePOCPOU, 
                                                               bray.weight.unrare.NMDS.PotableReusePOCPOU.ordu, 
                                                               type="sites", 
                                                               color="Includes_RO_UF",
                                                               shape="Sample_Local") + 
  geom_point(size=3) +
  theme(text = element_text(size=25, face = "bold"), 
        axis.text.x = element_text( hjust = 1)) +
  #ggtitle("Unweighted UniFrac Potable Source Water")+
  theme_ipsum()+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  labs(shape="Sample Type",col="Includes UF and/or RO")
PotableReusePOCPOU.Includes_RO_UF.bray.unrare.NMDS

##Includes_UV 
PotableReusePOCPOU.Includes_UV.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.PotableReusePOCPOU, 
                                                            bray.weight.unrare.NMDS.PotableReusePOCPOU.ordu, 
                                                            type="sites", 
                                                            color="Includes_UV",
                                                            shape="Sample_Local") + 
  geom_point(size=3) +
  theme(text = element_text(size=25, face = "bold"), 
        axis.text.x = element_text( hjust = 1)) +
  #ggtitle("Unweighted UniFrac Potable Source Water")+
  theme_ipsum()+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  labs(shape="Sample Type",col="Includes UV")
PotableReusePOCPOU.Includes_UV.bray.unrare.NMDS

##Type_DS 
PotableReusePOCPOU.Type_DS.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.PotableReusePOCPOU, 
                                                        bray.weight.unrare.NMDS.PotableReusePOCPOU.ordu, 
                                                        type="sites", 
                                                        color="Type_DS",
                                                        shape="Sample_Local") + 
  geom_point(size=3) +
  theme(text = element_text(size=25, face = "bold"), 
        axis.text.x = element_text( hjust = 1)) +
  #ggtitle("Unweighted UniFrac Potable Source Water")+
  theme_ipsum()+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  labs(shape="Sample Type",col="Type of POU")
PotableReusePOCPOU.Type_DS.bray.unrare.NMDS




##Classification_2v2 
PotableReusePOCPOU.Classification_2.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.PotableReusePOCPOU, 
                                                                 bray.weight.unrare.NMDS.PotableReusePOCPOU.ordu, 
                                                                 type="sites", 
                                                                 color="Classification_2v2",
                                                                 shape="Class1_Pot_NonPot") + 
  geom_point(size=3) +
  theme(text = element_text(size=25, face = "bold"), 
        axis.text.x = element_text( hjust = 1)) +
  #ggtitle("Unweighted UniFrac Potable Source Water")+
  theme_ipsum()+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  labs(shape="Potability",col="Classification_2")
PotableReusePOCPOU.Classification_2.bray.unrare.NMDS
# Classification_2

##Export
ggsave("PotableReusePOCPOU.wateruse.bray.unrare.NMDS.jpeg", 
       plot = PotableReusePOCPOU.wateruse.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableReuse/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("PotableReusePOCPOU.potnotpot.bray.unrare.NMDS.jpeg", 
       plot = PotableReusePOCPOU.potnotpot.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableReuse/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("PotableReusePOCPOU.climate.bray.unrare.NMDS.jpeg", 
       plot = PotableReusePOCPOU.climate.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableReuse/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("PotableReusePOCPOU.region.bray.unrare.NMDS.jpeg", 
       plot = PotableReusePOCPOU.region.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableReuse/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("PotableReusePOCPOU.disinfectionresidual.bray.unrare.NMDS.jpeg", 
       plot = PotableReusePOCPOU.disinfectionresidual.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableReuse/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("PotableReusePOCPOU.samplelocal.bray.unrare.NMDS.jpeg", 
       plot = PotableReusePOCPOU.samplelocal.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableReuse/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("PotableReusePOCPOU.Potable_Source_Water.bray.unrare.NMDS.jpeg", 
       plot = PotableReusePOCPOU.Potable_Source_Water.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableReuse/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("PotableReusePOCPOU.Potable_Disinfection_Primary.bray.unrare.NMDS.jpeg", 
       plot = PotableReusePOCPOU.Potable_Disinfection_Primary.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableReuse/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("PotableReusePOCPOU.Potable_Treatment_Type.bray.unrare.NMDS.jpeg", 
       plot = PotableReusePOCPOU.Potable_Treatment_Type.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableReuse/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("PotableReusePOCPOU.Disinfection_Total.bray.unrare.NMDS.jpeg", 
       plot = PotableReusePOCPOU.Disinfection_Total.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableReuse/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("PotableReusePOCPOU.Treatment_Total.bray.unrare.NMDS.jpeg", 
       plot = PotableReusePOCPOU.Treatment_Total.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableReuse/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("PotableReusePOCPOU.Treatment_Classification_Number.bray.unrare.NMDS.jpeg", 
       plot = PotableReusePOCPOU.Treatment_Classification_Number.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableReuse/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("PotableReusePOCPOU.Treatment_Classification_Cat.bray.unrare.NMDS.jpeg", 
       plot = PotableReusePOCPOU.Treatment_Classification_Cat.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableReuse/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("PotableReusePOCPOU.Mixed_Source_Water.bray.unrare.NMDS.jpeg", 
       plot = PotableReusePOCPOU.Mixed_Source_Water.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableReuse/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("PotableReusePOCPOU.Potable_Overlap.bray.unrare.NMDS.jpeg", 
       plot = PotableReusePOCPOU.Potable_Overlap.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableReuse/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("PotableReusePOCPOU.Includes_Ozone.bray.unrare.NMDS.jpeg", 
       plot = PotableReusePOCPOU.Includes_Ozone.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableReuse/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("PotableReusePOCPOU.Includes_RO_UF.bray.unrare.NMDS.jpeg", 
       plot = PotableReusePOCPOU.Includes_RO_UF.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableReuse/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("PotableReusePOCPOU.Includes_UV.bray.unrare.NMDS.jpeg", 
       plot = PotableReusePOCPOU.Includes_UV.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableReuse/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("PotableReusePOCPOU.Type_DS.bray.unrare.NMDS.jpeg", 
       plot = PotableReusePOCPOU.Type_DS.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableReuse/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)






ggsave("PotableReusePOCPOU.Classification_1_mod.bray.unrare.NMDS.jpeg", 
       plot = PotableReusePOCPOU.Classification_1_mod.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableReuse/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("PotableReusePOCPOU.Classification_2.bray.unrare.NMDS.jpeg", 
       plot = PotableReusePOCPOU.Classification_2.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableReuse/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("PotableReusePOCPOU.Climate_Region.bray.unrare.NMDS.jpeg", 
       plot = PotableReusePOCPOU.Climate_Region.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableReuse/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("PotableReusePOCPOU.Primers.bray.unrare.NMDS.jpeg", 
       plot = PotableReusePOCPOU.Primers.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableReuse/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("PotableReusePOCPOU.Trimmed.bray.unrare.NMDS.jpeg", 
       plot = PotableReusePOCPOU.Trimmed.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableReuse/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("PotableReusePOCPOU.Extraction_Type.bray.unrare.NMDS.jpeg", 
       plot = PotableReusePOCPOU.Extraction_Type.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableReuse/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
