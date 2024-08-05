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
            axis.text = element_text(size=rel(1)), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.box="vertical", 
            legend.margin=margin(),
            legend.direction = "horizontal",
            legend.key.size= unit(0.1, "cm"),
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

#resources to solve problems
# https://github.com/joey711/phyloseq/issues/936
# https://rdrr.io/github/vallenderlab/MicrobiomeR/man/root_phyloseq_tree.html
# https://forum.qiime2.org/t/qiime2r-trouble-with-phyloseq-and-rooted-tree/12052
#https://github.com/joey711/phyloseq/issues/936 


bray.weight.unrare.NMDS.PotableBothPOCPOU.ordu <- ordinate(physeq.combined.All.PotableBothPOCPOU, method="NMDS", distance="bray", weighted=TRUE, k=3, maxit = 5000, trymax = 2000 ,sfgrmin=1e-9, sratmax=0.999999999) #previous.best
bray.weight.unrare.NMDS.PotableBothPOCPOU.ordu

PotableBothPOCPOU.wateruse.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.PotableBothPOCPOU, 
                                                  bray.weight.unrare.NMDS.PotableBothPOCPOU.ordu, 
                                                  type="sites", 
                                                  color="Classification_1",
                                                  shape="Sample_Local") + 
  geom_point(size=3) +
  theme(text = element_text(size=25, face = "bold"), 
        axis.text.x = element_text( hjust = 1)) +
  #ggtitle("Unweighted UniFrac Potable Source Water")+
  theme_ipsum()+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  labs(shape="Sample Type",col="Intended Water Use")
PotableBothPOCPOU.wateruse.bray.unrare.NMDS



##Climate 
PotableBothPOCPOU.climate.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.PotableBothPOCPOU, 
                                                   bray.weight.unrare.NMDS.PotableBothPOCPOU.ordu, 
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
PotableBothPOCPOU.climate.bray.unrare.NMDS

##Region 
PotableBothPOCPOU.region.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.PotableBothPOCPOU, 
                                                 bray.weight.unrare.NMDS.PotableBothPOCPOU.ordu, 
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
PotableBothPOCPOU.region.bray.unrare.NMDS

##Final Disinfection Original 
PotableBothPOCPOU.disinfectionresidual.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.PotableBothPOCPOU, 
                                                bray.weight.unrare.NMDS.PotableBothPOCPOU.ordu, 
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
PotableBothPOCPOU.disinfectionresidual.bray.unrare.NMDS

##Sample Local 
PotableBothPOCPOU.samplelocal.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.PotableBothPOCPOU, 
                                                     bray.weight.unrare.NMDS.PotableBothPOCPOU.ordu, 
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
PotableBothPOCPOU.samplelocal.bray.unrare.NMDS

##Potable_Source_Water
PotableBothPOCPOU.Potable_Source_Water.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.PotableBothPOCPOU, 
                                                     bray.weight.unrare.NMDS.PotableBothPOCPOU.ordu, 
                                                     type="sites", 
                                                     color="Potable_Source_Water",
                                                     shape="Sample_Local") + 
  geom_point(size=3) +
  theme(text = element_text(size=25, face = "bold"), 
        axis.text.x = element_text( hjust = 1)) +
  #ggtitle("Unweighted UniFrac Potable Source Water")+
  theme_ipsum()+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  labs(shape="Sample Type",col="Treatment Plant Source Water")+ 
  guides(color = guide_legend(nrow = 3))
PotableBothPOCPOU.Potable_Source_Water.bray.unrare.NMDS


##Potable_Disinfection_Primary 
PotableBothPOCPOU.Potable_Disinfection_Primary.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.PotableBothPOCPOU, 
                                                     bray.weight.unrare.NMDS.PotableBothPOCPOU.ordu, 
                                                     type="sites", 
                                                     color="Potable_Disinfection_Primary",
                                                     shape="Sample_Local") + 
  geom_point(size=3) +
  theme(text = element_text(size=25, face = "bold"), 
        axis.text.x = element_text( hjust = 1)) +
  #ggtitle("Unweighted UniFrac Potable Source Water")+
  theme_ipsum()+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  labs(shape="Sample Type",col="Primary Disinfectant")+ 
  guides(color = guide_legend(nrow = 5))
PotableBothPOCPOU.Potable_Disinfection_Primary.bray.unrare.NMDS


##Potable_Treatment_Type 
PotableBothPOCPOU.Potable_Treatment_Type.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.PotableBothPOCPOU, 
                                                              bray.weight.unrare.NMDS.PotableBothPOCPOU.ordu, 
                                                              type="sites", 
                                                              color="Potable_Treatment_Type",
                                                              shape="Sample_Local") + 
  geom_point(size=3) +
  theme(text = element_text(size=25, face = "bold"), 
        axis.text.x = element_text( hjust = 1)) +
  #ggtitle("Unweighted UniFrac Potable Source Water")+
  theme_ipsum()+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  labs(shape="Sample Type",col="Treatment Type")+ 
  guides(color = guide_legend(nrow = 13))
PotableBothPOCPOU.Potable_Treatment_Type.bray.unrare.NMDS


##Disinfection_Total 
PotableBothPOCPOU.Disinfection_Total.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.PotableBothPOCPOU, 
                                                                bray.weight.unrare.NMDS.PotableBothPOCPOU.ordu, 
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
  guides(color = guide_legend(nrow = 12))
PotableBothPOCPOU.Disinfection_Total.bray.unrare.NMDS


##Treatment_Total 
PotableBothPOCPOU.Treatment_Total.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.PotableBothPOCPOU, 
                                                            bray.weight.unrare.NMDS.PotableBothPOCPOU.ordu, 
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
  guides(color = guide_legend(nrow = 12))
PotableBothPOCPOU.Treatment_Total.bray.unrare.NMDS





##Treatment_Classification_Number 
PotableBothPOCPOU.Treatment_Classification_Number.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.PotableBothPOCPOU, 
                                                                bray.weight.unrare.NMDS.PotableBothPOCPOU.ordu, 
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
PotableBothPOCPOU.Treatment_Classification_Number.bray.unrare.NMDS

##Treatment_Classification_Cat 
PotableBothPOCPOU.Treatment_Classification_Cat.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.PotableBothPOCPOU, 
                                                                bray.weight.unrare.NMDS.PotableBothPOCPOU.ordu, 
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
  labs(shape="Sample Type",col="Classified Treatment")+ 
  guides(color = guide_legend(nrow = 4))
PotableBothPOCPOU.Treatment_Classification_Cat.bray.unrare.NMDS

##Mixed_Source_Water 
PotableBothPOCPOU.Mixed_Source_Water.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.PotableBothPOCPOU, 
                                                                bray.weight.unrare.NMDS.PotableBothPOCPOU.ordu, 
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
PotableBothPOCPOU.Mixed_Source_Water.bray.unrare.NMDS

##Potable_Overlap 
PotableBothPOCPOU.Potable_Overlap.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.PotableBothPOCPOU, 
                                                                bray.weight.unrare.NMDS.PotableBothPOCPOU.ordu, 
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
  guides(color = guide_legend(nrow = 6))
PotableBothPOCPOU.Potable_Overlap.bray.unrare.NMDS

##Includes_Ozone 
PotableBothPOCPOU.Includes_Ozone.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.PotableBothPOCPOU, 
                                                                bray.weight.unrare.NMDS.PotableBothPOCPOU.ordu, 
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
PotableBothPOCPOU.Includes_Ozone.bray.unrare.NMDS

##Includes_RO_UF 
PotableBothPOCPOU.Includes_RO_UF.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.PotableBothPOCPOU, 
                                                                bray.weight.unrare.NMDS.PotableBothPOCPOU.ordu, 
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
PotableBothPOCPOU.Includes_RO_UF.bray.unrare.NMDS

##Includes_UV 
PotableBothPOCPOU.Includes_UV.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.PotableBothPOCPOU, 
                                                                bray.weight.unrare.NMDS.PotableBothPOCPOU.ordu, 
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
PotableBothPOCPOU.Includes_UV.bray.unrare.NMDS

##Export
ggsave("PotableBothPOCPOU.wateruse.bray.unrare.NMDS.jpeg", 
       plot = PotableBothPOCPOU.wateruse.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableBoth/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("PotableBothPOCPOU.potnotpot.bray.unrare.NMDS.jpeg", 
       plot = PotableBothPOCPOU.potnotpot.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableBoth/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("PotableBothPOCPOU.climate.bray.unrare.NMDS.jpeg", 
       plot = PotableBothPOCPOU.climate.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableBoth/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("PotableBothPOCPOU.region.bray.unrare.NMDS.jpeg", 
       plot = PotableBothPOCPOU.region.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableBoth/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("PotableBothPOCPOU.disinfectionresidual.bray.unrare.NMDS.jpeg", 
       plot = PotableBothPOCPOU.disinfectionresidual.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableBoth/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("PotableBothPOCPOU.samplelocal.bray.unrare.NMDS.jpeg", 
       plot = PotableBothPOCPOU.samplelocal.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableBoth/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("PotableBothPOCPOU.Potable_Source_Water.bray.unrare.NMDS.jpeg", 
       plot = PotableBothPOCPOU.Potable_Source_Water.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableBoth/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("PotableBothPOCPOU.Potable_Disinfection_Primary.bray.unrare.NMDS.jpeg", 
       plot = PotableBothPOCPOU.Potable_Disinfection_Primary.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableBoth/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("PotableBothPOCPOU.Potable_Treatment_Type.bray.unrare.NMDS.jpeg", 
       plot = PotableBothPOCPOU.Potable_Treatment_Type.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableBoth/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("PotableBothPOCPOU.Disinfection_Total.bray.unrare.NMDS.jpeg", 
       plot = PotableBothPOCPOU.Disinfection_Total.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableBoth/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("PotableBothPOCPOU.Treatment_Total.bray.unrare.NMDS.jpeg", 
       plot = PotableBothPOCPOU.Treatment_Total.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableBoth/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("PotableBothPOCPOU.Treatment_Classification_Number.bray.unrare.NMDS.jpeg", 
       plot = PotableBothPOCPOU.Treatment_Classification_Number.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableBoth/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("PotableBothPOCPOU.Treatment_Classification_Cat.bray.unrare.NMDS.jpeg", 
       plot = PotableBothPOCPOU.Treatment_Classification_Cat.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableBoth/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("PotableBothPOCPOU.Mixed_Source_Water.bray.unrare.NMDS.jpeg", 
       plot = PotableBothPOCPOU.Mixed_Source_Water.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableBoth/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("PotableBothPOCPOU.Potable_Overlap.bray.unrare.NMDS.jpeg", 
       plot = PotableBothPOCPOU.Potable_Overlap.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableBoth/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("PotableBothPOCPOU.Includes_Ozone.bray.unrare.NMDS.jpeg", 
       plot = PotableBothPOCPOU.Includes_Ozone.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableBoth/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("PotableBothPOCPOU.Includes_RO_UF.bray.unrare.NMDS.jpeg", 
       plot = PotableBothPOCPOU.Includes_RO_UF.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableBoth/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("PotableBothPOCPOU.Includes_UV.bray.unrare.NMDS.jpeg", 
       plot = PotableBothPOCPOU.Includes_UV.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableBoth/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("PotableBothPOCPOU.Type_DS.bray.unrare.NMDS.jpeg", 
       plot = PotableBothPOCPOU.Type_DS.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableBoth/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)


ggsave("PotableBothPOCPOU.Classification_1_mod.bray.unrare.NMDS.jpeg", 
       plot = PotableBothPOCPOU.Classification_1_mod.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableBoth/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("PotableBothPOCPOU.Classification_2.bray.unrare.NMDS.jpeg", 
       plot = PotableBothPOCPOU.Classification_2.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableBoth/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("PotableBothPOCPOU.Climate_Region.bray.unrare.NMDS.jpeg", 
       plot = PotableBothPOCPOU.Climate_Region.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableBoth/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("PotableBothPOCPOU.Primers.bray.unrare.NMDS.jpeg", 
       plot = PotableBothPOCPOU.Primers.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableBoth/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("PotableBothPOCPOU.Trimmed.bray.unrare.NMDS.jpeg", 
       plot = PotableBothPOCPOU.Trimmed.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableBoth/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("PotableBothPOCPOU.Extraction_Type.bray.unrare.NMDS.jpeg", 
       plot = PotableBothPOCPOU.Extraction_Type.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/PotableBoth/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)

