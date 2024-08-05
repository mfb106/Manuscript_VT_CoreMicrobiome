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
Metadata_POCPOU<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_DS_Special_V2', rowNames = TRUE)
#Metadata_POC<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_Final', rowNames = TRUE)
Metadata_POC<-Metadata_POCPOU[Metadata_POCPOU$Sample_Local == "POC",]
#Metadata_POU<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'DS_Final', rowNames = TRUE)
Metadata_POU<-Metadata_POCPOU[Metadata_POCPOU$Sample_Local == "POU",]

#Order Meta Data 
Metadata_POCPOU$Classification_1 = factor(Metadata_POCPOU$Classification_1, levels=c("Potable Conventional","Potable Reuse","Non-potable Reuse","Blank"))
Metadata_POC$Classification_1 = factor(Metadata_POC$Classification_1, levels=c("Potable Conventional","Potable Reuse","Non-potable Reuse","Blank"))
Metadata_POU$Classification_1 = factor(Metadata_POU$Classification_1, levels=c("Potable Conventional","Potable Reuse","Non-potable Reuse","Blank"))

#Bonus Figure
Metadata_POCPOU$Overlap_All = factor(Metadata_POCPOU$Overlap_All, levels=c("Potable Conventional","Potable Reuse","Non-potable Reuse","Blank","Potable Conventional with Limited Treatment and Mixed Source Water","Potable Conventional with Limited Treatment and No Mixed Source Water"))
#Test<-as.data.frame(Metadata_POCPOU$Overlap_All)

######All Samples, POCPOU
Meta.All.POCPOU<-Metadata_POCPOU
Meta.All.POCPOU<-Meta.All.POCPOU[Meta.All.POCPOU$Description != 'ConnerPrimer58',] #removed because of low read count
Meta.All.POCPOU<-Meta.All.POCPOU[Meta.All.POCPOU$Description != 'PR1WEEF',] #removed because of high read count 

Meta.All.POCPOU.Filter<-Meta.All.POCPOU[Meta.All.POCPOU$Classification_1 == "Potable Conventional"| Meta.All.POCPOU$Classification_1 =='Non-potable Reuse'|Meta.All.POCPOU$Classification_1 =='Potable Reuse'|Meta.All.POCPOU$Classification_1 =='',]#can add Blank
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

#resources to solve problems
# https://github.com/joey711/phyloseq/issues/936
# https://rdrr.io/github/vallenderlab/MicrobiomeR/man/root_phyloseq_tree.html
# https://forum.qiime2.org/t/qiime2r-trouble-with-phyloseq-and-rooted-tree/12052
#https://github.com/joey711/phyloseq/issues/936 


bray.weight.unrare.NMDS.POCPOU.ordu <- ordinate(physeq.combined.All.POCPOU, method="NMDS", distance="bray", weighted=TRUE, k=3, maxit = 5000, trymax = 2000 ,sfgrmin=1e-9, sratmax=0.999999999) #previous.best
bray.weight.unrare.NMDS.POCPOU.ordu

POCPOU.wateruse.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.POCPOU, 
                                                  bray.weight.unrare.NMDS.POCPOU.ordu, 
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
POCPOU.wateruse.bray.unrare.NMDS


###Overlap
POCPOU.overlapall.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.POCPOU, 
                                                  bray.weight.unrare.NMDS.POCPOU.ordu, 
                                                  type="sites", 
                                                  color="Overlap_All",
                                                  shape="Sample_Local") + 
  geom_point(size=3) +
  theme(text = element_text(size=25, face = "bold"), 
        axis.text.x = element_text( hjust = 1)) +
  #ggtitle("Unweighted UniFrac Potable Source Water")+
  theme_ipsum()+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  labs(shape="Sample Type",col="Intended Water Use and Overlap")+
  guides(color = guide_legend(nrow = 4))
POCPOU.overlapall.bray.unrare.NMDS





##Water Use, Classification Pot or Not Pot 
POCPOU.potnotpot.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.POCPOU, 
                                                  bray.weight.unrare.NMDS.POCPOU.ordu, 
                                                  type="sites", 
                                                  color="Class1_Pot_NonPot",
                                                  shape="Sample_Local") + 
  geom_point(size=3) +
  theme(text = element_text(size=25, face = "bold"), 
        axis.text.x = element_text( hjust = 1)) +
  #ggtitle("Unweighted UniFrac Potable Source Water")+
  theme_ipsum()+
  theme_Publication()+
  #scale_fill_Publication()+
  #scale_colour_Publication()+
  scale_color_manual(values=c("#fdb462","#386cb0"))+
  labs(shape="Sample Type",col="Potability")
POCPOU.potnotpot.bray.unrare.NMDS

##Climate 
POCPOU.climate.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.POCPOU, 
                                                   bray.weight.unrare.NMDS.POCPOU.ordu, 
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
POCPOU.climate.bray.unrare.NMDS

##Region 
POCPOU.region.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.POCPOU, 
                                                 bray.weight.unrare.NMDS.POCPOU.ordu, 
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
POCPOU.region.bray.unrare.NMDS

##Final Disinfection Original 
POCPOU.disinfectionresidual.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.POCPOU, 
                                                bray.weight.unrare.NMDS.POCPOU.ordu, 
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
POCPOU.disinfectionresidual.bray.unrare.NMDS

##Sample Local 
POCPOU.samplelocal.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.POCPOU, 
                                                     bray.weight.unrare.NMDS.POCPOU.ordu, 
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
POCPOU.samplelocal.bray.unrare.NMDS


##Classification_1_mod 
POCPOU.Classification_1_mod.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.POCPOU, 
                                                 bray.weight.unrare.NMDS.POCPOU.ordu, 
                                                 type="sites", 
                                                 color="Classification_1_mod",
                                                 shape="Class1_Pot_NonPot") + 
  geom_point(size=3) +
  theme(text = element_text(size=25, face = "bold"), 
        axis.text.x = element_text( hjust = 1)) +
  #ggtitle("Unweighted UniFrac Potable Source Water")+
  theme_ipsum()+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  labs(shape="Potability",col="Classification_1_mod")
POCPOU.Classification_1_mod.bray.unrare.NMDS
# Classification_1_mod

##Classification_2 
POCPOU.Classification_2.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.POCPOU, 
                                                 bray.weight.unrare.NMDS.POCPOU.ordu, 
                                                 type="sites", 
                                                 color="Classification_2",
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
POCPOU.Classification_2.bray.unrare.NMDS
# Classification_2

##Climate_Region 
POCPOU.Climate_Region.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.POCPOU, 
                                                 bray.weight.unrare.NMDS.POCPOU.ordu, 
                                                 type="sites", 
                                                 color="Climate_Region",
                                                 shape="Class1_Pot_NonPot") + 
  geom_point(size=3) +
  theme(text = element_text(size=25, face = "bold"), 
        axis.text.x = element_text( hjust = 1)) +
  #ggtitle("Unweighted UniFrac Potable Source Water")+
  theme_ipsum()+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  labs(shape="Potability",col="Climate_Region")
POCPOU.Climate_Region.bray.unrare.NMDS
# Climate_Region_Region

##Primers 
POCPOU.Primers.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.POCPOU, 
                                                 bray.weight.unrare.NMDS.POCPOU.ordu, 
                                                 type="sites", 
                                                 color="Primers",
                                                 shape="Class1_Pot_NonPot") + 
  geom_point(size=3) +
  theme(text = element_text(size=25, face = "bold"), 
        axis.text.x = element_text( hjust = 1)) +
  #ggtitle("Unweighted UniFrac Potable Source Water")+
  theme_ipsum()+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  labs(shape="Potability",col="Primers")
POCPOU.Primers.bray.unrare.NMDS
# Primers

##Trimmed 
POCPOU.Trimmed.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.POCPOU, 
                                                 bray.weight.unrare.NMDS.POCPOU.ordu, 
                                                 type="sites", 
                                                 color="Trimmed",
                                                 shape="Class1_Pot_NonPot") + 
  geom_point(size=3) +
  theme(text = element_text(size=25, face = "bold"), 
        axis.text.x = element_text( hjust = 1)) +
  #ggtitle("Unweighted UniFrac Potable Source Water")+
  theme_ipsum()+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  labs(shape="Potability",col="Trimmed")
POCPOU.Trimmed.bray.unrare.NMDS
# Trimmed

##Extraction_Type 
POCPOU.Extraction_Type.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.POCPOU, 
                                                 bray.weight.unrare.NMDS.POCPOU.ordu, 
                                                 type="sites", 
                                                 color="Extraction_Type",
                                                 shape="Class1_Pot_NonPot") + 
  geom_point(size=3) +
  theme(text = element_text(size=25, face = "bold"), 
        axis.text.x = element_text( hjust = 1)) +
  #ggtitle("Unweighted UniFrac Potable Source Water")+
  theme_ipsum()+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  labs(shape="Potability",col="Extraction_Type")
POCPOU.Extraction_Type.bray.unrare.NMDS
# Extraction_Type






##Export
ggsave("POCPOU.wateruse.bray.unrare.NMDS.jpeg", 
       plot = POCPOU.wateruse.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/All_Special/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("POCPOU.potnotpot.bray.unrare.NMDS.jpeg", 
       plot = POCPOU.potnotpot.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/All_Special/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("POCPOU.climate.bray.unrare.NMDS.jpeg", 
       plot = POCPOU.climate.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/All_Special/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("POCPOU.region.bray.unrare.NMDS.jpeg", 
       plot = POCPOU.region.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/All_Special/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("POCPOU.disinfectionresidual.bray.unrare.NMDS.jpeg", 
       plot = POCPOU.disinfectionresidual.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/All_Special/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("POCPOU.samplelocal.bray.unrare.NMDS.jpeg", 
       plot = POCPOU.samplelocal.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/All_Special/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)

ggsave("POCPOU.Classification_1_mod.bray.unrare.NMDS.jpeg", 
       plot = POCPOU.Classification_1_mod.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/All_Special/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("POCPOU.Classification_2.bray.unrare.NMDS.jpeg", 
       plot = POCPOU.Classification_2.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/All_Special/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("POCPOU.Climate_Region.bray.unrare.NMDS.jpeg", 
       plot = POCPOU.Climate_Region.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/All_Special/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("POCPOU.Primers.bray.unrare.NMDS.jpeg", 
       plot = POCPOU.Primers.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/All_Special/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("POCPOU.Trimmed.bray.unrare.NMDS.jpeg", 
       plot = POCPOU.Trimmed.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/All_Special/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("POCPOU.Extraction_Type.bray.unrare.NMDS.jpeg", 
       plot = POCPOU.Extraction_Type.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/All_Special/POCPOU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)


######All Samples, POC
Meta.All.POC<-Metadata_POC
Meta.All.POC<-Meta.All.POC[Meta.All.POC$Description != 'ConnerPrimer58',] #removed because of low read count
Meta.All.POC<-Meta.All.POC[Meta.All.POC$Description != 'PR1WEEF',] #removed because of high read count 

Meta.All.POC.Filter<-Meta.All.POC[Meta.All.POC$Classification_1 == "Potable Conventional"| Meta.All.POC$Classification_1 =='Non-potable Reuse'|Meta.All.POC$Classification_1 =='Potable Reuse'|Meta.All.POC$Classification_1 =='',]#can add Blank
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

#resources to solve problems
# https://github.com/joey711/phyloseq/issues/936
# https://rdrr.io/github/vallenderlab/MicrobiomeR/man/root_phyloseq_tree.html
# https://forum.qiime2.org/t/qiime2r-trouble-with-phyloseq-and-rooted-tree/12052
#https://github.com/joey711/phyloseq/issues/936 


bray.weight.unrare.NMDS.POC.ordu <- ordinate(physeq.combined.All.POC, method="NMDS", distance="bray", weighted=TRUE, k=3, maxit = 5000, trymax = 2000 ,sfgrmin=1e-9, sratmax=0.999999999)
bray.weight.unrare.NMDS.POC.ordu

POC.wateruse.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.POC, 
                                                  bray.weight.unrare.NMDS.POC.ordu, 
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
POC.wateruse.bray.unrare.NMDS


##Water Use, Classification Pot or Not Pot 
POC.potnotpot.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.POC, 
                                                   bray.weight.unrare.NMDS.POC.ordu, 
                                                   type="sites", 
                                                   color="Class1_Pot_NonPot",
                                                   shape="Sample_Local") + 
  geom_point(size=3) +
  theme(text = element_text(size=25, face = "bold"), 
        axis.text.x = element_text( hjust = 1)) +
  #ggtitle("Unweighted UniFrac Potable Source Water")+
  theme_ipsum()+
  theme_Publication()+
  #scale_fill_Publication()+
  #scale_colour_Publication()+
  scale_color_manual(values=c("#fdb462","#386cb0"))+
  labs(shape="Sample Type",col="Potability")
POC.potnotpot.bray.unrare.NMDS

##Climate 
POC.climate.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.POC, 
                                                 bray.weight.unrare.NMDS.POC.ordu, 
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
POC.climate.bray.unrare.NMDS

##Region 
POC.region.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.POC, 
                                                bray.weight.unrare.NMDS.POC.ordu, 
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
POC.region.bray.unrare.NMDS

##Final Disinfection Original 
POC.disinfectionresidual.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.POC, 
                                                              bray.weight.unrare.NMDS.POC.ordu, 
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
POC.disinfectionresidual.bray.unrare.NMDS

##Sample Local 
POC.samplelocal.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.POC, 
                                                     bray.weight.unrare.NMDS.POC.ordu, 
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
POC.samplelocal.bray.unrare.NMDS


##Classification_1_mod 
POC.Classification_1_mod.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.POC, 
                                                              bray.weight.unrare.NMDS.POC.ordu, 
                                                              type="sites", 
                                                              color="Classification_1_mod",
                                                              shape="Class1_Pot_NonPot") + 
  geom_point(size=3) +
  theme(text = element_text(size=25, face = "bold"), 
        axis.text.x = element_text( hjust = 1)) +
  #ggtitle("Unweighted UniFrac Potable Source Water")+
  theme_ipsum()+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  labs(shape="Potability",col="Classification_1_mod")
POC.Classification_1_mod.bray.unrare.NMDS
# Classification_1_mod

##Classification_2 
POC.Classification_2.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.POC, 
                                                          bray.weight.unrare.NMDS.POC.ordu, 
                                                          type="sites", 
                                                          color="Classification_2",
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
POC.Classification_2.bray.unrare.NMDS
# Classification_2

##Climate_Region 
POC.Climate_Region.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.POC, 
                                                        bray.weight.unrare.NMDS.POC.ordu, 
                                                        type="sites", 
                                                        color="Climate_Region",
                                                        shape="Class1_Pot_NonPot") + 
  geom_point(size=3) +
  theme(text = element_text(size=25, face = "bold"), 
        axis.text.x = element_text( hjust = 1)) +
  #ggtitle("Unweighted UniFrac Potable Source Water")+
  theme_ipsum()+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  labs(shape="Potability",col="Climate_Region")
POC.Climate_Region.bray.unrare.NMDS
# Climate_Region_Region

##Primers 
POC.Primers.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.POC, 
                                                 bray.weight.unrare.NMDS.POC.ordu, 
                                                 type="sites", 
                                                 color="Primers",
                                                 shape="Class1_Pot_NonPot") + 
  geom_point(size=3) +
  theme(text = element_text(size=25, face = "bold"), 
        axis.text.x = element_text( hjust = 1)) +
  #ggtitle("Unweighted UniFrac Potable Source Water")+
  theme_ipsum()+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  labs(shape="Potability",col="Primers")
POC.Primers.bray.unrare.NMDS
# Primers

##Trimmed 
POC.Trimmed.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.POC, 
                                                 bray.weight.unrare.NMDS.POC.ordu, 
                                                 type="sites", 
                                                 color="Trimmed",
                                                 shape="Class1_Pot_NonPot") + 
  geom_point(size=3) +
  theme(text = element_text(size=25, face = "bold"), 
        axis.text.x = element_text( hjust = 1)) +
  #ggtitle("Unweighted UniFrac Potable Source Water")+
  theme_ipsum()+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  labs(shape="Potability",col="Trimmed")
POC.Trimmed.bray.unrare.NMDS
# Trimmed

##Extraction_Type 
POC.Extraction_Type.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.POC, 
                                                         bray.weight.unrare.NMDS.POC.ordu, 
                                                         type="sites", 
                                                         color="Extraction_Type",
                                                         shape="Class1_Pot_NonPot") + 
  geom_point(size=3) +
  theme(text = element_text(size=25, face = "bold"), 
        axis.text.x = element_text( hjust = 1)) +
  #ggtitle("Unweighted UniFrac Potable Source Water")+
  theme_ipsum()+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  labs(shape="Potability",col="Extraction_Type")
POC.Extraction_Type.bray.unrare.NMDS
# Extraction_Type






##Export
ggsave("POC.wateruse.bray.unrare.NMDS.jpeg", 
       plot = POC.wateruse.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/All_Special/POC", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("POC.potnotpot.bray.unrare.NMDS.jpeg", 
       plot = POC.potnotpot.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/All_Special/POC", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("POC.climate.bray.unrare.NMDS.jpeg", 
       plot = POC.climate.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/All_Special/POC", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("POC.region.bray.unrare.NMDS.jpeg", 
       plot = POC.region.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/All_Special/POC", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("POC.disinfectionresidual.bray.unrare.NMDS.jpeg", 
       plot = POC.disinfectionresidual.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/All_Special/POC", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("POC.samplelocal.bray.unrare.NMDS.jpeg", 
       plot = POC.samplelocal.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/All_Special/POC", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)

ggsave("POC.Classification_1_mod.bray.unrare.NMDS.jpeg", 
       plot = POC.Classification_1_mod.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/All_Special/POC", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("POC.Classification_2.bray.unrare.NMDS.jpeg", 
       plot = POC.Classification_2.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/All_Special/POC", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("POC.Climate_Region.bray.unrare.NMDS.jpeg", 
       plot = POC.Climate_Region.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/All_Special/POC", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("POC.Primers.bray.unrare.NMDS.jpeg", 
       plot = POC.Primers.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/All_Special/POC", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("POC.Trimmed.bray.unrare.NMDS.jpeg", 
       plot = POC.Trimmed.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/All_Special/POC", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("POC.Extraction_Type.bray.unrare.NMDS.jpeg", 
       plot = POC.Extraction_Type.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/All_Special/POC", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)


######All Samples, POU
Meta.All.POU<-Metadata_POU
Meta.All.POU<-Meta.All.POU[Meta.All.POU$Description != 'ConnerPrimer58',] #removed because of low read count
Meta.All.POU<-Meta.All.POU[Meta.All.POU$Description != 'PR1WEEF',] #removed because of high read count 

Meta.All.POU.Filter<-Meta.All.POU[Meta.All.POU$Classification_1 == "Potable Conventional"| Meta.All.POU$Classification_1 =='Non-potable Reuse'|Meta.All.POU$Classification_1 =='Potable Reuse'|Meta.All.POU$Classification_1 =='',]#can add Blank
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

#resources to solve problems
# https://github.com/joey711/phyloseq/issues/936
# https://rdrr.io/github/vallenderlab/MicrobiomeR/man/root_phyloseq_tree.html
# https://forum.qiime2.org/t/qiime2r-trouble-with-phyloseq-and-rooted-tree/12052
#https://github.com/joey711/phyloseq/issues/936 


bray.weight.unrare.NMDS.POU.ordu <- ordinate(physeq.combined.All.POU, method="NMDS", distance="bray", weighted=TRUE, k=3, maxit = 5000, trymax = 2000 ,sfgrmin=1e-9, sratmax=0.999999999)
bray.weight.unrare.NMDS.POU.ordu

POU.wateruse.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.POU, 
                                                  bray.weight.unrare.NMDS.POU.ordu, 
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
POU.wateruse.bray.unrare.NMDS


##Water Use, Classification Pot or Not Pot 
POU.potnotpot.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.POU, 
                                                   bray.weight.unrare.NMDS.POU.ordu, 
                                                   type="sites", 
                                                   color="Class1_Pot_NonPot",
                                                   shape="Sample_Local") + 
  geom_point(size=3) +
  theme(text = element_text(size=25, face = "bold"), 
        axis.text.x = element_text( hjust = 1)) +
  #ggtitle("Unweighted UniFrac Potable Source Water")+
  theme_ipsum()+
  theme_Publication()+
  #scale_fill_Publication()+
  #scale_colour_Publication()+
  scale_color_manual(values=c("#fdb462","#386cb0"))+
  labs(shape="Sample Type",col="Potability")
POU.potnotpot.bray.unrare.NMDS

##Climate 
POU.climate.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.POU, 
                                                 bray.weight.unrare.NMDS.POU.ordu, 
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
POU.climate.bray.unrare.NMDS

##Region 
POU.region.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.POU, 
                                                bray.weight.unrare.NMDS.POU.ordu, 
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
POU.region.bray.unrare.NMDS

##Final Disinfection Original 
POU.disinfectionresidual.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.POU, 
                                                              bray.weight.unrare.NMDS.POU.ordu, 
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
POU.disinfectionresidual.bray.unrare.NMDS

##Sample Local 
POU.samplelocal.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.POU, 
                                                     bray.weight.unrare.NMDS.POU.ordu, 
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
POU.samplelocal.bray.unrare.NMDS


##Classification_1_mod 
POU.Classification_1_mod.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.POU, 
                                                              bray.weight.unrare.NMDS.POU.ordu, 
                                                              type="sites", 
                                                              color="Classification_1_mod",
                                                              shape="Class1_Pot_NonPot") + 
  geom_point(size=3) +
  theme(text = element_text(size=25, face = "bold"), 
        axis.text.x = element_text( hjust = 1)) +
  #ggtitle("Unweighted UniFrac Potable Source Water")+
  theme_ipsum()+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  labs(shape="Potability",col="Classification_1_mod")
POU.Classification_1_mod.bray.unrare.NMDS
# Classification_1_mod

##Classification_2 
POU.Classification_2.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.POU, 
                                                          bray.weight.unrare.NMDS.POU.ordu, 
                                                          type="sites", 
                                                          color="Classification_2",
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
POU.Classification_2.bray.unrare.NMDS
# Classification_2

##Climate_Region 
POU.Climate_Region.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.POU, 
                                                        bray.weight.unrare.NMDS.POU.ordu, 
                                                        type="sites", 
                                                        color="Climate_Region",
                                                        shape="Class1_Pot_NonPot") + 
  geom_point(size=3) +
  theme(text = element_text(size=25, face = "bold"), 
        axis.text.x = element_text( hjust = 1)) +
  #ggtitle("Unweighted UniFrac Potable Source Water")+
  theme_ipsum()+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  labs(shape="Potability",col="Climate_Region")
POU.Climate_Region.bray.unrare.NMDS
# Climate_Region_Region

##Primers 
POU.Primers.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.POU, 
                                                 bray.weight.unrare.NMDS.POU.ordu, 
                                                 type="sites", 
                                                 color="Primers",
                                                 shape="Class1_Pot_NonPot") + 
  geom_point(size=3) +
  theme(text = element_text(size=25, face = "bold"), 
        axis.text.x = element_text( hjust = 1)) +
  #ggtitle("Unweighted UniFrac Potable Source Water")+
  theme_ipsum()+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  labs(shape="Potability",col="Primers")
POU.Primers.bray.unrare.NMDS
# Primers

##Trimmed 
POU.Trimmed.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.POU, 
                                                 bray.weight.unrare.NMDS.POU.ordu, 
                                                 type="sites", 
                                                 color="Trimmed",
                                                 shape="Class1_Pot_NonPot") + 
  geom_point(size=3) +
  theme(text = element_text(size=25, face = "bold"), 
        axis.text.x = element_text( hjust = 1)) +
  #ggtitle("Unweighted UniFrac Potable Source Water")+
  theme_ipsum()+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  labs(shape="Potability",col="Trimmed")
POU.Trimmed.bray.unrare.NMDS
# Trimmed

##Extraction_Type 
POU.Extraction_Type.bray.unrare.NMDS<-plot_ordination(physeq.combined.All.POU, 
                                                         bray.weight.unrare.NMDS.POU.ordu, 
                                                         type="sites", 
                                                         color="Extraction_Type",
                                                         shape="Class1_Pot_NonPot") + 
  geom_point(size=3) +
  theme(text = element_text(size=25, face = "bold"), 
        axis.text.x = element_text( hjust = 1)) +
  #ggtitle("Unweighted UniFrac Potable Source Water")+
  theme_ipsum()+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  labs(shape="Potability",col="Extraction_Type")
POU.Extraction_Type.bray.unrare.NMDS
# Extraction_Type






##Export
ggsave("POU.wateruse.bray.unrare.NMDS.jpeg", 
       plot = POU.wateruse.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/All_Special/POU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("POU.potnotpot.bray.unrare.NMDS.jpeg", 
       plot = POU.potnotpot.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/All_Special/POU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("POU.climate.bray.unrare.NMDS.jpeg", 
       plot = POU.climate.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/All_Special/POU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("POU.region.bray.unrare.NMDS.jpeg", 
       plot = POU.region.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/All_Special/POU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("POU.disinfectionresidual.bray.unrare.NMDS.jpeg", 
       plot = POU.disinfectionresidual.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/All_Special/POU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("POU.samplelocal.bray.unrare.NMDS.jpeg", 
       plot = POU.samplelocal.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/All_Special/POU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)

ggsave("POU.Classification_1_mod.bray.unrare.NMDS.jpeg", 
       plot = POU.Classification_1_mod.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/All_Special/POU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("POU.Classification_2.bray.unrare.NMDS.jpeg", 
       plot = POU.Classification_2.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/All_Special/POU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("POU.Climate_Region.bray.unrare.NMDS.jpeg", 
       plot = POU.Climate_Region.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/All_Special/POU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("POU.Primers.bray.unrare.NMDS.jpeg", 
       plot = POU.Primers.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/All_Special/POU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("POU.Trimmed.bray.unrare.NMDS.jpeg", 
       plot = POU.Trimmed.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/All_Special/POU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("POU.Extraction_Type.bray.unrare.NMDS.jpeg", 
       plot = POU.Extraction_Type.bray.unrare.NMDS,
       path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Ordination/CA_Figures/UnRare/Bray/All_Special/POU", 
       scale = 1, width = 8.5, height = 8, units = "in", dpi = 300)







