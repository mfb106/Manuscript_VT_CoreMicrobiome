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

#rooted tree
comp_rooted_tree<-read_qza("C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/Qiime_Outputs_Rerun/Silva/Analysis/Comparison_Rerun_silva-rooted-tree.qza")

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

Meta.All.POCPOU.Filter<-Meta.All.POCPOU[Meta.All.POCPOU$Classification_1 == 'Potable'| Meta.All.POCPOU$Classification_1 =='Non-potable Reuse'|Meta.All.POCPOU$Classification_1 =='Potable Reuse'|Meta.All.POCPOU$Classification_1 =='Blank',]#can add Blank
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



##level1
physeq.All.POCPOU_level1 <- tax_glom(physeq.All.POCPOU, taxrank=rank_names(physeq.All.POCPOU)[1], NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))

#Pull out OTU Table 
OTU_level1<-as.data.frame(physeq.All.POCPOU_level1@otu_table)
OTU_level1_Trans<-as.data.frame(t(as.matrix(OTU_level1)))
#OTU_level1_Trans <- tibble::rownames_to_column(OTU_level1_Trans, "OTU")

Taxa_level1<-as.data.frame(physeq.All.POCPOU_level1@tax_table)
#Taxa_level1 <- tibble::rownames_to_column(Taxa_level1, "OTU")
Taxa_level1_Kingdom<-select(Taxa_level1, Kingdom)

level1_Collapsed_OTUs_data<-merge(OTU_level1_Trans,Taxa_level1_Kingdom, by=0)
rownames(level1_Collapsed_OTUs_data) <- level1_Collapsed_OTUs_data$Kingdom
level1_Collapsed_OTUs_data$Kingdom <- NULL
level1_Collapsed_OTUs_data$Row.names <- NULL


##level2
physeq.All.POCPOU_level2 <- tax_glom(physeq.All.POCPOU, taxrank=rank_names(physeq.All.POCPOU)[2], NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))

#Pull out OTU Table 
OTU_level2<-as.data.frame(physeq.All.POCPOU_level2@otu_table)
OTU_level2_Trans<-as.data.frame(t(as.matrix(OTU_level2)))
#OTU_level2_Trans <- tibble::rownames_to_column(OTU_level2_Trans, "OTU")

Taxa_level2<-as.data.frame(physeq.All.POCPOU_level2@tax_table)
#Taxa_level2 <- tibble::rownames_to_column(Taxa_level2, "OTU")
Taxa_level2_Phylum<-select(Taxa_level2, Phylum)

level2_Collapsed_OTUs_data<-merge(OTU_level2_Trans,Taxa_level2_Phylum, by=0)
rownames(level2_Collapsed_OTUs_data) <- level2_Collapsed_OTUs_data$Phylum
level2_Collapsed_OTUs_data$Phylum <- NULL
level2_Collapsed_OTUs_data$Row.names <- NULL

##level3
physeq.All.POCPOU_level3 <- tax_glom(physeq.All.POCPOU, taxrank=rank_names(physeq.All.POCPOU)[3], NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))

#Pull out OTU Table 
OTU_level3<-as.data.frame(physeq.All.POCPOU_level3@otu_table)
OTU_level3_Trans<-as.data.frame(t(as.matrix(OTU_level3)))
#OTU_level3_Trans <- tibble::rownames_to_column(OTU_level3_Trans, "OTU")

Taxa_level3<-as.data.frame(physeq.All.POCPOU_level3@tax_table)
#Taxa_level3 <- tibble::rownames_to_column(Taxa_level3, "OTU")
Taxa_level3_Class<-select(Taxa_level3, Class)

level3_Collapsed_OTUs_data<-merge(OTU_level3_Trans,Taxa_level3_Class, by=0)
level3_Collapsed_OTUs_data<-level3_Collapsed_OTUs_data[-(which(level3_Collapsed_OTUs_data$Class %in% "Incertae_Sedis")),]
level3_Collapsed_OTUs_data<-level3_Collapsed_OTUs_data[-(which(level3_Collapsed_OTUs_data$Class %in% "uncultured")),]

rownames(level3_Collapsed_OTUs_data) <- level3_Collapsed_OTUs_data$Class
level3_Collapsed_OTUs_data$Class <- NULL
level3_Collapsed_OTUs_data$Row.names <- NULL

##level4
physeq.All.POCPOU_level4 <- tax_glom(physeq.All.POCPOU, taxrank=rank_names(physeq.All.POCPOU)[4], NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))

#Pull out OTU Table 
OTU_level4<-as.data.frame(physeq.All.POCPOU_level4@otu_table)
OTU_level4_Trans<-as.data.frame(t(as.matrix(OTU_level4)))
#OTU_level4_Trans <- tibble::rownames_to_column(OTU_level4_Trans, "OTU")

Taxa_level4<-as.data.frame(physeq.All.POCPOU_level4@tax_table)
#Taxa_level4 <- tibble::rownames_to_column(Taxa_level4, "OTU")
Taxa_level4_Order<-select(Taxa_level4, Order)

level4_Collapsed_OTUs_data<-merge(OTU_level4_Trans,Taxa_level4_Order, by=0)
level4_Collapsed_OTUs_data<-level4_Collapsed_OTUs_data[-(which(level4_Collapsed_OTUs_data$Order %in% "Incertae_Sedis")),]
level4_Collapsed_OTUs_data<-level4_Collapsed_OTUs_data[-(which(level4_Collapsed_OTUs_data$Order %in% "uncultured")),]

rownames(level4_Collapsed_OTUs_data) <- level4_Collapsed_OTUs_data$Order
level4_Collapsed_OTUs_data$Order <- NULL
level4_Collapsed_OTUs_data$Row.names <- NULL

##level5
physeq.All.POCPOU_level5 <- tax_glom(physeq.All.POCPOU, taxrank=rank_names(physeq.All.POCPOU)[5], NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))

#Pull out OTU Table 
OTU_level5<-as.data.frame(physeq.All.POCPOU_level5@otu_table)
OTU_level5_Trans<-as.data.frame(t(as.matrix(OTU_level5)))
#OTU_level5_Trans <- tibble::rownames_to_column(OTU_level5_Trans, "OTU")

Taxa_level5<-as.data.frame(physeq.All.POCPOU_level5@tax_table)
#Taxa_level5 <- tibble::rownames_to_column(Taxa_level5, "OTU")
Taxa_level5_Family<-select(Taxa_level5, Family)

level5_Collapsed_OTUs_data<-merge(OTU_level5_Trans,Taxa_level5_Family, by=0)
level5_Collapsed_OTUs_data<-level5_Collapsed_OTUs_data[-(which(level5_Collapsed_OTUs_data$Family %in% "Incertae_Sedis")),]
level5_Collapsed_OTUs_data<-level5_Collapsed_OTUs_data[-(which(level5_Collapsed_OTUs_data$Family %in% "uncultured")),]
level5_Collapsed_OTUs_data<-level5_Collapsed_OTUs_data[-(which(level5_Collapsed_OTUs_data$Family %in% "Unknown_Family")),]

rownames(level5_Collapsed_OTUs_data) <- level5_Collapsed_OTUs_data$Family
level5_Collapsed_OTUs_data$Family <- NULL
level5_Collapsed_OTUs_data$Row.names <- NULL

##level6
physeq.All.POCPOU_level6 <- tax_glom(physeq.All.POCPOU, taxrank=rank_names(physeq.All.POCPOU)[6], NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))

#Pull out OTU Table 
OTU_level6<-as.data.frame(physeq.All.POCPOU_level6@otu_table)
OTU_level6_Trans<-as.data.frame(t(as.matrix(OTU_level6)))
#OTU_level6_Trans <- tibble::rownames_to_column(OTU_level6_Trans, "OTU")

Taxa_level6<-as.data.frame(physeq.All.POCPOU_level6@tax_table)
#Taxa_level6 <- tibble::rownames_to_column(Taxa_level6, "OTU")
Taxa_level6_Genus<-select(Taxa_level6, Genus)

level6_Collapsed_OTUs_data<-merge(OTU_level6_Trans,Taxa_level6_Genus, by=0)
level6_Collapsed_OTUs_data<-level6_Collapsed_OTUs_data[-(which(level6_Collapsed_OTUs_data$Genus %in% "Incertae_Sedis")),]
level6_Collapsed_OTUs_data<-level6_Collapsed_OTUs_data[-(which(level6_Collapsed_OTUs_data$Genus %in% "uncultured")),]
level6_Collapsed_OTUs_data<-level6_Collapsed_OTUs_data[-(which(level6_Collapsed_OTUs_data$Genus %in% "Unknown_Family")),]
rownames(level6_Collapsed_OTUs_data) <- level6_Collapsed_OTUs_data$Genus
level6_Collapsed_OTUs_data$Genus <- NULL
level6_Collapsed_OTUs_data$Row.names <- NULL


##Turn PA
level1_Collapsed_OTUs_data_PA<-as.data.frame(level1_Collapsed_OTUs_data)
level1_Collapsed_OTUs_data_PA[]<-as.integer((level1_Collapsed_OTUs_data_PA[]!=0))

level2_Collapsed_OTUs_data_PA<-as.data.frame(level2_Collapsed_OTUs_data)
level2_Collapsed_OTUs_data_PA[]<-as.integer((level2_Collapsed_OTUs_data_PA[]!=0))

level3_Collapsed_OTUs_data_PA<-as.data.frame(level3_Collapsed_OTUs_data)
level3_Collapsed_OTUs_data_PA[]<-as.integer((level3_Collapsed_OTUs_data_PA[]!=0))

level4_Collapsed_OTUs_data_PA<-as.data.frame(level4_Collapsed_OTUs_data)
level4_Collapsed_OTUs_data_PA[]<-as.integer((level4_Collapsed_OTUs_data_PA[]!=0))

level5_Collapsed_OTUs_data_PA<-as.data.frame(level5_Collapsed_OTUs_data)
level5_Collapsed_OTUs_data_PA[]<-as.integer((level5_Collapsed_OTUs_data_PA[]!=0))

level6_Collapsed_OTUs_data_PA<-as.data.frame(level6_Collapsed_OTUs_data)
level6_Collapsed_OTUs_data_PA[]<-as.integer((level6_Collapsed_OTUs_data_PA[]!=0))

#Metadata
#Read Meta Data
Metadata_POCPOU<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_DS_Special', rowNames = TRUE)

#Order Meta Data
Metadata_POCPOU$Classification_1 = factor(Metadata_POCPOU$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))

######All Samples, POCPOU
Meta.All.POCPOU<-Metadata_POCPOU
Meta.All.POCPOU<-Meta.All.POCPOU[Meta.All.POCPOU$Description != 'ConnerPrimer58',] #removed because of low read count
Meta.All.POCPOU<-Meta.All.POCPOU[Meta.All.POCPOU$Description != 'PR1WEEF',] #removed because of high read count

####Remove Samples that include limited treatment and mixed source waters 
Meta.All.POCPOU<-Meta.All.POCPOU[Meta.All.POCPOU$Potable_Overlap != 'Limited Treatment and Mixed Source Water',] #removed because of low read count


###trampose dataframe
Combined_Trans<-as.data.frame(t(as.matrix(table_all_data)))
Combined_Trans_PA<-as.data.frame(t(as.matrix(table_all_data_PA)))

level1_Collapsed_OTUs_Trans<-as.data.frame(t(as.matrix(level1_Collapsed_OTUs_data)))
level1_Collapsed_OTUs_Trans_PA<-as.data.frame(t(as.matrix(level1_Collapsed_OTUs_data_PA)))

level2_Collapsed_OTUs_Trans<-as.data.frame(t(as.matrix(level2_Collapsed_OTUs_data)))
level2_Collapsed_OTUs_Trans_PA<-as.data.frame(t(as.matrix(level2_Collapsed_OTUs_data_PA)))

level3_Collapsed_OTUs_Trans<-as.data.frame(t(as.matrix(level3_Collapsed_OTUs_data)))
level3_Collapsed_OTUs_Trans_PA<-as.data.frame(t(as.matrix(level3_Collapsed_OTUs_data_PA)))

level4_Collapsed_OTUs_Trans<-as.data.frame(t(as.matrix(level4_Collapsed_OTUs_data)))
level4_Collapsed_OTUs_Trans_PA<-as.data.frame(t(as.matrix(level4_Collapsed_OTUs_data_PA)))

level5_Collapsed_OTUs_Trans<-as.data.frame(t(as.matrix(level5_Collapsed_OTUs_data)))
level5_Collapsed_OTUs_Trans_PA<-as.data.frame(t(as.matrix(level5_Collapsed_OTUs_data_PA)))

level6_Collapsed_OTUs_Trans<-as.data.frame(t(as.matrix(level6_Collapsed_OTUs_data)))
level6_Collapsed_OTUs_Trans_PA<-as.data.frame(t(as.matrix(level6_Collapsed_OTUs_data_PA)))

###Merged tax data and sample metadata for final table 
Final_Mod_Merged<-transform(merge(Meta.All.POCPOU,Combined_Trans,by=0), row.names=Row.names, Row.names=NULL)
Final_Mod_Merged_PA<-transform(merge(Meta.All.POCPOU,Combined_Trans_PA,by=0), row.names=Row.names, Row.names=NULL)

Final_Mod_level1_Collapsed_OTUs_Trans_Merged<-transform(merge(Meta.All.POCPOU,level1_Collapsed_OTUs_Trans,by=0), row.names=Row.names, Row.names=NULL)
Final_Mod_level1_Collapsed_OTUs_Trans_Merged_PA<-transform(merge(Meta.All.POCPOU,level1_Collapsed_OTUs_Trans_PA,by=0), row.names=Row.names, Row.names=NULL)

Final_Mod_level2_Collapsed_OTUs_Trans_Merged<-transform(merge(Meta.All.POCPOU,level2_Collapsed_OTUs_Trans,by=0), row.names=Row.names, Row.names=NULL)
Final_Mod_level2_Collapsed_OTUs_Trans_Merged_PA<-transform(merge(Meta.All.POCPOU,level2_Collapsed_OTUs_Trans_PA,by=0), row.names=Row.names, Row.names=NULL)

Final_Mod_level3_Collapsed_OTUs_Trans_Merged<-transform(merge(Meta.All.POCPOU,level3_Collapsed_OTUs_Trans,by=0), row.names=Row.names, Row.names=NULL)
Final_Mod_level3_Collapsed_OTUs_Trans_Merged_PA<-transform(merge(Meta.All.POCPOU,level3_Collapsed_OTUs_Trans_PA,by=0), row.names=Row.names, Row.names=NULL)

Final_Mod_level4_Collapsed_OTUs_Trans_Merged<-transform(merge(Meta.All.POCPOU,level4_Collapsed_OTUs_Trans,by=0), row.names=Row.names, Row.names=NULL)
Final_Mod_level4_Collapsed_OTUs_Trans_Merged_PA<-transform(merge(Meta.All.POCPOU,level4_Collapsed_OTUs_Trans_PA,by=0), row.names=Row.names, Row.names=NULL)

Final_Mod_level5_Collapsed_OTUs_Trans_Merged<-transform(merge(Meta.All.POCPOU,level5_Collapsed_OTUs_Trans,by=0), row.names=Row.names, Row.names=NULL)
Final_Mod_level5_Collapsed_OTUs_Trans_Merged_PA<-transform(merge(Meta.All.POCPOU,level5_Collapsed_OTUs_Trans_PA,by=0), row.names=Row.names, Row.names=NULL)

Final_Mod_level6_Collapsed_OTUs_Trans_Merged<-transform(merge(Meta.All.POCPOU,level6_Collapsed_OTUs_Trans,by=0), row.names=Row.names, Row.names=NULL)
Final_Mod_level6_Collapsed_OTUs_Trans_Merged_PA<-transform(merge(Meta.All.POCPOU,level6_Collapsed_OTUs_Trans_PA,by=0), row.names=Row.names, Row.names=NULL)

###########Presence/Absence Based Analysis
row.names.remove <- c("SampleType",
                      "SampleEvent",
                      "Location",
                      "SourceWater",
                      "Matrix",
                      "Treatment",
                      "SourceData",
                      "Description",
                      "LocationAdjustedBetaDiversity",
                      "Include_Comparison",
                      "Paper",
                      "Classification_1",
                      "Classification_1_mod",
                      "Classification_2",
                      "Geolocation",
                      "Climate_Region",
                      "Region",
                      "Climate",
                      "Location_Type",
                      "Sample_Local",
                      "DS_Location.Time",
                      "Final_Mod_Disinfection_Residual",
                      "Potable_Source_Water",
                      "Potable_Disinfection_Primary",
                      "Potable_Treatment_Type",
                      "Potable_Scale",
                      "Reuse.Reclaimed_Disinfection",
                      "Reuse.Reclaimed_Treatment",
                      "Reuse.Reclaimed_Blending",
                      "Reuse.Reclaimed_Scale",
                      "PostBlend_Disinfection",
                      "PostBlend_Treatment",
                      "PostBlend_Scale",
                      "Final_Mod_Disinfection_Residual_Old",
                      "Class1_Pot_NonPot",
                      "Disinfection_Total",
                      "Treatment_Total",
                      "Treatment_Classification_Cat",
                      "Treatment_Classification_Number",
                      "Mixed_Source_Water",
                      "Potable_Overlap",
                      "Includes_Ozone",
                      "Includes_RO_UF",
                      "Includes_UV",
                      "Type_DS",
                      "Final_Disinfection_Residual_Old",
                      "Final_Disinfection_Residual")

#######
#######
#Filter All, WTP and DS, level2
#######
#######
PA_level2_All<-as.data.frame(Final_Mod_level2_Collapsed_OTUs_Trans_Merged_PA[Final_Mod_level2_Collapsed_OTUs_Trans_Merged_PA$Sample_Local == 'POC'| Final_Mod_level2_Collapsed_OTUs_Trans_Merged_PA$Sample_Local =='Bottled Water'| Final_Mod_level2_Collapsed_OTUs_Trans_Merged_PA$Sample_Local =='POU'| Final_Mod_level2_Collapsed_OTUs_Trans_Merged_PA$Sample_Local =='Blank',])

#Transpose, remove unwanted rows, leaving just SV  
PA_level2_All_Tran<-as.data.frame(t(as.matrix(PA_level2_All)))
PA_level2_All_Tran<-PA_level2_All_Tran[!(row.names(PA_level2_All_Tran) %in% row.names.remove), ]

Final_Mod_PA_level2_All_Collasped<-PA_level2_All_Tran

Final_Mod_PA_level2_All_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level2_All_Collasped,as.numeric))
Final_Mod_PA_level2_All_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level2_All_Collasped,na.rm=TRUE))


#######
#Filter Potable
#######
#All_Subset 
PA_level2_All_Pot<-as.data.frame(PA_level2_All[PA_level2_All$Classification_1 == 'Potable'| PA_level2_All$Classification_1 =='Taxonomy',])
#Transpose 
PA_level2_All_Pot_Tran<-as.data.frame(t(as.matrix(PA_level2_All_Pot)))
PA_level2_All_Pot_Tran<-as.data.frame(PA_level2_All_Pot_Tran[!(row.names(PA_level2_All_Pot_Tran) %in% row.names.remove), ])

Final_Mod_PA_level2_All_Pot_Collasped<-PA_level2_All_Pot_Tran

Final_Mod_PA_level2_All_Pot_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level2_All_Pot_Collasped,as.numeric))
Final_Mod_PA_level2_All_Pot_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level2_All_Pot_Collasped,na.rm=TRUE))

#######
#Filter Reuse
#######
#All_Subset 
PA_level2_All_Reuse<-as.data.frame(PA_level2_All[PA_level2_All$Classification_1 == 'Potable Reuse'| PA_level2_All$Classification_1 =='Taxonomy',])
#Transpose 
PA_level2_All_Reuse_Tran<-as.data.frame(t(as.matrix(PA_level2_All_Reuse)))
PA_level2_All_Reuse_Tran<-as.data.frame(PA_level2_All_Reuse_Tran[!(row.names(PA_level2_All_Reuse_Tran) %in% row.names.remove), ])

Final_Mod_PA_level2_All_Reuse_Collasped<-PA_level2_All_Reuse_Tran

Final_Mod_PA_level2_All_Reuse_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level2_All_Reuse_Collasped,as.numeric))
Final_Mod_PA_level2_All_Reuse_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level2_All_Reuse_Collasped,na.rm=TRUE))


#######
#Filter Reclaimed
#######
#All_Subset 
PA_level2_All_Reclaim<-as.data.frame(PA_level2_All[PA_level2_All$Classification_1 == 'Non-potable Reuse'| PA_level2_All$Classification_1 =='Taxonomy',])
#Transpose 
PA_level2_All_Reclaim_Tran<-as.data.frame(t(as.matrix(PA_level2_All_Reclaim)))
PA_level2_All_Reclaim_Tran<-as.data.frame(PA_level2_All_Reclaim_Tran[!(row.names(PA_level2_All_Reclaim_Tran) %in% row.names.remove), ])

Final_Mod_PA_level2_All_Reclaim_Collasped<-PA_level2_All_Reclaim_Tran

Final_Mod_PA_level2_All_Reclaim_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level2_All_Reclaim_Collasped,as.numeric))
Final_Mod_PA_level2_All_Reclaim_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level2_All_Reclaim_Collasped,na.rm=TRUE))


#######
#Filter Potable and Reuse 
#######
#All_Subset 
PA_level2_All_PotReuse<-as.data.frame(PA_level2_All[PA_level2_All$Classification_1 == 'Potable'| PA_level2_All$Classification_1 =='Potable Reuse'| PA_level2_All$Classification_1 =='Taxonomy',])
#Transpose 
PA_level2_All_PotReuse_Tran<-as.data.frame(t(as.matrix(PA_level2_All_PotReuse)))
PA_level2_All_PotReuse_Tran<-as.data.frame(PA_level2_All_PotReuse_Tran[!(row.names(PA_level2_All_PotReuse_Tran) %in% row.names.remove), ])

Final_Mod_PA_level2_All_PotReuse_Collasped<-PA_level2_All_PotReuse_Tran

Final_Mod_PA_level2_All_PotReuse_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level2_All_PotReuse_Collasped,as.numeric))
Final_Mod_PA_level2_All_PotReuse_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level2_All_PotReuse_Collasped,na.rm=TRUE))

#######
#Filter Bottled Water 
#######
#All_Subset 
PA_level2_All_BottledWater<-as.data.frame(PA_level2_All[PA_level2_All$Sample_Local == 'Bottled Water'| PA_level2_All$Classification_1 =='Taxonomy',])
#Transpose 
PA_level2_All_BottledWater_Tran<-as.data.frame(t(as.matrix(PA_level2_All_BottledWater)))
PA_level2_All_BottledWater_Tran<-as.data.frame(PA_level2_All_BottledWater_Tran[!(row.names(PA_level2_All_BottledWater_Tran) %in% row.names.remove), ])

Final_Mod_PA_level2_All_BottledWater_Collasped<-PA_level2_All_BottledWater_Tran

Final_Mod_PA_level2_All_BottledWater_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level2_All_BottledWater_Collasped,as.numeric))
Final_Mod_PA_level2_All_BottledWater_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level2_All_BottledWater_Collasped,na.rm=TRUE))

#######
#Filter Blank
#######
#All_Subset 
PA_level2_All_Blank<-as.data.frame(PA_level2_All[PA_level2_All$Classification_1 == 'Blank'| PA_level2_All$Classification_1 =='Taxonomy',])
#Transpose 
PA_level2_All_Blank_Tran<-as.data.frame(t(as.matrix(PA_level2_All_Blank)))
PA_level2_All_Blank_Tran<-as.data.frame(PA_level2_All_Blank_Tran[!(row.names(PA_level2_All_Blank_Tran) %in% row.names.remove), ])

Final_Mod_PA_level2_All_Blank_Collasped<-PA_level2_All_Blank_Tran

Final_Mod_PA_level2_All_Blank_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level2_All_Blank_Collasped,as.numeric))
Final_Mod_PA_level2_All_Blank_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level2_All_Blank_Collasped,na.rm=TRUE))

#Export Raw Files

list_PA_level2_All_Mod <- list(
  "All_Collasped" = Final_Mod_PA_level2_All_Collasped,
  "All_Pot_Collasped" = Final_Mod_PA_level2_All_Pot_Collasped ,
  "All_Reuse_Collasped" = Final_Mod_PA_level2_All_Reuse_Collasped,
  "All_Reclaim_Collasped" = Final_Mod_PA_level2_All_Reclaim_Collasped,
  "All_PotReuse_Collasped" = Final_Mod_PA_level2_All_PotReuse_Collasped,
  "Blanks"=Final_Mod_PA_level2_All_Blank_Collasped,
  "BottledWater"=Final_Mod_PA_level2_All_BottledWater_Collasped)

write.xlsx(list_PA_level2_All_Mod,colNames=TRUE,rowNames=TRUE, file = "C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Core/All/PA_level2_ALL_Collasped_Mod.xlsx")

#Merge and Export Averaged for Heatmap 
list_PA_level2_All_mean_Mod <- list(
  "All_Collasped_Mean" = Final_Mod_PA_level2_All_Collasped_Mean,
  "All_Pot_Collasped_Mean" = Final_Mod_PA_level2_All_Pot_Collasped_Mean ,
  "All_Reuse_Collasped_Mean" = Final_Mod_PA_level2_All_Reuse_Collasped_Mean,
  "All_Reclaim_Collasped_Mean" = Final_Mod_PA_level2_All_Reclaim_Collasped_Mean,
  "All_PotReuse_Collasped_Mean" = Final_Mod_PA_level2_All_PotReuse_Collasped_Mean,
  "Blanks_Mean"=Final_Mod_PA_level2_All_Blank_Collasped_Mean,
  "BottledWater_Mean"=Final_Mod_PA_level2_All_BottledWater_Collasped_Mean)

write.xlsx(list_PA_level2_All_mean_Mod,colNames=TRUE,rowNames=TRUE, file = "C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Core/All/PA_level2_ALL_Collasped_mean_Mod.xlsx")


list_rn_PA_level2_All_mean_Mod <- lapply(list_PA_level2_All_mean_Mod, rownames_to_column)
PA_level2_All_mean_Mod<- list_rn_PA_level2_All_mean_Mod %>% reduce(full_join, by= "rowname")

write.xlsx(PA_level2_All_mean_Mod,colNames=TRUE, file = "C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Core/All/PA_level2_All_mean_Mod.xlsx")

#######
#######
#Filter All, WTP and DS, level6
#######
#######
PA_level6_All<-as.data.frame(Final_Mod_level6_Collapsed_OTUs_Trans_Merged_PA[Final_Mod_level6_Collapsed_OTUs_Trans_Merged_PA$Sample_Local == 'POC'| Final_Mod_level6_Collapsed_OTUs_Trans_Merged_PA$Sample_Local =='Bottled Water'| Final_Mod_level6_Collapsed_OTUs_Trans_Merged_PA$Sample_Local =='POU'| Final_Mod_level6_Collapsed_OTUs_Trans_Merged_PA$Sample_Local =='Blank',])

#Transpose, remove unwanted rows, leaving just SV  
PA_level6_All_Tran<-as.data.frame(t(as.matrix(PA_level6_All)))
PA_level6_All_Tran<-PA_level6_All_Tran[!(row.names(PA_level6_All_Tran) %in% row.names.remove), ]

Final_Mod_PA_level6_All_Collasped<-PA_level6_All_Tran

Final_Mod_PA_level6_All_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level6_All_Collasped,as.numeric))
Final_Mod_PA_level6_All_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level6_All_Collasped,na.rm=TRUE))


#######
#Filter Potable
#######
#All_Subset 
PA_level6_All_Pot<-as.data.frame(PA_level6_All[PA_level6_All$Classification_1 == 'Potable'| PA_level6_All$Classification_1 =='Taxonomy',])
#Transpose 
PA_level6_All_Pot_Tran<-as.data.frame(t(as.matrix(PA_level6_All_Pot)))
PA_level6_All_Pot_Tran<-as.data.frame(PA_level6_All_Pot_Tran[!(row.names(PA_level6_All_Pot_Tran) %in% row.names.remove), ])

Final_Mod_PA_level6_All_Pot_Collasped<-PA_level6_All_Pot_Tran

Final_Mod_PA_level6_All_Pot_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level6_All_Pot_Collasped,as.numeric))
Final_Mod_PA_level6_All_Pot_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level6_All_Pot_Collasped,na.rm=TRUE))

#######
#Filter Reuse
#######
#All_Subset 
PA_level6_All_Reuse<-as.data.frame(PA_level6_All[PA_level6_All$Classification_1 == 'Potable Reuse'| PA_level6_All$Classification_1 =='Taxonomy',])
#Transpose 
PA_level6_All_Reuse_Tran<-as.data.frame(t(as.matrix(PA_level6_All_Reuse)))
PA_level6_All_Reuse_Tran<-as.data.frame(PA_level6_All_Reuse_Tran[!(row.names(PA_level6_All_Reuse_Tran) %in% row.names.remove), ])

Final_Mod_PA_level6_All_Reuse_Collasped<-PA_level6_All_Reuse_Tran

Final_Mod_PA_level6_All_Reuse_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level6_All_Reuse_Collasped,as.numeric))
Final_Mod_PA_level6_All_Reuse_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level6_All_Reuse_Collasped,na.rm=TRUE))


#######
#Filter Reclaimed
#######
#All_Subset 
PA_level6_All_Reclaim<-as.data.frame(PA_level6_All[PA_level6_All$Classification_1 == 'Non-potable Reuse'| PA_level6_All$Classification_1 =='Taxonomy',])
#Transpose 
PA_level6_All_Reclaim_Tran<-as.data.frame(t(as.matrix(PA_level6_All_Reclaim)))
PA_level6_All_Reclaim_Tran<-as.data.frame(PA_level6_All_Reclaim_Tran[!(row.names(PA_level6_All_Reclaim_Tran) %in% row.names.remove), ])

Final_Mod_PA_level6_All_Reclaim_Collasped<-PA_level6_All_Reclaim_Tran

Final_Mod_PA_level6_All_Reclaim_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level6_All_Reclaim_Collasped,as.numeric))
Final_Mod_PA_level6_All_Reclaim_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level6_All_Reclaim_Collasped,na.rm=TRUE))


#######
#Filter Potable and Reuse 
#######
#All_Subset 
PA_level6_All_PotReuse<-as.data.frame(PA_level6_All[PA_level6_All$Classification_1 == 'Potable'| PA_level6_All$Classification_1 =='Potable Reuse'| PA_level6_All$Classification_1 =='Taxonomy',])
#Transpose 
PA_level6_All_PotReuse_Tran<-as.data.frame(t(as.matrix(PA_level6_All_PotReuse)))
PA_level6_All_PotReuse_Tran<-as.data.frame(PA_level6_All_PotReuse_Tran[!(row.names(PA_level6_All_PotReuse_Tran) %in% row.names.remove), ])

Final_Mod_PA_level6_All_PotReuse_Collasped<-PA_level6_All_PotReuse_Tran

Final_Mod_PA_level6_All_PotReuse_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level6_All_PotReuse_Collasped,as.numeric))
Final_Mod_PA_level6_All_PotReuse_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level6_All_PotReuse_Collasped,na.rm=TRUE))

#######
#Filter Bottled Water 
#######
#All_Subset 
PA_level6_All_BottledWater<-as.data.frame(PA_level6_All[PA_level6_All$Sample_Local == 'Bottled Water'| PA_level6_All$Classification_1 =='Taxonomy',])
#Transpose 
PA_level6_All_BottledWater_Tran<-as.data.frame(t(as.matrix(PA_level6_All_BottledWater)))
PA_level6_All_BottledWater_Tran<-as.data.frame(PA_level6_All_BottledWater_Tran[!(row.names(PA_level6_All_BottledWater_Tran) %in% row.names.remove), ])

Final_Mod_PA_level6_All_BottledWater_Collasped<-PA_level6_All_BottledWater_Tran

Final_Mod_PA_level6_All_BottledWater_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level6_All_BottledWater_Collasped,as.numeric))
Final_Mod_PA_level6_All_BottledWater_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level6_All_BottledWater_Collasped,na.rm=TRUE))

#######
#Filter Blank
#######
#All_Subset 
PA_level6_All_Blank<-as.data.frame(PA_level6_All[PA_level6_All$Classification_1 == 'Blank'| PA_level6_All$Classification_1 =='Taxonomy',])
#Transpose 
PA_level6_All_Blank_Tran<-as.data.frame(t(as.matrix(PA_level6_All_Blank)))
PA_level6_All_Blank_Tran<-as.data.frame(PA_level6_All_Blank_Tran[!(row.names(PA_level6_All_Blank_Tran) %in% row.names.remove), ])

Final_Mod_PA_level6_All_Blank_Collasped<-PA_level6_All_Blank_Tran

Final_Mod_PA_level6_All_Blank_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level6_All_Blank_Collasped,as.numeric))
Final_Mod_PA_level6_All_Blank_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level6_All_Blank_Collasped,na.rm=TRUE))

#Export Raw Files

list_PA_level6_All_Mod <- list(
  "All_Collasped" = Final_Mod_PA_level6_All_Collasped,
  "All_Pot_Collasped" = Final_Mod_PA_level6_All_Pot_Collasped ,
  "All_Reuse_Collasped" = Final_Mod_PA_level6_All_Reuse_Collasped,
  "All_Reclaim_Collasped" = Final_Mod_PA_level6_All_Reclaim_Collasped,
  "All_PotReuse_Collasped" = Final_Mod_PA_level6_All_PotReuse_Collasped,
  "Blanks"=Final_Mod_PA_level6_All_Blank_Collasped,
  "BottledWater"=Final_Mod_PA_level6_All_BottledWater_Collasped)

write.xlsx(list_PA_level6_All_Mod,colNames=TRUE,rowNames=TRUE, file = "C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Core/All/PA_level6_ALL_Collasped_Mod.xlsx")

#Merge and Export Averaged for Heatmap 
list_PA_level6_All_mean_Mod <- list(
  "All_Collasped_Mean" = Final_Mod_PA_level6_All_Collasped_Mean,
  "All_Pot_Collasped_Mean" = Final_Mod_PA_level6_All_Pot_Collasped_Mean ,
  "All_Reuse_Collasped_Mean" = Final_Mod_PA_level6_All_Reuse_Collasped_Mean,
  "All_Reclaim_Collasped_Mean" = Final_Mod_PA_level6_All_Reclaim_Collasped_Mean,
  "All_PotReuse_Collasped_Mean" = Final_Mod_PA_level6_All_PotReuse_Collasped_Mean,
  "Blanks_Mean"=Final_Mod_PA_level6_All_Blank_Collasped_Mean,
  "BottledWater_Mean"=Final_Mod_PA_level6_All_BottledWater_Collasped_Mean)

write.xlsx(list_PA_level6_All_mean_Mod,colNames=TRUE,rowNames=TRUE, file = "C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Core/All/PA_level6_ALL_Collasped_mean_Mod.xlsx")


list_rn_PA_level6_All_mean_Mod <- lapply(list_PA_level6_All_mean_Mod, rownames_to_column)
PA_level6_All_mean_Mod<- list_rn_PA_level6_All_mean_Mod %>% reduce(full_join, by= "rowname")

write.xlsx(PA_level6_All_mean_Mod,colNames=TRUE, file = "C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Core/All/PA_level6_All_mean_Mod.xlsx")

#######
#######
#Filter All, WTP and DS, level2
#######
#######
PA_level2_All<-as.data.frame(Final_Mod_level2_Collapsed_OTUs_Trans_Merged_PA[Final_Mod_level2_Collapsed_OTUs_Trans_Merged_PA$Sample_Local == 'POC'| Final_Mod_level2_Collapsed_OTUs_Trans_Merged_PA$Sample_Local =='Bottled Water'| Final_Mod_level2_Collapsed_OTUs_Trans_Merged_PA$Sample_Local =='POU'| Final_Mod_level2_Collapsed_OTUs_Trans_Merged_PA$Sample_Local =='Blank',])

#Transpose, remove unwanted rows, leaving just SV  
PA_level2_All_Tran<-as.data.frame(t(as.matrix(PA_level2_All)))
PA_level2_All_Tran<-PA_level2_All_Tran[!(row.names(PA_level2_All_Tran) %in% row.names.remove), ]

Final_Mod_PA_level2_All_Collasped<-PA_level2_All_Tran

Final_Mod_PA_level2_All_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level2_All_Collasped,as.numeric))
Final_Mod_PA_level2_All_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level2_All_Collasped,na.rm=TRUE))


#######
#Filter Potable
#######
#All_Subset 
PA_level2_All_Pot<-as.data.frame(PA_level2_All[PA_level2_All$Classification_1 == 'Potable'| PA_level2_All$Classification_1 =='Taxonomy',])
#Transpose 
PA_level2_All_Pot_Tran<-as.data.frame(t(as.matrix(PA_level2_All_Pot)))
PA_level2_All_Pot_Tran<-as.data.frame(PA_level2_All_Pot_Tran[!(row.names(PA_level2_All_Pot_Tran) %in% row.names.remove), ])

Final_Mod_PA_level2_All_Pot_Collasped<-PA_level2_All_Pot_Tran

Final_Mod_PA_level2_All_Pot_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level2_All_Pot_Collasped,as.numeric))
Final_Mod_PA_level2_All_Pot_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level2_All_Pot_Collasped,na.rm=TRUE))

#######
#Filter Reuse
#######
#All_Subset 
PA_level2_All_Reuse<-as.data.frame(PA_level2_All[PA_level2_All$Classification_1 == 'Potable Reuse'| PA_level2_All$Classification_1 =='Taxonomy',])
#Transpose 
PA_level2_All_Reuse_Tran<-as.data.frame(t(as.matrix(PA_level2_All_Reuse)))
PA_level2_All_Reuse_Tran<-as.data.frame(PA_level2_All_Reuse_Tran[!(row.names(PA_level2_All_Reuse_Tran) %in% row.names.remove), ])

Final_Mod_PA_level2_All_Reuse_Collasped<-PA_level2_All_Reuse_Tran

Final_Mod_PA_level2_All_Reuse_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level2_All_Reuse_Collasped,as.numeric))
Final_Mod_PA_level2_All_Reuse_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level2_All_Reuse_Collasped,na.rm=TRUE))


#######
#Filter Reclaimed
#######
#All_Subset 
PA_level2_All_Reclaim<-as.data.frame(PA_level2_All[PA_level2_All$Classification_1 == 'Non-potable Reuse'| PA_level2_All$Classification_1 =='Taxonomy',])
#Transpose 
PA_level2_All_Reclaim_Tran<-as.data.frame(t(as.matrix(PA_level2_All_Reclaim)))
PA_level2_All_Reclaim_Tran<-as.data.frame(PA_level2_All_Reclaim_Tran[!(row.names(PA_level2_All_Reclaim_Tran) %in% row.names.remove), ])

Final_Mod_PA_level2_All_Reclaim_Collasped<-PA_level2_All_Reclaim_Tran

Final_Mod_PA_level2_All_Reclaim_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level2_All_Reclaim_Collasped,as.numeric))
Final_Mod_PA_level2_All_Reclaim_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level2_All_Reclaim_Collasped,na.rm=TRUE))


#######
#Filter Potable and Reuse 
#######
#All_Subset 
PA_level2_All_PotReuse<-as.data.frame(PA_level2_All[PA_level2_All$Classification_1 == 'Potable'| PA_level2_All$Classification_1 =='Potable Reuse'| PA_level2_All$Classification_1 =='Taxonomy',])
#Transpose 
PA_level2_All_PotReuse_Tran<-as.data.frame(t(as.matrix(PA_level2_All_PotReuse)))
PA_level2_All_PotReuse_Tran<-as.data.frame(PA_level2_All_PotReuse_Tran[!(row.names(PA_level2_All_PotReuse_Tran) %in% row.names.remove), ])

Final_Mod_PA_level2_All_PotReuse_Collasped<-PA_level2_All_PotReuse_Tran

Final_Mod_PA_level2_All_PotReuse_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level2_All_PotReuse_Collasped,as.numeric))
Final_Mod_PA_level2_All_PotReuse_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level2_All_PotReuse_Collasped,na.rm=TRUE))

#######
#Filter Bottled Water 
#######
#All_Subset 
PA_level2_All_BottledWater<-as.data.frame(PA_level2_All[PA_level2_All$Sample_Local == 'Bottled Water'| PA_level2_All$Classification_1 =='Taxonomy',])
#Transpose 
PA_level2_All_BottledWater_Tran<-as.data.frame(t(as.matrix(PA_level2_All_BottledWater)))
PA_level2_All_BottledWater_Tran<-as.data.frame(PA_level2_All_BottledWater_Tran[!(row.names(PA_level2_All_BottledWater_Tran) %in% row.names.remove), ])

Final_Mod_PA_level2_All_BottledWater_Collasped<-PA_level2_All_BottledWater_Tran

Final_Mod_PA_level2_All_BottledWater_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level2_All_BottledWater_Collasped,as.numeric))
Final_Mod_PA_level2_All_BottledWater_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level2_All_BottledWater_Collasped,na.rm=TRUE))

#######
#Filter Blank
#######
#All_Subset 
PA_level2_All_Blank<-as.data.frame(PA_level2_All[PA_level2_All$Classification_1 == 'Blank'| PA_level2_All$Classification_1 =='Taxonomy',])
#Transpose 
PA_level2_All_Blank_Tran<-as.data.frame(t(as.matrix(PA_level2_All_Blank)))
PA_level2_All_Blank_Tran<-as.data.frame(PA_level2_All_Blank_Tran[!(row.names(PA_level2_All_Blank_Tran) %in% row.names.remove), ])

Final_Mod_PA_level2_All_Blank_Collasped<-PA_level2_All_Blank_Tran

Final_Mod_PA_level2_All_Blank_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level2_All_Blank_Collasped,as.numeric))
Final_Mod_PA_level2_All_Blank_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level2_All_Blank_Collasped,na.rm=TRUE))

#Export Raw Files

list_PA_level2_All_Mod <- list(
  "All_Collasped" = Final_Mod_PA_level2_All_Collasped,
  "All_Pot_Collasped" = Final_Mod_PA_level2_All_Pot_Collasped ,
  "All_Reuse_Collasped" = Final_Mod_PA_level2_All_Reuse_Collasped,
  "All_Reclaim_Collasped" = Final_Mod_PA_level2_All_Reclaim_Collasped,
  "All_PotReuse_Collasped" = Final_Mod_PA_level2_All_PotReuse_Collasped,
  "Blanks"=Final_Mod_PA_level2_All_Blank_Collasped,
  "BottledWater"=Final_Mod_PA_level2_All_BottledWater_Collasped)

write.xlsx(list_PA_level2_All_Mod,colNames=TRUE,rowNames=TRUE, file = "C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Core/All/PA_level2_ALL_Collasped_Mod.xlsx")

#Merge and Export Averaged for Heatmap 
list_PA_level2_All_mean_Mod <- list(
  "All_Collasped_Mean" = Final_Mod_PA_level2_All_Collasped_Mean,
  "All_Pot_Collasped_Mean" = Final_Mod_PA_level2_All_Pot_Collasped_Mean ,
  "All_Reuse_Collasped_Mean" = Final_Mod_PA_level2_All_Reuse_Collasped_Mean,
  "All_Reclaim_Collasped_Mean" = Final_Mod_PA_level2_All_Reclaim_Collasped_Mean,
  "All_PotReuse_Collasped_Mean" = Final_Mod_PA_level2_All_PotReuse_Collasped_Mean,
  "Blanks_Mean"=Final_Mod_PA_level2_All_Blank_Collasped_Mean,
  "BottledWater_Mean"=Final_Mod_PA_level2_All_BottledWater_Collasped_Mean)

write.xlsx(list_PA_level2_All_mean_Mod,colNames=TRUE,rowNames=TRUE, file = "C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Core/All/PA_level2_ALL_Collasped_mean_Mod.xlsx")


list_rn_PA_level2_All_mean_Mod <- lapply(list_PA_level2_All_mean_Mod, rownames_to_column)
PA_level2_All_mean_Mod<- list_rn_PA_level2_All_mean_Mod %>% reduce(full_join, by= "rowname")

write.xlsx(PA_level2_All_mean_Mod,colNames=TRUE, file = "C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Core/All/PA_level2_All_mean_Mod.xlsx")

#######
#######
#Filter POC, WTP and DS, level6
#######
#######
PA_level6_POC<-as.data.frame(Final_Mod_level6_Collapsed_OTUs_Trans_Merged_PA[Final_Mod_level6_Collapsed_OTUs_Trans_Merged_PA$Sample_Local == 'POC'| Final_Mod_level6_Collapsed_OTUs_Trans_Merged_PA$Sample_Local =='Blank',])

#Transpose, remove unwanted rows, leaving just SV  
PA_level6_POC_Tran<-as.data.frame(t(as.matrix(PA_level6_POC)))
PA_level6_POC_Tran<-PA_level6_POC_Tran[!(row.names(PA_level6_POC_Tran) %in% row.names.remove), ]

Final_Mod_PA_level6_POC_Collasped<-PA_level6_POC_Tran

Final_Mod_PA_level6_POC_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level6_POC_Collasped,as.numeric))
Final_Mod_PA_level6_POC_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level6_POC_Collasped,na.rm=TRUE))


#######
#Filter Potable
#######
#POC_Subset 
PA_level6_POC_Pot<-as.data.frame(PA_level6_POC[PA_level6_POC$Classification_1 == 'Potable'| PA_level6_POC$Classification_1 =='Taxonomy',])
#Transpose 
PA_level6_POC_Pot_Tran<-as.data.frame(t(as.matrix(PA_level6_POC_Pot)))
PA_level6_POC_Pot_Tran<-as.data.frame(PA_level6_POC_Pot_Tran[!(row.names(PA_level6_POC_Pot_Tran) %in% row.names.remove), ])

Final_Mod_PA_level6_POC_Pot_Collasped<-PA_level6_POC_Pot_Tran

Final_Mod_PA_level6_POC_Pot_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level6_POC_Pot_Collasped,as.numeric))
Final_Mod_PA_level6_POC_Pot_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level6_POC_Pot_Collasped,na.rm=TRUE))

#######
#Filter Reuse
#######
#POC_Subset 
PA_level6_POC_Reuse<-as.data.frame(PA_level6_POC[PA_level6_POC$Classification_1 == 'Potable Reuse'| PA_level6_POC$Classification_1 =='Taxonomy',])
#Transpose 
PA_level6_POC_Reuse_Tran<-as.data.frame(t(as.matrix(PA_level6_POC_Reuse)))
PA_level6_POC_Reuse_Tran<-as.data.frame(PA_level6_POC_Reuse_Tran[!(row.names(PA_level6_POC_Reuse_Tran) %in% row.names.remove), ])

Final_Mod_PA_level6_POC_Reuse_Collasped<-PA_level6_POC_Reuse_Tran

Final_Mod_PA_level6_POC_Reuse_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level6_POC_Reuse_Collasped,as.numeric))
Final_Mod_PA_level6_POC_Reuse_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level6_POC_Reuse_Collasped,na.rm=TRUE))


#######
#Filter Reclaimed
#######
#POC_Subset 
PA_level6_POC_Reclaim<-as.data.frame(PA_level6_POC[PA_level6_POC$Classification_1 == 'Non-potable Reuse'| PA_level6_POC$Classification_1 =='Taxonomy',])
#Transpose 
PA_level6_POC_Reclaim_Tran<-as.data.frame(t(as.matrix(PA_level6_POC_Reclaim)))
PA_level6_POC_Reclaim_Tran<-as.data.frame(PA_level6_POC_Reclaim_Tran[!(row.names(PA_level6_POC_Reclaim_Tran) %in% row.names.remove), ])

Final_Mod_PA_level6_POC_Reclaim_Collasped<-PA_level6_POC_Reclaim_Tran

Final_Mod_PA_level6_POC_Reclaim_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level6_POC_Reclaim_Collasped,as.numeric))
Final_Mod_PA_level6_POC_Reclaim_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level6_POC_Reclaim_Collasped,na.rm=TRUE))


#######
#Filter Potable and Reuse 
#######
#POC_Subset 
PA_level6_POC_PotReuse<-as.data.frame(PA_level6_POC[PA_level6_POC$Classification_1 == 'Potable'| PA_level6_POC$Classification_1 =='Potable Reuse'| PA_level6_POC$Classification_1 =='Taxonomy',])
#Transpose 
PA_level6_POC_PotReuse_Tran<-as.data.frame(t(as.matrix(PA_level6_POC_PotReuse)))
PA_level6_POC_PotReuse_Tran<-as.data.frame(PA_level6_POC_PotReuse_Tran[!(row.names(PA_level6_POC_PotReuse_Tran) %in% row.names.remove), ])

Final_Mod_PA_level6_POC_PotReuse_Collasped<-PA_level6_POC_PotReuse_Tran

Final_Mod_PA_level6_POC_PotReuse_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level6_POC_PotReuse_Collasped,as.numeric))
Final_Mod_PA_level6_POC_PotReuse_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level6_POC_PotReuse_Collasped,na.rm=TRUE))

#######
#Filter Bottled Water 
#######
#POC_Subset 
PA_level6_POC_BottledWater<-as.data.frame(PA_level6_POC[PA_level6_POC$Sample_Local == 'Bottled Water'| PA_level6_POC$Classification_1 =='Taxonomy',])
#Transpose 
PA_level6_POC_BottledWater_Tran<-as.data.frame(t(as.matrix(PA_level6_POC_BottledWater)))
PA_level6_POC_BottledWater_Tran<-as.data.frame(PA_level6_POC_BottledWater_Tran[!(row.names(PA_level6_POC_BottledWater_Tran) %in% row.names.remove), ])

Final_Mod_PA_level6_POC_BottledWater_Collasped<-PA_level6_POC_BottledWater_Tran

Final_Mod_PA_level6_POC_BottledWater_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level6_POC_BottledWater_Collasped,as.numeric))
Final_Mod_PA_level6_POC_BottledWater_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level6_POC_BottledWater_Collasped,na.rm=TRUE))

#######
#Filter Blank
#######
#POC_Subset 
PA_level6_POC_Blank<-as.data.frame(PA_level6_POC[PA_level6_POC$Classification_1 == 'Blank'| PA_level6_POC$Classification_1 =='Taxonomy',])
#Transpose 
PA_level6_POC_Blank_Tran<-as.data.frame(t(as.matrix(PA_level6_POC_Blank)))
PA_level6_POC_Blank_Tran<-as.data.frame(PA_level6_POC_Blank_Tran[!(row.names(PA_level6_POC_Blank_Tran) %in% row.names.remove), ])

Final_Mod_PA_level6_POC_Blank_Collasped<-PA_level6_POC_Blank_Tran

Final_Mod_PA_level6_POC_Blank_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level6_POC_Blank_Collasped,as.numeric))
Final_Mod_PA_level6_POC_Blank_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level6_POC_Blank_Collasped,na.rm=TRUE))

#Export Raw Files

list_PA_level6_POC_Mod <- list(
  "POC_Collasped" = Final_Mod_PA_level6_POC_Collasped,
  "POC_Pot_Collasped" = Final_Mod_PA_level6_POC_Pot_Collasped ,
  "POC_Reuse_Collasped" = Final_Mod_PA_level6_POC_Reuse_Collasped,
  "POC_Reclaim_Collasped" = Final_Mod_PA_level6_POC_Reclaim_Collasped,
  "POC_PotReuse_Collasped" = Final_Mod_PA_level6_POC_PotReuse_Collasped,
  "Blanks"=Final_Mod_PA_level6_POC_Blank_Collasped,
  "BottledWater"=Final_Mod_PA_level6_POC_BottledWater_Collasped)

write.xlsx(list_PA_level6_POC_Mod,colNames=TRUE,rowNames=TRUE, file = "C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Core/POC/PA_level6_POC_Collasped_Mod.xlsx")

#Merge and Export Averaged for Heatmap 
list_PA_level6_POC_mean_Mod <- list(
  "POC_Collasped_Mean" = Final_Mod_PA_level6_POC_Collasped_Mean,
  "POC_Pot_Collasped_Mean" = Final_Mod_PA_level6_POC_Pot_Collasped_Mean ,
  "POC_Reuse_Collasped_Mean" = Final_Mod_PA_level6_POC_Reuse_Collasped_Mean,
  "POC_Reclaim_Collasped_Mean" = Final_Mod_PA_level6_POC_Reclaim_Collasped_Mean,
  "POC_PotReuse_Collasped_Mean" = Final_Mod_PA_level6_POC_PotReuse_Collasped_Mean,
  "Blanks_Mean"=Final_Mod_PA_level6_POC_Blank_Collasped_Mean,
  "BottledWater_Mean"=Final_Mod_PA_level6_POC_BottledWater_Collasped_Mean)

write.xlsx(list_PA_level6_POC_mean_Mod,colNames=TRUE,rowNames=TRUE, file = "C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Core/POC/PA_level6_POC_Collasped_mean_Mod.xlsx")


list_rn_PA_level6_POC_mean_Mod <- lapply(list_PA_level6_POC_mean_Mod, rownames_to_column)
PA_level6_POC_mean_Mod<- list_rn_PA_level6_POC_mean_Mod %>% reduce(full_join, by= "rowname")

write.xlsx(PA_level6_POC_mean_Mod,colNames=TRUE, file = "C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Core/POC/PA_level6_POC_mean_Mod.xlsx")

#######
#######
#Filter POC, WTP and DS, level2
#######
#######
PA_level2_POC<-as.data.frame(Final_Mod_level2_Collapsed_OTUs_Trans_Merged_PA[Final_Mod_level2_Collapsed_OTUs_Trans_Merged_PA$Sample_Local =='Bottled Water'| Final_Mod_level2_Collapsed_OTUs_Trans_Merged_PA$Sample_Local =='POC'| Final_Mod_level2_Collapsed_OTUs_Trans_Merged_PA$Sample_Local =='Blank',])

#Transpose, remove unwanted rows, leaving just SV  
PA_level2_POC_Tran<-as.data.frame(t(as.matrix(PA_level2_POC)))
PA_level2_POC_Tran<-PA_level2_POC_Tran[!(row.names(PA_level2_POC_Tran) %in% row.names.remove), ]

Final_Mod_PA_level2_POC_Collasped<-PA_level2_POC_Tran

Final_Mod_PA_level2_POC_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level2_POC_Collasped,as.numeric))
Final_Mod_PA_level2_POC_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level2_POC_Collasped,na.rm=TRUE))


#######
#Filter Potable
#######
#POC_Subset 
PA_level2_POC_Pot<-as.data.frame(PA_level2_POC[PA_level2_POC$Classification_1 == 'Potable'| PA_level2_POC$Classification_1 =='Taxonomy',])
#Transpose 
PA_level2_POC_Pot_Tran<-as.data.frame(t(as.matrix(PA_level2_POC_Pot)))
PA_level2_POC_Pot_Tran<-as.data.frame(PA_level2_POC_Pot_Tran[!(row.names(PA_level2_POC_Pot_Tran) %in% row.names.remove), ])

Final_Mod_PA_level2_POC_Pot_Collasped<-PA_level2_POC_Pot_Tran

Final_Mod_PA_level2_POC_Pot_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level2_POC_Pot_Collasped,as.numeric))
Final_Mod_PA_level2_POC_Pot_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level2_POC_Pot_Collasped,na.rm=TRUE))

#######
#Filter Reuse
#######
#POC_Subset 
PA_level2_POC_Reuse<-as.data.frame(PA_level2_POC[PA_level2_POC$Classification_1 == 'Potable Reuse'| PA_level2_POC$Classification_1 =='Taxonomy',])
#Transpose 
PA_level2_POC_Reuse_Tran<-as.data.frame(t(as.matrix(PA_level2_POC_Reuse)))
PA_level2_POC_Reuse_Tran<-as.data.frame(PA_level2_POC_Reuse_Tran[!(row.names(PA_level2_POC_Reuse_Tran) %in% row.names.remove), ])

Final_Mod_PA_level2_POC_Reuse_Collasped<-PA_level2_POC_Reuse_Tran

Final_Mod_PA_level2_POC_Reuse_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level2_POC_Reuse_Collasped,as.numeric))
Final_Mod_PA_level2_POC_Reuse_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level2_POC_Reuse_Collasped,na.rm=TRUE))


#######
#Filter Reclaimed
#######
#POC_Subset 
PA_level2_POC_Reclaim<-as.data.frame(PA_level2_POC[PA_level2_POC$Classification_1 == 'Non-potable Reuse'| PA_level2_POC$Classification_1 =='Taxonomy',])
#Transpose 
PA_level2_POC_Reclaim_Tran<-as.data.frame(t(as.matrix(PA_level2_POC_Reclaim)))
PA_level2_POC_Reclaim_Tran<-as.data.frame(PA_level2_POC_Reclaim_Tran[!(row.names(PA_level2_POC_Reclaim_Tran) %in% row.names.remove), ])

Final_Mod_PA_level2_POC_Reclaim_Collasped<-PA_level2_POC_Reclaim_Tran

Final_Mod_PA_level2_POC_Reclaim_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level2_POC_Reclaim_Collasped,as.numeric))
Final_Mod_PA_level2_POC_Reclaim_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level2_POC_Reclaim_Collasped,na.rm=TRUE))


#######
#Filter Potable and Reuse 
#######
#POC_Subset 
PA_level2_POC_PotReuse<-as.data.frame(PA_level2_POC[PA_level2_POC$Classification_1 == 'Potable'| PA_level2_POC$Classification_1 =='Potable Reuse'| PA_level2_POC$Classification_1 =='Taxonomy',])
#Transpose 
PA_level2_POC_PotReuse_Tran<-as.data.frame(t(as.matrix(PA_level2_POC_PotReuse)))
PA_level2_POC_PotReuse_Tran<-as.data.frame(PA_level2_POC_PotReuse_Tran[!(row.names(PA_level2_POC_PotReuse_Tran) %in% row.names.remove), ])

Final_Mod_PA_level2_POC_PotReuse_Collasped<-PA_level2_POC_PotReuse_Tran

Final_Mod_PA_level2_POC_PotReuse_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level2_POC_PotReuse_Collasped,as.numeric))
Final_Mod_PA_level2_POC_PotReuse_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level2_POC_PotReuse_Collasped,na.rm=TRUE))

#######
#Filter Bottled Water 
#######
#POC_Subset 
PA_level2_POC_BottledWater<-as.data.frame(PA_level2_POC[PA_level2_POC$Sample_Local == 'Bottled Water'| PA_level2_POC$Classification_1 =='Taxonomy',])
#Transpose 
PA_level2_POC_BottledWater_Tran<-as.data.frame(t(as.matrix(PA_level2_POC_BottledWater)))
PA_level2_POC_BottledWater_Tran<-as.data.frame(PA_level2_POC_BottledWater_Tran[!(row.names(PA_level2_POC_BottledWater_Tran) %in% row.names.remove), ])

Final_Mod_PA_level2_POC_BottledWater_Collasped<-PA_level2_POC_BottledWater_Tran

Final_Mod_PA_level2_POC_BottledWater_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level2_POC_BottledWater_Collasped,as.numeric))
Final_Mod_PA_level2_POC_BottledWater_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level2_POC_BottledWater_Collasped,na.rm=TRUE))

#######
#Filter Blank
#######
#POC_Subset 
PA_level2_POC_Blank<-as.data.frame(PA_level2_POC[PA_level2_POC$Classification_1 == 'Blank'| PA_level2_POC$Classification_1 =='Taxonomy',])
#Transpose 
PA_level2_POC_Blank_Tran<-as.data.frame(t(as.matrix(PA_level2_POC_Blank)))
PA_level2_POC_Blank_Tran<-as.data.frame(PA_level2_POC_Blank_Tran[!(row.names(PA_level2_POC_Blank_Tran) %in% row.names.remove), ])

Final_Mod_PA_level2_POC_Blank_Collasped<-PA_level2_POC_Blank_Tran

Final_Mod_PA_level2_POC_Blank_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level2_POC_Blank_Collasped,as.numeric))
Final_Mod_PA_level2_POC_Blank_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level2_POC_Blank_Collasped,na.rm=TRUE))

#Export Raw Files

list_PA_level2_POC_Mod <- list(
  "POC_Collasped" = Final_Mod_PA_level2_POC_Collasped,
  "POC_Pot_Collasped" = Final_Mod_PA_level2_POC_Pot_Collasped ,
  "POC_Reuse_Collasped" = Final_Mod_PA_level2_POC_Reuse_Collasped,
  "POC_Reclaim_Collasped" = Final_Mod_PA_level2_POC_Reclaim_Collasped,
  "POC_PotReuse_Collasped" = Final_Mod_PA_level2_POC_PotReuse_Collasped,
  "Blanks"=Final_Mod_PA_level2_POC_Blank_Collasped,
  "BottledWater"=Final_Mod_PA_level2_POC_BottledWater_Collasped)

write.xlsx(list_PA_level2_POC_Mod,colNames=TRUE,rowNames=TRUE, file = "C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Core/POC/PA_level2_POC_Collasped_Mod.xlsx")

#Merge and Export Averaged for Heatmap 
list_PA_level2_POC_mean_Mod <- list(
  "POC_Collasped_Mean" = Final_Mod_PA_level2_POC_Collasped_Mean,
  "POC_Pot_Collasped_Mean" = Final_Mod_PA_level2_POC_Pot_Collasped_Mean ,
  "POC_Reuse_Collasped_Mean" = Final_Mod_PA_level2_POC_Reuse_Collasped_Mean,
  "POC_Reclaim_Collasped_Mean" = Final_Mod_PA_level2_POC_Reclaim_Collasped_Mean,
  "POC_PotReuse_Collasped_Mean" = Final_Mod_PA_level2_POC_PotReuse_Collasped_Mean,
  "Blanks_Mean"=Final_Mod_PA_level2_POC_Blank_Collasped_Mean,
  "BottledWater_Mean"=Final_Mod_PA_level2_POC_BottledWater_Collasped_Mean)

write.xlsx(list_PA_level2_POC_mean_Mod,colNames=TRUE,rowNames=TRUE, file = "C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Core/POC/PA_level2_POC_Collasped_mean_Mod.xlsx")


list_rn_PA_level2_POC_mean_Mod <- lapply(list_PA_level2_POC_mean_Mod, rownames_to_column)
PA_level2_POC_mean_Mod<- list_rn_PA_level2_POC_mean_Mod %>% reduce(full_join, by= "rowname")

write.xlsx(PA_level2_POC_mean_Mod,colNames=TRUE, file = "C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Core/POC/PA_level2_POC_mean_Mod.xlsx")

#######
#######
#Filter POU, WTP and DS, level6
#######
#######
PA_level6_POU<-as.data.frame(Final_Mod_level6_Collapsed_OTUs_Trans_Merged_PA[Final_Mod_level6_Collapsed_OTUs_Trans_Merged_PA$Sample_Local == 'POC'| Final_Mod_level6_Collapsed_OTUs_Trans_Merged_PA$Sample_Local =='Bottled Water'| Final_Mod_level6_Collapsed_OTUs_Trans_Merged_PA$Sample_Local =='POU'| Final_Mod_level6_Collapsed_OTUs_Trans_Merged_PA$Sample_Local =='Blank',])

#Transpose, remove unwanted rows, leaving just SV  
PA_level6_POU_Tran<-as.data.frame(t(as.matrix(PA_level6_POU)))
PA_level6_POU_Tran<-PA_level6_POU_Tran[!(row.names(PA_level6_POU_Tran) %in% row.names.remove), ]

Final_Mod_PA_level6_POU_Collasped<-PA_level6_POU_Tran

Final_Mod_PA_level6_POU_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level6_POU_Collasped,as.numeric))
Final_Mod_PA_level6_POU_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level6_POU_Collasped,na.rm=TRUE))


#######
#Filter Potable
#######
#POU_Subset 
PA_level6_POU_Pot<-as.data.frame(PA_level6_POU[PA_level6_POU$Classification_1 == 'Potable'| PA_level6_POU$Classification_1 =='Taxonomy',])
#Transpose 
PA_level6_POU_Pot_Tran<-as.data.frame(t(as.matrix(PA_level6_POU_Pot)))
PA_level6_POU_Pot_Tran<-as.data.frame(PA_level6_POU_Pot_Tran[!(row.names(PA_level6_POU_Pot_Tran) %in% row.names.remove), ])

Final_Mod_PA_level6_POU_Pot_Collasped<-PA_level6_POU_Pot_Tran

Final_Mod_PA_level6_POU_Pot_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level6_POU_Pot_Collasped,as.numeric))
Final_Mod_PA_level6_POU_Pot_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level6_POU_Pot_Collasped,na.rm=TRUE))

#######
#Filter Reuse
#######
#POU_Subset 
PA_level6_POU_Reuse<-as.data.frame(PA_level6_POU[PA_level6_POU$Classification_1 == 'Potable Reuse'| PA_level6_POU$Classification_1 =='Taxonomy',])
#Transpose 
PA_level6_POU_Reuse_Tran<-as.data.frame(t(as.matrix(PA_level6_POU_Reuse)))
PA_level6_POU_Reuse_Tran<-as.data.frame(PA_level6_POU_Reuse_Tran[!(row.names(PA_level6_POU_Reuse_Tran) %in% row.names.remove), ])

Final_Mod_PA_level6_POU_Reuse_Collasped<-PA_level6_POU_Reuse_Tran

Final_Mod_PA_level6_POU_Reuse_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level6_POU_Reuse_Collasped,as.numeric))
Final_Mod_PA_level6_POU_Reuse_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level6_POU_Reuse_Collasped,na.rm=TRUE))


#######
#Filter Reclaimed
#######
#POU_Subset 
PA_level6_POU_Reclaim<-as.data.frame(PA_level6_POU[PA_level6_POU$Classification_1 == 'Non-potable Reuse'| PA_level6_POU$Classification_1 =='Taxonomy',])
#Transpose 
PA_level6_POU_Reclaim_Tran<-as.data.frame(t(as.matrix(PA_level6_POU_Reclaim)))
PA_level6_POU_Reclaim_Tran<-as.data.frame(PA_level6_POU_Reclaim_Tran[!(row.names(PA_level6_POU_Reclaim_Tran) %in% row.names.remove), ])

Final_Mod_PA_level6_POU_Reclaim_Collasped<-PA_level6_POU_Reclaim_Tran

Final_Mod_PA_level6_POU_Reclaim_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level6_POU_Reclaim_Collasped,as.numeric))
Final_Mod_PA_level6_POU_Reclaim_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level6_POU_Reclaim_Collasped,na.rm=TRUE))


#######
#Filter Potable and Reuse 
#######
#POU_Subset 
PA_level6_POU_PotReuse<-as.data.frame(PA_level6_POU[PA_level6_POU$Classification_1 == 'Potable'| PA_level6_POU$Classification_1 =='Potable Reuse'| PA_level6_POU$Classification_1 =='Taxonomy',])
#Transpose 
PA_level6_POU_PotReuse_Tran<-as.data.frame(t(as.matrix(PA_level6_POU_PotReuse)))
PA_level6_POU_PotReuse_Tran<-as.data.frame(PA_level6_POU_PotReuse_Tran[!(row.names(PA_level6_POU_PotReuse_Tran) %in% row.names.remove), ])

Final_Mod_PA_level6_POU_PotReuse_Collasped<-PA_level6_POU_PotReuse_Tran

Final_Mod_PA_level6_POU_PotReuse_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level6_POU_PotReuse_Collasped,as.numeric))
Final_Mod_PA_level6_POU_PotReuse_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level6_POU_PotReuse_Collasped,na.rm=TRUE))

#######
#Filter Bottled Water 
#######
#POU_Subset 
PA_level6_POU_BottledWater<-as.data.frame(PA_level6_POU[PA_level6_POU$Sample_Local == 'Bottled Water'| PA_level6_POU$Classification_1 =='Taxonomy',])
#Transpose 
PA_level6_POU_BottledWater_Tran<-as.data.frame(t(as.matrix(PA_level6_POU_BottledWater)))
PA_level6_POU_BottledWater_Tran<-as.data.frame(PA_level6_POU_BottledWater_Tran[!(row.names(PA_level6_POU_BottledWater_Tran) %in% row.names.remove), ])

Final_Mod_PA_level6_POU_BottledWater_Collasped<-PA_level6_POU_BottledWater_Tran

Final_Mod_PA_level6_POU_BottledWater_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level6_POU_BottledWater_Collasped,as.numeric))
Final_Mod_PA_level6_POU_BottledWater_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level6_POU_BottledWater_Collasped,na.rm=TRUE))

#######
#Filter Blank
#######
#POU_Subset 
PA_level6_POU_Blank<-as.data.frame(PA_level6_POU[PA_level6_POU$Classification_1 == 'Blank'| PA_level6_POU$Classification_1 =='Taxonomy',])
#Transpose 
PA_level6_POU_Blank_Tran<-as.data.frame(t(as.matrix(PA_level6_POU_Blank)))
PA_level6_POU_Blank_Tran<-as.data.frame(PA_level6_POU_Blank_Tran[!(row.names(PA_level6_POU_Blank_Tran) %in% row.names.remove), ])

Final_Mod_PA_level6_POU_Blank_Collasped<-PA_level6_POU_Blank_Tran

Final_Mod_PA_level6_POU_Blank_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level6_POU_Blank_Collasped,as.numeric))
Final_Mod_PA_level6_POU_Blank_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level6_POU_Blank_Collasped,na.rm=TRUE))

#Export Raw Files

list_PA_level6_POU_Mod <- list(
  "POU_Collasped" = Final_Mod_PA_level6_POU_Collasped,
  "POU_Pot_Collasped" = Final_Mod_PA_level6_POU_Pot_Collasped ,
  "POU_Reuse_Collasped" = Final_Mod_PA_level6_POU_Reuse_Collasped,
  "POU_Reclaim_Collasped" = Final_Mod_PA_level6_POU_Reclaim_Collasped,
  "POU_PotReuse_Collasped" = Final_Mod_PA_level6_POU_PotReuse_Collasped,
  "Blanks"=Final_Mod_PA_level6_POU_Blank_Collasped,
  "BottledWater"=Final_Mod_PA_level6_POU_BottledWater_Collasped)

write.xlsx(list_PA_level6_POU_Mod,colNames=TRUE,rowNames=TRUE, file = "C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Core/POU/PA_level6_POU_Collasped_Mod.xlsx")

#Merge and Export Averaged for Heatmap 
list_PA_level6_POU_mean_Mod <- list(
  "POU_Collasped_Mean" = Final_Mod_PA_level6_POU_Collasped_Mean,
  "POU_Pot_Collasped_Mean" = Final_Mod_PA_level6_POU_Pot_Collasped_Mean ,
  "POU_Reuse_Collasped_Mean" = Final_Mod_PA_level6_POU_Reuse_Collasped_Mean,
  "POU_Reclaim_Collasped_Mean" = Final_Mod_PA_level6_POU_Reclaim_Collasped_Mean,
  "POU_PotReuse_Collasped_Mean" = Final_Mod_PA_level6_POU_PotReuse_Collasped_Mean,
  "Blanks_Mean"=Final_Mod_PA_level6_POU_Blank_Collasped_Mean,
  "BottledWater_Mean"=Final_Mod_PA_level6_POU_BottledWater_Collasped_Mean)

write.xlsx(list_PA_level6_POU_mean_Mod,colNames=TRUE,rowNames=TRUE, file = "C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Core/POU/PA_level6_POU_Collasped_mean_Mod.xlsx")


list_rn_PA_level6_POU_mean_Mod <- lapply(list_PA_level6_POU_mean_Mod, rownames_to_column)
PA_level6_POU_mean_Mod<- list_rn_PA_level6_POU_mean_Mod %>% reduce(full_join, by= "rowname")

write.xlsx(PA_level6_POU_mean_Mod,colNames=TRUE, file = "C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Core/POU/PA_level6_POU_mean_Mod.xlsx")


#######
#######
#Filter POU, WTP and DS, level2
#######
#######
PA_level2_POU<-as.data.frame(Final_Mod_level2_Collapsed_OTUs_Trans_Merged_PA[Final_Mod_level2_Collapsed_OTUs_Trans_Merged_PA$Sample_Local == 'POC'| Final_Mod_level2_Collapsed_OTUs_Trans_Merged_PA$Sample_Local =='Bottled Water'| Final_Mod_level2_Collapsed_OTUs_Trans_Merged_PA$Sample_Local =='POU'| Final_Mod_level2_Collapsed_OTUs_Trans_Merged_PA$Sample_Local =='Blank',])

#Transpose, remove unwanted rows, leaving just SV  
PA_level2_POU_Tran<-as.data.frame(t(as.matrix(PA_level2_POU)))
PA_level2_POU_Tran<-PA_level2_POU_Tran[!(row.names(PA_level2_POU_Tran) %in% row.names.remove), ]

Final_Mod_PA_level2_POU_Collasped<-PA_level2_POU_Tran

Final_Mod_PA_level2_POU_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level2_POU_Collasped,as.numeric))
Final_Mod_PA_level2_POU_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level2_POU_Collasped,na.rm=TRUE))


#######
#Filter Potable
#######
#POU_Subset 
PA_level2_POU_Pot<-as.data.frame(PA_level2_POU[PA_level2_POU$Classification_1 == 'Potable'| PA_level2_POU$Classification_1 =='Taxonomy',])
#Transpose 
PA_level2_POU_Pot_Tran<-as.data.frame(t(as.matrix(PA_level2_POU_Pot)))
PA_level2_POU_Pot_Tran<-as.data.frame(PA_level2_POU_Pot_Tran[!(row.names(PA_level2_POU_Pot_Tran) %in% row.names.remove), ])

Final_Mod_PA_level2_POU_Pot_Collasped<-PA_level2_POU_Pot_Tran

Final_Mod_PA_level2_POU_Pot_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level2_POU_Pot_Collasped,as.numeric))
Final_Mod_PA_level2_POU_Pot_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level2_POU_Pot_Collasped,na.rm=TRUE))

#######
#Filter Reuse
#######
#POU_Subset 
PA_level2_POU_Reuse<-as.data.frame(PA_level2_POU[PA_level2_POU$Classification_1 == 'Potable Reuse'| PA_level2_POU$Classification_1 =='Taxonomy',])
#Transpose 
PA_level2_POU_Reuse_Tran<-as.data.frame(t(as.matrix(PA_level2_POU_Reuse)))
PA_level2_POU_Reuse_Tran<-as.data.frame(PA_level2_POU_Reuse_Tran[!(row.names(PA_level2_POU_Reuse_Tran) %in% row.names.remove), ])

Final_Mod_PA_level2_POU_Reuse_Collasped<-PA_level2_POU_Reuse_Tran

Final_Mod_PA_level2_POU_Reuse_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level2_POU_Reuse_Collasped,as.numeric))
Final_Mod_PA_level2_POU_Reuse_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level2_POU_Reuse_Collasped,na.rm=TRUE))


#######
#Filter Reclaimed
#######
#POU_Subset 
PA_level2_POU_Reclaim<-as.data.frame(PA_level2_POU[PA_level2_POU$Classification_1 == 'Non-potable Reuse'| PA_level2_POU$Classification_1 =='Taxonomy',])
#Transpose 
PA_level2_POU_Reclaim_Tran<-as.data.frame(t(as.matrix(PA_level2_POU_Reclaim)))
PA_level2_POU_Reclaim_Tran<-as.data.frame(PA_level2_POU_Reclaim_Tran[!(row.names(PA_level2_POU_Reclaim_Tran) %in% row.names.remove), ])

Final_Mod_PA_level2_POU_Reclaim_Collasped<-PA_level2_POU_Reclaim_Tran

Final_Mod_PA_level2_POU_Reclaim_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level2_POU_Reclaim_Collasped,as.numeric))
Final_Mod_PA_level2_POU_Reclaim_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level2_POU_Reclaim_Collasped,na.rm=TRUE))


#######
#Filter Potable and Reuse 
#######
#POU_Subset 
PA_level2_POU_PotReuse<-as.data.frame(PA_level2_POU[PA_level2_POU$Classification_1 == 'Potable'| PA_level2_POU$Classification_1 =='Potable Reuse'| PA_level2_POU$Classification_1 =='Taxonomy',])
#Transpose 
PA_level2_POU_PotReuse_Tran<-as.data.frame(t(as.matrix(PA_level2_POU_PotReuse)))
PA_level2_POU_PotReuse_Tran<-as.data.frame(PA_level2_POU_PotReuse_Tran[!(row.names(PA_level2_POU_PotReuse_Tran) %in% row.names.remove), ])

Final_Mod_PA_level2_POU_PotReuse_Collasped<-PA_level2_POU_PotReuse_Tran

Final_Mod_PA_level2_POU_PotReuse_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level2_POU_PotReuse_Collasped,as.numeric))
Final_Mod_PA_level2_POU_PotReuse_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level2_POU_PotReuse_Collasped,na.rm=TRUE))

#######
#Filter Bottled Water 
#######
#POU_Subset 
PA_level2_POU_BottledWater<-as.data.frame(PA_level2_POU[PA_level2_POU$Sample_Local == 'Bottled Water'| PA_level2_POU$Classification_1 =='Taxonomy',])
#Transpose 
PA_level2_POU_BottledWater_Tran<-as.data.frame(t(as.matrix(PA_level2_POU_BottledWater)))
PA_level2_POU_BottledWater_Tran<-as.data.frame(PA_level2_POU_BottledWater_Tran[!(row.names(PA_level2_POU_BottledWater_Tran) %in% row.names.remove), ])

Final_Mod_PA_level2_POU_BottledWater_Collasped<-PA_level2_POU_BottledWater_Tran

Final_Mod_PA_level2_POU_BottledWater_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level2_POU_BottledWater_Collasped,as.numeric))
Final_Mod_PA_level2_POU_BottledWater_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level2_POU_BottledWater_Collasped,na.rm=TRUE))

#######
#Filter Blank
#######
#POU_Subset 
PA_level2_POU_Blank<-as.data.frame(PA_level2_POU[PA_level2_POU$Classification_1 == 'Blank'| PA_level2_POU$Classification_1 =='Taxonomy',])
#Transpose 
PA_level2_POU_Blank_Tran<-as.data.frame(t(as.matrix(PA_level2_POU_Blank)))
PA_level2_POU_Blank_Tran<-as.data.frame(PA_level2_POU_Blank_Tran[!(row.names(PA_level2_POU_Blank_Tran) %in% row.names.remove), ])

Final_Mod_PA_level2_POU_Blank_Collasped<-PA_level2_POU_Blank_Tran

Final_Mod_PA_level2_POU_Blank_Collasped[]<-as.data.frame(lapply(Final_Mod_PA_level2_POU_Blank_Collasped,as.numeric))
Final_Mod_PA_level2_POU_Blank_Collasped_Mean<-as.data.frame(rowMeans(Final_Mod_PA_level2_POU_Blank_Collasped,na.rm=TRUE))

#Export Raw Files

list_PA_level2_POU_Mod <- list(
  "POU_Collasped" = Final_Mod_PA_level2_POU_Collasped,
  "POU_Pot_Collasped" = Final_Mod_PA_level2_POU_Pot_Collasped ,
  "POU_Reuse_Collasped" = Final_Mod_PA_level2_POU_Reuse_Collasped,
  "POU_Reclaim_Collasped" = Final_Mod_PA_level2_POU_Reclaim_Collasped,
  "POU_PotReuse_Collasped" = Final_Mod_PA_level2_POU_PotReuse_Collasped,
  "Blanks"=Final_Mod_PA_level2_POU_Blank_Collasped,
  "BottledWater"=Final_Mod_PA_level2_POU_BottledWater_Collasped)

write.xlsx(list_PA_level2_POU_Mod,colNames=TRUE,rowNames=TRUE, file = "C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Core/POU/PA_level2_POU_Collasped_Mod.xlsx")

#Merge and Export Averaged for Heatmap 
list_PA_level2_POU_mean_Mod <- list(
  "POU_Collasped_Mean" = Final_Mod_PA_level2_POU_Collasped_Mean,
  "POU_Pot_Collasped_Mean" = Final_Mod_PA_level2_POU_Pot_Collasped_Mean ,
  "POU_Reuse_Collasped_Mean" = Final_Mod_PA_level2_POU_Reuse_Collasped_Mean,
  "POU_Reclaim_Collasped_Mean" = Final_Mod_PA_level2_POU_Reclaim_Collasped_Mean,
  "POU_PotReuse_Collasped_Mean" = Final_Mod_PA_level2_POU_PotReuse_Collasped_Mean,
  "Blanks_Mean"=Final_Mod_PA_level2_POU_Blank_Collasped_Mean,
  "BottledWater_Mean"=Final_Mod_PA_level2_POU_BottledWater_Collasped_Mean)

write.xlsx(list_PA_level2_POU_mean_Mod,colNames=TRUE,rowNames=TRUE, file = "C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Core/POU/PA_level2_POU_Collasped_mean_Mod.xlsx")


list_rn_PA_level2_POU_mean_Mod <- lapply(list_PA_level2_POU_mean_Mod, rownames_to_column)
PA_level2_POU_mean_Mod<- list_rn_PA_level2_POU_mean_Mod %>% reduce(full_join, by= "rowname")

write.xlsx(PA_level2_POU_mean_Mod,colNames=TRUE, file = "C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Core/POU/PA_level2_POU_mean_Mod.xlsx")


