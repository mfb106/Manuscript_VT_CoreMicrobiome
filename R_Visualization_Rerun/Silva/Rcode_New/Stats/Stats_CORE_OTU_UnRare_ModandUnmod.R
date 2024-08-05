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
library(deseq2)
setwd("C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures")
getwd()
BiocManager::install("DESeq2")

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


#extract data version, one for RA and one for PA 
table_all_data<-as.data.frame(table_all$data)
taxonomy_all_data<-taxonomy_all$data

table_all_data_PA<-as.data.frame(table_all$data)
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

###Merged tax data and sample metadata for final table 
Final_Merged<-transform(merge(Meta.All.POCPOU,Combined_Trans,by=0), row.names=Row.names, Row.names=NULL)
Final_Merged_PA<-transform(merge(Meta.All.POCPOU,Combined_Trans_PA,by=0), row.names=Row.names, Row.names=NULL)

Final_OTU_Mod_Collapsed_OTUs_Trans_Merged_PA<-Final_Merged_PA

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

#Import Export For Core Mean
# PA_level2_All_mean<-read.xlsx("C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Core/All/PA_level2_All_mean.xlsx", sheet= 'Sheet 1', rowNames = TRUE)
# PA_level6_All_mean<-read.xlsx("C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Core/All/PA_level6_All_mean.xlsx", sheet= 'Sheet 1', rowNames = TRUE)
# 
# PA_level2_POU_mean<-read.xlsx("C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Core/POU/PA_level2_POU_mean.xlsx", sheet= 'Sheet 1', rowNames = TRUE)
# PA_level6_POU_mean<-read.xlsx("C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Core/POU/PA_level6_POU_mean.xlsx", sheet= 'Sheet 1', rowNames = TRUE)
# 
# PA_level2_POC_mean<-read.xlsx("C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Core/POC/PA_level2_POC_mean.xlsx", sheet= 'Sheet 1', rowNames = TRUE)
# PA_level6_POC_mean<-read.xlsx("C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Core/POC/PA_level6_POC_mean.xlsx", sheet= 'Sheet 1', rowNames = TRUE)

PA_means_OTU<-read.xlsx("C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Core/Combined_edited_OTU.xlsx", sheet= 'UnMod', rowNames = FALSE)
PA_means_OTU_Mod<-read.xlsx("C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Core/Combined_edited_OTU.xlsx", sheet= 'Mod', rowNames = FALSE)



#Set Rownames Name
rownames(PA_means_OTU)<- PA_means_OTU$OTU
rownames(PA_means_OTU_Mod) <- PA_means_OTU_Mod$OTU

#Drop Additional Taxa Columns
Final_PA_means_OTU <- subset(PA_means_OTU, select=-c(OTU,POC_Blank,POU_Blank,POC_BottledWater,POU_BottledWater))
Final_PA_means_OTU_Mod <- subset(PA_means_OTU_Mod, select=-c(OTU,POC_Blank,POU_Blank,POC_BottledWater,POU_BottledWater))

#Identify Max PA, can adjust .75|.80|.85 if needed  
Final_PA_means_OTU[,"MaxPA"]<-apply(Final_PA_means_OTU,1,max, na.rm = TRUE)
Final_PA_means_OTU<-filter(Final_PA_means_OTU, MaxPA > 0.75)
Final_PA_means_OTU <- tibble::rownames_to_column(Final_PA_means_OTU, "OTU")
Final_PA_means_OTU <- subset(Final_PA_means_OTU, select=-c(MaxPA,POC_POU_All ))

Final_PA_means_OTU_Long<-gather(Final_PA_means_OTU,Classifcation,Presence, 2:17, factor_key = TRUE)



Final_PA_means_OTU_Mod[,"MaxPA"]<-apply(Final_PA_means_OTU_Mod,1,max, na.rm = TRUE)
Final_PA_means_OTU_Mod<-filter(Final_PA_means_OTU_Mod, MaxPA > 0.75)
Final_PA_means_OTU_Mod <- tibble::rownames_to_column(Final_PA_means_OTU_Mod, "OTU_Mod")
Final_PA_means_OTU_Mod <- subset(Final_PA_means_OTU_Mod, select=-c(MaxPA,POC_POU_All))

Final_PA_means_OTU_Mod_Long<-gather(Final_PA_means_OTU_Mod,Classifcation,Presence, 2:17, factor_key = TRUE)


####Consider, Renaming X axis variables / column names. 
#Final_PA_means_OTU_Long$Classifcation <- gsub('POC_POU_All', 'All Water Uses, POC and POU', Final_PA_means_OTU_Long$Classifcation)
Final_PA_means_OTU_Long$Classifcation <- gsub('\\<POC_POU_Potable\\>', 'Potable Conventional, POC and POU', Final_PA_means_OTU_Long$Classifcation)
Final_PA_means_OTU_Long$Classifcation <- gsub('\\<POC_POU_NonPotableReuse\\>', 'Non-Potable Reuse, POC and POU', Final_PA_means_OTU_Long$Classifcation)
Final_PA_means_OTU_Long$Classifcation <- gsub('\\<POC_POU_PotableReuse\\>', 'Potable Reuse, POC and POU', Final_PA_means_OTU_Long$Classifcation)

Final_PA_means_OTU_Long$Classifcation <- gsub('\\<POC_POU_Blank\\>', 'Blanks', Final_PA_means_OTU_Long$Classifcation)
Final_PA_means_OTU_Long$Classifcation <- gsub('\\<POC_POU_BottledWater\\>', 'Bottled Water', Final_PA_means_OTU_Long$Classifcation)

Final_PA_means_OTU_Long$Classifcation <- gsub('\\<POC_All\\>', 'All Water Uses, POC', Final_PA_means_OTU_Long$Classifcation)
Final_PA_means_OTU_Long$Classifcation <- gsub('\\<POC_Potable\\>', 'Potable Conventional, POC', Final_PA_means_OTU_Long$Classifcation)
Final_PA_means_OTU_Long$Classifcation <- gsub('\\<POC_NonPotableReuse\\>', 'Non-Potable Reuse, POC', Final_PA_means_OTU_Long$Classifcation)
Final_PA_means_OTU_Long$Classifcation <- gsub('\\<POC_PotableReuse\\>', 'Potable Reuse, POC', Final_PA_means_OTU_Long$Classifcation)

Final_PA_means_OTU_Long$Classifcation <- gsub('\\<POU_All\\>', 'All Water Uses, POU', Final_PA_means_OTU_Long$Classifcation)
Final_PA_means_OTU_Long$Classifcation <- gsub('\\<POU_Potable\\>', 'Potable Conventional, POU', Final_PA_means_OTU_Long$Classifcation)
Final_PA_means_OTU_Long$Classifcation <- gsub('\\<POU_NonPotableReuse\\>', 'Non-Potable Reuse, POU', Final_PA_means_OTU_Long$Classifcation)
Final_PA_means_OTU_Long$Classifcation <- gsub('\\<POU_PotableReuse\\>', 'Potable Reuse, POU', Final_PA_means_OTU_Long$Classifcation)


Final_PA_means_OTU_Long$Classifcation <- factor(Final_PA_means_OTU_Long$Classifcation, levels=c(
  "Potable Conventional, POC and POU", 
  "Potable Conventional, POC", 
  "Potable Conventional, POU",
  "Potable Reuse, POC and POU", 
  "Potable Reuse, POC", 
  "Potable Reuse, POU",
  "Non-Potable Reuse, POC and POU", 
  "Non-Potable Reuse, POC",
  "Non-Potable Reuse, POU",
  "Bottled Water",
  "Blanks"))

#Drop categories that were not included above. 
Final_PA_means_OTU_Long<-Final_PA_means_OTU_Long %>% drop_na(Classifcation)

#Adjust Width of Axis
Final_PA_means_OTU_Long$Classifcationnew = str_wrap(Final_PA_means_OTU_Long$Classifcation, width = 10)

# Final_Join_OTU_Mod_Long$Classifcation <- gsub('Reuse', 'Potable Reuse', Final_Join_OTU_Mod_Long$Classifcation)
# Final_Join_OTU_Mod_Long$Classifcation <- gsub('Reclaimed', 'Non-Potable Reuse', Final_Join_OTU_Mod_Long$Classifcation)
# Final_Join_OTU_Mod_Long$Classifcation <- factor(Final_Join_OTU_Mod_Long$Classifcation, levels=c("All Potable Samples", "WTP Potable Samples", "DS Potable Samples",
#                                                                                             "All Potable Reuse Samples", "WTP Potable Reuse Samples", "DS Potable Reuse Samples",
#                                                                                             "All Non-Potable Reuse Samples", "WTP Non-Potable Reuse Samples", "DS Non-Potable Reuse Samples"))

#Heatmap, OTU 
Final_Join_OTU_Long_Norare_Mod<-ggplot(Final_PA_means_OTU_Long,aes(Classifcation,reorder(OTU,Presence),fill=Presence))+
  geom_tile()+
  scale_fill_gradient(low="white", high="dark blue") +
  #geom_text(aes(label = round(Presence, 2))) +
  theme_ipsum()+ 
    scale_colour_Publication()+ theme_Publication()+
  theme(legend.text=element_text(size=10), legend.title=element_text(size=10))+
  theme(legend.position="bottom")+
  theme(axis.text.x = element_text(angle = 0), axis.ticks = element_blank())+ #, vjust = 0, hjust = 0
  scale_x_discrete(labels = wrap_format(10))+
  ylab("Taxanomic Classification, OTU")
Final_Join_OTU_Long_Norare_Mod
leg_OTU <- get_legend(Final_Join_OTU_Long_Norare_Mod)
as_ggplot(leg_OTU)


####Consider, Renaming X axis variables / column names. 
#Final_PA_means_OTU_Mod_Long$Classifcation <- gsub('POC_POU_All', 'All Water Uses, POC and POU', Final_PA_means_OTU_Mod_Long$Classifcation)
Final_PA_means_OTU_Mod_Long$Classifcation <- gsub('\\<POC_POU_Potable\\>', 'Potable Conventional, POC and POU', Final_PA_means_OTU_Mod_Long$Classifcation)
Final_PA_means_OTU_Mod_Long$Classifcation <- gsub('\\<POC_POU_NonPotableReuse\\>', 'Non-Potable Reuse, POC and POU', Final_PA_means_OTU_Mod_Long$Classifcation)
Final_PA_means_OTU_Mod_Long$Classifcation <- gsub('\\<POC_POU_PotableReuse\\>', 'Potable Reuse, POC and POU', Final_PA_means_OTU_Mod_Long$Classifcation)

Final_PA_means_OTU_Mod_Long$Classifcation <- gsub('\\<POC_POU_Blank\\>', 'Blanks', Final_PA_means_OTU_Mod_Long$Classifcation)
Final_PA_means_OTU_Mod_Long$Classifcation <- gsub('\\<POC_POU_BottledWater\\>', 'Bottled Water', Final_PA_means_OTU_Mod_Long$Classifcation)

Final_PA_means_OTU_Mod_Long$Classifcation <- gsub('\\<POC_All\\>', 'All Water Uses, POC', Final_PA_means_OTU_Mod_Long$Classifcation)
Final_PA_means_OTU_Mod_Long$Classifcation <- gsub('\\<POC_Potable\\>', 'Potable Conventional, POC', Final_PA_means_OTU_Mod_Long$Classifcation)
Final_PA_means_OTU_Mod_Long$Classifcation <- gsub('\\<POC_NonPotableReuse\\>', 'Non-Potable Reuse, POC', Final_PA_means_OTU_Mod_Long$Classifcation)
Final_PA_means_OTU_Mod_Long$Classifcation <- gsub('\\<POC_PotableReuse\\>', 'Potable Reuse, POC', Final_PA_means_OTU_Mod_Long$Classifcation)

Final_PA_means_OTU_Mod_Long$Classifcation <- gsub('\\<POU_All\\>', 'All Water Uses, POU', Final_PA_means_OTU_Mod_Long$Classifcation)
Final_PA_means_OTU_Mod_Long$Classifcation <- gsub('\\<POU_Potable\\>', 'Potable Conventional, POU', Final_PA_means_OTU_Mod_Long$Classifcation)
Final_PA_means_OTU_Mod_Long$Classifcation <- gsub('\\<POU_NonPotableReuse\\>', 'Non-Potable Reuse, POU', Final_PA_means_OTU_Mod_Long$Classifcation)
Final_PA_means_OTU_Mod_Long$Classifcation <- gsub('\\<POU_PotableReuse\\>', 'Potable Reuse, POU', Final_PA_means_OTU_Mod_Long$Classifcation)


Final_PA_means_OTU_Mod_Long$Classifcation <- factor(Final_PA_means_OTU_Mod_Long$Classifcation, levels=c(
  "Potable Conventional, POC and POU", 
  "Potable Conventional, POC", 
  "Potable Conventional, POU",
  "Potable Reuse, POC and POU", 
  "Potable Reuse, POC", 
  "Potable Reuse, POU",
  "Non-Potable Reuse, POC and POU", 
  "Non-Potable Reuse, POC",
  "Non-Potable Reuse, POU",
  "Bottled Water",
  "Blanks"))

#Drop categories that were not included above. 
Final_PA_means_OTU_Mod_Long<-Final_PA_means_OTU_Mod_Long %>% drop_na(Classifcation)

#Adjust Width of Axis
Final_PA_means_OTU_Mod_Long$Classifcationnew = str_wrap(Final_PA_means_OTU_Mod_Long$Classifcation, width = 10)

# Final_Join_OTU_Mod_Long$Classifcation <- gsub('Reuse', 'Potable Reuse', Final_Join_OTU_Mod_Long$Classifcation)
# Final_Join_OTU_Mod_Long$Classifcation <- gsub('Reclaimed', 'Non-Potable Reuse', Final_Join_OTU_Mod_Long$Classifcation)
# Final_Join_OTU_Mod_Long$Classifcation <- factor(Final_Join_OTU_Mod_Long$Classifcation, levels=c("All Potable Samples", "WTP Potable Samples", "DS Potable Samples",
#                                                                                             "All Potable Reuse Samples", "WTP Potable Reuse Samples", "DS Potable Reuse Samples",
#                                                                                             "All Non-Potable Reuse Samples", "WTP Non-Potable Reuse Samples", "DS Non-Potable Reuse Samples"))

#Heatmap, OTU_Mod 
Final_Join_OTU_Mod_Long_Norare_Mod<-ggplot(Final_PA_means_OTU_Mod_Long,aes(Classifcation,reorder(OTU_Mod,Presence),fill=Presence))+
  geom_tile()+
  scale_fill_gradient(low="white", high="dark orange") +
  #geom_text(aes(label = round(Presence, 2))) +
  theme_ipsum()+ 
  scale_colour_Publication()+ theme_Publication()+
  theme(legend.text=element_text(size=10), legend.title=element_text(size=10))+
  theme(legend.position="bottom")+
  theme(axis.text.x = element_text(angle = 0), axis.ticks = element_blank())+ #, vjust = 0, hjust = 0
  scale_x_discrete(labels = wrap_format(10))+
  ylab("Taxanomic Classification, OTU_Mod")
Final_Join_OTU_Mod_Long_Norare_Mod
leg_OTU_Mod <- get_legend(Final_Join_OTU_Mod_Long_Norare_Mod)
as_ggplot(leg_OTU_Mod)



##Export
ggsave("Final_Join_OTU_Long_Norare_Mod.png",plot = Final_Join_OTU_Long_Norare_Mod,path = "C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Core/Figures/", width = 20, height = 16, units = "in", dpi = 300)
ggsave("Final_Join_OTU_Mod_Long_Norare_Mod.png",plot = Final_Join_OTU_Mod_Long_Norare_Mod,path = "C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Core/Figures/", width = 20, height = 16, units = "in", dpi = 300)

#Combo Figure
Manuscript_Combo <- ggarrange(Final_Join_OTU_Long_Norare_Mod,
                              Final_Join_OTU_Mod_Long_Norare_Mod,
                              labels = c("A", "B"),
                              font.label = list(size = 12, color = "black"),
                              ncol = 2, nrow = 1,
                              legend = "none")
Manuscript_Combo

ggsave("Manuscript_Combo_OTU.png",plot = Manuscript_Combo,path = "C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Core/Figures/", width = 8.5, height = 8, units = "in", dpi = 300)
ggsave("leg_OTU_Mod.png",plot = leg_OTU_Mod,path = "C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Core/Figures/", width = 6, height = 2, units = "in", dpi = 300)
ggsave("leg_OTU.png",plot = leg_OTU,path = "C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Core/Figures/", width = 6, height = 2, units = "in", dpi = 300)

#Export Excels 
write.xlsx(Final_PA_means_OTU,"C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Core/Figures/Final_PA_means_OTU_Mod_Norare.xlsx")
write.xlsx(Final_PA_means_OTU_Mod,"C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Core/Figures/Final_PA_means_OTU_Mod_Mod_Norare.xlsx")




