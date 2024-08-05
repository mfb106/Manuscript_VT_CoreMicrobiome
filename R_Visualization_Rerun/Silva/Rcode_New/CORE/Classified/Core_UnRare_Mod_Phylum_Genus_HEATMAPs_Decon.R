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



PA_means_Phylum<-read.xlsx("C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Core/Combined_edited_mod2.xlsx", sheet= 'Phylum', rowNames = FALSE)
PA_means_Genus<-read.xlsx("C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Core/Combined_edited_mod2.xlsx", sheet= 'Genus', rowNames = FALSE)


#Set Rownames Name
rownames(PA_means_Phylum)<- PA_means_Phylum$Phylum
rownames(PA_means_Genus) <- PA_means_Genus$Genus

#Drop Additional Taxa Columns
Final_PA_means_Phylum <- subset(PA_means_Phylum, select=-c(Phylum,POC_Blank,POU_Blank,POC_BottledWater,POU_BottledWater))
Final_PA_means_Genus <- subset(PA_means_Genus, select=-c(Genus,POC_Blank,POU_Blank,POC_BottledWater,POU_BottledWater))

#Identify Max PA, can adjust .75|.80|.85 if needed  
Final_PA_means_Phylum[,"MaxPA"]<-apply(Final_PA_means_Phylum,1,max, na.rm = TRUE)
Final_PA_means_Phylum<-filter(Final_PA_means_Phylum, MaxPA > 0.75)
Final_PA_means_Phylum <- tibble::rownames_to_column(Final_PA_means_Phylum, "Phylum")
Final_PA_means_Phylum <- subset(Final_PA_means_Phylum, select=-c(MaxPA,POC_POU_All ))

Final_PA_means_Phylum_Long<-gather(Final_PA_means_Phylum,Classifcation,Presence, 2:17, factor_key = TRUE)



Final_PA_means_Genus[,"MaxPA"]<-apply(Final_PA_means_Genus,1,max, na.rm = TRUE)
Final_PA_means_Genus<-filter(Final_PA_means_Genus, MaxPA > 0.80)
Final_PA_means_Genus <- tibble::rownames_to_column(Final_PA_means_Genus, "Genus")
Final_PA_means_Genus <- subset(Final_PA_means_Genus, select=-c(MaxPA,POC_POU_All))

Final_PA_means_Genus_Long<-gather(Final_PA_means_Genus,Classifcation,Presence, 2:17, factor_key = TRUE)


####Consider, Renaming X axis variables / column names. 
#Final_PA_means_Phylum_Long$Classifcation <- gsub('POC_POU_All', 'All Water Uses, POC and POU', Final_PA_means_Phylum_Long$Classifcation)
Final_PA_means_Phylum_Long$Classifcation <- gsub('\\<POC_POU_Potable\\>', 'Potable Conventional, POC and POU', Final_PA_means_Phylum_Long$Classifcation)
Final_PA_means_Phylum_Long$Classifcation <- gsub('\\<POC_POU_NonPotableReuse\\>', 'Non-Potable Reuse, POC and POU', Final_PA_means_Phylum_Long$Classifcation)
Final_PA_means_Phylum_Long$Classifcation <- gsub('\\<POC_POU_PotableReuse\\>', 'Potable Reuse, POC and POU', Final_PA_means_Phylum_Long$Classifcation)

Final_PA_means_Phylum_Long$Classifcation <- gsub('\\<POC_POU_Blank\\>', 'Blanks', Final_PA_means_Phylum_Long$Classifcation)
Final_PA_means_Phylum_Long$Classifcation <- gsub('\\<POC_POU_BottledWater\\>', 'Bottled Water', Final_PA_means_Phylum_Long$Classifcation)

Final_PA_means_Phylum_Long$Classifcation <- gsub('\\<POC_All\\>', 'All Water Uses, POC', Final_PA_means_Phylum_Long$Classifcation)
Final_PA_means_Phylum_Long$Classifcation <- gsub('\\<POC_Potable\\>', 'Potable Conventional, POC', Final_PA_means_Phylum_Long$Classifcation)
Final_PA_means_Phylum_Long$Classifcation <- gsub('\\<POC_NonPotableReuse\\>', 'Non-Potable Reuse, POC', Final_PA_means_Phylum_Long$Classifcation)
Final_PA_means_Phylum_Long$Classifcation <- gsub('\\<POC_PotableReuse\\>', 'Potable Reuse, POC', Final_PA_means_Phylum_Long$Classifcation)

Final_PA_means_Phylum_Long$Classifcation <- gsub('\\<POU_All\\>', 'All Water Uses, POU', Final_PA_means_Phylum_Long$Classifcation)
Final_PA_means_Phylum_Long$Classifcation <- gsub('\\<POU_Potable\\>', 'Potable Conventional, POU', Final_PA_means_Phylum_Long$Classifcation)
Final_PA_means_Phylum_Long$Classifcation <- gsub('\\<POU_NonPotableReuse\\>', 'Non-Potable Reuse, POU', Final_PA_means_Phylum_Long$Classifcation)
Final_PA_means_Phylum_Long$Classifcation <- gsub('\\<POU_PotableReuse\\>', 'Potable Reuse, POU', Final_PA_means_Phylum_Long$Classifcation)


Final_PA_means_Phylum_Long$Classifcation <- factor(Final_PA_means_Phylum_Long$Classifcation, levels=c(
  "Potable Conventional, POC and POU", 
  "Potable Conventional, POC", 
  "Potable Conventional, POU",
  "Potable Reuse, POC and POU", 
  "Potable Reuse, POC", 
  "Potable Reuse, POU",
  "Non-Potable Reuse, POC and POU", 
  "Non-Potable Reuse, POC",
  "Non-Potable Reuse, POU"))
# "Bottled Water",
# "Blanks"


#Drop categories that were not included above. 
Final_PA_means_Phylum_Long<-Final_PA_means_Phylum_Long %>% drop_na(Classifcation)

#Adjust Width of Axis
Final_PA_means_Phylum_Long$Classifcationnew = str_wrap(Final_PA_means_Phylum_Long$Classifcation, width = 10)


#Heatmap, Phylum 
Final_Join_Phylum_Long_Norare_Mod<-ggplot(Final_PA_means_Phylum_Long,aes(Classifcation,reorder(Phylum,Presence),fill=Presence))+
  geom_tile()+
  scale_fill_gradient(low="white", high="dark blue") +
  #geom_text(aes(label = round(Presence, 2))) +
  theme_ipsum()+ 
    scale_colour_Publication()+ theme_Publication()+
  theme(legend.text=element_text(size=10), legend.title=element_text(size=10))+
  theme(legend.position="bottom")+
  theme(axis.text.x = element_text(angle = 0), axis.ticks = element_blank())+ #, vjust = 0, hjust = 0
  scale_x_discrete(labels = wrap_format(10))+
  ylab("Taxanomic Classification, Phylum")
Final_Join_Phylum_Long_Norare_Mod
leg_phylum <- get_legend(Final_Join_Phylum_Long_Norare_Mod)
as_ggplot(leg_phylum)


####Consider, Renaming X axis variables / column names. 
#Final_PA_means_Genus_Long$Classifcation <- gsub('POC_POU_All', 'All Water Uses, POC and POU', Final_PA_means_Genus_Long$Classifcation)
Final_PA_means_Genus_Long$Classifcation <- gsub('\\<POC_POU_Potable\\>', 'Potable Conventional, POC and POU', Final_PA_means_Genus_Long$Classifcation)
Final_PA_means_Genus_Long$Classifcation <- gsub('\\<POC_POU_NonPotableReuse\\>', 'Non-Potable Reuse, POC and POU', Final_PA_means_Genus_Long$Classifcation)
Final_PA_means_Genus_Long$Classifcation <- gsub('\\<POC_POU_PotableReuse\\>', 'Potable Reuse, POC and POU', Final_PA_means_Genus_Long$Classifcation)

Final_PA_means_Genus_Long$Classifcation <- gsub('\\<POC_POU_Blank\\>', 'Blanks', Final_PA_means_Genus_Long$Classifcation)
Final_PA_means_Genus_Long$Classifcation <- gsub('\\<POC_POU_BottledWater\\>', 'Bottled Water', Final_PA_means_Genus_Long$Classifcation)

Final_PA_means_Genus_Long$Classifcation <- gsub('\\<POC_All\\>', 'All Water Uses, POC', Final_PA_means_Genus_Long$Classifcation)
Final_PA_means_Genus_Long$Classifcation <- gsub('\\<POC_Potable\\>', 'Potable Conventional, POC', Final_PA_means_Genus_Long$Classifcation)
Final_PA_means_Genus_Long$Classifcation <- gsub('\\<POC_NonPotableReuse\\>', 'Non-Potable Reuse, POC', Final_PA_means_Genus_Long$Classifcation)
Final_PA_means_Genus_Long$Classifcation <- gsub('\\<POC_PotableReuse\\>', 'Potable Reuse, POC', Final_PA_means_Genus_Long$Classifcation)

Final_PA_means_Genus_Long$Classifcation <- gsub('\\<POU_All\\>', 'All Water Uses, POU', Final_PA_means_Genus_Long$Classifcation)
Final_PA_means_Genus_Long$Classifcation <- gsub('\\<POU_Potable\\>', 'Potable Conventional, POU', Final_PA_means_Genus_Long$Classifcation)
Final_PA_means_Genus_Long$Classifcation <- gsub('\\<POU_NonPotableReuse\\>', 'Non-Potable Reuse, POU', Final_PA_means_Genus_Long$Classifcation)
Final_PA_means_Genus_Long$Classifcation <- gsub('\\<POU_PotableReuse\\>', 'Potable Reuse, POU', Final_PA_means_Genus_Long$Classifcation)


Final_PA_means_Genus_Long$Classifcation <- factor(Final_PA_means_Genus_Long$Classifcation, levels=c(
  "Potable Conventional, POC and POU", 
  "Potable Conventional, POC", 
  "Potable Conventional, POU",
  "Potable Reuse, POC and POU", 
  "Potable Reuse, POC", 
  "Potable Reuse, POU",
  "Non-Potable Reuse, POC and POU", 
  "Non-Potable Reuse, POC",
  "Non-Potable Reuse, POU"))
# "Bottled Water",
# "Blanks"

#Drop categories that were not included above. 
Final_PA_means_Genus_Long<-Final_PA_means_Genus_Long %>% drop_na(Classifcation)

#Adjust Width of Axis
Final_PA_means_Genus_Long$Classifcationnew = str_wrap(Final_PA_means_Genus_Long$Classifcation, width = 10)


#Heatmap, Genus 
Final_Join_Genus_Long_Norare_Mod<-ggplot(Final_PA_means_Genus_Long,aes(Classifcation,reorder(Genus,Presence),fill=Presence))+
  geom_tile()+
  scale_fill_gradient(low="white", high="dark orange") +
  #geom_text(aes(label = round(Presence, 2))) +
  theme_ipsum()+ 
  scale_colour_Publication()+ theme_Publication()+
  theme(legend.text=element_text(size=10), legend.title=element_text(size=10))+
  theme(legend.position="bottom")+
  theme(axis.text.x = element_text(angle = 0), axis.ticks = element_blank())+ #, vjust = 0, hjust = 0
  scale_x_discrete(labels = wrap_format(10))+
  ylab("Taxanomic Classification, Genus")
Final_Join_Genus_Long_Norare_Mod
leg_Genus <- get_legend(Final_Join_Genus_Long_Norare_Mod)
as_ggplot(leg_Genus)



##Export
ggsave("Final_Join_Phylum_Long_Norare_Mod.png",plot = Final_Join_Phylum_Long_Norare_Mod,path = "C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Core/Figures/", width = 20, height = 16, units = "in", dpi = 300)
ggsave("Final_Join_Genus_Long_Norare_Mod.png",plot = Final_Join_Genus_Long_Norare_Mod,path = "C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Core/Figures/", width = 20, height = 16, units = "in", dpi = 300)

#Combo Figure
Manuscript_Combo <- ggarrange(Final_Join_Phylum_Long_Norare_Mod,
                              Final_Join_Genus_Long_Norare_Mod,
                              labels = c("A", "B"),
                              font.label = list(size = 12, color = "black"),
                              ncol = 2, nrow = 1,
                              legend = "none")
Manuscript_Combo

ggsave("Manuscript_Combo_2mod.png",plot = Manuscript_Combo,path = "C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Core/Figures/", width = 18, height = 16, units = "in", dpi = 300)
ggsave("leg_genus.png",plot = leg_genus,path = "C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Core/Figures/", width = 6, height = 2, units = "in", dpi = 300)
ggsave("leg_phylum.png",plot = leg_phylum,path = "C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Core/Figures/", width = 6, height = 2, units = "in", dpi = 300)

#Export Excels 
write.xlsx(Final_PA_means_Phylum,"C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Core/Figures/Final_PA_means_Phylum_Mod_Norare.xlsx")
write.xlsx(Final_PA_means_Genus,"C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Core/Figures/Final_PA_means_Genus_Mod_Norare.xlsx")




