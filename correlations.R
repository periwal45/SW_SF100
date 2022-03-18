setwd("/Users/periwal/ShikiFactory/WP3/SW_SF100/Data")
library(dplyr)
library(tibble)
library(reshape2)
library(janitor) # removes empty rows and columns
library(ggfortify)
library(tidyr)
library(stringr)
library(RColorBrewer)
library(ggplot2)
library(ggsci) #fancy color palettes
library(Cairo) #high resolution images
library(ggforce) #facet wrapping
library(ggpubr) #reg equation
library(corrplot)
library(corrr)
library(gtools) #mixedsort
library(data.table)


########################### correlation analysis between technical replicates (fitted OD values)
fits<-read.table(file = "../Figures/fitted_values", header = TRUE, sep = '\t')
View(fits)
dim(fits)

#bug names
bug_desc<-read.table(file = "../Bug_desc", sep = "\t", header = TRUE)
head(bug_desc)
nrow(bug_desc) #26

fits$time<-as.character(fits$time)

correl_df<-data.frame(Bug_ID=character(),
                      Plate_no=character(),
                      compound=character(),
                      comparisons=character(),
                      value=double(),
                      Parameter=character(),
                      stringsAsFactors=FALSE)

bugs<-unique(fits$Bug_ID)
bugs

for(i in 1:length(bugs)){
  
  f<-fits %>% filter(Bug_ID == bugs[i])
  head(f)
  
  reps<-unique(f$Replicate_no)
  reps
  
  for(r in 1:length(reps)){
    
    r<-f %>% filter(Replicate_no == reps[r])
    r
    plates<-unique(r$Plate_no)
    plates
  
  for(j in 1:length(plates)){
    
    c<-r %>% filter(Plate_no == plates[j])
    head(c)
    c[ ,colSums(is.na(c)) == 0]
    
    #c<-na.omit(c)
    
    #print(nrow(c))
    Bug_ID<-c[1,1]
    Bug_ID
    Plate_no<-plates[j]
    Plate_no
    Replicate_no<-c[1,6]
    Replicate_no
    
    d<-c %>% arrange(col_name,time)
    head(d)
    
    compounds<-unique(c$compound)
    head(compounds)
    
    for(z in 1:length(compounds)){
      
      cmp<-d %>% filter(compound == compounds[z])
      head(cmp)
      cpd<-cmp[1,4]
      cpd
      
      mat_melt<-cmp %>% reshape2::melt() %>% reshape(idvar = c("Bug_ID","Plate_no","compound","Replicate_no","Phyla","Species","Sp_short","time","variable"), timevar = c("tech"), direction = "wide")
      head(mat_melt)
      
      mat<-mat_melt %>% filter(variable == 'fitted')
      head(mat)
      mat<-mat[,c(1:9,11,13)]
      colnames(mat)<-c("Bug_ID","Plate_no","compound","Replicate_no","time","Phyla","Species","Sp_short","variable","T2","T1")
      
      
      select_numeric_columns<-mat %>% select_if(., is.numeric)
      select_numeric_columns
      
      if(sum(select_numeric_columns$T2) != 0 & sum(select_numeric_columns$T1) != 0){
        
      correlations<-correlate(select_numeric_columns, method = "pearson")
      correlations
      corr<-reshape2::melt(correlations)
      corr<-na.omit(corr)
      nrow(corr)
      
      if(nrow(corr) > 0){
        
      corr$comp<-paste(corr$term,corr$variable,sep = "-")
      corr<-corr[,c(4,3)]
      
      combine_df<-cbind(data.frame(Bug_ID,Plate_no,Replicate_no,cpd), corr)
      print(combine_df)
      correl_df<-rbind(correl_df, data.frame(combine_df))
      
      }else{
        
        combine_df<-cbind(data.frame(Bug_ID,Plate_no,Replicate_no,cpd), comp = "T1-T2", value = 0)
      }
      }
    }
  }
  }
}


head(correl_df)
#
corr_values<-correl_df %>% filter(Bug_ID != 'NT5076')
head(corr_values)
nrow(corr_values) #34,714
corr_values<-corr_values %>% filter(!comp %in% c("T2-T1"))
nrow(corr_values) #17,357

head(corr_values)
write.table(corr_values, file = paste0("../Figures/tech_fittedOD_correl_pearson"), sep = "\t", quote = FALSE, row.names = FALSE)


# read fitted correaltions
corr_values<-read.table(file = paste0("../Figures/tech_fittedOD_correl_pearson"), sep = "\t", header = TRUE)

median_corr_bug_wise<-corr_values %>% dplyr::group_by(Bug_ID,Replicate_no,cpd) %>%
  dplyr::summarise(median_corr = median(value))

head(median_corr_bug_wise)

all_bugs_desc<-merge(median_corr_bug_wise,bug_desc, by="Bug_ID")
head(all_bugs_desc)

all_bugs_corr_desc<-merge(corr_values,bug_desc, by="Bug_ID")
View(all_bugs_corr_desc)

#plots colored histogram of all correlation values
CairoSVG(file=paste("../Figures/tech_fittedOD_correlation.svg", sep = ""), width = 7, height = 5, bg = "white")
ggplot(all_bugs_corr_desc, aes(value)) + geom_histogram(aes(fill=Sp_short), color="black", lwd=0.2) + 
  theme_bw() + th + scale_fill_igv() + scale_x_continuous(name = "correlation") +
  guides(fill=guide_legend(title="Species"))
dev.off()


# noisy wells
CairoSVG(file=paste("../Figures/technoise_fittedOD_correlation.svg", sep = ""), width = 7, height = 5, bg = "white")
View(all_bugs_corr_desc %>% filter(value < 0.8)) #%>% 
  ggplot(aes(value)) + geom_histogram(aes(fill=Sp_short), color="black", lwd=0.2) + 
  theme_bw() + th + scale_fill_uchicago() + scale_x_continuous(name = "correlation") +
  guides(fill=guide_legend(title="Species"))
dev.off()


########################### correlation analysis between batch replicates (fitted OD values)
fits<-read.table(file = "../Figures/fitted_values", header = TRUE, sep = '\t')
head(fits)
dim(fits)

fits$time<-as.character(fits$time)

mat_melt<-fits %>% reshape2::melt() %>% reshape(idvar = c("Bug_ID","Plate_no","col_name","compound","tech","Phyla","Species","Sp_short","time","variable"), timevar = c("Replicate_no"), direction = "wide")
head(mat_melt)
nrow(mat_melt)
colnames(mat_melt)<-c("Bug_ID","Plate_no","col_name","compound","time","tech","Phyla","Species","Sp_short","variable","Rep2","Rep3","Rep1")

#write.table(mat_melt, file = "../Figures/mat", sep = "\t", quote = FALSE, row.names = FALSE)

#Select variable to compute correlation
mat<-mat_melt %>% filter(variable == 'fitted')
head(mat)
nrow(mat) #232,272

mat<-mat %>% filter(Bug_ID != 'NT5076')
bugs<-unique(mat$Bug_ID)
bugs

correl_df_batch<-data.frame(Bug_ID=character(),
                      Plate_no=character(),
                      comparisons=character(),
                      value=double(),
                      Parameter=character(),
                      stringsAsFactors=FALSE)

for(i in 1:length(bugs)){
  
  f<-mat %>% filter(Bug_ID == bugs[i])
  head(f)
  
  plates<-unique(f$Plate_no)
  plates
  
  for(j in 1:length(plates)){
    
    c<-f %>% filter(Plate_no == plates[j])
    head(c)
    c[ ,colSums(is.na(c)) == 0]
    #c<-na.omit(c)
    
    #print(nrow(c))
    Bug_ID<-c[1,1]
    Bug_ID
    Plate_no<-plates[j]
    Plate_no
    
    d<-c %>% arrange(col_name,time)
    head(d)
    
    select_numeric_columns<-d %>% select_if(., is.numeric)
    select_numeric_columns
    
    correlations<-correlate(select_numeric_columns, method = "pearson")
    correlations
    corr<-reshape2::melt(correlations)
    corr<-na.omit(corr)
    corr
    corr$comp<-paste(corr$term,corr$variable,sep = "-")
    corr<-corr[,c(4,3)]
    
    combine_df<-cbind(data.frame(Bug_ID,Plate_no), corr)
    print(combine_df)
    correl_df_batch<-rbind(correl_df_batch, data.frame(combine_df))
    
  }
  
}

head(correl_df_batch)
#
nrow(correl_df_batch) #828
correl_df_batch<-correl_df_batch %>% filter(!comp %in% c("Rep3-Rep2","Rep2-Rep1","Rep3-Rep1"))
nrow(correl_df_batch) #414

head(correl_df_batch)
write.table(correl_df_batch, file = paste0("../Figures/batch_fittedOD_correl_pearson"), sep = "\t", quote = FALSE, row.names = FALSE)

median_corr_bug_wise<-correl_df_batch %>% dplyr::group_by(Bug_ID) %>%
  dplyr::summarise(median_corr = median(value))

View(median_corr_bug_wise)

#write.table(median_corr_bug_wise, file = paste0("/Users/vinitaperiwal/GrowthCurver/median_correl_pearson"), sep = "\t", quote = FALSE, row.names = FALSE)

all_bugs_desc<-merge(median_corr_bug_wise,bug_desc, by="Bug_ID")
head(all_bugs_desc)

all_bugs_corr_desc<-merge(correl_df_batch,bug_desc, by="Bug_ID")
head(all_bugs_corr_desc)

#plots colored histogram of all correlation values
CairoSVG(file=paste("../Figures/batch_fittedOD_correlation.svg", sep = ""), width = 7, height = 5, bg = "white")
ggplot(all_bugs_corr_desc, aes(value)) + geom_histogram(aes(fill=Sp_short), color="black", lwd=0.2) + 
  theme_bw() + th + scale_fill_igv() + scale_x_continuous(name = "correlation") +
  guides(fill=guide_legend(title="Species"))
dev.off()

#plots median correlation values of plates of each bug
CairoSVG(file=paste("../Figures/batch_median_correlation.svg", sep = ""), width = 4, height = 4, bg = "white")
all_bugs_desc %>% dplyr::group_by(Bug_ID,Sp_short) %>%
  ggplot(aes(x=Sp_short,y=median_corr)) + geom_bar(stat = "identity") + theme_bw() + th +
  theme(axis.text.x = element_text(angle = 90, size = 10, hjust = 1)) + scale_x_discrete(name = "Species")
dev.off()
