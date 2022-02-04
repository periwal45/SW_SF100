setwd("/Users/periwal/ShikiFactory/WP3/SF100-WP3/GrowthProfiler/")
library(dplyr)
library(tibble)
library(reshape2)
library(janitor) # removes empty rows and columns
library(ggfortify)
library(tidyr)
library(stringr)
library(RColorBrewer)
library(ggplot2)
library(lattice)
library(gridExtra)
library(ggsci) #fancy color palettes
library(Cairo) #high resolution images
library(ggforce) #facet wrapping
library(ggpubr) #reg equation
library(metRology) # t distribution
library(cluster) #autoplot cluster
library(broom)
library(car) #for Boxplots with labels
library(compare)
library(corrplot)
library(Hmisc)
library(pracma)
library(dunn.test)
library(purrr)
library(rcompanion) #cldList
library(rstatix) #KW test and others
library(corrr)
library(psych)
library(pROC) #auc function in growthcurver
library(growthcurver)
library(MASS) #fitdistr
library(gtools) #mixedsort
library(metap) #sumlog
library(esc) 
library(qvalue) #estimates FDR from a list of input p-values
library(fitdistrplus)
library(data.table)
library(BSDA)
library(bigutilsr) #tukey_mc_up outlier detection
library(forestplot)
library(meta)

# source the script with functions
source("/Users/periwal/ShikiFactory/WP3/SF100-WP3/Functions.R")

#read all files with .tab extension
file.names <- list.files(path = '/Users/periwal/ShikiFactory/WP3/SF100-WP3/GrowthProfiler/', recursive = TRUE, pattern = "\\.tab$") #recursive reads through all subfolders and files
file.names

annot<-data.frame(read.table(file = "/Users/periwal/ShikiFactory/WP3/SF100-WP3/plate_layout", sep = '\t', header = TRUE))
annot<-annot %>% reshape2::melt(measure.vars = c("T1","T2"), value.name = "technical_rep")
head(annot)

#loop for each .tab file, annotates
for(i in 1:length(file.names)){
  
  #read file/table
  input_file<-tools::file_path_sans_ext(file.names[i]) #read filename w/o extension
  input_file
  
  #create a new dataframe with additional annotations for each file
  
  annot <- tibble(Bug_ID = create_annot(input_file)[1], Replicate_no = create_annot(input_file)[2], Plate_no = create_annot(input_file)[3])
  annot
  X<-data.frame(read.table(file.names[i], header = TRUE, sep = '\t', stringsAsFactors = FALSE))
  
  head(X)
  colnames(X)[1]<-"time"
  X$time<-as.integer(X$time)
  
  annot_file<-merge(annot,X)
  head(annot_file)
  
  #write annotated files
  annot_outfile<-paste0(input_file,".annot")
  write.table(annot_file, file = annot_outfile, sep = "\t", quote = FALSE, row.names = FALSE)
  
}


#merge all replicates of a bug
# Merge files across all replicates of a strain (.annot files)
dir.names<-list.dirs(path = '/Users/periwal/ShikiFactory/WP3/SF100-WP3/GrowthProfiler', recursive = FALSE)
dir.names

for(d in 1:length(dir.names)){
  
  if(str_detect(dir.names[d], "/Users/periwal/ShikiFactory/WP3/SF100-WP3/GrowthProfiler/NT")){
    
    all_rep<-do.call(rbind, lapply(a<-list.files(path=dir.names[d], pattern="\\.annot$", full.names = TRUE), function(i){
      read.table(i,header=TRUE, sep="\t", stringsAsFactors = FALSE)}))
    write.table(all_rep, file = paste0(dir.names[d],"_all_replicates.merged"), sep = "\t", quote = FALSE, row.names = FALSE)
    
  }
  
}

#read all files with .merged extension
file.merged<-list.files(path = '/Users/periwal/ShikiFactory/WP3/SF100-WP3/GrowthProfiler', pattern = ".merged$") #recursive reads through all subfolders and files
file.merged

#create an empty data frame to store all raw points and fitted points
fitted_values<-data.frame(stringsAsFactors = FALSE)

#create an empty data frame to store all computed model parameters
gc_fit_params<-data.frame(stringsAsFactors = FALSE)

bug_min_reads<-data.frame(Bug_ID = character(), Replicate_no = character(), count_tp = double(), stringsAsFactors = FALSE)

#loop for each file, fits logistic curve, annotates
for(j in 1:length(file.merged)){
  
  #read file/table
  in_file<-tools::file_path_sans_ext(file.merged[j]) #read filename w/o extension
  in_file
  
  Y<-data.frame(read.table(file.merged[j], header = TRUE, sep = '\t', stringsAsFactors = FALSE))
  Y[,4]<-Y[,4]/60
  head(Y)
  
  #count time points in all replicates, all should have same time point reading else trim if needed
  
  time_points<-Y %>% dplyr::group_by(Bug_ID,Replicate_no) %>%
    summarise(count_tp = max(time))
  head(time_points)
  time_points<-time_points %>% mutate(min_tp = min(count_tp))
  min_tp<-min(time_points$count_tp)
  min_tp
  
  bug_min_reads<-rbind(bug_min_reads, data.frame(time_points))
  bug_min_reads
  
  Z<-Y %>% dplyr::group_by(Bug_ID,Replicate_no,Plate_no) %>%
    filter(time <= min_tp)
  
  head(Z)
  nrow(Z)
  
  write.table(Z, file = paste0("/Users/periwal/ShikiFactory/WP3/SF100-WP3/GrowthProfiler/trim_points/",in_file,".trim"), sep = "\t", quote = FALSE, row.names = FALSE)
  
  unique_rep<-unique(Z$Replicate_no)
  unique_rep
  
  unique_plate<-unique(Z$Plate_no)
  unique_plate
  
  #loop plate wise and then well wise (col_name) to fit growth curve on each well, later merge all data in a data frame
  for(j in 1:length(unique_plate)){
    
    dat<-Z %>% filter(Plate_no == unique_plate[j])
    dat
    
    if (dim(dat)[1] != 0){ #if missing plate
      
      for(k in 1:length(unique_rep)){
        
        plate_dat<-dat %>% filter(Replicate_no == unique_rep[k])
        plate_dat
        
        if (dim(plate_dat)[1] != 0){ #if missing replicate eg NT5022, plate16 has no rep1
          
          bug<-plate_dat[1,1]
          bug
          rep<-plate_dat[1,2]
          rep
          plate<-plate_dat[1,3]
          plate
          #drug<-plate_dat[1,4]
          #drug
          
          #plot_file<-paste0(ID,"_",rep,"_",plate)
          #plot_file
          
          data<-plate_dat[,-(1:3)]
          data

          #fit model to each well of a plate and save the fitted model
          for(col_name in names(data)){
            
            col_name
            
            if(col_name != "time"){
              
              current_well<-data[, c("time",col_name)]
              current_well
              min_OD<-min(current_well[, col_name])
              min_OD
              
              #do background correction using min OD of each well
              current_well[,col_name]<-current_well[,col_name] - min_OD
              current_well[,col_name]
              
              #each time create a new variable for each well
              gc_fit<-SummarizeGrowth(data_t = current_well[,"time"], data_n = current_well[,col_name])
              #saveRDS(gc_fit, file = paste0("/Users/vinitaperiwal/GrowthCurver/models/",ID,"_",rep,"_",plate,"_",drug,"_",col_name,".rds"))
              
              gc_fit
              
              #create a data frame of raw values and fitted values
              mod_t<-data.frame(matrix(unlist(gc_fit$data$t)))
              mod_t
              mod_N<-data.frame(matrix(unlist(gc_fit$data$N)))
              mod_N
              
              if(gc_fit$vals$k != 0 & gc_fit$vals$n0 != 0 & gc_fit$vals$r != 0){
                
                mod<-cbind(mod_t, mod_N, gc_fit$model$m$fitted()) #m is the model object with all fitted and residual values
                
              }else{
                
                mod<-cbind(mod_t, mod_N, gc_fit$vals$r)
                
              }
              
              colnames(mod)<-c("time","OD","fitted")
              mod
              
              #add annotation to each row
              annot_mod<-cbind(bug,rep,plate,col_name,mod)
              annot_mod
              
              fitted_values<-rbind(fitted_values, annot_mod)
              fitted_values
              
              fitted_param<-data.frame(cbind(bug,rep,plate,col_name,gc_fit$vals$k,gc_fit$vals$k_se,gc_fit$vals$k_p,gc_fit$vals$n0,gc_fit$vals$n0_se,gc_fit$vals$n0_p,gc_fit$vals$r,gc_fit$vals$r_se,gc_fit$vals$r_p,gc_fit$vals$sigma,gc_fit$vals$df,gc_fit$vals$t_mid,gc_fit$vals$t_gen,gc_fit$vals$auc_l,gc_fit$vals$auc_e,gc_fit$vals$note))
              colnames(fitted_param)<-c("Bug_ID","Replicate_no","Plate_no","well","k","k_se","k_p","n0","n0_se","n0_p","r","r_se","r_p","sigma","df","t_mid","t_gen","auc_l","auc_e","note")
              fitted_param
              
              gc_fit_params<-rbind(gc_fit_params, fitted_param)
              tail(gc_fit_params)
              
            }
          }
          
        }
        
      }
    }
    
  }
  
}


head(fitted_values)
nrow(fitted_values) #691,200

View(fitted_values)

fitted_values[fitted_values=="Plate1-control"]<-'Plate1-GALCH'
fitted_values[fitted_values=="Plate2-control"]<-'Plate2-Comm'


#mean_fit_GP<-fitted_values %>% filter(`file.names[i]` == "OD/SF100/SF100_all_bugs_growthtest_MTP02_OD_GPplate.tab") %>% dplyr::group_by(bug,strain,time) %>% summarise(avOD = mean(OD), avfit = mean(fitted))
#View(mean_fit_GP)

head(gc_fit_params)
nrow(gc_fit_params) #should be equal to number of wells: 96 for 1 plate

gc_fit_params[gc_fit_params=="Plate1-control"]<-'Plate1-GALCH'
gc_fit_params[gc_fit_params=="Plate2-control"]<-'Plate2-Comm'


write.table(fitted_values, file = "/Users/periwal/ShikiFactory/WP3/SF100-WP3/GrowthProfiler/GC_fitted_model_values", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gc_fit_params, file = "/Users/periwal/ShikiFactory/WP3/SF100-WP3/GrowthProfiler/GC_fit_params", sep = "\t", quote = FALSE, row.names = FALSE)


View(gc_fit_params)
############################### Read fitted values
#read fitted values of models
fits<-data.frame(read.table("/Users/periwal/ShikiFactory/WP3/SF100-WP3/GrowthProfiler/GC_fitted_model_values", header = TRUE, sep = '\t', stringsAsFactors = FALSE))
head(fits)
nrow(fits) #691,200

plate_layout<-data.frame(read.table("/Users/periwal/ShikiFactory/WP3/SF100-WP3/plate_layout", header = TRUE, sep = '\t', stringsAsFactors = FALSE))
head(plate_layout)
nrow(plate_layout) #4800
plate_layout<-plate_layout %>% reshape2::melt(measure.vars = c("T1","T2"))
head(plate_layout)
nrow(plate_layout) #9600
names(plate_layout)[names(plate_layout) == "variable"] <- "technical"
names(plate_layout)[names(plate_layout) == "value"] <- "col_name"

head(plate_layout)

merge_plate_layout_fits<-merge(plate_layout, fits, by=c("Bug_ID","Plate_no","col_name"))
nrow(merge_plate_layout_fits) #691,200

View(merge_plate_layout_fits)

# single well bug-only plot

head(fits)

CairoSVG(file=paste("/Users/periwal/ShikiFactory/WP3/SF100-WP3/Figures/NT5002.svg"), width = 3, height = 2, bg = "white")
fits %>% filter(Bug_ID == 'NT5002' & col_name == 'H4' & Replicate_no == 'Rep1' & Plate_no == 'Plate1-GALCH') %>%
  ggplot(aes(x=time, y=fitted,group=Replicate_no)) + geom_line() + 
  geom_point(aes(y=OD), size=0.001) + scale_color_aaas() + 
  ggtitle(paste('B. fragilis')) + scale_x_continuous(name = "Time (Hours)", limits = c(0,25)) + scale_y_continuous(name="bg corrected OD") + 
  th + theme_minimal()
dev.off()

# rename/remove wells modified during pipetting
# Plate2+C, NT5003, Rep1, well A7-A12, only DMSO
for(i in 1:nrow(merge_plate_layout_fits)){
  
  line = merge_plate_layout_fits[i,]
  #print(line)
  
  if(line$Bug_ID == 'NT5003' & line$Replicate_no == "Rep1" & line$Plate_no == "Plate2+C" & line$col_name %in% c('A7','A8','A9','A10','A11','A12')){
    
    #line$Compound<-'DMSO'
    merge_plate_layout_fits[i,]$Compound<-'DMSO'
    
  }
}

nrow(merge_plate_layout_fits) #691,200

#NT5002, Rep2, Plate2+D, Column 6 no bug - delete, column 7 has double bug

merge_plate_layout_fits<-merge_plate_layout_fits %>% filter(!(Bug_ID == 'NT5002' & Replicate_no == "Rep2" & Plate_no == "Plate2+D" & 
                                        col_name %in% c('A6','B6','C6','D6','E6','F6','G6','H6','A7','B7','C7','D7','E7','F7','G7','H7')))

nrow(merge_plate_layout_fits)  #690,816
head(merge_plate_layout_fits) 

#NT5002, Rep3, H1-H6 has Vanillin, H7-H12 is DMSO

for(i in 1:nrow(merge_plate_layout_fits)){
  
  line = merge_plate_layout_fits[i,]
  #print(line)
  
  if(line$Bug_ID == 'NT5002' & line$Replicate_no == "Rep3" & line$Plate_no %in% c("Plate1+V","Plate2+V") & line$col_name %in% c('H7','H8','H9','H10','H11','H12')){
    
    #line$Compound<-'DMSO'
    merge_plate_layout_fits[i,]$Compound<-'DMSO'
    
  }else if(line$Bug_ID == 'NT5002' & line$Replicate_no == "Rep3" & line$Plate_no %in% c("Plate1+V","Plate2+V") & line$col_name %in% c('H1','H2','H3','H4','H5','H6')){
    
    merge_plate_layout_fits[i,]$Compound<-'Vanillin'
    
  }
}

head(merge_plate_layout_fits)
nrow(merge_plate_layout_fits) #690,816

#NT5037, Rep2, Plate1-GALCH (row E and F are swapped during plate prep)

x<-merge_plate_layout_fits
head(x)
nrow(x) #690,816

for(i in 1:nrow(x)){
  
  line = x[i,]
  #print(line)
  
  if(line$Bug_ID == 'NT5037' & line$Replicate_no == "Rep2" & line$Plate_no == 'Plate1-GALCH'){
    
    x[i,]$col_name<-ifelse(x[i,]$col_name=="^E", "F", ifelse(x[i,]$col_name=="^F", "E", x[i,]$col_name))
    x[i,]$Compound<-ifelse(x[i,]$Compound=="MIN-2-073", "MIN-1-124", ifelse(x[i,]$Compound=="MIN-1-124", "MIN-2-073", x[i,]$Compound))
    x[i,]$Compound<-ifelse(x[i,]$Compound=="MIN-2-107", "MIN-2-013a", ifelse(x[i,]$Compound=="MIN-2-013a", "MIN-2-107", x[i,]$Compound))
    x[i,]$Compound<-ifelse(x[i,]$Compound=="MIN-3-076", "MIN-3-073", ifelse(x[i,]$Compound=="MIN-3-073", "MIN-3-076", x[i,]$Compound))
    x[i,]$Compound<-ifelse(x[i,]$Compound=="MIN-3-064", "MIN-3-066", ifelse(x[i,]$Compound=="MIN-3-066", "MIN-3-064", x[i,]$Compound))
    x[i,]$Compound<-ifelse(x[i,]$Compound=="MIN-3-029", "MIN-3-026", ifelse(x[i,]$Compound=="MIN-3-026", "MIN-3-029", x[i,]$Compound))
    x[i,]$Compound<-ifelse(x[i,]$Compound=="MIN-3-053", "MIN-2-014", ifelse(x[i,]$Compound=="MIN-2-014", "MIN-3-053", x[i,]$Compound))
    
  }
}

head(x)
nrow(x)

head(merge_plate_layout_fits)
merge_plate_layout_fits<-x

nrow(merge_plate_layout_fits) #690,816

#NT5037, Rep2, Plate1+A and plate2+A, controls swapped

for(i in 1:nrow(merge_plate_layout_fits)){
  
  line = merge_plate_layout_fits[i,]
  #print(line)
  
  if(line$Bug_ID == 'NT5037' & line$Replicate_no == "Rep2" & line$Plate_no %in% c("Plate1+A","Plate2+A") & line$col_name %in% c('H7','H8','H9','H10','H11','H12')){
    
    #line$Compound<-'DMSO'
    merge_plate_layout_fits[i,]$Compound<-'DMSO'
    
  }else if(line$Bug_ID == 'NT5037' & line$Replicate_no == "Rep2" & line$Plate_no %in% c("Plate1+A","Plate2+A") & line$col_name %in% c('H1','H2','H3','H4','H5','H6')){
    
    merge_plate_layout_fits[i,]$Compound<-'Aspartame'
    
  }
}

nrow(merge_plate_layout_fits) #690,816

#NT24007 and NT5019, Rep2 and 3, Plate1+A and plate2+A, controls swapped

for(i in 1:nrow(merge_plate_layout_fits)){
  
  line = merge_plate_layout_fits[i,]
  #print(line)
  
  if(line$Bug_ID == 'NT5037' & line$Replicate_no == "Rep2" & line$Plate_no %in% c("Plate1+A","Plate2+A") & line$col_name %in% c('H7','H8','H9','H10','H11','H12')){
    
    #line$Compound<-'DMSO'
    merge_plate_layout_fits[i,]$Compound<-'DMSO'
    
  }else if(line$Bug_ID == 'NT5037' & line$Replicate_no == "Rep2" & line$Plate_no %in% c("Plate1+A","Plate2+A") & line$col_name %in% c('H1','H2','H3','H4','H5','H6')){
    
    merge_plate_layout_fits[i,]$Compound<-'Aspartame'
    
  }
}

nrow(merge_plate_layout_fits) #690816

write.table(merge_plate_layout_fits, file = "/Users/periwal/ShikiFactory/WP3/SF100-WP3/GrowthProfiler/GC_fitted_model_values_filtered", sep = "\t", quote = FALSE, row.names = FALSE)

############### Merge technical replicates (plot raw and fitted ODs)
merge_plate_layout_fits<-data.frame(read.table("/Users/periwal/ShikiFactory/WP3/SF100-WP3/GrowthProfiler/GC_fitted_model_values_filtered", header = TRUE, sep = '\t', stringsAsFactors = FALSE))
head(merge_plate_layout_fits)

techs<-merge_plate_layout_fits %>% dplyr::group_by(Bug_ID,Plate_no,Replicate_no,Compound,time) %>%
  summarise(avg_OD = median(OD), avg_fit = median(fitted))
head(techs)
nrow(techs) #312,408

#bug names
bug_desc<-read.table(file = "/Users/periwal/ShikiFactory/WP3/SF100-WP3/Bug_desc", sep = "\t", header = TRUE)
head(bug_desc)
nrow(bug_desc) #10

#merge techs and bug desc
techs_bugs<-merge(techs,bug_desc, by="Bug_ID")
head(techs_bugs)

### plot plate wise replicates for each bug

bugs<-unique(techs_bugs$Bug_ID)
bugs

for(i in 1:length(bugs)){
  
  f<-techs_bugs %>% filter(Bug_ID == bugs[i])
  
  plates<-unique(f$Plate_no)
  plates
  
  for(j in 1:length(plates)){
    
    p<-f %>% filter(Plate_no == plates[j])
    head(p)
    
    CairoSVG(file=paste("/Users/periwal/ShikiFactory/WP3/SF100-WP3/Figures/", bugs[i], "_", plates[j], "_fil.svg", sep = ""), width = 12, height = 8, bg = "white")
    print(p %>% dplyr::group_by(Bug_ID,Plate_no,Compound,Replicate_no) %>%
            ggplot(aes(x=time, y=avg_fit,group=Replicate_no)) + geom_line(aes(color=Replicate_no)) + 
            geom_point(aes(y=avg_OD,color=Replicate_no), size=0.001) + facet_wrap("Compound", scales = "free", ncol = 6, nrow = 8) + scale_color_aaas() + 
            ggtitle(paste(bugs[i],p$Sp_short,plates[j],sep=":")) + scale_x_continuous(name = "Time (Hours)", limits = c(0,25)) + scale_y_continuous(name="bg corrected OD") + 
            th + theme_minimal())
    dev.off()
    
  }
}




########################### correlation analysis between replicates (fitted OD values)
merge_plate_layout_fits<-read.table(file = "/Users/periwal/ShikiFactory/WP3/SF100-WP3/GrowthProfiler/GC_fitted_model_values_filtered", header = TRUE, sep = '\t')
View(merge_plate_layout_fits)

#bug names
bug_desc<-read.table(file = "/Users/periwal/ShikiFactory/WP3/SF100-WP3/Bug_desc", sep = "\t", header = TRUE)
head(bug_desc)
nrow(bug_desc) #10

#merge techs and bug desc
fits<-merge(merge_plate_layout_fits,bug_desc, by="Bug_ID")
View(fits)

summary(fits)
fits<-fits[,c(1:4,6:12)]
fits$time<-as.character(fits$time)

mat_melt<-fits %>% reshape2::melt() %>% reshape(idvar = c("Bug_ID","Plate_no","col_name","Compound","Phyla","Species","Sp_short","time","variable"), timevar = c("Replicate_no"), direction = "wide")
head(mat_melt)
nrow(mat_melt)
colnames(mat_melt)<-c("Bug_ID","Plate_no","col_name","Compound","time","Phyla","Species","Sp_short","variable","Rep1","Rep2","Rep3")

write.table(mat_melt, file = "/Users/periwal/ShikiFactory/WP3/SF100-WP3/Figures/mat", sep = "\t", quote = FALSE, row.names = FALSE)

#Select variable to compute correlation
mat<-mat_melt %>% filter(variable == 'fitted')
head(mat)
nrow(mat) #232,272

bugs<-unique(mat$Bug_ID)
bugs

correl_df<-data.frame(Bug_ID=character(),
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
    c<-na.omit(c)
    
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
    correl_df<-rbind(correl_df, data.frame(combine_df))
    
  }
  
}

head(correl_df)
#
corr_values<-correl_df
head(corr_values)
nrow(corr_values) #600
corr_values<-corr_values %>% filter(!comp %in% c("Rep3-Rep2","Rep2-Rep1","Rep3-Rep1"))
nrow(corr_values) #300

head(corr_values)
write.table(corr_values, file = paste0("/Users/periwal/ShikiFactory/WP3/SF100-WP3/Figures/fittedOD_correl_pearson"), sep = "\t", quote = FALSE, row.names = FALSE)

median_corr_bug_wise<-corr_values %>% dplyr::group_by(Bug_ID) %>%
  dplyr::summarise(median_corr = median(value))

median_corr_bug_wise

#write.table(median_corr_bug_wise, file = paste0("/Users/vinitaperiwal/GrowthCurver/median_correl_pearson"), sep = "\t", quote = FALSE, row.names = FALSE)

#bug names
bug_desc<-read.table(file = "/Users/periwal/ShikiFactory/WP3/SF100-WP3/Bug_desc", sep = "\t", header = TRUE)
head(bug_desc)
nrow(bug_desc)
#
all_bugs_desc<-merge(median_corr_bug_wise,bug_desc, by="Bug_ID")
head(all_bugs_desc)

all_bugs_corr_desc<-merge(corr_values,bug_desc, by="Bug_ID")
head(all_bugs_corr_desc)

#plots colored histogram of all correlation values
CairoSVG(file=paste("/Users/periwal/ShikiFactory/WP3/SF100-WP3/Figures/fittedOD_correlation.svg", sep = ""), width = 6, height = 4, bg = "white")
ggplot(all_bugs_corr_desc, aes(value)) + geom_histogram(aes(fill=Sp_short), color="white") + theme_bw() + th + scale_fill_aaas()
dev.off()

#plots median correlation values of plates of each bug
CairoSVG(file=paste("/Users/periwal/ShikiFactory/WP3/SF100-WP3/Figures/median_correlation.svg", sep = ""), width = 4, height = 4, bg = "white")
all_bugs_desc %>% dplyr::group_by(Bug_ID,Sp_short) %>%
  ggplot(aes(x=Sp_short,y=median_corr)) + geom_bar(stat = "identity") + theme_bw() + th +
  theme(axis.text.x = element_text(angle = 90, size = 10, hjust = 1)) + scale_x_discrete(name = "Species")
dev.off()

#######read parameter values
params<-data.frame(read.table("/Users/periwal/ShikiFactory/WP3/SF100-WP3/GrowthProfiler/GC_fit_params", header = TRUE, sep = '\t', stringsAsFactors = FALSE))
head(params)

nrow(params) #28,800

#merge params and bug desc
params_bugs<-merge(params,bug_desc, by="Bug_ID")
head(params_bugs)
nrow(params_bugs)
#names(params_bugs)[names(params_bugs) == "well"] <- "col_name"

#merge plate_layout and params_bugs
params_bugs_cmp<-merge(plate_layout, params_bugs, by=c("Bug_ID","Plate_no","col_name"))
head(params_bugs_cmp)
nrow(params_bugs_cmp)

# rename/remove wells modified during pipetting
# Plate2+C, NT5003, Rep1, well A7-A12, only DMSO
for(i in 1:nrow(params_bugs_cmp)){
  
  line = params_bugs_cmp[i,]
  #print(line)
  
  if(line$Bug_ID == 'NT5003' & line$Replicate_no == "Rep1" & line$Plate_no == "Plate2+C" & line$col_name %in% c('A7','A8','A9','A10','A11','A12')){
    
    #line$Compound<-'DMSO'
    params_bugs_cmp[i,]$Compound<-'DMSO'
    
  }
}

nrow(params_bugs_cmp) #28800
#View(params_bugs_cmp)

#NT5002, Rep2, Plate2+D, Column 6 no bug - delete, column 7 has double bug

params_bugs_cmp<-params_bugs_cmp %>% filter(!(Bug_ID == 'NT5002' & Replicate_no == "Rep2" & Plate_no == "Plate2+D" & 
                                                                col_name %in% c('A6','B6','C6','D6','E6','F6','G6','H6','A7','B7','C7','D7','E7','F7','G7','H7')))

nrow(params_bugs_cmp)  #28,784

#NT5002, Rep3, H1-H6 has Vanillin, H7-H12 is DMSO

for(i in 1:nrow(params_bugs_cmp)){
  
  line = params_bugs_cmp[i,]
  #print(line)
  
  if(line$Bug_ID == 'NT5002' & line$Replicate_no == "Rep3" & line$Plate_no %in% c("Plate1+V","Plate2+V") & line$col_name %in% c('H7','H8','H9','H10','H11','H12')){
    
    #line$Compound<-'DMSO'
    params_bugs_cmp[i,]$Compound<-'DMSO'
    
  }else if(line$Bug_ID == 'NT5002' & line$Replicate_no == "Rep3" & line$Plate_no %in% c("Plate1+V","Plate2+V") & line$col_name %in% c('H1','H2','H3','H4','H5','H6')){
    
    params_bugs_cmp[i,]$Compound<-'Vanillin'
    
  }
}

nrow(params_bugs_cmp) #28,784

#NT5037, Rep2, Plate1-GALCH (row E and F are swapped during plate prep)

p<-params_bugs_cmp

for(i in 1:nrow(p)){
  
  line = p[i,]
  #print(line)
  
  if(line$Bug_ID == 'NT5037' & line$Replicate_no == "Rep2" & line$Plate_no == 'Plate1-GALCH'){
    
   p[i,]$col_name<-ifelse(p[i,]$col_name=="^E", "F", ifelse(p[i,]$col_name=="^F", "E", p[i,]$col_name))
   p[i,]$Compound<-ifelse(p[i,]$Compound=="MIN-2-073", "MIN-1-124", ifelse(p[i,]$Compound=="MIN-1-124", "MIN-2-073", p[i,]$Compound))
   p[i,]$Compound<-ifelse(p[i,]$Compound=="MIN-2-107", "MIN-2-013a", ifelse(p[i,]$Compound=="MIN-2-013a", "MIN-2-107", p[i,]$Compound))
   p[i,]$Compound<-ifelse(p[i,]$Compound=="MIN-3-076", "MIN-3-073", ifelse(p[i,]$Compound=="MIN-3-073", "MIN-3-076", p[i,]$Compound))
   p[i,]$Compound<-ifelse(p[i,]$Compound=="MIN-3-064", "MIN-3-066", ifelse(p[i,]$Compound=="MIN-3-066", "MIN-3-064", p[i,]$Compound))
   p[i,]$Compound<-ifelse(p[i,]$Compound=="MIN-3-029", "MIN-3-026", ifelse(p[i,]$Compound=="MIN-3-026", "MIN-3-029", p[i,]$Compound))
   p[i,]$Compound<-ifelse(p[i,]$Compound=="MIN-3-053", "MIN-2-014", ifelse(p[i,]$Compound=="MIN-2-014", "MIN-3-053", p[i,]$Compound))
  
  }
}

View(p)

params_bugs_cmp<-p

nrow(params_bugs_cmp) #28,784

#NT5037, Rep2, Plate1+A and plate2+A, controls swapped

for(i in 1:nrow(params_bugs_cmp)){
  
  line = params_bugs_cmp[i,]
  #print(line)
  
  if(line$Bug_ID == 'NT5037' & line$Replicate_no == "Rep2" & line$Plate_no %in% c("Plate1+A","Plate2+A") & line$col_name %in% c('H7','H8','H9','H10','H11','H12')){
    
    #line$Compound<-'DMSO'
    params_bugs_cmp[i,]$Compound<-'DMSO'
    
  }else if(line$Bug_ID == 'NT5037' & line$Replicate_no == "Rep2" & line$Plate_no %in% c("Plate1+A","Plate2+A") & line$col_name %in% c('H1','H2','H3','H4','H5','H6')){
    
    params_bugs_cmp[i,]$Compound<-'Aspartame'
    
  }
}

nrow(params_bugs_cmp) #28,784

write.table(params_bugs_cmp, file = "/Users/periwal/ShikiFactory/WP3/SF100-WP3/GrowthProfiler/GC_fit_params_filtered", sep = "\t", quote = FALSE, row.names = FALSE)

head(params_bugs_cmp)

nrow(params_bugs_cmp) #28,784

################### Quality control of fitted growth curves

params_bugs_cmp<-data.frame(read.table("/Users/periwal/ShikiFactory/WP3/SF100-WP3/GrowthProfiler/GC_fit_params_filtered", header = TRUE, sep = '\t', stringsAsFactors = FALSE))
head(params_bugs_cmp)

#### Plot histogram of sigma values and pca of all params
bugs<-unique(params_bugs_cmp$Bug_ID)
bugs

outlier_wells_pca<-data.frame(stringsAsFactors = FALSE)

for(i in 1:length(bugs)){
  
  f<-params_bugs_cmp %>% filter(Bug_ID == bugs[i])
  plates<-unique(f$Plate_no)
  plates
  
  for(j in 1:length(plates)){
    
    p<-f %>% filter(Plate_no == plates[j])
    head(p)
    
    CairoSVG(file=paste("/Users/periwal/ShikiFactory/WP3/SF100-WP3/Figures/", bugs[i], plates[j],"sigma.svg"), width = 6, height = 4, bg = "white")
    print(p %>% dplyr::group_by(Bug_ID,Plate_no,Compound,Replicate_no) %>%
      ggplot(aes(x=sigma)) + geom_histogram(aes(fill=Replicate_no),color="black", lwd=0.2) +
      scale_fill_nejm() +
      ggtitle(paste(bugs[i],p$Sp_short,plates[j],"sigma values",sep=":")) +
      th + theme_minimal())
    dev.off()
    
    reps<-unique(p$Replicate_no)
    reps
    
    for(a in 1:length(reps)){
      
      r<- p %>% filter(Replicate_no == reps[a])
      head(r)
      
      rownames(r)<-r$col_name
      
      pca.res <- prcomp(r %>% dplyr::select(k:sigma,t_mid:auc_e), center=TRUE, scale=TRUE)
      head(pca.res)
      summary(pca.res)
      
      U<-pca.res$x
      head(U)
      
      outliers_Ustd<-apply(U, 2, function(x) which( (abs(x - median(x)) / mad(x)) > 6 ))
      outliers_Ustd
      outs<-names(outliers_Ustd$PC1)
      outs
      
      if(length(outs) > 0){
        
        for(z in 1:length(outs)){
          
          # str<-str_split(outs[z], "-")
          # str
          
          outliers<-data.frame(Bug_ID = r$Bug_ID[1], Plate_no = r$Plate_no[1], 
                              col_name = outs[z],
                              Replicate_no = r$Replicate_no[1])
          head(outliers)

          outlier_wells_pca<-rbind(outlier_wells_pca,outliers)
          head(outlier_wells_pca)
        }
      }
    
      df_pca<-as_tibble(list(PC1=pca.res$x[,1],
                             PC2=pca.res$x[,2],
                             samples = r$col_name))
      df_pca
    
    
      CairoSVG(file=paste("/Users/periwal/ShikiFactory/WP3/SF100-WP3/Figures/", bugs[i], plates[j],reps[a],"pca.svg"), width = 6, height = 4, bg = "white")
      print(df_pca %>% ggplot(aes(x=PC1,y=PC2, label=samples)) +
              geom_text(size = 3, aes(color=r$Replicate_no)) + scale_color_jama())
      dev.off()
      
     
    }
    
  }
}

nrow(outlier_wells_pca) #553
head(outlier_wells_pca)
head(params_bugs_cmp)

# remove outliers (noisy) wells

head(params_bugs_cmp)
nrow(params_bugs_cmp)
head(outlier_wells_pca)

df1<-inner_join(params_bugs_cmp,outlier_wells_pca)
nrow(df1) #553
head(df1)

params_bugs_cmp_filtered<-anti_join(params_bugs_cmp,df1)
nrow(params_bugs_cmp_filtered) #28,231
head(params_bugs_cmp_filtered)

# filter questionable fits

params_bugs_cmp_filtered<-params_bugs_cmp_filtered %>% filter(note == "")
nrow(params_bugs_cmp_filtered) #28,166

params_bugs_cmp_filtered<-params_bugs_cmp_filtered[,c(1:4,6,20,23:25)]
head(params_bugs_cmp_filtered)
nrow(params_bugs_cmp_filtered) #28,166

write.table(params_bugs_cmp_filtered, file = "/Users/periwal/ShikiFactory/WP3/SF100-WP3/Figures/pca_filtered", sep = "\t", quote = FALSE, row.names = FALSE)

############ 1. Normalization

## control-based normalization

#mean AUC of control wells
control_wells<-params_bugs_cmp_filtered %>% filter(Compound == 'DMSO') %>% dplyr::group_by(Bug_ID,Replicate_no,Plate_no) %>%
  summarise(median_auc_control_wells = median(auc_l))
View(control_wells)
nrow(control_wells) #298 out of 300

##### normalizing each well by plate controls

# merge controls auc
params_bugs_cmp_auc<-merge(params_bugs_cmp_filtered,control_wells,by=c("Bug_ID","Replicate_no","Plate_no"))
nrow(params_bugs_cmp_auc) #27,999
head(params_bugs_cmp_auc)

# merge technical replicates
#techs_auc<-params_bugs_cmp_auc %>% dplyr::group_by(Bug_ID,Plate_no,Replicate_no,Compound,Species,Sp_short,Phyla,median_auc_control_wells) %>%
#  summarise(avg_auc = median(auc_l))

#View(techs_auc)
#nrow(techs_auc) #9046

#length(unique(techs_auc$Compound))

# normalize each well
normAUCl_control<-params_bugs_cmp_auc %>% dplyr::group_by(Bug_ID,Replicate_no,Plate_no,Compound,Phyla,Species,Sp_short,col_name) %>% 
  mutate(cnormAUC = round((auc_l/median_auc_control_wells), digits = 3))
head(normAUCl_control)
nrow(normAUCl_control) #27,999

write.table(normAUCl_control, file = "/Users/periwal/ShikiFactory/WP3/SF100-WP3/GrowthProfiler/cnormAUCs", sep = "\t", quote = FALSE, row.names = FALSE)

## read normAUCl
normAUCl_control<-data.frame(read.table("/Users/periwal/ShikiFactory/WP3/SF100-WP3/GrowthProfiler/cnormAUCs", header = TRUE, sep = '\t', stringsAsFactors = FALSE))

## all runs normAUC_control
#CairoSVG(file=paste("/Users/periwal/ShikiFactory/WP3/SF100-WP3/Figures/cnormAUC.svg", sep = ""), width = 10, height = 4, bg = "white")
#normAUCl_control %>% dplyr::group_by(Bug_ID,Replicate_no) %>%
#  ggplot(aes(x=Sp_short,y=cnormAUC)) + geom_boxplot(aes(fill = Plate_no), outlier.size = 0.01, lwd=0.15) +
#  theme_minimal() + th + scale_fill_npg() 
#dev.off()

CairoSVG(file=paste("/Users/periwal/ShikiFactory/WP3/SF100-WP3/Figures/cnormAUC_galch.svg", sep = ""), width = 8, height = 5, bg = "white")
normAUCl_control %>% filter(Plate_no == 'Plate1-GALCH') %>% dplyr::group_by(Bug_ID,Plate_no,Replicate_no) %>%
  ggplot(aes(cnormAUC)) + geom_density(aes(fill=Replicate_no,color=Replicate_no), alpha=0.3) +
  facet_wrap("Sp_short", scales = "free") +
  scale_fill_npg() + scale_color_npg() + th + theme_bw()
dev.off()

###################### boxplot of effect of comb drugs alone (compared with bug alone)
View(normAUCl_control)

normAUCl_control_asp<-normAUCl_control

for(i in 1:nrow(normAUCl_control_asp)){
  
  line = normAUCl_control[i,]
  #print(line)
  
  if(line$Compound == 'Aspartame' & line$col_name %in% c('H7','H8','H9','H10','H11','H12')){
    
    #line$Compound<-'DMSO'
    normAUCl_control_asp[i,]$Compound<-'Asp'
    
  }
}

View(normAUCl_control_asp)

SF_drugs_bugs<-normAUCl_control_asp %>% dplyr::group_by(Bug_ID,Plate_no,Replicate_no) %>% filter(Compound %in% c("Vanillin","Caffeine","Asp","Duloxetine") & !Plate_no %in% c('Plate1-GALCH','Plate2-Comm') & !col_name %in% c('D11','D12')) %>%
  dplyr::group_by(Bug_ID,Plate_no,Replicate_no,Species,Compound,Sp_short) %>% summarise(av_aucD = median(cnormAUC))
View(SF_drugs_bugs)
nrow(SF_drugs_bugs) #150

CairoSVG(file=paste("/Users/periwal/ShikiFactory/WP3/SF100-WP3/Figures/comb.svg", sep = ""), width = 13, height = 5, bg = "white")
SF_drugs_bugs %>% dplyr::group_by(Bug_ID) %>% 
  ggplot(aes(x=Sp_short,y=av_aucD)) + geom_boxplot(aes(fill = Plate_no),lwd=0.15) + geom_point(position=position_dodge2(width = 0.8),aes(fill = Plate_no), shape=21) +
  facet_grid(c("Compound","Bug_ID"), scales = "free") +  th + scale_fill_npg() + scale_color_npg() + scale_x_discrete(name = "Species") +
  scale_y_continuous(name = "normalized AUC")
dev.off()

##### plate-based normalization (robust z-score)
# params_bugs_cmp_filtered<-read.table("/Users/periwal/ShikiFactory/WP3/SF100-WP3/Figures/pca_filtered", header = TRUE, sep = '\t')
# head(params_bugs_cmp_filtered)
#  
# # #remove control wells from median calculations
# normAUCl_z<-params_bugs_cmp_filtered %>% dplyr::group_by(Bug_ID,Plate_no,Replicate_no) %>% filter(!Compound == 'DMSO') %>%
#    summarise(pl_med = median(auc_l), mean_ab_dev = mad(auc_l))
# head(normAUCl_z)
#  
# normAUCl_plate<-merge(params_bugs_cmp_filtered,normAUCl_z,by=c('Bug_ID','Plate_no','Replicate_no'))
# head(normAUCl_plate)
# 
# normAUCl_plate<-normAUCl_plate %>% dplyr::group_by(Bug_ID,Plate_no,Replicate_no,col_name) %>% filter(!Compound == 'DMSO') %>%
#   mutate(zscoreAUC = (auc_l-pl_med)/mean_ab_dev)
# head(normAUCl_plate)
# 
# CairoSVG(file=paste("/Users/periwal/ShikiFactory/WP3/SF100-WP3/Figures/zscoreAUC_Plate1.svg", sep = ""), width = 8, height = 5, bg = "white")
# normAUCl_plate %>% filter(Plate_no == 'Plate1-GALCH') %>% dplyr::group_by(Bug_ID,Plate_no,Replicate_no) %>%
#   ggplot(aes(zscoreAUC)) + geom_density(aes(fill=Replicate_no,color=Replicate_no), alpha=0.3) +
#   facet_wrap("Sp_short", scales = "free") +
#   scale_fill_npg() + scale_color_npg() + th + theme_bw()
# dev.off()
# 
#CairoSVG(file=paste("/Users/periwal/ShikiFactory/WP3/SF100-WP3/Figures/WP3_zscoreAUC.svg", sep = ""), width = 11, height = 3, bg = "white")
#normAUCl_plate %>% filter(Plate_no %in% c('Plate1-GALCH','Plate2-Comm')) %>% dplyr::group_by(Bug_ID,Plate_no,Replicate_no) %>%
#  ggplot(aes(zscoreAUC)) + geom_density(aes(fill=Replicate_no,color=Replicate_no), alpha=0.3) +
#  facet_grid("Plate_no ~ Sp_short", scales = "free") +
#  scale_fill_npg() + scale_color_npg() + th + theme_bw()
#dev.off()

# ###################### boxplot of effect of comb drugs alone (compared with bug alone)
# 
# SF_drugs_bugs<-normAUCl_plate %>% dplyr::group_by(Bug_ID,Plate_no,Replicate_no) %>% filter(Compound %in% c("Vanillin","Caffeine","Aspartame","Duloxetine") & !Plate_no %in% c('Plate1-GALCH','Plate2-Comm')) %>%
#   dplyr::group_by(Bug_ID,Replicate_no,Species,Compound,Sp_short) %>% summarise(av_aucD = median(zscoreAUC))
# View(SF_drugs_bugs)
# nrow(SF_drugs_bugs) #120
# 
# CairoSVG(file=paste("/Users/periwal/ShikiFactory/WP3/SF100-WP3/Figures/zscomb.svg", sep = ""), width = 10, height = 4, bg = "white")
# SF_drugs_bugs %>% dplyr::group_by(Bug_ID) %>% 
#   ggplot(aes(x=Bug_ID,y=av_aucD)) + geom_boxplot(aes(fill = Compound),lwd=0.15) + geom_point(position=position_dodge2(width = 0.8),aes(fill = Compound), shape=21) +
#   facet_wrap("Sp_short", scales = "free",ncol = 5) + theme_minimal() + th + scale_fill_npg() + scale_color_npg() + 
#   scale_y_continuous(name = "normalized AUC") + coord_flip() + theme(axis.text.y = element_blank(), axis.title.y = element_blank())
# dev.off()



######## 2. Hit selection (unpaired (unequal no of observations) t-test / welch two sample t-test)

View(normAUCl_plate)

#head(normAUCl_plate)
#var<-normAUCl_control %>% dplyr::group_by(Bug_ID,Plate_no,Replicate_no) %>% 
#  do(unpairedT_pval(.))

#View(var)

hits<-normAUCl_control %>% dplyr::group_by(Bug_ID,Plate_no,Replicate_no) %>% 
  do(unpairedT_pval(.))
View(hits)
nrow(hits) #12,171


CairoSVG(file=paste("/Users/periwal/ShikiFactory/WP3/SF100-WP3/Figures/pval_dist.svg", sep = ""), width = 3, height = 2, bg = "white")
hits %>% ggplot(aes(x=pv)) + geom_histogram(color="white",fill="#4d4d4d",bins = 30,lwd=0.1) + th
dev.off()

#### combining p values across replicates using fisher's method

combined_pval<-hits %>% dplyr::group_by(Bug_ID,Plate_no,Compound,Phyla,Species,Sp_short) %>%
  summarise(combined_pv = sumlog(pv)[["p"]], count = length(pv))
View(combined_pval)
nrow(combined_pval) #4247

sig_hits<-combined_pval %>% filter(combined_pv < 0.01 & count > 2 & Plate_no %in% c('Plate1-GALCH','Plate2-Comm'))
head(sig_hits)
nrow(sig_hits) #280

CairoSVG(file=paste("/Users/periwal/ShikiFactory/WP3/SF100-WP3/Figures/pval_combined.svg", sep = ""), width = 5, height = 2.3, bg = "white")
sig_hits %>% ggplot(aes(combined_pv)) + geom_histogram(aes(fill=Plate_no),color="white",lwd=0.1) + 
  th + scale_fill_lancet()
dev.off()

### multiple hypotheses testing: error correction

p_bh<-combined_pval %>% dplyr::group_by(Bug_ID,Plate_no) %>%
  mutate(p_bh = p.adjust(combined_pv,method = "BH"))
View(p_bh)
nrow(p_bh)

hits_bh<-p_bh %>% filter(p_bh < 0.01 & Plate_no %in% c('Plate1-GALCH','Plate2-Comm'))
nrow(hits_bh) #232
head(hits_bh)

CairoSVG(file=paste("/Users/periwal/ShikiFactory/WP3/SF100-WP3/Figures/pval_bh.svg", sep = ""), width = 5, height = 2, bg = "white")
hits_bh %>% ggplot(aes(p_bh)) + geom_histogram(aes(fill=Plate_no)) + 
  th + scale_fill_lancet()
dev.off()

#galch
CairoSVG(file=paste("/Users/periwal/ShikiFactory/WP3/SF100-WP3/Figures/p_bh_GALCH_heatmap.svg", sep = ""), width = 3.5, height = 7, bg = "white")
hits_bh %>% filter(Plate_no == "Plate1-GALCH") %>% ggplot(aes(x=Sp_short,y=Compound)) + geom_tile(aes(fill=as.character(count)), color = "white", lwd = 0.1) + theme_bw() +
  th + theme(axis.text.x = element_text(angle = 90,hjust = 1), axis.ticks = element_blank(), panel.grid = element_blank()) +
  scale_fill_nejm(name = "#Reps") + scale_y_discrete(name = "Compound") 
dev.off()

#commercial
CairoSVG(file=paste("/Users/periwal/ShikiFactory/WP3/SF100-WP3/Figures/p_bh_Comm_heatmap.svg", sep = ""), width = 5, height = 7, bg = "white")
hits_bh %>% filter(Plate_no == "Plate2-Comm") %>% ggplot(aes(x=Sp_short,y=Compound)) + geom_tile(aes(fill=as.character(count)), color = "white", lwd = 0.1) + theme_bw() +
  th + theme(axis.text.x = element_text(angle = 90,hjust = 1), axis.ticks = element_blank(), panel.grid = element_blank()) + 
  scale_fill_nejm(name = "#Reps") + scale_y_discrete(name = "Compound") 
dev.off()
  
#combination compound effect
CairoSVG(file=paste("/Users/periwal/ShikiFactory/WP3/SF100-WP3/Figures/comb_effect.svg", sep = ""), width = 5, height = 3, bg = "white")
p_bh %>% filter(p_bh < 0.05 & Compound %in% c("Vanillin","Caffeine","Asp","Duloxetine")) %>% ggplot(aes(x=Sp_short,y=Compound)) + 
  geom_point(aes(fill=Plate_no), shape=21, position = position_dodge2(width = 0.5)) + theme_bw() +
  th + theme(axis.text.x = element_text(angle = 90,hjust = 1), axis.ticks = element_blank()) + 
  scale_fill_aaas() + scale_y_discrete(name = "Compound") 
dev.off()


# combinations
CairoSVG(file=paste("/Users/periwal/ShikiFactory/WP3/SF100-WP3/Figures/comb_effect_galch.svg", sep = ""), width = 8.5, height = 7, bg = "white")
p_bh %>% filter(p_bh < 0.05 & Plate_no %in% c('Plate1-GALCH','Plate1+A','Plate1+C','Plate1+D','Plate1+V')) %>% ggplot(aes(x=Sp_short,y=Compound)) + 
  geom_point(aes(fill=as.character(count)), shape=21) + theme_bw() +
  th + theme(axis.text.x = element_text(angle = 90,hjust = 1), axis.ticks = element_blank(), axis.title.x = element_blank()) + 
  scale_fill_aaas(name = "#Reps") + scale_y_discrete(name = "Compound") +
  facet_grid(~Plate_no)
dev.off()


CairoSVG(file=paste("/Users/periwal/ShikiFactory/WP3/SF100-WP3/Figures/advan_effect.svg", sep = ""), width = 10, height = 4, bg = "white")
p_bh %>% filter(p_bh < 0.05 & Plate_no %in% c('Plate2-Comm','Plate2+A','Plate2+C','Plate2+D','Plate2+V') & !Compound %in% c("Vanillin","Caffeine","Asp","Duloxetine") & Compound %in% c('Advantame','Aspartame')) %>% ggplot(aes(x=Sp_short,y=Compound)) + 
  geom_point(aes(fill=as.character(count)), shape=21) + theme_bw() +
  th + theme(axis.text.x = element_text(angle = 90,hjust = 1), axis.ticks = element_blank(), axis.title.x = element_blank()) + 
  scale_fill_aaas(name = "#Reps") + scale_y_discrete(name = "Compound") +
  facet_grid(~Plate_no)
dev.off()




# pooled effect sizes (all replicates)

pooled_hits<-hits %>% dplyr::group_by(Bug_ID,Plate_no,Compound,Phyla,Species,Sp_short) %>% 
  do(pool_es(.))
head(pooled_hits)
View(pooled_hits)


pooled_hits %>% filter(Plate_no == "Plate1-GALCH") %>% dplyr::group_by(Bug_ID,Sp_short,Compound) %>%
  ggplot(aes(x=Sp_short,y=Compound)) + geom_tile(aes(fill=pooled_effect),color="#000000") +
  scale_fill_distiller(palette = "RdGy") + th + theme(axis.text.x = element_text(angle = 90,hjust = 1), axis.ticks = element_blank(), panel.grid = element_blank(), axis.title = element_blank())

sig_hits<-pooled_hits %>% filter(pval_es < 0.05)
View(sig_hits)
nrow(sig_hits) # 279

pooled_hits_monica_GALCH<-pooled_hits %>% filter(Plate_no == "Plate1-GALCH")
pooled_hits_monica_GALCH<-pooled_hits_monica_GALCH[,c(1:10)]
View(pooled_hits_monica_GALCH)
nrow(pooled_hits_monica_GALCH) #420 (42 compounds x 10 bugs)

write.csv(pooled_hits_monica_GALCH, file = "galch_es_10_bugs.csv", row.names = FALSE)


CairoSVG(file=paste("/Users/periwal/ShikiFactory/WP3/SF100-WP3/Figures/sig_GALCH.svg", sep = ""), width = 10, height = 6, bg = "white")
sig_hits %>% filter(Plate_no == "Plate1-GALCH") %>% dplyr::group_by(Bug_ID,Sp_short,Compound) %>%
  ggplot(aes(x=pooled_effect,y=Compound)) + geom_errorbar(aes(xmin = low.ci, xmax = high.ci, color=Bug_ID)) + geom_point(aes(color=Bug_ID)) +
  scale_color_npg() + th + theme(axis.text.x = element_text(angle = 90,hjust = 1), axis.ticks = element_blank(), panel.grid = element_blank(), axis.title = element_blank()) +
  facet_grid(~Sp_short, scales = "free")
dev.off()


CairoSVG(file=paste("/Users/periwal/ShikiFactory/WP3/SF100-WP3/Figures/sig_Comm.svg", sep = ""), width = 10, height = 6, bg = "white")
sig_hits %>% filter(Plate_no == "Plate2-Comm") %>% dplyr::group_by(Bug_ID,Sp_short,Compound) %>%
  ggplot(aes(x=pooled_effect,y=Compound)) + geom_errorbar(aes(xmin = low.ci, xmax = high.ci, color=Bug_ID)) + geom_point(aes(color=Bug_ID)) +
  scale_color_npg() + th + theme(axis.text.x = element_text(angle = 90,hjust = 1), axis.ticks = element_blank(), panel.grid = element_blank(), axis.title = element_blank()) +
  facet_grid(~Sp_short, scales = "free")
dev.off()


sig_hits %>% filter(Plate_no == "Plate1-GALCH") %>% dplyr::group_by(Bug_ID,Sp_short,Compound) %>%
  ggplot(aes(x=pooled_effect,y=Compound)) + geom_errorbar(aes(xmin = low.ci, xmax = high.ci, color=Bug_ID)) + geom_point(aes(color=Bug_ID)) +
  scale_color_npg() + th + theme(axis.text.x = element_text(angle = 90,hjust = 1), axis.ticks = element_blank(), panel.grid = element_blank(), axis.title = element_blank()) +
  facet_grid(~Sp_short, scales = "free")


sig_hits %>% filter(Plate_no == "Plate2+D") %>% dplyr::group_by(Bug_ID,Sp_short,Compound) %>%
  ggplot(aes(x=pooled_effect,y=Compound)) + geom_errorbar(aes(xmin = low.ci, xmax = high.ci, color=Bug_ID)) + geom_point(aes(color=Bug_ID)) +
  scale_color_npg() + th + theme(axis.text.x = element_text(angle = 90,hjust = 1), axis.ticks = element_blank(), panel.grid = element_blank(), axis.title = element_blank()) +
  facet_grid(~Sp_short, scales = "free")


#bug_aucs<-hits %>% filter(Bug_ID == 'NT24007', Plate_no == 'Plate1-GALCH',Compound == 'GAL-1-004')





###################### plot of compounds in control plates
#Galchimia compounds
CairoSVG(file=paste("/Users/periwal/ShikiFactory/WP3/SF100-WP3/Figures/plate1_heatmap.svg", sep = ""), width = 2.8, height = 7, bg = "white")
sig_hits %>% filter(Plate_no == "Plate1-GALCH") %>% dplyr::group_by(Bug_ID,Sp_short,Compound) %>%
  ggplot(aes(x=Sp_short,y=Compound)) + geom_tile(aes(fill=combined_pv),color="#000000") +
  scale_fill_distiller(palette = "RdGy") + th + theme(axis.text.x = element_text(angle = 90,hjust = 1), axis.ticks = element_blank(), panel.grid = element_blank(), axis.title = element_blank())
dev.off()


#Commercial sweeteners
CairoSVG(file=paste("/Users/periwal/ShikiFactory/WP3/SF100-WP3/Figures/plate2_heatmap.svg", sep = ""), width = 4, height = 7, bg = "white")
sig_hits %>% filter(Plate_no == "Plate2-Comm") %>% dplyr::group_by(Bug_ID,Sp_short,Compound) %>%
  ggplot(aes(x=Sp_short,y=Compound)) + geom_tile(aes(fill=combined_pv),color="#000000") +
  scale_fill_distiller(palette = "RdGy") + th + theme(axis.text.x = element_text(angle = 90,hjust = 1), axis.ticks = element_blank(), panel.grid = element_blank(), axis.title = element_blank())
dev.off()

#Commercial sweeteners + Sweeteners
CairoSVG(file=paste("/Users/periwal/ShikiFactory/WP3/SF100-WP3/Figures/plate1+A_heatmap.svg", sep = ""), width = 2.8, height = 7, bg = "white")
normAUCl %>% filter(Plate_no == "Plate2+A") %>% dplyr::group_by(Bug_ID,Replicate_no,Sp_short,Compound) %>% filter(!col_name %in% c('H6','H7','H8','H9','H10','H11','H12')) %>% summarise(av_auc = mean(normAUC)) %>%
  ggplot(aes(x=Sp_short,y=Compound)) + geom_tile(aes(fill=av_auc),color="#000000") +
  scale_fill_distiller(palette = "PiYG") + th + theme(axis.text.x = element_text(angle = 90,hjust = 1), axis.ticks = element_blank(), panel.grid = element_blank(), axis.title = element_blank())
dev.off()







###### fit distribution per bug plate wise (all replicates pooled)

# Plot normAUCs distribution
CairoSVG(file=paste("/Users/periwal/ShikiFactory/WP3/SF100-WP3/Figures/density_distribution.svg", sep = ""), width = 5, height = 2, bg = "white")
normAUCl %>% ggplot(aes(normAUC)) + geom_density(aes(fill=Species, color = Species), alpha = 0.5) + th + 
  scale_fill_npg() + scale_color_npg() + theme_minimal()
dev.off()

nrow(normAUCl) #11480
head(normAUCl)

# Effect of single compounds only
single_sweet<-normAUCl %>% filter(Plate_no %in% c("Plate1-GALCH","Plate2-Comm"))
head(single_sweet)

#plot normalized AUCs in control wells
CairoSVG(file=paste("/Users/periwal/ShikiFactory/WP3/SF100-WP3/Figures/control_wells_normAUC.svg", sep = ""), width = 5, height = 3, bg = "white")
single_sweet %>% filter(Compound == 'DMSO') %>% dplyr::group_by(Bug_ID,Replicate_no) %>% 
  ggplot(aes(x=Sp_short,y=normAUC)) + geom_boxplot(aes(fill = Replicate_no), lwd=0.1, position = "dodge2", outlier.size = 0.2) +
  th + scale_fill_jama() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

# Define reference wells
ref_wells<-single_sweet %>% filter(Compound == 'DMSO')
head(ref_wells)
ref_wells<-ref_wells[,c(1:6,23:27)]
ref_wells["label"]<-"reference"
nrow(ref_wells) #861
ncol(ref_wells) #12
View(ref_wells)

#descdist(ref_wells$normAUC, discrete = FALSE, boot = 1000)

sample_wells<-single_sweet %>% filter(!(Compound == 'DMSO'))
head(sample_wells)
sample_wells<-sample_wells[,c(1:6,23:27)]
sample_wells["label"]<-"sample"
nrow(sample_wells) #2014
ncol(sample_wells) #12
View(sample_wells)

total<-rbind(data.frame(ref_wells),data.frame(sample_wells))
head(total)
nrow(total) #2299

#compute z scores (first take mean of technical replicates)
z_scores<-sample_wells %>% dplyr::group_by(Bug_ID,Replicate_no,Plate_no,Species,Sp_short,Compound) %>% summarise(mean_normAUC = mean(normAUC)) %>% 
  mutate(plate_mean = mean(mean_normAUC), plate_sd = sd(mean_normAUC)) %>% dplyr::group_by(Bug_ID,Replicate_no,Plate_no,Species,Sp_short,Compound) %>%
  mutate(zscore = (mean_normAUC-plate_mean)/plate_sd)
View(z_scores)

# plot zscores
# use Z scores - to determine validity of hits
CairoSVG(file=paste("/Users/periwal/ShikiFactory/WP3/SF100-WP3/Figures/control_zscores.svg", sep = ""), width = 6, height = 3, bg = "white")
z_scores %>% ggplot(aes(zscore)) + geom_density(aes(fill=Species, color = Species), alpha = 0.5) + th + 
  scale_fill_npg() + scale_color_npg() + theme_minimal()
dev.off()

# filter zscores

significant_hits<-z_scores %>% filter(zscore < -1.5 | zscore > 1.5) %>% dplyr::group_by(Bug_ID,Plate_no,Compound) %>%
  mutate(count = length(Compound))
View(significant_hits)

# NON-PARAMETRIC TEST: two sample Welch test (different sample sizes)

# filter all sample wells
# hits<-total %>%  dplyr::group_by(Bug_ID,Replicate_no,Plate_no) %>% 
#   do(welch_pval(.))
# 
# View(hits)

#PROBLEM: many compounds don't have a paired tech replicate so can't be 
#processed with this method. P-values are strange 


#PARAMETRIC TEST: Requires distribution, can be used for single observation
set.seed(444)

#Compare different distributions (single compound plates)

# dist_compare<-total %>% filter(Plate_no %in% c("Plate1-GALCH","Plate2-Comm")) %>% dplyr::group_by(Bug_ID,Replicate_no,Plate_no) %>% 
#   do(fitCompare(.))
# names(dist_compare)[names(dist_compare) == "dist.i."] <- "distribution"
# names(dist_compare)[names(dist_compare) == "chisqpv"] <- "pv"
# View(dist_compare)
# 
# dist_compare %>% filter(pv > 0.05) %>% ggplot(aes(x=distribution,y=pv)) + geom_boxplot(aes(fill=Plate_no))
# 
# #Plot plate1 and plate 2 graphs
# 

distri<-total %>% filter(Plate_no %in% c("Plate1-GALCH","Plate2-Comm")) 
View(distri)
nrow(distri)

# # cullen n frey graph
CairoSVG(file=paste("/Users/periwal/ShikiFactory/WP3/SF100-WP3/Figures/distribution.svg", sep = ""), width = 5, height = 5, bg = "white")
descdist(distri$normAUC, discrete = FALSE, boot = 1000)
dev.off()

View(total)
nrow(total)

#Replicate and Plate wise

hits_rep_plate<-total %>% dplyr::group_by(Bug_ID,Replicate_no,Plate_no) %>% do(welch_pval(.))
View(hits_rep_plate)
nrow(hits_rep_plate) #11480

CairoSVG(file=paste("/Users/periwal/ShikiFactory/WP3/SF100-WP3/Figures/pval_nullhypothesis.svg", sep = ""), width = 3, height = 2, bg = "white")
hits %>% filter(label=='reference') %>% ggplot(aes(x=pv)) + geom_histogram(color="white",fill="#4d4d4d",bins = 30) + th
dev.off()

CairoSVG(file=paste("/Users/vinitaperiwal/GrowthCurver/Figures/pval_alternatehypothesis.svg", sep = ""), width = 3, height = 2, bg = "white")
wil_test %>% filter(label=='sample') %>% ggplot(aes(x=pv)) + geom_histogram(color="white",fill="#4d4d4d",bins = 30) + th
dev.off()

hits_rep_plate %>% ggplot(aes(x=pv)) + geom_histogram(aes(fill=label),color="white",bins = 30) + th + scale_fill_nejm()
