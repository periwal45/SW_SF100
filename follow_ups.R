setwd("/Users/periwal/ShikiFactory/WP3/SW_SF100/follow_ups")
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
library(car) #for Boxplots with labels
library(corrplot)
library(Hmisc)
library(pracma)
library(dunn.test)
library(purrr)
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
library(EnhancedVolcano)
library(wesanderson)

# source the script with functions
source("/Users/periwal/ShikiFactory/WP3/SW_SF100/functions_followup.R")

#bug names
bug_desc<-read.table(file = "../Bug_desc", sep = "\t", header = TRUE)
head(bug_desc)
nrow(bug_desc) #26

single_plate_layout<-data.frame(read.table("../single_cmp_layout", header = TRUE, sep = '\t'))
head(single_plate_layout)
nrow(single_plate_layout) #3552

CB_plate_layout<-data.frame(read.table("../CB_layout", header = TRUE, sep = '\t'))
head(CB_plate_layout)
nrow(CB_plate_layout) #384


##################################### READING RAW FILES AND FITTING CURVES
# Currently missing:
# NT5076 not growing in any replicate (kick-out)
# NT5023 Rep-3 not growing
# NT5025, 5081 has only 15 points for rep3 (not sufficient) 

#read all files with .tab extension
file.names <- list.files(path = '.', recursive = TRUE, pattern = "\\.tab$") #recursive reads through all subfolders and files
file.names

#loop for each .tab file, annotates
for(i in 1:length(file.names)){
  
  #read file/table
  input_file<-tools::file_path_sans_ext(file.names[i]) #read filename w/o extension
  input_file
  
  if(str_detect(input_file, "plate")){
  #create a new dataframe with additional annotations for each file
  
  annot <- tibble(Bug_ID = create_annot_single(input_file)[1], Replicate_no = create_annot_single(input_file)[2], Plate_no = create_annot_single(input_file)[3])
  annot
  X<-data.frame(read.table(file.names[i], header = TRUE, sep = '\t', stringsAsFactors = FALSE))
  
  head(X)
  colnames(X)[1]<-"time"
  #X$time<-as.numeric(X$time)
  
  annot_file<-merge(annot,X)
  head(annot_file)
  
  output_dir <- file.path('.', annot$Bug_ID)
  
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
  } else {
    print("Dir already exists!")
  }
  
  #write annotated files
  annot_outfile<-paste0(output_dir,"/",input_file,".annot")
  write.table(annot_file, file = annot_outfile, sep = "\t", quote = FALSE, row.names = FALSE)
  
  }
  
  if(str_detect(input_file, "CB")){
    
    annot <- tibble(Bug_ID = create_annot_CB(input_file)[1], Replicate_no = create_annot_CB(input_file)[2], Plate_no = create_annot_CB(input_file)[3])
    annot
    X<-data.frame(read.table(file.names[i], header = TRUE, sep = '\t', stringsAsFactors = FALSE))
    
    head(X)
    colnames(X)[1]<-"time"
    #X$time<-as.numeric(X$time)
    
    annot_file<-merge(annot,X)
    head(annot_file)
    
    output_dir <- file.path('.', annot$Plate_no)
    
    if (!dir.exists(output_dir)){
      dir.create(output_dir)
    } else {
      print("Dir already exists!")
    }
    
    #write annotated files
    annot_outfile<-paste0(output_dir,"/",input_file,".annot")
    write.table(annot_file, file = annot_outfile, sep = "\t", quote = FALSE, row.names = FALSE)
    
  }
  
}


#merge all replicates of a bug
# Merge files across all replicates of a strain (.annot files)
dir.names<-list.dirs(path = '.', recursive = FALSE)
dir.names

for(d in 1:length(dir.names)){
  
  if(str_detect(dir.names[d], "./NT") | str_detect(dir.names[d], "./plate")){
    
    all_rep<-do.call(rbind, lapply(a<-list.files(path=dir.names[d], pattern="\\.annot$", full.names = TRUE), function(i){
      read.table(i,header=TRUE, sep="\t", stringsAsFactors = FALSE)}))
    write.table(all_rep, file = paste0(dir.names[d],"_all_replicates.merged"), sep = "\t", quote = FALSE, row.names = FALSE)
    
  }
  
}


#read all files with .merged extension
file.merged<-list.files(path = '.', pattern = ".merged$") #recursive reads through all subfolders and files
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
  head(Y)
  lapply(Y, class)
  
  Y[, 4:ncol(Y)] <- lapply(4:ncol(Y), function(x) as.numeric(Y[[x]]))
  
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
  
  write.table(Z, file = paste0("trim_points/",in_file,".trim"), sep = "\t", quote = FALSE, row.names = FALSE)
  
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
              current_well
 
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
              annot_mod<-cbind(bug,rep,plate,col_name=col_name,mod)
              annot_mod
              
              fitted_values<-rbind(fitted_values, annot_mod)
              fitted_values
              
              fitted_param<-data.frame(cbind(bug,rep,plate,col_name,gc_fit$vals$k,gc_fit$vals$n0,gc_fit$vals$r,gc_fit$vals$sigma,gc_fit$vals$df,gc_fit$vals$t_mid,gc_fit$vals$auc_l,gc_fit$vals$auc_e,gc_fit$vals$note))
              colnames(fitted_param)<-c("Bug_ID","Replicate_no","Plate_no","well","k","n0","r","sigma","df","t_mid","auc_l","auc_e","note")
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
nrow(fitted_values) #depends on time points
View(fitted_values)

head(gc_fit_params)
nrow(gc_fit_params) #should be equal to number of wells
View(gc_fit_params)

write.table(fitted_values, file = "Follow_GC_fitted_model_values", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gc_fit_params, file = "Follow_GC_fit_params", sep = "\t", quote = FALSE, row.names = FALSE)

############################### Read fitted values
#read fitted values of models
fits<-data.frame(read.table("Follow_GC_fitted_model_values", header = TRUE, sep = '\t', stringsAsFactors = FALSE))
View(fits)
nrow(fits) #1,70,496

single_merge_plate_layout_fits<-merge(single_plate_layout, fits, by=c("Bug_ID","Plate_no","col_name"))
single_merge_plate_layout_fits<-merge(single_merge_plate_layout_fits,bug_desc, by="Bug_ID")
nrow(single_merge_plate_layout_fits) #115200

View(single_merge_plate_layout_fits)

CB_merge_plate_layout_fits<-merge(CB_plate_layout, fits, by=c("Bug_ID","Plate_no","col_name"))
CB_merge_plate_layout_fits<-merge(CB_merge_plate_layout_fits,bug_desc, by="Bug_ID")
nrow(CB_merge_plate_layout_fits) #55296

View(CB_merge_plate_layout_fits)

write.table(single_merge_plate_layout_fits, file = "Figures/single_fitted_values", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(CB_merge_plate_layout_fits, file = "Figures/CB_fitted_values", sep = "\t", quote = FALSE, row.names = FALSE)

############### plot raw and fitted ODs
merge_plate_layout_fits<-read.table(file = "Figures/single_fitted_values", sep = "\t", header = TRUE)
head(merge_plate_layout_fits)

### plot plate wise replicates for each bug

bugs<-unique(merge_plate_layout_fits$Bug_ID)
bugs

for(i in 1:length(bugs)){
  
  f<-merge_plate_layout_fits %>% filter(Bug_ID == bugs[i])
  
  plates<-unique(f$Plate_no)
  plates
  
  for(j in 1:length(plates)){
    
    p<-f %>% filter(Plate_no == plates[j])
    head(p)
    
    CairoSVG(file=paste("Figures/", bugs[i], "_", plates[j], ".svg", sep = ""), width = 8, height = 5, bg = "white")
    print(p %>% dplyr::group_by(Bug_ID,Plate_no,compound,Replicate_no,dose,col_name) %>%
            ggplot(aes(x=time, y=fitted, group=dose)) + geom_line(aes(color=as.character(dose)), size=0.2) + scale_color_d3(palette = "category20") +
            geom_point(aes(y=OD,color=as.character(dose)), size=0.000001) + facet_wrap(c("compound","Replicate_no"),scales = "free") +  
            ggtitle(paste(bugs[i],p$Sp_short,plates[j],sep=":")) + scale_x_continuous(name = "Time (Hours)", limits = c(0,25)) + scale_y_continuous(name="bg corrected OD") + 
            th + theme_minimal())
    dev.off()
    
  }
}

#CB
merge_plate_layout_fits<-read.table(file = "Figures/CB_fitted_values", sep = "\t", header = TRUE)
head(merge_plate_layout_fits)

merge_plate_layout_fits<-merge_plate_layout_fits %>% mutate(dose = paste0(dose1,'_',dose2), compound = paste0(cmp1,'_',cmp2))
  
### plot plate wise replicates for each bug

bugs<-unique(merge_plate_layout_fits$Bug_ID)
bugs

for(i in 1:length(bugs)){
  
  f<-merge_plate_layout_fits %>% filter(Bug_ID == bugs[i])
  
  plates<-unique(f$Plate_no)
  plates
  
  for(j in 1:length(plates)){
    
    p<-f %>% filter(Plate_no == plates[j])
    head(p)
    
    CairoSVG(file=paste("Figures/", bugs[i], "_", plates[j], ".svg", sep = ""), width = 15, height = 10, bg = "white")
    print(p %>% dplyr::group_by(Bug_ID,Plate_no,compound,Replicate_no) %>%
            ggplot(aes(x=time, y=fitted)) + geom_line(aes(color=Replicate_no), size=0.2) + scale_color_d3(palette = "category20") +
            geom_point(aes(y=OD,color=Replicate_no), size=0.000001) + facet_wrap(c("compound","dose"),scales = "free") +  
            ggtitle(paste(bugs[i],p$Sp_short,plates[j],sep=":")) + scale_x_continuous(name = "Time (Hours)", limits = c(0,25)) + scale_y_continuous(name="bg corrected OD") + 
            th + theme_minimal())
    dev.off()
    
  }
}

