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
library(EnhancedVolcano)
library(wesanderson)

# source the script with functions
source("/Users/periwal/ShikiFactory/WP3/SW_SF100/Functions.R")

#bug names
bug_desc<-read.table(file = "../Bug_desc", sep = "\t", header = TRUE)
head(bug_desc)
nrow(bug_desc) #26

plate_layout<-data.frame(read.table("../plate_layout", header = TRUE, sep = '\t'))
head(plate_layout)
nrow(plate_layout) #14,976

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
  
  #create a new dataframe with additional annotations for each file
  
  annot <- tibble(Bug_ID = create_annot(input_file)[1], Replicate_no = create_annot(input_file)[2], Plate_no = create_annot(input_file)[3])
  annot
  X<-data.frame(read.table(file.names[i], header = TRUE, sep = '\t', stringsAsFactors = FALSE))
  
  head(X)
  colnames(X)[1]<-"time"
  #X$time<-as.numeric(X$time)
  
  annot_file<-merge(annot,X)
  head(annot_file)
  
  #write annotated files
  annot_outfile<-paste0(input_file,".annot")
  write.table(annot_file, file = annot_outfile, sep = "\t", quote = FALSE, row.names = FALSE)
  
}


#merge all replicates of a bug
# Merge files across all replicates of a strain (.annot files)
dir.names<-list.dirs(path = '.', recursive = FALSE)
dir.names

for(d in 1:length(dir.names)){
  
  if(str_detect(dir.names[d], "./NT")){
    
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
              
              # # find zero OD and start reading from there
              # ind<-which(current_well[,col_name] == 0)
              # ind
              # length(ind)
              # 
              # if(ind %in% c('1','2','3','4')){
              #   
              #   m<-match(ind, c('1','2','3','4'))
              #   m<-max(na.omit(m))
              #   m
              #   current_well<-subset(current_well[m:nrow(current_well),])
              #   current_well$time<-c(0:(nrow(current_well)-1))
              #   current_well
              #   
              # }
              
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
nrow(fitted_values) #depends on time points
View(fitted_values)

fitted_values<-fitted_values %>% filter(Bug_ID != 'NT5076')
gc_fit_params<-gc_fit_params %>% filter(Bug_ID != 'NT5076')

head(gc_fit_params)
nrow(gc_fit_params) #should be equal to number of wells: 96 for 1 plate 41,472
View(gc_fit_params)

write.table(fitted_values, file = "GC_fitted_model_values", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gc_fit_params, file = "GC_fit_params", sep = "\t", quote = FALSE, row.names = FALSE)


############################### Read fitted values
#read fitted values of models
fits<-data.frame(read.table("GC_fitted_model_values", header = TRUE, sep = '\t', stringsAsFactors = FALSE))
head(fits)
nrow(fits) #97,5552

merge_plate_layout_fits<-merge(plate_layout, fits, by=c("Bug_ID","Plate_no","col_name"))
nrow(merge_plate_layout_fits) #97,5552

View(merge_plate_layout_fits)

merge_plate_layout_fits<-merge(merge_plate_layout_fits,bug_desc, by="Bug_ID")
write.table(merge_plate_layout_fits, file = "../Figures/fitted_values", sep = "\t", quote = FALSE, row.names = FALSE)

############### plot raw and fitted ODs
merge_plate_layout_fits<-read.table(file = "../Figures/fitted_values", sep = "\t", header = TRUE)
merge_plate_layout_fits$color<-paste0(merge_plate_layout_fits$Replicate_no,"_",merge_plate_layout_fits$tech)
View(merge_plate_layout_fits)

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
    
    CairoSVG(file=paste("../Figures/", bugs[i], "_", plates[j], ".svg", sep = ""), width = 15, height = 8, bg = "white")
    print(p %>% dplyr::group_by(Bug_ID,Plate_no,compound,Replicate_no) %>%
            ggplot(aes(x=time, y=fitted,group=color)) + geom_line(aes(color=Replicate_no), size=0.2) + scale_color_d3(palette = "category20") +
            geom_point(aes(y=OD,color=tech), size=0.000001) + facet_wrap("compound",scales = "free",ncol = 8, nrow = 6) +  
            ggtitle(paste(bugs[i],p$Sp_short,plates[j],sep=":")) + scale_x_continuous(name = "Time (Hours)", limits = c(0,25)) + scale_y_continuous(name="bg corrected OD") + 
            th + theme_minimal())
    dev.off()
    
  }
}

head(merge_plate_layout_fits)
# example single well plot
CairoSVG(file=paste("../Figures/NT5037_dmso_sacc_treh_p1.svg"), width = 6, height = 4, bg = "white")
merge_plate_layout_fits %>% dplyr::group_by(Bug_ID,Plate_no,compound,Replicate_no) %>% 
  filter(Bug_ID == 'NT5019' & compound %in% c('DMSO','Mannitol') & Plate_no %in% c('plate1')) %>%
  ggplot(aes(x=time, y=fitted,group=color)) + geom_line(aes(color=Replicate_no), size=0.5) + scale_color_d3(palette = "category20") +
  geom_point(aes(y=OD,color=tech), size=1) + facet_wrap("compound") +  
  scale_x_continuous(name = "Time (Hours)", limits = c(0,25)) + scale_y_continuous(name="bg corrected OD") + 
  th + theme_minimal()
dev.off()
#'H1','H2','H3','H4','H5','H6','H7','H8','H9','H10','H11','H12'


############ QC: REMOVE NOISY WELLS
# manually correct noise and try Re-fitting
# 1.NT5021 (plate1-all replicates) - kept, didn't remove
# 2.NT5022 (all)
# 3.NT5023

#read fitted values of models
fit_params<-data.frame(read.table("GC_fit_params", header = TRUE, sep = '\t', stringsAsFactors = FALSE))
head(fit_params)
nrow(fit_params) #44,928 removed one replicate from 3 bugs, that's why numbers differ now

fit_params_layout<-merge(plate_layout, fit_params, by=c("Bug_ID","Plate_no","col_name"))
head(fit_params_layout)

noisy_wells_handpicked<-data.frame(read.table("../noisy_wells", header = TRUE, sep = '\t', stringsAsFactors = FALSE))
head(noisy_wells_handpicked)
dim(noisy_wells_handpicked) #533

df1<-inner_join(fit_params_layout,noisy_wells_handpicked)
View(df1)
nrow(df1) #533

params_bugs_cmp_filtered<-anti_join(fit_params_layout,df1)
nrow(params_bugs_cmp_filtered) #40,939
head(params_bugs_cmp_filtered)

params_bugs_cmp_filtered<-merge(params_bugs_cmp_filtered,bug_desc,by='Bug_ID')
View(params_bugs_cmp_filtered)
nrow(params_bugs_cmp_filtered) #40,939

# plot of auc_l vs auc_e
params_bugs_cmp_filtered %>% group_by(Bug_ID,Sp_short,Replicate_no,Plate_no) %>% 
  ggplot(aes(x=auc_l,y=auc_e)) + geom_point(aes(color=Plate_no), size=1) + geom_abline(lwd=0.1) + stat_cor() + theme_bw() +
  scale_y_continuous(name = "Empirical AUC") + 
  scale_x_continuous(name = "Logistic AUC") + scale_color_npg() + 
  facet_wrap(~Sp_short, scales = "free") + th

############ 1. Normalization

##### plate-based normalization (robust z-score)

View(params_bugs_cmp_filtered)

params<-params_bugs_cmp_filtered[,c(1:6,23,25,20)]
dim(params) #40,939 X 9

cnormAUCl_z<-params %>% filter(compound == 'DMSO') %>% dplyr::group_by(Bug_ID,Plate_no,Replicate_no,Phyla,Sp_short) %>%
  summarise(ctrl_mean = mean(auc_l), std_ctrl = std(auc_l))
head(cnormAUCl_z)
dim(cnormAUCl_z)

normAUCl_z<-merge(params,cnormAUCl_z,by=c("Bug_ID","Replicate_no","Plate_no","Phyla","Sp_short"))
dim(normAUCl_z)
head(normAUCl_z)

znormAUCs<-normAUCl_z %>% dplyr::group_by(Bug_ID,Plate_no,Replicate_no,col_name) %>%
  mutate(czscoreAUC = (auc_l-ctrl_mean)/std_ctrl, FC = auc_l/ctrl_mean) #FC treatment vs control

View(znormAUCs %>% filter(Plate_no == 'plate1'))

write.table(znormAUCs, file = "../Figures/znormAUCs", sep = '\t', row.names = FALSE)

# plotting density of normalized wells
CairoSVG(file=paste("../Figures/zscoreAUC.svg", sep = ""), width = 8.5, height = 5.5, bg = "white")
znormAUCs %>%  dplyr::group_by(Bug_ID,Plate_no,Replicate_no) %>% filter(Plate_no != 'plate1') %>%
  ggplot(aes(czscoreAUC)) + geom_density(aes(fill=Replicate_no,color=Replicate_no), alpha=0.3, lwd=0.3) +
  facet_wrap("Sp_short", scales = "free") +
  scale_fill_jama() + scale_color_jama() + th + theme_bw()
dev.off()

# # plot of control and plate normalized AUCs
# znormAUCs %>% group_by(Bug_ID,Sp_short,Replicate_no,Plate_no) %>% 
#   ggplot(aes(x=czscoreAUC,y=pzscoreAUC)) + geom_point(aes(color=Plate_no), size=1) + geom_abline(lwd=0.1) + stat_cor() + 
#   scale_y_continuous(name = "ctrl norm AUC") + 
#   scale_x_continuous(name = "pl norm AUC") + scale_color_npg() + 
#   facet_grid(c('Plate_no','Sp_short'), scales = "free") + th + theme_bw()
# 

######## 2. Hit selection (unpaired (unequal no of observations) t-test / welch two sample t-test)
znormAUCs<-read.table(file = "../Figures/znormAUCs", header = TRUE, stringsAsFactors = FALSE)
dim(znormAUCs)
head(znormAUCs)

# t-test bug wise and replicate wise
hits_BR<-znormAUCs %>% filter(Plate_no == 'plate1') %>% dplyr::group_by(Bug_ID, Plate_no, Replicate_no) %>% 
  do(unpairedT_pval(.))
head(hits_BR)
nrow(hits_BR) #17,182

hits_BR_noDMSO<-hits_BR %>% filter(compound != 'DMSO')

CairoSVG(file=paste("../Figures/pval_dist.svg", sep = ""), width = 3, height = 2, bg = "white")
hits_BR_noDMSO %>% ggplot(aes(x=pv)) + geom_histogram(color="white",fill="#4d4d4d",bins = 30,lwd=0.1) + th
dev.off()

View(hits_BR_noDMSO)

# #### combining p values across replicates using fisher's (combined_pval)
pooled_es<-hits_BR_noDMSO %>% dplyr::group_by(Bug_ID,Plate_no,compound,Phyla,Sp_short) %>% do(pool_es(.))
head(pooled_es)
nrow(pooled_es)

combined_pval<-hits_BR_noDMSO %>% dplyr::group_by(Bug_ID,Plate_no,compound,Phyla,Sp_short) %>% 
  summarise(combined_pv = sumlog(pv)[["p"]], count = length(pv), l2FC = log2(median(avFC))) #log 2 of median FC of replicates

es_pv_l2fc<-merge(pooled_es, combined_pval, by=c('Bug_ID','Plate_no','compound','Phyla','Sp_short'))
nrow(es_pv_l2fc) #5,925
head(es_pv_l2fc)


CairoSVG(file=paste("../Figures/pval_comb.svg", sep = ""), width = 3, height = 2, bg = "white")
es_pv_l2fc %>% ggplot(aes(x=combined_pv)) + geom_histogram(color="white",fill="#4d4d4d",bins = 30,lwd=0.1) + th
dev.off()

nrow(es_pv_l2fc %>% filter(Plate_no == 'plate1' & combined_pv < 0.05)) #631

### multiple hypotheses testing: error correction
p_bh<-es_pv_l2fc %>% dplyr::group_by(Bug_ID,Plate_no) %>%
  mutate(pv_bh = p.adjust(combined_pv,method = "BH"))
head(p_bh)
nrow(p_bh) #5925

hits<-p_bh %>% filter(l2FC != '-Inf')
head(p_bh)
nrow(hits) #1,050

CairoSVG(file=paste("../Figures/pval_bh.svg", sep = ""), width = 3, height = 2, bg = "white")
hits %>% filter(Plate_no == 'plate1') %>% ggplot(aes(x=pv_bh)) + geom_histogram(color="white",fill="#4d4d4d",bins = 30,lwd=0.1) + th
dev.off()

nrow(hits %>% filter(pv_bh < 0.05)) #562

CairoSVG(file=paste("../Figures/volcano_reps.svg", sep = ""), width = 5, height = 5, bg = "white")
EnhancedVolcano(hits,
                lab = hits$compound,
                x = 'l2FC',
                y = 'pv_bh',
                pCutoff = 0.05,
                FCcutoff = 0.3,
                legendPosition = 'bottom',
                labSize = 0
 )
dev.off()

CairoSVG(file=paste("../Figures/heatmap_reps.svg", sep = ""), width = 3, height =5, bg = "white")
hits %>% filter(pv_bh < 0.05 & (l2FC > 0.3 | l2FC < -0.3)) %>% ggplot(aes(x=Sp_short,y=compound)) + geom_tile(aes(fill=l2FC), color = "white", lwd = 0.1) + theme_bw() +
     th + theme(axis.text.x = element_text(angle = 90,hjust = 1), axis.ticks = element_blank(), panel.grid = element_blank(), axis.text.y = element_blank()) +
     scale_fill_gradientn(colours = wes_palette("Zissou1", 100, type = "continuous"))
dev.off()

hits %>% ggplot(aes(x=l2FC)) + geom_density()
head(hits)

sig_p<-hits %>% filter(pv_bh < 0.05 & (l2FC > 0.3 | l2FC < -0.3))
nrow(sig_p)
head(sig_p)

#cumulative distribution frequency plot
ggplot(hits, aes(x=l2FC)) + 
  stat_ecdf(geom = 'point') + scale_x_continuous(breaks = seq(-1.5,1.5,by=0.05))


h<-hits %>% filter(l2FC > 0)
nrow(h) #31

## combination compounds
comb<-p_bh %>% filter(l2FC != '-Inf' & compound %in% c('Advantame-comb','Caffeine','Duloxetine','Vanillin','Cetrizine','Risperidone','Ethinylestradiol','Ranitidine','Montelukast','Acetaminophen','Aripirazole','Ibuprofen'))
View(comb)
EnhancedVolcano(comb,
                lab = comb$compound,
                x = 'l2FC',
                y = 'p_bh',
                pCutoff = 0.05,
                FCcutoff = 0.3,
                legendPosition = 'bottom',
                labSize = 4
)
comb %>% filter(p_bh < 0.05 & (l2FC > 0.3 | l2FC < -0.3)) %>% ggplot(aes(x=Sp_short,y=compound)) + geom_tile(aes(fill=l2FC), color = "white", lwd = 0.1) + theme_bw() +
  th + theme(axis.text.x = element_text(angle = 90,hjust = 1), axis.ticks = element_blank(), panel.grid = element_blank()) +
  scale_fill_gradientn(colours = wes_palette("Zissou1", 100, type = "continuous"))

# plot of control and plate normalized AUCs
p_bh %>% filter(combined_pv < 0.05) %>% group_by(Bug_ID,Sp_short,Plate_no) %>% 
  ggplot(aes(x=l2FC,y=pooled_effect)) + geom_point(aes(color=combined_pv), size=1) + geom_abline(lwd=0.1) + stat_cor() + 
  scale_y_continuous(name = "pooled effect") + 
  scale_x_continuous(name = "l2FC") + scale_color_continuous() + 
  facet_wrap(c('Sp_short'), scales = "free") + th + theme_bw()

############## Bliss interactions
znormAUCs<-read.table(file = "../Figures/znormAUCs", header = TRUE, stringsAsFactors = FALSE)
nrow(znormAUCs %>% filter(compound != 'DMSO')) #37,507
head(znormAUCs)

# take log of AUCs
znormAUCs<-znormAUCs[,c(1:7,9,10)] %>% filter(compound != 'DMSO') %>% dplyr::group_by(Bug_ID,Replicate_no,Plate_no,compound,Phyla,Sp_short,col_name) %>% 
  summarise(lAUC = log(auc_l), lctrl = log(ctrl_mean))
head(znormAUCs)
nrow(znormAUCs) #37,507

#w/o plate6
#normAUCs food compounds
SFa<-znormAUCs %>% dplyr::group_by(Bug_ID,Replicate_no,Plate_no,compound,col_name,Phyla,Sp_short) %>% filter(Plate_no == "plate1" & Plate_no != 'plate6') %>%
  summarise(SFa = lAUC-lctrl)
head(SFa)  
nrow(SFa) #6015

#normAUCs drugs
SFq<-znormAUCs %>% dplyr::group_by(Bug_ID,Replicate_no,Plate_no,compound,Phyla,Sp_short) %>% filter(compound %in% c('Advantame-comb','Caffeine','Duloxetine','Vanillin')) %>%
  summarise(SFq = median(lAUC)-median(lctrl))
head(SFq)
nrow(SFq) #285

############################################################
#normAUCs food+drug
SFaq<-znormAUCs %>% dplyr::group_by(Bug_ID,Replicate_no,Plate_no,compound,col_name,Phyla,Sp_short) %>% filter(!(compound %in% c('Advantame-comb','Caffeine','Duloxetine','Vanillin','DMSO')) & !(Plate_no %in% c('plate1','plate6'))) %>%
  summarise(SFaq = lAUC-median(lctrl))
head(SFaq)
nrow(SFaq) #23,773
View(SFaq) 

#merge Advantame-comb
p<-SFa
head(p)
p[p=="plate1"]<-"plate1+A"
nrow(p) #6015

nrow(SFaq %>% filter(Plate_no == 'plate1+A')) #6,032

merge_advantame<-merge(SFaq, p, by=c("Bug_ID","Plate_no","Replicate_no","Phyla","Sp_short","compound","col_name"))
head(merge_advantame)
nrow(merge_advantame) #6000

#colnames(merge_advantame)<-c("Bug_ID","Plate_no","Replicate_no","Phyla","Sp_short","compound","col_name",
#                             "SFaq","SFaqc","SFa","SFac")

#merge caffeine

l<-SFa
l[l=="plate1"]<-"plate1+C"
nrow(l) #6015
head(l)

merge_caff<-merge(SFaq, l, by=c("Bug_ID","Plate_no","Replicate_no","Phyla","Sp_short","compound","col_name"))
head(merge_caff)
nrow(merge_caff) #5784


#colnames(merge_caff)<-c("Bug_ID","Plate_no","Replicate_no","Phyla","Sp_short","compound","col_name",
#                        "SFaq","SFaqc","SFa","SFac")

#merge duloxetine

o<-SFa
o[o=="plate1"]<-"plate1+D"
nrow(o) #6015

merge_dulox<-merge(SFaq, o, by=c("Bug_ID","Plate_no","Replicate_no","Phyla","Sp_short","compound","col_name"))
head(merge_dulox)
nrow(merge_dulox) #5869

#colnames(merge_dulox)<-c("Bug_ID","Plate_no","Replicate_no","Phyla","Sp_short","compound","col_name",
#                         "SFaq","SFaqc","SFa","SFac")

#merge vanillin

f<-SFa
f[f=="plate1"]<-"plate1+V"
nrow(f) #6015

merge_vani<-merge(SFaq, f, by=c("Bug_ID","Plate_no","Replicate_no","Phyla","Sp_short","compound","col_name"))
head(merge_vani)
nrow(merge_vani) #5994

#colnames(merge_vani)<-c("Bug_ID","Plate_no","Replicate_no","Phyla","Sp_short","compound","col_name",
#                        "SFaq","SFaqc","SFa","SFac")

merge_SFaq_SFa<-rbind(merge_advantame,merge_caff,merge_dulox,merge_vani)

head(merge_SFaq_SFa)
nrow(merge_SFaq_SFa) #23,647

head(SFq)
merge_SFaq_SFa_SFq<-merge(merge_SFaq_SFa, SFq, by=c("Bug_ID","Plate_no","Replicate_no","Phyla","Sp_short"))
head(merge_SFaq_SFa_SFq)
nrow(merge_SFaq_SFa_SFq) #23,629

head(merge_SFaq_SFa_SFq)
colnames(merge_SFaq_SFa_SFq)<-c("Bug_ID","Plate_no","Replicate_no","Phyla","Sp_short","compound","col_name",
                                "SFaq","SFa","comb_comp","SFq")
head(merge_SFaq_SFa_SFq)
nrow(merge_SFaq_SFa_SFq) #23,629

#write.table(bliss, file = "../Figures/bliss_scores", sep = "\t", quote = FALSE, row.names = FALSE)

bliss<-merge_SFaq_SFa_SFq %>% mutate(bliss_score = SFa+SFq-SFaq, bliss_indp = exp(SFa+SFq-SFaq))
nrow(bliss) #23,629
View(bliss)  


############statistical determination of synergy/antagonism

A_bliss<-bliss %>% filter_if(~is.numeric(.), all_vars(!is.infinite(.)))
nrow(A_bliss) #23,533
View(A_bliss)

A<- A_bliss %>% dplyr::group_by(Bug_ID,Plate_no,compound,Phyla,Sp_short) %>% 
  do(anov(.))
head(A)
nrow(A) #4,173
AB<-merge(A_bliss,A,by=c('Bug_ID','Plate_no','compound','Phyla','Sp_short'))
AB<-AB %>% dplyr::group_by(Bug_ID,Plate_no,compound) %>% mutate(comp_count = length(compound))
head(AB)
nrow(AB) #4,173
View(AB)


# density to determine bliss score distribution
CairoSVG(file="/Users/periwal/GrowthCurver/Figures/bliss_density.svg", width = 9, height = 4, bg = "white")
bliss %>% ggplot(aes(bliss_score)) + geom_density(aes(fill=comb_comp, color=comb_comp), alpha=0.5) +
  facet_wrap("Sp_short", nrow = 3, scales = "free") + theme_bw() +
  scale_fill_jama() + scale_color_jama() + th
dev.off()

#statistical determination of synergy/antagonism

View(bliss)

# expected vs observed viability
CairoSVG(file="../Figures/bliss_obs_exp.svg", width = 9, height = 6, bg = "white")
AB %>% dplyr::group_by(Bug_ID,Sp_short,compound,Replicate_no,Plate_no) %>% 
  ggplot(aes(x=SFaq,y=SFa+SFq)) + geom_point(aes(color=comb_comp), size=0.0001) + geom_abline() + theme_bw() +
  scale_y_continuous(name = "Expected viability") + 
  scale_x_continuous(name = "Observed viability") + scale_color_npg() + 
  facet_wrap(~Sp_short, scales = "free") + th + theme(legend.position = 'none')
dev.off()

sig_bliss<-AB %>% filter(pv < 0.05)
nrow(sig_bliss) #491
head(sig_bliss)

sig_bliss %>% ggplot(aes(x=bliss_score)) + geom_density()

#take median of bliss scores
med_sig<-sig_bliss %>% dplyr::group_by(Bug_ID,Plate_no,compound,Phyla,Sp_short) %>% summarise(med_SFaq = median(SFaq), med_SFa = median(SFa), med_SFq = median(SFq),
                                                                                     med_bliss = median(bliss_score), med_blissInd = median(bliss_indp),
                                                                                     pv = median(pv), comp_count = median(comp_count))

View(med_sig)
nrow(med_sig) #1162
#heatmap
CairoSVG(file="../Figures/sig_bliss_cutoffs.svg", width = 5, height = 3, bg = "white")
med_sig %>% ggplot(aes(x=med_bliss, color = Plate_no)) + scale_color_npg() + th +
  stat_ecdf(geom = 'point', size = 0.001) + scale_x_continuous(breaks = seq(-3,3,by=0.2), name = "bliss independence") + scale_y_continuous(name = "% counts") +
  theme_minimal() + geom_vline(xintercept = c(0.2,-.2), lwd=0.2, color='red') + theme(legend.position = "none")
dev.off()

CairoSVG(file="../Figures/sig_bliss_interactions.svg", width = 5, height = 4, bg = "white")
med_sig %>% filter(med_bliss > 0.2 | med_bliss < -0.2) %>% ggplot(aes(x=Sp_short,y=compound)) + facet_grid(~Plate_no, scales = 'free') + geom_tile(aes(fill=med_bliss)) + theme_bw() +
  th + theme(axis.text.x = element_text(angle = 90,hjust = 1), axis.ticks = element_blank(), panel.grid = element_blank(), axis.text.y = element_blank()) + 
  scale_fill_gradient2() + scale_y_discrete(name = "compounds") 
dev.off()


nrow(med_sig %>% filter(med_bliss > 0.2 | med_bliss < -0.2)) #50

#case example - synergy
exam<-AB %>% filter(Bug_ID == 'NT5011', Plate_no == 'plate1+D', compound == 'Iso-Steviol') %>% mutate(expSFaq = exp(SFaq), expSFa = exp(SFa), expSFq = exp(SFq), expected = exp(SFa+SFq))
head(exam)

CairoSVG(file="../Figures/Isosteviol+D_synergy.svg", width = 4, height = 3, bg = "white")
exam %>% reshape2::melt() %>% dplyr::group_by(Replicate_no) %>% filter(variable %in% c('expSFaq','expSFa','expSFq','expected')) %>% ggplot(aes(x=variable,y=value)) + 
  geom_boxplot(lwd=0.1,aes(fill=variable),alpha=0.3) + scale_x_discrete(labels=c("expSFa" = "array cmp", "expSFq" = "query cmp",
                                             "expSFaq" = "observed", "expected" = "expected")) + th + scale_fill_npg() +
  theme(legend.position = 'none', axis.title.x = element_blank()) + geom_jitter() +
  scale_y_continuous(name = 'Surviving fraction')
dev.off()



#case example - antagonism
exam<-AB %>% filter(Bug_ID == 'NT5023', Plate_no == 'plate1+C', compound == 'Sorbitol') %>% mutate(expSFaq = exp(SFaq), expSFa = exp(SFa), expSFq = exp(SFq), expected = exp(SFa+SFq))
head(exam)
View(exam)

CairoSVG(file="../Figures/Sorbitol+C_antagonism.svg", width = 4, height = 3, bg = "white")
exam %>% reshape2::melt() %>% dplyr::group_by(Replicate_no) %>% filter(variable %in% c('expSFaq','expSFa','expSFq','expected')) %>% ggplot(aes(x=variable,y=value)) + 
  geom_boxplot(lwd=0.1,aes(fill=variable),alpha=0.3) + scale_x_discrete(labels=c("expSFa" = "array cmp", "expSFq" = "query cmp",
                                                                                 "expSFaq" = "observed", "expected" = "expected")) + th + scale_fill_npg() +
  theme(legend.position = 'none', axis.title.x = element_blank()) + geom_jitter() +
  scale_y_continuous(name = 'Surviving fraction')
dev.off()


