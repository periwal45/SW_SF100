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
library(paletteer)
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
              min_OD<-min((current_well[, col_name])[1,1])
              min_OD
              max_OD<-max(current_well[, col_name])
              max_OD
              
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
              annot_mod<-cbind(bug,rep,plate,col_name,mod)
              annot_mod
              
              fitted_values<-rbind(fitted_values, annot_mod)
              fitted_values
              
              fitted_param<-data.frame(cbind(bug,rep,plate,col_name,gc_fit$vals$k,gc_fit$vals$k_se,gc_fit$vals$k_p,gc_fit$vals$n0,gc_fit$vals$n0_se,gc_fit$vals$n0_p,gc_fit$vals$r,gc_fit$vals$r_se,gc_fit$vals$r_p,gc_fit$vals$sigma,gc_fit$vals$df,gc_fit$vals$t_mid,gc_fit$vals$t_gen,gc_fit$vals$auc_l,gc_fit$vals$auc_e,gc_fit$vals$note,max_OD,min_OD))
              colnames(fitted_param)<-c("Bug_ID","Replicate_no","Plate_no","well","k","k_se","k_p","n0","n0_se","n0_p","r","r_se","r_p","sigma","df","t_mid","t_gen","auc_l","auc_e","note","max_OD","min_OD")
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

head(gc_fit_params)
nrow(gc_fit_params) #should be equal to number of wells: 96 for 1 plate


write.table(fitted_values, file = "GC_fitted_model_values", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gc_fit_params, file = "GC_fit_params", sep = "\t", quote = FALSE, row.names = FALSE)


############################### Read fitted values
#read fitted values of models
fits<-data.frame(read.table("GC_fitted_model_values", header = TRUE, sep = '\t', stringsAsFactors = FALSE))
head(fits)
nrow(fits) #1,017,024

merge_plate_layout_fits<-merge(plate_layout, fits, by=c("Bug_ID","Plate_no","col_name"))
nrow(merge_plate_layout_fits) #1,017,024

View(merge_plate_layout_fits)

merge_plate_layout_fits<-merge(merge_plate_layout_fits,bug_desc, by="Bug_ID")
write.table(merge_plate_layout_fits, file = "../Figures/fitted_values", sep = "\t", quote = FALSE, row.names = FALSE)

############### plot raw and fitted ODs
merge_plate_layout_fits<-read.table(file = "../Figures/fitted_values", sep = "\t", header = TRUE)
merge_plate_layout_fits$color<-paste0(merge_plate_layout_fits$Replicate_no,"_",merge_plate_layout_fits$tech)
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
CairoSVG(file=paste("../Figures/NT24007_fig1.svg"), width = 4, height = 3, bg = "white")
merge_plate_layout_fits %>% dplyr::group_by(Bug_ID,Plate_no,compound,Replicate_no) %>% 
  filter(Bug_ID == 'NT24007' & compound %in% c('DMSO') & Plate_no %in% c('plate1') & tech == 'T1') %>%
  ggplot(aes(x=time, y=fitted)) + geom_line(aes(color=Replicate_no), size=0.4) + scale_color_uchicago() +
  geom_point(aes(y=OD, color=Replicate_no), size=0.01) + facet_wrap("compound", scales = "free") +  
  scale_x_continuous(name = "Time (Hours)", limits = c(0,25)) + scale_y_continuous(name="bg corrected OD") + 
  th + theme_minimal()
dev.off()
#'H1','H2','H3','H4','H5','H6','H7','H8','H9','H10','H11','H12'


############ QC: Data cleaning
# QC1: plate inspection after overlaying all replicates
# 1.NT5021 (plate1-all replicates) - kept, didn't remove
# 2.NT5022 (all)
# 3.NT5023

#read fitted values of models
fit_params<-data.frame(read.table("GC_fit_params", header = TRUE, sep = '\t', stringsAsFactors = FALSE))
head(fit_params)
nrow(fit_params) #43200

# removed not growing wells/replicates
# NT5076 not growing in any replicate (kick-out)
fit_params<-fit_params %>% filter(Bug_ID != 'NT5076')
nrow(fit_params) #41472
# already removed replicates
# NT5023 Rep-3 not growing
# NT5025, 5081 has only 15 points for rep3 (not sufficient) 

# QC2: individual noisy wells handpicked
noisy_wells_handpicked<-data.frame(read.table("../noisy_wells", header = TRUE, sep = '\t', stringsAsFactors = FALSE))
head(noisy_wells_handpicked)
dim(noisy_wells_handpicked) #533

fit_params_layout<-merge(plate_layout, fit_params, by=c("Bug_ID","Plate_no","col_name"))
head(fit_params_layout)

df1<-inner_join(fit_params_layout,noisy_wells_handpicked)
View(df1)
nrow(df1) #533

params_bugs_cmp_filtered<-anti_join(fit_params_layout,df1)
nrow(params_bugs_cmp_filtered) #40,939
head(params_bugs_cmp_filtered)

params_bugs_cmp_filtered<-merge(params_bugs_cmp_filtered,bug_desc,by='Bug_ID')
head(params_bugs_cmp_filtered)
nrow(params_bugs_cmp_filtered) #40,939

# plot of auc_l vs auc_e
params_bugs_cmp_filtered %>% group_by(Bug_ID,Sp_short,Replicate_no,Plate_no) %>% 
  ggplot(aes(x=auc_l,y=auc_e)) + geom_point(aes(color=Plate_no), size=1) + geom_abline(lwd=0.1) + stat_cor() + theme_bw() +
  scale_y_continuous(name = "Empirical AUC") + 
  scale_x_continuous(name = "Logistic AUC") + scale_color_npg() + 
  facet_wrap(~Sp_short, scales = "free") + th

# QC3: filter out questionable fits
head(params_bugs_cmp_filtered)

params<-params_bugs_cmp_filtered[,c(1:6,22,23:25,27,20)]
head(params)
dim(params) #40,939 X 12

nrow(params %>% filter(note != '')) #137

params<-params %>% filter(note == '') 
params<-params[,c(1:6,8:12)]
head(params)
nrow(params) #40802

write.table(params, file = '../Figures/cleaned_params', sep = '\t', row.names = FALSE)

############ 1. Normalization
params<-read.table(file = "../Figures/cleaned_params", header = TRUE, stringsAsFactors = FALSE)
dim(params) #40802 x 9
View(params)

##### plate-based normalization (median of controls)
PC<-params %>% filter(compound == 'DMSO') %>% dplyr::group_by(Bug_ID,Plate_no,Replicate_no,Phyla,Sp_short) %>%
  summarise(ctrl_med = round(median(auc_l),digits = 2))
View(PC)
dim(PC)

normAUCl_PC<-merge(params,PC,by=c("Bug_ID","Replicate_no","Plate_no","Phyla","Sp_short"))
dim(normAUCl_PC)
View(normAUCl_PC)

cnormAUCs<-normAUCl_PC %>% dplyr::group_by(Bug_ID,Plate_no,Replicate_no,col_name) %>%
  mutate(cAUC = round(auc_l/ctrl_med, digits = 2), FC = cAUC) #FC treatment vs control

head(cnormAUCs %>% filter(Plate_no == 'plate1'))

write.table(cnormAUCs, file = "../Figures/cnormAUCs", sep = '\t', row.names = FALSE)

# plotting density of normalized wells
CairoSVG(file=paste("../Figures/zscoreAUC.svg", sep = ""), width = 8.5, height = 5.5, bg = "white")
cnormAUCs %>%  dplyr::group_by(Bug_ID,Plate_no,Replicate_no) %>% filter(Plate_no != 'plate1') %>%
  ggplot(aes(cAUC)) + geom_density(aes(fill=Replicate_no,color=Replicate_no), alpha=0.3, lwd=0.3) +
  facet_wrap("Sp_short", scales = "free") +
  scale_fill_jama() + scale_color_jama() + th + theme_bw()
dev.off()

# plot of compound comparison with control
CairoSVG(file="../Figures/NT24007_fig1b.svg", width = 4, height = 3, bg = "white")
cnormAUCs %>% reshape2::melt() %>% filter(Bug_ID == 'NT5011' & compound %in% c('Duloxetine','Iso-Steviol') & variable == 'FC' & Plate_no %in% c('plate1','plate1+D')) %>% ggplot(aes(x=compound,y=value)) + 
  geom_boxplot(lwd=0.2,outlier.size = 0.1) + theme_minimal() +
  theme(axis.title.x = element_blank()) + geom_jitter(size = 0.1, aes(color=Replicate_no)) + scale_color_uchicago() +
  scale_y_continuous(name = 'Surviving fraction') + th
dev.off()


# # plot of control and plate normalized AUCs
# znormAUCs %>% group_by(Bug_ID,Sp_short,Replicate_no,Plate_no) %>% 
#   ggplot(aes(x=czscoreAUC,y=pzscoreAUC)) + geom_point(aes(color=Plate_no), size=1) + geom_abline(lwd=0.1) + stat_cor() + 
#   scale_y_continuous(name = "ctrl norm AUC") + 
#   scale_x_continuous(name = "pl norm AUC") + scale_color_npg() + 
#   facet_grid(c('Plate_no','Sp_short'), scales = "free") + th + theme_bw()
# 

######## 2. Hit selection (unpaired (unequal no of observations) t-test / welch two sample t-test)
cnormAUCs<-read.table(file = "../Figures/cnormAUCs", header = TRUE, stringsAsFactors = FALSE)
dim(cnormAUCs)
View(cnormAUCs)

# t-test bug wise and replicate wise (plate 1 is all sweeteners)
# hits_BR<-cnormAUCs %>% filter(Plate_no == 'plate1') %>% dplyr::group_by(Bug_ID, Plate_no, Replicate_no) %>% 
#   do(unpairedT_pval(.))
# View(hits_BR)
# nrow(hits_BR) #3051


# t-test all replicates pooled
hits_BR_norep<-cnormAUCs %>% filter(Plate_no == 'plate1') %>% dplyr::group_by(Bug_ID, Plate_no) %>% 
  do(unpairedT_pval(.))
View(hits_BR_norep)
nrow(hits_BR_norep) #1075



#hits_BR_noDMSO<-hits_BR_norep %>% filter(compound != 'DMSO')

CairoSVG(file=paste("../Figures/pval_dist.svg", sep = ""), width = 3, height = 2, bg = "white")
hits_BR_norep %>% ggplot(aes(x=pv)) + geom_histogram(color="white",fill="#4d4d4d",bins = 30,lwd=0.1) + th
dev.off()

# #### combining p values across replicates using fisher's (combined_pval)
# no need to combine effect sizes as it's already computed from pooled replicates
# pooled_es<-hits_BR_norep %>% dplyr::group_by(Bug_ID,Plate_no,compound,Phyla,Sp_short) %>% do(pool_es(.))
# head(pooled_es)
# nrow(pooled_es) #1075

# combined_pval<-hits_BR_norep %>% dplyr::group_by(Bug_ID,Plate_no,compound,Phyla,Sp_short) %>% 
#   summarise(combined_pv = sumlog(pv)[["p"]], count = length(pv), l2FC = log2(median(avFC))) #log 2 of median FC of replicates

# log2FC
log2FC<-hits_BR_norep %>% dplyr::group_by(Bug_ID,Plate_no,compound,Phyla,Sp_short) %>% 
  mutate(l2FC = round(log2(median(avFC)),digits = 2)) #log 2 of median FC of replicates
View(log2FC)
dim(log2FC)

nrow(log2FC %>% filter(Plate_no == 'plate1' & pv < 0.05)) #334

### multiple hypotheses testing: error correction
p_bh<-log2FC %>% dplyr::group_by(Bug_ID,Plate_no) %>%
  mutate(pv_bh = p.adjust(pv,method = "BH"))
View(p_bh)
nrow(p_bh) #1075

write.table(p_bh, file = "t-test_cnorms", sep = '\t', row.names = FALSE)

hits<-p_bh %>% filter(l2FC != '-Inf')
View(hits)
nrow(hits) #1,075

# CairoSVG(file=paste("../Figures/pval_bh.svg", sep = ""), width = 3, height = 2, bg = "white")
# hits %>% filter(Plate_no == 'plate1' & compound != 'DMSO') %>% ggplot(aes(x=pv_bh)) + geom_histogram(color="white",fill="#4d4d4d",bins = 30,lwd=0.1) + th
# dev.off()

filter_hits<-hits %>% filter(pv_bh < 0.05 & (l2FC >= 0.32 | l2FC <= -0.32)) #33

CairoSVG(file=paste("../Figures/volcano_reps.svg", sep = ""), width = 7, height = 5, bg = "white")
EnhancedVolcano(hits,
                lab = hits$compound,
                x = 'l2FC',
                y = 'pv_bh',
                pCutoff = 0.05,
                FCcutoff = 0.31,
                legendPosition = 'right',
                legendLabels=c('not sig.',bquote(log[2]~'fold change'),'p-value',
                               bquote('p-value &' ~log[2]~'fold change')),
                labSize = 4,
                pointSize = 0.5,
                col=c('#737373', '#02818a', '#fec44f', '#67001f'),
                colAlpha = 0.8,
                xlim = c(-1,1),
                ylim = c(0,12),
                legendLabSize = 12,
                legendIconSize = 3.0,
                axisLabSize = 12,
                title = "",
                subtitle = "",
                caption = ""
                )
dev.off()



CairoSVG(file=paste("../Figures/heatmap_reps.svg", sep = ""), width = 10, height = 3.5, bg = "white")
filter_hits %>% ggplot(aes(x=Sp_short,y=compound)) + geom_tile(aes(fill=l2FC), color = "white", lwd = 0.1) + theme_bw() +
     th + theme(axis.text.x = element_text(angle = 90,hjust = 1), axis.ticks = element_blank(), panel.grid = element_blank()) +
  scale_fill_gradientn(colors = wes_palette("Zissou1", 10, type = "continuous"), name = bquote(log[2]~'fold change')) + coord_flip()
dev.off()

hits %>% ggplot(aes(x=l2FC)) + geom_density()
head(hits)

#cumulative distribution frequency plot
ggplot(hits, aes(x=l2FC)) + 
  stat_ecdf(geom = 'point') + scale_x_continuous(breaks = seq(-1.5,1.5,by=0.05))

# plot of control and plate normalized AUCs
p_bh %>% filter(combined_pv < 0.05) %>% group_by(Bug_ID,Sp_short,Plate_no) %>% 
  ggplot(aes(x=l2FC,y=pooled_effect)) + geom_point(aes(color=combined_pv), size=1) + geom_abline(lwd=0.1) + stat_cor() + 
  scale_y_continuous(name = "pooled effect") + 
  scale_x_continuous(name = "l2FC") + scale_color_continuous() + 
  facet_wrap(c('Sp_short'), scales = "free") + th + theme_bw()




# # combination compound alone
# head(cnormAUCs)
# comb_cnormAUCs<-cnormAUCs
# 
# c_hits_BR<-comb_cnormAUCs %>% filter(compound %in% c('Advantame-comb','Caffeine','Duloxetine','Vanillin','DMSO') & Plate_no != 'plate1') %>% dplyr::group_by(Bug_ID, Plate_no, Replicate_no) %>% 
#   do(unpairedT_pval(.))
# head(c_hits_BR)
# nrow(c_hits_BR) #644
# 
# c_combined_FC<-c_hits_BR %>% dplyr::group_by(Bug_ID,Plate_no,compound,Phyla,Sp_short) %>% 
#   mutate(l2FC = round(log2(median(avFC)),digits = 2))
# c_p_bh<-c_combined_FC %>% dplyr::group_by(Bug_ID,Plate_no) %>%
#   mutate(pv_bh = p.adjust(pv,method = "BH"))
# 
# head(c_p_bh)
# nrow(c_p_bh)
# 
# cc<-merge(c_p_bh,c_hits_BR,by=c('Bug_ID','compound','Plate_no','Phyla','Sp_short'))
# 
# head(cc)
# 
# c_p_bh %>% filter(pv_bh < 0.05 | compound == 'DMSO') %>% dplyr::group_by(Bug_ID,Replicate_no,Plate_no) %>% ggplot(aes(x=compound, y=log2(avFC))) + geom_boxplot(aes(fill=compound), lwd=0.1, alpha=0.8) +
#   facet_wrap(~Bug_ID, scales = 'free') + scale_fill_nejm() + th + theme(axis.text.x = element_blank()) + geom_jitter(size=0.5) +
#   scale_color_nejm()
# 


############## Bliss interactions
cnormAUCs<-read.table(file = "../Figures/cnormAUCs", header = TRUE, stringsAsFactors = FALSE)
#nrow(znormAUCs %>% filter(compound != 'DMSO')) #37,507
View(cnormAUCs)

head(cnormAUCs)

# take log of AUCs
# blAUCs<-cnormAUCs %>% filter(compound != 'DMSO') %>% dplyr::group_by(Bug_ID,Replicate_no,Plate_no,compound,Phyla,Sp_short,col_name) %>%
#   mutate(lAUC = log(FC))
# View(blAUCs)
# nrow(blAUCs) #37,375

#w/o plate6
#normAUCs sweetener alone
SFa<-cnormAUCs %>% filter(Plate_no == "plate1" & compound != 'DMSO') %>% dplyr::group_by(Bug_ID,Replicate_no,Plate_no,compound,col_name,Phyla,Sp_short)  %>%
  summarise(SFa = FC)
head(SFa)  
nrow(SFa) #5995

#median of combination compound wells
SFq<-cnormAUCs %>% filter(compound %in% c('Advantame-comb','Caffeine','Duloxetine','Vanillin')) %>% dplyr::group_by(Bug_ID,Replicate_no,Plate_no,compound,Phyla,Sp_short) %>% 
  summarise(SFq = median(FC))
View(SFq)
nrow(SFq) #285

############################################################
#normAUCs sweetener+comb compound
SFaq<-cnormAUCs %>% filter(!(compound %in% c('Advantame-comb','Caffeine','Duloxetine','Vanillin','DMSO')) & !(Plate_no %in% c('plate1','plate6'))) %>% 
  dplyr::group_by(Bug_ID,Replicate_no,Plate_no,compound,col_name,Phyla,Sp_short) %>% 
  summarise(SFaq = FC)
head(SFaq)
nrow(SFaq) #23,677
View(SFaq) 

#merge Advantame-comb
p<-SFa
head(p)
p[p=="plate1"]<-"plate1+A"
nrow(p) #5995

nrow(SFaq %>% filter(Plate_no == 'plate1+A')) #6,016

merge_advantame<-merge(SFaq, p, by=c("Bug_ID","Plate_no","Replicate_no","Phyla","Sp_short","compound","col_name"))
head(merge_advantame)
nrow(merge_advantame) #5965


df1<-anti_join(SFaq %>% filter(Plate_no == 'plate1+A'),merge_advantame)
head(df1)
nrow(df1) #51


#colnames(merge_advantame)<-c("Bug_ID","Plate_no","Replicate_no","Phyla","Sp_short","compound","col_name",
#                             "SFaq","SFaqc","SFa","SFac")

#merge caffeine

l<-SFa
l[l=="plate1"]<-"plate1+C"
nrow(l) #5995
head(l)

merge_caff<-merge(SFaq, l, by=c("Bug_ID","Plate_no","Replicate_no","Phyla","Sp_short","compound","col_name"))
head(merge_caff)
nrow(merge_caff) #5742


#colnames(merge_caff)<-c("Bug_ID","Plate_no","Replicate_no","Phyla","Sp_short","compound","col_name",
#                        "SFaq","SFaqc","SFa","SFac")

#merge duloxetine

o<-SFa
o[o=="plate1"]<-"plate1+D"
nrow(o) #5995

merge_dulox<-merge(SFaq, o, by=c("Bug_ID","Plate_no","Replicate_no","Phyla","Sp_short","compound","col_name"))
head(merge_dulox)
nrow(merge_dulox) #5833

#colnames(merge_dulox)<-c("Bug_ID","Plate_no","Replicate_no","Phyla","Sp_short","compound","col_name",
#                         "SFaq","SFaqc","SFa","SFac")

#merge vanillin

f<-SFa
f[f=="plate1"]<-"plate1+V"
nrow(f) #5995

merge_vani<-merge(SFaq, f, by=c("Bug_ID","Plate_no","Replicate_no","Phyla","Sp_short","compound","col_name"))
head(merge_vani)
nrow(merge_vani) #5933

#colnames(merge_vani)<-c("Bug_ID","Plate_no","Replicate_no","Phyla","Sp_short","compound","col_name",
#                        "SFaq","SFaqc","SFa","SFac")

merge_SFaq_SFa<-rbind(merge_advantame,merge_caff,merge_dulox,merge_vani)

head(merge_SFaq_SFa)
nrow(merge_SFaq_SFa) #23,473

head(SFq)
merge_SFaq_SFa_SFq<-merge(merge_SFaq_SFa, SFq, by=c("Bug_ID","Plate_no","Replicate_no","Phyla","Sp_short"))
head(merge_SFaq_SFa_SFq)
nrow(merge_SFaq_SFa_SFq) #23,455

head(merge_SFaq_SFa_SFq)
colnames(merge_SFaq_SFa_SFq)<-c("Bug_ID","Plate_no","Replicate_no","Phyla","Sp_short","compound","col_name",
                                "SFaq","SFa","comb_comp","SFq")
View(merge_SFaq_SFa_SFq)
nrow(merge_SFaq_SFa_SFq) #23,455

write.table(bliss, file = "../Figures/bliss_scores", sep = "\t", quote = FALSE, row.names = FALSE)


bliss<-merge_SFaq_SFa_SFq %>% mutate(null_hyp = log(SFa)+log(SFq)-log(SFaq), exp_eff = SFa*SFq)
nrow(bliss) #23,455
head(bliss)
View(bliss)  

write.table(bliss, file = "../Figures/bliss_scores", sep = "\t", quote = FALSE, row.names = FALSE)

# density to determine bliss score distribution
CairoSVG(file="../Figures/bliss_density.svg", width = 9, height = 4, bg = "white")
bliss %>% ggplot(aes(null_hyp)) + geom_density(aes(fill=comb_comp, color=comb_comp), alpha=0.5) +
  facet_wrap("Sp_short", nrow = 3, scales = "free") + theme_bw() +
  scale_fill_uchicago() + scale_color_uchicago() + th
dev.off()


# density to determine bliss score distribution (single bug)
CairoSVG(file = "../Figures/bliss_density_NT24007.svg", width = 4, height = 2, bg = "white")
bliss %>% filter(Bug_ID == 'NT24007') %>% ggplot(aes(null_hyp)) + geom_density(aes(fill=comb_comp), lwd=0.1,alpha=0.6) +
  theme_bw() + scale_fill_manual(values = c('#67001f','#737373', '#02818a', '#fec44f')) + th
dev.off()

############statistical determination of synergy/antagonism

A<-bliss %>% dplyr::group_by(Bug_ID,Plate_no,compound,Phyla,Sp_short,comb_comp) %>% 
  do(t_distr(.))
View(A)
nrow(A) #4,200

write.table(A, file = "../Figures/Sig_bliss_scores", sep = "\t", quote = FALSE, row.names = FALSE)

# AB<-merge(bliss,A,by=c('Bug_ID','Plate_no','compound','Phyla','Sp_short'))
# AB<-AB %>% dplyr::group_by(Bug_ID,Plate_no,compound) %>% mutate(comp_count = length(compound))
head(A)
# nrow(AB) #4,173
# View(AB)

#statistical determination of synergy/antagonism

# expected vs observed viability
CairoSVG(file="../Figures/bliss_obs_exp.svg", width = 12, height = 8, bg = "white")
A %>% dplyr::group_by(Bug_ID,Sp_short,compound,Plate_no) %>% 
  ggplot(aes(x=mean_SFaq,y=mean_SFa+mean_SFq)) + geom_point(aes(color=comb_comp), size=0.0001) + geom_abline() + theme_bw() +
  scale_y_continuous(name = "Expected viability") + 
  scale_x_continuous(name = "Observed viability") + scale_color_npg() + 
  facet_wrap(~Sp_short, scales = "free") + th 
dev.off()

sig_bliss<-A %>% filter(pval < 0.05)
nrow(sig_bliss) #888
View(sig_bliss)

# # #take median of bliss scores
# med_sig<-sig_bliss %>% dplyr::group_by(Bug_ID,Plate_no,compound,Phyla,Sp_short) %>% summarise(med_SFaq = median(SFaq), med_SFa = median(SFa), med_SFq = median(SFq),
#                                                                                       med_bliss = median(bliss_null),
#                                                                                       pv = median(pv), comp_count = median(comp_count))
# 
# View(med_sig)
# nrow(med_sig) #1162
#heatmap
CairoSVG(file="../Figures/sig_bliss_cutoffs.svg", width = 5, height = 3, bg = "white")
sig_bliss %>% ggplot(aes(x=bliss_null, color = Plate_no)) + scale_color_npg() + th +
  stat_ecdf(geom = 'point', size = 0.001) + scale_x_continuous(breaks = seq(-3,3,by=0.2), name = "bliss independence") + scale_y_continuous(name = "% counts") +
  theme_minimal() + geom_vline(xintercept = c(0.2,-0.2), lwd=0.2, color='red') 
dev.off()

pal<-paletteer_c("ggthemes::Red-Blue Diverging", 100)
pal

CairoSVG(file="../Figures/sig_bliss_interactions.svg", width = 7, height = 7, bg = "white")
sig_bliss %>% filter(bliss_null > 0.2 | bliss_null < -0.2 & !(compound %in% c('MIN-3-035','MIN-0-60','MIN-3-100'))) %>% ggplot(aes(x=Sp_short,y=compound)) + facet_grid(~Plate_no, scales = 'free') + geom_tile(aes(fill=bliss_null)) + theme_bw() +
  th + theme(axis.text.x = element_text(angle = 90,hjust = 1), axis.ticks = element_blank(), panel.grid = element_blank()) + 
  scale_fill_gradientn(colors = pal) + scale_y_discrete(name = "compounds") 
dev.off()

wes_palette("FantasticFox1", 10, type = "continuous")

View(sig_bliss %>% filter(bliss_null > 0.2 | bliss_null < -0.2)) #50


#case example - synergy
exam<-bliss %>% filter(Bug_ID == 'NT5011', Plate_no == 'plate1+D', compound == 'Iso-Steviol')
View(exam)

write.table(exam, file = 'Synergy-I+D_NT5011', quote = FALSE, row.names = FALSE, sep = '\t')

CairoSVG(file="../Figures/Isosteviol+D_synergy.svg", width = 4, height = 3, bg = "white")
exam %>% reshape2::melt() %>% dplyr::group_by(Replicate_no) %>% filter(variable %in% c('SFaq','SFa','SFq','exp_eff')) %>% ggplot(aes(x=variable,y=value)) + 
  geom_boxplot(lwd=0.2, outlier.size = 0.1) + theme_minimal() + scale_x_discrete(labels=c("SFa" = "Iso-Setviol", "SFq" = "Duloxetine",
                                             "SFaq" = "observed", "exp_eff" = "expected")) + scale_color_uchicago() +
  theme(legend.position = "none",axis.title.x = element_blank()) + geom_jitter(size = 0.1, aes(color=Replicate_no)) + th + 
  scale_y_continuous(name = 'Surviving fraction')
dev.off()



#case example - antagonism
ant<-bliss %>% filter(Bug_ID == 'NT5022', Plate_no == 'plate1+V', compound == 'Saccharin') 
head(ant)
View(ant)

CairoSVG(file="../Figures/Saccharin+V_antagonism.svg", width = 6, height = 4, bg = "white")
ant %>% reshape2::melt() %>% dplyr::group_by(Replicate_no) %>% filter(variable %in% c('SFaq','SFa','SFq','exp_eff')) %>% ggplot(aes(x=variable,y=value)) + 
  geom_boxplot(lwd=0.2, outlier.size = 0.1) + theme_minimal() + scale_x_discrete(labels=c("SFa" = "Saccharin", "SFq" = "Vanillin",
                                                                        "SFaq" = "observed", "exp_eff" = "expected")) + scale_color_uchicago() +
  theme(legend.position = 'none', axis.title.x = element_blank()) + geom_jitter(size = 0.1, aes(color=Replicate_no)) + th +
  scale_y_continuous(name = 'Surviving fraction')
dev.off()


