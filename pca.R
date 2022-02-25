################### Quality check of fitted growth curves

View(fit_params_layout)
fit_params_layout<-merge(fit_params_layout,bug_desc, by="Bug_ID")
fit_params_layout$color<-paste0(fit_params_layout$Replicate_no,"_",fit_params_layout$tech)

#### Plot histogram of sigma values and pca of all params
bugs<-unique(fit_params_layout$Bug_ID)
bugs

outlier_wells_pca<-data.frame(stringsAsFactors = FALSE)

for(i in 1:length(bugs)){
  
  f<-fit_params_layout %>% filter(Bug_ID == bugs[i])
  plates<-unique(f$Plate_no)
  plates
  
  for(j in 1:length(plates)){
    
    p<-f %>% filter(Plate_no == plates[j])
    head(p)
    
    reps<-unique(p$Replicate_no)
    reps
    
    for(a in 1:length(reps)){
      
      r<- p %>% filter(Replicate_no == reps[a])
      head(r)
      
      rownames(r)<-r$col_name
      
      CairoSVG(file=paste("../Figures/hist/", bugs[i], plates[j],"sigma.svg"), width = 8, height = 4, bg = "white")
      print(p %>% dplyr::group_by(Bug_ID,Plate_no,compound,tech,Replicate_no) %>%
              ggplot(aes(x=sigma)) + geom_histogram(aes(fill=tech),color="black", lwd=0.2) +
              scale_fill_d3(palette = "category20") + facet_wrap("Replicate_no",scales = "free") +
              ggtitle(paste(bugs[i],p$Sp_short,plates[j],sep=":")) +
              th + theme_minimal())
      dev.off()
      
      
      pca.res <- prcomp(r %>% dplyr::select(k,n0,r,t_mid,t_gen,auc_l,auc_e,sigma), center=TRUE, scale=TRUE)
      head(pca.res)
      summary(pca.res)
      
      U<-pca.res$x
      head(U)
      
      outliers_Ustd<-apply(U, 2, function(x) which( (abs(x - median(x)) / mad(x)) > 6 ))
      outliers_Ustd
      
      if(length(outliers_Ustd) > 0){
        
        outs<-names(c(outliers_Ustd$PC1,outliers_Ustd$PC2))
        outs<-unique(outs)
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
        
      }
      
      df_pca<-as_tibble(list(PC1=pca.res$x[,1],
                             PC2=pca.res$x[,2],
                             samples = r$col_name))
      df_pca
      
      
      CairoSVG(file=paste("../Figures/pca/", bugs[i], plates[j],reps[a],"pca.svg"), width = 6, height = 4, bg = "white")
      print(df_pca %>% ggplot(aes(x=PC1,y=PC2, label=samples)) +
              geom_text(size = 3, aes(color=r$tech)) + scale_color_d3(palette = "category20") +
              ggtitle(paste(bugs[i],r$Sp_short,plates[j],r$Replicate_no,sep=":")))
      dev.off()
      
      
    }
    
  }
}

nrow(outlier_wells_pca) #994
View(outlier_wells_pca)

# remove outliers (noisy) wells (USE MANUAL CURATION)

# nrow(fit_params_layout) #41,472
# head(outlier_wells_pca)
# 
# df1<-inner_join(fit_params_layout,outlier_wells_pca)
# nrow(df1) #994
# View(df1)
# 
# fit_params_layout_filtered<-anti_join(fit_params_layout,df1)
# nrow(fit_params_layout_filtered) #40,478
# head(fit_params_layout_filtered)
# 
# # mark outlier wells
# fit_params_layout_filtered$outlier<-'no'
# 
# params_marked_outliers<-rbind(fit_params_layout_filtered,df1)
# params_marked_outliers<-params_marked_outliers[,c(1:6,13,16,20,23:27)]
# nrow(params_marked_outliers) #41,472
# View(params_marked_outliers)
# 
# write.table(params_marked_outliers, file = "../Figures/outliers_marked_params", sep = "\t", quote = FALSE, row.names = FALSE)
# 
# # group and filter outliers
# 
# params_marked_outliers<-read.table(file = "../Figures/outliers_marked_params", header = TRUE, sep = '\t')
# View(params_marked_outliers)
# dim(params_marked_outliers)
# 
# count(params_marked_outliers %>% filter(outlier == 'yes'))
# 
# params_outliers<-params_marked_outliers %>% dplyr::group_by(Bug_ID,Plate_no,compound,Replicate_no) %>% mutate(count = sum(outlier == 'yes'))
# dim(params_outliers)
# 
# #QC all oultier wells along with their technical and batch replicates
# params_outliers_filter<-params_outliers %>% filter(count >= 1)
# dim(params_outliers_filter)
# View(params_outliers_filter)
# 
# write.table(params_outliers_filter, file = "../Figures/params_pca_filtered.tab", sep = "\t", quote = FALSE, row.names = FALSE)
# 



