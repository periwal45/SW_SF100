############################################### Functions ##########################################################

# set global theme for all plots

th<-theme(plot.title = element_text(size = 12, face = "bold"),axis.title=element_text(size=12,color = "black"),
          axis.text.x = element_text(size=10, color = "black"),axis.text.y = element_text(size=10, color = "black"))

#melt dataframe of time and ODs (Used in Preprocess_Raw_Screens.R)
myfun <- function(x){
  melt(x,id.vars = "Time_h", variable.name = "wells", value.name= "OD")
}

# function for creating annotation from file name (Used in GC_processing.R)

create_annot <- function(input_file){
  
  elements<-strsplit(input_file, "_", perl = TRUE)
  ID<-elements[[1]][1]
  Replicate_no<-elements[[1]][2]
  Plate_no<-elements[[1]][3]
  sp<-strsplit(ID, "/", perl = TRUE)
  Bug_ID<-sp[[1]][1]
  Bug_ID
  annotation<-rbind(c(Bug_ID, Replicate_no, Plate_no))
  
  return(annotation)
  
}


# unpaired two samples welch t-test
unpairedT_pval<-function(bug_aucs){
 
  samples<-bug_aucs
  comp<-unique(samples$compound)
  comp
  
  ref_aucs<-(bug_aucs %>% filter(compound == 'DMSO'))$czscoreAUC
  ref_aucs
  
  out_df<-data.frame(stringsAsFactors = FALSE)
  
  for(i in 1:length(comp)){
    
    print(comp[i])
    print(samples$Bug_ID[1])
    print(samples$Plate_no[1])
    
    entry<-samples %>% filter(compound == comp[i])
    entry
    
    Bug_ID = entry$Bug_ID[1]
    Replicate_no = entry$Replicate_no[1]
    Plate_no = entry$Plate_no[1]
    compound = entry$compound[1]
    Phyla = entry$Phyla[1]
    Sp_short = entry$Sp_short[1]
    
    
    samp_aucs<-(samples %>% filter(compound == comp[i]))$czscoreAUC
    print(samp_aucs)
    print(length(samp_aucs))
    
    avFC = mean((samples %>% filter(compound == comp[i]))$FC) #mean for technical reps
    
    if(length(samp_aucs) > 1 & length(ref_aucs) > 1){
    
    test<-t.test(samp_aucs,ref_aucs,alternative = "two.sided", conf.level = 0.95)
    print(test)
    pv<-test$p.value
    conf.int.lower<-test$conf.int[1]
    conf.int.upper<-test$conf.int[2]
    sam<-length(samp_aucs)
    ref<-length(ref_aucs)
    t_stat<-test$statistic
    
    # estimate effect size from t-statistic (cohen's d, for correction and bias in small samples use: hedge's g)
    esize<-esc_t(test$statistic, grp1n = length(samp_aucs), grp2n = length(ref_aucs), es.type = "g")
    
    es<-esize$es
    gse<-esize$se
    lci<-esize$ci.lo
    uci<-esize$ci.hi
    
    z<-data.frame(Bug_ID,Replicate_no,Plate_no,compound,Phyla,Sp_short,pv,conf.int.lower,conf.int.upper,es,gse,lci,uci,ref,sam,t_stat,avFC)
    z
  
    out_df<-rbind(out_df,z)
    
    }else{
      
      print(paste0("Not enough observation:",comp[i]))
      next
    }
    
  }
  
  #View(out_df)
  return(out_df)
  
}


# pool effect sizes
pool_es<-function(bug_aucs){
  
  out_df<-data.frame(stringsAsFactors = FALSE)
  head(bug_aucs)    
  met<-metagen(data = bug_aucs, TE = es, seTE = gse, comb.fixed = FALSE, comb.random = TRUE, method.tau = "HE")
  #print(met)
  
  pval_es<-met$pval.random
  low.ci<-met$lower.random
  high.ci<-met$upper.random
  pooled_effect<-met$TE.random
  es_fixed<-met$TE.fixed
  studies_combined<-met$k
    
  z<-cbind(pval_es,low.ci,high.ci,pooled_effect,es_fixed,studies_combined)
  z
  
  #out_df<-rbind(out_df,z)
  return(z)
  
}

#adapted from https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0224137#sec020

t_distr<-function(y){ 
  
  n1 = length(y$SFa)
  n2 = length(y$SFq)
  n3 = length(y$SFaq)
  
  y1=mean(y$SFa)
  y2=mean(y$SFq)
  y3=mean(y$SFaq) 
  
  s1=var(y$SFa)*(n1-1)
  s2=var(y$SFq)*(n2-1)
  s3=var(y$SFaq)*(n3-1)
  
  sy=s1+s2+s3
  dft=n1+n2+n3-3
  denf=1/n1+1/n2+1/n3
  tss=(y1+y2-y3)/sqrt(sum(sy)/dft)/denf
  pv=2*(1-pt(abs(tss),df=dft))
  pvP=1-pt(tss,df=dft)
  
  #n - no of observations (replicates), sy - total sum of squares, dft,denf - df, tss - t-statistic, pv - p-value for Bliss independence hypothesis, pvP - One-sided p-value
  
  df<-data.frame("T-stat" = tss, "bliss_score" = y1+y2-y3, "bliss" = exp(y1+y2-y3), "pval" = pv)
  return(df)
  
}

anov<-function(y){

  res.aov<-aov(y$SFaq ~ y$SFa+y$SFq, data=y)
  pv<-summary(res.aov)[[1]][["Pr(>F)"]][1]
  df<-summary(res.aov)[[1]][["Df"]][3]
  f_val<-summary(res.aov)[[1]][["F value"]][1]
  
  if(length(pv) & length(df) & length(f_val) > 0){
    
    df<-data.frame("pv" = pv,"df" = df,"f_val" = f_val)
    
  }else{
    
    df<-data.frame("pv" = NULL,"df" = NULL,"f_val" = NULL)
    
  }

  return(df)
}
