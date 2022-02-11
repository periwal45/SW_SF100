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

## fit and plot different distributions (read paper fitdistrplus, 2018)
### function to compute p-value from t-distribution

compute_pval_syn<-function(total_list_bliss){
  
  head(total_list_bliss)
  nrow(total_list_bliss)
  scores<-NULL
  scores<-(total_list_bliss %>% filter(outcome == 'reference'))$bliss_score

  f <- fitStudent(abs(scores))
  f
  
  total_list_bliss$pv <- pstudent(abs(total_list_bliss$bliss_score), f$estimate["nu"], f$estimate["mean"], f$estimate["sigma"])
  
  df<-data.frame(total_list_bliss)
  return(df)

}

compute_pval_bliss_ant<-function(total_list_bliss){
  
  head(total_list_bliss)
  nrow(total_list_bliss)
  scores<-NULL
  scores<-(total_list_bliss %>% filter(outcome == 'reference'))$bliss_score
  
  f <- fitStudent_ant(scores)
  f
  
  total_list_bliss$pv <- 1-pstudent(total_list_bliss$bliss_score, f$estimate["nu"], f$estimate["mean"], f$estimate["sigma"])
  
  df<-data.frame(total_list_bliss)
  return(df)
  
}

fitStudent <- function(y) {
  
  f <- NULL
  f <- fitdistrplus::fitdist(y, "student", start=list(nu=30, mean=mean(y), sigma=sd(y)))
  f
}

fitLNorm <- function(y) {
  
  f <- fitdistrplus::fitdist(y, "lnorm")
  f
}

compute_pval<-function(bug_aucs){
  
  head(bug_aucs)
  nrow(bug_aucs)
  
  normAUCs<-(bug_aucs %>% filter(label == 'reference'))$normAUC
  head(normAUCs)
  
  f <- fitStudent(normAUCs)
  #if (is.na(f$estimate["nu"])) {
  
  #bug_aucs$pv <- plnorm(bug_aucs$normAUC, f$estimate["meanlog"], f$estimate["sdlog"])
  bug_aucs$pv <- pstudent(bug_aucs$normAUC, f$estimate["nu"], f$estimate["mean"], f$estimate["sigma"])
  #bug_aucs$pv_upper <- pnorm(bug_aucs$normAUC, f$estimate["mean"], f$estimate["sd"], lower.tail = F)
  
  #} else {
  
  #bug_aucs$pv <- pstudent(bug_aucs$normAUC, f$estimate["nu"], f$estimate["mean"], f$estimate["sigma"])
  
  #bug_aucs$pv_upper <- pstudent(bug_aucs$normAUC, f$estimate["nu"], f$estimate["mean"], f$estimate["sigma"], lower.tail = F)
  #}
  
  df<-data.frame(bug_aucs)
  return(df)
  
}

########
dlnorm<-function(x, meanlog, sdlog) dlnorm(x, meanlog, sdlog, log = FALSE)
plnorm<-function(q, meanlog, sdlog, lower.tail, log.p) plnorm(q, meanlog, sdlog, lower.tail = TRUE, log.p = FALSE)
qlnorm<-function(p, meanlog, sdlog, lower.tail, log.p) qlnorm(p, meanlog, sdlog, lower.tail = TRUE, log.p = FALSE)
rlnorm<-function(n, meanlog, sdlog) rlnorm(n, meanlog, sdlog)


dstudent <- function(x, nu, mean, sigma) dt((x-mean)/sigma, nu, log=F)/sigma
pstudent <- function(q, nu, mean, sigma) pt((q-mean)/sigma, nu)
qstudent <- function(p, nu, mean, sigma) (qt(p, nu)*sigma+mean)

# dnormal<-function(x, mean, sd) dnorm(x, mean, sd, log = FALSE)
# pnormal<-function(p, mean, sd) pnorm(q, mean, sd, lower.tail = TRUE, log.p = FALSE)
# qnormal<-function(q, mean, sd) qnorm(p, mean, sd, lower.tail = TRUE, log.p = FALSE)
# rnormal<-function(r, mean, sd) rnorm(n, mean, sd)


fitStudentOrLnormal <- function(y) {
  f <- NULL
  try({
   f <- fitdistrplus::fitdist(y, "student", start=list(nu=30, mean=mean(y), sigma=sd(y)))
  })
  if (is.null(f) || f$estimate[1] > 100) {
    f <- fitdistrplus::fitdist(y, "lnorm")
  }
  f
}


fitCompare<-function(bug_aucs){
  
  normAUCs<-(bug_aucs %>% filter(label == 'reference'))$normAUC
  print(normAUCs)
  
  gofs<-data.frame(stringsAsFactors = FALSE)
  dist <- c("norm", "lnorm", "gamma", "student", "weibull", "cauchy")
  
  # Loop through your list of distributions
  for(i in 1:length(dist)){
    
    print(dist[i])
    
    try({
    if(dist[i] == "student"){
      
      x<-fitdistrplus::fitdist(normAUCs, "student", start=list(nu=30, mean=mean(normAUCs), sigma=sd(normAUCs)))
      l<-gofstat(x)
      chisqpv<-l$chisqpvalue
      if(length(chisqpv) > 0){
        gofs<-rbind(gofs,data.frame(dist[i],chisqpv))
      }
      
      }else{
        
      x <- fitdistrplus::fitdist(normAUCs, dist[i])
      l<-gofstat(x)
      chisqpv<-l$chisqpvalue
      if(length(chisqpv) > 0){
        gofs<-rbind(gofs,data.frame(dist[i],chisqpv))
      }
  
      }
    })
    
  }

  return(gofs)

}

#shapiro wilk test for normality

shap_wilk<-function(bug_aucs){
  
  nrow(bug_aucs)
  #View(bug_aucs)
  aucs<-bug_aucs$zscoreAUC
  
  shap<-shapiro.test(aucs)
  head(shap)
  
  df<-data.frame(shap$p.value)
  
  return(df)
}


# unpaired two samples welch t-test
unpairedT_pval<-function(bug_aucs){
 
  samples<-bug_aucs %>% filter(Compound!='DMSO')
  comp<-unique(samples$Compound)
  comp
  
  ref_aucs<-(bug_aucs %>% filter(Compound == 'DMSO'))$cnormAUC
  ref_aucs
  
  out_df<-data.frame(stringsAsFactors = FALSE)
  
  for(i in 1:length(comp)){
    
    print(comp[i])
    print(samples$Bug_ID[1])
    print(samples$Plate_no[1])
    
    entry<-samples %>% filter(Compound == comp[i])
    entry
    
    Bug_ID = entry$Bug_ID[1]
    Replicate_no = entry$Replicate_no[1]
    Plate_no = entry$Plate_no[1]
    Compound = entry$Compound[1]
    Phyla = entry$Phyla[1]
    Species = entry$Species[1]
    Sp_short = entry$Sp_short[1]
    
    samp_aucs<-(samples %>% filter(Compound == comp[i]))$cnormAUC
    print(samp_aucs)
    print(length(samp_aucs))
    
    if(length(samp_aucs) > 1 & length(ref_aucs) > 1){
    
    test<-t.test(samp_aucs,ref_aucs,alternative = "two.sided", conf.level = 0.95)
    print(test)
    pv<-test$p.value
    conf.int.lower<-test$conf.int[1]
    conf.int.upper<-test$conf.int[2]
    
    # estimate effect size from t-statistic (cohen's d, for correction and bias in small samples use: hedge's g)
    esize<-esc_t(test$statistic, grp1n = length(samp_aucs), grp2n = length(ref_aucs), es.type = "g")
    
    es<-esize$es
    gse<-esize$se
    lci<-esize$ci.lo
    uci<-esize$ci.hi
    
    z<-data.frame(Bug_ID,Replicate_no,Plate_no,Compound,Phyla,Species,Sp_short,pv,conf.int.lower,conf.int.upper,es,gse,lci,uci)
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
  met<-metagen(data = bug_aucs, TE = es, seTE = gse, comb.fixed = FALSE, comb.random = TRUE, hakn = TRUE, method.tau = "REML")
  #print(met)
  
  pval_es<-met$pval.random
  low.ci<-met$lower.random
  high.ci<-met$upper.random
  pooled_effect<-met$TE.random
  es_fixed<-met$TE.fixed
    
  z<-data.frame(pval_es,low.ci,high.ci,pooled_effect,es_fixed)
  z
  
  out_df<-rbind(out_df,z)
  return(out_df)
  
}

# n<-length(sample_wells$well)
# n
# 
# total_wells<-rbind(ref_wells,sample_wells)
# nrow(total_wells)
# head(total_wells)
# 
# out_df<-data.frame(stringsAsFactors = FALSE)
# 
# for(i in 1:length(total_wells$well)){
#   
#   Bug_ID<-total_wells[i,1]
#   Replicate_no<-total_wells[i,2]
#   Plate_no<-total_wells[i,3]
#   Drug_name<-total_wells[i,4]
#   well<-total_wells[i,5]
#   auc_l<-total_wells[i,6]
#   normAUC<-total_wells[i,7]
#   Phyla<-total_wells[i,8]
#   Species<-total_wells[i,9]
#   Sp_short<-total_wells[i,10]
#   label<-total_wells[i,11]
#   
#   t.val<-(as.numeric(normAUC)-mean_ref_normAUC)/(sd_sample_normAUC/sqrt(n))
#   t.val
#   p.val<-pt(-abs(t.val),df=n-1)
#   p.val
#   
#   df<-data.frame(Bug_ID,Replicate_no,Plate_no,Drug_name,well,auc_l,normAUC,t.val,p.val,Phyla,Species,Sp_short,label)
#   out_df<-rbind(out_df,df)
#   
# }
# 
# return(out_df)
# }