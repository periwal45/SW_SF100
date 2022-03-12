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

