

params<-read.table(file = '../Figures/cleaned_params', header = TRUE, sep = '\t', stringsAsFactors = FALSE)
params<-params %>% filter(Plate_no == 'plate6')
head(params)
dim(params) #6,882 X 9

cnormAUCl_z<-params %>% filter(compound == 'DMSO') %>% dplyr::group_by(Bug_ID,Plate_no,Replicate_no,Phyla,Sp_short) %>%
  summarise(ctrl_mean = mean(auc_l), std_ctrl = std(auc_l))
head(cnormAUCl_z)
dim(cnormAUCl_z)

normAUCl_z<-merge(params,cnormAUCl_z,by=c("Bug_ID","Replicate_no","Plate_no","Phyla","Sp_short"))
dim(normAUCl_z)
head(normAUCl_z)

znormAUCs<-normAUCl_z %>% dplyr::group_by(Bug_ID,Plate_no,Replicate_no,col_name) %>%
  mutate(czscoreAUC = (auc_l-ctrl_mean)/std_ctrl, FC = auc_l/ctrl_mean) #FC treatment vs control

head(znormAUCs)

znormAUCs %>%  dplyr::group_by(Bug_ID,Plate_no,Replicate_no) %>% 
  ggplot(aes(czscoreAUC)) + geom_density(aes(fill=Replicate_no,color=Replicate_no), alpha=0.3, lwd=0.3) +
  facet_wrap("Sp_short", scales = "free") +
  scale_fill_jama() + scale_color_jama() + th + theme_bw()

View(znormAUCs)
# t-test bug wise and replicate wise
hits_BR<-znormAUCs %>% filter(compound %in% c('DMSO','Ace-K','Aspartame','Sucralose','Acetaminophen','Aripirazole','Cetrizine','Ethinylestradiol','Ibuprofen','Montelukast','Ranitidine','Risperidone')) %>% 
  dplyr::group_by(Bug_ID, Replicate_no) %>% 
  do(unpairedT_pval(.))
head(hits_BR)
nrow(hits_BR) #864

hits_BR_noDMSO<-hits_BR %>% filter(compound != 'DMSO')

CairoSVG(file=paste("../Figures/pval_dist.svg", sep = ""), width = 3, height = 2, bg = "white")
hits_BR_noDMSO %>% ggplot(aes(x=pv)) + geom_histogram(color="white",fill="#4d4d4d",bins = 30,lwd=0.1) + th
dev.off()

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

nrow(es_pv_l2fc %>% filter(combined_pv < 0.05)) #193

### multiple hypotheses testing: error correction
p_bh<-es_pv_l2fc %>% dplyr::group_by(Bug_ID) %>%
  mutate(pv_bh = p.adjust(combined_pv,method = "BH"))
head(p_bh)
nrow(p_bh) #275

hits<-p_bh %>% filter(l2FC != '-Inf')
head(p_bh)
nrow(hits) #275

CairoSVG(file=paste("../Figures/pval_bh.svg", sep = ""), width = 3, height = 2, bg = "white")
hits %>% ggplot(aes(x=pv_bh)) + geom_histogram(color="white",fill="#4d4d4d",bins = 30,lwd=0.1) + th
dev.off()

nrow(hits %>% filter(pv_bh < 0.05)) #185

CairoSVG(file=paste("../Figures/form_volcano_reps.svg", sep = ""), width = 7, height = 6, bg = "white")
EnhancedVolcano(hits,
                lab = hits$compound,
                x = 'l2FC',
                y = 'pv_bh',
                pCutoff = 0.05,
                FCcutoff = 0.3,
                legendPosition = 'bottom',
                labSize = 4
)
dev.off()

CairoSVG(file=paste("../Figures/form_heatmap_reps.svg", sep = ""), width = 6, height =4, bg = "white")
hits %>% filter(pv_bh < 0.05 & (l2FC > 0.3 | l2FC < -0.3)) %>% ggplot(aes(x=Sp_short,y=compound)) + geom_tile(aes(fill=l2FC), color = "white", lwd = 0.1) + theme_bw() +
  th + theme(axis.text.x = element_text(angle = 90,hjust = 1), axis.ticks = element_blank(), panel.grid = element_blank()) +
  scale_fill_gradientn(colours = wes_palette("Zissou1", 100, type = "continuous"))
dev.off()

View(hits_BR_noDMSO)

# plate6
#normAUCs food compounds
p6<-znormAUCs %>% dplyr::group_by(Bug_ID,Replicate_no,Plate_no,compound,col_name,Phyla,Sp_short) %>% filter(Plate_no == 'plate6')
View(p6)  

p6_SFa<-p6 %>% filter(compound %in% c('Ace-K','Aspartame','Sucralose')) %>% dplyr::group_by(Bug_ID,Replicate_no,Plate_no,compound,Phyla,Sp_short) %>% 
  summarise(SFa = median(lAUC)-median(lctrl))
head(p6_SFa)  
nrow(p6_SFa) #216

#normAUCs drugs
p6_SFq<-p6 %>% filter(compound %in% c('Acetaminophen','Aripirazole','Cetrizine','Ethinylestradiol','Ibuprofen','Montelukast','Ranitidine','Risperidone')) %>%
  dplyr::group_by(Bug_ID,Replicate_no,Plate_no,compound,Phyla,Sp_short) %>% 
  summarise(SFq = median(lAUC)-median(lctrl))
head(p6_SFq)
nrow(p6_SFq) #576

sd<-c('Ace-K','Aspartame','Sucralose','Acetaminophen','Aripirazole','Cetrizine','Ethinylestradiol','Ibuprofen','Montelukast','Ranitidine','Risperidone')
############################################################
#normAUCs food+drug
p6_SFaq<-p6 %>% filter(!compound %in% sd) %>% dplyr::group_by(Bug_ID,Replicate_no,Plate_no,compound,col_name,Phyla,Sp_short) %>%
  summarise(SFaq = lAUC-lctrl)
head(p6_SFaq)
nrow(p6_SFaq) #864
View(SFaq)

# combinations
#AceK+ACetaminophen
acek_acet<-merge((p6_SFaq %>% filter(compound == 'Ace-K+acetaminophen')),(p6_SFa %>% filter(compound == 'Ace-K')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
nrow(acek_acet)
head(acek_acet)
colnames(acek_acet)<-c("Bug_ID","Replicate_no","Phyla","Sp_short","Plate_no","comb_cmp","well",
                       "SFaq","cmp_a","SFa")

acek_acet<-merge(acek_acet,(p6_SFq %>% filter(compound == 'Acetaminophen')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
head(acek_acet)


# AceK+Cetrizine
acek_cet<-merge((p6_SFaq %>% filter(compound == 'Ace-K+cetrizine')),(p6_SFa %>% filter(compound == 'Ace-K')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
nrow(acek_cet)
head(acek_cet)
colnames(acek_cet)<-c("Bug_ID","Replicate_no","Phyla","Sp_short","Plate_no","comb_cmp","well",
                      "SFaq","cmp_a","SFa")

acek_cet<-merge(acek_cet,(p6_SFq %>% filter(compound == 'Cetrizine')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
head(acek_cet)

# AceK+Risperidone
acek_risp<-merge((p6_SFaq %>% filter(compound == 'Ace-K+Risperidone')),(p6_SFa %>% filter(compound == 'Ace-K')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
nrow(acek_risp)
head(acek_risp)
colnames(acek_risp)<-c("Bug_ID","Replicate_no","Phyla","Sp_short","Plate_no","comb_cmp","well",
                       "SFaq","cmp_a","SFa")

acek_risp<-merge(acek_risp,(p6_SFq %>% filter(compound == 'Risperidone')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
head(acek_risp)


# AceK+Ibuprofen
acek_ibu<-merge((p6_SFaq %>% filter(compound == 'Ace-K+Ibuprofen')),(p6_SFa %>% filter(compound == 'Ace-K')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
nrow(acek_ibu)
head(acek_ibu)
colnames(acek_ibu)<-c("Bug_ID","Replicate_no","Phyla","Sp_short","Plate_no","comb_cmp","well",
                      "SFaq","cmp_a","SFa")

acek_ibu<-merge(acek_ibu,(p6_SFq %>% filter(compound == 'Ibuprofen')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
head(acek_ibu)


# Asp+Aripiprazole
asp_ari<-merge((p6_SFaq %>% filter(compound == 'Asp+Aripiprazole')),(p6_SFa %>% filter(compound == 'Aspartame')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
nrow(asp_ari)
head(asp_ari)
colnames(asp_ari)<-c("Bug_ID","Replicate_no","Phyla","Sp_short","Plate_no","comb_cmp","well",
                     "SFaq","cmp_a","SFa")

asp_ari<-merge(asp_ari,(p6_SFq %>% filter(compound == 'Aripirazole')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
head(asp_ari)

# Asp+Cetrizine
asp_cet<-merge((p6_SFaq %>% filter(compound == 'Asp+cetrizine')),(p6_SFa %>% filter(compound == 'Aspartame')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
nrow(asp_cet)
head(asp_cet)
colnames(asp_cet)<-c("Bug_ID","Replicate_no","Phyla","Sp_short","Plate_no","comb_cmp","well",
                     "SFaq","cmp_a","SFa")

asp_cet<-merge(asp_cet,(p6_SFq %>% filter(compound == 'Cetrizine')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
head(asp_cet)


# Asp+Montelukast
asp_mont<-merge((p6_SFaq %>% filter(compound == 'Asp+Montelukast')),(p6_SFa %>% filter(compound == 'Aspartame')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
nrow(asp_mont)
head(asp_mont)
colnames(asp_mont)<-c("Bug_ID","Replicate_no","Phyla","Sp_short","Plate_no","comb_cmp","well",
                      "SFaq","cmp_a","SFa")

asp_mont<-merge(asp_mont,(p6_SFq %>% filter(compound == 'Montelukast')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
head(asp_mont)


# Asp+Risperidone
asp_risp<-merge((p6_SFaq %>% filter(compound == 'Asp+Risperidone')),(p6_SFa %>% filter(compound == 'Aspartame')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
nrow(asp_risp)
head(asp_risp)
colnames(asp_risp)<-c("Bug_ID","Replicate_no","Phyla","Sp_short","Plate_no","comb_cmp","well",
                      "SFaq","cmp_a","SFa")

asp_risp<-merge(asp_risp,(p6_SFq %>% filter(compound == 'Risperidone')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
head(asp_risp)


# Sucralose+EE
suc_ee<-merge((p6_SFaq %>% filter(compound == 'Sucralose+EE')),(p6_SFa %>% filter(compound == 'Sucralose')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
nrow(suc_ee)
head(suc_ee)
colnames(suc_ee)<-c("Bug_ID","Replicate_no","Phyla","Sp_short","Plate_no","comb_cmp","well",
                    "SFaq","cmp_a","SFa")

suc_ee<-merge(suc_ee,(p6_SFq %>% filter(compound == 'Ethinylestradiol')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
head(suc_ee)

# Sucralose+Ibuprofen
suc_ibu<-merge((p6_SFaq %>% filter(compound == 'Sucralose+Ibuprofen')),(p6_SFa %>% filter(compound == 'Sucralose')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
nrow(suc_ibu)
head(suc_ibu)
colnames(suc_ibu)<-c("Bug_ID","Replicate_no","Phyla","Sp_short","Plate_no","comb_cmp","well",
                     "SFaq","cmp_a","SFa")

suc_ibu<-merge(suc_ibu,(p6_SFq %>% filter(compound == 'Ibuprofen')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
head(suc_ibu)

# Sucralose+Montelukast
suc_mont<-merge((p6_SFaq %>% filter(compound == 'Sucralose+Montelukast')),(p6_SFa %>% filter(compound == 'Sucralose')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
nrow(suc_mont)
head(suc_mont)
colnames(suc_mont)<-c("Bug_ID","Replicate_no","Phyla","Sp_short","Plate_no","comb_cmp","well",
                      "SFaq","cmp_a","SFa")

suc_mont<-merge(suc_mont,(p6_SFq %>% filter(compound == 'Montelukast')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
head(suc_mont)

# Sucralose+Ranitidine
suc_rani<-merge((p6_SFaq %>% filter(compound == 'Sucralose+Ranitidine')),(p6_SFa %>% filter(compound == 'Sucralose')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
nrow(suc_rani)
head(suc_rani)
colnames(suc_rani)<-c("Bug_ID","Replicate_no","Phyla","Sp_short","Plate_no","comb_cmp","well",
                      "SFaq","cmp_a","SFa")

suc_rani<-merge(suc_rani,(p6_SFq %>% filter(compound == 'Ranitidine')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
head(suc_rani)

p6_merged<-rbind(acek_acet,acek_cet,acek_risp,acek_ibu,asp_ari,asp_cet,asp_mont,asp_risp,suc_ee,suc_ibu,suc_mont,suc_rani)
nrow(p6_merged) #864

View(p6_merged)

# bliss

p6_bliss<-p6_merged %>% filter_if(~is.numeric(.), all_vars(!is.infinite(.)))
nrow(p6_bliss) #858

p6_bliss<-p6_bliss %>% mutate(bliss_score = SFa+SFq-SFaq, bliss_indp = exp(SFa+SFq-SFaq))
nrow(p6_bliss) #858
head(p6_bliss)

comm_bliss<-p6_bliss %>% dplyr::group_by(Bug_ID,Plate_no,comb_cmp,Phyla,Sp_short) %>% 
  do(anov(.))
nrow(comm_bliss) #300
head(comm_bliss)

AQ<-merge(comm_bliss,p6_bliss,by=c('Bug_ID','Plate_no','comb_cmp','Phyla','Sp_short'))
AQ<-AQ %>% dplyr::group_by(Bug_ID,Plate_no,comb_cmp) %>% mutate(comp_count = length(comb_cmp))
head(AQ)
nrow(AQ) #573
View(AQ)

# expected vs observed viability
CairoSVG(file="../Figures/bliss_obs_exp.svg", width = 9, height = 6, bg = "white")
AQ %>% dplyr::group_by(Bug_ID,Sp_short,comb_cmp,Replicate_no,Plate_no) %>% 
  ggplot(aes(x=SFaq,y=SFa+SFq)) + geom_point(aes(color=comb_cmp), size=0.0001) + geom_abline() + theme_bw() +
  scale_y_continuous(name = "Expected viability") + 
  scale_x_continuous(name = "Observed viability") + scale_color_d3(palette = 'category20') + 
  facet_wrap(~Sp_short, scales = "free") + th
dev.off()

sig_bliss<-AQ %>% filter(pv < 0.05)
nrow(sig_bliss)
head(sig_bliss)


#take median of bliss scores
p6_med_sig<-sig_bliss %>% dplyr::group_by(Bug_ID,Plate_no,comb_cmp,Phyla,Sp_short) %>% summarise(med_SFaq = median(SFaq), med_SFa = median(SFa), med_SFq = median(SFq),
                                                                                              med_bliss = median(bliss_score), med_blissInd = median(bliss_indp),
                                                                                              pv = median(pv), comp_count = median(comp_count))

head(p6_med_sig)
nrow(p6_med_sig) #164


p6_med_sig %>% ggplot(aes(x=med_bliss, color = comb_cmp)) + scale_color_d3(palette = 'category20') + th +
  stat_ecdf() + scale_x_continuous(breaks = seq(-3,3,by=0.2), name = "bliss independence") + scale_y_continuous(name = "% counts") +
  theme_minimal() + geom_vline(xintercept = c(0.2,-.2), lwd=0.2, color='red') 



CairoSVG(file="../Figures/p6_sig_bliss_interactions.svg", width = 4, height = 3, bg = "white")
p6_med_sig %>% filter(med_bliss > 0.10 | med_bliss < -0.10) %>% ggplot(aes(x=Sp_short,y=comb_cmp)) + geom_tile(aes(fill=med_bliss)) + theme_bw() +
  th + theme(axis.text.x = element_text(angle = 90,hjust = 1), axis.ticks = element_blank(), panel.grid = element_blank()) + 
  scale_fill_gradient2() + scale_y_discrete(name = "compounds") 
dev.off()

#case example - synergy
exam<-AQ %>% filter(Bug_ID == 'NT5037', comb_cmp == 'Ace-K+Ibuprofen') %>% mutate(expSFaq = exp(SFaq), expSFa = exp(SFa), expSFq = exp(SFq), expected = exp(SFa+SFq))
View(exam)

CairoSVG(file="../Figures/Isosteviol+D_synergy.svg", width = 4, height = 3, bg = "white")
exam %>% reshape2::melt() %>% dplyr::group_by(Replicate_no) %>% filter(variable %in% c('expSFaq','expSFa','expSFq','expected')) %>% ggplot(aes(x=variable,y=value)) + 
  geom_boxplot(lwd=0.1,aes(fill=variable),alpha=0.3) + scale_x_discrete(labels=c("expSFa" = "sweetener", "expSFq" = "drug",
                                                                                 "expSFaq" = "observed", "expected" = "expected")) + th + scale_fill_npg() +
  theme(axis.title.x = element_blank()) + geom_jitter(aes(color=Replicate_no)) + scale_color_nejm()
  scale_y_continuous(name = 'Surviving fraction')
dev.off()


# sweetener alone


