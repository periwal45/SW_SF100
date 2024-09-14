setwd("/Users/periwal/ShikiFactory/WP3/SW_SF100")

params<-read.table(file = 'Figures/cleaned_params', header = TRUE, sep = '\t', stringsAsFactors = FALSE)
params<-params %>% filter(Plate_no == 'plate6')
head(params)
dim(params) #6,867 X 9

# source the script with functions
source("/Users/periwal/ShikiFactory/WP3/SW_SF100/p6_Functions.R")

cnormAUCl<-params %>% filter(compound == 'DMSO') %>% dplyr::group_by(Bug_ID,Plate_no,Replicate_no,Phyla,Sp_short) %>%
  summarise(ctrl_med = round(median(auc_l),digits = 2))
head(cnormAUCl)
dim(cnormAUCl)

normAUCl<-merge(params,cnormAUCl,by=c("Bug_ID","Replicate_no","Plate_no","Phyla","Sp_short"))
dim(normAUCl)
head(normAUCl)

cnormAUCs<-normAUCl %>% dplyr::group_by(Bug_ID,Plate_no,Replicate_no,col_name) %>%
  mutate(cAUC = round(auc_l/ctrl_med, digits = 2), FC = cAUC) #FC treatment vs control

head(cnormAUCs)

cnormAUCs %>%  dplyr::group_by(Bug_ID,Plate_no,Replicate_no) %>% 
  ggplot(aes(cAUC)) + geom_density(aes(fill=Replicate_no,color=Replicate_no), alpha=0.3, lwd=0.3) +
  facet_wrap("Sp_short", scales = "free") +
  scale_fill_jama() + scale_color_jama() + th + theme_bw()

View(cnormAUCs)


# t-test bug wise and replicate wise

hits_BR<-cnormAUCs %>% filter(compound %in% c('Ace-K','Aspartame','Sucralose','Acetaminophen','Aripirazole','Cetrizine','Ethinylestradiol','Ibuprofen','Montelukast','Ranitidine','Risperidone','DMSO')) %>%
  dplyr::group_by(Bug_ID, Plate_no) %>% 
  do(p6_unpairedT_pval(.))
View(hits_BR)
nrow(hits_BR) #300

hits_BR_noDMSO<-hits_BR %>% filter(compound != 'DMSO')

CairoSVG(file=paste("Figures/p6_pval_dist.svg", sep = ""), width = 3, height = 2, bg = "white")
hits_BR_noDMSO %>% ggplot(aes(x=pv)) + geom_histogram(color="white",fill="#4d4d4d",bins = 30,lwd=0.1) + th
dev.off()

log2FC<-hits_BR %>% dplyr::group_by(Bug_ID,Plate_no,compound,Phyla,Sp_short) %>% 
  mutate(l2FC = round(log2(median(avFC)),digits = 2)) #log 2 of median FC of replicates
head(log2FC)
dim(log2FC)

#nrow(log2FC %>% filter(pv < 0.05)) #136

### multiple hypotheses testing: error correction
p_bh<-log2FC %>% dplyr::group_by(Bug_ID) %>%
  mutate(pv_bh = p.adjust(pv,method = "BH"))
head(p_bh)
nrow(p_bh)

hits<-p_bh %>% filter(l2FC != '-Inf')
head(hits)
nrow(hits)


write.table(hits, file = 'Figures/plate6_aloneComp', sep = '\t', row.names = FALSE)


nrow(hits %>% filter(pv_bh < 0.05)) #119

CairoSVG(file=paste("Figures/p6_volcano_reps.svg", sep = ""), width = 10, height = 6, bg = "white")
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

CairoSVG(file=paste("Figures/p6_heatmap_reps.svg", sep = ""), width = 10, height = 4, bg = "white")
hits %>% filter(pv_bh < 0.05 & (l2FC >= 0.32 | l2FC <= -0.32)) %>% ggplot(aes(x=Sp_short,y=compound)) + geom_tile(aes(fill=l2FC), color = "white", lwd = 0.1) + theme_bw() +
  th + theme(axis.text.x = element_text(angle = 90,hjust = 1), axis.ticks = element_blank(), panel.grid = element_blank()) +
  scale_fill_gradientn(colors = wes_palette("BottleRocket1", 100, type = "continuous"), name = bquote(log[2]~'fold change'))
dev.off()


# plate6
#normAUCs food compounds
p6<-cnormAUCs
dim(p6)  
head(p6)

p6_SFa<-p6 %>% filter(compound %in% c('Ace-K','Aspartame','Sucralose')) %>% dplyr::group_by(Bug_ID,Replicate_no,Plate_no,compound,Phyla,Sp_short) %>% 
  summarise(SFa = log(median(cAUC)))
head(p6_SFa)  
nrow(p6_SFa) #216

#normAUCs drugs
p6_SFq<-p6 %>% filter(compound %in% c('Acetaminophen','Aripirazole','Cetrizine','Ethinylestradiol','Ibuprofen','Montelukast','Ranitidine','Risperidone')) %>%
  dplyr::group_by(Bug_ID,Replicate_no,Plate_no,compound,Phyla,Sp_short) %>% 
  summarise(SFq = log(median(cAUC)))
head(p6_SFq)
nrow(p6_SFq) #576

sd<-c('Ace-K','Aspartame','Sucralose','Acetaminophen','Aripirazole','Cetrizine','Ethinylestradiol','Ibuprofen','Montelukast','Ranitidine','Risperidone')
############################################################
#normAUCs food+drug
p6_SFaq<-p6 %>% filter(!compound %in% sd) %>% dplyr::group_by(Bug_ID,Replicate_no,Plate_no,compound,Phyla,Sp_short) %>%
  summarise(SFaq = log(median(cAUC)))
head(p6_SFaq)
nrow(p6_SFaq) #936
View(p6_SFaq)

# combinations
#AceK+Acetaminophen
acek_acet<-merge((p6_SFaq %>% filter(compound == 'Ace-K+acetaminophen')),(p6_SFa %>% filter(compound == 'Ace-K')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
nrow(acek_acet)
head(acek_acet)
colnames(acek_acet)<-c("Bug_ID","Replicate_no","Phyla","Sp_short","Plate_no","comb_cmp",
                       "SFaq","cmp_a","SFa")

acek_acet<-merge(acek_acet,(p6_SFq %>% filter(compound == 'Acetaminophen')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
head(acek_acet)


# AceK+Cetrizine
acek_cet<-merge((p6_SFaq %>% filter(compound == 'Ace-K+cetrizine')),(p6_SFa %>% filter(compound == 'Ace-K')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
nrow(acek_cet)
head(acek_cet)
colnames(acek_cet)<-c("Bug_ID","Replicate_no","Phyla","Sp_short","Plate_no","comb_cmp",
                      "SFaq","cmp_a","SFa")

acek_cet<-merge(acek_cet,(p6_SFq %>% filter(compound == 'Cetrizine')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
head(acek_cet)

# AceK+Risperidone
acek_risp<-merge((p6_SFaq %>% filter(compound == 'Ace-K+Risperidone')),(p6_SFa %>% filter(compound == 'Ace-K')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
nrow(acek_risp)
head(acek_risp)
colnames(acek_risp)<-c("Bug_ID","Replicate_no","Phyla","Sp_short","Plate_no","comb_cmp",
                       "SFaq","cmp_a","SFa")

acek_risp<-merge(acek_risp,(p6_SFq %>% filter(compound == 'Risperidone')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
head(acek_risp)


# AceK+Ibuprofen
acek_ibu<-merge((p6_SFaq %>% filter(compound == 'Ace-K+Ibuprofen')),(p6_SFa %>% filter(compound == 'Ace-K')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
nrow(acek_ibu)
head(acek_ibu)
colnames(acek_ibu)<-c("Bug_ID","Replicate_no","Phyla","Sp_short","Plate_no","comb_cmp",
                      "SFaq","cmp_a","SFa")

acek_ibu<-merge(acek_ibu,(p6_SFq %>% filter(compound == 'Ibuprofen')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
head(acek_ibu)


# Asp+Aripiprazole
asp_ari<-merge((p6_SFaq %>% filter(compound == 'Asp+Aripiprazole')),(p6_SFa %>% filter(compound == 'Aspartame')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
nrow(asp_ari)
head(asp_ari)
colnames(asp_ari)<-c("Bug_ID","Replicate_no","Phyla","Sp_short","Plate_no","comb_cmp",
                     "SFaq","cmp_a","SFa")

asp_ari<-merge(asp_ari,(p6_SFq %>% filter(compound == 'Aripirazole')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
head(asp_ari)

# Asp+Cetrizine
asp_cet<-merge((p6_SFaq %>% filter(compound == 'Asp+cetrizine')),(p6_SFa %>% filter(compound == 'Aspartame')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
nrow(asp_cet)
head(asp_cet)
colnames(asp_cet)<-c("Bug_ID","Replicate_no","Phyla","Sp_short","Plate_no","comb_cmp",
                     "SFaq","cmp_a","SFa")

asp_cet<-merge(asp_cet,(p6_SFq %>% filter(compound == 'Cetrizine')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
head(asp_cet)


# Asp+Montelukast
asp_mont<-merge((p6_SFaq %>% filter(compound == 'Asp+Montelukast')),(p6_SFa %>% filter(compound == 'Aspartame')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
nrow(asp_mont)
head(asp_mont)
colnames(asp_mont)<-c("Bug_ID","Replicate_no","Phyla","Sp_short","Plate_no","comb_cmp",
                      "SFaq","cmp_a","SFa")

asp_mont<-merge(asp_mont,(p6_SFq %>% filter(compound == 'Montelukast')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
head(asp_mont)


# Asp+Risperidone
asp_risp<-merge((p6_SFaq %>% filter(compound == 'Asp+Risperidone')),(p6_SFa %>% filter(compound == 'Aspartame')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
nrow(asp_risp)
head(asp_risp)
colnames(asp_risp)<-c("Bug_ID","Replicate_no","Phyla","Sp_short","Plate_no","comb_cmp",
                      "SFaq","cmp_a","SFa")

asp_risp<-merge(asp_risp,(p6_SFq %>% filter(compound == 'Risperidone')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
head(asp_risp)


# Sucralose+EE
suc_ee<-merge((p6_SFaq %>% filter(compound == 'Sucralose+EE')),(p6_SFa %>% filter(compound == 'Sucralose')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
nrow(suc_ee)
head(suc_ee)
colnames(suc_ee)<-c("Bug_ID","Replicate_no","Phyla","Sp_short","Plate_no","comb_cmp",
                    "SFaq","cmp_a","SFa")

suc_ee<-merge(suc_ee,(p6_SFq %>% filter(compound == 'Ethinylestradiol')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
head(suc_ee)

# Sucralose+Ibuprofen
suc_ibu<-merge((p6_SFaq %>% filter(compound == 'Sucralose+Ibuprofen')),(p6_SFa %>% filter(compound == 'Sucralose')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
nrow(suc_ibu)
head(suc_ibu)
colnames(suc_ibu)<-c("Bug_ID","Replicate_no","Phyla","Sp_short","Plate_no","comb_cmp",
                     "SFaq","cmp_a","SFa")

suc_ibu<-merge(suc_ibu,(p6_SFq %>% filter(compound == 'Ibuprofen')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
head(suc_ibu)

# Sucralose+Montelukast
suc_mont<-merge((p6_SFaq %>% filter(compound == 'Sucralose+Montelukast')),(p6_SFa %>% filter(compound == 'Sucralose')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
nrow(suc_mont)
head(suc_mont)
colnames(suc_mont)<-c("Bug_ID","Replicate_no","Phyla","Sp_short","Plate_no","comb_cmp",
                      "SFaq","cmp_a","SFa")

suc_mont<-merge(suc_mont,(p6_SFq %>% filter(compound == 'Montelukast')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
head(suc_mont)

# Sucralose+Ranitidine
suc_rani<-merge((p6_SFaq %>% filter(compound == 'Sucralose+Ranitidine')),(p6_SFa %>% filter(compound == 'Sucralose')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
nrow(suc_rani)
head(suc_rani)
colnames(suc_rani)<-c("Bug_ID","Replicate_no","Phyla","Sp_short","Plate_no","comb_cmp",
                      "SFaq","cmp_a","SFa")

suc_rani<-merge(suc_rani,(p6_SFq %>% filter(compound == 'Ranitidine')),by=c('Bug_ID','Replicate_no','Phyla','Sp_short','Plate_no'))
head(suc_rani)

p6_merged<-rbind(acek_acet,acek_cet,acek_risp,acek_ibu,asp_ari,asp_cet,asp_mont,asp_risp,suc_ee,suc_ibu,suc_mont,suc_rani)
nrow(p6_merged) #864

View(p6_merged)

# bliss
p6_bliss<-p6_merged %>% mutate(null_hyp = SFa+SFq-SFaq, exp_eff = SFa+SFq)
nrow(p6_bliss) #864
head(p6_bliss)

write.table(p6_bliss, file = "Figures/p6_bliss_scores", sep = "\t", quote = FALSE, row.names = FALSE)

p6_bliss %>% ggplot(aes(null_hyp)) + geom_density(aes(fill=comb_cmp, color=comb_cmp), alpha=0.5) +
  facet_wrap("Sp_short", nrow = 3, scales = "free") + theme_bw() +
  scale_fill_uchicago() + scale_color_uchicago() + th

############statistical determination of synergy/antagonism

p6_A<-p6_bliss %>% dplyr::group_by(Bug_ID,Plate_no,compound,Phyla,Sp_short,comb_cmp) %>% 
  do(t_distr(.))
View(p6_A)
nrow(p6_A) #300

write.table(p6_A, file = "Figures/p6_bliss_scores_sigTesting", sep = "\t", quote = FALSE, row.names = FALSE)

sig_bliss<-p6_A %>% filter(pval < 0.05)
nrow(sig_bliss) #6
head(sig_bliss)

# CairoSVG(file="../Figures/sig_bliss_cutoffs.svg", width = 5, height = 3, bg = "white")
# sig_bliss %>% ggplot(aes(x=bliss_null, color = Plate_no)) + scale_color_npg() + th +
#   stat_ecdf(geom = 'point', size = 0.001) + scale_x_continuous(breaks = seq(-3,3,by=0.2), name = "bliss independence") + scale_y_continuous(name = "% counts") +
#   theme_minimal() + geom_vline(xintercept = c(0.8,1.2), lwd=0.2, color='red') 
# dev.off()

pal<-paletteer_c("ggthemes::Red-Blue Diverging", 100)
pal

CairoSVG(file="Figures/p6_sig_bliss_interactions.svg", width = 7, height = 4, bg = "white")
sig_bliss %>% filter(pval < 0.05) %>% ggplot(aes(x=Sp_short,y=comb_cmp)) + geom_tile(aes(fill=bliss_null)) + theme_bw() +
  th + theme(axis.text.x = element_text(angle = 90,hjust = 1), axis.ticks = element_blank(), panel.grid = element_blank()) + 
  scale_fill_gradientn(colors = pal) + scale_y_discrete(name = "compounds") 
dev.off()

View(sig_bliss %>% filter(bliss_null > 0.1 | bliss_null < -0.1)) #50



