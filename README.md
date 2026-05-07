README\
Manuscript title: Impact of low-calorie sweeteners on gut bacteria is modulated by common xenobiotics
This repository contains code for processing, cleaning, and analysing data from primary screening microbial assays

Data:
The data folder contains raw OD values taken every hour for 24 hours for all microbes for all biological replicates. The folder also contains intermediary files generated such as fitted curves, fitted values, merged replicates and filtered hits. The raw OD files can be used to re-run the scripts mentioned below

Scripts:
1. primary_screen.R - used in Fig. 1\
     -reads raw OD values and fits growth curves on each plate for all replicates
     -plots all fitted growth curves
     -QC: does data cleaning on fitted curves
     -Normalizes: plate based normalization of area under the curve
     -Hit selection: two sample welch t-test
     -log2fold change calculation
     -Bliss interactions: compounds in combination and statistical determination of synergy & antagonism
   
2. Functions.R - used in Fig. 1\
     -common repeated functions called into primary_screen.R
     -function for two sample welch test
     -function to pool effect sizes
     -statistical determination of synergy/antagonism

3. p6.R - used in Appendix Figure S4\
     -additional xenobiotic interactions plate (sweeteners-commonly used drugs)
     -bliss interactions of all combinations of sweeteners and commonly used drugs

4. p6_Functions.R - used in Appendix Figure S4\
     -functions sourced in p6.R
