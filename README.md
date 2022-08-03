# FCStabilityOptimality_share

============================================/
Analysis code to recreate data for Functional connectome stability and optimality are predictors of cognitive performance

FC_Stab_Typ_Opt_Disc.m calculates individual and group-level connectome features for Datasets 1-3 and their relationship with behavioral performance on sustained attention and working memory tasks.

Network_StabTypOptDisc.m calculates individual and group-level connectome features within 10 canonical networks 

LesionNetworks_StabTyp.m calculates connectome features after computationally lesioning 10 canonical networks to quantify network contribution to prediction

LesionNetworks_NullDistribution_StabTyp.m computationlly lesions random nodes to create distribution of possible network contribution to prediction in order to determine significance of network contribution to prediction

CPM_model_comparison.m calculates sustained attention and working memory CPM network strength in each participant. Outputs a .txt file to be read by FC_StabTypOptDisc_Models.R to construct linear models

FC_StabTypOptDisc_Models.R runs model comparison for mixed-effects models with connectome feature predictors.


==============================================\
SOFTWARE VERSION INFO:

MATLAB 
version '9.11.0.1837725 (R2021b) Update 2'
date 'December 14, 2021'

R version 4.0.2 (2020-06-22)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Catalina 10.15.7

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib

Random number generation:
 RNG:     Mersenne-Twister 
 Normal:  Inversion 
 Sample:  Rounding 
 
locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] lmerTest_3.1-3    lme4_1.1-26       Matrix_1.2-18     data.table_1.14.0 ggplot2_3.3.2     apastats_0.3     

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.8.3        compiler_4.0.2      pillar_1.7.0        nloptr_1.2.2.2      tools_4.0.2        
 [6] boot_1.3-25         digest_0.6.25       statmod_1.4.35      evaluate_0.14       lifecycle_1.0.0    
[11] tibble_3.1.7        gtable_0.3.0        nlme_3.1-148        lattice_0.20-41     pkgconfig_2.0.3    
[16] rlang_1.0.2         rstudioapi_0.13     cli_3.2.0           yaml_2.2.1          xfun_0.30          
[21] fastmap_1.1.0       withr_2.4.2         dplyr_1.0.2         knitr_1.30          generics_0.0.2     
[26] vctrs_0.3.8         tidyselect_1.1.0    grid_4.0.2          glue_1.6.2          R6_2.4.1           
[31] fansi_0.4.1         rmarkdown_2.13      minqa_1.2.4         purrr_0.3.4         magrittr_1.5       
[36] scales_1.1.1        htmltools_0.5.2     ellipsis_0.3.2      MASS_7.3-51.6       splines_4.0.2      
[41] colorspace_1.4-1    numDeriv_2016.8-1.1 utf8_1.1.4          munsell_0.5.0       crayon_1.3.4    
![image](https://user-images.githubusercontent.com/73361771/182484176-dfde5336-14cb-4339-9662-c5dcf85a7587.png)




