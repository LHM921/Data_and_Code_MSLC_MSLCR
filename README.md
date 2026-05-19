***
Source R codes and data for the Mathematics and Computers in Simulation paper: 
"Multivariate Slash Censored Regression Approach with Multiple Detection Limits",
by Wan-Lun Wang, Hui-Min Li and Wei-Heng Huang
***

# Author responsible for the code #
For questions, comments or remarks about the code please contact responsible author, Wan-Lun Wang (wangwl@gs.ncku.edu.tw), Hui-Min Li (fa879727@gmail.com) and  Wei-Heng Huang (weihenghuang@gm.ntpu.edu.tw). 

# Configurations #
The code was written/evaluated in R with the following software versions:
R version 4.4.2 (2024-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26100)

Matrix products: default

locale:
[1] LC_COLLATE=Chinese (Traditional)_Taiwan.utf8 
[2] LC_CTYPE=Chinese (Traditional)_Taiwan.utf8   
[3] LC_MONETARY=Chinese (Traditional)_Taiwan.utf8
[4] LC_NUMERIC=C                                 
[5] LC_TIME=Chinese (Traditional)_Taiwan.utf8    

time zone: Asia/Taipei
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

loaded via a namespace (and not attached):
[1] compiler_4.4.2

# Descriptions of the codes #
Please extract the file **Data and Code.zip** to the current working directory of the R package.
The getwd() function shall determine an absolute pathname of the "current working directory".

Before running the codes **Table1.r**, **Table2.r**, **Table3.r**, **Table4.r**, **Table5.r**, **TableD1.r**, **TableD2.r**, **TableD3.r**, **TableE1.r**, **TableE2.r**, **Fig1.r**, **Fig2.r**, **Fig3.r**, **Fig4.r**, **Fig5.r**, **Fig6.r**, **Fig7.r** and **Fig8.r**, one needs to install the following R packages 

    install.packages("mvtnorm")  Version: 1.3-3
    install.packages("tmvtnorm") Version: 1.6
    install.packages("cubature") Version: 2.1.1
    install.packages("mclust")  Version: 6.1.1
    install.packages("PearsonDS") Version: 1.3.1
    install.packages("invgamma") Version: 1.1
    install.packages("Bessel")  Version: 0.6-1
    install.packages("scatterplot3d") Version: 0.3-44

 
R codes for the implementation of the proposed methodology are provided.

## Subfolder: `./function` ##
`./function`
contains the program (function) of
- (1) **TMSLmoment.r**: main script for evaluting the first two moments of truncated multivariate normal (MVN), multivariate-t (MVT) and multivariate slash (MSL) distributions;
- (2) **MSL.EM.r**: main script for implementing the MCECM algorithm for the MSL models;
- (3) **MSLR.EM.r**: main script for implementing the MCECM algorithm for the MSLR models;
- (4) **MSLC.EM.r**: main script for implementing the MCECM algorithm for the MSLC model;
- (5) **MSLCR.EM.r**: main script for implementing the MCECM algorithm for the MSLCR model;
- (6) **MVT.EM.r**: main script for implementing the MCECM algorithm for the MVT model;
- (7) **MVTR.EM.r**: main script for implementing the MCECM algorithm for the MVTR model;
- (8) **MVTC.EM.r**: main script for implementing the MCECM algorithm for the MVTC model;
- (9) **MVTCR.EM.r**: main script for implementing the MCECM algorithm for the MVTCR model.
- (10) **IMC.r**: main script for implementing the asymptotic standard errors for parameters.
- (11) **IMCR.r**: main script for implementing the asymptotic standard errors for regression coefficients.

## Subfolder: `./code` ##
`./code`
contains main scripts of 
- (1) **Fig1.r**: main script for reproducting Figure 1;
- (2) **Fig2.r**: main script for reproducting Figure 2;
- (3) **Fig3.r**: main script for reproducting Figure 3 (run **Table1.r** first or load **bioassay_sep.RData** and **bioassay_nov.RData*  directly, and then run **Fig3.r**);
- (4) **Fig4.r**: main script for reproducting Figure 4 (run **Table2.r** first or load **lpdata.RData** directly, and then run **Fig4.r**);
- (5) **Fig5.r**: main script for reproducting Figure 5 (run **Table3.r** first or load **lake.RData** directly, and then run **Fig5.r**);
- (6) **Fig6.r**: main script for reproducting Supplementary Figure 1; (read **Smse_all.csv** for 100 Monte Carlo datasets);
- (7) **Fig7.r**: main script for reproducting Supplementary Figure 2; (read **Smse_all.csv** for 100 Monte Carlo datasets);
- (8) **Fig8.r**: main script for reproducting Supplementary Figure 3; (read **Smse_all.csv** for 100 Monte Carlo datasets);
- (9) **Table1.r**: main script for reproducting all fitting results of MVN, MVT, MSL, MVNC, MVTC and MSLC models to the bioassay data;
- (10) **Table2.r**: main script for reproducting all fitting results of MVN, MVT, MSL, MVNC, MVTC and MSLC models to the LPdata;
- (11) **Table3.r**: main script for reproducting all fitting results of MVNR, MVTR, MSLR, MVNCR, MVTCR and MSLCR models to the Lake data.
- (12) **Table4.r**: main script for reproducting all fitting results of MVN, MVT, MSL, MVNC, MVTC and MSLC models;
- (13) **Table5.r**: main script for reproducting all fitting results of MVNR, MVTR, MSLR, MVNCR, MVTCR and MSLCR models;
- (14) **TableD1.r**: main script for reproducting all estimation results of MVN, MVT, MSL, MVNC, MVTC and MSLC models to the bioassay data;
- (15) **TableD2.r**: main script for reproducting all estimation results of MVN, MVT, MSL, MVNC, MVTC and MSLC models to the LPdata;
- (16) **TableD3.r**: main script for reproducting all estimation results of MVNR, MVTR, MSLR, MVNCR, MVTCR and MSLCR models to the Lake data.
- (17) **TableE1.r**: main script for reproducting all fitting results of MVN, MVT, MSL, MVNC, MVTC and MSLC models;
- (18) **TableE2.r**: main script for reproducting all fitting results of MVNR, MVTR, MSLR, MVNCR, MVTCR and MSLCR models;

## Subfolder: `./data` ##
`./data`
contains
- (1) **bioassay_cen.csv**: the photopigment concentrations data for bioassays dataset with censoring information;
- (2) **bioassay_dl.csv**: the photopigment concentrations data for bioassays dataset;
- (3) **lake_cen.csv**: the lake Michigan water chemistry dataset with censoring information;
- (4) **lake_dl.csv**: the lake Michigan water chemistry dataset replaced with detection limits;
- (5) **lake_ori.csv**: the original dataset of lake Michigan water chemistry data;
- (6) **LPdata.csv**: the La Paloma stream dataset.

## Subfolder: `./result` ##
`./result/SIM`
contains 
- (1) **Run_time_all.csv**: the run time of 100 Monte Carlo dataset from MVN, MVT, MSL, MVNC, MVTC and MSLC models;
- (2) **Sest_all.csv**: the estimation of parameters of 100 Monte Carlo dataset from MVN, MVT, MSL, MVNC, MVTC and MSLC models;
- (3) **Sfit_all.csv**: the model fitting results of 100 Monte Carlo dataset from MVN, MVT, MSL, MVNC, MVTC and MSLC models;
- (4) **Smse_all.csv**: the squared of mean squared errors results of 100 Monte Carlo dataset from MVN, MVT, MSL, MVNC, MVTC and MSLC models;
- (5) **SSE_all.csv**: the information-matrix based standard errors results of 100 Monte Carlo dataset from MVN, MVT, MSL, MVNC, MVTC and MSLC models.

`./result/RSIM`
contains 
- (1) **Run_time_all.csv**: the run time of 100 Monte Carlo dataset from MVNR, MVTR, MSLR, MVNCR, MVTCR and MSLCR models;
- (2) **Sest_all.csv**: the estimation of parameters of 100 Monte Carlo dataset from MVNR, MVTR, MSLR, MVNCR, MVTCR and MSLCR models;
- (3) **Sfit_all.csv**: the model fitting results of 100 Monte Carlo dataset from MVNR, MVTR, MSLR, MVNCR, MVTCR and MSLCR models;
- (4) **Smse_all.csv**: the squared of mean squared errors results of 100 Monte Carlo dataset from MVNR, MVTR, MSLR, MVNCR, MVTCR and MSLCR models;
- (5) **SSE_all.csv**: the information-matrix based standard errors results of 100 Monte Carlo dataset from MVNR, MVTR, MSLR, MVNCR, MVTCR and MSLCR models;
- (6) **SMSERC_all.csv**: the fitting responses and recovered censored values of mean squared errors of 100 Monte Carlo dataset from MVNR, MVTR, MSLR, MVNCR, MVTCR and MSLCR models.

`./result`
contains 
- (1) **bioassay_sep.RData**: the analysis results for the bioassays data (Application 1);
- (2) **bioassay_nov.RData**: the analysis results for the bioassays data (Application 1);
- (3) **lpdata.RData**: the analysis results for the La Paloma stream data (Additional Applications);
- (4) **lake.RData**: the analysis results for the lake Michigan water chemistry data (Application 2).
      
## Additional Remark ##
- Note (1): One can directly run each "source(.)" described in **master.r** file in the seperate R session to obtain the results.
- Note (2): Numerical results have been stored in "./result/", and the fitting results of the considered models obtained by running **Table1.r**, **Table2.r** and **Table3.r** have been stored in `./result/bioassay_sep.RData`, `./result/bioassay_nov.RData`, `./result/lpdata.RData` and `./result/lake.RData`.
- Note (3): To draw Figure 3, please load the **bioassay_sep.RData** and **bioassay_nov.RData* file in subfolder "./result/", and then run the **Fig3.r** script in subfolder `./code/`. 
- Note (4): To draw Supplementary Figure 4, please load the **lpdata.RData* file in subfolder "./result/", and then run **Fig4.r** scripts, respectively, in subfolder `./code/`.
- Note (5): To draw Figure 5, please load the **lake.RData* file in subfolder "./result/", and then run the **Fig5.r** script in subfolder `./code/`. 

