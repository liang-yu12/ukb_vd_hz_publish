# Analytic codes for "The Association between Vitamin D and Incident Herpes Zoster: A UK Biobank Study"

These are the codes for plotting the forest plot in the manuscript.

The hazard ratio of the Weibull regression model was obtained by using `survival` package (ver 3.2-11), and the regression output was organized using `tidyr()` function from the `Broom` package (ver 0.7.9). A forest plot organized these outputs was plotted using `ggplot2` package (ver 3.3.5).

The results of sensitivity analysis were summarized in Supplementary Figures.

| **Sensitivity analysis**                                     | Figure                                               |
| ------------------------------------------------------------ | ---------------------------------------------------- |
| Stop follow-up by 31 August 2013                             | Supplementary Figures 2<br />Supplementary Figures 3 |
| Using different covariates definitions                       | Supplementary Figures 4                              |
| Use Cox-regression to analyse the association  between exposure and outcomes | Supplementary Figures 5<br />Supplementary Figures 6 |
| Use Poisson regression to analyse the  association between exposure and outcomes | Supplementary Figures 7<br />Supplementary Figures 8 |



## Files Index:

### *Main analysis:*

-   [Figure 3: The association between vitamin D status and the risk of herpes zoster.](https://github.com/liang-yu12/ukb_vd_hz_publish/blob/main/Figure_3_vd_hz_weibull.R)

-   [Figure 4: The association between vitamin D intake and the risk of herpes zoster.](https://github.com/liang-yu12/ukb_vd_hz_publish/blob/main/Figure_4_supdrug_hz_weibull.R)

### *Supplementary Figures:*

-   [Supplementary Figure 2. Vitamin D status and the risk of herpes zoster excluding records after September 2013.](https://github.com/liang-yu12/ukb_vd_hz_publish/blob/main/Supplementary_Fig2_novaccine_vdhz_weibull.R)

-   [Supplementary Figure 3. Vitamin D intake and the risk of herpes zoster excluding records after September 2013.](https://github.com/liang-yu12/ukb_vd_hz_publish/blob/main/Supplementary_Fig3_novaccine_supdrug_hz_weibull.R)

-   [Supplementary Figure 4. Sensitivity analysis of using different definitions of clinical covariates.](https://github.com/liang-yu12/ukb_vd_hz_publish/blob/main/Supplementary_Fig4_defcov_vd_hz_weibull.R)

-   [Supplementary Figure 5. Sensitivity analysis of using stratified Cox regression to assess the association between vitamin D status and the hazards of incident herpes zoster before and after the vaccination program initiated.](https://github.com/liang-yu12/ukb_vd_hz_publish/blob/main/Supplementary_Fig5_vd_hz_cox.R)

-   [Supplementary Figure 6 Sensitivity analysis of using Cox proportional-hazards model to examine the association between vitaminD intake and the risk of herpes zoster.](https://github.com/liang-yu12/ukb_vd_hz_publish/blob/main/Supplementary_Fig6_supdrug_hz_cox.R)

-   [Supplementary Figure 7. Sensitivity analysis of using Poission regression model to examine the association between vitamin D status and the risk of herpes zoster](https://github.com/liang-yu12/ukb_vd_hz_publish/blob/main/Supplementary_Fig7_vd_hz_poisson.R)

-   [Supplementary Figure 8. Sensitivity analysis of using Poission regression model to examine the association between vitamin D supplementation and the risk of herpes zoster](https://github.com/liang-yu12/ukb_vd_hz_publish/blob/main/Supplementary_Fig8_supdrug_hz_poisson.R)

### *Code-lists*:
All diagnostic codes are in CTV3, SNOMED-CT, or read 2. Prescription codes are in BNF, DM+D and Read 2 codes

-   [Comorbidities](https://github.com/liang-yu12/ukb_vd_hz_publish/tree/main/code_lists/covariates_comorbidities): including asthma, COPD, Chronic kidney disease, rheumatoid arthritis, inflammatory bowel diseases, depression, 

-   [Immunosuppression](https://github.com/liang-yu12/ukb_vd_hz_publish/tree/main/code_lists/covariates_immunosuppression): including including organ transplantation, chemoradiotherapy, cell-mediated immunosuppression, HIV, blood cancers, chemotherapy (biological and non-biological agents), bone marrow transplantation and long-term oral steroid. 

-   [Vitamin D prescriptions](https://github.com/liang-yu12/ukb_vd_hz_publish/tree/main/code_lists/exposure_vitd_drug): vitamin D prescriptions codes

-   [Herpes zoster diagnosis](https://github.com/liang-yu12/ukb_vd_hz_publish/tree/main/code_lists/outcome_hz): diagnostic codes for diagnosing herpes zoster cases
