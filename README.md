# Analytic codes for "The Association between Vitamin D and Incident Herpes Zoster: A UK Biobank Study"

These are the codes for plotting the forest plot in the manuscript.

Rate ratio was obtained using the Poisson regression model, and the regression output was organized using `tidyr()` function from the Broom package. Then A forest plot of these organized outputs was plotted using `ggplot2` package.

In the sensitivity analysis reported in *Supplementary Fig4* and *Supplementary Fig5*, Cox regression model was used. In *Supplementary Figure 4*, because of the violation of proportional hazard assumption, a stratified Cox model was used. In *Supplementary Figure 5*, the proportional hazards assumption was not violated, so a Cox proportional hazard regression was used.

## Files Index:

*Main analysis:*

-   [Figure 3: The association between vitamin D status and the risk of herpes zoster.](https://github.com/liang-yu12/ukb_vd_hz_publish/blob/aca35f24a6c6b4015d3e69929bf3f79ffe4a48bb/figure_2_vd_hz.R)

-   [Figure 4: The association between vitamin D intake and the risk of herpes zoster.](https://github.com/liang-yu12/ukb_vd_hz_publish/blob/aca35f24a6c6b4015d3e69929bf3f79ffe4a48bb/figure_3_vdsupdrug_hz.R)

*Supplementary Figures:*

-   [Supplementary Figure 1. The association between vitamin D status and the risk of herpes zoster excluding records after September 2013.](https://github.com/liang-yu12/ukb_vd_hz_publish/blob/aca35f24a6c6b4015d3e69929bf3f79ffe4a48bb/supp_fig1_se_novaccine.R)

-   [Supplementary Figure 2. Vitamin D intake and the risk of herpes zoster excluding records after September 2013.](https://github.com/liang-yu12/ukb_vd_hz_publish/blob/aca35f24a6c6b4015d3e69929bf3f79ffe4a48bb/supp_fig2_se_novaccine.R)

-   [Supplementary Figure 3. Sensitivity analysis of using different definitions of clinical covariates.](https://github.com/liang-yu12/ukb_vd_hz_publish/blob/aca35f24a6c6b4015d3e69929bf3f79ffe4a48bb/supp_fig3_se_differentcovariates.R)

-   [Supplementary Figure 4. Sensitivity analysis of using stratified Cox regression to assess the association between vitamin D status and the hazards of incident herpes zoster before and after the vaccination program initiated.](https://github.com/liang-yu12/ukb_vd_hz_publish/blob/aca35f24a6c6b4015d3e69929bf3f79ffe4a48bb/supp_fig4_primary_surv.R)

-   [Supplementary Figure 5 Sensitivity analysis of using Cox proportional-hazards model to examine the association between vitaminD intake and the risk of herpes zoster.](https://github.com/liang-yu12/ukb_vd_hz_publish/blob/aca35f24a6c6b4015d3e69929bf3f79ffe4a48bb/supp_fig5_2ndexp_surv.R)
