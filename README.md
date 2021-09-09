# Analytic codes for "The Association between Vitamin D and Incident Herpes Zoster: A UK Biobank Study"

These are the codes for plotting the forest plot in the manuscript.

Rate ratio was obtained using the Poisson regression model, and the regression output was organized using `tidyr()` function from the Broom package. Then A forest plot of these organized outputs was plotted using `ggplot2` package.

In the sensitivity analysis reported in *Supplementary Fig4* and *Supplementary Fig5*, Cox regression model was used. In *Supplementary Figure 4*, because of the violation of proportional hazard assumption, a stratified Cox model was used. In *Supplementary Figure 5*, the proportional hazards assumption was not violated, so a Cox proportional hazard regression was used.
