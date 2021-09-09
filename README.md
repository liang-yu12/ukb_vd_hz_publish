# Analytic codes for "The Association between Vitamin D and Incident Herpes Zoster: A UK Biobank Study"

These are the codes for plotting the forest plot in the manuscript.

Rate ratio was obtained using the Poisson regression model, and the regression output was organized using `tidyr()` function from the Broom package. Then A forest plot of these organized outputs was plotted using `ggplot2` package. 
