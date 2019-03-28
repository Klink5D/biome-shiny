# meta-shiny
##  A Shiny R app for metagenomics analysis, built around the "microbiome" package

In essence, this app is a student project designed to create a semi or fully-automatic, easy to use pipeline, that could potentially be used by people with little knowledge of coding or bioinformatics.

It also serves as a way for me to further my own knowledge on the creation of pipelines, metagenomics analyses, the Shiny library for R and the R language as a whole.

*What it can do right now:*

Alpha Diversity displays a phyloseq summary of the input dataset (right now, one of the three sample sets available from the microbiome package). It displays a table containing the phyloseq file's metadata, as well as a series of diversity measures. Finally, it displays results in a visual form through a violin plot.

Beta Diversity displays a series of ordination plots (two split, one not split) based around user input, and allows the user to set a unique seed for outcome replication, if necessary.

Community Composition takes the user's input to display a variety of abundance plots: abundance by subject, relative abundance data and a taxa prevalence plot.

The app is in very early stages and thus is extremely buggy and unfit for any actual work. I would hope that by the end, it'll be a functional, user-friendly app that allows even those with little bioinformatics knowledge to perform a thorough metagenomics analysis without having to go look through heaps of documentation.

Citations:

Winston Chang, Joe Cheng, JJ Allaire, Yihui Xie and Jonathan McPherson (2018). shiny: Web Application Framework
for R. R package version 1.2.0. https://CRAN.R-project.org/package=shiny

Leo Lahti, Sudarshan Shetty et al. (2017). Tools for microbiome analysis in R. Version 1.5.28. URL: http://microbiome.github.com/microbiome. 

Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York. ISBN 978-3-319-24277-4, http://ggplot2.org.

Hadley Wickham (2007). Reshaping Data with the reshape Package. Journal of Statistical Software, 21(12), 1-20. URL
http://www.jstatsoft.org/v21/i12/.

Hadley Wickham, Romain François, Lionel Henry and Kirill Müller (2018). dplyr: A Grammar of Data Manipulation. R package version package version 0.8.0.1. https://CRAN.R-project.org/package=dplyr


Alboukadel Kassambara (2018). ggpubr: 'ggplot2' Based Publication Ready Plots. R package version 0.2.999.
http://www.sthda.com/english/rpkgs/ggpubr

Yihui Xie (2019). knitr: A General-Purpose Package for Dynamic Report Generation in R. R package version 1.22.

Bob Rudis (2019). hrbrthemes: Additional Themes, Theme Components and Utilities for 'ggplot2'. R package version
0.6.0. https://CRAN.R-project.org/package=hrbrthemes


