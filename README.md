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

 “Leo Lahti, Sudarshan Shetty et al. (2017). Tools for microbiome analysis in R. Version 1.5.28. URL: http://microbiome.github.com/microbiome. 

dplyr (Wickham, Francois, Henry, et al., 2017), ggplot2 (Wickham, 2009), phyloseq (McMurdie and Holmes, 2013), tidyr (Wickham, 2017), vegan (Oksanen, Blanchet, Friendly, et al., 2017)



