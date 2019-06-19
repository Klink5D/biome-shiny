# Biome-shiny
##  A Shiny R app for microbiome analysis, built around the "microbiome" package

Biome-shiny is an application that provides a simple, easy-to-use interface to visualize and analyze microbiomes through (mostly) interactive plots, and export the results to a file.
It was built using the Shiny library, and draws most of its features from the "microbiome" library.

Current features include:
	Viewing core microbiota
	Community composition	
	Alpha and Beta Diversity
	Community landscape
	DMM clustering
	PERMANOVA
	

The app accepts a phyloseq-style .biom file with metadata as an input, but can also use the sample datasets from the microbiome package for testing purposes.

Tutorial coming soon



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

  Martin Morgan (2019). DirichletMultinomial: Dirichlet-Multinomial Mixture Model Machine
  Learning for Microbiome Data. R package version 1.24.1.

