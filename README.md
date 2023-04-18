# OROV transmission risk models 
The current repository contains key scripts to recreate the crucial steps to develop the models and evaluation framework presented in the manuscript: "Transmission Risk of Oropouche fever across the Americas" published here: 
Specific steps that you will find within the script include: 

convexhull3d: Function developed by Luis Osorio-Olvera and modified by Daniel Romero-Alvarez to create convex hulls in the environmental space with more than two environmental dimensions. 

Hypervolume loop-post-loop: Script to use one-class support vector machine hypervolums and obtain a model with all the occurrences and calculate the inner performance using a bootstrap approach of 50 replicates. The loop output include performance statistics and descriptive rasters including the median, the 2.5 percentile of the boostrap, the 97.5 percentile of the boostrap, their difference (range), and their sum. 

Convex hull loop-post-loop: Script to use convex hulls. Same as previous description but focused on convex hulls. 

Occurrence contribution: This script calculates the proportion of geographical prediction changing by eliminating each occurrence in a Jackknife approach. The approach is the same but the algorithm use change. 

Environmental values per occurrence and visualization: This portion of the script calculates the environmental contribution across all the occurrence except those identified as impacting the geographical prediction in more than 10% as described in the main text of the manuscript. This analysis can be seen using a ggplot2 visualization as depicted in the publication. 

Randomization tests for vegetation loss: Here, randomization tests compare the observed values of vegetation across OROV occurrences versus the values of 1000 random draws of 35 points (the null distribution) to address statistical significant differences. The same approach is used for NVDI and EVI. 

Calculating the amount of people at risk: Rasters from the Worldpop database are used to calculate the amount of population pixels overlapping with the potential distribution of OROV from the raster selected as the best model. Also, here you can obtain a proxy of incidence by calculating the amount of population pixels from the OROV distribution divided by the total amount of pixels in the second administrative unit of countries of Latin America to present information as a choropleth map. There is a strategy to do this visualization in R using sf objects. 
