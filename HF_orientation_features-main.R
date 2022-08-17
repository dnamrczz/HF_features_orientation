library(CircStats)
library(fields)
library(circular)
library(tidyverse)
library(ggpubr)
library(readr)
library(grid)
library(colorspace)

###############################################
## Code to produce the bearing and plunge diagrams and related 
## statistics can be obtained from McPherron 2018 and Li et al. 2021
## link to the McPherron 2018's GitHub page: https://github.com/surf3s/Orientations
## link to the Li et al. 2021's GitHub page: https://github.com/lili0824/sdg_orientation/tree/main

# import source (from McPherron 2018)
source("orientations.R")

# change the name of my database in "lf" which has been used in the script
lf <- read.csv("HF_twoshots_features.csv", header = T)

#circular_statistics con min_sample = 30
circular_statistics(lf, level = lf$GH, min_sample = 30)
cs_lf_min_sample_30 = circular_statistics(lf, level = lf$GH, min_sample = 30)
str(circular_statistics(lf, level = lf$GH, min_sample = 30))
write.table(cs_lf_min_sample_30, file = "circular_statistics_min_sample_30_HF_orientation_features-main.csv", sep = ",")

#rose diagram for each feature with point for each feature
rose_diagram(plunge_and_bearing(lf %>% filter(GH=='7_IV_BEF6'))$bearing,          #in order to plot a different feature
             pts_on_edge = TRUE, pch = 19,                                    #change 7_IV_BEF6 with another feature name
             bg = "grey75", bar_col = 'grey50', main = "7_IV_BEF6")

#rose diagram plunge
rose_diagram_plunge(plunge_and_bearing(lf %>% filter(GH=='7_IV_BEF6'))$plunge,    #in order to plot a different feature
                    main="Plunge 7_IV_BEF6", bg = 'grey25', bar_col = 'grey75')     #change 7_IV_BEF6 with another feature name

###############################################
## how to make the Benn diagrams - from McPherron (2018) markdown file
# reference orientation data from Lenoble and Bertran (2004)
lenoble_and_bertran = readRDS('Lenoble_and_Bertran_2004.RDS') %>%
  filter(Type %in% c('Debris Flow','Runoff Shallow','Runoff Steep','Solifluction'))

# color for Lenoble and Bertran reference
formation_process_colors = adjustcolor(c('#404096','#57A3AD','#DEA73A','#D92120'), alpha = .3) 

# color for Lenoble and Bertran reference
lenoble_and_bertran_colors = factor(lenoble_and_bertran$Type, labels = formation_process_colors) 

p = .95 # p value for resampling - McPherron 2018

resampling = 10000	 # Number of resampling iterations

# Summary stats for isotropic ratio and elongation ratio for each layer
benn_lf = round(benn(xyz = lf, level = lf$GH, min_sample = 30), 3) 

layout(matrix(c(1,3,2,4), nrow = 2, ncol = 2)) # configure layout of the plot space

for (levelname in (unique(lf$GH))) {   # plot orientation ternary plot for each layer and draw area of bootstrapped zone.
  
  xyz_level = lf %>% filter(GH==levelname)
  
  benn_diagram(cbind(benn_lf[levelname,"EL"],benn_lf[levelname,"IS"]),
               main = paste("Layer",levelname),
               cex = .7,
               labels = ifelse(levelname == (unique(lf$GH))[1], 'outside', 'none'),
               dimnames_position = ifelse(levelname == (unique(lf$GH))[1], 'corner', 'none'),
               grid_labels = ifelse(levelname == (unique(lf$GH))[1], TRUE, FALSE))
  
  benn_diagram(lenoble_and_bertran, 
               drawhull = TRUE,
               new_page = FALSE,
               plot_points = FALSE, cex = .7,
               col = lenoble_and_bertran_colors,
               legend_names = switch((levelname == (unique(lf$GH))[1]) + 1,
                                     NULL, levels(factor(lenoble_and_bertran$Type))),
               legend_colors = formation_process_colors)
  
  resampling_contours = benn_resampling(xyz_level, resampling = resampling, p = p) 
  
  for (contour in resampling_contours)
    lines(benn_coords(cbind(elongation = contour$x, isotropy = contour$y))) }

###############################################
## Code to run the Kuiper's test

# run the Kuiper's test on all bearings as an example to show how
# the function works in R.
# call the plunge_and_bearing function from orientations.R to process the data
# and get bearings only

bearing_angle = plunge_and_bearing(lf)$bearing
bearing_angle = circular(bearing_angle, type = "angles", units = "degrees")

# according to the function's official documentation, alpha can be set to
# 0.15, 0.1, 0.05, 0.025, 0.01. as the significance level of the test.
# when an alpha level is omitted, a range for the p-value will be provided.

kuiper.test(bearing_angle, alpha = 0.05)
kuiper.test(bearing_angle)

# run the kuiper's test on each individual cultural layers
# GH7_IV_BEF6 data
GH7_IV_BEF6 = lf %>%
  filter(GH == '7_IV_BEF6')
GH7_IV_BEF6_bearing = plunge_and_bearing(GH7_IV_BEF6)$bearing
GH7_IV_BEF6_bearing_df = as.data.frame(GH7_IV_BEF6_bearing)
kuiper(circular(GH7_IV_BEF6_bearing, type = "angles", units = "degrees"), alpha = 0.05)

###############################################
# generate the histogram of bearings by GH
# GH7_IV_BEF6 histogram
GH7_IV_BEF6_bearing_plot = 
  ggplot(GH7_IV_BEF6_bearing_df, aes(GH7_IV_BEF6_bearing)) +
  geom_histogram(aes(y = ..density..), 
                 binwidth = 20, color = "grey30", fill = "grey87") +
  geom_density(alpha = .2, fill = "antiquewhite3") +
  theme_bw() +
  xlab("Bearing angle (degree)") +
  ylab("Density")
GH7_IV_BEF6_bearing_plot




