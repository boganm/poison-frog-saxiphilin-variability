#load libraries
library(tidyverse)
library(ggplot2)

#load datasheets and remove N/A's
binding_all_values_n1 <- as.data.frame(read.csv('R_FORMAT_N1_protopocket.csv'))
binding_all_values_c1 <- as.data.frame(read.csv('R_FORMAT_STX_binding_site.csv'))
binding_avgs_n1 <- as.data.frame(read.csv('R_FORMAT_N1_protopocket_avg.csv'))
binding_avgs_c1 <- as.data.frame(read.csv('R_FORMAT_STX_binding_site_avg.csv'))

#define color list
spec_colors <- c("brown", "blue", "orange","green")

#define theme for label sizes
my_theme <- theme(plot.title = element_text(size=40), 
                  axis.title.y = element_text(size = rel(3)), 
                  axis.title.x = element_text(size = rel(3)),
                  axis.text.x = element_text(size = rel(3.2)), 
                  axis.text.y = element_text(size = rel(3.2)),
                  legend.text = element_text(size = rel(2.3)))

#define fill colors
my_fill <- scale_fill_manual("legend", values = c("A. fem" = "tan", 
                                                  "D. tin" = "dark blue", 
                                                  "O. syl" = "dark orange", 
                                                  "R. cat" = "dark green"))

################################CODE TO EXECUTE################################

##Average H-bonds predicted at the N1-lobe proto-pocket
ggplot(binding_avgs_n1, aes(x = Compound, y = mean.h_bonds, fill = Species)) +
  geom_col(data = binding_avgs_n1, position = "dodge") +
  my_fill +
  ylab('H-bonds predicted') +
  ggtitle("N1-lobe proto-pocket H-bond prediction") + 
  my_theme

##Average H-bonds predicted at the C1-lobe STX-binding site
ggplot(binding_avgs_c1, aes(x = Compound, y = mean.h_bonds, fill = Species)) +
  geom_col(data = binding_avgs_c1, position = "dodge") +
  my_fill +
  ylab('H-bonds predicted') +
  ggtitle("STX-binding site H-bond prediction") + 
  my_theme

##Average binding score predicted at the STX-binding site
ggplot(binding_avgs_c1, aes(x = Compound, y = X.mean_score.,fill = Species)) +
  geom_col(data = binding_avgs_c1, position = "dodge") +
  my_fill +
  ylab('|mean binding score|') +
  ggtitle("STX-binding site scoring prediction") + 
  my_theme

##Average binding score predicted at the N1-lobe proto-pocket
ggplot(binding_avgs_n1, aes(x = Compound, y = X.mean_score., fill = Species)) +
  geom_col(data = binding_avgs_n1,position = "dodge") +
  my_fill +
  ylab('|mean binding score|') +
  ggtitle("N1-lobe proto-pocket scoring prediction") +
  my_theme