###############################################################################
####################___CMI_LAB_EMU_AUTOMATED_PLOTTING____######################
###############################################################################

### Created by Austin G. Marshall and Daniel T. Fuller ###

###############################################################################

### Load libraries for the packages used in this workflow ###
library(dplyr)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggsci)
library(janitor)
library(phyloseq)

###############################################################################

### Loading in combined abundance data from Emu ###
a = read.delim(file = "emu-combined-abundance-species.tsv")

###############################################################################
############___Change_the_order_around_to_plot_in_certain_order___#############
###############___Yes_this_is_ugly_but_who_cares_if_it_works___################
###############################################################################

### Order samples how you want them displayed in the plot ###
b = a %>%
  select(species, RT_30, RT_34, RT_36) %>%
  mutate_all(~replace(., is.na(.), 0))

### Grabbing the data without the blank taxa ###
without_blanks = b %>%
  na_if("") %>%
  na.omit()

### Grabbing the blank taxa only ###
with_blanks = b %>%
  filter(species=="") %>%
  na_if("")

### Summing the blanks together ###
blanksums = with_blanks %>%
  select(-species) %>%
  colSums() %>%
  t() %>%
  as.data.frame()

### Setting up the naming stuff for the filtering step ###
rownames(without_blanks)<-without_blanks$species
before_filt <- without_blanks %>%
  select(-species)

###############################################################################
################ Applying filter and binding rows together ####################
###############################################################################

abund_filt <- apply(before_filt, 2, FUN = function(x) (x/sum(x))<=0.01)
data.cleaned <- before_filt[-which(rowSums(abund_filt)==ncol(before_filt)),]
othercat <- colSums(before_filt[which(rowSums(abund_filt)==ncol(before_filt)),])
bind1 <- rbind(data.cleaned, othercat)
bind2 <- rbind(bind1, blanksums)
final_w_name = bind2 %>%
  rownames_to_column(var = "Sample")

###############################################################################
##########__Stop_here_and_change_what_samples_need_to_be_grabbed__#############
###############################################################################

grab_sums = bind2[c('10','1'),]
grab_blanksums = filter(final_w_name, final_w_name$Sample==1)
grab_taxasums = filter(final_w_name, final_w_name$Sample==10)

### Need to change these numbers to correspond to the number of samples ###
final_other = grab_taxasums[2:4] - grab_blanksums[2:4]

### Need to change the numbers to account for changing datasets ### ### Stuck here on big one ###
bind2_minusothers <- bind2[c(1:9),]

done <- rbind(bind2_minusothers, final_other) 
rownames(done)[length(done[,1])]="Other"

done2 = done %>%
  arrange(RT_30)

done3 = t(done2) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Sample")

###############################################################################
###### Gotta change the formatting once more for ggplot2 to like it ###########
###############################################################################

tbl_reshape <- melt(done3, id.vars = "Sample" , variable.name = "Species")

### Now We plot, change the title to whatever your experiment is ###

speciesplot <- ggplot(tbl_reshape, aes(x = fct_inorder(Sample), y = value, fill = Species)) +
  geom_bar(colour = "black", stat = "identity", position = "fill") +
  labs(x = "Sample", y = "Relative abundance", title = "const conc diff cycles arranged") +
  theme_bw() 

###############################################################################
#################### Change color patterns using ggsci ########################
###############################################################################

speciesplot_goodone = speciesplot + scale_fill_ucscgb(palette = "default", alpha = 1) 

### Final plotting of the figure ###

speciesplot_goodone

###############################################################################

bluh = read.delim("S1_rel-abundance-species.tsv")
