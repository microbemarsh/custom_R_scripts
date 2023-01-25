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
a = read.delim(file = "emu-combined-abundance-phylum.tsv")

###############################################################################
############___Change_the_order_around_to_plot_in_certain_order___#############
###############___Yes_this_is_ugly_but_who_cares_if_it_works___################
###############################################################################

### Order samples how you want them displayed in the plot ###
b = a %>%
  select(phylum, S1, S2, S3, S4, S5) %>%
  mutate_all(~replace(., is.na(.), 0))

### Grabbing the data without the blank taxa ###
without_blanks = b %>%
  na_if("") %>%
  na.omit()

### Grabbing the blank taxa only ###
with_blanks = b %>%
  filter(phylum=="") %>%
  na_if("")

### Summing the blanks together ###
blanksums = with_blanks %>%
  select(-phylum) %>%
  colSums() %>%
  t() %>%
  as.data.frame()

### Setting up the naming stuff for the filtering step ###
rownames(without_blanks)<-without_blanks$phylum
before_filt <- without_blanks %>%
  select(-phylum)

###############################################################################
################ Applying filter and binding rows together ####################
###############################################################################

abund_filt <- apply(before_filt, 2, FUN = function(x) (x/sum(x))<=0.005)
data.cleaned <- before_filt[-which(rowSums(abund_filt)==ncol(before_filt)),]
othercat <- colSums(before_filt[which(rowSums(abund_filt)==ncol(before_filt)),])
bind1 <- rbind(data.cleaned, othercat)
bind2 <- rbind(bind1, blanksums)
final_w_name = bind2 %>%
  rownames_to_column(var = "Sample")

###############################################################################
##########__Stop_here_and_change_what_samples_need_to_be_grabbed__#############
###############################################################################

grab_sums = bind2[c('7','1'),]
grab_blanksums = filter(final_w_name, final_w_name$Sample==7)
grab_taxasums = filter(final_w_name, final_w_name$Sample==1)

### Need to change these numbers to correspond to the number of samples ###
final_other = grab_taxasums[2:6] - grab_blanksums[2:6]

### Need to change the numbers to account for changing datasets ### ### Stuck here on big one ###
bind2_minusothers <- bind2[c(1:6),]

done <- rbind(bind2_minusothers, final_other) 
rownames(done)[length(done[,1])]="Other"

done2 = done %>%
  arrange(S1)

done3 = t(done2) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Sample")

###############################################################################
###### Gotta change the formatting once more for ggplot2 to like it ###########
###############################################################################

tbl_reshape <- melt(done3, id.vars = "Sample" , variable.name = "Phylum")

### Now We plot, change the title to whatever your experiment is ###

phylumplot <- ggplot(tbl_reshape, aes(x = fct_inorder(Sample), y = value, fill = Phylum)) +
  geom_bar(colour = "black", stat = "identity", position = "fill") +
  labs(x = "Sample", y = "Relative abundance", title = "Varying conc. and cycles RDP") +
  theme_bw() 

###############################################################################
#################### Change color patterns using ggsci ########################
###############################################################################

phylumplot_goodone = phylumplot + scale_fill_ucscgb(palette = "default", alpha = 1) 

### Final plotting of the figure ###

phylumplot_goodone
