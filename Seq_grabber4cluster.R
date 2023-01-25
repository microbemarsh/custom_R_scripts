###########################################################################################
##########__Austin Marshall_PhD_student_Dept_of_Biology_Clarkson_University__##############
###########################################################################################


### This script is used to read in the bam file after samtools conversion from the sam file produced by Emu minimap2 step

### We load in the bam file and using a modified script from https://gist.github.com/davetang/6460320 convert to dataframe

### We then use the Emu database reference number to find the species hit of interest and select only those sequences to 
### provide a better clustering with NGSpeciesID (less noise allows for higher accuracy)

library(Rsamtools)
library(dplyr)
library(data.table)
library(Biostrings)
library(tgutil)

#read in entire BAM file
bam <- scanBam("S1_emu_aln_sorted.bam")

#names of the BAM fields
names(bam[[1]])
# [1] "qname"  "flag"   "rname"  "strand" "pos"    "qwidth" "mapq"   "cigar"
# [9] "mrnm"   "mpos"   "isize"  "seq"    "qual"

#distribution of BAM flags
table(bam[[1]]$flag)

#      0       4      16 
#1472261  775200 1652949

#function for collapsing the list of lists into a single list
#as per the Rsamtools vignette
.unlist <- function (x){
  ## do.call(c, ...) coerces factor to integer, which is undesired
  x1 <- x[[1L]]
  if (is.factor(x1)){
    structure(unlist(x), class = "factor", levels = levels(x1))
  } else {
    do.call(c, x)
  }
}

#store names of BAM fields
bam_field <- names(bam[[1]])

#go through each BAM field and unlist
list <- lapply(bam_field, function(y) .unlist(lapply(bam, "[[", y)))

#store as data frame
bam_df <- do.call("DataFrame", list)
names(bam_df) <- bam_field

dim(bam_df)
#[1] 3900410      13

S1 = data.frame(bam_df)

filt_S1 = S1 %>%
  select(qname, rname, seq, qual)%>%
  na_if("") %>%
  na.omit()

main_ral <- filt_S1[filt_S1$rname %like% "305:emu_db:16544", ] 
colnames(main_ral)[1] = "id"

almost_fastq = main_ral %>%
  select(id, seq, qual)

shiz = gsub("^", "@", as.character(almost_fastq$id))

shiz2 = shiz %>%
  as.data.frame()
colnames(shiz2)[1] = "id_good"

soclose = cbind(shiz2, almost_fastq) %>%
  select(id_good, seq, qual)
colnames(soclose)[1] = "id"

### Need to add the @ symbol to beginning of the read ids before writing

### This final writing of the fastq is what will be loaded into NGSpeciesID

write_fastq(soclose, "S1_main_ralstonia.fastq")

#################################################################

