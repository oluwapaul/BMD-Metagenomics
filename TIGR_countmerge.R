#install.packages('stringr', repos="http://ftp.ussg.iu.edu/CRAN/")
#install.packages('tidyr', repos="http://ftp.ussg.iu.edu/CRAN/")
#install.packages('dplyr', repos="http://ftp.ussg.iu.edu/CRAN/")
library("stringr")
library('tidyr')
library("dplyr")

setwd("./TIGR_counts/")
sample_list <- list.files()
sample_list_short <- str_replace(sample_list, "_mapped_sorted.bam.TIGR.counts", "")
n_genes <- 653739
count <- 0
rm(TIGR_count_table)
for (s in sample_list){
  count <- count + 1
  temp <- read.table(s, sep = "\t", nrows = n_genes, stringsAsFactors = FALSE)
  temp <- temp[,c(1,7)]
  colnames(temp) <- c("geneID", sample_list_short[count])
  print(sample_list[count])
  if(exists("TIGR_count_table")){
    print("yes")
    TIGR_count_table <- merge(TIGR_count_table, temp, by.x = "geneID", by.y = "geneID")
    print(count)
  } else {
    print("no")
    print(count)
    TIGR_count_table <- temp
  }
}

write.table(TIGR_count_table, "../OutputFiles/TIGR_count_table.txt", sep = "\t", quote = F, row.names = F)