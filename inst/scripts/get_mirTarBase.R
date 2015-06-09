# The following script downloads and constructs the mirTarBase dataset,
# included in the mirIntegrator package. 
# This dataset was published by Hsu SD et al. (2014) Nucleic acids research.
# This dataset is included for demostration purposes and it should be used 
# according to its licensing terms:
# http://mirtarbase.mbc.nctu.edu.tw/cache/download/LICENSE

# 1. Download the file from mirTarBase website.

download.file(url = "http://mirtarbase.mbc.nctu.edu.tw/cache/download/4.5/hsa_MTI.xls", 
              destfile = "~/hsa_MTI.xls")
# 2. Load the table
library("gdata")  # gdata needs perl
mirTarBase <- read.xls ("~/hsa_MTI.xls", sheet = 1, header = TRUE)

# 3. Rename columns
col.nam <- names(mirTarBase) 
col.nam[5] <- "Target.ID"
names(mirTarBase) <- col.nam


# 4. Delete the empty rows
targets_db <- mirTarBase[,c("miRNA","Target.ID")]
mirTarBase <- mirTarBase[!(is.na(targets_db[,2])),]

# 5. Save files
save(mirTarBase, file = "mirTarBase.rda")
library(tools)
resaveRdaFiles("mirTarBase.rda")

# 6. Now you can delete the downloaded file from your home directory
if (file.exists("~/hsa_MTI.xls")) file.remove("~/hsa_MTI.xls")
