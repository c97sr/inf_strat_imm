rm(list=ls(all=TRUE))
load('serological.RDA') 
write.table(serological, "E:/Documents/Github/mrc_nihr_phe/source/isltr/main/extract_antibody_titres/serological_11.21.csv", sep=",") 