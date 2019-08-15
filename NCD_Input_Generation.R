###########################################################################
#		Generate NCD Input from VCF Files - Unphased										
###########################################################################
## Date: 08/08/2019
## Involved: Miles
## Task: Parse an unphased VCF file for it's genotype information to generate a table of allele frequencies
## Note: This has not been tested on Phased genotype data. It should work, but check the format of the gt field
## 	 in vcfR objects to make sure there are no other numbers present.

library(vcfR)
library(dplyr)
library(stringr)
library(data.table)

vcf_dir <- "~/QM_PhD/Chapter1/VCFs/Attempt3/All_Ind_VCFs/"
setwd("~/NCD/NCD_Input")

## Create a population dataframe, could also use a reference panel. What you need ultimately is a vector with the column 
## names in the vcf for your vcf for each population. Alternatively, you could split your VCFs prior to this step. 
Population_list <- data.frame(  CA_L=c("BS44", "BS46", "BS48", "BS50", "BS52", "BS54"), 
								CA_R=c("BS43", "BS45", "BS47", "BS49", "BS51", "BS53"), 
								DE_M=c("BS25", "BS26", "BS27", "BS28", "BS29", "BS30"), 
								G1_L=c("BS2", "BS4", "BS6", "BS8", "BS10", "BS12"), 
								G1_R=c("BS1", "BS3", "BS5", "BS7", "BS9", "BS11"), 
								G2_L=c("BS14", "BS16", "BS18", "BS20", "BS22", "BS24"), 
								G2_R=c("BS13", "BS15", "BS17", "BS19", "BS21", "BS23"), 
								NO_L=c("BS56", "BS58", "BS60", "BS62", "BS64", "BS66"), 
								NO_R=c("BS55", "BS57", "BS59", "BS61", "BS63", "BS65"), 
								US_L=c("BS32", "BS34", "BS36", "BS38", "BS40", "BS42"), 
								US_R=c("BS31", "BS33", "BS35", "BS37", "BS39", "BS41"))

## A function that uses the piped input in two separate functions, since dplyr struggled to do this even using the (.) signifier. 
## The GT_separator is what separates gt info, I think this changes to | in phased VCFs. 
Freq_Fun <- function(input, value){
	str_count(string = input, pattern = value)/str_count(string = input, pattern = "[0-9]")
}

## A vectorised function that writes into data.table format by parsing the "fix" and "gt" fields in the vcfR object.
## AFx are calculated by only taking the information is the gt field that designates the nucleotides at that position for each individual
## then collpasing all individuals into a single string and counting the number of 0s, 1s, 2s, 3s, and 4s+ for AF1-5, respectively.
## It then divides the number by the number of numbers in the gt field string. This accounts for individuals with missing information
## at that site, expressed as "./.".
## AF1-5 fields are replaced with NA's if the number is 0. Otherwise MAF will always be 0, even if, for example, AF3 and AF4 are 0.25 and 0.75
## because the AF1 is 0.00. Thus, the default MAF for an invariant site is now 1.00. This is handled fine in NCD1s code.   
NCD_Gen <- function(vcf_obj, population_ID){
	output <- data.table(CHR = vcf_obj@fix[,"CHROM"], POS = vcf_obj@fix[,"POS"] %>% as.integer(), 
			     				ID = paste(vcf_obj@fix[,"CHROM"], vcf_obj@fix[,"POS"], sep = "|"),
							REF = vcf_obj@fix[,"REF"], 
			     				ALT = vcf_obj@fix[,"ALT"],
							AF1 = gsub(pattern = ":.*", replacement =  "", vcf_obj@gt[,as.character(Population_list[,population_ID])]) %>% as.data.frame() %>%
									apply(MARGIN = 1, FUN = paste, collapse = "|") %>% Freq_Fun("0") %>% gsub(pattern = "^0$", replacement = NA) %>% as.numeric(),
							AF2 = gsub(pattern = ":.*", replacement =  "", vcf_obj@gt[,as.character(Population_list[,population_ID])]) %>% as.data.frame() %>%
									apply(MARGIN = 1, FUN = paste, collapse = "|") %>% Freq_Fun("1") %>% gsub(pattern = "^0$", replacement = NA) %>% as.numeric(),
							AF3 = gsub(pattern = ":.*", replacement =  "", vcf_obj@gt[,as.character(Population_list[,population_ID])]) %>% as.data.frame() %>%
									apply(MARGIN = 1, FUN = paste, collapse = "|") %>% Freq_Fun("2") %>% gsub(pattern = "^0$", replacement = NA) %>% as.numeric(),
							AF4 = gsub(pattern = ":.*", replacement =  "", vcf_obj@gt[,as.character(Population_list[,population_ID])]) %>% as.data.frame() %>%
									apply(MARGIN = 1, FUN = paste, collapse = "|") %>% Freq_Fun("3") %>% gsub(pattern = "^0$", replacement = NA) %>% as.numeric(),
							AF5 = gsub(pattern = ":.*", replacement =  "", vcf_obj@gt[,as.character(Population_list[,population_ID])]) %>% as.data.frame() %>%
									apply(MARGIN = 1, FUN = paste, collapse = "|") %>% Freq_Fun("[4-9]") %>% gsub(pattern = "^0$", replacement = NA) %>% as.numeric())
	## Finally, once the AFs are calculated, this creates the MAF column by taking the min values of the AF columns. 
	output %>% mutate(MAF = pmin(AF1, AF2, AF3, AF4, AF5, na.rm = TRUE)) %>% as.data.table()
}

## Loop to calculate and write output as Rdata object. It makes most sense to load one chromosome at a time, rather than loading and unloading chromosome VCFs for each population.
## Depending on size of output, it might then be worthwhile outputting, then another loop to create autosomal list for each population. 
## Note: vcfR can read in gunzipped VCFs
## If you are doing this on a computing cluster. I would recommend reading in all chromosomes, and using lapply or foreach to parallelise this step. It you're on a 64-bit operating
## system, set memory to as high as it can go if you have several populations you want to run at the same time. The loop is only because it would take too much space to load in 
## all the VCFs. 
for(chrom in c(1:22)){
	start_time <- Sys.time(); if(chrom == 1){ult_start <- Sys.time()}
	if(chrom != 22){
		temp_chrom <- read.vcfR(paste(vcf_dir,"All_chr", as.roman(chrom),".filtered.selected.SNPs.vcf.gz", sep = ""))
	} else if(chrom == 22){
		temp_chrom <- read.vcfR(paste(vcf_dir,"All_chrUn.filtered.selected.SNPs.vcf.gz", sep = ""))
	}
	checkpoint <- Sys.time(); print(checkpoint - start_time)
	for(pops in names(Population_list)){
		if(pops != "DE_M"){
			print(paste("Starting extraction for", pops))
			## This is my attempt to not load and reloading the list objects. If the datasets are really large, this may cause R to crash by loading to much into the environment.
			## Alternatively, you could save each output as an RData file and then reload it, but it's a lot of waster memory
			## allocation. 
			if(chrom == 1){
				temp_list <- list()
				assign(pops, temp_list)
			}	else	{
				temp_list <- get(pops)
			}
			temp_list[[chrom]] <- NCD_Gen(temp_chrom, pops)
			checkpoint2 <- Sys.time(); print(checkpoint2 - checkpoint); print(paste("from start of Chr", chrom, sep = ""))
			assign(pops, temp_list)
			## saving the output
			if(chrom == 22){
				save(temp_list, file = paste("All_autosome", pops, "NCD_Input_Table_AllAFs.RData", sep = "_"))
				print(paste("Saved output for", pops))
			}
		}
	}
## Just a timer to see how long the entire data set took to process
if(chrom == 22){ult_end <- Sys.time(); print(ult_end - ult_start)}
}

