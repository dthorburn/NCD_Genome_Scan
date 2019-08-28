# NCD Input Generation and NCD1

The scripts published here are based on the NCD statistic from the publication [Bitarello et al. 2019 in GBE](https://academic.oup.com/gbe/article/10/3/939/4938688). The original scripts upon which mine are based are published [here](https://github.com/bbitarello/NCD-Statistics). 

Additionally, all my scripts are heavily annotated and written for use in R. 

# Input Generation

The ```NCD_Gen``` function has some alterations compared to the one used by the original authors. The main alterations are:
* Allele frequencies are calculated for all possible variants at a site. I have added columns for sites where up to 3 alternate variants (or 5 if you include wildcards like GATK's ```*```).
* ```MAF``` is the lowest value between all 5 ```AF``` columns.

My script works using the ```vcfR``` package to access the genotype field in VCFs. I have only tested this on unphased data, where the genotype field is separated by ```/```, rather than ```|```, but both should be handled correctly. 

As it currently stands, the ```NCD_Gen``` function takes ~0.6 minutes to process a single chromosome, and is easily parallelised across multiple populations or across multiple chromosomes using ```mapply``` or ```foreach```. The problem would be in memory allocation when loading multiple VCFs. Hence, the example I put up is a loop where only one population at a time is calculated for each chromosome.  

# NCD1

There are some minor changes I am making to the NCD1 code too. The ealier changes should increase the accuracy of NCD1 due to a more accurate calculation of ```MAF``` being achieved. The major changes I've made are:
* Windows now start at position 1 of the chromosome, rather than the first SNP. This means your first informative site can now be included in multiple windows, rather than just the first - window size dependent. 
* I've added extra columns for ```Chr```, ```Start```, and ```End``` to the final output of NCD1 for ease of use in downstream analyses. As well as re-ordering the output by the ```Start``` column.  
