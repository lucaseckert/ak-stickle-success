## README: raw-data

This folder contains all raw genotyping and environmental data, organized into the following folders:

-   `env-data/`: Contains all environmental data and ecological data, mostly water chemistry, macroinvertebrates, and zooplankton. The `env-data-cleaning.R` script is used to compile and clean this data, for further analysis.

-   `snp-array-data/`: Contains all raw genotyping data from the SNP arrays, including the raw genotyping calls and information on the set of SNPs. The `ancestry-inference.R` script is used to infer proportional ancestry of each fish from these data.
