# ProstateCancer_TCGA_VCF
This project aims at processing and analyzing genetic sequence variants in VCF formats.
The op_vcf.py transforms each sample in VCF files into a corresponding csv file. It requires 4 parameters as follows:
(1) vcf directory: full path to all vcf files to be processed.
(2) output directory: folder to store csv files
(3) gtf file: full path to the refflat file, including the file name. This folder should also contain reference sequence files split by chromosome.
(4) number of cores: number of cores to be used, usually one less than total cores.
