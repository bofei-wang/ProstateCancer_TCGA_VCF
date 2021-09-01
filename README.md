# ProstateCancer_TCGA_VCF
This project aims at processing and analyzing genetic sequence variants in VCF formats.
The op_vcf.py transforms each sample in VCF files into a corresponding csv file. It requires 4 parameters as follows:
(1) vcf directory: full path to all vcf files to be processed.
(2) output directory: folder to store csv files
(3) gtf file: full path to the refflat file, including the file name. This folder should also contain reference sequence files split by chromosome.
(4) number of cores: number of cores to be used, usually one less than total cores.
VCF2CSV.py, parse_VCF.py, miscVar.py need to be put at the same folder as op_vcf.py.

fathmm_submission.py takes an OMI file, submit it to FATHMM for predictions and download and parse results. It requires 2 parameters as follows:
(1) OMI file: full path to the OMI file, including the file name.
(2) Output directory: folder to store FATHMM result file. It could be the same directory. The output file is named as OMI filename+"_fathmm.csv".

gsv_count_subject.py takes a list of OMI files and count for each variants, the number of files containg that variant, and the distribution of this variant among all files. It requires 3 parameters as follows:
(1) tumor directory:full path to all tumor OMI files.
(2) tumor directory:full path to all normal OMI files. It could be the same as tumor directory.
(3) output directory: folder to store result file. The file is named as sub_count_dist.csv
