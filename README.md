# ProstateCancer_TCGA_VCF
This project aims at processing and analyzing genetic sequence variants in VCF formats.
The op_vcf.py transforms each sample in VCF files into a corresponding csv file. It requires 4 parameters as follows:
(1) -i --input_dir: full path to all vcf files to be processed.
(2) -o --out_dir: folder to store csv files
(3) -f --ref_flat: full path to the refflat file, including the file name. This folder should also contain reference sequence files split by chromosome.
(4) -c --num_cores: optional, number of cores to be used.
VCF2CSV.py, parse_VCF.py, miscVar.py need to be put at the same folder as op_vcf.py.

fathmm_submission.py takes an OMI file, submit it to FATHMM for predictions and download and parse results. It requires 2 required and 2 optional parameters as follows:
(1) -i --input_file: full path to the input OMI file, including the file name.
(2) -o --output_dir: folder to store FATHMM result file. It could be the same directory. The output file is named as OMI filename+"_fathmm.csv".
(3) -r --reference build: optional, specify reference genome build to be used
(4) -t --runtime: optional, specify sleeptime for website running

gsv_count_subject.py takes a list of OMI files and count for each variants, the number of files containg that variant, and the distribution of this variant among all files. It requires 3 parameters as follows:
(1) -t --input_dir:full path to all tumor OMI files.
(2) -n --input_dir:full path to all normal OMI files. It could be the same as tumor directory.
(3) -o --out_dir: folder to store result file. The file is named as sub_count_dist.csv

provean.py takes nonsynonymous variants and generate sequences files and variation files for PROVEAN predictions. It takes 5 parameters as follows:
(1) -i --input_file: CSV file for containing all nonsynonymous variants.
(2) -s --outseq_dir: folder to store all amino acid sequences.
(3) -v --outvar_dir: folder to store all amino acid variations.
(4) -r --ref_sequence: reference sequence file used to translate aa,including path.
(5) -f --ref_flat: reference genomic field used to identify gene region,including path.
