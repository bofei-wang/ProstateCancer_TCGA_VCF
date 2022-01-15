from multiprocessing import Pool,Manager
import os  
from math import fabs         
import sys   
from miscVar import AADict
import time
import csv
import gc

MAX_DIST = 5000
#function for processing vcf, returns dict
def vcf_func(vcffile,gtf,genome):     
    chr_start = time.time()
    vcffile = open('%s'%(vcffile), 'r')
    vcf = vcffile.readlines()
    vcffile.close()
    gtffile = gtf    
    f = open(gtffile,'r')
    flat = f.readlines()
    f.close()
    seqfile=open(genome,'r')
    next(seqfile)
    sequence=seqfile.read().replace('\n','')
    seqfile.close()
    
    dictionary= {} 
    position = 0
    chrom=''
    #genestart = 0
    #genestop = 0
    start_time = time.time()
    for l in range(0, len(vcf)): #line in vcf:
        if vcf[l].startswith('chr'):
            line = vcf[l].strip()
            vl = line.split() #split on all white space
            chrom=vl[0]
            position = int(vl[1])
            ref=vl[3]
            alt=vl[4]
            gtfCounter = -1 #point to each line in gtf file
            genestart = 0
            genestop = 0
            if position not in dictionary:
                boo = 0
                where = []
                genename = []
                change_type=[]
                aachange=[]
                counter=[]
                strand_gene=[]
                while boo == 0 and gtfCounter < len(flat):
                    if position >genestart and position <= genestop and gtfCounter < (len(flat)-1): #GSV within gene
                        counter.append(gtfCounter)
                        strand_gene.append(flat[gtfCounter].split()[3])
                        #if (gtfCounter < (len(flat)-1)):
                        gtfCounter += 1
                        pgstop = genestop
                        genestart = int(flat[gtfCounter].split()[4]) #transcript start
                        genestop = int(flat[gtfCounter].split()[5]) #transcript end
                    elif position > genestop and gtfCounter < (len(flat)-1):  #Moves to the next gene
                        gtfCounter += 1
                        pgstop = genestop #prev gene stop
                        genestart = int(flat[gtfCounter].split()[4]) #transcript start
                        genestop = int(flat[gtfCounter].split()[5]) #transcript end
                    elif position < genestart and position > pgstop:
                        boo = 1                       
                    else:
                        boo = 1
                if len(strand_gene)==0:
                    strand_gene.append('')
                strand = flat[gtfCounter].split()[3] #ref strand direction
                run_time = time.time() - start_time
                if len(counter)==0:
                    if gtfCounter == len(flat)-1 and position > genestop:
                        if (position - genestop) > MAX_DIST: #if variant position is far from the end of transcript
                            genename.append('NoName')
                            where.append('Not close to gene')
                            change_type.append('')
                            aachange.append('')
                        else:
                            genename.append(flat[gtfCounter].split()[0]) #else whether close to 3' or 5'
                            if strand == '-':
                                where.append('5\' untranscribed')
                                change_type.append('')
                                aachange.append('')
                            else:
                                where.append('3\' untranscribed')
                                change_type.append('')
                                aachange.append('')
                    elif position < genestart and position > pgstop:
                        # check if closer to start or end of the one before or current
                        diff1 = fabs(position-genestart) #current
                        diff2 = fabs(position-pgstop) #previous
                        if diff1 <= diff2 and diff1 <= MAX_DIST: #if closer to current "new " gene and within cutoff distance
                            genename.append(flat[gtfCounter].split()[0])
                            if strand == '-':
                                where.append('3\' untranscribed')
                                change_type.append('')
                                aachange.append('')
                            else:
                                where.append('5\' untranscribed')
                                change_type.append('')
                                aachange.append('')
                        elif diff2 <= MAX_DIST: #if closer to prev gene and within cutoff distance
                            genename.append(flat[gtfCounter-1].split()[0])
                            strand = flat[gtfCounter-1].split()[3]
                            if strand == '-':
                                where.append('5\' untranscribed')
                                change_type.append('')
                                aachange.append('')
                            else:
                                where.append('3\' untranscribed')
                                change_type.append('')
                                aachange.append('')
                        else:
                            genename.append('NoName')
                            where.append('Not close to gene')
                            change_type.append('')
                            aachange.append('')
                    else:
                        continue
                else:
                    for i in range(len(counter)):
                        genename.append(flat[counter[i]].split()[0])
                        #tr_start = int(flat[counter[i]].split()[4]) #transcript  start; check
                        #tr_end = int(flat[counter[i]].split()[5]) #transcript  end; check
                        cd_start=int(flat[counter[i]].split()[6]) #CDS start
                        cd_end=int(flat[counter[i]].split()[7]) #CDS end
                        ex_starts = []
                        ex_ends=[]
                        x=0
                        tmp = 'Unknown'
                        ex_starts = [int(i) for i in flat[counter[i]].split()[-2].split(',')[0:-1]]
                        ex_ends = [int(i) for i in flat[counter[i]].split()[-1].split(',')[0:-1]]

                        while tmp == 'Unknown' and x < len(ex_starts):
                            if position > ex_starts[x] and position <= ex_ends[x]: #if position in exon
                                tmp = 'Exon'
                            elif position <= ex_starts[x]:
                                tmp = 'Intron'
                            else:
                                x+=1
                        if tmp == 'Intron' or tmp=='Unknown':
                            where.append('Intron')
                            change_type.append('')
                            aachange.append('')
                        #see if position is in translated region:
                        elif tmp == 'Exon':
                            if position > cd_start and position <= cd_end:
                                tmp = 'Translated'
                            else: #not translated
                                tmp = 'Not translated'
                            if tmp == 'Not translated':
                                if cd_start == cd_end: #ncRNAs
                                    where.append('Not translated, ncRNA')
                                    change_type.append('')
                                    aachange.append('')
                                else: #utrs
                                    if position  <= cd_start:
                                        if strand_gene[i] == '-':
                                            where.append('Not translated,3\' utr')
                                            change_type.append('')
                                            aachange.append('')
                                        else:
                                            where.append('Not translated,5\' utr')
                                            change_type.append('')
                                            aachange.append('')
                                    else:
                                        if strand_gene[i] == '-':
                                            where.append('Not translated,5\' utr')
                                            change_type.append('')
                                            aachange.append('')
                                        else:
                                            where.append('Not translated,3\' utr')
                                            change_type.append('')
                                            aachange.append('')
                            elif tmp == 'Translated':
                                where.append('CDS')
                                ref_codon=''
                                var_codon=''
                                aalength=0
                                m=0                                
                                toRev=0
                                #m = (position - cd_start) % 3
                                if strand_gene[i]=="+":
                                    ex=0
                                    #pre_ex=0
                                    while ex != (len(ex_ends)-1) and position > ex_ends[ex]:
                                        if cd_start > ex_ends[ex]:
                                            ex+=1
                                        elif cd_start > ex_starts[ex] and cd_start < ex_ends[ex] and position > ex_ends[ex]:
                                            #pre_ex=ex
                                            m = ((ex_ends[ex] - cd_start )+m)%3
                                            aalength+=(ex_ends[ex] - cd_start)
                                            ex+=1
                                        else:
                                           # pre_ex=ex
                                            m = ((ex_ends[ex]-ex_starts[ex])+m)%3
                                            aalength+=(ex_ends[ex] - ex_starts[ex])
                                            ex+=1
                                    #print("ex is %d" %ex)
                                    if cd_start > ex_starts[ex]:
                                        m=(position-cd_start -1+m )%3
                                        aalength+=(position-cd_start) 
                                    else:
                                        m=(position-ex_starts[ex] -1+m)%3
                                        aalength+=(position-ex_starts[ex])
                                    if aalength%3==0:
                                        aapos=aalength//3
                                    else:
                                        aapos=aalength//3 + 1
                    
                                #if len(ref)==len(alt)and position<len(sequence)-2: 
                                    if m==0 :
                                        if len(ref)==1:
                                            ref_codon=ref+sequence[position]+sequence[position+1]
                                            var_codon=alt+sequence[position]+sequence[position+1]
                                        elif len(ref)==2:
                                            ref_codon=ref+sequence[position]
                                            var_codon=alt+sequence[position]
                                        else:
                                            ref_codon=ref
                                            var_codon=alt
                                    elif m==1:
                                        if len(ref)==1:
                                            ref_codon = sequence[position-2]+ref+ sequence[position]
                                            var_codon = sequence[position - 2]+alt+sequence[position]
                                        elif len(ref)==2:
                                            ref_codon = sequence[position-2]+ref
                                            var_codon = sequence[position - 2]+alt
                                        else:
                                            ref_codon = sequence[position-2]+ref[0:2]
                                            var_codon = sequence[position - 2]+alt[0:2]
                                    else:
                                        if len(ref)==1:
                                            ref_codon = sequence[position-3]+sequence[position-2]+ref
                                            var_codon = sequence[position-3]+sequence[position-2]+alt
                                        elif len(ref)==2: 
                                            ref_codon = ref[1]+sequence[position]+sequence[position+1]
                                            var_codon = alt[1]+sequence[position]+sequence[position+1]  
                                        else:
                                            ref_codon = ref[1:]+sequence[position]
                                            var_codon = alt[1:]+sequence[position]                           
                                else:                        
                                    ex = len(ex_starts)-1
                                    while ex != 0 and position < int(ex_starts[ex]): 
                                        if cd_end < ex_starts[ex]:
                                            ex-=1
                                        elif cd_end < ex_ends[ex] and cd_end > ex_starts[ex] and position < ex_starts[ex]:
                                            m = ((cd_end - ex_starts[ex])+m)%3
                                            aalength=(cd_end-ex_starts[ex])
                                            ex-=1
                                        else:
                                            m = ((ex_ends[ex]-ex_starts[ex])+m)%3
                                            aalength+=(ex_ends[ex]-ex_starts[ex]) 
                                            ex-=1
                                    if cd_end > ex_ends[ex]:
                                        m=(ex_ends[ex] - position +m)%3
                                        aalength+=(ex_ends[ex]-position+1)
                                    else:
                                        m=(cd_end - position + m)%3
                                        aalength+=(cd_end-position+1)
                                    if aalength%3==0:
                                        aapos=aalength//3
                                    else:
                                        aapos=aalength//3 + 1
            
                                    toRev=1
                                    if m==0 :
                                        if len(ref)==1:
                                            ref_codon=ref+sequence[position-2]+sequence[position-3]
                                            var_codon=alt+sequence[position-2]+sequence[position-3]
                                        elif len(ref)==2:
                                            ref_codon=ref+sequence[position-2]
                                            var_codon=alt+sequence[position-2]
                                        else:
                                            ref_codon=ref
                                            var_codon=alt
                                    elif m==1:
                                        if len(ref)==1:
                                            ref_codon = sequence[position]+ref+ sequence[position-2]
                                            var_codon = sequence[position]+alt+sequence[position-2]
                                        elif len(ref)==1:
                                            ref_codon = sequence[position]+ref
                                            var_codon = sequence[position]+alt
                                        else:
                                            ref_codon = sequence[position]+ref[0:2]
                                            var_codon = sequence[position]+alt[0:2]
                                    else:
                                        if len(ref)==1:
                                            ref_codon = sequence[position+1]+sequence[position]+ref
                                            var_codon = sequence[position+1]+sequence[position]+alt
                                        elif len(ref)==1:
                                            ref_codon=ref[1]+sequence[position-2]+sequence[position-3]
                                            var_codon=alt[1]+sequence[position-2]+sequence[position-3]
                                        else:
                                            ref_codon=ref[1:]+sequence[position-2]
                                            var_codon=alt[1:]+sequence[position-2]
                                ref_codon=ref_codon.upper()
                                var_codon=var_codon.upper()
                                if toRev==1:
                                    intab = "ATCG"
                                    outtab = "TAGC"
                                    trantab = str.maketrans(intab, outtab)
                                    var_codon = var_codon.translate(trantab)
                                    ref_codon = ref_codon.translate(trantab)
                                    
                                var_codon = var_codon[0:3]
                                ref_codon = ref_codon[0:3]
                                ref_aa=''
                                var_aa=''
                                if var_codon in AADict:
                                    var_aa = AADict[var_codon] #translate
                                if ref_codon in AADict:
                                    ref_aa = AADict[ref_codon]
                                aachange.append(ref_aa+str(aapos)+var_aa)
                                if var_aa == ref_aa:
                                    change_type.append('Synonymous')
                                else:
                                    change_type.append('Non-synonymous')
                            else:  #To make sure to fill change type in on the ignored ones.
                                where.append('Ignore')
                                change_type.append('')
                        else:
                            where.append('Ignore')
                            change_type.append('')   
                dictionary[position]=[chrom,genename,where,change_type,aachange]   #include strand_gene for SIFT format 
    chr_end =  time.time() - chr_start
    del vcffile
    del seqfile
    gc.collect()    
    return dictionary
