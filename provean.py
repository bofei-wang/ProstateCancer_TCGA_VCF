# -*- coding: utf-8 -*-
"""
Created on Tue Aug 31 22:49:21 2021

@author: Bofei
"""

import re
from miscVar import AADict
MAX_DIST = 5000
import bisect
import pandas as pd
import os

def get_same_element_index(ob_list, word):
    return [i for (i, v) in enumerate(ob_list) if v == word]

intab = "ATCG"
outtab = "TAGC"
trantab = str.maketrans(intab, outtab)

def annotation(unify,n=3):
    dictionary= {} 
    position = 0
    chrom=''
    prechrom=0
    
    for l in range(unify.shape[0]): #line in vcf:
        chrom=unify["chrom"][l]
        position = unify["left"][l]
        ref=unify['ref_seq'][l]
        alt=unify['alt'][l]
        
        gtfCounter = -1 #point to each line in gtf file
        genestart = 0
        genestop = 0
        #if position not in dictionary:
        boo = 0
        where = []
        genename = []
        change_type=[]
        counter=[]
        strand_gene=[]
        aavar=[]
        dnasequence=[]
        codon_pos=[]
        ref_trip=[]
        var_trip=[]
        if chrom != prechrom:
            f = open(r"D:\Personal materials\PhD\Oncominer\OP_program\ref_genome\\"+str(chrom)+'_hg38refFlat_2020.txt','r')
            flat = f.readlines()
            f.close()
            seqfile=open(r"D:\Personal materials\PhD\Oncominer\OP_program\ref_genome\genome_"+str(chrom)+".fa",'r')
            next(seqfile)
            sequence=seqfile.read().replace('\n','')
            seqfile.close()
        prechrom=chrom
        
        while boo == 0 and gtfCounter < len(flat):
            if position >= genestart and position <= genestop and gtfCounter < (len(flat)-1): #GSV within gene
                counter.append(gtfCounter)
                strand_gene.append(flat[gtfCounter].split()[3])
                
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
       
        for i in range(len(counter)):
            genename.append(flat[counter[i]].split()[0]+'_'+str(counter[i]))
            #tr_start = int(flat[counter[i]].split()[4]) #transcript  start; check
            #tr_end = int(flat[counter[i]].split()[5]) #transcript  end; check
            cd_start=int(flat[counter[i]].split()[6]) #CDS start
            cd_end=int(flat[counter[i]].split()[7]) #CDS end
            ex_starts = []
            ex_ends=[]
            dnaseq=''
            x=0
            tmp = 'Unknown'
            ex_starts = [int(i) for i in flat[counter[i]].split()[-2].split(',')[0:-1]]
            ex_ends = [int(i) for i in flat[counter[i]].split()[-1].split(',')[0:-1]]

            while tmp == 'Unknown' and x < len(ex_starts):
                if position > ex_starts[x] and position <= ex_ends[x]: #if position in exon
                    tmp = 'Exon'
                elif position < ex_starts[x]:
                    tmp = 'Intron'
                else:
                    x+=1
            if tmp == 'Intron' or tmp=='Unknown':
                where.append('Intron')
                change_type.append('')
            #see if position is in translated region:
            elif tmp == 'Exon':
                if position >  cd_start and position <= cd_end:
                    tmp = 'Translated'
                else: #not translated
                    tmp = 'Not translated'
                if tmp == 'Not translated':
                    if cd_start == cd_end: #ncRNAs
                        where.append('Not translated, ncRNA')
                        change_type.append('')
                    else: #utrs
                        if position  < cd_start:
                            if strand_gene[i] == '-':
                                where.append('Not translated,3\' utr')
                                change_type.append('')
                            else:
                                where.append('Not translated,5\' utr')
                                change_type.append('')
                        else:
                            if strand_gene[i] == '-':
                                where.append('Not translated,5\' utr')
                                change_type.append('')
                            else:
                                where.append('Not translated,3\' utr')
                                change_type.append('')
                elif tmp == 'Translated':
                    where.append('CDS')
                    ref_codon=''
                    var_codon=''
                    m=0
                    toRev=0
                    aapos=0
                    aaseq=''
                    #m = (position - cd_start) % 3
                    if strand_gene[i]=="+":
                        ex=0
                        #pre_ex=0
                        while ex != (len(ex_ends)-1) and position > ex_ends[ex]:
                            if cd_start > ex_ends[ex]:
                                ex+=1
                            elif cd_start > ex_starts[ex] and cd_start < ex_ends[ex] and position > ex_ends[ex]:
                                #pre_ex=ex
                                m = ((ex_ends[ex] - cd_start)+m)%3
                                aapos+=(ex_ends[ex] - cd_start)
                                ex+=1
                            else:
                               # pre_ex=ex
                                m = ((ex_ends[ex]-ex_starts[ex])+m)%3
                                aapos+=(ex_ends[ex]-ex_starts[ex])
                                ex+=1
                        #print("ex is %d" %ex)
                        if cd_start > ex_starts[ex]:
                            m=(position-cd_start -1+m )%3
                            aapos+=(position-cd_start) 
                        else:
                            m=(position-ex_starts[ex] -1+m)%3
                            aapos+=(position-ex_starts[ex]) 
                        if aapos%3==0:
                            pos=aapos//3
                        else:
                            pos=aapos//3 + 1
                        firstex=bisect.bisect(ex_starts, cd_start)-1
                        lastex=bisect.bisect(ex_ends,cd_end)
                        if lastex>=len(ex_ends):
                            lastex=len(ex_ends)-1
                        if ex_ends[firstex]>cd_start: 
                            dnaseq = dnaseq+sequence[cd_start:ex_ends[firstex]]
                        firstex+=1
                        while firstex < lastex:
                            dnaseq=dnaseq+sequence[ex_starts[firstex]:ex_ends[firstex]]
                            firstex+=1
                        if cd_end > ex_starts[lastex]:
                            dnaseq=dnaseq+sequence[ex_starts[lastex]:cd_end]
#                            aapos=(ex_ends[firstex]-cd_start)/3
#                            for idx in range(firstex+1,ex):
#                                aapos+=(ex_ends[idx]-ex_starts[idx])/3 
#                            aapos+=(position-ex_starts[ex])//3 + 1 
#                            aapos=int(aapos)
#                            aapos=(position-cd_start)//3 + 1                               
                            
                    #if len(ref)==len(alt)and position<len(sequence)-2: 
                        if m==0 :
                            ref_codon=ref+sequence[position]+sequence[position+1]
                            var_codon=alt+sequence[position]+sequence[position+1]
                            codon_pos.append(1)
                        elif m==1:
                            ref_codon = sequence[position-2]+ref+ sequence[position]
                            var_codon = sequence[position - 2]+alt+sequence[position]
                            codon_pos.append(2)
                        else:
                            ref_codon = sequence[position-3]+sequence[position-2]+ref
                            var_codon = sequence[position-3]+sequence[position-2]+alt  
                            codon_pos.append(3)                           
                    else:                        
                        ex = len(ex_starts)-1
                        while ex != 0 and position < int(ex_starts[ex]): 
                            if cd_end < ex_starts[ex]:
                                ex-=1
                            elif cd_end < ex_ends[ex] and cd_end > ex_starts[ex] and position < ex_starts[ex]:
                                m = ((cd_end - ex_starts[ex] )+m)%3
                                aapos=(cd_end-ex_starts[ex])
                                ex-=1
                            else:
                                m = ((ex_ends[ex]-ex_starts[ex] )+m)%3
                                aapos+=(ex_ends[ex]-ex_starts[ex]) 
                                ex-=1
                        if cd_end > ex_ends[ex]:
                            m=(ex_ends[ex] - position + m)%3
                            aapos+=(ex_ends[ex]-position)                          
                        else:
                            m=(cd_end - position + m)%3
                            aapos+=(cd_end-position)
                            
                        if aapos%3==0:
                            pos=aapos//3
                        else:
                            pos=aapos//3 + 1
                        firstex=bisect.bisect(ex_starts, cd_start)-1
                        lastex=bisect.bisect(ex_ends,cd_end)
                        if lastex>=len(ex_ends):
                            lastex=len(ex_ends)-1
                        if ex_ends[firstex]>cd_start: 
                            dnaseq = dnaseq+sequence[cd_start:ex_ends[firstex]]
                        firstex+=1
                        while firstex < lastex:
                            dnaseq=dnaseq+sequence[ex_starts[firstex]:ex_ends[firstex]]
                            firstex+=1
                        if cd_end > ex_starts[lastex]:
                            dnaseq=dnaseq+sequence[ex_starts[lastex]:cd_end]
                        dnaseq=dnaseq[::-1]    
                        
                        toRev=1
                        if m==0 :
                            ref_codon=ref+sequence[position-2]+sequence[position-3]
                            var_codon=alt+sequence[position-2]+sequence[position-3]
                            codon_pos.append(1)
                        elif m==1:
                            ref_codon = sequence[position]+ref+ sequence[position-2]
                            var_codon = sequence[position]+alt+sequence[position-2]
                            codon_pos.append(2)
                        else:
                            ref_codon = sequence[position+1]+sequence[position]+ref
                            var_codon = sequence[position+1]+sequence[position]+alt
                            codon_pos.append(3)
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
                    if var_aa == ref_aa:
                        change_type.append('Synonymous')
                    else:
                        change_type.append('Non-synonymous')
                    aavar.append(ref_aa+str(int(pos))+var_aa)
                    ref_trip.append(ref_codon)
                    var_trip.append(var_codon)
                    dnasequence.append(dnaseq)
#                    for k in range(0,len(dnaseq),3):
#                        aaseq=aaseq+AADict[dnaseq[k:k+3].upper()]
#                    f=open(path+genename[-1]+".fasta","w")
#                    f.write(aaseq)
#                    f.close()    
                else:  #To make sure to fill change type in on the ignored ones.
                    where.append('Ignore')
                    change_type.append('')
            else:
                where.append('Ignore')
                change_type.append('')
        dictionary[l]=[chrom,position,strand_gene,genename,where,change_type,dnasequence] 
    return dictionary   

if __name__=='__main__': 
    nonsynpath=input("nonsynonymous variants file:")
    aaseqpath=input("directory to store amino acid sequence:")
    aavarpath=input("directory to store amino acid variation:")
    
    nonsyn=pd.read_csv(nonsynpath)
    seqname=[]
    result=annotation(nonsyn)
    for j in range(len(result)):
        nonsynidx=get_same_element_index(result[j][5],'Non-synonymous')
        for n in range(len(nonsynidx)):
            varfile=open(aavarpath+result[j][3][nonsynidx[n]]+".var", "a")     
            #varfile.write(result[j][6][n]+"\n")
            varfile.write(nonsyn['AAchange'][j][nonsynidx[n]])
            varfile.close() 
            if result[j][3][nonsynidx[n]] not in seqname:
                seqname.append(result[j][3][nonsynidx[n]])
                if result[j][2][nonsynidx[n]]=='-':               
                    geneseq=result[j][6][n].upper().translate(trantab)   
                else:
                    geneseq=result[j][6][n].upper()
                aaseq=''
                for i in range(0,len(geneseq),3):
                    if len(geneseq)>=i+3 and AADict[geneseq[i:i+3]]!='*':
                        aaseq=aaseq+AADict[geneseq[i:i+3]]
                    else:
                        break
                aaseq=re.sub("(.{60})","\\1\n",aaseq)
                output=open(aaseqpath+result[j][3][nonsynidx[n]]+".fasta", "w")
                output.write("> sp|" + result[j][3][nonsynidx[n]] + "\n")
                output.write(aaseq)
                output.close()
       
    output=open(aavarpath+"name.txt",'w')
    for ele in seqname:
        output.write(ele +"\n")    
    output.close()        
    #### modify the last character of sequence to make the stop codon as *                   
    for i in range(len(seqname)):        
        var=pd.read_csv(aavarpath+seqname[i]+".var",header=None)
        f=open(aaseqpath+seqname[i]+".fasta",'r')         
        sequence=f.readlines()
        f.close()  
        if sequence[-1][-1]=="+":
            sequence[-1]=sequence[-1][0:-1]+'*'
        elif sequence[-1][-1]=="*":
            continue
        else:
            sequence[-1]=sequence[-1]+'*'
    
        for k in range(var.shape[0]):
            loc=int(var.iloc[k,0][1:-1])
            if sequence[loc-1]==var.iloc[k,0][0]:
                continue
            elif sequence[loc-2]==var.iloc[k,0][0]:
                var.iloc[k,0]=var.iloc[k,0][0]+str(loc-1)+var.iloc[k,0][-1]
            elif sequence[loc]==var.iloc[k,0][0]:
                var.iloc[k,0]=var.iloc[k,0][0]+str(loc+1)+var.iloc[k,0][-1]
            else:
                print(seqname[i])
        var.to_csv(aavarpath+seqname[i]+".var",index=False,header=False)    
        sequence=re.sub("(.{60})","\\1\n",sequence)
        fw=open(aaseqpath+seqname[i]+".fasta",'w') 
        fw.write("> sp|" + seqname[i] + "\n")
        fw.write(sequence)
        fw.close()   

        
        
        
        
        