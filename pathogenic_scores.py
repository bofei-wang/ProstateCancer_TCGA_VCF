# -*- coding: utf-8 -*-
"""
Created on Fri Oct  8 20:27:07 2021

@author: Bofei
"""

import pandas as pd
from collections import Counter
import math
import argparse

def get_same_element_index(ob_list, word):
    return [i for (i, v) in enumerate(ob_list) if v == word]

def collect_gene(var_df):
    singlegene=[]
    provean_min=[]
    fat_combine=[]
    for i in range(var_df.shape[0]):
        word=var_df['genename'][i].replace('\'', '').replace('[', '').replace(']', '').replace(' ', '').split(',')
        anno=var_df['where'][i].replace('\'', '').replace('[', '').replace(']', '').replace(' ', '').split(',')
        anno=list(filter(lambda x:x!='"Nottranslated' and x!='Nottranslated',anno))
        idx=get_same_element_index(anno,'CDS')
        word1=[]
        for ele in idx:
            word1.append(word[ele])
        singlegene.append(Counter(word1).most_common(1)[0][0])
        provean=var_df['provean score'][i].split(",")
        provean=list(map(lambda x:float(x),provean))
        provean_min.append(min(provean))
        if str(var_df['Coding Score'][4])=='nan':
            fat_combine.append(float(var_df['Non-Coding Score'][i]))
        else:
            fat_combine.append(float(var_df['Coding Score'][i]))
            
    var_df.insert(var_df.shape[1],'single gene',singlegene) 
    var_df.insert(var_df.shape[1],'provean min',provean_min)
    var_df.insert(var_df.shape[1],'fat_combine',fat_combine)
    var_df.insert(var_df.shape[1],'provean rank',var_df['provean min'].rank())   
    var_df.insert(var_df.shape[1],'fathmm rank',var_df['fat_combine'].rank(ascending=False))
    
    var_df['rank sum']=var_df['fathmm rank']+var_df['provean rank']
    var_df['rank-score']=1-var_df['rank sum']/(2*var_df.shape[0])
    var_df=var_df.sort_values('rank sum',axis=0)
    return var_df

def select_cds(curr):
    maxium=0
    for i in range(curr.shape[0]):
        ex=0   
        m=0                  
        while ex != (len(curr.iloc[i,10])):
            if curr.iloc[i,6] > curr.iloc[i,10][ex]:
                ex+=1
            elif curr.iloc[i,6] > curr.iloc[i,9][ex] and curr.iloc[i,6] < curr.iloc[i,10][ex]:
                m = (curr.iloc[i,10][ex] - curr.iloc[i,6])+m
                ex+=1
            else:
                m = (curr.iloc[i,10][ex]-curr.iloc[i,9][ex])+m
                ex+=1
        if m> maxium:
            maxium=m
    return m

def pathogenic_score(var_df,gtf,out):
    gene_counts=var_df['single gene'].value_counts()        
    patho_gene=gene_counts.reset_index()
    patho_gene.rename(columns={'index':'genename','single gene':'distinct pathogenic nsGSV'},inplace=True)
    cds_length=[]
    for i in range(patho_gene.shape[0]):
        curr=gtf[gtf[0]==patho_gene.iloc[i,0]]
        cds_length.append(select_cds(curr))
    patho_gene.insert(patho_gene.shape[1],'CDS length',cds_length)
    plg=[]
    ptg=[]
    for i in range(patho_gene.shape[0]):
        curr=var_df[var_df['single gene']==patho_gene.iloc[i,0]]
        plg.append(sum((curr['SubT1N0']-curr['SubT0N1'])*curr['rank-score'])/math.log(patho_gene.iloc[i,2]))
        #curr=curr[(curr['SubT1N0']+curr['SubT1N1']>0) & (curr['SubT0N1']+curr['SubT1N1']==0)]
        ptg.append(sum((curr['SubT1N0']-curr['SubT0N1'])*curr['rank-score']))
     
    patho_gene.insert(patho_gene.shape[1],'Pt(g)',ptg)    
    patho_gene.insert(patho_gene.shape[1],'Pl(g)',plg) 
    patho_gene.insert(patho_gene.shape[1],'Pt(g) rank',patho_gene['Pt(g)'].rank(ascending=False))    
    patho_gene.insert(patho_gene.shape[1],'Pl(g) rank',patho_gene['Pl(g)'].rank(ascending=False))  
    patho_gene.to_csv(out,index=False)
   
if __name__=='__main__':    
    parser = argparse.ArgumentParser(description='count subjects for each variant')
    parser.add_argument('-t', '--input_dir', metavar="nsGSV", dest='nonsyn', required=True, help='Required: nsGSV files.')
    parser.add_argument('-f', '--ref_flat', metavar="REFFLAT", dest='ref_flat', required=True, help='Required: Reference genomic field file.')
    parser.add_argument('-o', '--out_dir', metavar="DIRNAME", dest='outfile', required=True, help='Required: Path and name of output file.')
    args = parser.parse_args()

    var_df=pd.read_csv(args.nonsyn)
    gtf = pd.read_table(args.ref_flat)
    outdir = args.outdir
    var_df=collect_gene(var_df)
    pathogenic_score(var_df,gtf,outdir)
    
    
    
    
    
    
    
    
    
    
    
    
    
    