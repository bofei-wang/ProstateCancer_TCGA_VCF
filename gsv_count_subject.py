# -*- coding: utf-8 -*-
"""
Created on Wed May 19 09:48:00 2021

@author: Bofei
"""
import numpy as np
import os
import pandas as pd
import argparse

def subject_count(tumordir,normaldir,outdir):
    #tumordir="D:\\Personal materials\\PhD\\Oncominer\\Data\\PRODOMI\\OMI_PROD"
    #normaldir="D:\\Personal materials\\PhD\\Oncominer\\Data\\PRODOMI\\OMI_NORMAL"
    tumors=os.listdir(tumordir)
    normals=os.listdir(normaldir)
    tumor=[value for value in tumors if value.endswith('tumor.csv')]
    normal=[value for value in normals if value.endswith('normal.csv')]
    tumorname=[value.split("_")[0] for value in tumor]
    normalname=[value.split("_")[0] for value in normal]
    del tumors,normals,tumor,normal
    
    tumorvar=[]
    normalvar=[]
    commonvar=[]
    counttumor=[]
    countnormal=[]
    countcommon=[]
    numtum=0
    numnor=0
    numcom=0
    
    if(normalname==tumorname):
        disease_var=pd.DataFrame()
        normal_var=pd.DataFrame()
        common_var=pd.DataFrame()
        for name in tumorname:
            tf=pd.read_csv(tumordir+name+"_tumor.csv")
            nf=pd.read_csv(normaldir+name+"_normal.csv")
            norepeat_tf = tf.drop_duplicates(subset=['chrom', 'left','right','ref_seq','var_seq1',
                                                     'var_seq2'], keep='first')
            norepeat_nf = nf.drop_duplicates(subset=['chrom', 'left','right','ref_seq','var_seq1',
                                                     'var_seq2'], keep='first')
        
            tidx=[]
            nidx=[]
            for i in range(norepeat_tf.shape[0]):
                for j in range(norepeat_nf.shape[0]):
                    if (norepeat_tf.iloc[i,1:7].values.tolist() == norepeat_nf.iloc[j,1:7].values.tolist()):
                        tidx.append(i)
                        nidx.append(j)
                        break
            
            a=pd.DataFrame()
            norepeat_tf=a.append(norepeat_tf,ignore_index=True)
            norepeat_nf=a.append(norepeat_nf,ignore_index=True)
            
            common=a.append(norepeat_tf.loc[tidx],ignore_index=True)
            tempcom=a.append(norepeat_tf.loc[tidx],ignore_index=True)
            norepeat_tf=norepeat_tf.drop(tidx,axis=0) 
            norepeat_nf=norepeat_nf.drop(nidx,axis=0) 
            temptf=a.append(norepeat_tf,ignore_index=True)
            tempnf=a.append(norepeat_nf,ignore_index=True)
            
            
            for j in range(common.shape[0]):
                if (list(common.iloc[j,1:7]) in commonvar):
                    countcommon[commonvar.index(list(common.iloc[j,1:7]))]+=1
                    tempcom=tempcom.drop(j,axis=0)
                else:
                    commonvar.append(list(common.iloc[j,1:7]))
                    countcommon.append(1)
            numcom=numcom+tempcom.shape[0]
            common_var=common_var.append(tempcom,ignore_index=True)
    
            for j in range(norepeat_tf.shape[0]):
                if (list(norepeat_tf.iloc[j,1:7]) in tumorvar):
                    counttumor[tumorvar.index(list(norepeat_tf.iloc[j,1:7]))]+=1
                    temptf=temptf.drop(j,axis=0)
                else:
                    tumorvar.append(list(norepeat_tf.iloc[j,1:7]))
                    counttumor.append(1)
            numtum=numtum+temptf.shape[0]
            disease_var=disease_var.append(temptf,ignore_index=True)
            
            
            for j in range(norepeat_nf.shape[0]):
                if (list(norepeat_nf.iloc[j,1:7]) in normalvar):
                    countnormal[normalvar.index(list(norepeat_nf.iloc[j,1:7]))]+=1
                    tempnf=tempnf.drop(j,axis=0)
                else:
                    normalvar.append(list(norepeat_nf.iloc[j,1:7]))
                    countnormal.append(1)
            numnor=numnor+tempnf.shape[0]
            normal_var=normal_var.append(tempnf,ignore_index=True)
    ### keep count in each group
    disease_var.insert(disease_var.shape[1],"count",counttumor)
    normal_var.insert(normal_var.shape[1],"count",countnormal)
    common_var.insert(common_var.shape[1],"count",countcommon)
    transcript=disease_var.append(normal_var)
    transcript=transcript.append(common_var)
    transcript=transcript.drop_duplicates(subset=['chrom', 'left','right','ref_seq','var_seq1',
                                                     'var_seq2'], keep='first')
    
    t0n0=[0]*transcript.shape[0]
    t0n1=[0]*transcript.shape[0]
    t1n0=[0]*transcript.shape[0]
    t1n1=[0]*transcript.shape[0]
    for i in range(transcript.shape[0]):
        var=transcript.iloc[i,0:7].tolist()
        for name in tumorname:
            tf=pd.read_csv(tumordir+name+"_tumor.csv")
            nf=pd.read_csv(normaldir+name+"_normal.csv")
            norepeat_tf = tf.drop_duplicates(subset=['chrom', 'left','right','ref_seq','var_seq1',
                                                     'var_seq2'], keep='first')
            norepeat_nf = nf.drop_duplicates(subset=['chrom', 'left','right','ref_seq','var_seq1',
                                                     'var_seq2'], keep='first')
            tfile=np.array(norepeat_tf.iloc[:,1:7]).tolist()
            nfile=np.array(norepeat_nf.iloc[:,1:7]).tolist()
            if var in tfile and var in nfile:
                t1n1[i] += 1
            elif var in tfile and var not in nfile:
                t1n0[i] += 1
            elif var not in tfile and var in nfile:    
                t0n1[i]+=1
            else:
                t0n0[i]+=1
    transcript.insert(transcript.shape[1],"subt1n1",t1n1)
    transcript.insert(transcript.shape[1],"subt1n0",t1n0)
    transcript.insert(transcript.shape[1],"subt0n1",t0n1)
    transcript.insert(transcript.shape[1],"subt0n0",t0n0)
    
    alt=[]
    for k in range(transcript.shape[0]):
        if transcript.iloc[k,4]==transcript.iloc[k,5]:
            alt.append(transcript.iloc[k,6])
        else:
            alt.append(transcript.iloc[k,5])
    transcript.insert(5,'alt',alt)
    transcript.drop(columns=['var_index','var_seq1', 'var_seq2'],inplace=True)    
    
        
    '''
    transcript variants distribution combining tumor, normal and common
    '''                
    file_total_tumor=[0]*transcript.shape[0]
    file_total_normal=[0]*transcript.shape[0]
    for name in tumorname:
        tf=pd.read_csv(tumordir+name+"_tumor.csv")
        nf=pd.read_csv(normaldir+name+"_normal.csv")
        norepeat_tf = tf.drop_duplicates(subset=['chrom', 'left','right','ref_seq','var_seq1',
                                                'var_seq2'], keep='first')
        norepeat_nf = nf.drop_duplicates(subset=['chrom', 'left','right','ref_seq','var_seq1',
                                                'var_seq2'], keep='first')
        
        alt=[]
        for k in range(norepeat_tf.shape[0]):
            if norepeat_tf.iloc[k,4]==norepeat_tf.iloc[k,5]:
                alt.append(norepeat_tf.iloc[k,6])
            else:
                alt.append(norepeat_tf.iloc[k,5])
        norepeat_tf.insert(5,'alt',alt)
        norepeat_tf.drop(columns=['var_seq1', 'var_seq2'],inplace=True)
        
        alt=[]
        for k in range(norepeat_nf.shape[0]):
            if norepeat_nf.iloc[k,4]==norepeat_nf.iloc[k,5]:
                alt.append(norepeat_nf.iloc[k,6])
            else:
                alt.append(norepeat_nf.iloc[k,5])
        norepeat_nf.insert(5,'alt',alt)
        norepeat_nf.drop(columns=['var_seq1', 'var_seq2'],inplace=True)
        
        per_file_tumor=[0]*transcript.shape[0]
        per_file_normal=[0]*transcript.shape[0]
        
        merge_tumor=pd.merge(transcript,norepeat_tf,on=['chrom', 'left', 'right', 'ref_seq', 'alt'],how='left')
        include_tumor=merge_tumor[~merge_tumor['count1_y'].isnull()]
        index_tumor=include_tumor._stat_axis.tolist()
        for ele in index_tumor:
            per_file_tumor[ele]=1
            file_total_tumor[ele]+=1
        merge_normal=pd.merge(transcript,norepeat_nf,on=['chrom', 'left', 'right', 'ref_seq', 'alt'],how='left')
        include_normal=merge_normal[~merge_normal['count1_y'].isnull()]
        index_normal=include_normal._stat_axis.tolist()
        for ele in index_normal:
            per_file_normal[ele]=1
            file_total_normal[ele]+=1
    
    
        transcript.insert(transcript.shape[1],name+'_tumor',per_file_tumor)
        transcript.insert(transcript.shape[1],name+'_normal',per_file_normal)
    
    #return transcript
    transcript.to_csv(outdir+"sub_count_dist.csv",index=False)


if __name__=='__main__':    
    #tumordir=input("directory of tumor OMI file:")
    #normaldir=input("directory of normal OMI file:")
    #outdir=input("directory of output file:")
    parser = argparse.ArgumentParser(description='count subjects for each variant')
    parser.add_argument('-t', '--input_dir', metavar="TUMORDIR", dest='tumor', required=True, help='Required: Directory to all tumor OMI files.')
    parser.add_argument('-o', '--out_dir', metavar="DIRNAME", dest='outdir', required=True, help='Required: Name of output directory to store result file.')
    parser.add_argument('-n', '--input_dir', metavar="NORMALDIR", dest='normal', required=True, help='Required: Directory to all normal OMI files.')
    args = parser.parse_args()

    tumordir = args.tumor
    normaldir = args.normal
    outdir = args.outdir
    subject_count(tumordir,normaldir,outdir)
    
