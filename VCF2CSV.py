import vcf
import csv
import time
import timeit
import gc

def writecsv(inputvcf):    
    chr_start = time.time()
    vcf_reader = vcf.Reader(open(inputvcf,'r'))
    outfile_name=[]
    col_name=[['var_index','chrom','left','right','ref_seq','var_seq1','var_seq2','count1','count2','genotype','var_score','gene_name','where_in_transcript']]#,'ID'
    index=[]
    chrom=[]
    left=[]
    right=[]
    ref=[]
    alt=[]
    alt2=[]
    tags=[]
    names=locals()
    sample_names=[]
    sample_temps=[]
    base_pool=[['A'],['C'],['G'],['T']]

    record=next(vcf_reader)
    #dynamically names output file for each sample
    for i in range(0,len(record.samples)):
        names['%s' %''.join([char for char in record.samples[i].sample if char.isalnum()])]=[]
        names['%s' %''.join([char for char in record.samples[i].sample.lower() if char.isalnum()])]=[]
        names['count1%s' %i]=[]
        names['count2%s' %i]=[]
        names['GT%s' %i]=[]  ################################     new
        sample_temps.append(''.join([char for char in record.samples[i].sample if char.isalnum()]))
        sample_names.append(''.join([char for char in record.samples[i].sample.lower() if char.isalnum()]))
    vcf_reader = vcf.Reader(open(inputvcf,'r'))
    for record in vcf_reader:         
            chrom.append(record.CHROM.split())
            right.append(record.POS)
            left.append(str(record.POS).split())
            ref.append(record.REF.split())
            alt.append(record.ALT)
            tags.append(record.FORMAT.split(':'))
            for j in range(0,len(record.samples)):
                eval(sample_temps[j]).append(record.samples[j])    
    for m in range(0,len(tags)):
        for k in range(0,len(record.samples)):
            names['temp%s' %k]=[]
        for n in range(0,len(tags[0])):
            if(tags[0][n] not in tags[m][n]):
                for j in range(0,len(record.samples)):
                    eval('temp'+str(j)).append(None)
                tags[m].insert(n,tags[0][n])
            else:
                for j in range(0,len(record.samples)):
                    eval('temp'+str(j)).append(eval(sample_temps[j])[m][tags[0][n]])
        for j in range(0,len(record.samples)):
            eval(sample_names[j]).append(eval('temp'+str(j)))
    def get_index(input_list):
        AD=-1
        BCOUNT=-1                         
        for j in range(0,len(input_list)):
            if (input_list[j]=='AD'):
                AD=j
                return AD
            elif (input_list[j]=='BCOUNT'):
                BCOUNT=j
                return BCOUNT
            else:
                continue
        return
    index_count=get_index(tags[0])                    
    if (len(chrom)==len(left)==len(ref)==len(alt)):
            for j in range(0,len(record.samples)):
                names['list0%s' %j]=[[]]*len(chrom)
                names['temp%s' %j]=[[]]*len(chrom)
            index=[num for num in range(1,(len(chrom)+1))]
            score=[['28']]*len(chrom)
            alt2=[['']]*len(chrom)
            for i in range(0,len(chrom)):
                    if(len(alt[i])>1):
                        alt2[i]=[str(alt[i][1])]
                        alt[i]=[str(alt[i][0])]
                    else:
                        if(ref[i][0]<str(alt[i][0])):
                            alt2[i]=[str(alt[i][0])]
                            alt[i]=ref[i]
                        else:
                            alt2[i]=ref[i]
                            alt[i]=[str(alt[i][0])]
                    
                    index[i]=[str(index[i])]
                    right[i]=[str(right[i]+len(ref[i]))]

                    for j in range(0,len(record.samples)):  ##############################################   This is new
                        eval('GT'+str(j)).append([",".join(eval(sample_names[j])[i][0].split('/'))])
                    
                    if('BCOUNT' in tags[0]):
                        pos1=base_pool.index(alt[i])
                        pos2=base_pool.index(alt2[i])
                        for k in range(0,len(record.samples)):
                            eval('count1'+str(k)).append([eval(sample_names[k])[i][index_count][pos1]])
                            eval('count2'+str(k)).append([eval(sample_names[k])[i][index_count][pos2]])
                    elif('AD' in tags[0]):
                        if(alt[i]==ref[i]):
                            for k in range(0,len(record.samples)):
                                eval('count1'+str(k)).append([eval(sample_names[k])[i][index_count][0]])
                                eval('count2'+str(k)).append([eval(sample_names[k])[i][index_count][1]])    
                        elif(alt2[i]==ref[i]):
                            for k in range(0,len(record.samples)):
                                eval('count1'+str(k)).append([eval(sample_names[k])[i][index_count][1]])
                                eval('count2'+str(k)).append([eval(sample_names[k])[i][index_count][0]])
                        else:
                            for k in range(0,len(record.samples)):
                                eval('count1'+str(k)).append([eval(sample_names[k])[i][index_count][0]])
                                eval('count2'+str(k)).append([eval(sample_names[k])[i][index_count][1:len(eval(sample_names[k])[i][index_count])]])
                    else:
                        break                   
                    for j in range(0,len(record.samples)):
                        eval('list0'+str(j))[i]= index[i]+chrom[i]+ left[i]+ right[i]+ ref[i]+ alt[i]+alt2[i]+eval('count1'+str(j))[i]+eval('count2'+str(j))[i]+eval('GT'+str(j))[i]+score[i]#+normal[i]+ID[i]
                        #eval('temp'+str(j))[i]=eval('list0'+str(j))[i]
                        eval('temp'+str(j))[i]= index[i]+chrom[i]+ left[i]+ right[i]+ ref[i]+ alt[i]+alt2[i]+eval('count1'+str(j))[i]+eval('count2'+str(j))[i]+eval('GT'+str(j))[i]+score[i]

    for p in range(0,len(record.samples)):  
        for ele in eval('temp'+str(p)):
            if (type(ele[8])==int):
                if ((ele[7]+ele[8])<=5 or ele[7]/(ele[7]+ele[8])<=0.05 or ele[8]/(ele[7]+ele[8])<=0.05):
                    eval('list0'+str(p)).remove(ele)
            else:
                pos=eval('list0'+str(p)).index(ele)
                eval('list0'+str(p)).insert(pos+1,ele)
                if (ele[4]<ele[5]):
                    eval('list0'+str(p))[pos][8]=ele[8][0]
                    eval('list0'+str(p))[pos][5]=ele[4]
                    eval('list0'+str(p))[pos][6]=ele[5]
                else:
                    eval('list0'+str(p))[pos][7]=ele[8][0]
                    eval('list0'+str(p))[pos][5]=ele[5]
                    eval('list0'+str(p))[pos][6]=ele[4]
                    eval('list0'+str(p))[pos][8]=ele[7]
                if (ele[4]<ele[6]):
                    eval('list0'+str(p))[pos+1][8]=ele[8][1]
                    eval('list0'+str(p))[pos+1][5]=ele[4]
                    eval('list0'+str(p))[pos+1][6]=ele[6]
                else:
                    eval('list0'+str(p))[pos+1][7]=ele[8][1]
                    eval('list0'+str(p))[pos+1][5]=ele[6]
                    eval('list0'+str(p))[pos+1][6]=ele[4]
                    eval('list0'+str(p))[pos+1][8]=ele[7]
                if ((eval('list0'+str(p))[pos+1][7]+eval('list0'+str(p))[pos+1][8])<=5 or eval('list0'+str(p))[pos+1][7]/(eval('list0'+str(p))[pos+1][7]+eval('list0'+str(p))[pos+1][8])<=0.05
                    or eval('list0'+str(p))[pos+1][8]/(eval('list0'+str(p))[pos+1][7]+eval('list0'+str(p))[pos+1][7])<=0.05):
                    eval('list0'+str(p)).pop(pos+1)
                if ((eval('list0'+str(p))[pos][7]+eval('list0'+str(p))[pos][8])<=5 or eval('list0'+str(p))[pos][7]/(eval('list0'+str(p))[pos][7]+eval('list0'+str(p))[pos][8])<=0.05
                    or eval('list0'+str(p))[pos][8]/(eval('list0'+str(p))[pos][7]+eval('list0'+str(p))[pos][7])<=0.05):
                    eval('list0'+str(p)).pop(pos)
        index1=[num for num in range(1,(len(eval('list0'+str(p)))+1))]
        for q in range(len(eval('list0'+str(p)))):
            eval('list0'+str(p))[q][0]=index1[q]     
    for h in range(0,len(record.samples)):
        names['file1%s' %h]=open(sample_names[h]+'.csv',"a",newline='')
        outfile_name.append(sample_names[h]+'.csv')       
        with eval('file1'+str(h)):
                writer=csv.writer(eval('file1'+str(h)))
                writer.writerows(col_name)
                writer.writerows(eval('list0'+str(h)))                
    for h in range(0,len(record.samples)):
        eval('file1'+str(h)).close()
    end = time.time()-chr_start
    del index
    del chrom
    del ref
    del alt
    del alt2
    for key in list(locals()):
        if (key.startswith('list0') or key.startswith('temp')):
            del locals()[key]
    del vcf_reader
    gc.collect()
    return outfile_name

#if __name__=='__main__':   
#    writecsv(“ffbe667a-e317-4fd7-bba8-3caf626d7f14.vcf”): 
        
        
        
        
        
        
