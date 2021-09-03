from parse_VCF import vcf_func
from VCF2CSV import writecsv
#from functools import partial
import csv
import vcf
import gc
import time
#import sys
import os
import xlwt
from xlutils.copy import copy
from xlrd import open_workbook
import multiprocessing
import argparse
#import resource


#Function combines all information and writes an OMI file
def writeOMI(infile,outfile,newfile):    
    start=time.time()
    f1=open(outfile,'r')
    rows=csv.reader(f1)
    f2=open(newfile,'a',newline='')
    writer=csv.writer(f2)
    for row in rows:
        for key in infile:
            if (row[1]==infile[key][0] and row[2]==str(key)):               
                row.append(infile[key][1])
                row.append(infile[key][2])
                row.append(infile[key][3])
                row.append(infile[key][4])
                writer.writerow(row)
                break                                        
            
    f1.close()
    f2.close()
    end=time.time()-start
    del rows
    gc.collect()

#Function that splits VCF file by chromosome for parallel processing
def splitvcf(vcffile):
    start=time.time()
    for i in range(1,23):
        vcf_reader = vcf.Reader(open('%s' %(vcffile), 'r'))
        outvcf=open('chr'+str(i)+'.vcf','w')#names['chr%s' %i]=[]
        my_out=vcf.Writer(outvcf,vcf_reader)
        for record in vcf_reader:
            if (record.CHROM== 'chr'+str(i)):
                my_out.write_record(record)
        outvcf.close()
    vcf_reader = vcf.Reader(open('%s' %(vcffile), 'r'))   
    outvcf1=open('chrX.vcf','w')
    outvcf2=open('chrY.vcf','w')
    my_out1=vcf.Writer(outvcf1,vcf_reader)
    my_out2=vcf.Writer(outvcf2,vcf_reader)

    for record in vcf_reader:
        if (record.CHROM== 'chrX'):
            my_out1.write_record(record)
        elif(record.CHROM== 'chrY'):
            my_out2.write_record(record)
        else:
            continue
    outvcf1.close()
    outvcf2.close()
    end=time.time()-start
    del vcf_reader
    gc.collect()
    print('Time for generating vcf files:'+str(round(end))+'s')
    return end

#Function that splits reference file by chromosome for parallel processing
def splitgtf(gtfdir):
    start=time.time()
    filegtf=open(gtfdir,'r')
    flat=filegtf.readlines()
    filegtf.close()
    resultx=[]
    resulty=[]
    for i in range(1,23):
        myout=open('chr'+str(i)+'_ref.txt','a')
        result=[]
        for line in flat:
            if (line.split()[2]== 'chr'+str(i)):
                result.append(line)     
            else:
                continue
        result=sorted(result, key=lambda x:int(x.split()[4]))
        myout.writelines(result)
        myout.close()       

    for line in flat:
        if (line.split()[2]== 'chrX'):
            resultx.append(line)
        elif (line.split()[2]== 'chrY'):
            resulty.append(line)
        else:
            continue
    myoutx=open('chrX_ref.txt','a')
    myouty=open('chrY_ref.txt','a')    
    resultx.sort()
    resulty.sort()
    myoutx.writelines(resultx)
    myouty.writelines(resulty)
    myoutx.close()
    myouty.close()
    

    end=time.time()
    del flat
    del result
    del resultx
    del resulty
    gc.collect()
    print('Time for generating gtf files:'+str(round(end-start))+'s')

#Function that counts total number of variants, chromosomes and variants per chromosome    
def count_chrom(vcfdir):
    vcffile=open('%s'%(vcfdir), 'r')
    vcfall = vcffile.readlines()
    chrom=[]
    var_chr=[]
    var_total=0
    count=0
    chrpool=['chrX','chrY']
    for i in range(1,23):
        chrpool.insert(i-1,'chr'+str(i))
    for k in range(len(vcfall)):
        line = vcfall[k].strip()
        vl = line.split()
        if(vl[0] in chrpool):
            var_total+=1
            if(vl[0] not in chrom):
                chrom.append(vl[0])
                var_chr.append(count)
                count=1
            else:
                count+=1 
    var_chr.append(count)

    del var_chr[0]
    info=[chrom,var_chr,var_total]
    vcffile.close()
    del vcfall
    del chrpool
    gc.collect()
    return info
    
if __name__=='__main__':    
    #vcfdir = input("vcf directory:")
    #outdir = input("output directory:")
    #gtfdir = input("gtf file:")
    #cores = int(input("number of cores:"))
    parser = argparse.ArgumentParser(description='Process VCF files to generate OMI files')
    parser.add_argument('-i', '--input_dir', metavar="DIRNAME", dest='vcf_dir', required=True, help='Required: Input directory for all VCF files.')
    parser.add_argument('-o', '--out_dir', metavar="DIRNAME", dest='outdir', required=True, help='Required: Name of output directory to store aa sequence.')
    parser.add_argument('-f', '--ref_flat', metavar="ref_flat", dest='ref_flat', required=True, help='Required: Reference genomic field used to identify gene region,including path.')
    parser.add_argument('-c', '--num_cores', metavar="NUMCORES", dest='cores', type=int,default=3)
    args = parser.parse_args()

    vcfdir = args.vcf_dir
    outdir = args.outdir
    cores = args.cores
    gtfdir = args.ref_flat
    
    vcfs=os.listdir(vcfdir)
    vcfs=[value for value in vcfs if value.endswith('.vcf')] #remove non-vcf files under vcfdir
    (refpath, gtfname) = os.path.split(gtfdir)
########## create the excel file to record runtimes if file does not exist
    if(not os.path.exists('runtime.xls')):
        workbook=xlwt.Workbook()
        sheet1=workbook.add_sheet('sheet1',cell_overwrite_ok=True)
        colname=['filename','Size','# GSVs','# chroms','# of samples']
        for i in range(len(colname)):
            sheet1.write(0,i,colname[i])
        workbook.save('runtime.xls')
########## deal with each file 
    for file in vcfs:
        #st_mem=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024 #memory usage at the beginning
        starting = time.time()
        time0=splitvcf(vcfdir+'/'+file)
        splitgtf(gtfdir)
        
        info=count_chrom(vcfdir+'/'+file)
        (filename,extension)=os.path.splitext(file)
        filesize=os.path.getsize(vcfdir+'/'+file)/(1024*1024)
        if ('chrM' in info[0]):
            info[0].remove('chrM') 
        vcf_by_chrom=[b+'.vcf' for b in info[0]]
        time1=time.time()
    ######## create multiprocessing pool for copying information
        pool=multiprocessing.Pool(processes=cores,maxtasksperchild=1)
        pool_list1=[]
        pool_list1=[pool.apply_async(writecsv,(vcf_by_chrom[i],)) for i in range(len(info[0]))]
        result1=[c.get() for c in pool_list1]
        csv_name=result1[0]

        pool.close()
        pool.join()
        time2=time.time()-time1
        print("Time for writing all csv:%.2f seconds" %(time2))
    ######## create final OMI files and write the header
        for i in range(0,len(csv_name)):   
            f1=open(csv_name[i],'r')
            header=f1.readline().strip('\n')
            header1=header.split(',')
            #header1.append('gene_name')
            header1.append('change_type1')
            header1.append('AAchange\n')
            header1=",".join(header1)
            f2=open(outdir+filename+'_'+csv_name[i],'a')
            f2.write(header1)
            f1.close()
            f2.close()
        #cores=multiprocessing.cpu_count()-1
        pool=multiprocessing.Pool(processes=cores,maxtasksperchild=1)
    ######## prepare lists containing names of all reference files 
        ref=[b+'_ref.txt' for b in info[0]]
        genome=[refpath+'/'+'genome_'+b+'.fa' for b in info[0]]
        pool_list=[]
    ######## create multiprocessing pool for parsing
        time3=time.time()
        pool_list=[pool.apply_async(vcf_func,(vcf_by_chrom[i],ref[i],genome[i],)) for i in range(len(info[0]))]
        result2=[c.get() for c in pool_list]
        time4=time.time()-time3
        print("Time for parsing whole vcf: %.2f seconds" %(time4))
        
        time5=time.time()
        if (len(csv_name) < cores):
            for j in range(0,len(csv_name)):
                pool_list=[pool.apply_async(writeOMI,(result2[i],csv_name[j],outdir+filename+'_'+csv_name[j],)) for i in range(len(info[0]))] 
        else:
            result=result2[0].copy()
            for j in range(1,len(result2)):
                result.update(result2[j])
            pool_list=[pool.apply_async(writeOMI,(result,csv_name[k],outdir+filename+'_'+csv_name[k],)) for k in range(len(csv_name))] 
        pool.close()
        pool.join()
        time6=time.time()-time5
        print("Time for writing all OMI: %.3f reconds" %(time6))
    ######## delete intermediate files
        chrpool=['chrX','chrY']
        for i in range(1,23):
            chrpool.insert(i-1,'chr'+str(i))
        for i in chrpool:
            os.remove(i+'_ref.txt')
        for k in chrpool:
            os.remove(k+'.vcf')
        for j in csv_name:
            os.remove(j)
        #en_mem=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024#maximum memory usage during the process 
        times=time.time()-starting
    ######## delete intermediate variables to release memory
        del result1
        del result2
        del pool_list
        del pool_list1
        gc.collect()
    ######## record runtime information to excel file
        rexcel=open_workbook('runtime.xls')
        filelist=rexcel.sheet_by_name('sheet1').col_values(0)
        header=rexcel.sheet_by_name('sheet1').row_values(0)
        rows=rexcel.sheets()[0].nrows
        excel=copy(rexcel)
        table=excel.get_sheet(0)
        #if filename already exists, append runtime information
        if(filename in filelist):
            file_index=filelist.index(filename)
            if (cores in header):
                core_index=header.index(cores)
                table.write(file_index,core_index,times)
                table.write(file_index,core_index+1,time0)
                table.write(file_index,core_index+2,time2)
                table.write(file_index,core_index+3,time4)
                table.write(file_index,core_index+4,time6)
            else:
                table.write(0,len(header),cores)
                table.write(0,len(header)+1,'generating vcf')
                table.write(0,len(header)+2,'copy csv')  
                table.write(0,len(header)+3,'parsing')  
                table.write(0,len(header)+4,'write OMI')  
                table.write(file_index,len(header),times)  
                table.write(file_index,len(header)+1,time0)  
                table.write(file_index,len(header)+2,time2)  
                table.write(file_index,len(header)+3,time4)
                table.write(file_index,len(header)+4,time6)
        #if file does not exist, first write file information then runtime 
        else:
            new_row=[filename,filesize,info[2],len(info[0]),len(csv_name)]
            for i in range(len(new_row)):
                table.write(rows,i,new_row[i])
            if (cores in header):
                core_index=header.index(cores)
                table.write(rows,core_index,times)
                table.write(rows,core_index+1,time0)
                table.write(rows,core_index+2,time2)
                table.write(rows,core_index+3,time4)
                table.write(rows,core_index+4,time6)
            else:
                table.write(0,len(header),cores)
                table.write(0,len(header)+1,'splitvcf')  
                table.write(0,len(header)+2,'copy csv')  
                table.write(0,len(header)+3,'parsing')
                table.write(0,len(header)+4,'write OMI')  
                table.write(rows,len(header),times)  
                table.write(rows,len(header)+1,time0)  
                table.write(rows,len(header)+2,time2)  
                table.write(rows,len(header)+3,time4)
                table.write(rows,len(header)+4,time6)
        excel.save('runtime.xls')
    ######## print out filename and total time used            
        print('Time elapsed:'+str(round(times,2))+' seconds')
        print('Number of cores used:%d' %(cores))
        print('VCF filename:%s' %(file))
        #print(st_mem)
        #print(en_mem)
