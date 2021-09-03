# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import os
import pandas as pd
#import mechanize  import urllib
import mechanicalsoup
#from urllib.request import Request
from urllib.request import urlopen
from bs4 import BeautifulSoup
import time
import csv
import argparse

#variants=[value for value in allfiles if value.endswith('.txt')] 


####---------------------------------------python3.6--------------------------------------------------
#browser = mechanicalsoup.Browser()
#page=browser.get("http://fathmm.biocompute.org.uk/fathmm-xf/")
#submission = page.soup.select("#myForm")[0]
#print(submission)
#submission["batch"] = "2,15424278,C,G"
#submission["hg38"]="hg38"
#response = browser.submit(submission,page.url)
#content = response.content.decode("UTF-8", "ignore")

####-------------------------------------another version---------------------------------------------
def fathmm(omifile,fatout,build='hg38',sleeptime=60):
    omi = pd.read_csv(omifile)
    directory,name=os.path.split(omifile)
    name=name.split(".")[0]
    del directory
    temp=omi[['left', 'ref_seq']]
    alt=[]
    for i in range(omi.shape[0]):
        if omi.iloc[i,5]==omi.iloc[i,6]:
            alt.append(omi.iloc[i,7])
        else:
            alt.append(omi.iloc[i,6])
    chrom=omi['chrom'].apply(lambda x : x[3:])
    temp.insert(0,'chrom',chrom)
    temp.insert(3,'alt',alt)
    temp['left']=temp['left'].apply(lambda x:str(x))
    temp=temp.drop_duplicates()
    fathmm=temp.apply(lambda x:x['chrom']+','+ x['left']+','+ x['ref_seq']+','+ x['alt'],axis=1) 
    del temp,i,alt,chrom,omi
    
    
    browser = mechanicalsoup.StatefulBrowser()
    browser.open("http://fathmm.biocompute.org.uk/fathmm-xf/")
    browser.get_url()
    browser.get_current_page()
    browser.select_form('form[name="myForm"]')
    #browser.get_current_form().print_summary()
    #browser["batch"]=lines
    browser["batch"]='\n'.join(fathmm.tolist())
    if build=='hg38':
        browser["hg38"]="hg38"
    
    response = browser.submit_selected()
    time.sleep(sleeptime)
    #print(response.text)
    browser.get_url()
    temp=response.text.split("session=")
    session=temp[1].split("\'")[0]
    resulturl="http://fathmm.biocompute.org.uk/fathmm-xf/cgi-bin/results.cgi?session="+session
    #resultpage=urllib.request.urlopen(resulturl).read()
    resultpage=urlopen(resulturl).read()
    #import requests
    #requests.get(url)
    soup = BeautifulSoup(resultpage, 'html.parser')
    table = soup.find('table', attrs={'class': 'table table-striped table-bordered'})
    headings = [th.get_text() for th in table.find("tr").find_all("th")]
    datasets = []
    for row in table.find_all("tr")[1:]:
        dataset = [td.get_text() for td in row.find_all("td")]
        datasets.append(dataset)
    
    output=open(fatout+name+"_fathmm.csv",'a',newline='')
    writer=csv.writer(output)
    writer.writerow(headings)
    for ele in datasets:
        writer.writerow(ele)
    output.close()
    return
   
if __name__=='__main__':   
    parser = argparse.ArgumentParser(description='prepare for FATHMM submission and parse results')
    parser.add_argument('-i', '--input_file', metavar="FILENAME", dest='omi', required=True, help='Required: Input OMI file.')
    parser.add_argument('-o', '--output_dir', metavar="DIRNAME", dest='outdir', required=True, help='Required: Name of output directory to store results.')
    parser.add_argument('-r', '--reference_build', metavar="ref_build", dest='ref', default='hg38', help='Choose reference genome going to be used')
    parser.add_argument('-t', '--runtime', metavar="runtime", dest='time', type=int, default=60)
    args = parser.parse_args()

    omifile = args.omi
    outdir = args.outdir
    build=args.ref
    sleeptime=args.time
##    omifile = input("OMI file:")
##    outdir = input("Output directory:")
    fathmm(omifile,outdir,build,sleeptime)  
    
    
    
    









   
    
