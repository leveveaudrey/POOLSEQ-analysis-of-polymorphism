#### Polymorphism summary file for VCF poolsed ####
#### Le Veve Audrey ##############################
#### 08/10/2025 ##################################

''' 
We create a csv file that summarize, based on
vcf and annotation files, the polymorphism
pi and the tajima'D of each population 

'''

# module required
import os
import sys
import argparse
import math
import random
from function_for_poolseq import *

def poolseq_Genetic_population_program(ArgsVal):
    global args
    description = """ Study of genetic variation accross population of poolseq data """
    parser = argparse.ArgumentParser(prog="Set of analysis",description=description)
    parser.add_argument("-bed", "--bed_file", default="candidate.bed")
    parser.add_argument("-pop", "--pop_files", default="pop.txt")
    parser.add_argument("-annotNS", "--annotation_fileNS", default="lyrata_filter_NS.csv")
    parser.add_argument("-annotS", "--annotation_fileS", default="lyrata_filter_S.csv")
    parser.add_argument("-directory_vcf", "--directory_vcf", default="")
    parser.add_argument("-low_cov", "--low_cov", default="25")
    parser.add_argument("-high_cov", "--high_cov", default="500")
    parser.add_argument("-qual", "--qual", default="10")
    parser.add_argument("-nb_ind", "--nb_ind", default="25")
    parser.add_argument("-p", "--ploidy", default="2")
    
    args = parser.parse_args(ArgsVal)
    bed_file=str(args.bed_file)
    pop_files=str(args.pop_files)
    directory_vcf=str(args.directory_vcf)
    file_annotation_NS=str(args.annotation_fileNS)
    file_annotation_S=str(args.annotation_fileS)
    low_cov=int(args.low_cov)
    high_cov=int(args.high_cov)
    qual=int(args.qual)
    nb_ind=int(args.nb_ind)
    ploidy=int(args.ploidy)
    
    x=nb_ind-1
    nb_comp=x
    while (x > 0): 
        x=x-1
        nb_comp=nb_comp+x
        

    pop_files=open(str(pop_files),"r")
    pop_files=pop_files.read()
    pop_files=pop_files.split('\n')
    
    BED=open(str(bed_file),"r")
    BED=BED.read()
    BED=BED.split('\n')
    
    
    ### fichier avec le pi a chaque pos + annotation 
    for pop in pop_files:
        if os.path.exists(str(pop)+"_analysis_polymorphism_all.csv")==False:
            polymorphism(file_annotation_NS,file_annotation_S,pop,low_cov,high_cov,nb_comp,directory_vcf,qual)
            
    ### dictionaries of pi and allele
    pi_NS={}
    pi_S={}
    pi={}
    allele={}
    for pop in pop_files:
        input_file=open(str(pop)+"_analysis_polymorphism_all.csv", "r")
        input_file=input_file.read()
        input_file = input_file.split("\n")
        
        for i in input_file[1:]:
            if i!="":
                line=i.split(";")
                ID=str(pop)+"_"+str(line[0])+"_"+str(line[1])
                pi[ID]=float(line[3])
                allele[ID]=line
                if str(line[2])=="NS":pi_NS[ID]=float(line[3])
                if str(line[2])=="S":pi_S[ID]=float(line[3])
    
    ## analyse of pi by gene
    pi_gene(pop_files,BED,pi_NS,pi_S,pi)
    
    ## Analysis of tajima'D
    Tajima(BED,pop_files,nb_ind,pi,nb_comp,ploidy)
    
    ## FST analysis
    if len(pop_files)>1:FST(pop_files,allele,BED,nb_ind)
        


if __name__=='__main__':
    
    ArgsVal = sys.argv[1:]
    poolseq_Genetic_population_program(ArgsVal)


