#### list of function for polymorphism summary file for VCF poolsed ####
#### Le Veve Audrey ##############################
#### 08/10/2025 ##################################


# module require
import math

# filtration and annotation of studied positions based on mean depth

def polymorphism(file_annotation_NS,file_annotation_S,pop,low_cov,high_cov,nb_comp,directory_vcf,qual):
    
    annotation=open(str(file_annotation_NS), "r")
    annotation=annotation.read()
    annotation = annotation.split("\n")
    dict_annotation={}
    for line in annotation:
        line_splitted = line.split(";")
        dict_annotation[line_splitted[0]+";"+line_splitted[1]]=line_splitted
    annotation2=open(str(file_annotation_S), "r")
    annotation2=annotation2.read()
    annotation2 = annotation2.split("\n")
    dict_annotation2={}
    for line in annotation2:
        line_splitted = line.split(";")
        if len(line_splitted)>=2:dict_annotation2[line_splitted[0]+";"+line_splitted[1]]=line_splitted       
    depth = open(str(pop)+"_depth.csv" , "r")
    depth=depth.read()
    depth = depth.split("\n")
    x2=str(pop)+"_depth_filtrered.csv" 
    fichier=open(str(x2), "a")
    fichier.write("chrom;pos;depth;annotation")
    k=0
    while k<len(depth) and depth[k]!="" :
        x=depth[k]
        x=x.split("\t")
        if len(x)>=3 and int(x[2])>=low_cov and int(x[2])<=high_cov: 
            search_annotation = str(x[0])+";"+str(x[1])
            if search_annotation in dict_annotation:
                x="\n"+str(x[0])+";"+str(x[1])+";"+str(x[2])+";NS"
                fichier.write(x)
            elif search_annotation in dict_annotation2:
                x="\n"+str(x[0])+";"+str(x[1])+";"+str(x[2])+";S"
                fichier.write(x)
            else :
                x="\n"+str(x[0])+";"+str(x[1])+";"+str(x[2])+";other"
                fichier.write(x)
        k=k+1
    fichier.close()
    print("step1_ok")
    
    # analyse of vcf file
    depth = open(str(pop)+"_depth_filtrered.csv" , "r")
    depth=depth.read()
    depth = depth.split("\n")
    file_sync = open("GATK_HC_"+str(directory_vcf)+str(pop)+".vcf", "r")
    file_sync=file_sync.read()
    file_sync = file_sync.split("\n")
    deb_vcf=0
    i=0
    deb_vcf2=file_sync[i]
    while deb_vcf2[0]=="#":
        i=i+1
        deb_vcf2=file_sync[i]
    deb_vcf=i
    dict_annotation={}
    for line in file_sync[i:]:
        line_splitted = line.split("\t")
        if len(line_splitted)>1 :
          if float(line_splitted[5])>=qual:
            ind = line_splitted[9]
            ind=ind.split(":")
            AD=str(ind[1])
            AD=AD.split(",")
            AD= list(map(int, AD))
            AD.sort()
            while 0 in AD:del AD[AD.index(0)]
            DP=int(ind[2])
            if len(AD)<=2 and DP>=low_cov and DP<=high_cov:
              MAF=min(AD)*1.0/DP
              pi=(round(MAF*50,0)*round(50-(MAF*50),0)*1.0)/nb_comp
              allele_all=["R",line_splitted[4]]
              allelemin=line_splitted[4]
              if min(AD)==AD[0]:allelemin="R"
              if MAF==1:
                MAF=0.0
                pi=0.0
              dict_annotation[line_splitted[0]+";"+line_splitted[1]]=[pi,MAF,allele_all,allelemin]
    x2=str(pop)+"_analysis_polymorphism_all.csv" 
    fichier=open(str(x2), "a")
    fichier.write("Chrom;pos;type;pi;MAF;allele_all;allele_minor")
    # filtration of files
    for seq in depth[1:]:
        seq=seq.split(";")
        chrom=str(seq[0])
        pos=int(seq[1])
        annot=str(seq[3])
        search_annotation = chrom+";"+str(pos)
        following_line="\n"+str(chrom)+";"+str(pos)+";"+str(annot)+(";0;0;R;R")
        if search_annotation in dict_annotation:
            line=dict_annotation[search_annotation]
            following_line="\n"+str(chrom)+";"+str(pos)+";"+str(annot)+";"+str(line[0])+";"+str(line[1])+";"+str(line[2])+";"+str(line[3])
        fichier.write(str(following_line))
    fichier.close()
    print("step_1_ok")
    
    
## analyse of pi by gene

def pi_gene(pop_files,BED,pi_NS,pi_S,pi):

    x2="Analysis_polymorphism_bed_filtration_candidate.csv" 
    fichier=open(str(x2), "w")
    fichier.write("chrom;start;end;size;nb_site;pi;piNS;piS;ratio;Gene_names;Gene_cluster;Pop")
    
    
    for seq in BED[1:]:
        seq=seq.split("\t")
        chrom=str(seq[0])
        start=int(seq[1])
        end=int(seq[2])
        size=end-start
        GENE=str(seq[3])
        Gene_cluster=str(seq[4])
        
        for pop in pop_files:
            pi2=[]
            piNS=[]
            piS=[]
            for i in range(start,end+1):
                ID=str(pop)+"_"+str(chrom)+"_"+str(i)
                if ID in pi.keys():
                    pi2.append(float(pi[ID]))
                    if ID in pi_NS.keys():piNS.append(float(pi_NS[ID]))
                    if ID in pi_S.keys():piS.append(float(pi_S[ID]))
            
            nb_site=len(pi2)
            if len(pi2)!=0:pi2=sum(pi2)/len(pi2)
            if len(pi2)==0:pi2="NA"
            if len(piNS)!=0:piNS=sum(piNS)/len(piNS)
            if len(piNS)==0:piNS="NA"
            if len(piS)!=0:
                piS=sum(piS)/len(piS)
                if piS!=0:ratio=float(piNS*1.0/piS)
                if piS==0:ratio="Inf"
            if len(piS)==0:
                piS="NA"
                ratio="NA"
            following_line="\n"+str(chrom)+";"+str(start)+";"+str(end)+";"+str(size)+";"+str(nb_site)+";"+str(pi2)+";"+str(piNS)+";"+str(piS)+";"+str(ratio)+";"+str(GENE)+";"+str(Gene_cluster)+";"+str(pop)
            fichier.write(str(following_line))
    
    fichier.close()

    
## Analysis of tajima'D
    
'''
we estimate the tajima D by gene following:
- a1= sum(1/i) for i between n-1 haplotypes
- a2= sum(1/i*i) for i between n-1 haplotypes
- S= number of polymorphic sites
- theta = S/a
- pi2= mean  pi
- d=pi2-theta
- b1=(n+1)/3*(n+1)
- b2=(2*(n*n+n+3))/(9*n*(n-1))
- c1=b1-(1/a1)
- c2=b2-((n+2)/(a1*n-1))+(a2/(a1*a1))
- e1=c1/a1
- e2=c2/(a1*a1+a2)
- V=e1​*S+e2*​S*(S−1)
- D= (d)- sqrt(V)

'''


def Tajima(BED,pop_files,nb_ind,pi,nb_comp,ploidy):

    # estimation of a

    a1=0
    a2=0
    for i in range(1,nb_ind*ploidy):
        a1=a1+(1/i)
        a2=a2+(1/i*i)
    n=nb_ind*ploidy
    b1=(n+1)/3*(n+1)
    b2=(2*(n*n+n+3))/(9*n*(n-1))
    c1=b1-(1/a1)
    c2=b2-((n+2)/(a1*n-1))+(a2/(a1*a1))
    e1=c1/a1
    e2=c2/(a1*a1+a2)
    
    x2="TajimaD.csv" 
    fichier=open(str(x2), "a")
    fichier.write("chrom;start;end;gene_size;Gene_name;Gene_cluster;S;Tajima;pop")
    
    
    for seq in BED[1:]:
        seq=seq.split("\t")
        chrom=str(seq[0])
        start=int(seq[1])
        end=int(seq[2])
        GENE=str(seq[3])
        Gene_cluster=str(seq[4])
        
        for pop in pop_files:
            pi2=[]
            for i in range(start,end+1):
                ID=str(pop)+"_"+str(chrom)+"_"+str(i)
                S=0
                if ID in pi.keys():
                    k=float(pi[ID])
                    pi2.append(float(k))

                    
            theta=S/a
            pi2=sum(pi2)/len(pi2)
            D=pi2-theta
            V=(e1*S)+(e2*S*(S-1))
            V=sqrt(V)
            if V!=0: 
                D= (d)- V
            elif V==0:D="NA"
            following_line="\n"+str(chrom)+";"+str(start)+";"+str(end)+";"+str(end-start)+";"+str(GENE)+";"+str(Gene_cluster)+";"+str(S)+";"+str(D)+";"+str(pop)
            fichier.write(str(following_line))
    
    fichier.close()
    print("step_4_ok")
    

## FST analysis


''' We create a csv file that summarize, based on
a csv files created by step1, the FST between pairwise of pop
of each position
'''
def FST(pop_files,allele,BED,nb_ind):
    pop = pop_files
    pop.remove("")
    size_pop=nb_ind*2
    nb_comp_inter=0
    i=size_pop*2
    while i>0:
        nb_comp_inter=nb_comp_inter+i
        i=i-1
    
    x2="FST_polymorphism_North_America.csv"
    fichier=open(str(x2), "a")
    fichier.write("chrom;start;end;Gene_name;Gene_cluster;pop1;pop2;FST")
    
    for seq in BED[1:]:
        seq=seq.split("\t")
        chrom=str(seq[0])
        start=int(seq[1])
        end=int(seq[2])
        GENE=str(seq[3])
        Gene_cluster=str(seq[4])
        i=0
        for Pop in pop_files[:-1]:
            for Pop2 in pop_files[i+1:]:
                FST=[]
                nb_site_fst=0
                for i in range(start,end+1):
                    ID1=str(Pop)+"_"+str(chrom)+"_"+str(i)
                    ID2=str(Pop2)+"_"+str(chrom)+"_"+str(i)
                    if ID1 in allele.keys() and ID2 in allele.keys():
                        pop1=allele[ID1]
                        pop2=allele[ID2]
                        allelepop1=pop1[5]
                        allelepop1=allelepop1.replace("[","")
                        allelepop1=allelepop1.replace("]","")
                        allelepop1=allelepop1.replace("'","")
                        allelepop1=allelepop1.replace(" ","")
                        allelepop1=allelepop1.split(",")
                        if allelepop1==["R"]:pop1[4]=1
                        allelepop2=pop2[5]
                        allelepop2=allelepop2.replace("[","")
                        allelepop2=allelepop2.replace("]","")
                        allelepop2=allelepop2.replace("'","")
                        allelepop2=allelepop2.replace(" ","")
                        allelepop2=allelepop2.split(",")
                        if allelepop2==["R"]:pop2[4]=1
                        if len(allelepop1)<3 and len(allelepop2)<3:
                            if len(allelepop2)>=1 or len(allelepop1)>=1 :
                                nb_site_fst=nb_site_fst+1
                                # on calcul pi intra moyen
                                pi_intra=(float(pop1[3])+float(pop2[3]))/2
                                # on liste les variants presents dans l'ensemble des pop
                                allele_all=allelepop1+allelepop2
                                allele_all=set(allele_all)
                                allele_all=list(allele_all)
                                if len(allele_all)==1:FST.append(0)
                                elif len(allele_all)>1:
                                    # on compte le nombre de copies de chaque allele dans l'ensemble
                                    FREQ=[]
                                    for i2 in allele_all:
                                        x=0
                                        if i2 in allelepop1 and pop1[6]==i2:x=x+float(pop1[4])*size_pop
                                        if i2 in allelepop1 and pop1[6]!=i2:x=x+(1-float(pop1[4]))*size_pop
                                        if i2 in allelepop2 and pop2[6]==i2:x=x+float(pop2[4])*size_pop
                                        if i2 in allelepop2 and pop12[6]!=i2:x=x+(1-float(pop2[4]))*size_pop
                                        FREQ.append(x)
                                    i2=0
                                    pi_inter=0
                                    while i2<len(FREQ):
                                        pi_inter=pi_inter+(sum(FREQ[i2+1:])*FREQ[i2])
                                        i2=i2+1
                                    pi_inter=pi_inter/nb_comp_inter
                                    #on calcul FST
                                    fst=(pi_inter-pi_intra)/pi_inter
                                    if fst<0:fst=0
                                    FST.append(fst)
                if len(FST)>0:FST=sum(FST)/len(FST)
                elif len(FST)==0:FST="NA"
                following_line="\n"+str(chrom)+";"+str(start)+";"+str(end)+";"+str(GENE)+";"+str(Gene_cluster)+";"+str(Pop)+";"+str(Pop2)+";"+str(FST)
                fichier.write(str(following_line))
            i=i+1


    

    fichier.close()


