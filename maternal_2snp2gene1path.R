###########################################################################################
library(rio)
library(trio)
library(sim1000G)
############################### pathway information########################################

# Cytosolic Metabolism Pathway including genes as follow####
path1 = c("MTHFR","DHFR","CTH","MTR","MTRR","BHMT2",
          "BHMT","GNMT","CBS","DMGDH","MTHFD1",
          "MTHFS","SHMT1","FTCD","TYMS","ALDH1L1") 

#row range for all the genes on pathway
#  1  12  13  37  38  56  63 145 146 233 234 253 254 258 259 267 268 279 280 283 450 462
# 463 486 487 492 499 512 513 529 545 559

# Mitochondrial Metabolism Pathway including genes as follow####
path2 = c("MTHFD2","FPGS","SHMT2","MTHFD1L","ALDH1L2")

#row range for all the genes on pathway 2
#57:62, 284:387, 395:401, 425:427, 428:449


# Cross-membrane Transport Pathway including genes as follow###
path3 = c("FOLH1","FOLR3","FOLR1","FOLR2","SLC19A1","SLC46A1","SLC25A32")
#row range for all the genes on pathway 3
#388:394, 402:412, 413:418, 419:420, 421:424, 493:498, 530:544

##########
m1= matrix(ncol=2,nrow=sum(length(path1),length(path2),length(path3)))
m1[,1:2]=c(c(rep(1,length(path1)),rep(2,length(path2)),rep(3,length(path3))),c(path1,path2,path3))

trans=matrix(c(1,12 ,13 ,37,38,56,63,145,146,233, 234, 253, 254, 258, 259, 267, 268, 279, 280, 283,
               450, 462,463, 486, 487, 492, 499, 512, 513, 529, 545, 559,57,62, 284,387, 395,401, 425,427, 428,449,
               388,394, 402,412, 413,418, 419,420, 421,424, 493,498, 530,544),byrow = T,ncol = 2 )
pinfo= cbind(m1,trans)
colnames(pinfo)=c("pathway","genename","row","range")
pinfo=as.data.frame(pinfo)
#write.csv(pinfo,"/Users/Shelly/Documents/Project_R_code/pathway/pathway-information.csv")

############################################################################################
############################################################################################






###########################       -SIMULATION-          ###################################

jobindex = 221
################find minimum p values for each gene########################
#set directory
#setwd("/Users/Shelly/Documents/Project_R_code/gene_min-p")
#vcfdir="/global/home/hpc4429/XINYOU/DATA/subset chr vcf/" 
vcfdir="/Users/Shelly/Documents/Fall_Project/subset chr vcf/" 

#Simulation settings (p0 is get from another r script)
numtrios=1000
p0=0.35
numsim=3

chroms=c(1,2,3,5,6,8,9,11,12,14,15,17,18,21)
# Store results in a matrix and write out at the end
# newmat=matrix(nrow=numsim,ncol=1)
# colnames(newmat)=c("Cytosolic Metabolism pathway")
geno_list=list()
for(k in 1:3){
  ######################  generate genotypes  #######################################
  genomat=NULL
  for (j in chroms){
    SIM$reset()
    vcf_file=paste(vcfdir,"MengcheSNPs_1000G_Euro_chr",j,".vcf",sep="")
    vcf = readVCF( vcf_file , maxNumberOfVariants = 600 , min_maf = NA , max_maf = NA ) 
    readGeneticMap( chromosome = j )
    startSimulation( vcf, totalNumberOfIndividuals = 3*numtrios )
    
    allfamilies=NULL
    for (p in 1:numtrios){
      fam = newFamilyWithOffspring(paste("FID-",p,sep=""),1)
      allfamilies=rbind(allfamilies,fam)
    }
    genotype=retrieveGenotypes(allfamilies$gtindex)
    genomat=cbind(genomat,genotype)
  }
  geno_list[[k]]=genomat
}

###################  generate environment variable  ##############################
envvar=rbinom(numtrios,1,1-p0)

######################  generate pedfile  #######################################
pid_pattern=c(1,2,3)
pid_newcol=rep(pid_pattern,numtrios)
my_pattern= paste(allfamilies[,1],"_", pid_newcol,sep="")
colnames(genomat)=paste("SNP",1:559,sep="")
rownames(genomat)=my_pattern

############################################################################################ 
############################################################################################

snp=vector()
snp=samplegenomes(nsnp=2,ngene=2,npathway=1) 
father_index = seq(1,3000,3) 
mather_index = seq(1,3000,3) + 1
child_index=seq(1,3000,3)+2
#maternal_genotype=genomat[mather_index, snp]
b1=2.5;b2=2.5;b3=-5;b4=-4 #this is the exp(logit),all the factors,b4 is the interaction between
#maternal genotype and enviromental, I gave it 1.
y=vector()
odds=vector()
prob_logit=vector()
disease=vector()

maternal_genotype=list()
disease_genomat=list()
disease_genomat_final=list()
maternalgeno1=matrix(nrow = 1000,ncol = 1)
maternalgeno2=matrix(nrow = 1000,ncol = 1)
maternalgeno3=matrix(nrow = 1000,ncol = 1)
add.maternal.geno=list(maternalgeno1,maternalgeno2,maternalgeno3)
for(k in 1:3){
  maternal_genotype[[k]]=as.matrix(geno_list[[k]][mather_index, c(snp[1],snp[2],snp[3],snp[4])])
  for(q in 1:1000){
    add.maternal.geno[[k]][q]=sum(maternal_genotype[[k]][q,])
  }
  
  for(i in 1:length(add.maternal.geno[[k]])){
    y[i]=add.maternal.geno[[k]][i]*envvar[i]
    odds[i]=exp(b1+b2*add.maternal.geno[[k]][i]+b3*envvar[i]+b4*y[i])
    prob_logit[i]=odds[i]/(1+odds[i])
    disease[i]=rbinom(1,1,prob_logit[i])
  }
  x=which(disease==0)
  f=father_index[x]
  m=mather_index[x]
  c=child_index[x]
  d.index=c(f,m,c)
  #length(d.index)
  disease_genomat[[k]]=geno_list[[k]][-d.index,]
  
}
for(l in 1:3){
  disease_genomat_final[[l]]=disease_genomat[[l]][-c(1001:nrow(disease_genomat[[l]])),]
  
}
mm1=as.data.frame(disease_genomat_final[[1]])
mm2=as.data.frame(disease_genomat_final[[2]])
mm3=as.data.frame(disease_genomat_final[[3]])
disease_genomat_finalfinal=rbind(mm1,mm2,mm3)
pid_pattern=c(1,2,3)
pid_newcol=rep(pid_pattern,numtrios)
my_pattern= paste(allfamilies[,1],"_", pid_newcol,sep="")
colnames(disease_genomat_finalfinal)=paste("SNP",1:559,sep="")
rownames(disease_genomat_finalfinal)=my_pattern

write.table(disease_genomat_finalfinal,"/Users/Shelly/Documents/Project_R_code/MATERNAL/2snp2gene1path/disease_geno221.ped",
            quote=F,col.names = T,row.names = T)



disease.data=read.table("/Users/Shelly/Documents/Project_R_code/MATERNAL/2snp2gene1path/disease_geno221.ped",header = T)
en.index=which(envvar==0)
length(en.index)
paternal.gene=disease.data[father_index,]
maternal.gene=disease.data[mather_index,]
child.gene=disease.data[child_index,]
row1=father_index[en.index]
row2=mather_index[en.index]
therow=sort(c(row1,row2))
thedata=disease.data[therow,]
write.csv(thedata,"/Users/Shelly/Documents/Project_R_code/MATERNAL/2snp2gene1path/diseaseparents_e=0_221.csv")
dad.index=vector()
mom.index=vector()
dad.index=seq(1,length(en.index),2)
mom.index=dad.index+1
count0=0
count1=0
count2=0
count3=0
count4=0
count5=0
for(i in 1:length(en.index)){
  for(j in 1:559){
    if(thedata[dad.index[i],j]==1&&thedata[mom.index[i],j]==0){
      count0=count0+1
    } 
    
    if(thedata[dad.index[i],j]==0&&thedata[mom.index[i],j]==1){
      count1=count1+1
      
    }
    if(thedata[dad.index[i],j]==0&&thedata[mom.index[i],j]==2){
      count2=count2+1
    }
    if(thedata[dad.index[i],j]==2&&thedata[mom.index[i],j]==0){
      count3=count3+1
    }
    if(thedata[dad.index[i],j]==1&&thedata[mom.index[i],j]==2){
      count4=count4+1
    }
    if(thedata[dad.index[i],j]==2&&thedata[mom.index[i],j]==1){
      count5=count5+1
    }else{next}
  }
  return(count0)
  return(count1)
  return(count2)
  return(count3)
  return(count4)
  return(count5)
}
# count0: 69, count1: 85, count2: 11, count3: 14, count4: 68, count5: 30

p01=count0/length(en.index)
p10=count1/length(en.index)
p02=count2/length(en.index)
p20=count3/length(en.index)
p12=count4/length(en.index)
p21=count5/length(en.index)
p01;p10;p02;p20;p12;p21
#[1] 0.2664756
#[1] 0.2722063
#[1] 0.06303725
#[1] 0.09169054
#[1] 0.09455587
#[1] 0.06876791 this is the probabilities when environment=0;situation is 2 snps,2 genes, 1pathway.

en.index=which(envvar==1)
length(en.index)
paternal.gene=disease.data[father_index,]
maternal.gene=disease.data[mather_index,]
child.gene=disease.data[child_index,]
row1=father_index[en.index]
row2=mather_index[en.index]
therow=sort(c(row1,row2))
thedata=disease.data[therow,]
write.csv(thedata,"/Users/Shelly/Documents/Project_R_code/MATERNAL/2snp2gene1path/diseaseparents_e=1_221.csv")
dad.index=vector()
mom.index=vector()
dad.index=seq(1,length(en.index),2)
mom.index=dad.index+1
count0=0
count1=0
count2=0
count3=0
count4=0
count5=0
for(i in 1:length(en.index)){
  for(j in 1:559){
    if(thedata[dad.index[i],j]==1&&thedata[mom.index[i],j]==0){
      count0=count0+1
    } 
    
    if(thedata[dad.index[i],j]==0&&thedata[mom.index[i],j]==1){
      count1=count1+1
      
    }
    if(thedata[dad.index[i],j]==0&&thedata[mom.index[i],j]==2){
      count2=count2+1
    }
    if(thedata[dad.index[i],j]==2&&thedata[mom.index[i],j]==0){
      count3=count3+1
    }
    if(thedata[dad.index[i],j]==1&&thedata[mom.index[i],j]==2){
      count4=count4+1
    }
    if(thedata[dad.index[i],j]==2&&thedata[mom.index[i],j]==1){
      count5=count5+1
    }else{next}
  }
  return(count0)
  return(count1)
  return(count2)
  return(count3)
  return(count4)
  return(count5)
}


p01=count0/length(en.index)
p10=count1/length(en.index)
p02=count2/length(en.index)
p20=count3/length(en.index)
p12=count4/length(en.index)
p21=count5/length(en.index)
p01;p10;p02;p20;p12;p21
#[1] 0.1182796
#[1] 0.1198157
#[1] 0.04608295
#[1] 0.04454685
#[1] 0.04454685
#[1] 0.07526882 this is the probability when environment=1, 2snp 2genes 1pathway.


##########################             -SAMPLE DATA-          #####################
pathmat=pinfo
#       -------random selecting function------           #
# there are three pathways in our study, t is the pathway number, means pathway1, pathway2,
#or pathway3, you can check the define of three pathways from the top of this scripts
# x is the number of genes, so it's either one or two.
gfun=function(t,x){
  if (t==1){
    y=pathmat[sample(1:16,x),]
    return(y)
  } else if (t==2){
    y=pathmat[sample(17:21,x),]
    return(y)
  } else {
    y=pathmat[sample(22:28,x),]
    return(y)
  }
  
}

# -------------------------------------------------------------------------------
#random sample function for snp level gene level pathway level. nsnp represents the number of snps we want, ngene
#represents the number of genes we want, and npathway represents the number of the pathways we want.
lsnp=NULL
samplegenomes=function(nsnp,ngene,npathway){
  if (nsnp == 1&&ngene==1&&npathway==1){
    t=sample(c(1:3),size=1)
    geno=gfun(t,1) #randomly choose one gene
    lsnp=sample(x=as.numeric(as.character(geno[1,3])):as.numeric(as.character(geno[1,4])),size=1) # sample snps
    return(lsnp)
  }
  if (nsnp ==2&&ngene==1&npathway==1){
    t=sample(c(1:3),size=1)
    geno=gfun(t,1) #randomly choose one gene
    lsnp=sample(x=as.numeric(as.character(geno[1,3])):as.numeric(as.character(geno[1,4])),size=2) # sample snps
    return(lsnp)
  }
  if (nsnp==1&&ngene==2&&npathway==1){
    t=sample(c(1:3),size=1)
    geno=gfun(t,2) #randomly choose one gene
    lsnp1=sample(as.numeric(as.character(geno[1,3])):as.numeric(as.character(geno[1,4])),size=1) # sample snps
    lsnp2=sample(as.numeric(as.character(geno[2,3])):as.numeric(as.character(geno[2,4])),size=1) # sample snps
    step=rbind(lsnp1,lsnp2)
    return(step)
  }
  if (nsnp==2&&ngene==2&&npathway==1){
    geno = list()
    step = NULL
    lsnp1=data.frame()
    t=sample(c(1:3),1)
    geno=gfun(t,2)
    lsnp1=sample(as.numeric(as.character(geno[1,3])):as.numeric(as.character(geno[1,4])),size=2)
    lsnp2=sample(as.numeric(as.character(geno[2,3])):as.numeric(as.character(geno[2,4])),size=2)
    step=rbind(lsnp1,lsnp2)
    return(step)
  }
  if (nsnp==1&&ngene==1&&npathway==2){
    geno=list()
    step=NULL
    t=sample(c(1:3),size = 2,replace = FALSE)
    geno1=gfun(t[[1]][1],1)
    geno2=gfun(t[[2]][1],1)
    lsnp1=sample(x=as.numeric(as.character(geno1[1,3])):as.numeric(as.character(geno1[1,4])),size=1) # sample snps
    lsnp2=sample(x=as.numeric(as.character(geno2[1,3])):as.numeric(as.character(geno2[1,4])),size=1) # sample snps
    step=rbind(lsnp1,lsnp2)
    return(step)
  }
  if (nsnp==1&&ngene==2&&npathway==2){
    geno=list()
    step=NULL
    t=sample(c(1:3),size = 2,replace = FALSE)
    geno1=gfun(t[[1]][1],2)
    geno2=gfun(t[[2]][1],2)
    lsnp1=sample(x=as.numeric(as.character(geno1[1,3])):as.numeric(as.character(geno1[1,4])),size=1) # sample snps
    lsnp2=sample(x=as.numeric(as.character(geno1[2,3])):as.numeric(as.character(geno1[2,4])),size=1) # sample snps
    lsnp3=sample(x=as.numeric(as.character(geno2[1,3])):as.numeric(as.character(geno2[1,4])),size=1) # sample snps
    lsnp4=sample(x=as.numeric(as.character(geno2[2,3])):as.numeric(as.character(geno2[2,4])),size=1)
    step=rbind(lsnp1,lsnp2,lsnp3,lsnp4)
    return(step)
  }
  if (nsnp==2&&ngene==1&&npathway==2){
    geno=list()
    step=NULL
    t=sample(c(1:3),size = 2,replace = FALSE)
    geno1=gfun(t[[1]][1],1)
    geno2=gfun(t[[2]][1],1)
    lsnp1=sample(x=as.numeric(as.character(geno1[1,3])):as.numeric(as.character(geno1[1,4])),size=2) # sample snps
    lsnp2=sample(x=as.numeric(as.character(geno2[1,3])):as.numeric(as.character(geno2[1,4])),size=2) # sample snps
    step=rbind(lsnp1,lsnp2)
    return(step)
  }
  if (nsnp==2&&ngene==2&&npathway==2){
    geno=list()
    step=NULL
    t=sample(c(1:3),size = 2,replace = FALSE)
    geno1=gfun(t[[1]][1],2)
    geno2=gfun(t[[2]][1],2)
    lsnp1=sample(x=as.numeric(as.character(geno1[1,3])):as.numeric(as.character(geno1[1,4])),size=2) # sample snps
    lsnp2=sample(x=as.numeric(as.character(geno1[2,3])):as.numeric(as.character(geno1[2,4])),size=2)
    lsnp3=sample(x=as.numeric(as.character(geno2[1,3])):as.numeric(as.character(geno2[1,4])),size=2) # sample snps
    lsnp4=sample(x=as.numeric(as.character(geno2[2,3])):as.numeric(as.character(geno2[2,4])),size=2)
    step=rbind(lsnp1,lsnp2,lsnp3,lsnp4)
    return(step)
  }
}  





#--------------------------------------------------------------------------------------------------------




