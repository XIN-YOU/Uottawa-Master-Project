jobindex=1000
######################## settings ######################### 
library(trio)
library(sim1000G)
################find minimum p values for each gene########################
#set directory
#setwd("/Users/Shelly/Documents/Project_R_code/gene_min-p")
vcfdir="/global/home/hpc4429/XINYOU/DATA/subset chr vcf/" 
#Simulation settings (p0 is get from another r script)
numtrios=1000
p0=0.35
numsim=1000

chroms=c(1,2,3,5,6,8,9,11,12,14,15,17,18,21)
# Store results in a matrix and write out at the end
newmat=matrix(nrow=numsim,ncol=1)
colnames(newmat)=c("MTHFD1L-min-p")

for(k in 1:numsim){
  
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
  
  ###################  generate environment variable  ##############################
  envvar=rbinom(numtrios,1,1-p0)
  
  ######################  generate pedfile  #######################################
  pid_pattern=c(1,2,3)
  pid_newcol=rep(pid_pattern,numtrios)
  my_pattern= paste(allfamilies[,1],"_", pid_newcol,sep="")
  colnames(genomat)=paste("SNP",1:559,sep="")
  rownames(genomat)=my_pattern
  ######################  TDT test  #######################################
  gxe.add <- colGxE(genomat,envvar, model = "additive")
  minindex<-min(gxe.add$pval[284:387,2])
  newmat[k,]<-c(minindex)
  
  
}

write.table(newmat, paste("MTHFD1L",jobindex,"_","min-p",".txt",sep=""),quote=F,row=F,col.names = T,append = F)

#THIS is G(SHMT2)XE p-value

