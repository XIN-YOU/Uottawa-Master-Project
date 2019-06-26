# When run this on the cluster (as an array job), the jobindex will 
# be set based on the array index in job submission script 
jobindex=1000

# load an R library 
library(trio)
library(sim1000G)

#set directory
#setwd("/Users/Shelly/Documents/Project_R_code/1000")
vcfdir="/Users/Shelly/Documents/Fall_Project/subset chr vcf/" 
#vcfdir="/global/home/hpc4429/XINYOU/DATA/subset chr vcf/" 

# Simulation settings (p0 is get from another r script)
numtrios=10
p0=0.35
numsim=1000


# Chromosomes with variant data  
chroms=c(1,2,3,5,6,8,9,11,12,14,15,17,18,21)

# Store results in a matrix and write out at the end
resmat=matrix(nrow=numsim,ncol=2)
colnames(resmat)=c("SNP","GXE_min_p")
newmat=matrix(nrow = numsim,ncol = 2)
colnames(newmat)=c("SNP","main_gene_effect(min p)")
#Loop for number of simulations
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
  g.main <- colTDT(genomat)
  minrow=which.min(g.main$pval)
  minval=min(g.main$pval)
  gxe.add <- colGxE(genomat,envvar, model = "additive")
  minindex<-which.min(gxe.add$pval[,2])
  resmat[k,]<-c(minindex,gxe.add$pval[minindex,2])
  newmat[k,]<-c(minrow,minval)

  
}

write.table(resmat, paste("GXE","-",jobindex,".txt",sep=""),quote=F,row=F,col.names = T,append = F)
write.table(newmat, paste("MainGene","-",jobindex,".txt",sep=""),quote=F,row=F,col.names = T,append = F)
# for each simulation, keep only the relevant information. In our case, it's the minimum p-value from a simulation

