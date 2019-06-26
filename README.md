An Algorithm for Simulating Data for the Investigation of Maternal-Gene-Environment and Pathway Effects

author: Xin You
Project submitted in partial fulﬁllment of the requirements for the degree of Master of Science Mathematics and Statistics
Department:  Mathematics and Statistics Faculty of Science, University of Ottawa
Date:2019

This poject was supervised by Dr.Kelly Burkett, Faculty of Science, University of Ottawa and Marie-Helene RoyGagnon, Faculty of Medicine, University of Ottawa. 


World-wide, there is approximately one case of oral cleft disease for every seven hundred births, making it a fairly common birth defect. Researchers are, therefore interested in discovering how genetic and environmental risk factors contribute to this disease risk. Birth defects may be linked to the mother’s genetic variants in conjunction with environmental factors. These factors may contribute to the risk of disease while the baby is in utero. Luckily, pathway-based association analyses have become popular approaches for discovering how genetic variants contribute to diseases, such as cleft palate.


In this project, we develop an approach to simulate data on trios where the probability of disease is increased by maternal genetic effects and maternal gene-environment effects on genes in selected pathways. Development of the simulation algorithm was motivated by a need to compare the performance of pathway-based statistical approaches in detecting maternal gene-environment effects on a cleft palate data set collected from trios. We illustrate the use of the simulation algorithm by investigating the type 1 error rate when there are no true genetic effects and by examining mating type asymmetry when there are true maternal gene-environment effects. Our algorithm will be of use for investigating the power and type 1 error

Keywords: Orofacial cleft, maternal effects, simulation, genetic effects, gene-environment interactions, causal variants.

We built the model under the Null hypothesis and Alternative hypothesis. Null hypothesis: there are no genetic environment interactions between genotypes and phenotypes. Alternative hypothesis : there are genetic environment interactions between genotypes and phenotypes.

To test the Null hypothesis, we did the analysis in SNP level, gene level and pathway level. For testing the SNPs associations, the codes are in job1new.R, with a shell script job1.sh. For gene level experiment, we chose 4 genes as targets. The names of the 4 genes are CTH, MTHFR, MTHFD1L, and SHMT2. The codes for analyzing these genes are in CTH.r, MTHFR.r, MTHFD1L.r AND SHMT2.r, with 4 shell scripts as well.
For the pathway level analysis, we chose Cytosolic Metabolism Pathway.

For the alternative hypothesis, which we think there is associations. To test the maternal effects, we built a function for randomly choosing variants within genes within pathways. Under each senario, we extract the proportion of the mother's genotype and father's genotype, and check the mating type asymmetry.


1. Randomly choose 1 pathway from the three, randomly choose 1 gene from n genes on that pathway, and randomly choose 1 SNP as the causal variant with a gene within the pathway. We only have 1 causal variant in this scenario.

2. Randomly choose 1 pathway from the three, randomly choose 2 genes from n genes on that pathway, randomly choose 1 SNP as the causal variant from each genes we selected within the pathway. 

3. Randomly choose 1 pathway from the three, randomly choose 1 gene within that pathway, and randomly choose 2 SNPs as causal variants within that gene within that pathway. In this scenario, we also have 2 causal variants in total.

4. Randomly choose 1 pathway from the three, randomly choose 2 genes within that pathway, and randomly choose 2 SNPs as causal variants within that gene and pathway. There will be 4 causal variants in total for this case.

5. Randomly choose 2 pathways from the three, randomly choose 1 gene from each pathways we selected, and randomly choose 1 SNP from each gene within the pathway. There will be 2 causal variants in total from two diﬀerent pathways and diﬀerent genes within that pathway.

6. Randomly choose 2 pathways from the three variables, randomly choose 1 gene within each pathways, and randomly choose 2 SNPs as causal variants within each gene. There will be 4 causal variants in total.

7. Randomly choose 2 pathways from the three, on each pathway, randomly choose 2 genes, and randomly choose 1 SNP from each gene within the pathway. There will be 4 causal variants in tota.

8. Randomly choose 2 pathways, randomly choose 2 genes within each pathway, and randomly choose 2 SNPs within each gene within that pathway. In this scenario, we will have 8 causal variants in total. 
