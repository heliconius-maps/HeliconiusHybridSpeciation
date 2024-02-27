#!/usr/bin/env Rscript 

#./qtl_pheromones.r &> out.log &
	
library(qtl)

r_qtl_F2s<-read.cross("csvr",file="qtl_dat_pheromones_ordinations_PCA_F2s_only_27-09-22_trim_tips.csv",genotypes=c("EE","PE","PP"),alleles=c("E","P"))
r_qtl_F2s<-jittermap(r_qtl_F2s)
r_qtl_F2s<-calc.genoprob(r_qtl_F2s,step=1, off.end=0.0, error.prob=0.001,map.function="haldane",stepwidth="fixed")


#markers for confidence intervals

list_all_markers<-read.csv("C:/Users/NER334/Work/Harvard/Genomics/RAD_2018/elev_x_pard/QTL_F2s/qtl/Fst/list_all_markers_Fst_trim_tips_reset_cM.csv",header=T,sep=",")
list_all_markers$X<-NULL
list_all_markers$Fst<-NULL
list_all_markers$LG[which(list_all_markers$LG=="21")]<-"X"

# list_all_markers<-read.table("C:/Users/NER334/Work/Harvard/Genomics/RAD_2018/elev_x_pard/linkage_maps/all_markers")
# colnames(list_all_markers)<-list_all_markers[1,]
# list_all_markers<-list_all_markers[-1,]
# list_all_markers$cM<-as.numeric(list_all_markers$cM)
# list_all_markers$LG[which(list_all_markers$LG=="21")]<-"X"

n.perm<-1000

#Use these when you plot a separate sex chrom sig threshold
#These x values are specific to this particular map and not generalisable
# autosomes_x_coords_4_plot<-c(0,2165)
# sex_x_coords_4_plot<-c(2190,2244)
#for trimmed datafile
autosomes_x_coords_4_plot<-c(0,1518)
sex_x_coords_4_plot<-c(1519,2000)


for (i in 10:160){
	pheno_name<-names(r_qtl_F2s$pheno)[i]
	png(filename = paste("results/pheromones_ordination/r_qtl_approach_PCA_27-09-22_trim_tips/",pheno_name,".png",sep=""), width = 1000, height = 300)
	#par(mfrow=c(2,1))
		
	r_qtl_F2s.scan<-scanone(r_qtl_F2s,pheno.col=i,model="normal",method="hk")
	r_qtl_F2s.perm<-scanone(r_qtl_F2s,pheno.col=i,n.perm=n.perm,model="normal",method="hk",perm.Xsp=FALSE)	
	plot(r_qtl_F2s.scan,main="F2s_only")
	abline(h=summary(r_qtl_F2s.perm, alpha=c(0.10, 0.05, 0.01,0.001))[1],col="green")	
	abline(h=summary(r_qtl_F2s.perm, alpha=c(0.10, 0.05, 0.01,0.001))[2],col="red")
		
			
	dev.off()	
}           


###############################################################################################
#batch generate results
#remove outfile if previously created
file.remove("results/pheromones_ordination/r_qtl_approach_PCA_27-09-22_trim_tips/pheromones_ordination_r_qtl_approach_PCA_dat_27-09-22_trim_tips.txt")
header<-paste("Cross","Compound", "LOD_max","chrom_max", "physical_peak", "cM_max", "cM_limits", "physical_limits", "BO","B1","B2","B3_pgm","B1_B3_pgm_int","R2")
write(header,file="results/pheromones_ordination/r_qtl_approach_PCA_27-09-22_trim_tips/pheromones_ordination_r_qtl_approach_PCA_dat_27-09-22_trim_tips.txt",append=TRUE)

cross<-"r_qtl_F2s"
datset_to_analyse<-r_qtl_F2s

for (i in 10:160){
	pheno_name<-names(datset_to_analyse$pheno)[i]
	datset_to_analyse.scan<-scanone(datset_to_analyse,pheno.col=i,model="normal",method="hk")
	datset_to_analyse.perm<-scanone(datset_to_analyse,pheno.col=i,n.perm=n.perm,model="normal",method="hk") #perm.Xsp=TRUE	
	threshold<-summary(datset_to_analyse.perm, alpha=c(0.10, 0.05, 0.01,0.001))[2]
	#you could add separate threshold for sex here if you wanted. 
	for (chrom in unique(datset_to_analyse.scan$chr)){
		temp_chrom<-datset_to_analyse.scan[datset_to_analyse.scan$chr==chrom,]
		if(sum(temp_chrom[3])>0){		
			if(chrom!="X" & max(temp_chrom)[3]>threshold | chrom=="X" & max(temp_chrom)[3]>threshold){
				chrom_QTL<-as.character(max(temp_chrom)[1,1])
				
				#First get the genetic confidence intervals
				confint_<-print(bayesint(datset_to_analyse.scan,chrom))			
				LOD_max<-round(confint_[2,3],2)
				chrom_max<-confint_[2,1]
				cM_max<-round(confint_[2,2],2)
				cM_limits<-paste(round(confint_[1,2],2),"-",round(confint_[3,2],2),sep="")
				
				#now remove the inferred positions and get the physical confidence intervals
				temp_chrom_physical_only<-rownames(temp_chrom)
				temp_chrom_physical_only<-temp_chrom[!grepl("loc", rownames(temp_chrom), fixed=TRUE),]
				confint_<-print(bayesint(temp_chrom_physical_only))	
				mname1 <- find.marker(datset_to_analyse, confint_[2,1], confint_[2,2])
				png(filename = paste("results/pheromones_ordination/r_qtl_approach_PCA_27-09-22_trim_tips/significant_geno_pheno_plot/",cross,"_",pheno_name," chrom = ",chrom,".png",sep=""))
				plotPXG(datset_to_analyse,pheno.col=i,c(mname1))
				dev.off()
				lower_interval<-list_all_markers[list_all_markers$LG==confint_[1,1] & round(list_all_markers$cM,3)==round(confint_[1,2],3),] 
				peak<-list_all_markers[list_all_markers$LG==confint_[2,1] & round(list_all_markers$cM,3)==round(confint_[2,2],3),] 
				if(nrow(peak)==1){
					median_peak<-peak$pos
					median_scaff<-peak[which(as.numeric(peak$pos)==median_peak),2]
				}else{
					median_peak<-peak$pos[length(peak$pos)/2]
					median_scaff<-peak[which(as.numeric(peak$pos)==median_peak),2]
				}	
				upper_interval<-list_all_markers[list_all_markers$LG==confint_[3,1] & round(list_all_markers$cM,3)==round(confint_[3,2],3),] 
				physical_peak<-paste(median_scaff,":",median_peak,sep="")
				physical_limits<-paste(lower_interval[1,2],":",lower_interval[1,3],"-",upper_interval[nrow(upper_interval),2],":",upper_interval[nrow(upper_interval),3],sep="")
				
				#now model and extract coefficients
				pheno_test<-i
				pheno_for_lm<-pull.pheno(datset_to_analyse, pheno_test)
				covariate_pgm<-pull.pheno(datset_to_analyse, which(names(r_qtl_F2s$pheno)=="pgm"))
				covariate_sex<-pull.pheno(datset_to_analyse, which(names(r_qtl_F2s$pheno)=="Sex"))				
				geno_for_lm<-pull.genoprob(datset_to_analyse, chrom_QTL, omit.first.prob=FALSE, include.pos.info=TRUE, rotate=TRUE)
				geno_for_lm_1<-t(geno_for_lm[geno_for_lm$marker==mname1,][5:ncol(geno_for_lm)])
				geno_for_lm_1<-geno_for_lm_1[,-1] #sets EE as as the intercept
				if(chrom!="X"){				
					summary(lm(pheno_for_lm~geno_for_lm_1))
					model<-lm(pheno_for_lm~geno_for_lm_1)
					out<-data.frame(round(coef(summary(model))[, "Estimate"],2),round(coef(summary(model))[, "Std. Error"],2),round(coef(summary(model))[, "Pr(>|t|)"],3),round(coef(summary(model))[, "Pr(>|t|)"],4))
					colnames(out)<-c("coef","stderr","p")
					out$stars<-NA
					out$stars[which(out$p<0.1)]<-"^"
					out$stars[which(out$p<0.05)]<-"*"
					out$stars[which(out$p<0.01)]<-"**"
					out$stars[which(out$p<0.001)]<-"***"
					out$stars[which(out$p>0.1)]<-""
					out$p<-out$stars
					#print QTL
					output_2_write<-paste(cross,names(r_qtl_F2s$pheno)[i],LOD_max," ",chrom_max," ",physical_peak," ",cM_max," ",cM_limits," ",physical_limits,noquote(paste(out[1,1],"±",out[1,2],out[1,3]," ",out[2,1],"±",out[2,2],out[2,3]," ",out[3,1],"±",out[3,2],out[3,3]," ","NA"," ","NA"," ",round(summary(model)$r.squared,2),sep="")))
					write(output_2_write,file="results/pheromones_ordination/r_qtl_approach_PCA_27-09-22_trim_tips/pheromones_ordination_r_qtl_approach_PCA_dat_27-09-22_trim_tips.txt",append=TRUE)
				}else if(chrom=="X"){
					summary(lm(pheno_for_lm~geno_for_lm_1 * covariate_pgm))
					model<-lm(pheno_for_lm~geno_for_lm_1 * covariate_pgm)
					out<-data.frame(round(coef(summary(model))[, "Estimate"],2),round(coef(summary(model))[, "Std. Error"],2),round(coef(summary(model))[, "Pr(>|t|)"],3),round(coef(summary(model))[, "Pr(>|t|)"],4))
					colnames(out)<-c("coef","stderr","p")
					out$stars<-NA
					out$stars[which(out$p<0.1)]<-"^"
					out$stars[which(out$p<0.05)]<-"*"
					out$stars[which(out$p<0.01)]<-"**"
					out$stars[which(out$p<0.001)]<-"***"
					out$stars[which(out$p>0.1)]<-""
					out$p<-out$stars
					#print QTL
					output_2_write<-paste(cross,names(r_qtl_F2s$pheno)[i],LOD_max," ",chrom_max," ",physical_peak," ",cM_max," ",cM_limits," ",physical_limits,noquote(paste(out[1,1],"±",out[1,2],out[1,3]," ",out[2,1],"±",out[2,2],out[2,3]," ","NA"," ",out[3,1],"±",out[3,2],out[3,3]," ",out[4,1],"±",out[4,2],out[4,3]," ",round(summary(model)$r.squared,2),sep="")))					
					write(output_2_write,file="results/pheromones_ordination/r_qtl_approach_PCA_27-09-22_trim_tips/pheromones_ordination_r_qtl_approach_PCA_dat_27-09-22_trim_tips.txt",append=TRUE)
				}				
			}
		}	
	}           
}






