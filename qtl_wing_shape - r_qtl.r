#!/usr/bin/env Rscript 

#./qtl_wing_shape_r_qtl.r &  #> qtl_wing_shape_r_qtl_OUT.log &
	
library(qtl)
options(bitmapType='cairo')
print(version)


r_qtl_F2s<-read.cross("csvr",file="qtl_dat_wing_shape_F2s_trim_tips.csv",genotypes=c("EE","PE","PP"),alleles=c("E","P"))
r_qtl_F2s<-jittermap(r_qtl_F2s)
r_qtl_F2s<-calc.genoprob(r_qtl_F2s,step=1, off.end=0.0, error.prob=0.001,map.function="haldane",stepwidth="fixed")



#markers for confidence intervals
list_all_markers<-read.csv("C:/Users/NER334/Work/Harvard/Genomics/RAD_2018/elev_x_pard/QTL_F2s/qtl/Fst/list_all_markers_Fst_trim_tips_reset_cM.csv",header=T,sep=",")
list_all_markers$X<-NULL
list_all_markers$Fst<-NULL
list_all_markers$LG[which(list_all_markers$LG=="21")]<-"X"


# list_all_markers<-read.table("C:/Users/NER334/Work/Harvard/Genomics/RAD_2018/elev_x_pard/linkage_maps/all_markers")
# colnames(list_all_markers)<-as.character(unlist(list_all_markers[1,]))
# list_all_markers<-list_all_markers[-1,]
# list_all_markers$cM<-as.numeric(as.character(list_all_markers$cM))
# list_all_markers$pos<-as.numeric(as.character(list_all_markers$pos))
# list_all_markers$scaff<-as.character(list_all_markers$scaff)
# list_all_markers$LG<-as.character(list_all_markers$LG)
# list_all_markers$LG[which(list_all_markers$LG=="21")]<-"X"

n.perm<-1000



#Need to do the FW and HW separately cos they take a different covariate (centroid size)

which(names(r_qtl_F2s$pheno)=="FW_1")
which(names(r_qtl_F2s$pheno)=="FW_rot_36")
which(names(r_qtl_F2s$pheno)=="HW_1")
which(names(r_qtl_F2s$pheno)=="HW_rot_26")

#Use these when you plot a separate sex chrom sig threshold
#These x values are specific to this particular map and not generalisable
# autosomes_x_coords_4_plot<-c(0,2165)
# sex_x_coords_4_plot<-c(2190,2244)
#for trimmed datafile
autosomes_x_coords_4_plot<-c(0,1518)
sex_x_coords_4_plot<-c(1519,2000)


for (i in 11:82){
	pheno_name<-names(r_qtl_F2s$pheno)[i]
	png(filename = paste("results/wing_shape/r_qtl_approach_trimmed/",pheno_name,".png",sep=""), width = 1000, height = 300)
	#par(mfrow=c(4,1))
	
	r_qtl_F2s.scan<-scanone(r_qtl_F2s,pheno.col=i,model="normal",method="hk",addcovar=r_qtl_F2s$pheno[,which(names(r_qtl_F2s$pheno)=="FW_csize")])
	r_qtl_F2s.perm<-scanone(r_qtl_F2s,pheno.col=i,n.perm=n.perm,model="normal",method="hk",addcovar=r_qtl_F2s$pheno[,which(names(r_qtl_F2s$pheno)=="FW_csize")],perm.Xsp=TRUE)
	plot(r_qtl_F2s.scan,main="r_qtl_F2s")
	auto_thresh_10<-summary(r_qtl_F2s.perm, alpha=c(0.10, 0.05, 0.01,0.001))$A[1]
	auto_thresh_5<-summary(r_qtl_F2s.perm, alpha=c(0.10, 0.05, 0.01,0.001))$A[2]
	sex_thresh_10<-summary(r_qtl_F2s.perm, alpha=c(0.10, 0.05, 0.01,0.001))$X[1]
	sex_thresh_5<-summary(r_qtl_F2s.perm, alpha=c(0.10, 0.05, 0.01,0.001))$X[2]	
	segments(autosomes_x_coords_4_plot[1], auto_thresh_10, autosomes_x_coords_4_plot[2], auto_thresh_10, col= 'green')
	segments(autosomes_x_coords_4_plot[1], auto_thresh_5, autosomes_x_coords_4_plot[2], auto_thresh_5, col= 'red')
	segments(sex_x_coords_4_plot[1], sex_thresh_10, sex_x_coords_4_plot[2], sex_thresh_10, col= 'green')
	segments(sex_x_coords_4_plot[1], sex_thresh_5, sex_x_coords_4_plot[2], sex_thresh_5, col= 'red')
		

	dev.off()	
}           

for (i in 83:134){
	pheno_name<-names(r_qtl_F2s$pheno)[i]
	png(filename = paste("results/wing_shape/r_qtl_approach_trimmed/",pheno_name,".png",sep=""), width = 1000, height = 300)
	#par(mfrow=c(4,1))
	
		
	r_qtl_F2s.scan<-scanone(r_qtl_F2s,pheno.col=i,model="normal",method="hk",addcovar=r_qtl_F2s$pheno[,which(names(r_qtl_F2s$pheno)=="HW_csize")])
	r_qtl_F2s.perm<-scanone(r_qtl_F2s,pheno.col=i,n.perm=n.perm,model="normal",method="hk",addcovar=r_qtl_F2s$pheno[,which(names(r_qtl_F2s$pheno)=="HW_csize")],perm.Xsp=TRUE)
	plot(r_qtl_F2s.scan,main="r_qtl_F2s")
	auto_thresh_10<-summary(r_qtl_F2s.perm, alpha=c(0.10, 0.05, 0.01,0.001))$A[1]
	auto_thresh_5<-summary(r_qtl_F2s.perm, alpha=c(0.10, 0.05, 0.01,0.001))$A[2]
	sex_thresh_10<-summary(r_qtl_F2s.perm, alpha=c(0.10, 0.05, 0.01,0.001))$X[1]
	sex_thresh_5<-summary(r_qtl_F2s.perm, alpha=c(0.10, 0.05, 0.01,0.001))$X[2]	
	segments(autosomes_x_coords_4_plot[1], auto_thresh_10, autosomes_x_coords_4_plot[2], auto_thresh_10, col= 'green')
	segments(autosomes_x_coords_4_plot[1], auto_thresh_5, autosomes_x_coords_4_plot[2], auto_thresh_5, col= 'red')
	segments(sex_x_coords_4_plot[1], sex_thresh_10, sex_x_coords_4_plot[2], sex_thresh_10, col= 'green')
	segments(sex_x_coords_4_plot[1], sex_thresh_5, sex_x_coords_4_plot[2], sex_thresh_5, col= 'red')

	dev.off()	
}     




###############################################################################################
#batch generate results
#remove outfile if previously created
file.remove("results/wing_shape/r_qtl_approach_trimmed/output_wingshape_trimmed.txt")
header<-paste("Cross","Compound", "LOD_max","chrom_max", "physical_peak", "cM_max", "cM_limits", "physical_limits", "BO","B1","B2","B3_csize","B4_sex","B4_pgm","B5_q1_sex","B5_q1_pgm","R2")
write(header,file="results/wing_shape/r_qtl_approach_trimmed/output_wingshape_trimmed.txt",append=TRUE)

cross<-"r_qtl_F2s"
datset_to_analyse<-r_qtl_F2s

#FW
for (i in 11:82){
	pheno_name<-names(datset_to_analyse$pheno)[i]
	datset_to_analyse.scan<-scanone(datset_to_analyse,pheno.col=i,model="normal",method="hk",addcovar=datset_to_analyse$pheno[,which(names(datset_to_analyse$pheno)=="FW_csize")])
	datset_to_analyse.perm<-scanone(datset_to_analyse,pheno.col=i,n.perm=n.perm,model="normal",method="hk",addcovar=datset_to_analyse$pheno[,which(names(datset_to_analyse$pheno)=="FW_csize")],perm.Xsp=TRUE)	
	threshold_autosomes<-summary(datset_to_analyse.perm, alpha=c(0.10, 0.05, 0.01,0.001))$A[2]
	threshold_sex<-summary(datset_to_analyse.perm, alpha=c(0.10, 0.05, 0.01,0.001))$X[2]
	for (chrom in unique(datset_to_analyse.scan$chr)){
		temp_chrom<-datset_to_analyse.scan[datset_to_analyse.scan$chr==chrom,]
		if(sum(temp_chrom[3])>0){
			if(chrom!="X" & max(temp_chrom)[3]>threshold_autosomes | chrom=="X" & max(temp_chrom)[3]>threshold_sex){
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
				png(filename = paste("results/wing_shape/r_qtl_approach_trimmed/significant_geno_pheno_plot/",pheno_name," chrom = ",chrom,".png",sep=""))
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
				covariate_csize<-pull.pheno(datset_to_analyse, which(names(datset_to_analyse$pheno)=="FW_csize"))
				covariate_pgm<-pull.pheno(datset_to_analyse, which(names(datset_to_analyse$pheno)=="pgm"))
				covariate_sex<-pull.pheno(datset_to_analyse, which(names(datset_to_analyse$pheno)=="Sex"))
				geno_for_lm<-pull.genoprob(datset_to_analyse, chrom_QTL, omit.first.prob=FALSE, include.pos.info=TRUE, rotate=TRUE)
				geno_for_lm_1<-t(geno_for_lm[geno_for_lm$marker==mname1,][5:ncol(geno_for_lm)])
				geno_for_lm_1<-geno_for_lm_1[,-1] #sets EE as as the intercept
				if(chrom!="X"){				
					summary(lm(pheno_for_lm~geno_for_lm_1 + covariate_csize))
					model<-lm(pheno_for_lm~geno_for_lm_1 + covariate_csize)
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
					output_2_write<-paste(cross,names(datset_to_analyse$pheno)[i],LOD_max," ",chrom_max," ",physical_peak," ",cM_max," ",cM_limits," ",physical_limits,noquote(paste(out[1,1],"±",out[1,2],out[1,3]," ",out[2,1],"±",out[2,2],out[2,3]," ",out[3,1],"±",out[3,2],out[3,3]," ",out[4,1],"±",out[4,2],out[4,3]," ","NA"," ","NA"," ","NA"," ","NA"," ",round(summary(model)$r.squared,2),sep="")))
					write(output_2_write,file="results/wing_shape/r_qtl_approach_trimmed/output_wingshape_trimmed.txt",append=TRUE)
				}else if(chrom=="X"){
					summary(lm(pheno_for_lm~geno_for_lm_1 * covariate_sex + geno_for_lm_1 * covariate_pgm + covariate_csize))
					model<-lm(pheno_for_lm~geno_for_lm_1 * covariate_sex + geno_for_lm_1 * covariate_pgm + covariate_csize)
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
					output_2_write<-paste(cross,names(datset_to_analyse$pheno)[i],LOD_max," ",chrom_max," ",physical_peak," ",cM_max," ",cM_limits," ",physical_limits,noquote(paste(out[1,1],"±",out[1,2],out[1,3]," ",out[2,1],"±",out[2,2],out[2,3]," ","NA"," ",out[5,1],"±",out[5,2],out[5,3]," ",out[3,1],"±",out[3,2],out[3,3]," ",out[4,1],"±",out[4,2],out[4,3]," ",out[6,1],"±",out[6,2],out[6,3]," ",out[7,1],"±",out[7,2],out[7,3]," ",round(summary(model)$r.squared,2),sep="")))
					write(output_2_write,file="results/wing_shape/r_qtl_approach_trimmed/output_wingshape_trimmed.txt",append=TRUE)
				}					
			}
		}	
	}           
}
#HW
for (i in 83:134){
	pheno_name<-names(datset_to_analyse$pheno)[i]
	datset_to_analyse.scan<-scanone(datset_to_analyse,pheno.col=i,model="normal",method="hk",addcovar=datset_to_analyse$pheno[,which(names(datset_to_analyse$pheno)=="HW_csize")])
	datset_to_analyse.perm<-scanone(datset_to_analyse,pheno.col=i,n.perm=n.perm,model="normal",method="hk",addcovar=datset_to_analyse$pheno[,which(names(datset_to_analyse$pheno)=="HW_csize")],perm.Xsp=TRUE)	
	threshold_autosomes<-summary(datset_to_analyse.perm, alpha=c(0.10, 0.05, 0.01,0.001))$A[2]
	threshold_sex<-summary(datset_to_analyse.perm, alpha=c(0.10, 0.05, 0.01,0.001))$X[2]
	for (chrom in unique(datset_to_analyse.scan$chr)){
		temp_chrom<-datset_to_analyse.scan[datset_to_analyse.scan$chr==chrom,]
		if(chrom!="X" & max(temp_chrom)[3]>threshold_autosomes | chrom=="X" & max(temp_chrom)[3]>threshold_sex){
			if(max(temp_chrom)[3]>threshold_autosomes){
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
				png(filename = paste("results/wing_shape/r_qtl_approach_trimmed/significant_geno_pheno_plot/",pheno_name," chrom = ",chrom,".png",sep=""))
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
				covariate_csize<-pull.pheno(datset_to_analyse, which(names(datset_to_analyse$pheno)=="HW_csize"))
				covariate_pgm<-pull.pheno(datset_to_analyse, which(names(datset_to_analyse$pheno)=="pgm"))
				covariate_sex<-pull.pheno(datset_to_analyse, which(names(datset_to_analyse$pheno)=="Sex"))
				geno_for_lm<-pull.genoprob(datset_to_analyse, chrom_QTL, omit.first.prob=FALSE, include.pos.info=TRUE, rotate=TRUE)
				geno_for_lm_1<-t(geno_for_lm[geno_for_lm$marker==mname1,][5:ncol(geno_for_lm)])
				geno_for_lm_1<-geno_for_lm_1[,-1] #sets EE as as the intercept
				if(chrom!="X"){				
					summary(lm(pheno_for_lm~geno_for_lm_1 + covariate_csize))
					model<-lm(pheno_for_lm~geno_for_lm_1 + covariate_csize)
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
					output_2_write<-paste(cross,names(datset_to_analyse$pheno)[i],LOD_max," ",chrom_max," ",physical_peak," ",cM_max," ",cM_limits," ",physical_limits,noquote(paste(out[1,1],"±",out[1,2],out[1,3]," ",out[2,1],"±",out[2,2],out[2,3]," ",out[3,1],"±",out[3,2],out[3,3]," ",out[4,1],"±",out[4,2],out[4,3]," ","NA"," ","NA"," ","NA"," ","NA"," ",round(summary(model)$r.squared,2),sep="")))
					write(output_2_write,file="results/wing_shape/r_qtl_approach_trimmed/output_wingshape_trimmed.txt",append=TRUE)
				}else if(chrom=="X"){
					summary(lm(pheno_for_lm~geno_for_lm_1 * covariate_sex + geno_for_lm_1 * covariate_pgm + covariate_csize))
					model<-lm(pheno_for_lm~geno_for_lm_1 * covariate_sex + geno_for_lm_1 * covariate_pgm + covariate_csize)
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
					output_2_write<-paste(cross,names(datset_to_analyse$pheno)[i],LOD_max," ",chrom_max," ",physical_peak," ",cM_max," ",cM_limits," ",physical_limits,noquote(paste(out[1,1],"±",out[1,2],out[1,3]," ",out[2,1],"±",out[2,2],out[2,3]," ","NA"," ",out[5,1],"±",out[5,2],out[5,3]," ",out[3,1],"±",out[3,2],out[3,3]," ",out[4,1],"±",out[4,2],out[4,3]," ",out[6,1],"±",out[6,2],out[6,3]," ",out[7,1],"±",out[7,2],out[7,3]," ",round(summary(model)$r.squared,2),sep="")))
					write(output_2_write,file="results/wing_shape/r_qtl_approach_trimmed/output_wingshape_trimmed.txt",append=TRUE)
				}	
			}
		}	
	}           
}





