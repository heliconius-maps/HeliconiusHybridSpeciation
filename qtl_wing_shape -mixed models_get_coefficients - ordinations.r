
#This is to get the coefficients for the mixed effects models in a format suitable for pasting into excel
#Because they are slow to run, best is probably to use the cluster script to get the plots and see which compounds are significant, and then do them manually. 
#you will need to set the permutation threshold manually below, as well as the phenotype you want to extract details from


library(qtl)
library(doParallel)
library(lme4)

wingshape_custom<-read.csv(file="qtl_dat_wing_shape_lme_trim_tips.csv",header=FALSE)

FW_rotation_rows<-c(48:83)
HW_rotation_rows<-c(110:135)

row_start_genotypes<-136

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


#now output results

header<-paste("Cross","Pheno", "LOD_max","chrom_max", "physical_peak", "cM_max", "cM_limits", "physical_limits", "BO","B1","B2","B3_csize","R2","Var_ran_effect_CROSS","Var_ran_effect_SEX")

write(header,file="results/Feb 2023/wingshape_custom_cross_type_1000_perms_rotated/output_wingshape_lme_cross_type.txt",append=TRUE)

cross<-"wingshape_custom"
datset_to_analyse<-wingshape_custom


#load LOD significance thresholds
LOD_sig_thresholds<-read.table("C:/Users/NER334/Work/Harvard/Genomics/RAD_2018/elev_x_pard/QTL_F2s/qtl/results/Feb 2023/wingshape_custom_cross_type_1000_perms_rotated/sig_thresholds.txt",header=T)

for (i in FW_rotation_rows){	
	phenotype_row<-i
	pheno_name<-wingshape_custom[phenotype_row,1]
	threshold<-LOD_sig_thresholds[which(LOD_sig_thresholds[,1]==pheno_name),2]	
	
	#png(filename = paste("plots/cross_type_10000_perms/",pheno_name,".png",sep=""), width = 1000, height = 300)
	#par(mfrow=c(1,1))

	#run custom scripts
	all_lumped<-wingshape_custom[,-which(wingshape_custom[phenotype_row,]=="-")]
	pedigree<-all_lumped[1:row_start_genotypes-1,]
	genotypes<-all_lumped[row_start_genotypes:nrow(all_lumped),4:ncol(all_lumped)]
	markers<-all_lumped[row_start_genotypes:nrow(all_lumped),1:3]
	autosome_genotypes<-genotypes[which(markers$V2!="X"),]
	sex_genotypes<-genotypes[which(markers$V2=="X"),]	
	phenotype<-as.numeric(paste(pedigree[phenotype_row,4:ncol(pedigree)]))
	pgm<-as.factor(paste(pedigree[6,4:ncol(pedigree)]))
	cross_type<-as.factor(paste(pedigree[7,4:ncol(pedigree)]))	
	cross_type2<-as.factor(paste(pedigree[8,4:ncol(pedigree)]))
	families<-as.factor(paste(pedigree[1,4:ncol(pedigree)]))
	FW_csize<-as.numeric(paste(pedigree[10,4:ncol(pedigree)]))
	sex<-as.numeric(paste(pedigree[5,4:ncol(pedigree)]))
	lod_null_auto<-logLik(lmer(phenotype ~ FW_csize + (1|cross_type),REML=FALSE))	
	lod_null_sex <-logLik(lmer(phenotype ~ FW_csize + (1|cross_type) + (1|sex),REML=FALSE))		
	lod_auto <- apply(autosome_genotypes, 1, function(a) (logLik(lmer(phenotype ~ a + FW_csize + (1|cross_type),REML=FALSE))-lod_null_auto)/log(10))
	lod_sex <- apply(sex_genotypes, 1, function(a) (logLik(lmer(phenotype ~ a + FW_csize + (1|cross_type) + (1|sex),REML=FALSE))-lod_null_sex)/log(10))
	lod<-c(lod_auto,lod_sex)
	out<-cbind(markers[,2:3],lod)
	rownames(out)<-markers[,1]
	colnames(out)<-c("chr","pos","lod")

	class(out) <- c("scanone", "data.frame")

	datset_to_analyse.scan<-out

	for (chrom in unique(datset_to_analyse.scan$chr)){
			temp_chrom<-datset_to_analyse.scan[datset_to_analyse.scan$chr==chrom,]
			if(sum(temp_chrom[3])>0){
				if(max(temp_chrom)[3]>threshold){
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
					#mname1 <- find.marker(datset_to_analyse, confint_[2,1], confint_[2,2])
					#png(filename = paste("results/wing_shape/significant_geno_pheno_plot/lme",pheno_name," chrom = ",chrom,".png",sep=""))
					#plotPXG(datset_to_analyse,pheno.col=i,c(mname1))
					#dev.off()
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
					genotypes_chrom<-genotypes[datset_to_analyse.scan$chr==chrom,]						
					genotypes_for_model<-unlist(genotypes_chrom[which(temp_chrom$lod==max(temp_chrom$lod))[1],])
					if(chrom!="X"){
						model_summary<-summary(lmer(phenotype ~ genotypes_for_model + FW_csize + (1|cross_type),REML=FALSE))				
						intercept<-round(model_summary$coefficients[1,1],2)
						intercept_std_err<-round(model_summary$coefficients[1,2],2)
						intercept_t_val<-round(model_summary$coefficients[1,3],2)

						B1<-round(model_summary$coefficients[2,1],2)
						B1_std_err<-round(model_summary$coefficients[2,2],2)
						B1_t_val<-round(model_summary$coefficients[2,3],2)

						B2<-round(model_summary$coefficients[3,1],2)
						B2_std_err<-round(model_summary$coefficients[3,2],2)
						B2_t_val<-round(model_summary$coefficients[3,3],2)

						B3<-round(model_summary$coefficients[4,1],2)
						B3_std_err<-round(model_summary$coefficients[4,2],2)
						B3_t_val<-round(model_summary$coefficients[4,3],2)

						var_r.effect<-round(data.frame(model_summary$varcor)[4][1,1],2)
						#print QTL
						output_2_write<-print(paste("lme_r.crosstype_additive",pheno_name,LOD_max," ",chrom_max," ",physical_peak," ",cM_max," ",cM_limits," ",physical_limits,noquote(paste(intercept,"±",intercept_std_err," ",B1,"±",B1_std_err," ",B2,"±",B2_std_err," ",B3,"±",B3_std_err," ","NA"," ",var_r.effect," NA",sep=""))))
						write(output_2_write,file="results/Feb 2023/wingshape_custom_cross_type_1000_perms_rotated/output_wingshape_lme_cross_type.txt",append=TRUE)
					}else{
						model_summary<-summary(lmer(phenotype ~ genotypes_for_model + FW_csize + (1|cross_type) + (1|sex),REML=FALSE))				
						intercept<-round(model_summary$coefficients[1,1],2)
						intercept_std_err<-round(model_summary$coefficients[1,2],2)
						intercept_t_val<-round(model_summary$coefficients[1,3],2)

						B1<-round(model_summary$coefficients[2,1],2)
						B1_std_err<-round(model_summary$coefficients[2,2],2)
						B1_t_val<-round(model_summary$coefficients[2,3],2)

						B2<-round(model_summary$coefficients[3,1],2)
						B2_std_err<-round(model_summary$coefficients[3,2],2)
						B2_t_val<-round(model_summary$coefficients[3,3],2)

						B3<-round(model_summary$coefficients[4,1],2)
						B3_std_err<-round(model_summary$coefficients[4,2],2)
						B3_t_val<-round(model_summary$coefficients[4,3],2)
						
						var_r.effect1<-round(data.frame(model_summary$varcor)[4][1,1],2)
						var_r.effect2<-round(data.frame(model_summary$varcor)[4][2,1],2)
						#print QTL
						output_2_write<-print(paste("lme_r.crosstype_additive",pheno_name,LOD_max," ",chrom_max," ",physical_peak," ",cM_max," ",cM_limits," ",physical_limits,noquote(paste(intercept,"±",intercept_std_err," ",B1,"±",B1_std_err," ",B2,"±",B2_std_err," ",B3,"±",B3_std_err," ","NA"," ",var_r.effect1," ",var_r.effect2,sep=""))))
						write(output_2_write,file="results/Feb 2023/wingshape_custom_cross_type_1000_perms_rotated/output_wingshape_lme_cross_type.txt",append=TRUE)				
					}
				}
			}	
		}          
	}


for (i in HW_rotation_rows){	
	phenotype_row<-i
	phenotype_row<-i
	pheno_name<-wingshape_custom[phenotype_row,1]
	threshold<-LOD_sig_thresholds[which(LOD_sig_thresholds[,1]==pheno_name),2]	
	
	#png(filename = paste("plots/cross_type_10000_perms/",pheno_name,".png",sep=""), width = 1000, height = 300)
	#par(mfrow=c(1,1))

	#run custom scripts
	all_lumped<-wingshape_custom[,-which(wingshape_custom[phenotype_row,]=="-")]
	pedigree<-all_lumped[1:row_start_genotypes-1,]
	genotypes<-all_lumped[row_start_genotypes:nrow(all_lumped),4:ncol(all_lumped)]
	markers<-all_lumped[row_start_genotypes:nrow(all_lumped),1:3]
	autosome_genotypes<-genotypes[which(markers$V2!="X"),]
	sex_genotypes<-genotypes[which(markers$V2=="X"),]	
	phenotype<-as.numeric(paste(pedigree[phenotype_row,4:ncol(pedigree)]))
	pgm<-as.factor(paste(pedigree[6,4:ncol(pedigree)]))
	cross_type<-as.factor(paste(pedigree[7,4:ncol(pedigree)]))	
	cross_type2<-as.factor(paste(pedigree[8,4:ncol(pedigree)]))
	families<-as.factor(paste(pedigree[1,4:ncol(pedigree)]))
	HW_csize<-as.numeric(paste(pedigree[11,4:ncol(pedigree)]))
	sex<-as.numeric(paste(pedigree[5,4:ncol(pedigree)]))
	lod_null_auto<-logLik(lmer(phenotype ~ HW_csize + (1|cross_type),REML=FALSE))	
	lod_null_sex <-logLik(lmer(phenotype ~ HW_csize + (1|cross_type) + (1|sex),REML=FALSE))		
	lod_auto <- apply(autosome_genotypes, 1, function(a) (logLik(lmer(phenotype ~ a + HW_csize + (1|cross_type),REML=FALSE))-lod_null_auto)/log(10))
	lod_sex <- apply(sex_genotypes, 1, function(a) (logLik(lmer(phenotype ~ a + HW_csize + (1|cross_type) + (1|sex),REML=FALSE))-lod_null_sex)/log(10))
	lod<-c(lod_auto,lod_sex)
	out<-cbind(markers[,2:3],lod)
	rownames(out)<-markers[,1]
	colnames(out)<-c("chr","pos","lod")

	class(out) <- c("scanone", "data.frame")

	datset_to_analyse.scan<-out
	for (chrom in unique(datset_to_analyse.scan$chr)){
			temp_chrom<-datset_to_analyse.scan[datset_to_analyse.scan$chr==chrom,]
			if(sum(temp_chrom[3])>0){
				if(max(temp_chrom)[3]>threshold){
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
					#mname1 <- find.marker(datset_to_analyse, confint_[2,1], confint_[2,2])
					#png(filename = paste("results/wing_shape/significant_geno_pheno_plot/lme",pheno_name," chrom = ",chrom,".png",sep=""))
					#plotPXG(datset_to_analyse,pheno.col=i,c(mname1))
					#dev.off()
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
					genotypes_chrom<-genotypes[datset_to_analyse.scan$chr==chrom,]
						
					genotypes_for_model<-unlist(genotypes_chrom[which(temp_chrom$lod==max(temp_chrom$lod))[1],])
					if(chrom!="X"){
						model_summary<-summary(lmer(phenotype ~ genotypes_for_model + HW_csize + (1|cross_type),REML=FALSE))				
						intercept<-round(model_summary$coefficients[1,1],2)
						intercept_std_err<-round(model_summary$coefficients[1,2],2)
						intercept_t_val<-round(model_summary$coefficients[1,3],2)

						B1<-round(model_summary$coefficients[2,1],2)
						B1_std_err<-round(model_summary$coefficients[2,2],2)
						B1_t_val<-round(model_summary$coefficients[2,3],2)

						B2<-round(model_summary$coefficients[3,1],2)
						B2_std_err<-round(model_summary$coefficients[3,2],2)
						B2_t_val<-round(model_summary$coefficients[3,3],2)

						B3<-round(model_summary$coefficients[4,1],2)
						B3_std_err<-round(model_summary$coefficients[4,2],2)
						B3_t_val<-round(model_summary$coefficients[4,3],2)

						var_r.effect<-round(data.frame(model_summary$varcor)[4][1,1],2)
						#print QTL
						output_2_write<-print(paste("lme_r.crosstype_additive",pheno_name,LOD_max," ",chrom_max," ",physical_peak," ",cM_max," ",cM_limits," ",physical_limits,noquote(paste(intercept,"±",intercept_std_err," ",B1,"±",B1_std_err," ",B2,"±",B2_std_err," ",B3,"±",B3_std_err," ","NA"," ",var_r.effect," NA",sep=""))))
						write(output_2_write,file="results/Feb 2023/wingshape_custom_cross_type_1000_perms_rotated/output_wingshape_lme_cross_type.txt",append=TRUE)
					}else{
						model_summary<-summary(lmer(phenotype ~ genotypes_for_model + HW_csize + (1|cross_type) + (1|sex),REML=FALSE))				
						intercept<-round(model_summary$coefficients[1,1],2)
						intercept_std_err<-round(model_summary$coefficients[1,2],2)
						intercept_t_val<-round(model_summary$coefficients[1,3],2)

						B1<-round(model_summary$coefficients[2,1],2)
						B1_std_err<-round(model_summary$coefficients[2,2],2)
						B1_t_val<-round(model_summary$coefficients[2,3],2)

						B2<-round(model_summary$coefficients[3,1],2)
						B2_std_err<-round(model_summary$coefficients[3,2],2)
						B2_t_val<-round(model_summary$coefficients[3,3],2)

						B3<-round(model_summary$coefficients[4,1],2)
						B3_std_err<-round(model_summary$coefficients[4,2],2)
						B3_t_val<-round(model_summary$coefficients[4,3],2)
						
						var_r.effect1<-round(data.frame(model_summary$varcor)[4][1,1],2)
						var_r.effect2<-round(data.frame(model_summary$varcor)[4][2,1],2)
						#print QTL
						output_2_write<-print(paste("lme_r.crosstype_additive",pheno_name,LOD_max," ",chrom_max," ",physical_peak," ",cM_max," ",cM_limits," ",physical_limits,noquote(paste(intercept,"±",intercept_std_err," ",B1,"±",B1_std_err," ",B2,"±",B2_std_err," ",B3,"±",B3_std_err," ","NA"," ",var_r.effect1," ",var_r.effect2,sep=""))))
						write(output_2_write,file="results/Feb 2023/wingshape_custom_cross_type_1000_perms_rotated/output_wingshape_lme_cross_type.txt",append=TRUE)														
					}	
				}
			}	
		}          
	}

