#!/usr/bin/env Rscript 

#./qtl_colour_pattern_PCA_lme_FW_cluster.r & #Better not to put an out file because it will be huge (lme spits a lot of messages during all the perms)
	
library(qtl)
library(doParallel)
library(lme4)

options(bitmapType='cairo')
print(version)
#permutations for single locus scan
scanone.lm.perm <-
	function(genotypes, phenotype,cross_type,n.perm=10,sex,sex.perms=TRUE,thresholds=c(0.1,0.05,0.01)){				
		result <- NULL
		registerDoParallel(cores=4)			
		perm.list<-foreach(perm.it=1:n.perm) %dopar%  {
			result.temp<-NULL
			print(paste("permutation = ", perm.it,sep=""))	
			rand_index<-sample(1:length(phenotype))	
			pheno.randomised<-phenotype[rand_index]
			#pgm.randomised<-pgm[rand_index]
			sex.randomised<-sex[rand_index]			
			cross_type.randomised<-cross_type[rand_index]
			library(lme4)
			if (sex.perms==FALSE){
				lod_null<-logLik(lmer(pheno.randomised ~ (1|cross_type.randomised),REML=FALSE))								
				lod <- apply(genotypes, 1, function(a) (logLik(lmer(pheno.randomised ~ a + (1|cross_type.randomised),REML=FALSE))-lod_null)/log(10))
			}else{
				lod_null_sex<-logLik(lmer(pheno.randomised ~ (1|sex.randomised) + (1|cross_type.randomised),REML=FALSE))
				lod <- apply(genotypes, 1, function(a) (logLik(lmer(pheno.randomised ~ a + (1|sex.randomised) + (1|cross_type.randomised),REML=FALSE))-lod_null_sex)/log(10))
			}			
			max(lod)
		}
		stopImplicitCluster()
		result<-unlist(perm.list)	
		result<-sort(result)			
		result_p0.1<-result[round((1-thresholds[1])*length(result))]
		result_p0.05<-result[round((1-thresholds[2])*length(result))]
		result_p0.01<-result[round((1-thresholds[3])*length(result))]
		result_out<-data.frame(cbind(result_p0.1,result_p0.05,result_p0.01))
		result_out
}	



qtl_dat<-read.csv(file="qtl_dat_NR_all_FWs_rm_bg_aligned_NR16_29_white_95%_unconst_PCA_FINAL.csv",header=FALSE)

#this is every PC until they stop explaining >0.5% of the variance. 
phenotype_rows<-c(17:49)

row_start_genotypes<-330


#Permutations
n.perm<-1000
chom_lengths<-c()
for (i in unique(qtl_dat$V2[row_start_genotypes:nrow(qtl_dat)])){
	chom_lengths<-append(chom_lengths,max(qtl_dat[qtl_dat$V2==i,]$V3))
}

n.permX<-round(n.perm*(sum(chom_lengths[1:20])/chom_lengths[21])) #For the X chromosome, n.perm × L_A/L_X permutations are performed 
#false positive alpha correction for separate sex chromosome permutations 
L_A<-sum(chom_lengths[1:20])
L_X<-sum(chom_lengths[21])
L<-L_A+L_X
LOD_A_0.1<-1-(1-0.1)^(L_A/L)
LOD_X_0.1<-1-(1-0.1)^(L_X/L)
LOD_A_0.05<-1-(1-0.05)^(L_A/L)
LOD_X_0.05<-1-(1-0.05)^(L_X/L)
LOD_A_0.01<-1-(1-0.01)^(L_A/L)
LOD_X_0.01<-1-(1-0.01)^(L_X/L)

#Use these when you plot a separate sex chrom sig threshold
#These x values are specific to this particular map and not generalisable
# autosomes_x_coords_4_plot<-c(0,2165)
# sex_x_coords_4_plot<-c(2190,2244)
#for trimmed datafile
autosomes_x_coords_4_plot<-c(0,1518)
sex_x_coords_4_plot<-c(1519,2000)


header<-paste("phenotype","auto_5%_threshold","sex_5%_threshold")
write(header,file="results/colour_pattern_PCA/NR_all_FWs_rm_bg_aligned_NR16_29_white_95_unconst_PCA_FINAL_lme/sig_thresholds_custom_lme_colour_PCA_FW.txt",append=TRUE)

sig_thresholds<-c()
for (i in phenotype_rows){
	phenotype_row<-i
	pheno_name<-print(qtl_dat[phenotype_row,1])
	png(filename = paste("results/colour_pattern_PCA/NR_all_FWs_rm_bg_aligned_NR16_29_white_95_unconst_PCA_FINAL_lme/",pheno_name,".png",sep=""), width = 1000, height = 300)
	par(mfrow=c(1,1))
	#run standard r/qtl
	#run custom scripts
	all_lumped<-qtl_dat[,-which(qtl_dat[phenotype_row,]=="-")]
	pedigree<-all_lumped[1:row_start_genotypes-1,]
	chromosomes<-all_lumped[row_start_genotypes:nrow(all_lumped),2]
	genotypes<-all_lumped[row_start_genotypes:nrow(all_lumped),4:ncol(all_lumped)]
	genotypes_autosomes<-genotypes[which(chromosomes!="X"),]
	genotypes_sex<-genotypes[which(chromosomes=="X"),]
	markers<-all_lumped[row_start_genotypes:nrow(all_lumped),1:3]
	phenotype<-as.numeric(paste(pedigree[phenotype_row,4:ncol(pedigree)]))
	pgm<-as.factor(paste(pedigree[6,4:ncol(pedigree)]))
	cross_type<-as.factor(paste(pedigree[7,4:ncol(pedigree)]))	
	families<-as.factor(paste(pedigree[1,4:ncol(pedigree)]))
	sex<-as.numeric(paste(pedigree[5,4:ncol(pedigree)]))		
	lod_null<-logLik(lmer(phenotype ~ (1|cross_type),REML=FALSE))
	lod_null_sex<-logLik(lmer(phenotype ~ (1|sex) + (1|cross_type),REML=FALSE)) 
	lod_auto <- apply(genotypes_autosomes, 1, function(a) (logLik(lmer(phenotype ~ a + (1|cross_type),REML=FALSE))-lod_null)/log(10))		
	lod_sex <- apply(genotypes_sex, 1, function(a) (logLik(lmer(phenotype ~ a + (1|sex) + (1|cross_type),REML=FALSE))-lod_null_sex)/log(10))
	lod<-c(lod_auto,lod_sex)
	out<-cbind(markers[,2:3],lod)
	rownames(out)<-markers[,1]
	colnames(out)<-c("chr","pos","lod")
	#bind sex markers with no values so all plots line up for easy comparison
	#out<-rbind(out,dummy_sex)
	class(out) <- c("scanone", "data.frame")
	threshold_autosomes<-scanone.lm.perm(genotypes=genotypes_autosomes,phenotype=phenotype,cross_type=cross_type,n.perm=n.perm,sex=sex,sex.perms=FALSE,thresholds=c(LOD_A_0.1,LOD_A_0.05,LOD_A_0.01))
	threshold_sex<-scanone.lm.perm(genotypes_sex,phenotype=phenotype,cross_type=cross_type,n.perm=n.permX,sex=sex,sex.perms=TRUE,thresholds=c(LOD_X_0.1,LOD_X_0.05,LOD_X_0.01))
	plot(out,main=paste("custom with cross_type, 5% sig thresh auto = ",round(threshold_autosomes[1,2],3), "5% sig thresh sex = ",round(threshold_sex[1,2],3),sep=""))
	segments(autosomes_x_coords_4_plot[1], threshold_autosomes[1,1], autosomes_x_coords_4_plot[2], threshold_autosomes[1,1], col= 'green')
	segments(autosomes_x_coords_4_plot[1], threshold_autosomes[1,2], autosomes_x_coords_4_plot[2], threshold_autosomes[1,2], col= 'red')
	segments(sex_x_coords_4_plot[1], threshold_sex[1,1], sex_x_coords_4_plot[2], threshold_sex[1,1], col= 'green')
	segments(sex_x_coords_4_plot[1], threshold_sex[1,2], sex_x_coords_4_plot[2], threshold_sex[1,2], col= 'red')
	dev.off()
	sig_thresholds<-rbind(sig_thresholds,c(pheno_name,threshold_autosomes[1,2],threshold_sex[1,2]))
	sig_thresholds2<-paste(pheno_name,paste(threshold_autosomes[1,2],threshold_sex[1,2]))
	write(sig_thresholds2,file="results/colour_pattern_PCA/NR_all_FWs_rm_bg_aligned_NR16_29_white_95_unconst_PCA_FINAL_lme/sig_thresholds_custom_lme_colour_PCA_FW.txt",append=TRUE)
}      

     
colnames(sig_thresholds)<-c("phenotype","5%_threshold_auto","5%_threshold_sex")
write.csv(sig_thresholds,"results/colour_pattern_PCA/colour_pattern_PCA/NR_all_FWs_rm_bg_aligned_NR16_29_white_95_unconst_PCA_FINAL_lme/sig_thresholds_custom_lme_colour_PCA_FW.csv")



