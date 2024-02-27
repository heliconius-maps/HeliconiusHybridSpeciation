#!/usr/bin/env Rscript 

#./qtl_pheromones_lme_cluster_ordinations_PCA.r & #Better not to put an out file because it will be huge (lme spits a lot of messages during all the perms)


	
library(qtl)
library(doParallel)
library(lme4)

options(bitmapType='cairo')
print(version)
#permutations for single locus scan
scanone.lm.perm <-
	function(genotypes, phenotype,cross_type,n.perm=10){				
		result <- NULL
		registerDoParallel(cores=30)
			
		perm.list<-foreach(perm.it=1:n.perm) %dopar%  {
			result.temp<-NULL
			print(paste("permutation = ", perm.it,sep=""))	
			rand_index<-sample(1:length(phenotype))				
			pheno.randomised<-phenotype[rand_index]
			cross_type.randomised<-cross_type[rand_index]
			library(lme4)
			lod_null<-logLik(lmer(pheno.randomised ~ (1|cross_type.randomised),REML=FALSE))
			lod <- apply(genotypes, 1, function(a) (logLik(lmer(pheno.randomised ~ a + (1|cross_type.randomised),REML=FALSE))-lod_null)/log(10))
			result.temp <- append(result.temp, as.vector(lod))		
			max(result.temp)
		}
		stopImplicitCluster()
		result<-unlist(perm.list)	
		result<-sort(result)			
		result_p0.1<-result[0.9*length(result)]
		result_p0.05<-result[0.95*length(result)]
		result_p0.01<-result[0.99*length(result)]
		result_out<-data.frame(cbind(result_p0.1,result_p0.05,result_p0.01))
		result_out
}	



pheromones_custom<-read.csv(file="qtl_dat_pheromones_ordinations_PCA_F2s_BCs_27-09-22_trim_tips.csv",header=FALSE)

#If you want to analyse only autosomes you can use this, but be sure to also uncomment the bit in the loop below (#out<-rbind(out,dummy_sex))
#remove sex chromosome but first create dummy scanone for it (see below)
# dummy_sex<-pheromones_custom[which(pheromones_custom$V2=="X"),1:3]
# dummy_sex$lod<-rep(0,nrow(dummy_sex))
# rownames(dummy_sex)<-dummy_sex[,1]
# dummy_sex[,1]<-NULL
# colnames(dummy_sex)<-c("chr","pos","lod")
# pheromones_custom<-pheromones_custom[which(pheromones_custom$V2!="X"),]


row_start_phenotypes<-161
row_stop_phenotypes<-185	
row_start_genotypes<-186


header<-paste("phenotype","5%_threshold")
write(header,file="qtl_dat_pheromones_ordinations_PCA_F2s_BCs_27-09-22_trim_tips/sig_thresholds_trimmed.txt",append=TRUE)


sig_thresholds<-c()
for (i in row_start_phenotypes:row_stop_phenotypes){
	phenotype_row<-i
	pheno_name<-pheromones_custom[phenotype_row,1]
	png(filename = paste("qtl_dat_pheromones_ordinations_PCA_F2s_BCs_27-09-22_trim_tips/",pheno_name,"_trimmed.png",sep=""), width = 1000, height = 300)
	par(mfrow=c(1,1))
	#run standard r/qtl
	#run custom scripts
	all_lumped<-pheromones_custom[,-which(pheromones_custom[phenotype_row,]=="-")]
	pedigree<-all_lumped[1:row_start_genotypes-1,]
	genotypes<-all_lumped[row_start_genotypes:nrow(all_lumped),4:ncol(all_lumped)]
	markers<-all_lumped[row_start_genotypes:nrow(all_lumped),1:3]
	phenotype<-as.numeric(paste(pedigree[phenotype_row,4:ncol(pedigree)]))
	pgm<-as.factor(paste(pedigree[6,4:ncol(pedigree)]))
	cross_type<-as.factor(paste(pedigree[7,4:ncol(pedigree)]))	
	cross_type2<-as.factor(paste(pedigree[8,4:ncol(pedigree)]))
	families<-as.factor(paste(pedigree[1,4:ncol(pedigree)]))
	
	lod_null<-logLik(lmer(phenotype ~ (1|cross_type),REML=FALSE))	
	lod <- apply(genotypes, 1, function(a) (logLik(lmer(phenotype ~ a + (1|cross_type),REML=FALSE))-lod_null)/log(10))		
	out<-cbind(markers[,2:3],lod)
	rownames(out)<-markers[,1]
	colnames(out)<-c("chr","pos","lod")
	#bind sex markers with no values so all plots line up for easy comparison
	#out<-rbind(out,dummy_sex)
	class(out) <- c("scanone", "data.frame")
	out.perm<-scanone.lm.perm(genotypes,phenotype,cross_type,1000)
	plot(out,main=paste("custom with cross_type, 5% sig thresh = ",round(out.perm[2],3),sep=""))
	abline(h=out.perm[1],col="green")
	abline(h=out.perm[2],col="red")	
	abline(h=out.perm[3],col="red",lty=3)
	dev.off()
	sig_thresholds<-paste(pheno_name,paste(out.perm[2]))
	write(sig_thresholds,file="qtl_dat_pheromones_ordinations_PCA_F2s_BCs_27-09-22_trim_tips/sig_thresholds_trimmed.txt",append=TRUE)
}           
