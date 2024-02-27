#!/usr/bin/env Rscript 

#./qtl_male_preference_glmer_cluster.r &> qtl_male_preference_glmer_cluster.log & 


#I had a lot of problems running glmer in parallel. It seems that it does some internal parallellising of some kind, and basically distributed everything over as many nodes as it can, making the cluster hopelessly overloaded. 
#In order to force only using 1 core per R, you use RhpcBLASctl and blas_set_num_threads (below). See this here: https://github.com/lme4/lme4/issues/492
# I guess for running the initial scan it would be better to let it use more, but when you are parellelising in the perms better to force 1, and then you can control how many perms / cores run in parallel. But just stuck with 1 throughout as the initial scan doesn't take long anyway. 



#In some case better not to put an out file because it will be huge if glmm spits lots of errors
	
library(qtl)
library(doParallel)
library(lme4)
library(MASS)

library(RhpcBLASctl)
blas_set_num_threads(1) #

options(bitmapType='cairo')
print(version)

#Single locus scan using binomial glmm. Returns scanone object that can be plotted.
scanone.glm<- function(cross=r_qtl_F2s, pheno.col=c(15,16), cross_direction=c("both"),chromosomes=c(1:21)){
		pheno <- pull.pheno(cross, pheno.col)
		cross <- subset(cross, ind=!is.na(pheno)[,1])
		pheno <- pheno[which(!is.na(pheno[,1])),]
		pgm <- pull.pheno(cross, 6)
		id <- pull.pheno(cross, 2)
		chr <- names(cross$geno)
		result <- NULL
		for(i in chromosomes) {
			if(!("prob" %in% names(cross$geno[[i]]))) {
					warning("First running calc.genoprob.")
					cross <- calc.genoprob(cross)
					}
			p <- cross$geno[[i]]$prob
			# pull out map; drop last column of probabilities
			if(i!=21){
				M0<-logLik(glmer(cbind(pheno[,1],pheno[,2]) ~ 1 + (1|id),family="binomial",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000))))
			}else{ 
				if(cross_direction=="both"){
					M0_sex<-logLik(glmer(cbind(pheno[,1],pheno[,2]) ~ pgm + (1|id),family="binomial",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000))))
				}else if (cross_direction=="one_way"){
					M0_sex<-logLik(glmer(cbind(pheno[,1],pheno[,2]) ~ (1|id),family="binomial",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000))))
				}
			}	
			map <- attr(p, "map")
			p <- p[,,-dim(p)[3],drop=FALSE]
			if(i!=21){
				lod <- apply(p, 2, function(a,b)	
				(logLik(glmer(cbind(pheno[,1],pheno[,2]) ~ a + (1|id),family="binomial",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000))))-M0)/log(10))					
				z <- data.frame(chr=chr[i], pos=map, lod=lod)
			}else{
				if(cross_direction=="both"){
					lod <- apply(p, 2, function(a,b)	
					(logLik(glmer(cbind(pheno[,1],pheno[,2]) ~ a*pgm + (1|id),family="binomial",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000))))-M0_sex)/log(10))					
					z <- data.frame(chr=chr[i], pos=map, lod=lod)
				}else if (cross_direction=="one_way"){
					lod <- apply(p, 2, function(a,b)	
					(logLik(glmer(cbind(pheno[,1],pheno[,2]) ~ a + (1|id),family="binomial",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000))))-M0_sex)/log(10))					
					z <- data.frame(chr=chr[i], pos=map, lod=lod)
				}					
			}
			# special names for rows
			w <- names(map)
			o <- grep("^loc-*[0-9]+", w)
			if(length(o) > 0) # locations cited as "c*.loc*"
			w[o] <- paste("c",names(cross$geno)[i],".",w[o],sep="")
			rownames(z) <- w
			result <- rbind(result, z)
		}	
		class(result) <- c("scanone", "data.frame")
		return(result)	
}	



#permutations for single locus scan using binomial glmm
#default runs a single genome-wide test with thresholds=c(0.1,0.05,0.01)
#If cross_direction =one_way, run autosomes and sex perms separately, and specify appropriate thresholds for false positive alpha correction (below)
scanone.glm.perm <-
	function(cross=r_qtl_F2s,  pheno.col=c(12,13),n.perm=10,cross_direction=c("both"),chromosomes=c(1:21),thresholds=c(0.1,0.05,0.01)){
		pheno <- pull.pheno(cross, pheno.col)
		cross <- subset(cross, ind=!is.na(pheno)[,1])
		pheno <- pheno[which(!is.na(pheno[,1])),]
		pgm <- pull.pheno(cross, 6)
		id <- pull.pheno(cross, 2)
		chr <- names(cross$geno)
		result <- NULL
		registerDoParallel(cores=48)		
		perm.list<-foreach(perm.it=1:n.perm,.packages="lme4") %dopar%  {		
			result.temp<-NULL
			print(paste("permutation = ", perm.it,sep=""))
			rand_index<-sample(1:nrow(pheno))	
			pheno.randomised<-pheno[rand_index,]
			pgm.randomised<-pgm[rand_index]
			id.randomised<-id[rand_index]
			for(i in chromosomes) {
				if(!("prob" %in% names(cross$geno[[i]]))) {
						warning("First running calc.genoprob.")
						cross <- calc.genoprob(cross)
						}
				p <- cross$geno[[i]]$prob
				M0<-logLik(glmer(cbind(pheno.randomised[,1],pheno.randomised[,2]) ~ 1 + (1|id.randomised),family="binomial",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000))))				
				if(cross_direction=="both"){
					M0_sex<-logLik(glmer(cbind(pheno.randomised[,1],pheno.randomised[,2]) ~ pgm.randomised + (1|id.randomised),family="binomial",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000))))
				}else if (cross_direction=="one_way"){
					M0_sex<-logLik(glmer(cbind(pheno.randomised[,1],pheno.randomised[,2]) ~ (1|id.randomised),family="binomial",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000))))
				}					
				map <- attr(p, "map")
				p <- p[,,-dim(p)[3],drop=FALSE]
				lod<-c()
				if(i!=21){			
					for (j in 1:ncol(p)){
						genotypes<-p[,j,]
						lod_score<-tryCatch((logLik(glmer(cbind(pheno.randomised[,1],pheno.randomised[,2]) ~ genotypes + (1|id.randomised),family="binomial",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000))))-M0)/log(10), error=function(e) NA)
						lod<-append(lod,lod_score)
					}
				}else{
					for (j in 1:ncol(p)){
						genotypes<-p[,j,]
						if(cross_direction=="both"){						
							lod_score<-tryCatch((logLik(glmer(cbind(pheno.randomised[,1],pheno.randomised[,2]) ~ genotypes*pgm.randomised + (1|id.randomised),family="binomial",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000))))-M0_sex)/log(10), error=function(e) NA)
						}else if (cross_direction=="one_way"){	
							lod_score<-tryCatch((logLik(glmer(cbind(pheno.randomised[,1],pheno.randomised[,2]) ~ genotypes + (1|id.randomised),family="binomial",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000))))-M0_sex)/log(10), error=function(e) NA)
						}												
						lod<-append(lod,lod_score)
					}
				}	
				#z <- data.frame(chr=chr[i], pos=map, lod=lod)
				# special names for rows
				w <- names(map)
				o <- grep("^loc-*[0-9]+", w)
				if(length(o) > 0) # locations cited as "c*.loc*"
				w[o] <- paste("c",names(cross$geno)[i],".",w[o],sep="")
				#rownames(z) <- w
				result.temp <- append(result.temp, as.vector(lod))
			}
			max(result.temp)
		}
	result<-unlist(perm.list)	
	result<-sort(result)			
	result_p0.1<-result[round((1-thresholds[1])*length(result))]
	result_p0.05<-result[round((1-thresholds[2])*length(result))]
	result_p0.01<-result[round((1-thresholds[3])*length(result))]
	result_out<-data.frame(cbind(result_p0.1,result_p0.05,result_p0.01))
	return(result_out)
	stopImplicitCluster()
}	


r_qtl_F2s<-read.cross("csvr",file="qtl_dat_male_preference_trim_tips.csv",genotypes=c("EE","PE","PP"),alleles=c("E","P"))
r_qtl_F2s<-jittermap(r_qtl_F2s)
r_qtl_F2s<-calc.genoprob(r_qtl_F2s,step=1, off.end=0.0, error.prob=0.001,map.function="haldane",stepwidth="fixed")




bc_individual<-which(r_qtl_F2s$pheno[2]=="NR16-316")
r_qtl_F2s_no_BC<-r_qtl_F2s

for (pheno_index in c(1:35)){
	r_qtl_F2s_no_BC$pheno[pheno_index][bc_individual,]<-NA
}

# #set up single direction crosses
# pgm0<-which(r_qtl_F2s$pheno[6]==0)
# pgm1<-which(r_qtl_F2s$pheno[6]==1)

# r_qtl_F2s_pgm0<-r_qtl_F2s
# r_qtl_F2s_pgm1<-r_qtl_F2s

# for (pheno_index in c(7:length(names(r_qtl_F2s_pgm0$pheno)))){
	# r_qtl_F2s_pgm0$pheno[pheno_index][pgm1,]<-NA
# }

# for (pheno_index in c(7:length(names(r_qtl_F2s_pgm1$pheno)))){
	# r_qtl_F2s_pgm1$pheno[pheno_index][pgm0,]<-NA
# }



#markers for confidence intervals
list_all_markers<-read.csv("list_all_markers_Fst_trim_tips_reset_cM.csv",header=T,sep=",")
list_all_markers$X<-NULL
list_all_markers$Fst<-NULL
list_all_markers$LG[which(list_all_markers$LG=="21")]<-"X"

# list_all_markers<-read.table("C:/Users/NER334/Work/Harvard/Genomics/RAD_2018/elev_x_pard/linkage_maps/all_markers")
# colnames(list_all_markers)<-list_all_markers[1,]
# list_all_markers<-list_all_markers[-1,]
# list_all_markers$cM<-as.numeric(list_all_markers$cM)
# list_all_markers$LG[which(list_all_markers$LG=="21")]<-"X"

#Permutations
n.perm<-10

# n.permX<-round(n.perm*(sum(chrlen(r_qtl_F2s)[1:20])/chrlen(r_qtl_F2s)[21])) #For the X chromosome, n.perm × L_A/L_X permutations are performed 

# #false positive alpha correction for separate sex chromosome permutations 
# L_A<-sum(chrlen(r_qtl_F2s_pgm1)[1:20])
# L_X<-sum(chrlen(r_qtl_F2s_pgm1)[21])
# L<-L_A+L_X
# LOD_A_0.1<-1-(1-0.1)^(L_A/L)
# LOD_X_0.1<-1-(1-0.1)^(L_X/L)
# LOD_A_0.05<-1-(1-0.05)^(L_A/L)
# LOD_X_0.05<-1-(1-0.05)^(L_X/L)
# LOD_A_0.01<-1-(1-0.01)^(L_A/L)
# LOD_X_0.01<-1-(1-0.01)^(L_X/L)

#Use these when you plot a separate sex chrom sig threshold
#These x values are specific to this particular map and not generalisable
autosomes_x_coords_4_plot<-c(0,2165)
sex_x_coords_4_plot<-c(2190,2244)



###################################################################
courtship_col_indices<-c(9,10)
pheno_name<-"Approach"
print(pheno_name)
png(filename = paste("results/male_preference/glm_binomial/",pheno_name,".png",sep=""), width = 1000, height = 600)
par(mfrow=c(2,1))


r_qtl_F2s.scan<-scanone.glm(cross=r_qtl_F2s,pheno.col=courtship_col_indices,cross_direction=c("both"),chromosomes=c(1:21))
r_qtl_F2s.perm<-scanone.glm.perm(cross=r_qtl_F2s,pheno.col=courtship_col_indices,n.perm=n.perm,cross_direction=c("both"),chromosomes=c(1:21))	
plot(r_qtl_F2s.scan,main=paste("r_qtl_F2s_",pheno_name,"_0.1_sig=",round(r_qtl_F2s.perm[1],4),"_0.05_sig=",round(r_qtl_F2s.perm[2],4),"_0.01_sig=",round(r_qtl_F2s.perm[3],4),sep=""))
abline(h=r_qtl_F2s.perm[1],col='green')
abline(h=r_qtl_F2s.perm[2],col='red')

r_qtl_F2s_no_BC.scan<-scanone.glm(cross=r_qtl_F2s_no_BC,pheno.col=courtship_col_indices,cross_direction=c("both"),chromosomes=c(1:21))
r_qtl_F2s_no_BC.perm<-scanone.glm.perm(cross=r_qtl_F2s_no_BC,pheno.col=courtship_col_indices,n.perm=n.perm,cross_direction=c("both"),chromosomes=c(1:21))	
plot(r_qtl_F2s_no_BC.scan,main=paste("r_qtl_F2s_no_BC_",pheno_name,"_0.1_sig=",round(r_qtl_F2s_no_BC.perm[1],4),"_0.05_sig=",round(r_qtl_F2s_no_BC.perm[2],4),"_0.01_sig=",round(r_qtl_F2s_no_BC.perm[3],4),sep=""))
abline(h=r_qtl_F2s_no_BC.perm[1],col='green')
abline(h=r_qtl_F2s_no_BC.perm[2],col='red')	
	
	
dev.off()	
##########################################################################################
courtship_col_indices<-c(11,12)
pheno_name<-"Hover"
print(pheno_name)
png(filename = paste("results/male_preference/glm_binomial/",pheno_name,".png",sep=""), width = 1000, height = 600)
par(mfrow=c(2,1))


r_qtl_F2s.scan<-scanone.glm(cross=r_qtl_F2s,pheno.col=courtship_col_indices,cross_direction=c("both"),chromosomes=c(1:21))
r_qtl_F2s.perm<-scanone.glm.perm(cross=r_qtl_F2s,pheno.col=courtship_col_indices,n.perm=n.perm,cross_direction=c("both"),chromosomes=c(1:21))	
plot(r_qtl_F2s.scan,main=paste("r_qtl_F2s_",pheno_name,"_0.1_sig=",round(r_qtl_F2s.perm[1],4),"_0.05_sig=",round(r_qtl_F2s.perm[2],4),"_0.01_sig=",round(r_qtl_F2s.perm[3],4),sep=""))
abline(h=r_qtl_F2s.perm[1],col='green')
abline(h=r_qtl_F2s.perm[2],col='red')

r_qtl_F2s_no_BC.scan<-scanone.glm(cross=r_qtl_F2s_no_BC,pheno.col=courtship_col_indices,cross_direction=c("both"),chromosomes=c(1:21))
r_qtl_F2s_no_BC.perm<-scanone.glm.perm(cross=r_qtl_F2s_no_BC,pheno.col=courtship_col_indices,n.perm=n.perm,cross_direction=c("both"),chromosomes=c(1:21))	
plot(r_qtl_F2s_no_BC.scan,main=paste("r_qtl_F2s_no_BC_",pheno_name,"_0.1_sig=",round(r_qtl_F2s_no_BC.perm[1],4),"_0.05_sig=",round(r_qtl_F2s_no_BC.perm[2],4),"_0.01_sig=",round(r_qtl_F2s_no_BC.perm[3],4),sep=""))
abline(h=r_qtl_F2s_no_BC.perm[1],col='green')
abline(h=r_qtl_F2s_no_BC.perm[2],col='red')	

	
dev.off()	
##########################################################################################
courtship_col_indices<-c(13,14)
pheno_name<-"Alight"
print(pheno_name)
png(filename = paste("results/male_preference/glm_binomial/",pheno_name,".png",sep=""), width = 1000, height = 600)
par(mfrow=c(2,1))


r_qtl_F2s.scan<-scanone.glm(cross=r_qtl_F2s,pheno.col=courtship_col_indices,cross_direction=c("both"),chromosomes=c(1:21))
r_qtl_F2s.perm<-scanone.glm.perm(cross=r_qtl_F2s,pheno.col=courtship_col_indices,n.perm=n.perm,cross_direction=c("both"),chromosomes=c(1:21))	
plot(r_qtl_F2s.scan,main=paste("r_qtl_F2s_",pheno_name,"_0.1_sig=",round(r_qtl_F2s.perm[1],4),"_0.05_sig=",round(r_qtl_F2s.perm[2],4),"_0.01_sig=",round(r_qtl_F2s.perm[3],4),sep=""))
abline(h=r_qtl_F2s.perm[1],col='green')
abline(h=r_qtl_F2s.perm[2],col='red')
	
r_qtl_F2s_no_BC.scan<-scanone.glm(cross=r_qtl_F2s_no_BC,pheno.col=courtship_col_indices,cross_direction=c("both"),chromosomes=c(1:21))
r_qtl_F2s_no_BC.perm<-scanone.glm.perm(cross=r_qtl_F2s_no_BC,pheno.col=courtship_col_indices,n.perm=n.perm,cross_direction=c("both"),chromosomes=c(1:21))	
plot(r_qtl_F2s_no_BC.scan,main=paste("r_qtl_F2s_no_BC_",pheno_name,"_0.1_sig=",round(r_qtl_F2s_no_BC.perm[1],4),"_0.05_sig=",round(r_qtl_F2s_no_BC.perm[2],4),"_0.01_sig=",round(r_qtl_F2s_no_BC.perm[3],4),sep=""))
abline(h=r_qtl_F2s_no_BC.perm[1],col='green')
abline(h=r_qtl_F2s_no_BC.perm[2],col='red')		
	
dev.off()	
##########################################################################################
courtship_col_indices<-c(15,16)
pheno_name<-"Total_courtship"
print(pheno_name)
png(filename = paste("results/male_preference/glm_binomial/",pheno_name,".png",sep=""), width = 1000, height = 600)
par(mfrow=c(2,1))


r_qtl_F2s.scan<-scanone.glm(cross=r_qtl_F2s,pheno.col=courtship_col_indices,cross_direction=c("both"),chromosomes=c(1:21))
r_qtl_F2s.perm<-scanone.glm.perm(cross=r_qtl_F2s,pheno.col=courtship_col_indices,n.perm=n.perm,cross_direction=c("both"),chromosomes=c(1:21))	
plot(r_qtl_F2s.scan,main=paste("r_qtl_F2s_",pheno_name,"_0.1_sig=",round(r_qtl_F2s.perm[1],4),"_0.05_sig=",round(r_qtl_F2s.perm[2],4),"_0.01_sig=",round(r_qtl_F2s.perm[3],4),sep=""))
abline(h=r_qtl_F2s.perm[1],col='green')
abline(h=r_qtl_F2s.perm[2],col='red')
	
r_qtl_F2s_no_BC.scan<-scanone.glm(cross=r_qtl_F2s_no_BC,pheno.col=courtship_col_indices,cross_direction=c("both"),chromosomes=c(1:21))
r_qtl_F2s_no_BC.perm<-scanone.glm.perm(cross=r_qtl_F2s_no_BC,pheno.col=courtship_col_indices,n.perm=n.perm,cross_direction=c("both"),chromosomes=c(1:21))	
plot(r_qtl_F2s_no_BC.scan,main=paste("r_qtl_F2s_no_BC_",pheno_name,"_0.1_sig=",round(r_qtl_F2s_no_BC.perm[1],4),"_0.05_sig=",round(r_qtl_F2s_no_BC.perm[2],4),"_0.01_sig=",round(r_qtl_F2s_no_BC.perm[3],4),sep=""))
abline(h=r_qtl_F2s_no_BC.perm[1],col='green')
abline(h=r_qtl_F2s_no_BC.perm[2],col='red')		
	
dev.off()	
##########################################################################################
courtship_col_indices<-c(28,29)
pheno_name<-"approach_hover"
print(pheno_name)
png(filename = paste("results/male_preference/glm_binomial/",pheno_name,".png",sep=""), width = 1000, height = 600)
par(mfrow=c(2,1))

r_qtl_F2s.scan<-scanone.glm(cross=r_qtl_F2s,pheno.col=courtship_col_indices,cross_direction=c("both"),chromosomes=c(1:21))
r_qtl_F2s.perm<-scanone.glm.perm(cross=r_qtl_F2s,pheno.col=courtship_col_indices,n.perm=n.perm,cross_direction=c("both"),chromosomes=c(1:21))	
plot(r_qtl_F2s.scan,main=paste("r_qtl_F2s_",pheno_name,"_0.1_sig=",round(r_qtl_F2s.perm[1],4),"_0.05_sig=",round(r_qtl_F2s.perm[2],4),"_0.01_sig=",round(r_qtl_F2s.perm[3],4),sep=""))
abline(h=r_qtl_F2s.perm[1],col='green')
abline(h=r_qtl_F2s.perm[2],col='red')

r_qtl_F2s_no_BC.scan<-scanone.glm(cross=r_qtl_F2s_no_BC,pheno.col=courtship_col_indices,cross_direction=c("both"),chromosomes=c(1:21))
r_qtl_F2s_no_BC.perm<-scanone.glm.perm(cross=r_qtl_F2s_no_BC,pheno.col=courtship_col_indices,n.perm=n.perm,cross_direction=c("both"),chromosomes=c(1:21))	
plot(r_qtl_F2s_no_BC.scan,main=paste("r_qtl_F2s_no_BC_",pheno_name,"_0.1_sig=",round(r_qtl_F2s_no_BC.perm[1],4),"_0.05_sig=",round(r_qtl_F2s_no_BC.perm[2],4),"_0.01_sig=",round(r_qtl_F2s_no_BC.perm[3],4),sep=""))
abline(h=r_qtl_F2s_no_BC.perm[1],col='green')
abline(h=r_qtl_F2s_no_BC.perm[2],col='red')	
	
dev.off()	
##########################################################################################