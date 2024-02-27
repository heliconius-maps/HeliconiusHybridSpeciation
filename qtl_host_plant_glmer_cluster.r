
	
library(qtl)
library(doParallel)
library(lme4)
library(MASS)

library(RhpcBLASctl)
#blas_set_num_threads(1) #

options(bitmapType='cairo')
print(version)

#Single locus scan using binomial glmm. Returns scanone object that can be plotted.
scanone.glm<- function(cross=r_qtl_F2s, pheno.col=c(9,10)){
		pheno <- pull.pheno(cross, pheno.col)
		cross <- subset(cross, ind=!is.na(pheno)[,1])
		pheno <- pheno[which(!is.na(pheno[,1])),]
		pgm <- pull.pheno(cross, 6)
		id <- pull.pheno(cross, 2)
		chr <- names(cross$geno)
		result <- NULL
		for(i in 1:nchr(cross)) {
			if(!("prob" %in% names(cross$geno[[i]]))) {
					warning("First running calc.genoprob.")
					cross <- calc.genoprob(cross)
					}
			p <- cross$geno[[i]]$prob
			M0<-logLik(glmer(cbind(pheno[,1],pheno[,2]) ~ 1 + (1|id),family="binomial",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000))))
			M0_sex<-logLik(glmer(cbind(pheno[,1],pheno[,2]) ~ pgm + (1|id),family="binomial",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000))))
			map <- attr(p, "map")
			p <- p[,,-dim(p)[3],drop=FALSE]
			if(i!=21){
				lod <- apply(p, 2, function(a,b)	
				(logLik(glmer(cbind(pheno[,1],pheno[,2]) ~ a + (1|id),family="binomial",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000))))-M0)/log(10))					
				z <- data.frame(chr=chr[i], pos=map, lod=lod)
			}else{
				lod <- apply(p, 2, function(a,b)	
				(logLik(glmer(cbind(pheno[,1],pheno[,2]) ~ a*pgm + (1|id),family="binomial",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000))))-M0_sex)/log(10))					
				z <- data.frame(chr=chr[i], pos=map, lod=lod)
			}
			w <- names(map)
			o <- grep("^loc-*[0-9]+", w)
			if(length(o) > 0)
			w[o] <- paste("c",names(cross$geno)[i],".",w[o],sep="")
			rownames(z) <- w
			result <- rbind(result, z)
		}	
		class(result) <- c("scanone", "data.frame")
		return(result)	
}	



#permutations for single locus scan using binomial glmm
scanone.glm.perm <-
	function(cross=r_qtl_F2s, pheno.col=c(9,10),n.perm=10){
		pheno <- pull.pheno(cross, pheno.col)
		cross <- subset(cross, ind=!is.na(pheno)[,1])
		pheno <- pheno[which(!is.na(pheno[,1])),]
		pgm <- pull.pheno(cross, 6)
		id <- pull.pheno(cross, 2)
		chr <- names(cross$geno)
		result <- NULL
		registerDoParallel(cores=25)		
		perm.list<-foreach(perm.it=1:n.perm) %dopar%  {		
			result.temp<-NULL
			print(paste("permutation = ", perm.it,sep=""))
			rand_index<-sample(1:nrow(pheno))				
			pheno.randomised<-pheno[rand_index,]
			pgm.randomised<-pgm[rand_index]
			id.randomised<-id[rand_index]
			for(i in 1:21) {
				if(!("prob" %in% names(cross$geno[[i]]))) {
						warning("First running calc.genoprob.")
						cross <- calc.genoprob(cross)
						}
				p <- cross$geno[[i]]$prob
				M0<-logLik(glmer(cbind(pheno.randomised[,1],pheno.randomised[,2]) ~ 1 + (1|id.randomised),family="binomial",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000))))
				M0_sex<-logLik(glmer(cbind(pheno[,1],pheno[,2]) ~ pgm.randomised + (1|id.randomised),family="binomial",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000))))
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
						lod_score<-tryCatch((logLik(glmer(cbind(pheno.randomised[,1],pheno.randomised[,2]) ~ genotypes*pgm.randomised + (1|id.randomised),family="binomial",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000))))-M0_sex)/log(10), error=function(e) NA)
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
	result_p0.1<-result[0.9*length(result)]
	result_p0.05<-result[0.95*length(result)]
	result_p0.01<-result[0.99*length(result)]
	result_out<-data.frame(cbind(result_p0.1,result_p0.05,result_p0.01))
	return(result_out)
	stopImplicitCluster()
}	


r_qtl_F2s<-read.cross("csvr",file="qtl_dat_host_plants_trim_tips.csv",genotypes=c("EE","PE","PP"),alleles=c("E","P"))
r_qtl_F2s<-jittermap(r_qtl_F2s)
r_qtl_F2s<-calc.genoprob(r_qtl_F2s,step=1, off.end=0.0, error.prob=0.001,map.function="haldane",stepwidth="fixed")


courtship_col_indices<-as.numeric(c("10","9")) #10,9 for prop laid on nitida, 9,10 for prop laid on riparia

pheno_name<-"prop_nititda"

png(filename = paste("host_plants/glmer_binomial/",pheno_name,".png",sep=""), width = 1000, height = 300)
print("running scan")
r_qtl_F2s.scan<-scanone.glm(cross=r_qtl_F2s,pheno.col=courtship_col_indices)
print("running permutations")
r_qtl_F2s.perm<-scanone.glm.perm(cross=r_qtl_F2s,pheno.col=courtship_col_indices,n.perm=10)
plot(r_qtl_F2s.scan,main=paste("r_qtl_F2s_host_plants_0.1_sig=",round(r_qtl_F2s.perm[1],4),"_0.05_sig=",round(r_qtl_F2s.perm[2],4),"_0.01_sig=",round(r_qtl_F2s.perm[3],4),sep=""))
abline(h=r_qtl_F2s.perm[1],col="green")
abline(h=r_qtl_F2s.perm[2],col="red")
abline(h=r_qtl_F2s.perm[3],col="red")
dev.off()	
          


