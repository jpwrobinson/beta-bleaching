
## function for comparing shifted and recovering sites
## using betapart package, compute Bray-Curtis estimates for all pairwise comparisons between sits in dat1 and dat2

## We use this function to estimate spatial beta diversity among 1) recovering and 2) shifted sites in each survey year

beta.multi.man<-function(dat1, dat2){
	
	## empty matrix for filling with beta estimates
	returns<-matrix(NA, nrow=dim(dat1)[1], ncol=1)
	rownames(returns)<- rownames(dat1)

	for(j in 1:dim(dat1)[1]) {
	
	betas<-numeric()
	
	for(i in 1:dim(dat2)[1]){
		temp<-rbind(dat1[j,], dat2[i,])
		beta<-betapart::beta.pair.abund(temp)$beta.bray
		betas<-rbind(betas, beta)
		## drop comparison between the same site
		if(i == j){betas<-betas[-i]}
	}

	## average bray-curtis across sites
	site.beta<-mean(betas)
	returns[j, 1]<-site.beta

}	
	return(returns)
}
