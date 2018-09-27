
## Adapted Stegen et al. (2013) Appendix S1 
## to use betapart::beta.pair.abund to make temporal comparison across sites

null.fe<-function(mat, mat.baseline=NULL, reshuffle=1){

	## round abundances to ensure SADs are maintained in reshuffling 
	mat <- round(mat, 0)
	## set up empty data frame for beta vals
	null.beta <- data.frame(Location = NULL, beta.bray.null = NULL, rep=NULL)

	null.alpha.comp = numeric()

	for (i in 1:reshuffle) {  
		
		null.dist = mat;
		## for each species in regional pool
		for (species in 1:ncol(null.dist)) {
			## calculate total sp. abundance
			tot.abund = sum(null.dist[,species])

			## set 0 abundances for species absent across all sites
			if(tot.abund == 0){ null.dist[,species]=0 } else {

			# reset sp. abundance to 0 for null
			null.dist[,species] = 0;
			## for each individual fish in total abundance
			for (individual in 1:tot.abund) {
				## pull out a random site (1 row)
				sampled.site = sample(c(1:nrow(mat)),1);
				## fill null distribution with 1 individual at 1 site and 1 species
				null.dist[sampled.site,species] = null.dist[sampled.site,species] + 1;
			}}}


		## set up empty vector for beta estimates
		beta.est<-numeric()

		## estimate beta abundance across all site pairs
		bt <- betapart::beta.pair.abund(rbind(mat.baseline, null.dist), index.family="bray")
		## fill in beta estimate data frame
		## 3 indexes Bray-Curtis value; 
		## nrow(mat)*j loops through matrix to pull out same site comparisons
			for(j in 1:nrow(mat)){
				beta.est[j]<-as.matrix(bt[[3]])[nrow(mat)+j, j]
			}

		beta.est<-data.frame(Location = rownames(mat.baseline), 
				beta.bray.null=beta.est, rep=i)

		null.beta<-rbind(null.beta, beta.est) 

	} ## end randomize loop

	return(null.beta)

	}
