#########SCRIPT GOALS
#investigating permutations and cutoffs:
#permuting the permutation...
#generating figures of basic results
#applying permuted results to conclusions
##########
##Significance of SNPs in the real dataset is assessed based beta values estimated for that SNP in datasets where phenotypes and genotypes are permuted
###If a SNP betavalue in the real dataset falls beyond the 95% interval (or matching interval for FP rates other than 5% ) it is significant.
####We use the permutations to verify that this method works as expected for SNPs with reasonable variance
#####We then repeat for tests about patterns of significant SNPs with respect to origin, and estimated effects on different traits.
##This script considers a number of different false positive rates for some analyses, and includes more detailed analyses than are included in the associated manuscript
#########
#main parts of snps are denoted with the following notationL
##########
#Section goal/title
##########


#a shorcut function for reading in files, used at various points.
extract_fun <- function(files){
  traits <- read.table(files,header=T,sep="\t") # read in the file
    return(traits)
}  



##########
#Processing permutations
##########

#get filenames of permutations
files <- list.files(path ="Run_GEMMA_files/Perm_runs/", pattern="perm_comb", full.names = T)

# ## measuring the false positives in the permutations is a function
# fprates <- function(filename,ncols=2:1001) { #, prblim = c(0.00007,0.99993)){
# 	randgwas <- read.csv(filename,header=T,sep="\t") #read in the file
# 	fp <- sapply(1:length(ncols), function(z)  sapply(1:nrow(randgwas), function(x) # repeat over each column (each permuation), and repeat over each row (for each snp)
# 			 findInterval(randgwas[x,ncols[z]],sort(randgwas[x,(ncols)[-z]]))/length(ncols)    #is the beta for the snp outside the quantiles defined above, across the betas estimated in all permutuations (witholding the current one, so 999 permutations)
# 			)   )
# 	return(fp)	
# }
# #run the function on the files
# allfprates <- lapply(files, function(z) fprates(filename = z,ncols=2:1001) ) #, prblim = c(0.00007,0.99993) ) ) #run for every permutation file
# 
# save(allfprates,file="~/allfprates.Rdata")
####The above was run on a computing cluster, further split into four parts by running only some traits at a time (traits are independent of each other)
######It is commented out because it takes too much time to run on a personal laptop. output files are available for those who want to start from here
#######The script continues below assuming access to our output files from our 4-batch version of the above

library(lawstat)

load("../GWAS/Data/allfprates1-3.Rdata") # because of file naming, this corresponds to traits 1, 10, 11
load("../GWAS/Data/allfprates4-6.Rdata")# because of file naming, this corresponds to traits 12, 2, 3
load("../GWAS/Data/allfprates7-9.Rdata")# because of file naming, this corresponds to traits 4-6
load("../GWAS/Data/allfprates10-12.Rdata")# because of file naming, this corresponds to traits 7-9
allfprates <- c(allfprates1.3,allfprates4.6,allfprates7.9,allfprates10.12)[c(1,5:12,2:4)] 
traitnames <- c("beta.shoot267", "beta.shoot270", "beta.shoot276", "beta.shoot279", "beta.shoot313",
                    "beta.shootall", "beta.nod267", "beta.nod270", "beta.nod276", "beta.nod279", "beta.nod313",
                    "beta.nodall") # traits here correspond, in order, to traits 1-12, so are resorted
files <- list.files(path ="Run_GEMMA_files/Perm_runs/", pattern="perm_comb", full.names = T) [c(1,5:12,2:4)] #these are put in the same order
files #check arranged as expected


##########
#what is the actually number of "false positives" that we observe given this cutoff, using rates of signficant association in each permutation
## -- e.g. verify it is the rate we expect it to be
##########

#average false positives per trait
prbs <- list(c(0.995,0.005),c(0.99,0.01),c(0.985,0.015),c(0.98,0.02),c(0.975,0.025)) # 1 to 5% cutoffs
intout <- c(0.01,0.02,0.03,0.04,0.05)
snpsoutbyint <- lapply(1:length(prbs), function(y)  sapply(1:12, function(x)
		 mean(sapply(1:1000, function(z) sum(allfprates[[x]][,z]>prbs[[y]][1]  | allfprates[[x]][,z]<prbs[[y]][2]) )
		 ) ) )
snpsoutbyquant <- lapply(1:length(prbs), function(y)  sapply(1:12, function(x)
		 mean(sapply(1:1000, function(z) sum(allfprates[[x]][,z]>prbs[[y]][1]  | allfprates[[x]][,z]<prbs[[y]][2]) )
		 ) ) )
#verify that observed false positive rates in permuted data only match those expected based on the interval used
pdf("Figures/sigsnpsFPs1000perm_1-5.pdf",h=3,w=3) 
par(mfrow=c(1,1))
par(mar=c(4,4,1,1))
plot(c(0.01,0.02,0.03,0.04,0.05)~intout,pch=NA,ylim=c(0,0.06),xlim=c(0,0.06),ylab="avg proportion false positive snps",xlab="expected proportion false positive")
for(i in 1:(length(intout))){
	points(snpsoutbyint[[i]]/363~rep(intout[i],times=12))
}
abline(a=0, b=1)
dev.off()

###Here we examine variation among permuted datasets. Using the hold-one-out method for a permuted dataset, there is therefore a distribution of false positive rates -- 1 per permuted dataset
###E.g. if the FP rate is 5%, what range of actual FP percents are observed for an individually sampled random dataset? 
####out of 1000 random datasets, some will be more, some will be less. This evaluates the risk of actually having a different FP rate than expected given a cutoff
#distribution of false positives per trait
pdf("Figures/sigsnpsFP1000pertrait1percent.pdf",width=8,height=4)
par(mfrow=c(2,6))
par(oma = c(3,1,0,0))
sapply(1:12, function(x) 
	hist(sapply(1:1000, function(z) 
		100*sum(allfprates[[x]][,z]>0.995 | allfprates[[x]][,z]<0.005) /363 ),
	ylab=ifelse(x == 1 | x==7,"# permutations",""),xlab="%FP",main=traitnames[x],ylim=c(0,400))
	) 
mtext("at 99.5-0.5% cutoff",side=1, adj = 4,line=6)
dev.off()
pdf("Figures/sigsnpsFP1000pertrait5percent.pdf",width=8,height=4)
par(mfrow=c(2,6))
par(oma = c(3,1,0,0))
sapply(1:12, function(x) 
	hist(sapply(1:1000, function(z) 
		100*sum(allfprates[[x]][,z]>0.975 | allfprates[[x]][,z]<0.025) /363 ),
	ylab=ifelse(x == 1 | x==7,"# permutations",""),xlab="%FP",main=traitnames[x],ylim=c(0,400))
	) 
mtext("at 97.5-2.5% cutoff",side=1, adj = 4,line=6)
dev.off()
pdf("Figures/sigsnpsFP1000pertrait4percent.pdf",width=8,height=4)
par(mfrow=c(2,6))
par(oma = c(3,1,0,0))
sapply(1:12, function(x) 
	hist(sapply(1:1000, function(z) 
		100*sum(allfprates[[x]][,z]>0.98 | allfprates[[x]][,z]<0.02) /363 ),
	ylab=ifelse(x == 1 | x==7,"# permutations",""),xlab="%FP",main=traitnames[x],ylim=c(0,400))
	) 
mtext("at 98.0-2.0% cutoff",side=1, adj = 4,line=6)
dev.off()
pdf("Figures/sigsnpsFP1000pertrait3percent.pdf",width=8,height=4)
par(mfrow=c(2,6))
par(oma = c(3,1,0,0))
sapply(1:12, function(x) 
	hist(sapply(1:1000, function(z) 
		100*sum(allfprates[[x]][,z]>0.985 | allfprates[[x]][,z]<0.015) /363 ),
	ylab=ifelse(x == 1 | x==7,"# permutations",""),xlab="%FP",main=traitnames[x],ylim=c(0,400))
	) 
mtext("at 98.5-1.5% cutoff",side=1, adj = 4,line=6)
dev.off()
pdf("Figures/sigsnpsFP1000pertrait2percent.pdf",width=8,height=4)
par(mfrow=c(2,6))
par(oma = c(3,1,0,0))
sapply(1:12, function(x) 
	hist(sapply(1:1000, function(z) 
		100*sum(allfprates[[x]][,z]>0.99 | allfprates[[x]][,z]<0.01) /363 ),
	ylab=ifelse(x == 1 | x==7,"# permutations",""),xlab="%FP",main=traitnames[x],ylim=c(0,400))
	) 
mtext("at 99.0-1.0% cutoff",side=1, adj = 4,line=6)
dev.off()

##########
#is there beta bias -- e.g a prevalence of positive/negative SNPs. this is permuted snps only.\
###note this would not be expected given the nature of GWAS linear models, but is a useful sanity check to examine random beta distributions
##########

pdf("Figures/positivebetas1000perm.pdf",height=8,width=8)
par(mfrow=c(3,4))
for(i in 1:12){
	permut <- read.table(files[i],header=T,sep="\t") #read in the file
	hist(rowMeans(permut[,2:1001]),breaks=14, main = traitnames[[i]] ,xlab="snp avg betascores" ) #
	abline(v=0,col="green")
}
dev.off()
#####
#there is no bias in the permuted snps.
#####


##########
#are some SNPs outside the distribution? 
##########

#by SNP for each trait, FP1%: sum the times counted as significant in the 1000 permutations, find the range, double check means
sapply(1:12, function(x) range(rowSums(allfprates[[x]]>0.995 | allfprates[[x]]<0.005) )) 
sapply(1:12, function(x) mean(rowSums(allfprates[[x]]>0.995 | allfprates[[x]]<0.005) )) 
#means are just under 10, i.e. just below expected. But some snps might be high
#again, with FP 2%
sapply(1:12, function(x) range(rowSums(allfprates[[x]]>0.99 | allfprates[[x]]<0.01) )) 
sapply(1:12, function(x) mean(rowSums(allfprates[[x]]>0.99 | allfprates[[x]]<0.01) )) 
#mean is 20,  some variation beyond
#again, FP 5%
sapply(1:12, function(x) range(rowSums(allfprates[[x]]>0.975 | allfprates[[x]]<0.025) )) 
sapply(1:12, function(x) mean(rowSums(allfprates[[x]]>0.975 | allfprates[[x]]<0.025) )) 
#mean is 50,  some variation beyond.
######Now in more detail.
####and for that, it calculates false positive counts for each SNP in each trait here
#read example permutation file for positions
f1 <- read.table(files[1],header=T,sep="\t")
#FP 1%
tagged <- matrix(NA,nrow=nrow(f1),ncol=12)
for(i in 1:12){
tagged[,i] <- sapply(1:363, function(z) sum((allfprates[[i]][z,] > 0.995 | allfprates[[i]][z,] < 0.005)))
}
#FP2%
tagged2 <- matrix(NA,nrow=nrow(f1),ncol=12)
for(i in 1:12){
tagged2[,i] <- sapply(1:363, function(z) sum((allfprates[[i]][z,] > 0.99 | allfprates[[i]][z,] < 0.01)))
}
#FP5%
tagged5 <- matrix(NA,nrow=nrow(f1),ncol=12)
for(i in 1:12){
tagged5[,i] <- sapply(1:363, function(z) sum((allfprates[[i]][z,] > 0.975 | allfprates[[i]][z,] < 0.025)))
}
# for how many traits*permutations does a SNP exceed the threshold?
pdf("Figures/persnpbias1000.pdf")
plot(rowSums(tagged)~f1$ps,ylab="Times FP",xlab="genomic position",ylim=c(min(rowSums(tagged)),1000*0.01*12*1.05)) 
abline(h=1000*0.01*12)# based on max false positive rate, expected is here
abline(h=9*12,lty=2) # based on ~mean observed FP rate 
##(see above, was about 9 times in 1000 permutations, or 9*12 across 1000 permutations in 12 traits)
dev.off()
#a few SNPs do seem to be tagged because of data structure, but generally do not exceed our expected FP rate
#Again at FP 2%, and then 5%
pdf("Figures/persnpbias1000_2percent.pdf")
plot(rowSums(tagged2)~f1$ps,ylab="Times FP",xlab="genomic position",ylim=c(min(rowSums(tagged2)),1000*0.02*12*1.05)) 
abline(h=1000*0.02*12)
dev.off()
pdf("Figures/persnpbias1000_5percent.pdf")
plot(rowSums(tagged5)~f1$ps,ylab="Times FP",xlab="genomic position",ylim=c(min(rowSums(tagged5)),1000*0.05*12*1.05)) # for how many traits*permutations does a SNP exceed the threshold? expected is 120
abline(h=1000*0.05*12)
dev.off()

##########
##Basic comparisons to real data -- calling significant SNPS, checking rates, distribution of betas
##########

files <- list.files(path ="Run_GEMMA_files/Perm_runs/", pattern="perm_comb", full.names = T)[c(1,5:12,2:4)]
 # files are not in order of plant lines and need to be rearranged
files #files defined a second time (the same way), repeated for code chunk clarity.

#real files now
files2 <- list.files(path ="Run_GEMMA_files/Real_runs/", pattern=".assoc", full.names = T)[c(1,5:12,2:4)] 
files2
traits <- lapply(files2, extract_fun) #extract_fun is defined at the beginning of this script, it is a read.csv() call
realres.ps <- cbind(traits[[1]][,c("ps","beta")], traits[[2]][,c("beta")], traits[[3]][,c("beta")], 
                traits[[4]][,c("beta")], traits[[5]][,c("beta")], traits[[6]][,c("beta")], traits[[7]][,c("beta")],
                traits[[8]][,c("beta")], traits[[9]][,c("beta")], traits[[10]][,c("beta")], traits[[11]][,c("beta")],
                traits[[12]][,c("beta")])  # extract SNP position column from one, beta columns from each
 names(realres.ps) <- c("ps","beta.shoot267", "beta.shoot270", "beta.shoot276", "beta.shoot279", "beta.shoot313",
                    "beta.shootall", "beta.nod267", "beta.nod270", "beta.nod276", "beta.nod279", "beta.nod313",
                    "beta.nodall") 
#in contrast with the gwas_res file later, Real_runs/*.assoc positions are indeed sorted the same as permuted output.                     
sum(f1$ps==realres.ps$ps)/length(realres.ps$ps) #check of identical sorting.

betacols.realres <- c( 2:13 ) #first shoot columns then nod columns
colnames(realres.ps[betacols.realres])#double check this
realsnpsBETAS <- realres.ps[,betacols.realres]
realintFull <- matrix(NA,ncol=12,nrow=363) #set empty matrix
for(i in 1:12){
	permut <- read.table(files[i],header=T,sep="\t") #read in the file
	realintFull[,i] <- sapply(1:nrow(realres.ps), function(z) 
			findInterval(realsnpsBETAS[ z,i ],sort(permut[z,-1]))/1000  
			 #Find the interval where real SNPs fall relative to the permuted SNPs, same method used as when holding out one permutation abbove
			)
}	
realsigFull <- (realintFull > 0.995 | realintFull < 0.005) # call  signicant snps for FP1%
realsigFull.2 <- (realintFull > 0.99 | realintFull < 0.01) #2%
realsigFull.3 <- (realintFull > 0.985 | realintFull < 0.015) #3%
realsigFull.4 <- (realintFull > 0.98 | realintFull < 0.02) #4%
realsigFull.5 <- (realintFull > 0.975 | realintFull < 0.025) #5%
realsigvecFull <- c(sum(realsigFull ),sum(realsigFull.2),sum(realsigFull.3),sum(realsigFull.4),sum(realsigFull.5))
plot(100*realsigvecFull/(363*12)~c(1:5),ylab="% Real SNPs significant",xlab="False Positive %") # visualize the number of significant SNPs at each threshold, at just the FP rate would be: 363*12*FP; all of these are clearly greater.
abline(a=0,b=1) # compare this plot to sigsnpsFPs1000perm_1-5.pdf, if interested
#inspect for individual trait categories, for a few FP rates
data.frame(trait =colnames(realres.ps[betacols.realres]) , snps = colSums(realintFull > 0.995 | realintFull < 0.005))
data.frame(trait =colnames(realres.ps[betacols.realres]) , snps = colSums(realintFull > 0.99 | realintFull < 0.01))
data.frame(trait =colnames(realres.ps[betacols.realres]) , snps = colSums(realintFull > 0.975 | realintFull < 0.025))

dfFIsig <- data.frame(realsigFull) #format as dataframe object for export
colnames(dfFIsig) <- paste(colnames(realres.ps[betacols.realres]),".FIsig",sep="") # add column names
FIsig <- cbind(realres.ps,dfFIsig)  # add SNP positions
write.csv(FIsig,"Output/sigsnps_findInterval.csv") #write file; previously known as annasversionsigsnps_findInterval.csv, and may be labeled as such in gitrepo

pdf("Figures/positivebetas1000real.pdf",height=8,width=8) #inspect distribution of beta scores in real snps
par(mfrow=c(3,4))
for(i in 1:12){
	hist(realsnpsBETAS[,i],breaks=14, main = traitnames[[i]] ,xlab="snp avg betascores")
	abline(v=0,col="green")
}
# as for permuted SNPs (above), there is no real pattern expected or observed, given that GWAS models ask how SNPs correlate w/ deviations from a mean phenotype
dev.off()

##########
###SNP EFFECTS AND ORIGINS
###Score and tally rates of cooperative and defective SNPs, origins of SNPs, in real data and permutations, evaluate significance
##########

#SNP categories for permutations
#tally for each trait (here, trait = plant host line) shoot-nod host environment pair, which snps are in which categories, 
##defined by 4 possible states of sign of effects on plants (shoot weight) and bacteria (nodule number)
snpOTCM5 <- array(NA, dim=c(363,1000,5)) 
# up to 5, do not want to include nodall and shootall (averages across hosts), as the home away of below doesn't make sense for averages
for(i in 1:5) {
	nodtmp <- 	read.table(files[i],header=T,sep="\t")[,-1]
	shoottmp <- read.table(files[i+6],header=T,sep="\t") [,-1] #plus 6, again to avoid including nodall and shootall
 	nodPtmp <-  allfprates[[i]] > 0.975 | allfprates[[i]] < 0.025
 	shootPtmp <- allfprates[[i+6]] > 0.975 | allfprates[[i+6]] < 0.025
 	for(k in 1:363) {  #which snps fit each category, where the effect on at least one of plants or bacteria is significant at the FP rate (5% here)
 		snpOTCM5[ k, which( nodtmp[k,] > 0 & shoottmp[k,] > 0  & (shootPtmp[k,] | nodPtmp[k,])  ), i] <- "cooperative" 
 		snpOTCM5[ k, which( nodtmp[k,] < 0 & shoottmp[k,] < 0  & (shootPtmp[k,] | nodPtmp[k,])  ), i] <- "defective"
 		snpOTCM5[ k, which( nodtmp[k,] > 0 & shoottmp[k,] < 0  & (shootPtmp[k,] | nodPtmp[k,])  ), i] <- "bactwins"
 		snpOTCM5[ k, which( nodtmp[k,] < 0 & shoottmp[k,] > 0  & (shootPtmp[k,] | nodPtmp[k,])  ), i] <- "hostwins"
 		}
} #hosts in order of 267 270 276 279 313
hostorder <- c(267, 270, 276, 279, 313)

#tally SNP origins 
realres.nbu <- read.csv("./Data/gwas_res_11Nov19.csv",header=T)
 #unsorted GWAS results; this file is a curated version of the individual gwas files read in above for each trait
 #Importantly, this one includes host origin information for each SNP
sims.ps <- read.csv("./Run_GEMMA_files/Perm_runs/perm_comb_trait9.tsv",sep="\t",header=T)$ps #sims and the above gwas_res file are sorted differently, but both have ps info
#above is read as an example file for sorting order of SNP positions
sortindex <- sapply(sims.ps,function(z) which(realres.nbu$ps == z))
realres.nb <- realres.nbu[sortindex,] # now realres object sorted the same way as sims
origincols <- grep("origin",colnames(realres.nb)) #extract origin information
#for each host-environment, does the SNP occur in a bacterial line on that host? -- ancestrally segregating SNPS are not assumed yes
snpORIGIN <- matrix(NA, nrow=363,ncol=5) # only individual pairs. otherwise home vs away doesn't mean anything
for(i in 1:5){
	snpORIGIN[,i] <- (1:363)%in%grep(hostorder[i], realres.nb[,origincols])
}
#(NB permutations and real data are the same SNPs) 

#merge origin (home vs away) and outcome categories to significant GWAS SNPS from the 1000 permutations
permCD5 <- list()
for(i in 1:5){   # only individual pairs. otherwise home vs away doesn't mean anything
	tmp <- data.frame(
		Hmhostwins = sapply(1:1000, function(z) table(snpOTCM5[snpORIGIN[,i]==T,z,i])["hostwins"] ),
		Awhostwins = sapply(1:1000, function(z) table(snpOTCM5[snpORIGIN[,i]==F,z,i])["hostwins"] ),
		Hmcoop = sapply(1:1000, function(z) table(snpOTCM5[snpORIGIN[,i]==T,z,i])["cooperative"] ),
		Awcoop = sapply(1:1000, function(z) table(snpOTCM5[snpORIGIN[,i]==F,z,i])["cooperative"] ),
		Hmdef = sapply(1:1000, function(z) table(snpOTCM5[snpORIGIN[,i]==T,z,i])["defective"] ),
		Awdef = sapply(1:1000, function(z) table(snpOTCM5[snpORIGIN[,i]==F,z,i])["defective"] ),
		Hmbactwins = sapply(1:1000, function(z) table(snpOTCM5[snpORIGIN[,i]==T,z,i])["bactwins"] ),
		Awbactwins = sapply(1:1000, function(z) table(snpOTCM5[snpORIGIN[,i]==F,z,i])["bactwins"] ) 
	)
	tmp[is.na(tmp)] <- 0  #categories with no SNPs receive an NA, replace with 0
	permCD5[[i]] <- tmp
}

#same for real results
realsnpOTCM.5 <- matrix(NA,ncol=5,nrow=363) #categorize real signifcant snps for effect type
for(i in 1:5) {
 		realsnpOTCM.5[  which( realsnpsBETAS[,i] > 0 & realsnpsBETAS[,i+6] > 0  & (realsigFull.5[,i] | realsigFull.5[,i+6])  ), i] <- "cooperative"
 		realsnpOTCM.5[  which( realsnpsBETAS[,i] < 0 & realsnpsBETAS[,i+6] < 0  & (realsigFull.5[,i] | realsigFull.5[,i+6])  ), i] <- "defective"
 		realsnpOTCM.5[  which( realsnpsBETAS[,i] > 0 & realsnpsBETAS[,i+6] < 0  & (realsigFull.5[,i] | realsigFull.5[,i+6])  ), i] <- "hostwins"
 		realsnpOTCM.5[  which( realsnpsBETAS[,i] < 0 & realsnpsBETAS[,i+6] > 0  & (realsigFull.5[,i] | realsigFull.5[,i+6])  ), i] <- "bactwins"
 		}
realsnpCD5 <- data.frame( #categorize real signifcant snps across effect types for origin; e.g. as for permutation results, merge origin and outcome categoriees
		Hmhostwins = sapply(1:5, function(i) table(realsnpOTCM.5[snpORIGIN[,i]==T,i])["hostwins"] ),
		Awhostwins = sapply(1:5, function(i)table(realsnpOTCM.5[snpORIGIN[,i]==F,i])["hostwins"] ),
		Hmcoop = sapply(1:5, function(i)table(realsnpOTCM.5[snpORIGIN[,i]==T,i])["cooperative"] ),
		Awcoop =sapply(1:5, function(i) table(realsnpOTCM.5[snpORIGIN[,i]==F,i])["cooperative"]) ,
		Hmdef = sapply(1:5, function(i)table(realsnpOTCM.5[snpORIGIN[,i]==T,i])["defective"] ),
		Awdef = sapply(1:5, function(i)table(realsnpOTCM.5[snpORIGIN[,i]==F,i])["defective"] ),
		Hmbactwins = sapply(1:5, function(i)table(realsnpOTCM.5[snpORIGIN[,i]==T,i])["bactwins"]) ,
		Awbactwins =sapply(1:5, function(i) table(realsnpOTCM.5[snpORIGIN[,i]==F,i])["bactwins"] ) 
	)
	realsnpCD5[is.na(realsnpCD5)] <- 0  #categories with no SNPs receive an NA, replace with 0
rownames(realsnpCD5) <- hostorder

#evaluate
fIntCATsnps5hi <- matrix(unlist(lapply(1:5, function(z) sapply(1:ncol(permCD5[[z]]), function(k)   findInterval(realsnpCD5[z,k],sort(permCD5[[z]][,k]),left.open=T) )/ 1000 	)),ncol=8,byrow=T) #left open puts burden on counts to be higher in case of ties, especially useful when variable is small number/small variance counts
fIntCATsnps5lo <- matrix(unlist(lapply(1:5, function(z) sapply(1:ncol(permCD5[[z]]), function(k)   findInterval(realsnpCD5[z,k],sort(permCD5[[z]][,k]),left.open=F) )/ 1000 	)),ncol=8,byrow=T) #left closed puts burden on counts to be lower in case of ties
sigCATsnps5 <- matrix("ns", ncol=ncol(fIntCATsnps5lo),nrow=nrow(fIntCATsnps5lo))
sigCATsnps5[ fIntCATsnps5hi > 0.975] <- "sigHi"
 sigCATsnps5[fIntCATsnps5lo < 0.025] <- "sigLo"
colnames(sigCATsnps5) <- colnames(permCD5[[1]])
rownames(sigCATsnps5) <- hostorder
realsnpCD5
sigCATsnps5
#some categories are significantly higher or lower than the distribution of tallies observed in the permuted datasets

##a different approach to the same question -- note the conclusions remain the same.
#similar analysis as above, but as a sum of the betas, for significant SNPs.
#real data
HomeEffPlnt.sumbeta <- sapply(c("cooperative","defective","hostwins","bactwins"), function(x) sapply(1:5, function(z) sum(realsnpsBETAS[which(realsnpOTCM.5[,z]==x & snpORIGIN[,z]==TRUE), z])   ) )
AwayEffPlnt.sumbeta <- sapply(c("cooperative","defective","hostwins","bactwins"), function(x) sapply(1:5, function(z) sum(realsnpsBETAS[which(realsnpOTCM.5[,z]==x & snpORIGIN[,z]==FALSE),z])   ) )
HomeEffBact.sumbeta <- sapply(c("cooperative","defective","hostwins","bactwins"), function(x) sapply(1:5, function(z) sum(realsnpsBETAS[which(realsnpOTCM.5[,z]==x & snpORIGIN[,z]==TRUE), z+6]) ) ) #plus 6 to go to nodule number traits & skip average traits
AwayEffBact.sumbeta <- sapply(c("cooperative","defective","hostwins","bactwins"), function(x) sapply(1:5, function(z) sum(realsnpsBETAS[which(realsnpOTCM.5[,z]==x & snpORIGIN[,z]==FALSE),z+6]) ) )#as above
rownames(HomeEffPlnt.sumbeta) <- hostorder; rownames(AwayEffPlnt.sumbeta) <- hostorder; rownames(HomeEffBact.sumbeta) <- hostorder; rownames(AwayEffBact.sumbeta) <- hostorder
#permutations
HomeEffPlnt.sumbetaSIM <- array(NA, dim=c(1000,4,5))
for(i in 1:5){
	nodtmp <- 	read.table(files[i],header=T,sep="\t")[,-1]
	shoottmp <- read.table(files[i+6],header=T,sep="\t") [,-1]
 	HomeEffPlnt.sumbetaSIM[,,i] <- sapply(c("cooperative","defective","hostwins","bactwins"), function(x) sapply(1:1000, function(z) sum(shoottmp[which(snpOTCM5[,z,i]==x & snpORIGIN[,i]==TRUE), z])   ) )
}
AwayEffPlnt.sumbetaSIM <- array(NA, dim=c(1000,4,5))
for(i in 1:5){
	nodtmp <- 	read.table(files[i],header=T,sep="\t")[,-1]
	shoottmp <- read.table(files[i+6],header=T,sep="\t") [,-1]
 	AwayEffPlnt.sumbetaSIM[,,i] <- sapply(c("cooperative","defective","hostwins","bactwins"), function(x) sapply(1:1000, function(z) sum(shoottmp[which(snpOTCM5[,z,i]==x & snpORIGIN[,i]==FALSE), z])   ) )
}
HomeEffBact.sumbetaSIM <- array(NA, dim=c(1000,4,5))
for(i in 1:5){
	nodtmp <- 	read.table(files[i],header=T,sep="\t")[,-1]
	shoottmp <- read.table(files[i+6],header=T,sep="\t") [,-1]
 	HomeEffBact.sumbetaSIM[,,i] <- sapply(c("cooperative","defective","hostwins","bactwins"), function(x) sapply(1:1000, function(z) sum(nodtmp[which(snpOTCM5[,z,i]==x & snpORIGIN[,i]==TRUE), z])   ) )
}
AwayEffBact.sumbetaSIM <- array(NA, dim=c(1000,4,5))
for(i in 1:5){
	nodtmp <- 	read.table(files[i],header=T,sep="\t")[,-1]
	shoottmp <- read.table(files[i+6],header=T,sep="\t") [,-1]
 	AwayEffBact.sumbetaSIM[,,i] <- sapply(c("cooperative","defective","hostwins","bactwins"), function(x) sapply(1:1000, function(z) sum(nodtmp[which(snpOTCM5[,z,i]==x & snpORIGIN[,i]==FALSE), z])   ) )
}
#Testing whether these are higher or lower sums of betas than expected.
#betas are summed by either the effect on plants or on bacteria, regardless of which effect is significant (or whether both effects are)
#and we are interested in whether these are higher or lower than expected, for which we have to use different sorts of intervals (right or left sides of distributions open) so as to avoid nonsensical conclusions like 0 is less than 0
betasumsig.HomeEffPlnthi <-t( sapply(1:5, function(z) sapply(1:4, function(k) findInterval(HomeEffPlnt.sumbeta[z,k],sort(HomeEffPlnt.sumbetaSIM[,k,z]),left.open=T)/1000)  ))
betasumsig.AwayEffPlnthi <-t( sapply(1:5, function(z) sapply(1:4, function(k) findInterval(AwayEffPlnt.sumbeta[z,k],sort(AwayEffPlnt.sumbetaSIM[,k,z]),left.open=T)/1000)  ))
betasumsig.HomeEffBacthi <-t( sapply(1:5, function(z) sapply(1:4, function(k) findInterval(HomeEffBact.sumbeta[z,k],sort(HomeEffBact.sumbetaSIM[,k,z]),left.open=T)/1000)  ))
betasumsig.AwayEffBacthi <-t( sapply(1:5, function(z) sapply(1:4, function(k) findInterval(AwayEffBact.sumbeta[z,k],sort(AwayEffBact.sumbetaSIM[,k,z]),left.open=T)/1000)  ))
betasumsig.HomeEffPlntlo <-t( sapply(1:5, function(z) sapply(1:4, function(k) findInterval(HomeEffPlnt.sumbeta[z,k],sort(HomeEffPlnt.sumbetaSIM[,k,z]),left.open=F)/1000)  ))
betasumsig.AwayEffPlntlo <-t( sapply(1:5, function(z) sapply(1:4, function(k) findInterval(AwayEffPlnt.sumbeta[z,k],sort(AwayEffPlnt.sumbetaSIM[,k,z]),left.open=F)/1000)  ))
betasumsig.HomeEffBactlo <-t( sapply(1:5, function(z) sapply(1:4, function(k) findInterval(HomeEffBact.sumbeta[z,k],sort(HomeEffBact.sumbetaSIM[,k,z]),left.open=F)/1000)  ))
betasumsig.AwayEffBactlo <-t( sapply(1:5, function(z) sapply(1:4, function(k) findInterval(AwayEffBact.sumbeta[z,k],sort(AwayEffBact.sumbetaSIM[,k,z]),left.open=F)/1000)  ))
betasumsig.HomeEffPlnt <- matrix("ns", ncol=ncol(betasumsig.HomeEffPlnthi),nrow=nrow(betasumsig.HomeEffPlnthi))
betasumsig.HomeEffPlnt[ betasumsig.HomeEffPlnthi > 0.975] <- "sigHi"
 betasumsig.HomeEffPlnt[betasumsig.HomeEffPlntlo < 0.025] <- "sigLo"
betasumsig.AwayEffPlnt <- matrix("ns", ncol=ncol(betasumsig.HomeEffPlnthi),nrow=nrow(betasumsig.HomeEffPlnthi))
betasumsig.AwayEffPlnt[ betasumsig.AwayEffPlnthi > 0.975] <- "sigHi"
 betasumsig.AwayEffPlnt[betasumsig.AwayEffPlntlo < 0.025] <- "sigLo"
betasumsig.HomeEffBact <- matrix("ns", ncol=ncol(betasumsig.HomeEffPlnthi),nrow=nrow(betasumsig.HomeEffPlnthi))
betasumsig.HomeEffBact[ betasumsig.HomeEffBacthi > 0.975] <- "sigHi"
 betasumsig.HomeEffBact[betasumsig.HomeEffBactlo < 0.025] <- "sigLo"
betasumsig.AwayEffBact <- matrix("ns", ncol=ncol(betasumsig.HomeEffPlnthi),nrow=nrow(betasumsig.HomeEffPlnthi))
betasumsig.AwayEffBact[ betasumsig.AwayEffBacthi > 0.975] <- "sigHi"
 betasumsig.AwayEffBact[betasumsig.AwayEffBactlo < 0.025] <- "sigLo"
colnames(betasumsig.HomeEffPlnt)<- c("cooperative","defective","hostwins","bactwins");colnames(betasumsig.AwayEffPlnt)<- c("cooperative","defective","hostwins","bactwins");colnames(betasumsig.HomeEffBact)<- c("cooperative","defective","hostwins","bactwins");colnames(betasumsig.AwayEffBact)<- c("cooperative","defective","hostwins","bactwins")
rownames(betasumsig.HomeEffPlnt)<- hostorder;rownames(betasumsig.AwayEffPlnt)<- hostorder;rownames(betasumsig.HomeEffBact)<- hostorder;rownames(betasumsig.AwayEffBact)<- hostorder
HomeEffPlnt.sumbeta
 betasumsig.HomeEffPlnt
AwayEffPlnt.sumbeta
 betasumsig.AwayEffPlnt
HomeEffBact.sumbeta
 betasumsig.HomeEffBact
AwayEffBact.sumbeta
 betasumsig.AwayEffBact
#CONCLUSIONS are the same as for SNP counting. 

#a slightly different take
####### TRADEOFFS #####
##Are snps often coop in one host and def in another?
#in the real data
SNPtype.H <- sapply(c("cooperative","defective","hostwins","bactwins"), function(x) sapply(1:363, function(z) length(which(realsnpOTCM.5[z,]==x & snpORIGIN[z,]==TRUE)))  )
 #SNPS CAN OVERLAP bt home and away, i.e. same snp can be home and significant in one of these categories and also significant in another or the same category when away
SNPtype.A <- sapply(c("cooperative","defective","hostwins","bactwins"), function(x) sapply(1:363, function(z) length(which(realsnpOTCM.5[z,]==x & snpORIGIN[z,]==FALSE))) )
 #SNPS CAN OVERLAP bt home and away. 
realCH.XA <- sum(SNPtype.H[,1] >0 & rowSums(SNPtype.A[,2:4]) > 0)
 #only one SNP is cooperative in effect at home, and also ever in a different category (defective, hostwins or bactwins) away.
#***THIS SUGGESTS LOCAL ADAPTATION IS NOT ABOUT TRADING OFF.
realCH.CA <- sum(SNPtype.H[,1] >0 & (SNPtype.A[,1]) > 0) 
# but also only 6 instances (could include repeats of snps on various away hosts) are cooperative at home and cooperative on any other away host
#so not really much sharing of snps
##repeat in the PERMS
SNPtype.Hp <- lapply(1:1000, function(p) sapply(c("cooperative","defective","hostwins","bactwins"), function(x) sapply(1:363, function(z) length(which(snpOTCM5[z,p,]==x & snpORIGIN[z,]==TRUE)))  ) ) #
SNPtype.Ap <- lapply(1:1000, function(p) sapply(c("cooperative","defective","hostwins","bactwins"), function(x) sapply(1:363, function(z) length(which(snpOTCM5[z,p,]==x & snpORIGIN[z,]==FALSE))) ) ) #
pCH.CA <- unlist(lapply(1:1000, function(z) sum(SNPtype.Hp[[z]][,1] >0 &        (SNPtype.Ap[[z]][,1]  ) > 0) ))
pCH.XA <- unlist(lapply(1:1000, function(z) sum(SNPtype.Hp[[z]][,1] >0 & rowSums(SNPtype.Ap[[z]][,2:4]) > 0) ))
sigCH.CA <- c( "sighi"=findInterval(realCH.CA,sort(pCH.CA),left.open=T)/1000 > 0.975  , "lo"=findInterval(realCH.CA,sort(pCH.CA),left.open=F)/1000 < 0.025 )
sigCH.XA <- c( "sighi"=findInterval(realCH.XA,sort(pCH.CA),left.open=T)/1000 > 0.975 , "lo"=findInterval(realCH.XA,sort(pCH.CA),left.open=F)/1000  < 0.025 )
sigCH.XA #n.s.
sigCH.CA ##significantly more than in the permutations, but not necessarily interesting, as is a small number.
#as a follow-up, are betascores of SNPs similar across hosts?
cormat <- cor(realsnpsBETAS[,-c(6,12)])
cormat[upper.tri(cormat)] <- NA
rb <- colorRampPalette(colors=c(rgb(0,0,1),rgb(1,1,1),rgb(1,0,0)) ) #low is blue, high is red
par(mar=c(8,8,1,1))
 image(cormat,zlim=c(-1,1),col=rb(100)  ,xaxt="n",yaxt="n",ylab="",xlab="")
 axis(side=1, at=seq(from=0, to = 1, length.out=10), labels=colnames(realsnpsBETAS[-c(6,12)]) ,las=2)
 axis(side=2, at=seq(from=0, to = 1, length.out=10), labels=colnames(realsnpsBETAS)[-c(6,12)] ,las=2)
#shoot-nod scores within a host are positively correlated, but across hosts no 
##regardless of whether comparison is for shoot-shoot, nod-nod, or shoot-nod, beta scores are only weakly correlated when host lines differ reflecting results above 
##this also reflects the raw phenotype data for isolate quality and isolate fitness, in manuscript Figure S4.




###########################
#NOVEL SNPS 
#1. are there more novel snps with significant effects on traits than those from standing variation?  (NO, numbers are are similar given prevalence in dataset)
#2. do novel snps have stronger effects than those that were standing variation? -- i.e. relaxation of selection on mutuations with tradeoffs across hosts (maybe)
#3. Do denovo mutations have more negative effects on non local hosts than standing variation does? (suggests yes)
#####

#re-specifying trait names for clarity here. the same variable as earlier.
traitnames <- c("beta.shoot267", "beta.shoot270", "beta.shoot276", "beta.shoot279", "beta.shoot313",
                    "beta.shootall", "beta.nod267", "beta.nod270", "beta.nod276", "beta.nod279", "beta.nod313",
                    "beta.nodall") # 

#1. compare for real and permuted snps for rates of denovo (maf_anc==0) and standing variants (maf_anc>0)
perm.dnv <- list() # get rates of significant denovo snps in each trait and permutation
for(i in 1:12) {
	permsig <- allfprates[[i]]
	perm.dnv[[i]] <- sapply(1:1000, function(z) sum(realres.nb$maf_anc==0 & (permsig[,z]>0.975 | permsig[,z]<0.025)  ) ) 
}
real.dnv <- sapply(1:12, function(z) sum(realres.nb$maf_anc==0 & realsigFull.5[,z])  ) # now for snps tagged by gwas with real data
names(real.dnv) <- traitnames 
perm.stand <- list()
for(i in 1:12) {  # now for standing variants in permuted dataset gwases
	permsig <- allfprates[[i]]
	perm.stand[[i]] <- sapply(1:1000, function(z) sum(realres.nb$maf_anc>0 & (permsig[,z]>0.975 | permsig[,z]<0.025)  ) ) 
}
real.stand <- sapply(1:12, function(z) sum(realres.nb$maf_anc>0 & realsigFull.5[,z])  ) # and standing variants in gwas with real data
names(real.stand) <- traitnames 
perm.standVdenovo <- lapply(1:12, function(i) perm.dnv[[i]]/perm.stand[[i]] ) # take ratio of denovo to standing variants taggeed in permuted
real.standVdenovo <- real.dnv/real.stand #same for real
#evaulate whether real ratios are different than ratios found for permuted datasets.
fint.standVdenovo <- lapply(1:12, function(z) findInterval(real.standVdenovo[z], sort(perm.standVdenovo[[z]]),left.open=T)/1000) #these are proportions, not whole numbers, and so ties are not an issue 
names(fint.standVdenovo) <- names(real.standVdenovo)
real.stand
real.dnv
real.standVdenovo
unlist(fint.standVdenovo)
#the proportion of significant snps that are standing or denovo is not different between real data and permuted datasets

#2. Do denovo snps underlying traits have larger effect sizes?
#first some plots of real datasets to visually assess
pdf("Figures/betas_dnvVstand_bytraitD.pdf",height=4,width=6)
 #for each trait, plot a density histogram of beta values for the snps tagged that were standing variation and snps tagged that are denovo (overlay)
par(mfrow=c(3,4))
par(mar=c(2,4,2,2))
par(oma = c(0,0,0,4))
for(i in 1:12){
	hist((realsnpsBETAS[,i])[(realres.nb$maf_anc != 0)] ,freq=F,col=rgb(0,0,0,alpha=0.25),xlab="" ,main=traitnames[i],breaks=seq(from=-1.5,to=1.5,length.out=12))
	hist((realsnpsBETAS[,i])[(realres.nb$maf_anc == 0)] ,freq=F,add=T,col=rgb(0,0,1,alpha=0.25),breaks=seq(from=-1.5,to=1.5,length.out=12))
	}
 par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
 plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend(x=0.75,y=0.75,c("Standing","De novo"),fill=c(rgb(0,0,0,alpha=0.25),rgb(0,0,1,alpha=0.25)),bty="n",xpd=TRUE)
dev.off()
###visually looks like denovo are spanning a wider range of beta values
#subset betas to standing (271) and denovo (92) snps for each trait.
realbetasdnv <- lapply(c(1:5,7:11), function(i) realsnpsBETAS[,i][(realres.nb$maf_anc == 0)]) #note this is just now for non-averaged traits
realbetasstand <- lapply(c(1:5,7:11), function(i) realsnpsBETAS[,i][(realres.nb$maf_anc > 0)])
#again visually assess differences, this time pooling by trait and checking both density and frequency
pdf("Figures/betas_dnvVstandD.pdf",height=3,width=8)
par(mfrow=c(1,2))
par(mar=c(4,4,2,2))
par(oma = c(0,0,0,4))
betasSD <- hist(unlist(realbetasstand),col=rgb(0,0,0,alpha=0.5),main="standing vs denovo snp effects",xlab = "betas",freq=F, breaks = seq(from=-1.5,to=1.5,length.out=15)) 
hist(unlist(realbetasdnv),col=rgb(0,0,1,alpha=0.5),add=T,freq=F,breaks=betasSD$breaks) 
betasSD <- hist(unlist(realbetasstand),col=rgb(0,0,0,alpha=0.5),main="standing vs denovo snp effects",xlab = "betas",freq=T, breaks = seq(from=-1.5,to=1.5,length.out=15)) 
hist(unlist(realbetasdnv),col=rgb(0,0,1,alpha=0.5),add=T,freq=T,breaks=betasSD$breaks) 
 par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
 plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend(x=0.75,y=0.75,c("Standing","De novo"),fill=c(rgb(0,0,0,alpha=0.25),rgb(0,0,1,alpha=0.25)),bty="n",xpd=TRUE)
dev.off()
#
#Now aquire means and variances of these betas, here pooled across traits (shoot weight and nodule number and both across hosts)
realbetasdnvM <- mean(unlist(realbetasdnv))
realbetasdnvV <- var(unlist(realbetasdnv))
realbetasstandV <- var(unlist(realbetasstand))
realbetasstandM <- var(unlist(realbetasstand))
#
##collect the same information for permutations
relevanttraits <- c(1:5,7:11)
permbetasdnv <- list()
permbetasstand <- list()
for(i in 1:length(relevanttraits)){ 
#get betas by trait (list of 10), and then by permutation within trait (list of 1000), then betas for all snps in the category (92 denovo, 271 standing)
	betatmp <- 	read.table(files[relevanttraits[i]],header=T,sep="\t") [,-1]
	permbetasdnv[[i]] <- lapply(1:1000, function(p) betatmp[,p][(realres.nb$maf_anc == 0)])
	permbetasstand[[i]] <- lapply(1:1000, function(p) betatmp[,p][(realres.nb$maf_anc > 0)])
}
##for each permutation, pool across traits and get means and variances of slopes by category
permbetasdnvM <-   sapply(1:1000, function(p) mean(unlist(lapply(1:10, function(z) permbetasdnv[[z]][[p]]  ) ) ) ) 
permbetasdnvV <-   sapply(1:1000, function(p)  var(unlist(lapply(1:10, function(z) permbetasdnv[[z]][[p]]  ) ) ) )
permbetasstandM <- sapply(1:1000, function(p) mean(unlist(lapply(1:10, function(z) permbetasstand[[z]][[p]]) ) ) )
permbetasstandV <- sapply(1:1000, function(p)  var(unlist(lapply(1:10, function(z) permbetasstand[[z]][[p]]) ) ) )
permbetasMdiff <-   permbetasdnvM - permbetasstandM #difference in bettavalue?
permbetasVdiff <-   permbetasdnvV - permbetasstandV
#evaluate
findInterval(realbetasdnvM - realbetasstandM,sort(permbetasMdiff))/1000 #this asks if one category more positive or more negative betas
#means not significantly different, denovo variants not significantly more negative in effects than standing variants (pooled across all hosts without respect to SNP origin)
findInterval(realbetasdnvV - realbetasstandV,sort(permbetasVdiff))/1000 #this asks if variance of betas wider in one category than the other
#variance difference between real denovo and standing betas is significant
#
#Similar, but now not pooling by trait for first real, and then permuted datasets
realbetasdnvMt <- lapply(realbetasdnv,mean)
realbetasdnvVt <- lapply(realbetasdnv,var)
realbetasstandVt <- lapply(realbetasstand,var)
realbetasstandMt <- lapply(realbetasstand,mean)
permbetasdnvMt <-     sapply(1:1000, function(p) sapply(1:10, function(z) mean(permbetasdnv[[z]][[p]]  )   )  )
permbetasdnvVt <-     sapply(1:1000, function(p) sapply(1:10, function(z)  var(permbetasdnv[[z]][[p]]  )   )  )
permbetasstandMt <-   sapply(1:1000, function(p) sapply(1:10, function(z) mean(permbetasstand[[z]][[p]])   )  )
permbetasstandVt <-   sapply(1:1000, function(p) sapply(1:10, function(z) var(permbetasstand[[z]][[p]])   )  )
sapply(1:10, function(z) realbetasdnvMt[[z]] - realbetasstandMt[[z]]) # all traits have lower mean beta for denovo snps
sapply(1:10, function(z) findInterval(realbetasdnvMt[[z]] - realbetasstandMt[[z]],sort(permbetasdnvMt[z,]-permbetasstandMt[z,]))/1000 )#
#but the difference falls within the distributions observed in permuted datasets (for all traits individually)
sapply(1:10, function(z) realbetasdnvVt[[z]] - realbetasstandVt[[z]]) # the variance of betas is higher for denovo snps in each trait in real data
sapply(1:10, function(z) findInterval(realbetasdnvVt[[z]] - realbetasstandVt[[z]],sort(permbetasdnvVt[z,]-permbetasstandVt[z,]))/1000 )
## N.S. or marginal (2 traits) with respect to the distributions of differences in variances observed.
#brief conclusion for method: 
#significance difference in variance depends on pooling across traits by this test, denovo snps are not more likely to be associated with negative effects on isolate quality or fitness
#
#
#now considering effect strength
#compare sizes of absolute values of betas for significant snps between those that were standing and those that are denovo
# -- e.g. just another way of asking if denovo effects are drawn from a wider distribution  (averaged traits included here, though we are not interested in them)
perm.dnvBetas <- list()
perm.standBetas <- list()
for(i in 1:12){
	betatmp <- 	read.table(files[i],header=T,sep="\t") [,-1]
	permsig <- allfprates[[i]]
 	perm.dnvBetas[[i]] <- sapply(1:1000, function(z) mean(abs(betatmp[ realres.nb$maf_anc==0 & (permsig[,z]>0.975 | permsig[,z]<0.025), z] ) ) ) 
	perm.standBetas[[i]] <- sapply(1:1000, function(z) mean(abs(betatmp[ realres.nb$maf_anc>0 & (permsig[,z]>0.975 | permsig[,z]<0.025), z] ) ) ) 
}	 #for each trait and permutation, extract mean of the absolute values of betas
#evaluate: compare permuted and real datasets for difference between denovo and standing variant average absolute value of the betas (by trait)
##also compare means of absolute values of betas of denovos in permuted and real datasets, and same for standing variants.
perm.dnvVstandBetas <- lapply(1:12,function(z) perm.dnvBetas[[z]] - perm.standBetas[[z]])
real.dnvBetas <- sapply(1:12, function(z) mean(abs(realsnpsBETAS[ realres.nb$maf_anc==0 & realsigFull.5[,z], z ]  )) )
real.standBetas <- sapply(1:12, function(z) mean(abs(realsnpsBETAS[ realres.nb$maf_anc>0 & realsigFull.5[,z], z]  )) )
real.dnvVstandBetas <- real.dnvBetas - real.standBetas
fint.dnvBetas <- lapply(1:12, function(z) findInterval(real.dnvBetas[z] ,sort(perm.dnvBetas[[z]]) )/1000)
fint.standBetas <- lapply(1:12, function(z) findInterval(real.standBetas[z] ,sort(perm.standBetas[[z]]) )/1000)
fint.dnvVstandBetas  <- lapply(1:12, function(z) findInterval(real.dnvVstandBetas[z] ,sort(perm.dnvVstandBetas[[z]]) )/1000)
names(real.dnvBetas) <- traitnames; names(fint.dnvBetas) <- traitnames
names(real.standBetas) <- traitnames; names(fint.standBetas) <- traitnames
names(real.dnvVstandBetas) <- traitnames; names(fint.dnvVstandBetas) <- traitnames
real.dnvBetas
real.standBetas
real.dnvVstandBetas
unlist(fint.dnvVstandBetas)
# in the real datasets, significant denovo snps have stronger estimated trait effects than significant standing snps, but this is not significant
unlist(fint.dnvBetas)
unlist(fint.standBetas)
#largely n.s. or marginal.

#3. Do denovo mutations have more negative effects on non local hosts ? is this more significant that standing differences
##here a similar analysis, but now paying attention to whether SNPs occur in bacteria from the host line on which the trait is assessed
##estimated by trait, then differences pooled across traits at the final step, otherwise number of snps sampled is small by these categories
##still comparing distributions of beta values, not numbers of snps
#
#need denovo betas split by: on local denovo, and on nonlocal;  repeat widening beyond significant snps,; same for standing; then repeat over permuted.
real.dnvBetasHs <- lapply(c(1:5,7:11), function(z) (realsnpsBETAS[ realres.nb$maf_anc==0 & realsigFull.5[,z]  & snpORIGIN[,ifelse(z>5,z-6,z)]  , z ]  )) #sig only
real.dnvBetasAs <- lapply(c(1:5,7:11), function(z) (realsnpsBETAS[ realres.nb$maf_anc==0 & realsigFull.5[,z]  & !snpORIGIN[,ifelse(z>5,z-6,z)]  , z ]  )) #sig only
real.dnvBetasH <- lapply(c(1:5,7:11), function(z) (realsnpsBETAS[ realres.nb$maf_anc==0   & snpORIGIN[,ifelse(z>5,z-6,z)]  , z ]  ))  #all
real.dnvBetasA <- lapply(c(1:5,7:11), function(z) (realsnpsBETAS[ realres.nb$maf_anc==0   & !snpORIGIN[,ifelse(z>5,z-6,z)]  , z ]  )) #all
#these compare denovo snps home effects for snps from a focal hos, to denovo snps from other all other hosts on the focal host. repeated across all hosts as focal.
# STANDING VARIANTS
real.standBetasH <- lapply(c(1:5,7:11), function(z) (realsnpsBETAS[ realres.nb$maf_anc>0   & snpORIGIN[,ifelse(z>5,z-6,z)]  , z ]  )) #all
real.standBetasA <- lapply(c(1:5,7:11), function(z) (realsnpsBETAS[ realres.nb$maf_anc>0   & !snpORIGIN[,ifelse(z>5,z-6,z)]  , z ]  )) 
real.standBetasHs <- lapply(c(1:5,7:11), function(z) (realsnpsBETAS[ realres.nb$maf_anc>0 & realsigFull.5[,z]   & snpORIGIN[,ifelse(z>5,z-6,z)]  , z ]  )) #sig only
real.standBetasAs <- lapply(c(1:5,7:11), function(z) (realsnpsBETAS[ realres.nb$maf_anc>0 & realsigFull.5[,z]   & !snpORIGIN[,ifelse(z>5,z-6,z)]  , z ]  )) 
#DIFFERENCES of mean and variance of betas between hom and away cases. for denovo and standing, all snps and significant snps only.
#denovo: betas mean(home) - mean(away); var(Home) - var(Away)
meand.dnvreal <- mean(unlist(real.dnvBetasH)) -  mean(unlist(real.dnvBetasA)) #all; denovo
vard.dnvreal <- var(unlist(real.dnvBetasH)) -  var(unlist(real.dnvBetasA)) 
meand.standreal <- mean(unlist(real.standBetasH)) -  mean(unlist(real.standBetasA)) #standing
vard.standreal <- var(unlist(real.standBetasH)) -  var(unlist(real.standBetasA)) 
meand.dnvreals <- mean(unlist(real.dnvBetasHs)) -  mean(unlist(real.dnvBetasAs)) #sig; denovo
vard.dnvreals <- var(unlist(real.dnvBetasHs)) -  var(unlist(real.dnvBetasAs)) 
meand.standreals <- mean(unlist(real.standBetasHs)) -  mean(unlist(real.standBetasAs)) #standing
vard.standreals <- var(unlist(real.standBetasHs)) -  var(unlist(real.standBetasAs)) 
#SAME FOR PERMUTED datasets
perm.dnvBetasH <- list()
perm.dnvBetasA <- list() #note that these lists go to 11, with an empty entry at 6, they are unlisted across traits before calculating means and variances, however, so this doesn't matter
for(i in c(1:5,7:11)){ #denovo, all 
	betatmp <- 	read.table(files[i],header=T,sep="\t")[,-1]
 	perm.dnvBetasH[[i]] <- lapply(1:1000, function(z) betatmp[ realres.nb$maf_anc==0   & snpORIGIN[,ifelse(i>5,i-6,i)]  , z ] ) 
 	perm.dnvBetasA[[i]] <- lapply(1:1000, function(z) betatmp[ realres.nb$maf_anc==0   & !snpORIGIN[,ifelse(i>5,i-6,i)]  , z ] ) 
}
perm.dnvBetasHs <- list()
perm.dnvBetasAs <- list()
for(i in c(1:5,7:11)){ #denovo, sig only
	permsig <- allfprates[[i]]
	betatmp <- 	read.table(files[i],header=T,sep="\t")[,-1]
 	perm.dnvBetasHs[[i]] <- lapply(1:1000, function(z) betatmp[ realres.nb$maf_anc==0   & snpORIGIN[,ifelse(i>5,i-6,i)]  &  (permsig[,z]>0.975 | permsig[,z]<0.025) , z ] ) 
 	perm.dnvBetasAs[[i]] <- lapply(1:1000, function(z) betatmp[ realres.nb$maf_anc==0   & !snpORIGIN[,ifelse(i>5,i-6,i)] &  (permsig[,z]>0.975 | permsig[,z]<0.025) , z ] ) 
}
perm.standBetasH <- list()
perm.standBetasA <- list()
for(i in c(1:5,7:11)){ #standing, all
	betatmp <- 	read.table(files[i],header=T,sep="\t")[,-1]
 	perm.standBetasH[[i]] <- lapply(1:1000, function(z) betatmp[ realres.nb$maf_anc > 0   & snpORIGIN[,ifelse(i>5,i-6,i)]  , z ] ) 
 	perm.standBetasA[[i]] <- lapply(1:1000, function(z) betatmp[ realres.nb$maf_anc > 0   & !snpORIGIN[,ifelse(i>5,i-6,i)]  , z ] ) 
}
perm.standBetasHs <- list()
perm.standBetasAs <- list()
for(i in c(1:5,7:11)){#standing, sig only
	betatmp <- 	read.table(files[i],header=T,sep="\t")[,-1]
 	perm.standBetasHs[[i]] <- lapply(1:1000, function(z) betatmp[ realres.nb$maf_anc > 0   & snpORIGIN[,ifelse(i>5,i-6,i)] &  (permsig[,z]>0.975 | permsig[,z]<0.025)   , z ] ) 
 	perm.standBetasAs[[i]] <- lapply(1:1000, function(z) betatmp[ realres.nb$maf_anc > 0   & !snpORIGIN[,ifelse(i>5,i-6,i)] &  (permsig[,z]>0.975 | permsig[,z]<0.025)  , z ] ) 
}
#DIFFERENCES of mean and variance of betas between hom and away cases (permutations). for denovo and standing, all snps and significant snps only.
#denovo.
meand.dnvperm <-  sapply(1:1000, function(z)  mean(unlist(lapply(perm.dnvBetasH, function(k) k[[z]]))) - mean(unlist(lapply(perm.dnvBetasA, function(k) k[[z]]))) )
vard.dnvperm <-  sapply(1:1000, function(z)  var(unlist(lapply(perm.dnvBetasH, function(k) k[[z]]))) - var(unlist(lapply(perm.dnvBetasA, function(k) k[[z]]))) )
meand.dnvperms <-  sapply(1:1000, function(z)  mean(unlist(lapply(perm.dnvBetasHs, function(k) k[[z]]))) - mean(unlist(lapply(perm.dnvBetasAs, function(k) k[[z]]))) )
vard.dnvperms <-  sapply(1:1000, function(z)  var(unlist(lapply(perm.dnvBetasHs, function(k) k[[z]]))) - var(unlist(lapply(perm.dnvBetasAs, function(k) k[[z]]))) )
#standing 
meand.standperm <-  sapply(1:1000, function(z)  mean(unlist(lapply(perm.standBetasH, function(k) k[[z]]))) - mean(unlist(lapply(perm.standBetasA, function(k) k[[z]]))) )
vard.standperm <-  sapply(1:1000, function(z)  var(unlist(lapply(perm.standBetasH, function(k) k[[z]]))) - var(unlist(lapply(perm.standBetasA, function(k) k[[z]]))) )
meand.standperms <-  sapply(1:1000, function(z)  mean(unlist(lapply(perm.standBetasHs, function(k) k[[z]]))) - mean(unlist(lapply(perm.standBetasAs, function(k) k[[z]]))) )
vard.standperms <-  sapply(1:1000, function(z)  var(unlist(lapply(perm.standBetasHs, function(k) k[[z]]))) - var(unlist(lapply(perm.standBetasAs, function(k) k[[z]]))) )
#compare real to  than differences in random permutations ; denovo and standing, all snps and sig snps only
meand.dnvreal; findInterval(meand.dnvreal , sort(meand.dnvperm))/1000	 #DENOVO home higher than DENOVO away for all snps; and this value sig > than permuted datasets
meand.dnvreals;findInterval(meand.dnvreals , sort(meand.dnvperms))/1000	 ##for sig snps only, DENOVO home > denovo away betas, but MARGINALLY sig higher permuted datasets
vard.dnvreal; findInterval(vard.dnvreal , sort(vard.dnvperm))/1000	 #variance lower in betas of denovo snps at home than away, but this is n.s., not different than permutations
vard.dnvreals; findInterval(vard.dnvreals , sort(vard.dnvperms))/1000	 ##sig snps only, variance of denovos again is lower for home than away betas, but again n.s.; pattern not different than permuted datasets
meand.standreal; findInterval(meand.standreal , sort(meand.standperm))/1000	 #means for standing variant snps lower at home than away, but difference is not outside what is observed in permutations
meand.standreals; findInterval(meand.standreals , sort(meand.standperms))/1000	 ##sig snps only, means now higher for standing variants snps in home vs away contexts; but difference within range observed in permutations 
vard.standreal; findInterval(vard.standreal , sort(vard.standperm),left.open=T)/1000   #variance lower for standing variant betas at home vs in away hosts; but this reduction in variance is significantly less than reductions in permuted data. variances are MORE EQUAL than expected by chance
vard.standreals; findInterval(vard.standreals , sort(vard.standperms))/1000	 ##sig snps only, variance of standing var home effects again LOWER for away effects; but now this reduction is STRONGER than we'd expect by chance.
#NOTE: there are especially few standing variants that were not retained by bacteria on some hosts AND are significant, see bottom left of freq. histogram below, but this is criteria for "away" standing variants
#visualize results
pdf("Figures/betas_dnvVstandHvA_dens.pdf",height=5,width=6)
par(mfrow=c(2,2))
	hist(unlist(real.dnvBetasA),col=rgb(0,0,0,alpha=0.5),freq=F,ylim=c(0,2.5),main = "effects for all de novo variants",xlab="SNP effect") #denovo, all
	hist(unlist(real.dnvBetasH),col=rgb(1,0,0,alpha=0.5),add=T,freq=F)
	hist(unlist(real.dnvBetasAs),col=rgb(0,0,0,alpha=0.5),freq=F,ylim=c(0,2.5),main = "significant de novo variants only",xlab="SNP effect")#sig
	hist(unlist(real.dnvBetasHs),col=rgb(1,0,0,alpha=0.5),add=T,freq=F)
	hist(unlist(real.standBetasA),col=rgb(0,0,0,alpha=0.5),freq=F,ylim=c(0,2.5),main = "effects for all standing variants",xlab="SNP effect")#standing,all
	hist(unlist(real.standBetasH),col=rgb(1,0,0,alpha=0.5),add=T,freq=F)
	hist(unlist(real.standBetasAs),col=rgb(0,0,0,alpha=0.5),freq=F,ylim=c(0,2.5),main = "significant standing variants only",xlab="SNP effect")#sig
	hist(unlist(real.standBetasHs),col=rgb(1,0,0,alpha=0.5),add=T,freq=F)
dev.off()
# same as above, but now frequency histograms, this highlights smaller sample sizes in denovo snps; and especially few standing variants that are not retained in some lines and can occur in an away context 
pdf("Figures/betas_dnvVstandHvA_freq.pdf",height=5,width=6)
par(mfrow=c(2,2))
	hist(unlist(real.dnvBetasA),col=rgb(0,0,0,alpha=0.5),freq=T,ylim=c(0,500),main = "effects for all de novo variants",xlab="SNP effect") #denovo,all
	hist(unlist(real.dnvBetasH),col=rgb(1,0,0,alpha=0.5),add=T,freq=T)
	hist(unlist(real.dnvBetasAs),col=rgb(0,0,0,alpha=0.5),freq=T,ylim=c(0,100),main = "significant de novo variants only",xlab="SNP effect")#sig
	hist(unlist(real.dnvBetasHs),col=rgb(1,0,0,alpha=0.5),add=T,freq=T)
	hist(unlist(real.standBetasA),col=rgb(0,0,0,alpha=0.5),	  freq=T,ylim=c(0,500),main = "effects for all standing variants",xlab="SNP effect")#standing, all
	hist(unlist(real.standBetasH),col=rgb(1,0,0,alpha=0.5),add=T,freq=T)
	hist(unlist(real.standBetasAs),col=rgb(0,0,0,alpha=0.5),freq=T,ylim=c(0,100),main = "significant standing variants only",xlab="SNP effect")#sig
	hist(unlist(real.standBetasHs),col=rgb(1,0,0,alpha=0.5),add=T,freq=T)
dev.off()

##########
#P-values of GEMMA. how do they relate to our significance measures?
#########

# files2 <- list.files(path ="./GWAS/Run_GEMMA_files/Real_runs/", pattern=".assoc", full.names = T)[c(1,5:12,2:4)] 
# files2
# traits <- lapply(files2, extract_fun)  # these are the same as the above files2 objects, and are only reproduced here for ease of reading this chunk;
###uncomment the file read in above and run if you are only using this chunk and the object "traits" does not exist
realres.pvals <- cbind(traits[[1]][,c("ps","p_lrt")], traits[[2]][,c("p_lrt")], traits[[3]][,c("p_lrt")], 
                traits[[4]][,c("p_lrt")], traits[[5]][,c("p_lrt")], traits[[6]][,c("p_lrt")], traits[[7]][,c("p_lrt")],
                traits[[8]][,c("p_lrt")], traits[[9]][,c("p_lrt")], traits[[10]][,c("p_lrt")], traits[[11]][,c("p_lrt")],
                traits[[12]][,c("p_lrt")])
 names(realres.pvals) <- c("ps","p_lrt.shoot267", "p_lrt.shoot270", "p_lrt.shoot276", "p_lrt.shoot279", "p_lrt.shoot313",
                    "p_lrt.shootall", "p_lrt.nod267", "p_lrt.nod270", "p_lrt.nod276", "p_lrt.nod279", "p_lrt.nod313",
                    "p_lrt.nodall")
#in contrast with the gwas_res file later, Real_runs/*.assoc files are indeed sorted the same as permuted output.                     
realsnpsPV <- realres.pvals[,2:13]
realsnpsPVsig <-realsnpsPV < 0.05
par(mfrow=c(3,4))
tabledYpYb <- list()
for(i in 1:12){
tabledYpYb[[i]] <- table(paste(realsnpsPVsig[,i],realsigFull.5[,i]))
}
tabledYpYb #counting SNPS in different categories; by trait
#comparison of select rates
lapply(tabledYpYb, function(z) z/363)
lapply(tabledYpYb, function(z) c(agree=(z/363)[1] +  (z/363)[4] , disagree= (z/363)[2] + (z/363)[3])  )
#some disagreement between the methods over which SNPs are TRUE. -- sums of FALSE/TRUE or TRUE/FALSE (disagree) often exceed TRUE/TRUE (agree)
