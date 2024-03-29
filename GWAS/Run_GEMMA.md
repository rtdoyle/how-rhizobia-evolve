Running GEMMA
================
Rebecca Batstone
2020-05-27

``` r
knitr::opts_chunk$set(echo = TRUE)

# packages
library("tidyverse") #includes ggplot2, dplyr, readr, stringr
library("cowplot") # paneled graphs
library("reshape2") # dcast function
library("magick") ## include images
library("rlist") ## bind vectors within a list
library("wesanderson") ## color palettes
```

## Recap: Variant filtering prior to GEMMA (see variant\_filtering.md for full details)

  - Filtered out most singletons (kept 43 representatives)
  - Linkage groups: picked one representative SNP for each group (kept
    31 representatives)

## Get files ready for for GEMMA:

### Convert vcf to plink format (GEMMA input) using plink

``` bash
# using plink, in GEMMA directory
$PLINK --vcf GEMMA_363.ann.vcf --recode --make-bed --out "KH35c_final_plink" --allow-extra-chr --allow-no-sex --maf 0.001
## no need to include pheno file at this point, it only takes the first column
```

## Format the phenotype file (to make .fam file)

First, we formatted the phenotypes file:

  - NO HEADERS
  - column 1 = strain identifiers
  - column 2 = phenotype associated with each strain

IMPORTANT: Ensure that strain identifiers generated from the vcf to
plink conversion ended up being in the same order as the strain
identifiers in the .fam
file.

``` r
load("../Single-inoculation experiment/Output/stand_means_for_GEMMA.Rdata") ## load stand_means.w
  
# add in vcf sample names
isolate_info <- read_csv("../Single-inoculation experiment/Data/isolate_info.csv")
```

    ## Parsed with column specification:
    ## cols(
    ##   isolate = col_character(),
    ##   batch = col_character(),
    ##   state = col_character(),
    ##   line_origin = col_character(),
    ##   plant_origin = col_character(),
    ##   strain = col_character(),
    ##   isolate_no = col_double(),
    ##   extraction_ID = col_double(),
    ##   notes = col_character(),
    ##   CAGEF_ID = col_character(),
    ##   tree_code = col_character(),
    ##   VCF_New = col_character(),
    ##   morph = col_character(),
    ##   vcf_order = col_double()
    ## )

``` r
# merge phenotypes and info by isolate
fam_file <- merge(isolate_info, stand_means.w, by = "isolate")

# order rows by vcf order
fam_file.o <- fam_file[order(fam_file$vcf_order),]

# keep only particular cols
fam_file.ok <- fam_file.o[c(12,15:26)]

# save without headers
write.table(fam_file.ok, file="./Run_GEMMA_files/pheno.txt", sep="\t",  
            col.names=FALSE, row.names=FALSE)
```

### Merge automatically generated .fam file (empty phenotype cols) with the actual phenotypes (pheno.txt) on the server

``` bash

# first, need to remove missing phenotype (field 6)
awk '{$6=""; print $0}' KH35c_final_plink.fam > KH35c_final_plink.fam2
# then, add the columns from the pheno.txt file
paste KH35c_final_plink.fam2 <(cut -f 2-13 pheno.txt) > KH35c_final_plink.fam
```

## GEMMA: Relatedness matrix

  - “-gk 1” calculates the centered relatedness matrix while
  - “-gk 2” calculates the standardized relatedness matrix;

Used the gk1 option, as recommended:

``` bash
$GEMMA -bfile KH35c_final_plink -gk 1 -miss 1 -o KH35c_final_kmat
## make sure you include all sites
```

## LMM option in GEMMA

If you include multiple phenotypes in the same .fam file, you can use
this loop -n <number> specifies column corresponding to phenotype in
.fam file. Epstein used -lmm option 4, for all three tests, but only
reports the likelihood ratio test (-lmm 2)

``` bash
# paste into bash script (KH35c_final_ulmm.sh):
for i in {1..12}
do
  $GEMMA -bfile KH35c_final_plink -k ./output/KH35c_final_kmat.cXX.txt -miss 1 -lmm 4 -n $i -o KH35c_final_ulmm_trait_$i
done
## run in the bg: nohup bash KH35c_final_ulmm.sh > KH35c_final_ulmm.out &

# extract pve's and se's from out file:
grep "pve estimate" KH35c_final_ulmm.out > pve_ulmm.txt
grep "se(pve)" KH35c_final_ulmm.out > se_ulmm.txt
paste pve_ulmm.txt se_ulmm.txt > pve_se_ulmm.txt
```

## Permutation test

### Randomize phenotypes

Use same phenotypes file as before, but shuffle values within each
phenotype for each phenotype, and produce 1000 shuffled sets.

``` r
# randomize .fam file
phenotypes <- read.delim("./Run_GEMMA_files/pheno.txt", header=FALSE)

# create a function to randomize phenotype file
random_phenos <- function(df, x){
  
    # randomize each col, excluding first (sample_IDs):
    rnd <- lapply(phenotypes[,-1], sample) 
    # bind randomize cols
    rnd.b <- as.data.frame(list.cbind(rnd)) 
    # add three columns of 0's:
    phenotypes$V14 <- 0 
    phenotypes$V15 <- 0 
    phenotypes$V16 <- 0 
    # format for GEMMA (.fam file), and save:
    phenotypes_rnd <- cbind(phenotypes[,c(1,1,14,15,16)],rnd.b) 
    write.table(phenotypes_rnd,file=paste('./Run_GEMMA_files/Random_phenos/','phenotype_rnd', x, '.fam', sep=''), 
                 quote=FALSE, row.names=FALSE,col.names=FALSE,sep="\t")
}
# create 1000 randomized phenotype files:
random <- lapply(1:1000, random_phenos, df = phenotypes) 
```

### GEMMA permutations (on server)

Need names of all file inputs to match:

``` bash
# executed directly from the Random_phenos directory (make_beds.sh)
for i in {1..1000} 
do 
   cp ../KH35c_final_plink.bed phenotype_rnd$i.bed 
   cp ../KH35c_final_plink.bim phenotype_rnd$i.bim 
done

# executed from GEMMA directory: nano gemma_random.sh
for (( i = 1; i <= 1000; i++ ));      ### Outer for loop, rnd phenotypes ###
do

    for (( j = 1 ; j <= 12; j++ )) ### Inner for loop, traits ###
    do
      $GEMMA -bfile ./Random_phenos/phenotype_rnd$i -k ./output/KH35c_final_kmat.cXX.txt -miss 1 -lmm 4 -n $j -o /perm/phenotype_rnd${i}_lmm_trait${j}
    done
done

## nohup bash gemma_random.sh > gemma_random.out &
## make sure this code is the exact same used for the actual data above
## Once it's fully run (expected 12000 assoc. files), extract pve and se estimates
grep "pve estimate" gemma_random.out > pve_ulmm_random.txt
grep "se(pve)" gemma_random.out > se_ulmm_random.txt
paste pve_ulmm_random.txt se_ulmm_random.txt > pve_se_ulmm_random.txt

# Move perms into trait directories: nano trait_directories.sh
for i in {1..12}
do
        mv *trait${i}.* trait_${i}/
done
```

### Merge outputs from perm test (on server)

``` bash
# nano merge_betas.sh
i=0
for x in {1..12}
do
   cd trait_${x}/
   cut -f 3 phenotype_rnd1_lmm_trait${x}.assoc.txt > delim
for file in *.assoc.txt
do
   i=$(($i+1))
   cut -f 8 ${file} > ${file}_${i}.temp
done
paste -d\\t delim *.temp > perm_comb_trait${x}.tsv
rm *.temp
cd ../
mv ./trait_${x}/perm_comb_trait${x}.tsv ./
echo "merged perms for trait ${x}"
done
```

Result: permutation files with first column as position locations, and
cols 2-1001 as beta scores for randomized
phenotypes

## findinterval method for checking whether betas fall outside the distribution

``` r
# First, need to combine real res into one spreadsheet:
## Import GEMMA association files from 'real data':
extract_fun <- function(files){
  traits <- read.table(files,header=T,sep="\t") # read in the file
    return(traits)
}  

files <- list.files(path ="./Run_GEMMA_files/Real_runs/", pattern=".assoc", full.names = T)[c(1,5:12,2:4)] 
## file order needs to be specified, it will read in 1, 10, 11, 12... if not specified

traits <- lapply(files, extract_fun)

# extract positions and beta scores
realres <- cbind(traits[[1]][,c("ps","beta")], traits[[2]][,c("beta")], traits[[3]][,c("beta")], 
                traits[[4]][,c("beta")], traits[[5]][,c("beta")], traits[[6]][,c("beta")], traits[[7]][,c("beta")],
                traits[[8]][,c("beta")], traits[[9]][,c("beta")], traits[[10]][,c("beta")], traits[[11]][,c("beta")],
                traits[[12]][,c("beta")])

# rename cols
names(realres) <- c("ps","beta.shoot267", "beta.shoot270", "beta.shoot276", "beta.shoot279", "beta.shoot313",
                    "beta.shootall", "beta.nod267", "beta.nod270", "beta.nod276", "beta.nod279", "beta.nod313",
                    "beta.nodall")

# make sure realres object is sorted the same way as sims
sims.ps <- read.csv("./Run_GEMMA_files/Perm_runs/perm_comb_trait1.tsv",sep="\t",header=T)$ps
sortindex <- sapply(sims.ps,function(z) which(realres$ps == z))
realres.ps <- realres[sortindex,] 

# Import GEMMA association files from permutations:
files <- list.files(path ="./Run_GEMMA_files/Perm_runs/", pattern="perm_comb_", full.names = T)[c(1,5:12,2:4)] ## file order needs to be specified

betacols.realres <- c(2:13) ## first shoot columns then nod columns
colnames(realres.ps[betacols.realres]) ## double check this
realsnpsBETAS <- realres.ps[,betacols.realres]
realintFull <- matrix(NA,ncol=12,nrow=363)
for(i in 1:12){
    permut <- read.table(files[i],header=T,sep="\t") #read in the file
    realintFull[,i] <- sapply(1:nrow(realres.ps), function(z) 
            findInterval(realsnpsBETAS[ z,i ],sort(permut[z,-1]))/1000 
            )
}   

# Does real data exceed 5% cuttoff?
realsigFull <- (realintFull > 0.975 | realintFull < 0.025)
colSums(realsigFull)
data.frame(trait =colnames(realres.ps[betacols.realres]) , snps = colSums(realintFull > 0.975 | realintFull < 0.025))

dfFIsig <- data.frame(realsigFull)
colnames(dfFIsig) <- paste(colnames(realres.ps[betacols.realres]),".FIsig",sep="")
realres.psFIsig <- cbind(realres.ps,dfFIsig)
save(realres.psFIsig, file = "./Run_GEMMA_files/Real_runs/betasigs.Rdata")
```

## Merge variant info

### First, import vcf info file, and merge with GEMMA output file

``` r
# Import VCF info file, merge
load(file = "../Variant discovery/Output/GEMMA_363_sum.Rdata") ## load GEMMA.vars.f

# load beta sig file
load(file = "./Run_GEMMA_files/Real_runs/betasigs.Rdata") ## load realres.psFIsig

# check that vars match
setdiff(GEMMA_vars.f$ps, realres.psFIsig$ps) ## no diffs
```

    ## integer(0)

``` r
setdiff(realres.psFIsig$ps,GEMMA_vars.f$ps) ## no diffs
```

    ## integer(0)

``` r
# merge with real data to get sig vars
GEMMA_vars_sig <- merge(GEMMA_vars.f, realres.psFIsig, by = "ps") # all 363 overlapping vars

# save file
save(GEMMA_vars_sig, file = "./Output/GEMMA_363_sig.Rdata")
```

## Figures

### Figure 3

Plot shoot and nod beta scores one each host and averaged over all hosts

``` r
# general settings
shapes=c(22,21,24)
palette <- "Royal1"

# load df
load(file = "./Output/GEMMA_363_sig.Rdata") # loads GEMMA_vars_sig

# Make "All" figure
sig_all_env <- subset(GEMMA_vars_sig, 
              beta.shootall.FIsig == "TRUE" | beta.nodall.FIsig == "TRUE") #Significant SNPs, either trait
sig_all_env$SIG <- ifelse(sig_all_env$beta.shootall.FIsig == "TRUE" & 
                        sig_all_env$beta.nodall.FIsig == "TRUE", "Both traits",
                  ifelse(sig_all_env$beta.shootall.FIsig == "TRUE", "Shoot biomass", 
                         "Nodule number"))
sig_all_env$ORIG <- "NA"
sig_all_env$effect <- ifelse(sig_all_env$beta.nodall > 0 & sig_all_env$beta.shootall > 0, "Cooperative", ifelse(
  sig_all_env$beta.nodall < 0 & sig_all_env$beta.shootall > 0, "Altruistic", ifelse(
    sig_all_env$beta.nodall > 0 & sig_all_env$beta.shootall < 0, "Cheater", "Defective"
)))

save(sig_all_env, file = "./Output/sig_all_env.Rdata")

# plot
gwas_res_all <- GEMMA_vars_sig
gwas_res_all$line <- "All"

p.all <- ggplot(data = gwas_res_all, aes(x=beta.nodall, y=beta.shootall)) +
  geom_rect(aes(xmin=0, xmax=Inf, ymin=-Inf, ymax=0), fill="grey90", alpha=0.5) + 
  geom_point(size = 2, alpha = 0.2) + 
  geom_smooth(method=lm, se=FALSE, colour="black") + 
  geom_point(aes(x=beta.nodall, y=beta.shootall, shape=Novel), 
             data=sig_all_env, size = 4, fill="black") + 
  geom_hline(aes(yintercept=0), linetype="dashed")+
  geom_vline(aes(xintercept=0), linetype="dashed")+
  scale_shape_manual(values = shapes, name="Variant type", labels=c("Standing", "De novo")) +
  facet_grid(. ~ line) +
  theme_bw() +
  theme(axis.text.y = element_text(size=12), 
        legend.position=c(0.15, 0.8), 
        legend.background = element_blank(),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(size=12), 
        plot.title = element_text(size=16, face = "bold"),
        panel.background = element_rect(fill=wes_palette(palette)[[3]]),
        strip.background = element_rect(fill=wes_palette(palette)[[4]]),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
p.all
```

    ## `geom_smooth()` using formula 'y ~ x'

![](Run_GEMMA_files/figure-gfm/Fig3-1.png)<!-- -->

``` r
Fig3_fun <- function(df, line){
  
  beta.shoot <- paste0("beta.shoot", line)
  perm.shoot <- paste0("beta.shoot", line, ".FIsig")
  beta.nod <- paste0("beta.nod", line)
  perm.nod <- paste0("beta.nod", line, ".FIsig")
  
## Significant SNPs, either trait
tmp <- subset(df, get(perm.shoot) == "TRUE" | get(perm.nod) == "TRUE") 
## add in sig traits (shoot, nod or both)
tmp$SIG <- ifelse(tmp[[perm.shoot]] == "TRUE" & tmp[[perm.nod]] == "TRUE", "Both traits",
                  ifelse(tmp[[perm.shoot]] == "TRUE", "Shoot biomass", "Nodule number"))
## add in shared history
tmp$ORIG <- ifelse(grepl(paste(line), tmp$origin, fixed=TRUE), TRUE, FALSE)
tmp$ORIG <- factor(tmp$ORIG, levels = c(TRUE, FALSE))
levels(tmp$ORIG) <- c("Yes", "No")
# summarize SNP effects
tmp$quadrant <- ifelse(tmp[[beta.nod]] > 0 & tmp[[beta.shoot]] > 0, "up_R", ifelse(
  tmp[[beta.nod]] < 0 & tmp[[beta.shoot]] > 0, "up_L", ifelse(
    tmp[[beta.nod]] > 0 & tmp[[beta.shoot]] < 0, "bt_R", "bt_L"
)))
tmp$effect <- ifelse(tmp[[beta.nod]] > 0 & tmp[[beta.shoot]] > 0, "Cooperative", ifelse(
  tmp[[beta.nod]] < 0 & tmp[[beta.shoot]] > 0, "Altruistic", ifelse(
    tmp[[beta.nod]] > 0 & tmp[[beta.shoot]] < 0, "Cheater", "Defective"
)))
# separate into x or y axis for plotting
tmp$axis <- tmp$quadrant
tmp <- separate(data = tmp, col = axis, into = c("y","x"), sep = "\\_")

# summarize number of vars within each quadrant
tmp2 <- tmp %>% 
  group_by(quadrant, effect, y, x, ORIG, .drop = FALSE) %>% 
  summarize(n=n())
tmp2$quadrant <- as.factor(tmp2$quadrant)
tmp2$effect <- as.factor(tmp2$effect)
tmp2$y <- as.factor(tmp2$y)
tmp2$x <- as.factor(tmp2$x)

# aggregate numbers for each quadrant
tmp3 <- aggregate(n ~ quadrant + effect + x + y, data=tmp2, paste, collapse = ", ")
tmp4 <- aggregate(n ~ quadrant + effect + x + y, data=tmp2, sum)
## inset figure
tmp2$tempvar <- paste(line)
inset.p <- ggplot(data=tmp2)+
  geom_tile(data=tmp4, aes(x=x, y=y, fill=as.numeric(n)), color="black")+
  geom_text(data=tmp3, aes(x=x, y=y, label=n), size=6)+
  geom_text(data=tmp3, aes(x=x, y=y, label=effect), size=5, nudge_y = 0.25)+
  ylab(NULL)+
  xlab(NULL)+
  facet_grid(. ~ tempvar) +
  scale_fill_gradient(low = "white", high = wes_palette("Royal1")[[4]])+
  guides(fill=FALSE)+
  theme_bw() +
  theme(axis.line = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title.y = element_blank(),
        axis.text = element_blank(),
        strip.background = element_rect(fill=wes_palette(palette)[[4]]),
        strip.text = element_text(size = 14, face = "bold"),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
# Main figure
df$tempvar <- paste(line)

p <- ggplot(data = df, aes(x=get(beta.nod), y=get(beta.shoot))) +
  geom_rect(aes(xmin=0, xmax=Inf, ymin=-Inf, ymax=0), fill="grey90", alpha=0.5) + 
  geom_point(size = 2, alpha = 0.2) + 
  geom_smooth(method=lm, se=FALSE, colour="black") + 
  geom_point(aes(x=get(beta.nod), y=get(beta.shoot), shape=Novel, fill=ORIG), 
             data=tmp, size = 4, alpha=0.8) +
  xlab("Effect size (nodule number)") + 
  ylab("Effect size (shoot biomass)") +
  facet_grid(. ~ tempvar) +
  geom_hline(aes(yintercept=0), linetype="dashed")+
  geom_vline(aes(xintercept=0), linetype="dashed")+
  scale_shape_manual(values=shapes)+ 
  scale_fill_manual(values=c(wes_palette(palette)[[2]], 
                             wes_palette(palette)[[1]]), name="Local variant?")+
  guides(shape = FALSE)+ 
  guides(fill = guide_legend(override.aes=list(shape=21))) +
  theme_bw() +
  theme(axis.text.y = element_text(size=12), 
        legend.position="none",  # c(0.15, 0.8) for 270 panel
        legend.background = element_blank(),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(size=12), 
        plot.title = element_text(size=16, face = "bold"),
        panel.background = element_rect(fill=wes_palette(palette)[[3]]),
        strip.background = element_rect(fill=wes_palette(palette)[[4]]),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
    
  return(list(inset.p, p, tmp, tmp4))
  
}

# specify vars to loop over
line.list <- c("270", "276","279", "313", "267")

fig3_plots <- lapply(line.list, Fig3_fun, df = GEMMA_vars_sig)
save(fig3_plots, file = "./Output/fig3_plots.Rdata")

# cowplot function
fig3a_base <- plot_grid(p.all, 
                       fig3_plots[[1]][[2]] + theme(legend.position=c(0.15, 0.8)), 
                       fig3_plots[[2]][[2]], 
                       fig3_plots[[3]][[2]], 
                       fig3_plots[[4]][[2]], 
                       fig3_plots[[5]][[2]], 
                       ncol = 3, nrow = 2, align="hv", labels=NULL)
```

    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'

``` r
fig3a_xaxis <- add_sub(fig3a_base, "Beta score (nodule number)", size=20, fontface = "bold")
fig3a_yaxis <- ggdraw() + draw_label("Beta score (shoot biomass)", size=20, fontface = "bold", angle=90)
fig3a <- plot_grid(fig3a_yaxis, fig3a_xaxis, ncol=2, rel_widths=c(0.05, 1))

# add in zeros for altruistic and cheaters
header <- c("quadrant","effect","y","x","n")
row_add1 <- c("up_L", "Altruistic", "up", "L", 0)
row_add2 <- c("bt_R", "Cheater", "bt", "R", 0)

missing <- rbind(row_add1, row_add2)
colnames(missing) <- header
missing.df <- as.data.frame(missing)

# add in whether significantly above or below
row_add1 <- c("up_L", "Altruistic", "up", "L", 0)
row_add2 <- c("bt_R", "Cheater", "bt", "R", 0)
row_add3 <- c("up_R", "Cooperative", "up", "R", 0)
row_add4 <- c("bt_L", "Defective", "bt", "L", 0)

# number of vars in each category
fig3b_base <- plot_grid(
                   # 270 inset    
                   fig3_plots[[1]][[1]] +
                     geom_tile(data=missing.df, aes(x=x, y=y), fill="white", color="black") +
                     geom_text(data=missing.df, aes(x=x, y=y, label="0, 0"), size=6) +
                     geom_text(data=missing.df, aes(x=x, y=y, label=effect), size=5, 
                               nudge_y = 0.25),
                   # 276 inset  
                   fig3_plots[[2]][[1]] + 
                      geom_tile(data=missing.df, aes(x=x, y=y), fill="white", color="black") +
                      geom_text(data=missing.df, aes(x=x, y=y, label="0, 0"), size=6) +
                      geom_text(data=missing.df, aes(x=x, y=y, label=effect), size=5, 
                                nudge_y = 0.25),
                   # 279 inset  
                   fig3_plots[[3]][[1]] + 
                      geom_tile(data=missing.df %>%
                                  filter(effect == "Altruistic"), 
                                aes(x=x, y=y), fill="white", color="black") +
                      geom_text(data=missing.df %>%
                                  filter(effect == "Altruistic"), 
                                aes(x=x, y=y, label="0, 0"), size=6) +
                      geom_text(data=missing.df %>%
                                  filter(effect == "Altruistic"), 
                                aes(x=x, y=y, label=effect), size=5, nudge_y = 0.25), 
                   # 313 inset  
                   fig3_plots[[4]][[1]] + 
                      geom_tile(data=missing.df, aes(x=x, y=y), fill="white", color="black") +
                      geom_text(data=missing.df, aes(x=x, y=y, label="0, 0"), size=6) +
                      geom_text(data=missing.df, aes(x=x, y=y, label=effect), size=5, 
                                nudge_y = 0.25), 
                   # 267 inset  
                   fig3_plots[[5]][[1]],
                       nrow=1, 
                       ncol=5,
                       align = "hv")

fig3b <- ggdraw() + 
  draw_plot(fig3b_base) +
  theme(plot.margin=unit(c(0,1,0,1),"cm")) +
  # 270
   geom_text(aes(x=0.055, y = 0.7, label="-"), size=6) +        ## up_L
   geom_text(aes(x=0.15, y = 0.7, label="+"), size=6) +         ## up_R
   geom_text(aes(x=0.145, y = 0.3, label="-"), size=6) +        ## bt_R
  # 276
   geom_text(aes(x=0.35, y = 0.7, label="+"), size=6) +        ## up_R
  # 279
   geom_text(aes(x=0.455, y = 0.7, label="-"), size=6) +       ## up_L
   geom_text(aes(x=0.55, y = 0.7, label="+"), size=6) +        ## up_R
   geom_text(aes(x=0.545, y = 0.3, label="-"), size=6) +       ## bt_R
  # 313
   geom_text(aes(x=0.745, y = 0.7, label="+"), size=6) +        ## up_R
  # 267
   geom_text(aes(x=0.855, y = 0.7, label="-"), size=6) +   ## up_L
   geom_text(aes(x=0.95, y = 0.7, label="+"), size=6)    ## up_R
  
  
# put them all together
fig3 <- plot_grid(fig3a, fig3b,
                  nrow=2, 
                  ncol=1,
                  align="hv",
                  rel_heights = c(1, 0.4),
                  labels=c("A", "B"))

fig3
```

![](Run_GEMMA_files/figure-gfm/Fig3-2.png)<!-- -->

``` r
save_plot("./Figures/Figure3.pdf", fig3, dpi= 1000,
          ncol = 1, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_height=5, 
          base_width=14)
```

## Table

This code generates the information in supplemental table S3 in the
paper. It summarizes variants by combining dataframes generated to make
Figure 3

``` r
load(file = "./Output/sig_all_env.Rdata") # loads sig_all_env
load(file = "./Output/fig3_plots.Rdata") # fig3_plots

sig_all <- sig_all_env %>%
  select(-starts_with("maf_"),-starts_with("beta"))
sig_all$line <- "All" 
sig_all$ORIG <- "All" 

sig_270 <- as.data.frame(fig3_plots[[1]][[3]])
sig_270 <- sig_270 %>%
  select(-quadrant, -x, -y, -starts_with("maf_"),-starts_with("beta"))
sig_270$line <- "270" 

sig_276 <- as.data.frame(fig3_plots[[2]][[3]])
sig_276 <- sig_276 %>%
  select(-quadrant, -x, -y, -starts_with("maf_"),-starts_with("beta"))
sig_276$line <- "276" 

sig_279 <- as.data.frame(fig3_plots[[3]][[3]])
sig_279 <- sig_279 %>%
  select(-quadrant, -x, -y, -starts_with("maf_"),-starts_with("beta"))
sig_279$line <- "279" 

sig_313 <- as.data.frame(fig3_plots[[4]][[3]])
sig_313 <- sig_313 %>%
  select(-quadrant, -x, -y, -starts_with("maf_"),-starts_with("beta"))
sig_313$line <- "313" 

sig_267 <- as.data.frame(fig3_plots[[5]][[3]])
sig_267 <- sig_267 %>%
  select(-quadrant, -x, -y, -starts_with("maf_"),-starts_with("beta"))
sig_267$line <- "267" 

sig_comb <- rbind(sig_all, sig_270, sig_276, sig_279, sig_313, sig_267)
sig_comb$ORIG <- recode_factor(sig_comb$ORIG, Yes = "L", No = "N")

# rename effects
sig_comb$effect <- ifelse(sig_comb$effect == "Cooperative", "Co",
                          ifelse(sig_comb$effect == "Defective", "D",
                                 ifelse(sig_comb$effect == "Cheater","Ch", "Alt")))

# replace with all when present in all hosts
sig_comb$origin <- recode_factor(sig_comb$origin, 'anc, 270, 276, 279, 313, 267' = "all")

sig_comb_sum <- sig_comb %>%
  filter(ORIG != "All") %>%
  group_by(ORIG, effect) %>%
  summarize(count = n())

#ggplot(sig_comb_sum, aes(x = ORIG, y = count, fill = effect)) +
#  geom_bar(stat="identity", position = position_dodge())

# Summarize at the position-level
sig_ps_sum <- sig_comb %>%
  mutate(effect_cat = paste(ORIG, effect, sep = '-'),
         effect_cat2 = paste(line, ORIG, effect, sep = '-')) %>%
  spread(key = line, value = effect_cat) %>%
  rename(
    "L270" = '270',
    "L276" = '276',
    "L279" = '279',
    "L313" = '313',
    "L267" = '267'
    ) %>%
  group_by(CHROM, ps, gene, RefSeq_ID, ncbi_func, var_type, origin, Novel, maf) %>%
  summarize(traits = paste(unique(SIG), collapse = ", "),
    effect_cats = paste(unique(effect_cat2), collapse = ", "),
    no_effects = n(),
    effect_All = paste(unique(All[!is.na(All)]), collapse = ", "),
    effect_270 = paste(unique(L270[!is.na(L270)]), collapse = ", "),
    effect_276 = paste(unique(L276[!is.na(L276)]), collapse = ", "),
    effect_279 = paste(unique(L279[!is.na(L279)]), collapse = ", "),
    effect_313 = paste(unique(L313[!is.na(L313)]), collapse = ", "),
    effect_267 = paste(unique(L267[!is.na(L267)]), collapse = ", ")) %>%
  as.data.frame(.)

sig_ps_sum$var_cat <- ifelse(grepl("L-Co", sig_ps_sum$effect_cats) &
                                  grepl("N-C", sig_ps_sum$effect_cats), "BC",
                             ifelse(grepl("L-D", sig_ps_sum$effect_cats) &
                                  grepl("N-D", sig_ps_sum$effect_cats), "BD",
                                ifelse(grepl("L-Co", sig_ps_sum$effect_cats), "LC", 
                                        ifelse(grepl("N-Co", sig_ps_sum$effect_cats), "NC",
                                               ifelse(grepl("N-D", sig_ps_sum$effect_cats), "ND",
                                                      ifelse(grepl("L-D", 
                                                                   sig_ps_sum$effect_cats), "LD",
                                                             "O"))))))

# rename traits
sig_ps_sum$traits <- ifelse(sig_ps_sum$traits == "Shoot biomass", "Quality",
                                      ifelse(sig_ps_sum$traits == "Nodule number","Fitness",
                                             "Both"))

write.csv(sig_ps_sum, "./Output/sig_ps_sum.csv", row.names = FALSE)
save(sig_ps_sum, file = "./Output/sig_ps_sum.Rdata")

sig_ps_sum_counts <- sig_ps_sum %>%
  group_by(var_cat) %>%
  summarize(no_vars = n())

# add var_cat to sig_comb
sig_comb_cat <- merge(sig_comb, sig_ps_sum[,c("ps","var_cat")], by = "ps", all.x = TRUE)

# rename NA in ncbi function
sig_comb_cat$ncbi_func[is.na(sig_comb_cat$ncbi_func)] <- "Intergenic"

# Summarize at the gene-level
sig_gene_sum1 <- sig_comb_cat %>%
   mutate(gene_ID = paste(CHROM, gene, ncbi_func, sep = "-")) %>%
   group_by(CHROM, gene, ncbi_func, gene_ID) %>%
   summarize(no_vars = n_distinct(ps),
            min_ps = min(ps),
            max_ps = max(ps),
            var_types = paste(unique(var_type), collapse = ", "),
            Novel = paste(unique(Novel), collapse = ", "),
            traits = paste(unique(SIG), collapse = ", "),
            no_effects = n()) %>%
  arrange(CHROM, min_ps) %>%
  as.data.frame(.)

sig_gene_sum2 <- sig_comb_cat %>%
  mutate(
    gene_ID = paste(CHROM, gene, ncbi_func, sep = "-"),
    effect_cat = paste(ORIG, effect, sep = "-")
    ) %>%
  group_by(CHROM, gene, ncbi_func, gene_ID, line, effect_cat) %>%
  mutate(no_effects = n(),
         effect_cat_no = paste0(effect_cat, "(", no_effects, ")")) %>%
  ungroup(.) %>%
  spread(key = line, value = effect_cat_no) %>%
   rename(
     "L270" = '270',
     "L276" = '276',
     "L279" = '279',
     "L313" = '313',
     "L267" = '267'
     ) %>%
  group_by(gene_ID) %>%
  summarize(
    effect_All = paste(unique(All[!is.na(All)]), collapse = ", "),
    effect_270 = paste(unique(L270[!is.na(L270)]), collapse = ", "),
    effect_276 = paste(unique(L276[!is.na(L276)]), collapse = ", "),
    effect_279 = paste(unique(L279[!is.na(L279)]), collapse = ", "),
    effect_313 = paste(unique(L313[!is.na(L313)]), collapse = ", "),
    effect_267 = paste(unique(L267[!is.na(L267)]), collapse = ", ")
    ) %>%
    as.data.frame(.)

sig_gene_sum <- merge(x = sig_gene_sum1, y = sig_gene_sum2, by = "gene_ID", all.x=TRUE)

# format position range
sig_gene_sum$ps_range <- ifelse(sig_gene_sum$min_ps == sig_gene_sum$max_ps, 
                                          sig_gene_sum$min_ps, 
                                          paste(sig_gene_sum$min_ps,
                                                sig_gene_sum$max_ps, sep = "-"))

# rename traits
sig_gene_sum$traits <- ifelse(sig_gene_sum$traits == "Shoot biomass", "Quality",
                                      ifelse(sig_gene_sum$traits == "Nodule number",
                                             "Fitness","Both"))

# save wide_gene df
write.csv(sig_gene_sum, "./Output/sig_gene_sum.csv", row.names = FALSE)
```
