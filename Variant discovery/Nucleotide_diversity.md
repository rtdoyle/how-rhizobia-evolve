Nucleotide diversity
================
Rebecca Batstone
2020-05-27

## Load packages

``` r
library(tidyverse)
library(knitr)
```

## Calculating nucelotide diversity

Script available at Dryad, <https://doi.org/10.5061/dryad.77k64>
(<doi:10.5061/dryad.77k64>)

copied (cp) KH35c\_align.filtered.recode.vcf to python\_scripts, and
renamed to KH35c\_WSM1022.vcf nohup python diversity2.py
KH35c\_WSM1022.vcf \> nuc\_div.out & Notes: - to get samples names in
correct format, used excel, then vimmed into python script - Need to
change vcf file in the script

Need to get count of total SNPs in each component (chromosome, plasA,
plasB) Divide average no. of differences by total

Used this to extract info for each region: grep ‘^REGION’ file \>
new\_file Type Name RefSeq INSDC Chr - NZ\_CP021825.1 CP021825.1 Plsm
accessoryA NZ\_CP021826.1 CP021826.1 Plsm psymA NZ\_CP021827.1
CP021827.1 Plsm psymB NZ\_CP021828.1 CP021828.1

Notes: use actual INSDC accession numbers from NCBI (Chromosome isn’t in
the code)

## Import results

``` r
# full dataset
nuc_div_chrom <- read_csv("./Nucleotide_diversity_files/nuc_div_chrom.csv")
```

    ## Parsed with column specification:
    ## cols(
    ##   chrom = col_character(),
    ##   samp1 = col_character(),
    ##   samp2 = col_character(),
    ##   num_sites_where_both_havedata = col_double(),
    ##   var_sites = col_double(),
    ##   N_both_het = col_double(),
    ##   N_hom_het = col_double(),
    ##   N_hom_alt_hom_ref = col_double(),
    ##   N_hom_hom = col_double(),
    ##   pairwise_geno = col_double()
    ## )

``` r
nuc_div_pSymA <- read_csv("./Nucleotide_diversity_files/nuc_div_pSymA.csv")
```

    ## Parsed with column specification:
    ## cols(
    ##   chrom = col_character(),
    ##   samp1 = col_character(),
    ##   samp2 = col_character(),
    ##   num_sites_where_both_havedata = col_double(),
    ##   var_sites = col_double(),
    ##   N_both_het = col_double(),
    ##   N_hom_het = col_double(),
    ##   N_hom_alt_hom_ref = col_double(),
    ##   N_hom_hom = col_double(),
    ##   pairwise_geno = col_double()
    ## )

``` r
nuc_div_pSymB <- read_csv("./Nucleotide_diversity_files/nuc_div_pSymB.csv")
```

    ## Parsed with column specification:
    ## cols(
    ##   chrom = col_character(),
    ##   samp1 = col_character(),
    ##   samp2 = col_character(),
    ##   num_sites_where_both_havedata = col_double(),
    ##   var_sites = col_double(),
    ##   N_both_het = col_double(),
    ##   N_hom_het = col_double(),
    ##   N_hom_alt_hom_ref = col_double(),
    ##   N_hom_hom = col_double(),
    ##   pairwise_geno = col_double()
    ## )

``` r
nuc_div <- rbind(nuc_div_chrom, nuc_div_pSymA, nuc_div_pSymB)

# calculate nucleotide differences
# We averaged all pairwise nucleotide differences across strains to
# obtain π, and divided it by the number of loci (variant and nonvariant)
# called by GATK to obtain per site values.

# rename genome regions
nuc_div$region <- recode_factor(nuc_div$chrom, 'CP021825.1'="Chromosome", 'CP021827.1'="pSymA", 'CP021828.1'="pSymB")

# add in total sites (both variant and invariant)
nuc_div$tot_sites <- ifelse(nuc_div$region == "Chromosome", 3667614, ifelse(nuc_div$region == "pSymA", 1279311, 1623934))
nuc_div$tot_sites <- as.numeric(nuc_div$tot_sites)

# calculate per site nuc div
nuc_div$nuc_diff <- (nuc_div$pairwise_geno)/(nuc_div$tot_sites)
nuc_div$nuc_diff.t <- (nuc_div$nuc_diff)*10000

# load isolate info file
iso_info <- read_csv("./Nucleotide_diversity_files/isolate_info_KH35c.csv", 
    col_types = cols(line_origin = col_factor(levels = c("anc1022", 
        "270", "276", "279", "313", "267")), 
        morph = col_factor(levels = c("EPS-", 
            "EPS+")),
        state = col_factor(levels = c("anc", "der"))))

# include population defs (state) in diversity dataset
nuc_div$state1 <- iso_info$state[match(nuc_div$samp1,iso_info$vcf_names)]
nuc_div$state2 <- iso_info$state[match(nuc_div$samp2,iso_info$vcf_names)]

# concatenate pops to get comp_type
nuc_div$state_comp <- do.call(paste, c(nuc_div[c("state1","state2")], sep = "vs"))

# include host_env in nuc_div dataset
nuc_div$host1 <- iso_info$line_origin[match(nuc_div$samp1,iso_info$vcf_names)]
nuc_div$host2 <- iso_info$line_origin[match(nuc_div$samp2,iso_info$vcf_names)]

# concatenate hosts to get host_comp_type
nuc_div$host_comp <- do.call(paste, c(nuc_div[c("host1","host2")], sep = "vs"))
```

## Compare pw nucleotide diffs between ancestral and derived isolates

``` r
## across all genomic regions:
div_sum_state <- nuc_div %>%
  group_by(state_comp) %>%
  summarise(mean_div = mean(nuc_diff.t), sd_div = sd(nuc_diff.t),
            mean_pw = mean(pairwise_geno), SE_pw = (sd(pairwise_geno)/sqrt(length(pairwise_geno))))
  
div_sum_state$state_comp <- recode_factor(div_sum_state$state_comp,'ancvsanc'="Ancestral",'dervsder'="Derived")

kable(div_sum_state)
```

| state\_comp | mean\_div |   sd\_div | mean\_pw |    SE\_pw |
| :---------- | --------: | --------: | -------: | --------: |
| Ancestral   | 0.2105536 | 0.1963617 | 30.81481 | 2.3000951 |
| ancvsder    | 0.2103440 | 0.1947775 | 31.51029 | 0.7670309 |
| dervsanc    | 0.2177059 | 0.1924949 | 32.64815 | 2.2365445 |
| Derived     | 0.2162197 | 0.1996412 | 32.65812 | 0.5139515 |
