# Human generation times across the past 250,000 years

The generation times of our recent ancestors can tell us about both the biology and social organization of prehistoric humans, placing human evolution on an absolute timescale. We implement a method for predicting historical male and female generation times based on changes in the mutation spectrum.

Our method combines data from two different types of studies:
- **Mutations from pedigree studies** <br>
We apply a Dirichlet-multinomial regression to mutation count data to capture the relationship between the underlying mutation spectrum and parental ages. 
- **Variants from population genetic studies** <br>
We use human variants from the 1000 Genomes Project with allele ages estimated from the Genealogical Estimation of Variant Age (GEVA) approach.

## Summary of analysis workflow:
1. Preprocess 1000 Genomes variant data with ages from GEVA
2. Count binned variants, including for each continental population
3. Load mutation data and build Dirichlet-multinomial model
4. Estimate best-fit parental ages for variant spectrum in each bin

## Brief descriptions for folders and files in top-level of repository:
### folders
- bootstraps/<br>
Recalculate estimates for each 100x100 double-bootstrap of model and variants
- neanderthal_masked/<br> 
Reanalysis masking genomic tracts with potential Neanderthal introgression
- resample_alleleage/<br>
Reanalysis after drawing new allele ages based on 95% CI from GEVA
- var_count/<br>
Preprocess variant data, bin variants, and count each mutation class
### files
* age_modeling.R<br>
Loads mutation data and builds the probabilistic model for estimating parental ages
* analyze_main.R<br>
Analysis script for main plots, depends on age_modeling.R and plot_helper.R
* analyze_populations.R<br>
Analysis script for separate continental human populations
* calculate_SSE.R<br>
Calculate the sum of squared error (SSE) for generation time estimates
* cross_validation.R<br>
Short script for calculating sample variance SSE
* plot_helper.R<br>
Auxillary scripts for shaping output and plotting
* recombination_analysis.R<br>
Investigation of connection between recombination rate and mutation spectrum
* sim_famvariance.R<br>
Simulate variance in parental ages and calculate SSE
