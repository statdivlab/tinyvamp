# Using dilution series to remove contamination

## Scope and purpose

In this vignette, we will walk through how to use `tinyvamp` to estimate
detection efficiencies and contamination using dilution series.

For this analysis, we are going to consider a single mock community that
was sequenced 9 times at different dilutions. The composition of the
mock community is known to be 8 taxa at 12.5% relative abundance each.
248 taxa were detected, however – the remaining 240 were contaminants.

Specifically, we will fit the following model:
``` math
\begin{align}
  \text{expected counts for taxon $j$ in sample $i$} =  e^{\gamma_i} \left( p_{i j}e^{\beta_j} + \tilde{p}_j e^{\tilde{\gamma}_1 }3^{d_i} \right)
\end{align}
```

where

- $`p_{ij}`$ is the known composition of the mock community
- $`\gamma_i`$ reflects how deeply sequenced sample $`i`$ was
- $`\beta_j`$ is the detectability of taxon $`j`$, relative to the
  reference taxon *L. fermentum*. (we need to set one element equal to
  zero for identifiability, and we arbitrarily chose *L. fermentum*).
  Because detectabilities are not identifiable for contaminant taxa, we
  fix them to be zero for these taxa.
- $`\tilde{p}_j`$ is the unknown composition of the contaminants. The
  contamination profile affects all samples, but to differing degrees.
- $`d_i = 0, 1, \ldots, 8`$ is the number of 3-fold dilutions undertaken
  by sample $`i`$.
- $`\tilde{\gamma}_1`$ reflects the contamination intensity for the
  undiluted sample. Notice how the contribution of $`\tilde{p}_j`$ is
  proportional to $`e^{\tilde{\gamma}_1 }3^{d_i}`$. That is, greater and
  greater contribution for more dilute samples.

This model reflects the assumption that the ratio of expected
contaminant reads to expected non-contaminant reads is proportional to
$`3^{d_i}`$, as we saw empirically.

## Setup

We will start by loading the relevant packages. We recommend using
`speedyseq` instead of `phyloseq`. It’s awesome, and you can get it with
`remotes::install_github("mikemc/speedyseq")`.

``` r
library(tidyverse)
library(tinyvamp)
library(phyloseq) 
# library(speedyseq)
```

Now let’s load the relevant data and inspect it. This data is from
[Karstens et
al. (2019)](https://journals.asm.org/doi/full/10.1128/mSystems.00290-19?rfr_dat=cr_pub++0pubmed&url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org),
who generated 9 samples via three-fold dilutions of a synthetic
community containing 8 distinct strains of bacteria which each account
for 12.5% of the DNA in the community. Despite only 8 strains being
present in the synthetic community, 248 total strains were identified
using DADA2 (see [our paper](https://arxiv.org/pdf/2204.12733.pdf) for
details on data processing).

``` r
data("karstens_phyloseq")
karstens_phyloseq
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 248 taxa and 9 samples ]
#> sample_data() Sample Data:       [ 9 samples by 6 sample variables ]
#> tax_table()   Taxonomy Table:    [ 248 taxa by 7 taxonomic ranks ]
```

We can see the DNA concentrations and what-fold dilution each sample is
as follows.

``` r
karstens_phyloseq %>% sample_data
#>    X.SampleID SampleDescription DNA_conc    SampleType Description is.neg
#> D1         D1            MC_1:3    109.2 MockCommunity        Mock  FALSE
#> D0         D0           MC_Neat    138.8 MockCommunity        Mock  FALSE
#> D3         D3           MC_1:27     60.6 MockCommunity        Mock  FALSE
#> D2         D2            MC_1:9    105.8 MockCommunity        Mock  FALSE
#> D5         D5          MC_1:243     22.6 MockCommunity        Mock  FALSE
#> D4         D4           MC_1:81     34.9 MockCommunity        Mock  FALSE
#> D7         D7         MC_1:2187     11.4 MockCommunity        Mock  FALSE
#> D6         D6          MC_1:729     20.0 MockCommunity        Mock  FALSE
#> D8         D8         MC_1:6561     13.8 MockCommunity        Mock  FALSE
```

We know which genera correspond to taxa in the mock, so let’s create a
vector telling us which rows of the ASV table contain data on the mock
taxa

``` r
genera_data <- karstens_phyloseq %>% tax_table %>% as.data.frame %>% as_tibble %>% select("Genus") %>% pull
Pseudomonas <- genera_data %>% str_detect("Pseudomonas") %>% which
Escherichia <- genera_data %>% str_detect("Escherichia") %>% which
Salmonella <- genera_data %>% str_detect("Salmonella") %>% which
Limosilactobacillus <- genera_data %>% str_detect("Limosilactobacillus") %>% which
Enterococcus <- genera_data %>% str_detect("Enterococcus") %>% which
Staphylococcus <- genera_data %>% str_detect("Staphylococcus") %>% which
Listeria <- genera_data %>% str_detect("Listeria") %>% which
Bacillus <- genera_data %>% str_detect("Bacillus") %>% which
mock_taxa <- c(Pseudomonas,
               Escherichia,
               Salmonella[1],
               Enterococcus,
               Staphylococcus,
               Listeria,
               Bacillus,
               Limosilactobacillus)
```

Note that there were some complexities regarding Salmonella, with two
strains being detected: one strain as a likely mutant of synthetic
community member *S. enterica*. For this reason we just consider the
more prevalent Salmonella strain (hence the `Salmonella[1]` in the
above).

Now we can construct the count data matrix $`\mathbf{W}`$. We will
reorder its columns so the taxa in the mock come last, and reorder rows
so dilutions are in increasing order – this just makes our lives easier
later!

``` r
W <- karstens_phyloseq %>% otu_table %>% as.matrix
jj <- ncol(W)
W <- cbind(W[ , !(1:jj %in% mock_taxa)],
           W[ , mock_taxa])
W <- W[order(rownames(W)), ]
```

## Model assumptions and data

In this vignette, we’re going to fit a model to all the data we have,
and inspect the estimated detection efficiencies and composition of the
contaminant profile. In contrast, if you want to evaluate the
performance of `tinyvamp`’s fit, you might want to split your sample and
do cross-validation, like we did in the paper. You can check out code
for that
[here](https://github.com/statdivlab/tinyvamp_supplementary/blob/main/karstens_cv.R),
but note that it is more complicated!

To do this analysis, we are going to tell `tinyvamp` that

- the same specimen is analyzed in all nine samples (this is
  $`\mathbf{Z}`$)
- we are not interested in estimating detection efficiencies for
  contaminants (this is $`\tilde{\mathbf{X}}`$)
- there are eight taxa that are present in 12.5% abundance each, and
  that every other taxon is not present in the true composition of the
  sample (this is $`\mathbf{p}`$)

and we are going to ask `tinyvamp` to estimate

- the contamination relative abundance profile ($`\tilde{\mathbf{p}}`$)
- the detection efficiencies for mock taxa ($`\mathbf{\beta}`$)

We’re going to assume that

- every sample has the same contamination profile (this is
  $`\tilde{\mathbf{Z}}`$)
- the ratio of expected contaminant reads to expected mock reads is
  proportional to the dilution number of the sample (you can see support
  for this in Figure 2, upper right, of [the
  paper](https://arxiv.org/pdf/2204.12733.pdf)). Mo’ dilution = less
  biomass = mo’ problems!

Ok, let’s dive in!

## Parameterizing the model

The same specimen is analyzed in all nine samples. The rows of
$`\mathbf{Z}`$ represent the samples (of which there are nine), and the
columns represent the distinct specimens (of which there is one).

``` r
Z <- matrix(1, ncol = 1, nrow = 9)
```

If you had a study where you had, let’s say, two distinct specimens of
different compositions to analyze, and 3 samples for each, your
$`\mathbf{Z}`$ might look like
$`\begin{pmatrix} 1 & 0 \\ 1 & 0 \\ 1 & 0 \\ 0 & 1 \\ 0 & 1 \\ 0 & 1 \end{pmatrix}`$.

The true composition of the specimen is that there are eight taxa that
are present in 12.5% abundance each, and we ordered them last.
Everything else has 0% abundance.

``` r
P <- matrix(c(rep(0, 240), rep(1/8, 8)), nrow = 1, ncol = jj)
```

and we know this to be true ($`\mathbf{p}`$ is fixed):

``` r
P_fixed_indices <- matrix(TRUE, nrow = 1, ncol = jj)
```

Note that `P_fixed_indices`, `B_fixed_indices`, etc., are how we tell
`tinyvamp` if we *do* versus *do not* know the value of these
parameters. So by specifying `P` and `P_fixed_indices` as above, we say
“keep $`\mathbf{p}`$ fixed at the values that we tell you; you don’t
need to estimate it yourself!”

In contrast, we do need to estimate *some* of the $`\beta`$’s, a.k.a.
the detection efficiencies. We want to estimate the detection
efficiencies of the mock taxa, but not for the contaminants.

``` r
B_fixed_indices <- matrix(TRUE, ncol = jj, nrow = 1)
B_fixed_indices[1, 241:247] <- FALSE
```

(`B_fixed_indices` = FALSE means “do estimate these beta’s!”)

Hang on, why is $`\beta_{248}`$ fixed? Well, we always need some taxon
to compare detection efficiencies to – otherwise there’s no difference
between efficiencies 1, 2, and 4 versus 8, 16 and 32! (The fancy term
for this is “identifiability constraint”.) So here we’re deciding to set
*Limosilactobacillus fermentum* as the “baseline” taxon, which
corresponds to the taxon in the last column of $`\mathbf{W}`$.

Remember that the $`\beta`$’s give the detection efficiencies on the
log-scale, so let’s set:

``` r
B <- matrix(0, ncol = jj, nrow = 1)
```

Note that this analysis treats the contaminant taxa as having the same
detection efficiency as *Limosilactobacillus fermentum*.

**If a parameter is** *not fixed* (e.g., $`\beta_{241}`$ in the above),
**the corresponding value is where is** *initialized* **for the
estimation algorithm to proceed.**

The rows of `B` correspond to sequencing protocols, and here we have
only one protocol – no data on batches, for example – hence one row for
`B`. You can check out the “Compare experiments” vignette for an example
showing how to compare three different experimental protocols!

$`\mathbf{X}`$ connects the samples ($`n`$ rows) to the sequencing
protocols ($`p`$ columns). Here we have only one protocol (no data on
batches, for example), so one column. We’ll just fill this matrix with
1’s to say “estimate efficiencies for our single protocol”.

``` r
X <- matrix(1, nrow = nrow(W), ncol = 1)
```

In contrast, if you had two sequencing protocols on the same mock
community, you might set
$`\begin{pmatrix} 1 & 0 \\ \vdots & \vdots \\ 1 & 0 \\ 0 & 1 \\ \vdots & \vdots \\ 0 & 1 \end{pmatrix}`$.

Ok, phew, that’s most of the hard stuff!

We *do* need to estimate the sampling intensity parameters (this will
usually be the case). Let’s initialize them at the log-total-read count
(plus a bit of noise).

``` r
gammas_fixed_indices <- rep(FALSE, nrow(W))
gammas <- log(apply(W, 1, sum)) + rnorm(nrow(W))
```

Now let’s talk contamination! We assume a single contamination profile,
so $`\tilde{K}=1`$ and the number of columns of $`\tilde{\mathbf{Z}}`$
is one. We want to assume that the ratio of expected contaminant reads
to expected mock reads is proportional to $`3^{d}`$, where $`d`$ is the
dilution number of the sample ($`d_i=0`$ corresponds to the undiluted
sample). To do this, we set $`\tilde{\mathbf{Z}}_i`$ proportional to
$`\exp(\gamma_i) \times 3^{d_i}`$ for sample $`i`$. To address the
$`3^{d_i}`$ piece, we set

``` r
Z_tilde <- matrix(3^(0:8)/exp(mean(log(3^(0:8)))), ncol = 1)
```

and to address the $`\exp(\gamma_i)`$ piece, we set

``` r
Z_tilde_gamma_cols <- 1 
```

We do want to estimate the contamination intensity, so let’s initialize
it at zero (on the log-scale) and tell `tinyvamp` to estimate it:

``` r
gamma_tilde <- matrix(0, nrow = 1, ncol = 1)
gamma_tilde_fixed_indices <- matrix(FALSE, nrow = 1, ncol = 1)
```

Finally, we want to tell `tinyvamp` to estimate the contamination
profile, which we will initialize uniformly over all taxa. Because there
is only one specimen in this study (sampled 9 times), we do need to
choose one taxon that we believe is *not* a contaminant. We’re going to
choose *L. fermentum*, the 248th taxon, but we could have chosen any
taxon in the mock (or, if we didn’t have a mock community, any taxon
more common in the undiluted samples).

``` r
P_tilde <- matrix(c(rep(1/(jj - 1), jj - 1), 0), ncol = jj, nrow = 1)
P_tilde_fixed_indices <- matrix(c(rep(FALSE, jj - 1), TRUE), ncol = jj, nrow = 1)
```

Finally, we are not interested in estimating detection efficiencies for
the contaminants, so we set $`\tilde{\mathbf{X}} = \mathbf{0}`$:

``` r
X_tilde <- matrix(0, nrow = 1, ncol = 1) 
```

## Fitting the model

With all of that hard work parameterizing the model out of the way, we
can fit the model! The function we use is `estimate_parameters`, and we
give it all of the information we just put together. Also, we’re going
to use the reweighted estimator, which we generally recommend because it
has lower variance when the noise in our data is not Poisson-like… that
is, always!

``` r
full_karstens_model <- estimate_parameters(W = W,
                                           X = X,
                                           Z = Z,
                                           Z_tilde = Z_tilde,
                                           Z_tilde_gamma_cols = Z_tilde_gamma_cols,
                                           gammas = gammas,
                                           gammas_fixed_indices = gammas_fixed_indices,
                                           P = P,
                                           P_fixed_indices = P_fixed_indices,
                                           B = B,
                                           B_fixed_indices = B_fixed_indices,
                                           X_tilde = X_tilde,
                                           P_tilde = P_tilde,
                                           P_tilde_fixed_indices = P_tilde_fixed_indices,
                                           gamma_tilde = gamma_tilde,
                                           gamma_tilde_fixed_indices = gamma_tilde_fixed_indices,
                                           criterion = "reweighted_Poisson") 
```

Wooohooo! That took only a minute or two on my computer, which is pretty
impressive given how many parameters there are in this model. The output
is a list of named elements, which you can see as follows:

``` r
full_karstens_model %>% names
#>  [1] "optimization_status"       "objective"                
#>  [3] "varying"                   "fixed"                    
#>  [5] "W"                         "X"                        
#>  [7] "Z"                         "Z_tilde"                  
#>  [9] "Z_tilde_gamma_cols"        "gammas"                   
#> [11] "gammas_fixed_indices"      "P"                        
#> [13] "P_fixed_indices"           "B"                        
#> [15] "B_fixed_indices"           "X_tilde"                  
#> [17] "P_tilde"                   "P_tilde_fixed_indices"    
#> [19] "gamma_tilde"               "gamma_tilde_fixed_indices"
#> [21] "alpha_tilde"               "Z_tilde_list"             
#> [23] "criterion"                 "weights"                  
#> [25] "variance_function"
```

Let’s check everything worked okay:

``` r
full_karstens_model$optimization_status
#> [1] "Converged"
```

Terrific! The estimated parameter values are in the object `varying`:

``` r
full_karstens_model$varying %>% as_tibble
#> # A tibble: 264 × 4
#>       value param       k     j
#>       <dbl> <chr>   <dbl> <dbl>
#>  1 0.0123   P_tilde     1     1
#>  2 0.0333   P_tilde     1     2
#>  3 0.0603   P_tilde     1     3
#>  4 0.0109   P_tilde     1     4
#>  5 0.0220   P_tilde     1     5
#>  6 0.00670  P_tilde     1     6
#>  7 0.00502  P_tilde     1     7
#>  8 0.00350  P_tilde     1     8
#>  9 0.00731  P_tilde     1     9
#> 10 0.000752 P_tilde     1    10
#> # ℹ 254 more rows
```

Wow, 264 parameters were estimated pretty quickly! That’s thanks to *a
lot* of hard work and cleverness on David’s part. Lets look specifically
at the estimated efficiencies, and I will join it to the taxon names
data to make it easier to interpret

``` r
estd_bs <- full_karstens_model$varying %>% as_tibble %>% dplyr::filter(param == "B")
estd_bs %>%
  inner_join(tibble(j=241:247, name = genera_data[mock_taxa][-8]), by="j")
#> # A tibble: 7 × 5
#>   value param     k     j name                
#>   <dbl> <chr> <dbl> <dbl> <chr>               
#> 1 -3.72 B         1   241 Pseudomonas         
#> 2 -2.78 B         1   242 Escherichia-Shigella
#> 3 -2.90 B         1   243 Salmonella          
#> 4 -4.96 B         1   244 Enterococcus        
#> 5 -5.47 B         1   245 Staphylococcus      
#> 6 -4.40 B         1   246 Listeria            
#> 7 -4.08 B         1   247 Bacillus
```

On the basis of this model, we estimate that in an equal mixture of
*Pseudomonas aeruginosa* and our reference taxon, *L. fermentum*, we
expect on average to observe exp(-3.72) = 0.02 *P. aeruginosa* reads for
each *L. fermentum* read. Since all of the estimated $`\beta`$’s were
negative, this tells us that *L. fermentum* is the most easily detected
taxon in the mock community.

Aside: The estimated values in this table are pretty similar to those
from our paper, with differences being due to different subsets of the
data being used for fitting. Here we used all nine samples.

You can look at which other parameters were estimated as follows:

``` r
full_karstens_model$varying %>% pull(param) %>% unique
#> [1] "P_tilde"     "B"           "gamma"       "gamma_tilde"
```

To obtain 95% confidence intervals for the detection efficiencies, you
can run the following. This step can take a while (proportional to
`n_boot`), but can be parallelized across your computer’s cores using
the package `parallel`. Using only one core, 10 bootstrap iterations
took \<6 minutes on my 1.4 GHz Dual-Core Intel Core i7 from 2017. So, on
a more modern machine using parallelization, even 100 iterations should
take even less than this. Also, it’s good for you to go for a walk and
stretch your neck occasionally – enjoy!

``` r
full_cis <- bootstrap_ci(W,
                         alpha = 0.05, # to get 95% CIs
                         fitted_model = full_karstens_model,
                         n_boot = 500, # 100 would probably be fine
                         parallelize = TRUE,
                         ncores = 6,
                         seed  = 4324)
```

You can use the following to pull out the confidence intervals

``` r
full_cis$ci %>%
  filter(param == "B")
full_cis$ci %>%
  filter(param == "P_tilde") # %>% tail(20)
```

Please note that confidence intervals (even clever subsampled
bootstrapped CIs like this one) inherently rely a large sample to have
correct coverage. So while this is a useful uncertainty quantification
tool, I wouldn’t put much store in the 95-percent-ness of these
intervals, because they are only based on 9 samples.

Side note: If you look at
`full_karstens_model$varying %>% filter(param == "P_tilde")`, you’ll
notice that all of the mock taxa had some non-zero estimated relative
abundance in the contamination profile. This is not surprising, as point
estimates for maximum likelihood estimates are rarely on the boundary of
the parameter space. However, if you run 500 bootstraps with
`seed = 4324` you’ll find that 5 out of the 7 of these taxa had 95%
confidence intervals that contain zero. Go tinyvamp!

## Wrap-up

I hope you found this to be helpful in setting up the `tinyvamp` model
to analyze your dilution series! Please let me (Amy) know if you have
any questions or would like further clarification. You can do this via a
GitHub issue (preferred), or via email.

Happy estimating and experiment-designing!

## References

- Karstens, L., Asquith, M., Davin, S. et al. Controlling for
  Contaminants in Low-Biomass 16S rRNA Gene Sequencing Experiments.
  *mSystems* **4**(4):e00290-19 (2019). doi: 10.1128/mSystems.00290-19.
- Clausen, D.S. and Willis, A.D. Modeling complex measurement error in
  microbiome experiments. *arxiv* 2204.12733.
  <https://arxiv.org/abs/2204.12733>.
