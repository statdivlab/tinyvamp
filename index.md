# tinyvamp

[![Coverage
status](https://codecov.io/gh/statdivlab/tinyvamp/branch/main/graph/badge.svg)](https://codecov.io/github/statdivlab/tinyvamp?branch=main)

`tinyvamp` is a package for estimation and removal of measurement error
in high-throughput sequencing data.

There are many use cases for `tinyvamp`; here are a few: - Estimating
the detectability of different species using mock communities -
Comparing experimental protocols (or batches) with respect to how they
well they detect different species - Estimating “true” relative
abundances after adjusting for differential detectability - Estimating
“true” relative abundances after removing contamination relative
abundance profiles - Estimating contamination relative abundances

The online documentation is available
[here](https://statdivlab.github.io/tinyvamp).

Full details on the model and estimation method are available in the
[preprint](https://arxiv.org/abs/2204.12733).

Need something? File an
[issue](https://github.com/statdivlab/tinyvamp/issues) or [send Amy an
email](http://statisticaldiversitylab.com/team/amy-willis)!

## Installation

You can install the development version of tinyvamp from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("statdivlab/tinyvamp")
```

Are you using `tinyvamp`? We’d love to hear about it!

## Humans

Authors: [David
Clausen](https://www.biostat.washington.edu/people/david-clausen) and
[Amy Willis](http://statisticaldiversitylab.com)

Do you have a request for us? Let us know! We want folks to use
`tinyvamp` and will try to make it as easy as possible.

Do you have a question? Check out the above documentation list, then
shoot us an email or open a
[Discussion](https://github.com/adw96/breakaway/discussions). We receive
a lot of emails from users, so we try to answer questions on the Wiki
rather than responding to everyone individually.
