# Summary
The scripts in this repository can be used to generate summary statistics to use in Mendelian Randomization simulations

There are many different statistical methods available to perform Mendelian Randomization (MR). When these methods are introduced in the literature, simulations are performed to evaluate their performance. There is currently no standard for performing these simulations, hence some methods perform well under some settings but poorly under others.

The code in this repo is open source, meaning **you** can modify it but initiating a [pull request](https://github.com/noahlorinczcomi/simmr/pulls). We encourage the community to make changes to our code with the intention that the community will eventually settle on an accepted standard for how simulations in MR should be performed.

# Who is this for?
This software is primarily for researchers evaluating the performance of their new or existing MR estimator. This software can also be used by reviewers and others to validate the performance of pre-published manuscripts reporting MR simulations using GWAS summary statistics.

# What is in the software?
This software generates GWAS summary statistics under different scenarios of
- Exposure and outcome GWAS sample overlap (ref)[https://doi.org/10.1101/2021.06.28.21259622]
- Uncorrelated horizontal pleiotropy (ref)[https://doi.org/10.1093/ije/dyv080]
- Correlated horizontal pleiotropy (ref)[https://doi.org/10.1038/s41588-020-0631-4]
- Weak instruments (ref)[https://doi.org/10.1101/2023.01.10.523480]
- Winner's curse (ref)[https://doi.org/10.1101/2021.06.28.21259622]
- Linkage disequilibrium between instruments (ref)[https://doi.org/10.1002/gepi.22506]

that researchers may want to use in the evaluation of univariable or multivariable Mendelian Randomization methods.

# How do I use it?

# What settings can I change?

# Can I add ____ to the simulation?

