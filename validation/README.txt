Reproduce the comparison after installing ididPolicy and the paper dependencies:

    Rscript validation/compare_legacy.R 1000 1
    Rscript validation/compare_legacy.R 1500 1

The script uses replication 1 at each sample size and t = 1 and 4. It
independently fits the current two-part simulation and the package on the
same training and evaluation samples, with the same nuisance settings,
seeds, optimizer settings, and frozen reference. It compares IPW1,
IPW2 (DS in the simulation code), AIPW, and the logistic/Cox IPW1
comparator. For AIPW, it checks both evaluation of the known fixed optimum
and independent evaluation of the learned rule. The simulation code and
synthetic reference under legacy/ are used for reproduction.

The CSV files give the largest absolute difference by metric and
configuration. The package check also exercises fractional times,
arbitrary covariate names, custom learner wrappers, both survival engines,
and the IPW2 calculation at a censoring jump. These are implementation
checks on selected runs. They do not recompute the study's Monte Carlo
coverage or bias summaries.

Verified using R 4.5.1 with SuperLearner 2.0.29, survival 3.8.3,
randomForestSRC 3.4.5, lightgbm 4.6.0 and rgenoud 5.9.0.11.
