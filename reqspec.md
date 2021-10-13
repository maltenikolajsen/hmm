## Tentative contents of the package

### Functions

- Forward algorithm for HMM's
- Backward algorithm for HMM's
- EM-algorithm for non-stationary HMM's, where users can either choose from standard distributions (Poisson, Normal etc.) or provide custom density and estimation functions
- Direct likelihood optimization for HMM's
- Functions for prediction: Decoding, outlier detection, forecasting etc.

### Classes

Mainly one `hmm` class similar to that of `lm`, including many of the same methods such as `plot`, `summary`, `predict`, `simulate` etc.

### Data

Different simulated data-sets using different marginal distributions.

Perhaps also fitting real-world data, if any such is freely available.