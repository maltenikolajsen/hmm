# HMM

### Purpose and goal

The purpose of this package is to implement methods for handling Hidden Markov Models in a sensible way.

The goal is to define a HMM as a sensible S3 or S4 object, and define various methods for such objects.
Possible methods include computing likelihood, random sampling, fitting data as a HMM, decoding hidden states, forecasting, predicting.
These methods should be implemented using classic methods from HMM theory, such as the forwards-, backwards- and viterbi-algorithm for general distribution.

The interface should be somewhat similar to that of lm.

### Authors

This project is made by:
- [https://github.com/asbjornholk](Asbjørn Holk Thomsen)
- [https://github.com/maltenikolajsen](Malte Nikolajsen)
- [https://github.com/LeanderTilsted](Leander Tilsted Kristensen)
