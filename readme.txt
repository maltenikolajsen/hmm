The purpose of this package is to implement methods for handling Hidden Markov Models in a sensible way.

The goal is to define a HMM as a sensible S3 or S4 object, and define various methods for such objects.
Possible methods include computing likelihood, random sampling, fitting data as a HMM, decoding hidden states, forecasting, predicting.
To implement these methods we must, it is necessary to implement classic methods such as the forwards-, backwards- and viterbi-algorithm for general distribution.

The interface should be somewhat similar to that of lm.

This project is made by:
Asbjørn Holk Thomsen     
Malte Nikolajsen
Leander Tilsted Kristensen