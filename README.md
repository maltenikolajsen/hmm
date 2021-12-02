# HMM

### Purpose and goal

The purpose of this package is to implement methods for handling Hidden Markov Models (HHMs) in a sensible way.

The package considers HHMs as S3 object.
An HMM class is created using the function `hmm`.
If emission observations are given, an HMM is fitted using the EM-algorithm.
Methods for the class include 
- `plot`
- `simulate`
- `summary`
Furthermore, the function `hmm` provides classical decoding of hidden states using posterior probabilities and the Viterbi algorithm.
What sets apart the package from other packages dealing with HMMs is the ease at which custom emission distributions can be added. 

### Authors

This project is made by:
- [Malte Nikolajsen](https://github.com/maltenikolajsen) (Project Manager)
- [Asbjørn Holk Thomsen](https://github.com/asbjornholk) (Quality Manager)
- [Leander Tilsted Kristensen](https://github.com/LeanderTilsted) (Documentation Manager)
