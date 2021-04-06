An implementation of the simulated annealing algorithm for clustering, written in C++ for use in R. This code is for demonstration purposes.

The simulated annealing algorithm gives you more guarantee that you will end up with the global optimum. However, its complexity is high. If speed is your concern, you're better off running the k-means algorithm several times with different initializations.

# Description of the files

`sa.R` contains the simulated annealing algorithm, written in R.

`sa.cpp` contains the simulated annealing algorithm, written in C++.

`sample.R` contains sample code, comparing the above methods with the original implementation of the k-means algorithm.

# Dependencies

`Rcpp` and `RcppArmadillo`

# References

Sanghamitra Bandyopadhyay, Ujjwal Maulik and Malay Kumar Pakhira.
“Clustering using simulated annealing with probabilistic redistribution”.
International Journal of Pattern Recognition and Artificial Intelligence 15.02 (2001), p. 269–285.
https://doi.org/10.1142/S0218001401000927

Zülal Güngör and Alper Ünler.
“K-harmonic means data clustering with simulated annealing heuristic”.
Applied Mathematics and Computation 184.2 (2007), p. 199–209.
https://doi.org/10.1016/j.amc.2006.05.166

Moon-Won Park et Yeong-Dae Kim. “A systematic procedure for setting parameters in simulated annealing algorithms”.
Computers & Operations Research 25.3 (1998), p. 207–217.
https://doi.org/10.1016/S0305-0548(97)00054-3

Shokri Z. Selim et K. Alsultan.
“A simulated annealing algorithm for the clustering problem”.
Pattern Recognition 24.10 (1991).
https://doi.org/10.1016/0031-3203(91)90097-O