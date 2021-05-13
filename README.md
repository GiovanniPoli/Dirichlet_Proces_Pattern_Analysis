# Dirichlet process mixture model for Bayesian pattern analysis

Experiments designed with the aim of identifying bio-markers in a biologicalframework are usually characterized by a high number of measures and low sample size.This usually leads to inference based on separate multiple tests which imply an high risk of first-type errors. We propose a new approach that identifies patterns based on latent clusters among the measures to reduce that risk. This approach is based on a statisticalmodel composed of a parametric part for the known effects (due to the structure of the experiment) and a non-parametric part for the patterns. Non-parametric effects are defined through a Dirichlet process mixture model. 

This repository stores all the experimental results obtained using the new semi-parametric approach in simulation studies. 

# Dirichlet Process mixture models and pattern analysis

Let $\Theta$ be a finite-dimensional parameter space and $G\simDP(M,G_0)$ such that $M\in\mathhb{R}$ and $G_0$ be a probability measure defined on $\Theta$. A Dirichlet Process Mixture Model is p.d.f. described as:

\begin{equation}\label{DPMM}
f_G(y)=\int f(y|\theta)\ dG(\theta)
\end{equation}
