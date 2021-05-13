# Dirichlet process mixture model for Bayesian pattern analysis

Experiments designed with the aim of identifying bio-markers in a biologicalframework are usually characterized by a high number of measures and low sample size.This usually leads to inference based on separate multiple tests which imply an high risk of first-type errors. We propose a new approach that identifies patterns based on latent clusters among the measures to reduce that risk. This approach is based on a statisticalmodel composed of a parametric part for the known effects (due to the structure of the experiment) and a non-parametric part for the patterns. Non-parametric effects are defined through a Dirichlet process mixture model. 

This repository stores all the experimental results obtained using the new semi-parametric approach in simulation studies. 

# Dirichlet Process mixture models and pattern analysis

Let $\Theta$ be a finite-dimensional parameter space and G~DP(M,G_0) such that $M\in\mathhb{R}$ and $G_0$ be a probability measure defined on $\Theta$. A Dirichlet Process Mixture Model is p.d.f. described as:

![\begin{equation}
f_G(y)=\int f(y\gv\theta)\ dG(\theta)
\end{equation}
](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Cbegin%7Bequation%7D%0Af_G%28y%29%3D%5Cint+f%28y%5Cgv%5Ctheta%29%5C+dG%28%5Ctheta%29%0A%5Cend%7Bequation%7D%0A)


<img src=
"https://render.githubusercontent.com/render/math?math=%5Ctextstyle++%5Cbegin%7Bequation%7D%0AG+%7EDP%28M%2CG_0%29%0A%5Cend%7Bequation%7D" 
alt=" \begin{equation}
G ~DP(M,G_0)
\end{equation}">
