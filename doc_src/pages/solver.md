title: PCG Solver

The following is an overview of the observation-space preconditioned conjugate gradient solver, the design of which was borrowed heavly from the Navy's NCODA system [^1][^2]


The 3DVar formulation used here is the observation space formulation:

\( \mathbf{x}_a - \mathbf{x}_b = [\mathbf{HBH}^T+\mathbf{R}]^{-1}\mathbf{H}^T
  \mathbf{R}^{-1}(\mathbf{y}-H(\mathbf{x}_b)) \)

which is solved in two step. First is to solve via iteration the following for \(\mathbf{z}\)

\( [\mathbf{HBH}^T+\mathbf{R}]\mathbf{z} = \mathbf{y}-H(\mathbf{x}_b) \)

then finding the analysis increments by:

\( \mathbf{x}_a - \mathbf{x}_b = \mathbf{BH}^T\mathbf{z} \)


This is all handled by several modules in the code. The main conjugate gradient solver is in the [[g3dv_pcg]] module, the background covariance model is handled by [[g3dv_bgcov]] and the separation of observations in to blocks by [[g3dv_obs]].

# Preconditioned conjugate gradient solver

with \( \mathbf{A} = [\mathbf{HBH}^T+\mathbf{R}] \) and \(\mathbf{d} = \mathbf{y}-H(\mathbf{x}_b) \) iteratively solve for \(\mathbf{z}\) as follows until either a specified number of iterations have occurred or the residual, \(\mathbf{r}_k\) has reduced from the initial value, \( \mathbf{r}_0 \) by a specified factor. The formation of the preconditioning matrix, \(\mathbf{A}^{*}\) is described later.

*initialize the following variables:*

> \( \mathbf{z}_0 = 0\)

> \( \mathbf{r}_0 = \mathbf{y}-H(\mathbf{x}_b) \)

> \( \mathbf{s}_0 = \mathbf{A}^{*-1}\mathbf{r}_0 \)

> \( \mathbf{p}_1 = \mathbf{s}_0 \)

*begin iterations until solution converges:*

*if \(k > 1\):*

>\( \beta_k = \cfrac{ \mathbf{r}_{k-1}^T \mathbf{s}_{k-1} }{ \mathbf{r}_{k-2}^T \mathbf{s}_{k-2} }\)

>\( \mathbf{p}_k = \mathbf{s}_{k-1} + \beta_k\mathbf{p}_{k-1} \)

*then for each iteration with \(k>0\):*

> \( \mathbf{q}_k = \mathbf{Ap}_k \)

> \(\alpha_k = \cfrac{ \mathbf{s}_{k-1}^T\mathbf{r}_{k-1} }{\mathbf{p}_{k}^T\mathbf{q}_{k} } \)

> \(\mathbf{z}_k = \mathbf{z}_{k-1} + \alpha_k\mathbf{p}_k \)

> \(\mathbf{r}_k = \mathbf{r}_{k-1} + \alpha_k\mathbf{q}_k \)

> \( \mathbf{s}_k = \mathbf{A}^{*-1}\mathbf{r}_k \)

It should be noted that in the actual code only the most recent vectors for \(\mathbf{r}_k\) and \(\mathbf{s}_k\) are stored, and so only the final dot products of \(\mathbf{r}^T\mathbf{s}\) are kept from 2 steps ago for calculating \(\beta\).

# Preconditioning

TODO


[^1]: Daley, R., & Barker, E. (2001). NAVDAS Source Book 2001. NRL Publication NRL. PU/7530—01-441, 161pp.

[^2]: Cummings, J. A., & Smedstad, O. M. (2013). Variational data assimilation for the global ocean. In Data Assimilation for Atmospheric, Oceanic and Hydrologic Applications (Vol. II) (pp. 303-343). Springer Berlin Heidelberg.

