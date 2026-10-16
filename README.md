# Solving PDEs with Random Walks

A small, self-contained notebook that approximates the solution of the heat equation with **Monte Carlo simulation of random walks**, and compares it against the exact solution given by convolution with the heat kernel.

It was originally built as a teaching demo for an undergraduate PDE course, as a concrete first look at the bridge between **partial differential equations and probability**.

![animation_N](https://github.com/user-attachments/assets/724171df-9008-47c9-b606-479bfd8cf768)

## The problem

We consider the Cauchy problem for the heat equation on the whole real line:

$$
\begin{cases}
u_t = \kappa\, u_{xx}, & x \in \mathbb{R},\ t > 0, \\
u(x,0) = g(x), & x \in \mathbb{R}.
\end{cases}
$$

Its solution is the convolution of the initial datum with the heat kernel:

$$
u(x,t) = \int_{\mathbb{R}} \frac{1}{\sqrt{4\pi\kappa t}}\, e^{-\frac{(x-y)^2}{4\kappa t}}\, g(y)\, dy .
$$

## The probabilistic idea

The heat kernel is exactly the density of a Gaussian with variance $2\kappa t$. So, if $B_t$ is a standard Brownian motion,

$$
u(x,t) = \mathbb{E}\left[\, g\big(x + \sqrt{2\kappa}\, B_t\big) \right].
$$

This is the simplest instance of the **Feynman–Kac formula**: solving the PDE amounts to computing an expectation over random paths.

To turn this into an algorithm:

1. **Replace Brownian motion by a simple random walk.** With $K$ time steps of size $\Delta t = T/K$, each particle jumps $\pm\sigma$ with probability $1/2$, where $\sigma = \sqrt{2\kappa\,\Delta t}$. After $k$ steps its displacement has variance $2\kappa t_k$, and by the central limit theorem (Donsker's principle) the walk converges to $\sqrt{2\kappa}\,B_t$.
2. **Replace the expectation by an empirical mean.** Launch $N$ independent walks from each grid point $x_i$ and average:

$$
u(x_i, t_k) \approx \frac{1}{N} \sum_{j=1}^{N} g\big(X^{(j)}_k\big), \qquad X^{(j)}_0 = x_i .
$$

By the law of large numbers this converges as $N \to \infty$, and by the CLT the statistical error decays like $N^{-1/2}$.

Note that no spatial discretization of the derivatives is involved: the spatial grid is only the set of points where we evaluate $u$, and each point is computed independently of the others.

## What's in the notebook

`Heat_MonteCarlo.ipynb` contains:

- **`solve_heat_mc(g, ...)`** — vectorized Monte Carlo solver. Returns the space grid, the time grid and the matrix `U_mc[i, k] ≈ u(x_i, t_k)`.
- **`solve_heat_exact(g, x_grid, t_grid, ...)`** — reference solution via discrete convolution with the heat kernel on an extended mesh (the domain is padded by three standard deviations of the kernel).
- An example initial condition (a sum of two Gaussians), which you can replace with any vectorized function `g`.
- Comparisons between both solutions:
  - space–time heatmaps of the Monte Carlo solution, the exact solution and the absolute error;
  - mean absolute error as a function of time;
  - experiments with $N = 10, 100, 1000, 10000$ walks per point;
  - animations of $u(\cdot, t)$, Monte Carlo vs exact.

### Main parameters of `solve_heat_mc`

| Parameter | Meaning | Default |
|---|---|---|
| `g` | initial condition (vectorized function) | — |
| `x_min`, `x_max` | spatial window where $u$ is evaluated | `-10`, `10` |
| `t_max` | final time $T$ | `1` |
| `K` | number of time steps | `100` |
| `M` | number of spatial grid points | `100` |
| `N` | number of random walks per grid point | `1000` |
| `kappa` | diffusion coefficient $\kappa$ | `0.5` |
| `seed` | random seed (for reproducibility) | `None` |

## Quick start

```bash
git clone https://github.com/cesaramr314/Solving-PDEs-with-Random-Walks.git
cd Solving-PDEs-with-Random-Walks
pip install numpy matplotlib jupyter
jupyter notebook Heat_MonteCarlo.ipynb
```

Minimal usage inside the notebook:

```python
import numpy as np

g = lambda x: 2*np.exp(-(x + 2)**2) + 0.5*np.exp(-(x - 1)**2)

x_grid, t_grid, U_mc = solve_heat_mc(g, N=10_000, seed=42)
U_exact = solve_heat_exact(g, x_grid, t_grid)

print("MAE:", np.mean(np.abs(U_mc - U_exact)))
```

> **Note:** GitHub does not render animations created with `matplotlib.animation.FuncAnimation`. Run the notebook locally (or in Google Colab) to see them.

## Results

The Monte Carlo solution converges to the exact one as the number of walks grows, with error of order

$$
\text{error} \sim N^{-1/2}.
$$

Multiplying $N$ by 100 reduces the error by roughly a factor of 10. This is the characteristic Monte Carlo rate: slow, but independent of the spatial dimension, which is why probabilistic methods become attractive for high-dimensional problems where grid-based schemes are too expensive.

## Going further

This notebook is a first step into a broad area connecting PDEs and probability: the Feynman–Kac formula for equations with potentials and source terms, diffusions with variable coefficients, boundary value problems via stopped processes, and hydrodynamic limits where macroscopic PDEs emerge from microscopic particle systems.

## License

See [LICENSE](LICENSE).



