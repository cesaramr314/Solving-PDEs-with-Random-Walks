# Solving PDEs with Random Walks!

This project implements a Monte Carlo method to approximate solutions of the heat equation via Simple Random Walks and compares it with the analytical solution via convolution with the heat kernel.

## Key Ideas

- Connection between Brownian Motion and the Heat Equation is well known:
  
  u(x,t) = E[g(x + B_t)]

- Monte Carlo approximation of the expected value using simple random walks
- Comparison with exact solution (Gaussian convolution)
- Visualization of convergence as N increases

## Features

- Interactive notebook
- Space-time heatmaps
- Animated solution evolution
- Convergence plots (MAE vs N)

## Results

Monte Carlo approximations converge to the exact solution at rate:

    error ~ N^{-1/2}

## Example

![animation_N](https://github.com/user-attachments/assets/724171df-9008-47c9-b606-479bfd8cf768)


```bash
pip install -r requirements.txt
jupyter notebook
