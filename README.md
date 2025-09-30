# Basket-Options-Pricing-Variance-Reduction-Techniques
This project investigates simulation-based methods for pricing vanilla basket options and barrier options, where no closed-form solutions exist and Monte Carlo techniques are therefore essential. We model the underlying assets in the basket as correlated geometric Brownian motions and apply Cholesky decomposition to generate correlated path trajectories. The option price is then estimated as the discounted expected payoff, computed via the sample mean across simulated paths.

Beyond the basic Monte Carlo approach, we implement and compare more advanced techniques to improve efficiency. These include variance reduction strategies such as antithetic variables and specialized control variates (including Pezzirelli-type control variates discussed in research literature, where we identified an overlooked flaw). When possible, we also combine multiple methods to further enhance performance, achieving comparable accuracy with significantly fewer simulations.

The project also explores exotic derivatives, specifically barrier options. Here, we analyze the convergence behavior of the Euler–Maruyama scheme, which exhibits relatively slow convergence (order 1/2). We then apply insights from Mike Giles’ Oxford lectures on Brownian Bridge methods to accelerate convergence to order 1. This part of the work deepened our understanding of numerical schemes for stochastic differential equations (SDEs), particularly the dual sources of error: discretization bias from approximating the SDE, and statistical error from Monte Carlo estimation.

For detailed problem statements, derivations, and complete solutions, please refer to the accompanying PDF in this repository.

<img width="1687" height="865" alt="PezzirelliVsMonteCarloEfficiency" src="https://github.com/user-attachments/assets/d05e847b-9dff-4cee-9fa7-09f8d448d150" />
<img width="875" height="656" alt="results_table" src="https://github.com/user-attachments/assets/274e1f7b-2a55-4568-a979-c99f3f4020cd" />

![barriererror](https://github.com/user-attachments/assets/d85cb536-ca74-4afd-9d75-25a13dac6b31)
