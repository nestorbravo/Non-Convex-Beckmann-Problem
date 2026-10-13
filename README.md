# Dynamic Beckmann problem with two interacting flows

Numerical experiments for two disjoint roads coupled through a congestion cost. Each road has an initial mass that must be discharged by a terminal time. The main quadratic experiments study the nonconvex regime, `gamma > 1`.

## Installation and execution

Install the dependencies and open the notebook:

```bash
python -m pip install numpy scipy matplotlib pyscipopt jupyterlab
jupyter lab Nonconvex_Beckmann_comentada.ipynb
```

Run the cells in order. To run an individual experiment, first execute the imports and the function definitions it uses. Some experiments reuse `potential`, so set it explicitly when rerunning a section out of order.

## Model and notation

| Code | Meaning |
|---|---|
| `mu_1`, `mu_2` | Initial masses |
| `eps_1`, `eps_2` | Linear holding coefficients |
| `T`, `N` | Terminal time and number of time cells |
| `delta` | Time step, `T / N` |
| `gamma` | Interaction coefficient |
| `j_1`, `j_2` | Discharge rates on the `N` cells |
| `m_1`, `m_2` | Pending masses at the `N + 1` time nodes |
| `w` | Total discharge rate, `j_1 + j_2` |

The discrete constraints are

$$
j_{\alpha,k}=\frac{m_{\alpha,k-1}-m_{\alpha,k}}{\delta}\geq0,
\qquad m_{\alpha,0}=\mu_\alpha,
\qquad m_{\alpha,N}=0.
$$

For the quadratic potential, the discrete energy is

$$
\mathscr H_N^{\mathrm{disc}}(m)
=\delta\sum_{k=1}^N\left[
\frac12(j_{1,k}^2+j_{2,k}^2+2\gamma j_{1,k}j_{2,k})
+\varepsilon_1m_{1,k}+\varepsilon_2m_{2,k}
\right].
$$

## Methods and example

- `solve_two_flows`: solves the discharge-rate formulation with SCIP. `time_limit` and `tolerance` control its stopping limits; `polish=True` adds a local SLSQP refinement.
- `solve_chang_zapeta`: implements successive local minimizations in the pending masses for the quadratic potential. `initial` accepts `"linear"`, `"early"`, or `"late"`; `tolerance` and `max_sweeps` control termination.
- `plot_solution`: displays discharge rates and pending masses, and prints energy and overlap diagnostics.

After executing the function definitions, compare both methods on the same parameters:

```python
gamma = 2.0
params = dict(mu_1=0.7, mu_2=0.3, T=1.0, N=40,
              eps_1=2.0, eps_2=1.5)

scip_solution = solve_two_flows(
    potential=partial(quadratic, gamma=gamma), **params,
)
cz_solution = solve_chang_zapeta(
    gamma=gamma, **params, initial="linear",
)

plot_solution(scip_solution)
plot_solution(cz_solution)
```

The supplied `potential(j_1, j_2, t)` contains only the congestion term; `solve_two_flows` adds the holding cost. A custom potential must support both NumPy arrays and SCIP symbolic variables.

## Experiments

The notebook includes an individual example, time-grid refinement, variation of `gamma`, a Chang-Zapeta example, a random comparison of both methods, and a nonquadratic power-law potential.

In the random comparison cell, change `m` to select the number of experiments. The seed `42` reproduces the parameter samples, and both methods receive the same parameters in each experiment. `results` stores scalar summaries; `solutions` stores the full solution dictionaries. Output is printed without pandas. Execution times can vary between runs.

## Reading the results

| Result key | Meaning |
|---|---|
| `discrete_energy`, `value` | Computed discrete energy |
| `energy` | Continuous action of the cellwise constant rates for the supplied time-independent congestion potentials |
| `correction` | Boundary correction between the two energies |
| `overlap` | `delta * sum(j_1 * j_2)` |
| `max_j_1_j_2_product` | Largest cellwise rate product |
| `overlap_cells`, `overlap_fraction` | Cells above `overlap_tolerance`, and their proportion |
| `mass_error`, `positivity_violation` | Numerical feasibility diagnostics |
| `seconds` | Elapsed time for the solver call |

For the supplied potentials, the exact energy relation is

$$
\mathscr H(j_N)=\mathscr H_N^{\mathrm{disc}}(m)
+\frac{\delta}{2}(\varepsilon_1\mu_1+\varepsilon_2\mu_2).
$$

Compare `discrete_energy` across methods using identical parameters. Zero overlap alone does not establish global optimality. SCIP returns its best feasible profile within the configured stopping limits; the current wrapper does not return its final optimality gap. For Chang-Zapeta, `converged` indicates stabilization of the energy, while `sweeps`, `max_change`, and `energy_history` describe the iteration.

Saved figures retain the labels from their original runs. Rerun the plotting cells to refresh them.

## Solver documentation

- [SCIP](https://www.scipopt.org/)
- [PySCIPOpt installation guide](https://pyscipopt.readthedocs.io/en/stable/install.html)
