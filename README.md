# NSYSU CSE608 - Collective Intelligence: Particle Swarm Optimization Final Project

## Overview
This repository contains a baseline Particle Swarm Optimization (PSO) implementation and an improved variant with a simple local exploitation mechanism to mitigate premature convergence. Both variants are benchmarked on a suite of standard continuous optimization test functions (Ackley, Rastrigin, HappyCat, Rosenbrock, Zakharov, Michalewicz, Schwefel, BentCigar, DropWave, Step) across multiple dimensionalities. Generated logs include convergence traces, best-so-far fitness values, timing, and per‑iteration global best fitness for downstream plotting.

## Build Requirements
- C++11 compatible compiler (e.g., g++)
- Unix-like shell (tested on macOS / Linux)

## Quick Start
Build and run the baseline PSO on (example) Ackley 30D with 50 particles, 1 trial:
```bash
cd PSO
./run.sh Ackley 1 30 1 1.5 1.5 50
```
Arguments:
1. function_type  : One of Ackley|Rastrigin|HappyCat|Rosenbrock|Zakharov|Michalewicz|Schwefel|BentCigar|DropWave|Step
2. run            : Number of independent runs (integer)
3. dimension      : Dimensionality (int, controls iterations = dim * 10000)
4. k              : Velocity scaling factor for vMin/vMax = ± k * (xMax - xMin)/2
5. c1             : Cognitive coefficient
6. c2             : Social coefficient
7. numParticle    : Swarm size

Improved variant (adds adaptive local perturbation when global best stagnates):
```bash
cd PSO_improve
./run.sh Ackley 1 30 1 1.5 1.5 50
```

Batch all benchmark functions (default settings in `experiment.sh`):
```bash
cd PSO
./experiment.sh
```
(Repeat in `PSO_improve/` for the enhanced version.)

## Output Files
Each run writes per-iteration global best fitness values to a function-specific text file under a result subfolder. Examples:
- `PSO/result/coverage/Ackley_30D.txt`
- `PSO_improve/result/improve_coverage/Ackley_30D.txt`

Timing logs append to `run_time.txt` (baseline currently enabled; improve variant can be toggled on by uncommenting the line in its `main.cpp`).

`statistics.csv` (if populated) can aggregate summary metrics (min / mean / std) across runs or particle sizes (files like `stats_50P.csv` in nested `plot/` directory reflect processed statistics for plotting).

## Benchmark Functions & Domains
| Function | Domain | Notes |
|----------|--------|-------|
| Ackley | [-32.768, 32.768]^d | Many local minima |
| Rastrigin | [-5.12, 5.12]^d | Highly multimodal |
| HappyCat | [-20, 20]^d | Shift/rot invariant variant |
| Rosenbrock | [-10, 10]^d | Narrow curved valley |
| Zakharov | [-10, 10]^d | Continuous, convex-like |
| Michalewicz | [0, π]^d | Deceptive, steep valleys |
| Schwefel | [-500, 500]^d | Many distant minima |
| BentCigar | [-100, 100]^d | Ill-conditioned axis scaling |
| DropWave | [-5.12, 5.12]^d | Only meaningful in low D; still generalized here |
| Step | [-100, 100]^d | Flat plateaus |

## PSO Algorithm (Baseline)
For each iteration (T = dimension * 10000):
1. Update inertia weight linearly from wMax to wMin.
2. For each particle: update velocity with cognitive & social terms (clamped to vMin/vMax) then position (clamped to bounds).
3. Evaluate fitness; update personal and global bests.
4. Record current global best fitness to `gBest_list`.

## Improvements Implemented
In `PSO_improve/pso_improve.h` an additional stagnation detector tracks whether global best fitness has improved. If no improvement for `_threshold` (default 100) iterations, a lightweight local perturbation phase samples random particles and randomly perturbs selected dimensions (half of dimensions across a small loop). If a new personal best improves, it's adopted, potentially advancing the global best and resetting the stagnation counter.

## Parameter Tuning Tips
- Increase `numParticle` for harder multimodal landscapes (Ackley, Rastrigin, Schwefel).
- Adjust `k` if velocities saturate (frequent clipping) or search is too slow (small steps).
- Balance `c1` vs `c2`: higher `c1` encourages exploration around personal experiences; higher `c2` accelerates convergence to global best (risk of premature convergence).
- For high dimensions, consider reducing implicit iteration budget (currently `d * 10000`) to manage runtime or modify `_maxIter` in the constructor.


## Run Experiments
1. Run baseline and improved batch scripts.
2. Collect convergence files & timings.
3. Compute statistics (extend `cal_result.py` or add a new analysis script if needed).
4. Plot comparative curves (log-scale often helpful for clarity).

## Extending
To add a new benchmark function:
1. Implement `double NewFunc(const double* x, const int d)` in `pso.h` (and `pso_improve.h` if desired).
2. Add a domain branch in `main.cpp` mapping function name to pointer and bounds.
3. Invoke via `./run.sh NewFunc 1 30 1 1.5 1.5 50`.

To change stagnation behavior in improved variant, adjust `_threshold` or the local search strategy (e.g., perturbation magnitude, number of sampled particles, selection strategy).

## Known Limitations / Future Ideas
- Memory management: baseline uses raw pointers without deallocation (acceptable for short-lived process; could migrate to `std::vector` / RAII for safety).
- Missing seed control: add CLI param to set deterministic RNG seed.
- Local search step size currently uniform in [0,1]; could scale by domain size or anneal.
- Iteration budget fixed at `d * 10000`; expose as CLI arg for flexibility.
- Timing logging disabled in improved variant (line commented) — re-enable if needed.

