HASTS 416 – Tutorial 1 Results

## A1 – 5-State Markov Chain

**State Classification**
- Recurrent class: {S1} (absorbing)
- Transient classes: {S2}, {S3, S4, S5}

**Periods**
- S1: 1
- S2: transient (no return)
- S3, S4, S5: 3

**Steady-State Distribution**
- All probability mass eventually concentrates in S1: (1, 0, 0, 0, 0)

**Simulated Trajectories** (30 steps)
- Starting in S1 stays in S1.
- Starting in S2 moves to S1 and stays.
- Starting in S2 (another run) follows S2 → S5 → S4 → S3 → S5 → … and never leaves the transient class {S3,S4,S5}.

**Convergence of Unconditional Probabilities**
- Probability of S1 → 1; all others → 0 by n ≈ 20.


## A2 – 7-State Markov Chain

**State Classification**
- Recurrent classes: {S1, S2} (period 2), {S4, S5, S6, S7} (aperiodic)
- Transient class: {S3}

**Periods**
- S1, S2: 2
- S3: transient
- S4–S7: 1

**Stationary Distributions**
- For class {S1, S2}: (0.5, 0.5, 0, 0, 0, 0, 0)
- The chain is not ergodic because two recurrent classes exist.

**Simulated Trajectories** (50 steps)
- Starting in S2: oscillates between S1 and S2 (period 2).
- Starting in S5: wanders in {S4,S5,S6,S7} for a while, then eventually enters {S1,S2} and begins alternating.

## A3 – Time-Inhomogeneous Traffic Chain

**Analytical Distribution at 6 PM** (starting from Light at 1 PM)
Monte Carlo Simulation (N = 10,000):
Light   Heavy   Jammed
0.0147  0.1283  0.8570
