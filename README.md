Directed Rod Simulation for Drosophila Class IV Dendritic Arbor Morphogenesis

This repository contains a C implementation of a simulation model that recapitulates the mean-field picture of Drosophila Class IV dendritic arbor morphogenesis using *directed rod dynamics.

Model Description
The simulation models the steady-state tip density of growing dendritic arbors as a collection of directed rods within a square box of size L = 200 μm with periodic boundary conditions.

Each rod represents a dendritic segment with:
A base located randomly and uniformly within the simulation box.
A tip that grows in a random direction θ ∈ [0, 2π].

Tip Dynamics
Each rod's tip stochastically transitions among three dynamic states: Growing  (G), Paused  (P) and Shrinking (S)

Transitions follow experimentally measured rates see Table 1 in Ouyang et al., 2025). The simulation uses a Monte Carlo scheme for state transitions
Upon contact, rods disappear instantly, mimicking experimental observations that dendritic branches shrink after collision (Shree et al., 2022).
Lateral branching is implemented as birth of a new rod at a random location and orientation within the box.

For each existing tip, a branching probability is calculated: 𝑃𝑖_b = 1 − 𝑒−𝑙𝑖(𝑡)𝑘𝑏∆

Compilation
To compile the code, use a standard C compiler. For example:
gcc -O2 -o simulate_rods simulate_rods.c -lm

Once compiled, you can run the simulation as:
./simulate_rods

References
Ouyang et al. (2025). Title of preprint. bioRxiv. https://doi.org/10.1101/2025.02.24.639873
Shree et al. (2022). Experimental data and dynamics of dendritic tip shrinkage.

This is an examplery simulation snapshot:
![SnapShot_24hr](https://github.com/user-attachments/assets/5a8e83cc-6931-4822-87ee-7f52497f9275)
