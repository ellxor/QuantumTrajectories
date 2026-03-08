## Weakly Permutation-Symmetric Quantum Trajectories
> [!WARNING]
> This library is work-in-progress and is not yet suited for research work.

The idea of this library is to provide an efficient way to study weakly permutation symmetric states of two-level systems (spin-1/2) in an optical cavity with O(N) scaling for N systems. This is done using the results of [Barberena](https://doi.org/10.48550/arXiv.2508.05751) to reformulate the Lindblad equation in terms of purely collective operators.

The systems that are currently provided are:
- Dicke Model: $H_D = \omega_c a^\dagger a + \omega_0 \sum_i \sigma^z_i + g \sum_i \sigma^x_i (a^\dagger + a)$
- Tavis-Cummings Model: $H_{TC} = \omega_c a^\dagger a + \omega_0 \sum_i \sigma^z_i + g \sum_i (\sigma^-_i a^\dagger + \sigma^+_i a)$

with Lindblad $\dot{\rho} = -i[H,\rho] + \kappa D[a] + \gamma_\phi \sum_i D[\sigma^z] + \gamma_\downarrow D[\sigma^-]$,
where $D[x] = x \rho x^\dagger - \frac{1}{2} (x^\dagger x \rho + \rho x^\dagger x )$

It uses a mixture of quantum jumps and quantum state diffusion to preserve and break symmetries respectively, where useful. It has been tested on, and able to reproduce both:
- [Kirton paper](https://doi.org/10.1103/PhysRevLett.118.123602) on supressing and restoring Dicke superradiance from 2017.
- [Freter paper](https://doi.org/10.1515/nanoph-2025-0427) on organic superradiance from 2025.

**Goals:**
- [x] Scale linearly with number of atoms.
- [x] Study large system size, to the order of $N^5$.
- [x] Compile quickly (under a second) as future wrapper libraries might require on-the-fly code compilation for performance.
- [ ] Add Python/Julia wrappers for this library to make it easier for others to use.

**Build:**

Requires a C23 compliant compiler to build (GCC works, clang does not). Use `src/parameters.h` to set model and trajectory options. Recommended flags to build and run:
```
gcc -std=c23 -o main src/main.c -lm -O3 -ffast-math -march=native -Wall -Wno-psabi
./main
```

**TODO:**
- [x] Add both Dicke and Tavis-Cummings models.
- [ ] Add optimisation for Tavis-Cummings model starting in symmetry-unbroken state.
- [ ] Add graphs of latest results to README.
- [ ] Refactor code to make symmetry-breaking cleaner.
- [ ] Refactor main function and add more detailed logging / output.
- [ ] Explore using GPU acceleration to speed up code.
- [ ] Add timeout variable to stop simulation at a certain time theshold (useful for university servers).
	- [ ] Dump data to file after each round as a backup, to prevent data loss if process is killed.
- [ ] Add Python/Julia wrapper libraries.
	- [ ] Include a way for user to define their own jumps and Hamiltonian using operator objects.
	- [ ] Publish!
