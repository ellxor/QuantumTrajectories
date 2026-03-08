/*
 * Parameters for Dicke and Tavis-Cummings models.
 *
 * Hamiltonian:
 *  - OmegaC:    for energy of photon mode [(a dagger)a]
 *  - Omega0:    for on-site energy of spins [Jz]
 *  - GCoupling: interaction strength of atom/cavity
 *               Dicke Model: [Jx((a dagger) + a)]
 *               Tavis-Cummings Model: [(J+)(a) + (J-)(a dagger)]
 *
 * Lindblad superoperator:
 *  - Kappa:     rate of photon loss
 *  - GammaDown: rate of spin loss
 *  - GammaPhi:  rate of spin dephasing
 *
 * and various other integration parameters below.
 */
#define _GNU_SOURCE
#include <math.h> // sqrtf

enum Model {
	Dicke,
	TavisCummings,
};

constexpr enum Model model = Dicke;

constexpr bool UseDisplacementTransform = true; // reduce cavity limit
constexpr bool UseQuantumStateDiffusion = true; // break photon symmetry (necessary to use displacement transform on symmetric state)

constexpr int N = 100;		// number of atoms in optical cavity
constexpr int CavityLimit = 4;  // truncation limit of displaced photon fluctuations

constexpr float SbAngle = M_PI_4f; // angle of symmetry breaking for initial state (Tavis-Cummings model only!)

constexpr float OmegaC = 1.0f;
constexpr float Omega0 = 0.5f;
constexpr float GCoupling = 0.9f / sqrtf(N);

constexpr float Kappa     = 1.0f;
constexpr float GammaDown = 0.2f;
constexpr float GammaPhi  = 0.0f;

constexpr int ThreadCount = 12;
constexpr int TrajectoryCount = 20 * ThreadCount;

constexpr int SimulationIterations = 40'000;	// number of integration steps
constexpr int LoggingStep = 32;			// number of integration steps to skip when logging data
constexpr float TimeStep = 0.001f;		// time step of integration
constexpr char OutputDataDirectory[] = "data";	// where to write data to
