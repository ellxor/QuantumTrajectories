#define _GNU_SOURCE
#include <complex.h>
#include <math.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdio.h>
#include <threads.h>
#include <time.h>
#include <sys/sysinfo.h>

#include "random.h"
#include "parameters.h"
#include "utility.h"

static_assert(N % 2 == 0, "Number of atoms in optical cavity must be even.");
static_assert(SimulationIterations % LoggingStep == 0,
    "The logging step size must wholly divide the total number of simulation iterations");
static_assert(TrajectoryCount % ThreadCount == 0,
	"The total number of trajectories must be an integer multiple of the number of allocated threads.");

// Structure to hold auxiliary information about the current wave state.
typedef struct WaveInfo WaveInfo;

struct WaveInfo {
	int j_sector;		// J sector that state resides in.
	complex float alpha;	// parameter of photon displacement
	complex float jx;	// expectation of Jx operator
	complex float xi;	// complex Wiener fluctuation for quantum state diffusion
};

// Stores a wave vector for a quantum state, with dimension N+1 (for magnetisation M = -J to J) by
// CavityLimit which is the truncation for photons. It is padded with an extra cavity limit at each
// end to reduce the number of necessary boundary checks.
typedef complex float WaveVector[(N + 3)*CavityLimit];

// Starting index for WaveVector in a given J sector.
size_t begin(int j_sector) { return (N/2 - j_sector + 1) * CavityLimit; }


enum {
	JUMP_DEPHASING_SAME_J,	// Jz
	JUMP_DEPHASING_LOWER_J,	// Lz
	JUMP_DEPHASING_UPPER_J,	// Kz
	JUMP_SPIN_LOSS_SAME_J,	// (iJx + Jy)
	JUMP_SPIN_LOSS_LOWER_J,	// (iLx + Ly)
	JUMP_SPIN_LOSS_UPPER_J, // (iKx + Ky)
	JUMP_PHOTON_LOSS,	// a

	JUMP_COUNT,
	EFFECTIVE_HAMILTONIAN = JUMP_COUNT,
};

// Collective J operators
float precompute_E(int J) { return (float)(N/2 + 1)     / (2*J*(J + 1));         }
float precompute_F(int J) { return (float)(N/2 + J + 2) / (2*(J + 1)*(2*J + 3)); }
float precompute_G(int J) { return (float)(N/2 - J + 1) / (2*J*(2*J - 1));       }


// Selects a jump (or effective Hamiltonian) at random using probabilities in the jump table.
int select_random_jump(float jump_table[]) {
	float r = random_uniform();
	float seen_probability = 0;

	for (int choice = 0; choice < JUMP_COUNT; ++choice) {
		seen_probability += jump_table[choice];
		if (r < seen_probability) return choice;
	}

	// If no jump has been selected then we evolve using the effective Hamiltonian.
	return EFFECTIVE_HAMILTONIAN;
}


/*
 * Implementation of Effective Hamiltonian
 */

// Compute linear evolution of effective Hamiltonian using shifts rather than a matrix.
// [-i Heff dt + L dxi]
void linear_hamiltonian_integration_step(WaveVector src, WaveVector dst, WaveInfo info) {
	auto index = begin(info.j_sector);
	auto final = index + (2*info.j_sector + 1)*CavityLimit;

	// Important! zero destination wave first as operations only add (due to overlapping).
	fill(dst + index, dst + final, 0);

	for (int m = -info.j_sector; m <= info.j_sector; ++m) {
		for (int a = 0; a < CavityLimit; ++a) {
			// First the T dagger T jump terms
			dst[index] -= src[index] * 0.5f*TimeStep * (Kappa*a + GammaPhi*N + GammaDown*(m + N/2));

			// Hamiltonian terms
			auto coeff = I * TimeStep * src[index];
			dst[index] -= coeff * (OmegaC * (a + cnormf(info.alpha)) + 2*Omega0 * m);

			dst[index + CavityLimit] -= coeff * GCoupling *       info.alpha  * sqrtf((info.j_sector + m + 1) * (info.j_sector - m));
			dst[index - CavityLimit] -= coeff * GCoupling * conjf(info.alpha) * sqrtf((info.j_sector - m + 1) * (info.j_sector + m));

			dst[index + CavityLimit] -= coeff * GCoupling * conjf(info.alpha) * sqrtf((info.j_sector + m + 1) * (info.j_sector - m)) * (model == Dicke);
			dst[index - CavityLimit] -= coeff * GCoupling *       info.alpha  * sqrtf((info.j_sector - m + 1) * (info.j_sector + m)) * (model == Dicke);

			if (a) {
				dst[index + CavityLimit - 1] -= coeff * GCoupling * sqrtf((info.j_sector + m + 1) * (info.j_sector - m) * a);
				dst[index - CavityLimit - 1] -= coeff * GCoupling * sqrtf((info.j_sector - m + 1) * (info.j_sector + m) * a) * (model == Dicke);
				dst[index - 1] += coeff * GCoupling * conjf(info.jx) * sqrtf(a);

				if (UseQuantumStateDiffusion) {
					dst[index - 1] += info.xi * src[index] * sqrtf(Kappa * a);
				}
			}

			if (a < CavityLimit - 1) {
				dst[index - CavityLimit + 1] -= coeff * GCoupling * sqrtf((info.j_sector - m + 1) * (info.j_sector + m) * (a + 1));
				dst[index + CavityLimit + 1] -= coeff * GCoupling * sqrtf((info.j_sector + m + 1) * (info.j_sector - m) * (a + 1)) * (model == Dicke);
				dst[index + 1] += coeff * GCoupling * info.jx * sqrtf(a + 1);
			}

			++index;
		}
	}
}


// Compute higher order integration step of evolution of effective Hamiltonian.
// [-i Heff dt + L dxi]
void hamiltonian_integration_step(WaveVector wave, WaveInfo info) {
	constexpr int RUNGE_KUTTA_POLY = 4; // order of integration step

	if (UseQuantumStateDiffusion) {
		info.xi = random_complex_gaussian(TimeStep);
	}

	// In this case of an exponential and linear Hamiltonian, the Runge-Kutta method
	// is identical to a Taylor series expansion, so this is performed for efficiency.

	static thread_local WaveVector _a, _b; // Create two temporary wave vectors as a double-buffering technique.
	auto a = wave; // Controlled by pointers which are cheap to swap.
	auto b = _b;

	int factorial = 1;

	for (int i = 1; i <= RUNGE_KUTTA_POLY; ++i) {
		linear_hamiltonian_integration_step(a, b, info); // b now contains -i Heff dt a
		factorial *= i;

		// accumulate Taylor series expansion
		auto index = begin(info.j_sector);
		float factor = 1.0 / factorial;

		for (int m = -info.j_sector; m <= info.j_sector; ++m) {
			for (int a = 0; a < CavityLimit; ++a) {
				wave[index] += factor * b[index];
				++index;
			}
		}

		if (i == 1) a = _a; // a is temporarily set to wave for first iteration to avoid a copy
		swap(a, b); // perform double-buffering
	}
}


/*
 * Implementation of Quantum Jumps:
 * The coefficients of jumps are ignored as the states will be renormalised afterwards.
 * The phase is also ignored as it doesn't determine the dynamics (Lindblad jumps are invariant under phase rotations).
 */

void jump_dephasing_same_j(WaveVector wave, WaveInfo *info) {
	auto index = begin(info->j_sector);

	for (int m = -info->j_sector; m <= info->j_sector; ++m) {
		for (int a = 0; a < CavityLimit; ++a) {
			wave[index] *= m;
			++index;
		}
	}
}

void jump_dephasing_lower_j(WaveVector wave, WaveInfo *info) {
	auto index = begin(info->j_sector);

	for (int m = -info->j_sector; m <= info->j_sector; ++m) {
		for (int a = 0; a < CavityLimit; ++a) {
			wave[index] *= sqrtf((info->j_sector - m)*(info->j_sector + m));
			++index;
		}
	}

	info->j_sector -= 1;
}

void jump_dephasing_upper_j(WaveVector wave, WaveInfo *info) {
	auto index = begin(info->j_sector);

	for (int m = -info->j_sector; m <= info->j_sector; ++m) {
		for (int a = 0; a < CavityLimit; ++a) {
			wave[index] *= sqrtf((info->j_sector - m + 1)*(info->j_sector + m + 1));
			++index;
		}
	}

	info->j_sector += 1;

	// Zero out new regions that the wave state has expanded into.
	auto start = begin(info->j_sector);
	fill(wave + start, wave + start + CavityLimit, 0); // zero out M = -J
	fill(wave + index, wave + index + CavityLimit, 0); // zero out M = +J
}

void jump_spin_loss_same_j(WaveVector wave, WaveInfo *info) {
	auto index = begin(info->j_sector - 1);

	for (int m = -info->j_sector + 1; m <= info->j_sector; ++m) {
		for (int a = 0; a < CavityLimit; ++a) {
			wave[index - CavityLimit] = sqrtf((info->j_sector - m + 1)*(info->j_sector + m)) * wave[index];
			++index;
		}
	}

	// Zero out M = +J, (no states can be lowered into here).
	fill(wave + index - CavityLimit, wave + index, 0);
}

void jump_spin_loss_lower_j(WaveVector wave, WaveInfo *info) {
	auto index = begin(info->j_sector - 2);

	for (int m = -info->j_sector + 2; m <= info->j_sector; ++m) {
		for (int a = 0; a < CavityLimit; ++a) {
			wave[index - CavityLimit] = sqrtf((info->j_sector + m - 1)*(info->j_sector + m)) * wave[index];
			++index;
		}
	}

	info->j_sector -= 1;
}

void jump_spin_loss_upper_j(WaveVector wave, WaveInfo *info) {
	auto index = begin(info->j_sector);

	for (int m = -info->j_sector; m <= info->j_sector; ++m) {
		for (int a = 0; a < CavityLimit; ++a) {
			wave[index - CavityLimit] = sqrtf((info->j_sector - m + 1)*(info->j_sector - m + 2)) * wave[index];
			++index;
		}
	}

	info->j_sector += 1;
	// Zero out upper two sectors (no states can be lowered to M = +J AND expanding into M = J + 1)
	fill(wave + index - CavityLimit, wave + index + CavityLimit, 0);
}

void jump_photon_loss(WaveVector wave, WaveInfo *info) {
	auto index = begin(info->j_sector);

	for (int m = -info->j_sector; m <= info->j_sector; ++m) {
		for (int a = 0; a < CavityLimit; ++a) {
			wave[index - 1] = wave[index] * sqrtf(a);
			++index;
		}
	}

	wave[index] = 0; // last element must be 0
}


/*
 * Quantum Trajectory Implementation
 */

constexpr int LogFileLines = SimulationIterations / LoggingStep;
static atomic(float) expectation_of_a[LogFileLines] = {};
static atomic(float) expectation_of_b[LogFileLines] = {};
static atomic(float) expectation_of_m[LogFileLines] = {};
static atomic(float) expectation_of_j[LogFileLines] = {};
static atomic(float) expectation_of_m2[LogFileLines] = {};


complex float compute_expectation_of_jx(WaveVector wave, WaveInfo info) {
	auto index = begin(info.j_sector);
	complex float jx = 0;

	for (int m = -info.j_sector; m <= info.j_sector; ++m) {
		for (int a = 0; a < CavityLimit; ++a) {
			jx += conjf(wave[index - CavityLimit]) * wave[index] * sqrtf((info.j_sector - m + 1) * (info.j_sector + m));
			jx += conjf(wave[index + CavityLimit]) * wave[index] * sqrtf((info.j_sector + m + 1) * (info.j_sector - m)) * (model == Dicke);
			++index;
		}
	}

	return jx;
}


void run_simulation() {
	// Start simulation in empty state [J = Jmax, M = -J, a = 0]
	WaveInfo info = {
		.j_sector = N/2,
		.alpha = 0,
	};

	static thread_local WaveVector wave = {};

	if (model == Dicke) {
		wave[begin(N/2)] = 1.0f;
	}

	if (model == TavisCummings) {
		constexpr float LogCos = logf(cosf(0.5f * SbAngle));
		constexpr float LogSin = logf(sinf(0.5f * SbAngle));

		for (int m = -info.j_sector; m <= info.j_sector; ++m) {
			auto index = begin(0) + m*CavityLimit;
			auto k = m + info.j_sector;

			auto logp = 0.5f * ( lgammaf(N + 1) - lgammaf(k + 1) - lgammaf(N - k + 1) ) + k * LogCos + (N - k) * LogSin;
			complex float p = expf(logp);

			switch ((N - k) % 4) {
				case 1: p *= -I; break;
				case 2: p *= -1; break;
				case 3: p *=  I; break;
			}

			wave[index] = p;
		}
	}

	float factor_of_E = precompute_E(info.j_sector);
	float factor_of_F = precompute_F(info.j_sector - 1);
	float factor_of_G = precompute_G(info.j_sector + 1);

	// Start with empty jump table for simplicity (i.e. choosing Heff for first step).
	float jump_table[JUMP_COUNT] = {};

	for (int step = 0; step < SimulationIterations; ++step) {
		if (UseDisplacementTransform) info.jx = compute_expectation_of_jx(wave, info);
		int choice = select_random_jump(jump_table);
		bool update_prefactors = false;

		switch (choice) {
			case JUMP_DEPHASING_SAME_J:  jump_dephasing_same_j(wave, &info); break;
			case JUMP_DEPHASING_LOWER_J: jump_dephasing_lower_j(wave, &info); update_prefactors = true; break;
			case JUMP_DEPHASING_UPPER_J: jump_dephasing_upper_j(wave, &info); update_prefactors = true; break;
			case JUMP_SPIN_LOSS_SAME_J:  jump_spin_loss_same_j(wave, &info); break;
			case JUMP_SPIN_LOSS_LOWER_J: jump_spin_loss_lower_j(wave, &info); update_prefactors = true; break;
			case JUMP_SPIN_LOSS_UPPER_J: jump_spin_loss_upper_j(wave, &info); update_prefactors = true; break;
			case JUMP_PHOTON_LOSS:
				if (UseQuantumStateDiffusion) { jump_photon_loss(wave, &info); break; }
				[[fallthrough]]; // if using purely jumps
			case EFFECTIVE_HAMILTONIAN:  hamiltonian_integration_step(wave, info); break;
		}

		if (update_prefactors) {
			factor_of_E = (info.j_sector == 0) ? 0.0f : precompute_E(info.j_sector);
			factor_of_F = (info.j_sector == 0) ? 0.0f : precompute_F(info.j_sector - 1);
			factor_of_G = precompute_G(info.j_sector + 1);
		}

		if (UseDisplacementTransform) {
			auto alpha_dot = -(I*OmegaC + 0.5f*Kappa) * info.alpha - I * GCoupling * info.jx;
			info.alpha += TimeStep * alpha_dot;
		}

		// Recompute wave norm and jump table for next iteration.
		for (int i = 0; i < JUMP_COUNT; ++i) jump_table[i] = 0;
		float inner_product = 0;
		float m_expectation = 0;

		auto index = begin(info.j_sector);

		for (int m = -info.j_sector; m <= info.j_sector; ++m) {
			for (int a = 0; a < CavityLimit; ++a) {
				auto norm = cnormf(wave[index]);
				++index;

				inner_product += norm;
				m_expectation += norm * m;
				jump_table[JUMP_DEPHASING_SAME_J]  += norm * m*m;
				jump_table[JUMP_DEPHASING_LOWER_J] += norm * (info.j_sector - m) * (info.j_sector + m);
				jump_table[JUMP_DEPHASING_UPPER_J] += norm * (info.j_sector - m + 1)*(info.j_sector + m + 1);
				jump_table[JUMP_SPIN_LOSS_SAME_J]  += norm * (info.j_sector - m + 1)*(info.j_sector + m);
				jump_table[JUMP_SPIN_LOSS_LOWER_J] += norm * (info.j_sector + m - 1)*(info.j_sector + m);
				jump_table[JUMP_SPIN_LOSS_UPPER_J] += norm * (info.j_sector - m + 1)*(info.j_sector - m + 2);
				jump_table[JUMP_PHOTON_LOSS] += norm * a;
			}
		}

		// Renormalise state and jump table.
		float scale = 1.0f / sqrtf(inner_product);
		index = begin(info.j_sector);

		for (int m = -info.j_sector; m <= info.j_sector; ++m) {
			for (int a = 0; a < CavityLimit; ++a) {
				wave[index] *= scale;
				++index;
			}
		}

		float scale_squared = scale*scale;
		for (int i = 0; i < JUMP_COUNT; ++i) jump_table[i] *= scale_squared;

		// Log out the expectation of a whilst it is equal to the non-rated probability of photon annihilation.
		if (step % LoggingStep == LoggingStep - 1) {
			auto write_index = step / LoggingStep;
			m_expectation *= scale_squared;

			expectation_of_a[write_index] += jump_table[JUMP_PHOTON_LOSS] + cnormf(info.alpha);
			expectation_of_b[write_index] += jump_table[JUMP_PHOTON_LOSS];
			expectation_of_m[write_index] += m_expectation;
			expectation_of_j[write_index] += info.j_sector;
			expectation_of_m2[write_index] += jump_table[JUMP_DEPHASING_SAME_J];
		}

		// Scale jump table by process rates.
		jump_table[JUMP_DEPHASING_SAME_J]  *= TimeStep * 4.0f*GammaPhi * factor_of_E;
		jump_table[JUMP_DEPHASING_LOWER_J] *= TimeStep * 4.0f*GammaPhi * factor_of_F;
		jump_table[JUMP_DEPHASING_UPPER_J] *= TimeStep * 4.0f*GammaPhi * factor_of_G;
		jump_table[JUMP_SPIN_LOSS_SAME_J]  *= TimeStep * GammaDown * factor_of_E;
		jump_table[JUMP_SPIN_LOSS_LOWER_J] *= TimeStep * GammaDown * factor_of_F;
		jump_table[JUMP_SPIN_LOSS_UPPER_J] *= TimeStep * GammaDown * factor_of_G;
		jump_table[JUMP_PHOTON_LOSS] *= TimeStep * Kappa;
	}
}


int start_simulation_thread(void *thread_index) {
	set_random_seed((uint64_t)thread_index);
	run_simulation();
	return 0;
}


double get_time_from_os() {
	struct timespec timestamp;
	clock_gettime(CLOCK_REALTIME, &timestamp);

	return timestamp.tv_sec + 1.0e-9 * timestamp.tv_nsec;
}


int main() {
	constexpr int Repeats = TrajectoryCount / ThreadCount;
	size_t core_count = get_nprocs();
	fprintf(stdout, "[INFO] %zu cpu cores detected.\n", core_count);

	if (ThreadCount > core_count) {
		fprintf(stdout, "[WARNING] running with more threads than cores available.\n");
	}

	thrd_t threads[ThreadCount];
	size_t index = 0;

	auto t1 = get_time_from_os();
	for (int i = 0; i < Repeats; ++i) {
		for (int j = 0; j < ThreadCount; ++j) thrd_create(&threads[j], start_simulation_thread, (void *)index++);
		for (int j = 0; j < ThreadCount; ++j) thrd_join(threads[j], nullptr);
		fprintf(stdout, "Round %d/%d complete.\n", i + 1, Repeats);
	}
	auto t2 = get_time_from_os();

	fprintf(stdout, "%.0f iterations per second per thread.\n", Repeats * SimulationIterations / (t2 - t1));
	fprintf(stdout, "%.0f iterations per second.\n", Repeats * ThreadCount * SimulationIterations / (t2 - t1));

	char filename[100];
	snprintf(filename, sizeof filename, "%s/log-N-%d.txt", OutputDataDirectory, N);
	fprintf(stdout, "Writing data to %s...\n", filename);

	auto log = fopen(filename, "wb");

	for (int i = 0; i < LogFileLines; ++i) {
		fprintf(log, "%g,%g,%g,%g,%g\n",
			expectation_of_a[i] / TrajectoryCount,
			expectation_of_b[i] / TrajectoryCount,
			expectation_of_m[i] / TrajectoryCount,
			expectation_of_j[i] / TrajectoryCount,
			expectation_of_m2[i] / TrajectoryCount
		);
	}

	fclose(log);
}
