#include <algorithm>
#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <threads.h>

using namespace std::complex_literals;

//=================================================================================================
// COMMAND LINE OPTIONS
//=================================================================================================
#ifndef OPTION_N
#  define OPTION_N 32
#endif
//=================================================================================================


//=================================================================================================
// MAIN OPTIONS
//=================================================================================================
constexpr int N = OPTION_N;
constexpr int CAVITY_LIMIT = std::max(3*N/8, 32);
static_assert(N % 2 == 0, "N must be even!");

constexpr float KAPPA = 1.0;
constexpr float GAMMA_PHI = 0.1;
constexpr float GAMMA_DOWN = 0.2;
constexpr float OMEGA_C = 1.0;
constexpr float OMEGA_0 = 0.5;
constexpr float G_COUPLING = 0.9 / sqrtf(N);

constexpr float TIME_DELTA = 1e-3;
constexpr int SIMULATION_ITERATIONS = 10'000;
constexpr int LOGGING_ITERATIONS = 16; // log every this number of steps (ideally power of 2)

static_assert(SIMULATION_ITERATIONS % LOGGING_ITERATIONS == 0,
	"total iterations must be an integer multiple of the logging iteration steps");
//=================================================================================================

// An (M,a) vector constrained to a single J sector.
// The size of the Hilbert space is (N + 1) * CAVITY_LIMIT, however padding is added to
// either side of the vector of coefficients. This is done for performance so that shifts
// by M +- 1 (and a +- 1) will always fall inside a valid memory region. This eliminates the
// need for if statements to check if shifted elements are in-bounds. This is also aligned
// to a 64-byte boundary for SIMD friendliness.
struct WaveVector {
	alignas(64) std::complex<float> coeffs[(N + 3)*CAVITY_LIMIT];
};

constexpr size_t begin(int j_sector) {
	return (N/2 - j_sector + 1)*CAVITY_LIMIT;
}

// Compute Effective Hamiltonian efficiently by using shifts rather than a matrix.
// This term is also multiplied by -i dt.
void compute_effective_hamiltonian_step(const WaveVector &src, WaveVector &dst, int j_sector) {
	auto index = begin(j_sector);
	auto final = index + (2*j_sector + 1)*CAVITY_LIMIT;

	// Important! zero destination wave first as operations only add (due to overlapping).
	std::fill(dst.coeffs + index, dst.coeffs + final, 0);

	for (int m = -j_sector; m <= j_sector; ++m) {
		for (int a = 0; a < CAVITY_LIMIT; ++a) {
			// First the T dagger T jump terms: these summed over produce a neat
			// closed form, verified using sympy, coded below.
			dst.coeffs[index] -= src.coeffs[index] * (TIME_DELTA/2)
			                  * (KAPPA*a + GAMMA_PHI*N + 4*GAMMA_DOWN*(m + N/2));

			// Next the Hamiltonian terms which result in a 5-way branching
			auto coeff = 1.0if * TIME_DELTA * src.coeffs[index];
			dst.coeffs[index] -= coeff * (OMEGA_C * a + 2*OMEGA_0 * m);

			if (a) {
				dst.coeffs[index + CAVITY_LIMIT - 1] -= coeff * G_COUPLING * sqrtf((j_sector + m + 1) * (j_sector - m) * a);
				dst.coeffs[index - CAVITY_LIMIT - 1] -= coeff * G_COUPLING * sqrtf((j_sector - m + 1) * (j_sector + m) * a);
			}

			if (a < CAVITY_LIMIT - 1) {
				dst.coeffs[index + CAVITY_LIMIT + 1] -= coeff * G_COUPLING * sqrtf((j_sector + m + 1) * (j_sector - m) * (a + 1));
				dst.coeffs[index - CAVITY_LIMIT + 1] -= coeff * G_COUPLING * sqrtf((j_sector - m + 1) * (j_sector + m) * (a + 1));
			}

			++index;
		}
	}
}


void evolve_under_effective_hamiltonian(WaveVector &wave, int j_sector) {
	constexpr int RUNGE_KUTTA_POLY = 4; // order of integration step

	// In this case of an exponential and linear Hamiltonian, the Runge-Kutta method
	// is identical to a Taylor series expansion, so this is performed for efficiency.

	WaveVector _a, _b; // create two temporary wave vectors as a double-buffering technique.
	WaveVector *a = &wave, *b = &_b; // these are controlled by pointers which are cheap to swap.

	int factorial = 1;

	for (int i = 1; i <= RUNGE_KUTTA_POLY; ++i) {
		// The branch in the argument should be eliminated by loop unrolling
		compute_effective_hamiltonian_step(*a, *b, j_sector); // b now contains -i Heff dt a
		factorial *= i;

		// accumulate Taylor series expansion
		auto index = begin(j_sector);
		float factor = 1.0 / factorial;

		for (int m = -j_sector; m <= j_sector; ++m) {
			for (int a = 0; a < CAVITY_LIMIT; ++a) {
				wave.coeffs[index] += factor * b->coeffs[index];
				++index;
			}
		}

		if (i == 1) a = &_a; // a is temporarily set to &wave for first iteration to avoid a copy
		std::swap(a, b); // perform float buffering
	}
}


// Collective J operators
constexpr float precompute_E(int J) { return (float)(N/2 + 1)     / (2*J*(J + 1));         }
constexpr float precompute_F(int J) { return (float)(N/2 + J + 2) / (2*(J + 1)*(2*J + 3)); }
constexpr float precompute_G(int J) { return (float)(N/2 - J + 1) / (2*J*(2*J - 1));       }


void run_simulation(size_t thread_index) {
	// xorshift-* peudo-random 64-bit number algorithm
	constexpr uint64_t MAGIC = 0x2545f4914f6cdd1d;
	uint64_t seed = thread_index ^ MAGIC;

	auto rand = [&seed]() {
		seed ^= seed >> 12;
		seed ^= seed << 25;
		seed ^= seed >> 27;
		return seed * MAGIC;
	};

	char filename[100];
	std::snprintf(filename, sizeof filename, "data/test-log-%d-thread-%zu.txt", N, thread_index);
	auto *log = fopen(filename, "wb");

	// start with maximum allowed J (free choice)
	int j_sector = N/2;

	// start in state M = 0, a = max
	WaveVector wave = {};
	wave.coeffs[begin(j_sector) + CAVITY_LIMIT - 1] = 1;

	enum {
		JUMP_PHOTON_ANNIHILATION,
		JUMP_DEPHASING_SAME_SPIN,
		JUMP_DEPHASING_LOWER_SPIN,
		JUMP_DEPHASING_UPPER_SPIN,
		JUMP_SPIN_LOSS_SAME_SPIN,
		JUMP_SPIN_LOSS_LOWER_SPIN,
		JUMP_SPIN_LOSS_UPPER_SPIN,

		JUMP_COUNT,
		EFFECTIVE_HAMILTONIAN = JUMP_COUNT,
	};

	float factor_of_E = precompute_E(j_sector);
	float factor_of_F = precompute_F(j_sector - 1);
	float factor_of_G = precompute_G(j_sector + 1);
	float jump_table[JUMP_COUNT] = {};

	for (size_t i = 0; i < SIMULATION_ITERATIONS; ++i) {
		auto r = rand() / (float)UINT64_MAX;

		float total_probability = 0;
		size_t choice = 0;

		for (float probability : jump_table) {
			if (r < (total_probability += probability)) break;
			++choice;
		}

		switch (choice) {
		case JUMP_PHOTON_ANNIHILATION: {
			auto index = begin(j_sector);

			for (int m = -j_sector; m <= j_sector; ++m) {
				for (int a = 0; a < CAVITY_LIMIT; ++a) {
					wave.coeffs[index - 1] = wave.coeffs[index] * sqrtf(a);
					++index;
				}
			}

			wave.coeffs[index] = 0; // last element must be 0
			break;
		}

		case JUMP_DEPHASING_SAME_SPIN: {
			auto index = begin(j_sector);

			for (int m = -j_sector; m <= j_sector; ++m) {
				for (int a = 0; a < CAVITY_LIMIT; ++a) {
					wave.coeffs[index] *= m;
					++index;
				}
			}

			break;
		}

		case JUMP_DEPHASING_LOWER_SPIN: {
			auto index = begin(j_sector);

			for (int m = -j_sector; m <= j_sector; ++m) {
				for (int a = 0; a < CAVITY_LIMIT; ++a) {
					wave.coeffs[index] *= sqrtf((j_sector - m)*(j_sector + m));
					++index;
				}
			}

			j_sector -= 1;
			factor_of_E = j_sector == 0 ? 0 : precompute_E(j_sector);
			factor_of_F = j_sector == 0 ? 0 : precompute_F(j_sector - 1);
			factor_of_G = precompute_G(j_sector + 1);
			break;
		}

		case JUMP_DEPHASING_UPPER_SPIN: {
			auto index = begin(j_sector);

			for (int m = -j_sector; m <= j_sector; ++m) {
				for (int a = 0; a < CAVITY_LIMIT; ++a) {
					wave.coeffs[index] *= sqrtf((j_sector - m + 1)*(j_sector + m + 1));
					++index;
				}
			}

			j_sector += 1;
			factor_of_E = precompute_E(j_sector);
			factor_of_F = precompute_F(j_sector - 1);
			factor_of_G = precompute_G(j_sector + 1);

			// zero out new region
			auto start = begin(j_sector);
			std::fill(wave.coeffs + start, wave.coeffs + start + CAVITY_LIMIT, 0);
			std::fill(wave.coeffs + index, wave.coeffs + index + CAVITY_LIMIT, 0);
			break;
		}

		case JUMP_SPIN_LOSS_SAME_SPIN: {
			auto index = begin(j_sector - 1);

			for (int m = -j_sector + 1; m <= j_sector; ++m) {
				for (int a = 0; a < CAVITY_LIMIT; ++a) {
					wave.coeffs[index - CAVITY_LIMIT] = sqrtf((j_sector - m + 1)*(j_sector + m)) * wave.coeffs[index];
					++index;
				}
			}

			std::fill(wave.coeffs + index - CAVITY_LIMIT, wave.coeffs + index, 0);
			break;
		}

		case JUMP_SPIN_LOSS_LOWER_SPIN: {
			auto index = begin(j_sector - 2);

			for (int m = -j_sector + 2; m <= j_sector; ++m) {
				for (int a = 0; a < CAVITY_LIMIT; ++a) {
					wave.coeffs[index - CAVITY_LIMIT] = sqrtf((j_sector + m - 1)*(j_sector + m)) * wave.coeffs[index];
					++index;
				}
			}

			j_sector -= 1;
			factor_of_E = j_sector == 0 ? 0 : precompute_E(j_sector);
			factor_of_F = j_sector == 0 ? 0 : precompute_F(j_sector - 1);
			factor_of_G = precompute_G(j_sector + 1);
			break;
		}

		case JUMP_SPIN_LOSS_UPPER_SPIN: {
			auto index = begin(j_sector);

			for (int m = -j_sector; m <= j_sector; ++m) {
				for (int a = 0; a < CAVITY_LIMIT; ++a) {
					wave.coeffs[index - CAVITY_LIMIT] = sqrtf((j_sector - m + 1)*(j_sector - m + 2)) * wave.coeffs[index];
					++index;
				}
			}

			// zero out 2 upper sectors!
			std::fill(wave.coeffs + index - CAVITY_LIMIT, wave.coeffs + index + CAVITY_LIMIT, 0);

			j_sector += 1;
			factor_of_E = precompute_E(j_sector);
			factor_of_F = precompute_F(j_sector - 1);
			factor_of_G = precompute_G(j_sector + 1);
			break;
		}

		case EFFECTIVE_HAMILTONIAN:
			evolve_under_effective_hamiltonian(wave, j_sector);
			break;

		default:
			__builtin_unreachable();
		}

		// Recompute wave norm and jump table for next iteration.
		std::fill(jump_table, jump_table + JUMP_COUNT, 0);

		auto index = begin(j_sector);
		float inner_product = 0;

		for (int m = -j_sector; m <= j_sector; ++m) {
			for (int a = 0; a < CAVITY_LIMIT; ++a) {
				auto norm = std::norm(wave.coeffs[index]);
				++index;

				inner_product += norm;
				jump_table[JUMP_PHOTON_ANNIHILATION] += norm * a;
				jump_table[JUMP_DEPHASING_SAME_SPIN] += norm * m*m;
				jump_table[JUMP_DEPHASING_LOWER_SPIN] += norm * (j_sector - m)    *(j_sector + m);
				jump_table[JUMP_DEPHASING_UPPER_SPIN] += norm * (j_sector - m + 1)*(j_sector + m + 1);
				jump_table[JUMP_SPIN_LOSS_SAME_SPIN]  += norm * (j_sector - m + 1)*(j_sector + m);
				jump_table[JUMP_SPIN_LOSS_LOWER_SPIN] += norm * (j_sector + m - 1)*(j_sector + m);
				jump_table[JUMP_SPIN_LOSS_UPPER_SPIN] += norm * (j_sector - m + 1)*(j_sector - m + 2);
			}
		}

		// Renormalise state and jump table.
		auto scale = 1.0 / sqrtf(inner_product);
		index = begin(j_sector);

		for (int m = -j_sector; m <= j_sector; ++m) {
			for (int a = 0; a < CAVITY_LIMIT; ++a) {
				wave.coeffs[index] *= scale;
				++index;
			}
		}

		for (float &probability : jump_table) {
			probability *= scale*scale;
		}

		// Log out the expectation of a whilst it is equal to the non-rated probability of photon annihilation.
		if (i % LOGGING_ITERATIONS == LOGGING_ITERATIONS - 1) {
			std::fprintf(log, "%f\n", jump_table[JUMP_PHOTON_ANNIHILATION] / N);
		}

		// Scale jump table by process rates.
		jump_table[JUMP_PHOTON_ANNIHILATION]  *= TIME_DELTA * KAPPA;
		jump_table[JUMP_DEPHASING_SAME_SPIN]  *= TIME_DELTA * 4*GAMMA_PHI  * factor_of_E;
		jump_table[JUMP_DEPHASING_LOWER_SPIN] *= TIME_DELTA * 4*GAMMA_PHI  * factor_of_F;
		jump_table[JUMP_DEPHASING_UPPER_SPIN] *= TIME_DELTA * 4*GAMMA_PHI  * factor_of_G;
		jump_table[JUMP_SPIN_LOSS_SAME_SPIN]  *= TIME_DELTA * 4*GAMMA_DOWN * factor_of_E;
		jump_table[JUMP_SPIN_LOSS_LOWER_SPIN] *= TIME_DELTA * 4*GAMMA_DOWN * factor_of_F;
		jump_table[JUMP_SPIN_LOSS_UPPER_SPIN] *= TIME_DELTA * 4*GAMMA_DOWN * factor_of_G;
	}

	fclose(log);
}


int start_simulation_thread(void *thread_index) {
	run_simulation((size_t)thread_index);
	return 0;
}


double get_time_from_os() {
	timespec timestamp;
	clock_gettime(CLOCK_REALTIME, &timestamp);

	return timestamp.tv_sec + 1.0e-9 * timestamp.tv_nsec;
}


int main() {
	constexpr size_t REPEATS = 50;
	constexpr size_t THREAD_COUNT = 12;

	thrd_t threads[THREAD_COUNT];
	size_t index = 0;

	auto t1 = get_time_from_os();
	for (size_t i = 0; i < REPEATS; ++i) {
		for (auto& thread : threads) thrd_create(&thread, start_simulation_thread, (void *)index++);
		for (auto& thread : threads) thrd_join(thread, nullptr);
		std::fprintf(stderr, "Round %zu/%zu complete.\n", i + 1, REPEATS);
	}
	auto t2 = get_time_from_os();

	std::fprintf(stderr, "%.0f iterations per second per thread.\n", REPEATS * SIMULATION_ITERATIONS / (t2 - t1));
	std::fprintf(stderr, "%.0f iterations per second.\n", REPEATS * THREAD_COUNT * SIMULATION_ITERATIONS / (t2 - t1));
}
