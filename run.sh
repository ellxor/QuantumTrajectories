set -e

CC="gcc"
CFLAGS="-std=c++20 -ffast-math -O3 -march=native -Wall -Wextra"

echo "Starting green tests..."

for N in 2 4 6 8 10 16 32 64 128 256 512 1024; do
    echo "Running test for N = ${N}..."
    ${CC} -o main main.cpp ${CFLAGS} \
        -DOPTION_N=${N} \
        -DOPTION_PHI=0.1 \
        -DOPTION_DOWN=0.2 \
        -DOPTION_SIMULATION_ITERATIONS=40000 \
        -DOPTION_CAVITY_LIMIT="3*N/8" \
        -DOUTPUT_PREFIX=\"data/green-\"
    ./main
done

echo "Starting blue tests..."

for N in 2 4 6 8 10 16 32 64 128 256 512 1024; do
    echo "Running test for N = ${N}..."
    ${CC} -o main main.cpp ${CFLAGS} \
        -DOPTION_N=${N} \
        -DOPTION_PHI=0 \
        -DOPTION_DOWN=0.2 \
        -DOPTION_SIMULATION_ITERATIONS=40000 \
        -DOPTION_CAVITY_LIMIT=N/2 \
        -DOUTPUT_PREFIX=\"data/blue-\"
    ./main
done

echo "Starting red tests..."

for N in 2 4 6 8 10 16 32 64 128 256 512 1024; do
    echo "Running test for N = ${N}..."
    ${CC} -o main main.cpp ${CFLAGS} \
        -DOPTION_N=${N} \
        -DOPTION_PHI=0.01 \
        -DOPTION_DOWN=0 \
        -DOPTION_SIMULATION_ITERATIONS=120000 \
        -DOPTION_CAVITY_LIMIT=32 \
        -DOUTPUT_PREFIX=\"data/red-\"
    ./main
done

echo "Starting black tests..."

for N in 2 4 6 8 10 16 32 64 128 256 512 1024; do
    echo "Running test for N = ${N}..."
    ${CC} -o main main.cpp ${CFLAGS} \
        -DOPTION_N=${N}   \
        -DOPTION_PHI=0  \
        -DOPTION_DOWN=0 \
        -DOPTION_SIMULATION_ITERATIONS=40000 \
        -DOPTION_CAVITY_LIMIT=N \
        -DOUTPUT_PREFIX=\"data/black-\"
    ./main
done
