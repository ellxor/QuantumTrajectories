set -xe

CC="gcc"
CFLAGS="-ffast-math -O3 -march=native -Wall -Wextra"

for N in 2 4 6 8 10 12 14 16 32 64 128; do
    echo "Running test for N = ${N}..."
    ${CC} -o main main.cpp ${CFLAGS} -DOPTION_N=${N}
    ./main
done
