CFLAGS = -O3 -g -march=skylake-avx512 -ftree-vectorize -fopt-info-vec-optimized -fopt-info-loop-optimized

shcentroids: shcentroids.o shfunctions.o

shfunctions.o: shfunctions.c sh.h

shcentroids.o: shfunctions.c sh.h


#	cc -O3 -o shcentroids shcentroids.c
#	cc -O3 -g -march=skylake-avx512 -ftree-vectorize -fopt-info-vec-optimized -fopt-info-loop-optimized -o shcentroids shcentroids.c
#	gcc -O3 -g -march=native -fopt-info-loop-optimized -o shcentroids shcentroids.c


vectest: vectest.c
#	cc -O3 -o vectest vectest.c
	cc -O3 -g -march=skylake-avx512 -ftree-vectorize -fopt-info-vec-optimized -fopt-info-loop-optimized -o vectest vectest.c

clean:
	rm shcentroids vectest

prof:
	sudo operf ./shcentroids
	opannotate --source --output-dir=annotated
	opannotate --source --assembly > assembly.txt
	less annotated/home/bmcleod/src/agwsbench/shcentroids.c

