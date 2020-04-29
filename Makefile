CFLAGS = -O3 -g -march=skylake-avx512 -ftree-vectorize -fopt-info-vec-optimized -fopt-info-loop-optimized
CFLAGS =  -g 

all: shcentroids dfs_timing

shcentroids: shcentroids.o shfunctions.o sh-fft.o
	cc  shcentroids.o shfunctions.o sh-fft.o -lfftw3 -o shcentroids 

shfunctions.o: shfunctions.c sh.h

sh-fft.o: sh-fft.c sh.h

shcentroids.o: shcentroids.c sh.h

fields.txt: agwsvalid-makefields.py
	python agwsvalid-makefields.py > fields.txt

dfs_timing: dfs_computeffts.o dfs_timing.o
	g++ $(CFLAGS) dfs_computeffts.o dfs_timing.o -lfftw3 -o dfs_timing

dfs_preprocess.o: dfs_preprocess.cpp
	g++ -c -I /home/bmcleod/include $(CFLAGS) dfs_preprocess.cpp -o dfs_preprocess.o

dfs_timing.o: dfs_timing.c
	g++ -c -I /home/bmcleod/include $(CFLAGS) dfs_timing.c -o dfs_timing.o

dfs_computeffts.o: dfs_computeffts.c
	g++ -c -I /home/bmcleod/include $(CFLAGS) dfs_computeffts.c -o dfs_computeffts.o

#dfs_preprocess: dfs_preprocess.cpp
#	g++ -I /home/bmcleod/include -O3 -mavx512f -mavx512cd -fopt-info-vec-all -march=native dfs_preprocess.cpp -o dfs_preprocess

test:	fields.txt
	./agwsvalid-time

#	cc -O3 -o shcentroids shcentroids.c
#	cc -O3 -g -march=skylake-avx512 -ftree-vectorize -fopt-info-vec-optimized -fopt-info-loop-optimized -o shcentroids shcentroids.c
#	gcc -O3 -g -march=native -fopt-info-loop-optimized -o shcentroids shcentroids.c


vectest: vectest.c
#	cc -O3 -o vectest vectest.c
	cc -O3 -g -march=skylake-avx512 -ftree-vectorize -fopt-info-vec-optimized -fopt-info-loop-optimized -o vectest vectest.c

clean:
	rm shcentroids vectest dfs_preprocess *.o

prof:
	sudo operf ./shcentroids
	opannotate --source --output-dir=annotated
	opannotate --source --assembly > assembly.txt
	less annotated/home/bmcleod/src/agwsbench/shcentroids.c

