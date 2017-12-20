SRC=mpi_spike.cpp

pfp: $(SRC)
	mpic++ -std=c++11 -fopenmp -Wall -fno-strict-overflow -ftrapv -fwrapv -fno-delete-null-pointer-checks -o pfp $(SRC)

clean:
	rm pfp
