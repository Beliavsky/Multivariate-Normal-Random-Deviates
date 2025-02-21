exec   = mvt_gfort.exe
obj    = kind.o constants.o linear_algebra.o statistics.o random.o mv_normal.o xmv_normal_test.o
FC     = gfortran
FFLAGS = -O0 -Wall -Werror=unused-parameter -Werror=unused-variable -Werror=unused-function -Wno-maybe-uninitialized -Wno-surprising -fbounds-check -static -g -fmodule-private

all: $(exec)

# Compile .f90 to .o
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

$(exec): $(obj)
	$(FC) -o $(exec) $(obj) $(FFLAGS)

run: $(exec)
	./$(exec)

clean:
	rm -f $(exec) $(obj)

