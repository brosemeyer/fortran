# LOCALLY
gfortran -o sus -fopenmp main.f03 SelfConsistent.f03 Energy.f03 Susceptibility.f03 Utils.f03 -llapack -lblas -lfftw3 && ./sus

# HYALITE
gfortran -o sus -fopenmp main.f03 SelfConsistent.f03 Energy.f03 Susceptibility.f03 Utils.f03 -L$LAPACK_DIR -L$FFTWDIR -L$BLASDIR -llapack -lfftw3 -lblas

ssh -X k58w857@hyalite.rcg.montana.edu

scp local_file.f03 k58w857@hyalite.rcg.montana.edu:/mnt/lustrefs/scratch/benjamin.rosemeyer/test/.
