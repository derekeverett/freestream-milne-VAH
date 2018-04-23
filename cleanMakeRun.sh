#export OMP_NUM_THREADS=40 #set this to the number of threads on CPU
rm RunFreestreamMilne
rm -R output
mkdir output
cp freestream_input output/freestream_input
#make clean
make
#./RunFreestreamMilne
