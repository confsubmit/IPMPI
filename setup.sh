PROFILER=ipmpi
mkdir $PWD/lib
mkdir -p ${HOME}/${PROFILER}_profiles/
echo "export LD_PRELOAD=$PWD/lib/lib${PROFILER}.so" > .${PROFILER}
cp .${PROFILER} ${HOME}
make clean
make install
echo "Initial Setup Finished"