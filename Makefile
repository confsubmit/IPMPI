CC      = mpiicpc
COPT    = -O3 -shared -fPIC
CFLAGS  = -Wall -fPIC -no-pie -std=c++11
LD      = $(CC)
LDFLAGS = $(COPT)

PROFILER=ipmpi

MPI_DIR = $(PWD)/lib/

all: lib$(PROFILER).so

lib$(PROFILER).so: $(PROFILER).o helper.o
	$(CC) $(COPT) -o $@ $^

$(PROFILER).o: $(PROFILER).cpp
	$(CC) $(CFLAGS) -c $< -o $@

helper.o: helper.cpp
	$(CC) $(CFLAGS) -c $< -o $@

install: lib$(PROFILER).so
	cp lib$(PROFILER).so $(MPI_DIR)

.PHONY: all clean install

clean:
	rm -rf *.o *.a
