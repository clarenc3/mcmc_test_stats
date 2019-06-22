# Have some ROOT dependencies
ROOT_LIB := $(shell root-config --libs)
ROOT_CFLAGS := $(shell root-config --cflags)

MCStats: MCStats.cpp
	g++ -g ${ROOT_CFLAGS} -o $@ $^ ${ROOT_LIB}

clean:
	rm MCStats
