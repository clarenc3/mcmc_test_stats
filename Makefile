# Have some ROOT dependencies
ROOT_LIB := $(shell root-config --libs)
ROOT_CFLAGS := $(shell root-config --cflags)

all: MCStats MCStats_Minuit

MCStats: MCStats.cpp
	g++ -g ${ROOT_CFLAGS} -std=c++11 -o $@ $^ ${ROOT_LIB} 

MCStats_Minuit: MCStats_Minuit.cpp
	g++ -g ${ROOT_CFLAGS} -std=c++11 -o $@ $^ -Wl,--start-group ${ROOT_LIB} -Wl,--end-group

clean:
	rm MCStats
	rm MCStats_Minuit
