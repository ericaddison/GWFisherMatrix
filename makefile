
CC=g++
CFLAGS = -g -W -Wall `gsl-config --cflags --libs`  
sources = GwaveTest.cpp GravWave.cpp fisher.cpp gslWrappers.cpp

all: $(sources)
	$(CC) $(sources) -o GTest $(CFLAGS)   

run: $(sources)
	$(CC)  $(sources) -o GTest $(CFLAGS)
	./GTest


