all: main.cpp
	c++ -o fvs.o -std=c++11 -stdlib=libc++ -I ../lemon-1.3.1/ -I ../lemon-1.3.1/build/ main.cpp reductions.cpp -DLEMON_ONLY_TEMPLATES -DNEWCEIL

clean:
	rm fvs.o