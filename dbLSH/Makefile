CC=g++ -O3 -std=c++11 -Wall -O3 -ffast-math -march=native
SRCS=$(wildcard ./src/*.cpp)
#INCLUDE = -I/home/qrstu/yty/wd/PM-LSH

dblsh:$(SRCS)
	rm -rf wlsh
	$(CC) $(SRCS) -o dblsh 

clean:
	rm -rf dblsh
	rm -rf $(wildcard *.txt)

cleancodes:
	rm -rf $(wildcard *.cpp) $(wildcard *.h) $(wildcard ./src/*.cpp) $(wildcard ./src/*.p)
