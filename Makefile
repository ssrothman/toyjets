test: test.o gaus.o gen.o
	g++ test.o gaus.o gen.o -o test -larmadillo

%.o: %.cc
	g++ -c -o $@ $< -I/usr/local/include/Minuit2 

clean:
	rm -f *.o test
