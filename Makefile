all:
	h5pcc e-bench.c -o e-bench -g

clean:
	rm -f *.o e-bench
