all:
	h5pcc data_reorg.c -o data_reorg
	h5pcc e-bench.c    -o e-bench

debug:
	h5pcc data_reorg.c -o data_reorg -g
	h5pcc e-bench.c    -o e-bench    -g

clean:
	rm -f *.o e-bench data_reorg
