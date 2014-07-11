all:
	h5pcc data_reorg.c -o data_reorg
	h5pcc e-bench.c    -o e-bench
	h5pcc data_expand.c -o data_expand

debug:
	h5pcc data_reorg.c  -o data_reorg  -g
	h5pcc e-bench.c     -o e-bench     -g
	h5pcc data_expand.c -o data_expand -g

clean:
	rm -f *.o e-bench data_reorg data_expand
