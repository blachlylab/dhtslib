all: htslib/libhts.a

htslib/hts.c:
	git submodule init
	git submodule update

htslib/libhts.a: htslib/hts.c
	cd htslib;make

clean:
	rm bgzfreader tabix_gffreader vcfwriter output.vcf
