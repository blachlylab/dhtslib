all: deps
deps:
	cd htslib;autoreconf;./configure --disable-bz2 --disable-lzma;make -j8;cd ..;
