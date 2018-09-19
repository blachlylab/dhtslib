all: deps
deps:
	cd htslib;autoreconf;./configure CFLAGS="-DBGZF_MT" CPPFLAGS="-DBGZF_MT" --disable-bz2 --disable-lzma;make -j8;cd ..;
