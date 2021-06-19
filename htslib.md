# htslib
### Install
#### Prequisites
`dhtslib` relies on `htslib`. `htslib` relies on a handful of compression and web-access libraries: `zlib, bzip2, lzma, curl, and ssl`. Technically `htslib` can be built without some of these libraries, though to get all the benefits of `dhtslib` we recommend installing all of them. For more information you can visit the [htslib](https://github.com/samtools/htslib) repository. 

To intall `htslib` dependencies:
```
Debian / Ubuntu
---------------

sudo apt-get update  # Ensure the package list is up to date
sudo apt-get install autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev

Note: libcurl4-openssl-dev can be used as an alternative to libcurl4-gnutls-dev.

RedHat / CentOS / Amazon Linux
---------------

sudo yum install autoconf automake make gcc perl-Data-Dumper zlib-devel bzip2 bzip2-devel xz-devel curl-devel openssl-devel

Alpine Linux
------------

sudo apk update  # Ensure the package list is up to date
sudo apk add autoconf automake make gcc musl-dev perl bash zlib-dev bzip2-dev xz-dev curl-dev libressl-dev

OpenSUSE
--------

sudo zypper install autoconf automake make gcc perl zlib-devel libbz2-devel xz-devel libcurl-devel libopenssl-devel

MacOS
-----

brew install xz autoconf automake
````
`libdeflate` can also be installed for faster (de)compression but is optional.

#### htslib build and install
You will then need to download htslib (latest release [here](https://github.com/samtools/htslib/releases/latest)). As of now we support versions >=1.10. To install htslib:

```
curl https://github.com/samtools/htslib/releases/download/1.12/htslib-1.12.tar.bz2 -o htslib-1.12.tar.bz2
tar -xjf htslib-1.12.tar.bz2
cd htslib-1.12
./configure
make 
sudo make install
```
htslib is also linked as a submodule in the dhtslib repo and potentially could be used as it reflects the currently supported version (you must git clone recursively).
```
cd htslib
autoreconf -i
./configure
make 
sudo make install
```

### Troubleshooting
#### dub linking issues
By default `htslib` is dynamically linked to `dhtslib` via `dub`.Sometimes it is neccessary to specify the `LIBRARY_PATH` or `LD_LIBRARY_PATH` environment variables in order for `dub` and the D compilers to find your `htslib` installation. This could also apply if your `htslib` isn't a standard install under `/usr/local/lib`.

```
LIBRARY_PATH=/usr/local/lib dub build
LD_LIBRARY_PATH=/usr/local/lib LIBRARY_PATH=/usr/local/lib dub build
```


#### Statically linking to htslib
`libhts.a` needs to be added to your project's source files.
Remember to link to all dynamic libraries configured when htslib was built. This may
include bz2, lzma, zlib, defalate, crypto, pthreads, curl.
Finally, if statically linking, the `-lhts` flag needs to be removed from compilation
by selecting the dub configuration `source-static` as the dub configuration type for dhtslib
within your own project's dub configuration file:

```
"subConfigurations": {
    "dhtslib": "source-static"
},
```
