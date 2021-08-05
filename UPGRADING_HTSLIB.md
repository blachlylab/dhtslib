Instructions for upgrading dhtslib's API based on changes in htslib
===================================================================

```
cd htslib
git pull

# insert new version here 
git checkout 1.13
git submodule update

# insert old version here
git --no-pager diff 1.12 -- '*.h' > htslib.1.13.diff
```
Look over the output diff and pull over any changes in header files under the htslib folder that we have a corresponding D file for. Pull over any and all changes including documentation and license changes

For new header files, you will need to use dstep to generate bindings. New bindings must be manually inspected and any functions designated as `static inline` will have to be ported by hand from the C header.
