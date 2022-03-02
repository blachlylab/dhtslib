SHELL=bash

# commandline from features I was working on
# adding to dstep
# DSTEP=../../dstep/bin/dstep \
	--package htslib \
	--space-after-function-name=false \
	--roots ../htslib/htslib \
	--global-attribute @system \
	--global-attribute nothrow \
	--global-attribute @nogc

DSTEP_PATH=dstep
DSTEP=${DSTEP_PATH} --package htslib --space-after-function-name=false 
outputfiles = \
	output/bgzf.d \
	output/cram.d \
	output/faidx.d \
	output/hfile.d \
	output/hts.d \
	output/hts_endian.d \
	output/hts_expr.d \
	output/hts_log.d \
	output/kbitset.d \
	output/kfunc.d \
	output/knetfile.d \
	output/kroundup.d \
	output/kstring.d \
	output/regidx.d \
	output/sam.d \
	output/synced_bcf_reader.d \
	output/tbx.d \
	output/thread_pool.d \
	output/vcf.d \
	output/vcf_sweep.d \
	output/vcfutils.d

patchFiles = \
	patches/bgzf.diff \
	patches/cram.diff \
	patches/faidx.diff \
	patches/hfile.diff \
	patches/hts.diff \
	patches/hts_endian.diff \
	patches/hts_expr.diff \
	patches/hts_log.diff \
	patches/kbitset.diff \
	patches/kfunc.diff \
	patches/knetfile.diff \
	patches/kroundup.diff \
	patches/kstring.diff \
	patches/regidx.diff \
	patches/sam.diff \
	patches/synced_bcf_reader.diff \
	patches/tbx.diff \
	patches/thread_pool.diff \
	patches/vcf.diff \
	patches/vcf_sweep.diff \
	patches/vcfutils.diff

changePatchFiles = \
	d_changes/bgzf.diff \
	d_changes/cram.diff \
	d_changes/faidx.diff \
	d_changes/hfile.diff \
	d_changes/hts.diff \
	d_changes/hts_endian.diff \
	d_changes/hts_expr.diff \
	d_changes/hts_log.diff \
	d_changes/kbitset.diff \
	d_changes/kfunc.diff \
	d_changes/knetfile.diff \
	d_changes/kroundup.diff \
	d_changes/kstring.diff \
	d_changes/regidx.diff \
	d_changes/sam.diff \
	d_changes/synced_bcf_reader.diff \
	d_changes/tbx.diff \
	d_changes/thread_pool.diff \
	d_changes/vcf.diff \
	d_changes/vcf_sweep.diff \
	d_changes/vcfutils.diff

headerPatchFiles = \
	htslib_changes/bgzf.diff \
	htslib_changes/cram.diff \
	htslib_changes/faidx.diff \
	htslib_changes/hfile.diff \
	htslib_changes/hts.diff \
	htslib_changes/hts_endian.diff \
	htslib_changes/hts_expr.diff \
	htslib_changes/hts_log.diff \
	htslib_changes/kbitset.diff \
	htslib_changes/kfunc.diff \
	htslib_changes/knetfile.diff \
	htslib_changes/kroundup.diff \
	htslib_changes/kstring.diff \
	htslib_changes/regidx.diff \
	htslib_changes/sam.diff \
	htslib_changes/synced_bcf_reader.diff \
	htslib_changes/tbx.diff \
	htslib_changes/thread_pool.diff \
	htslib_changes/vcf.diff \
	htslib_changes/vcf_sweep.diff \
	htslib_changes/vcfutils.diff

all: $(outputfiles) $(patchFiles) $(headerPatchFiles) $(changePatchFiles)
.PHONY: all

patches: $(patchFiles)
.PHONY: patches

changes: $(changePatchFiles) $(headerPatchFiles) 
.PHONY: changes

# clean out files
clean:
	rm patches/*
	rm output/*
	rm /tmp/dstep/*

apply:
	cp output/* ../source/htslib/
	cp include/* ../source/htslib/include

vpath %.h include
vpath %.diff patches

# perform dstep conversion of latest
# htslib headers in sub-repo
/tmp/dstep/%.d: %.h
	$(DSTEP) -o /tmp/dstep/$(notdir $@) $<

# create patch files for patching
# converted headers to current bindings 
$(patchFiles): patches/%.diff: /tmp/dstep/%.d %.h
	diff $< ../source/htslib/$*.d > $@; [ $$? -eq 1 ]


# create patch files for viewing
# significant changes by ignoring
# changes in comment style and whitespace
$(changePatchFiles): d_changes/%.diff: /tmp/dstep/%.d %.h
	diff -w \
	<(sed -E 's;///.*;;g; s;//.*;;g; s;/\*\*;/*;g;' $<) \
	<(sed -E 's;///.*;;g; s;//.*;;g; s;/\*\*;/*;g;' ../source/htslib/$*.d) \
	> $@; [ $$? -eq 1 ]

# create patch files for viewing
# significant changes by ignoring
# changes in comment style and whitespace
$(headerPatchFiles): htslib_changes/%.diff: %.h
	mkdir -p htslib_changes
	diff $< ../source/htslib/include/$*.h > $@; [ $$? -eq 1 ]


# patch converted to match bindings
$(outputfiles): output/%.d: /tmp/dstep/%.d %.h %.diff
	patch $< -o $@ < patches/$*.diff