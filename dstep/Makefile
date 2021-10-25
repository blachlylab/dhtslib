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


DSTEP=dstep --package htslib --space-after-function-name=false 
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
	changes/bgzf.diff \
	changes/cram.diff \
	changes/faidx.diff \
	changes/hfile.diff \
	changes/hts.diff \
	changes/hts_endian.diff \
	changes/hts_expr.diff \
	changes/hts_log.diff \
	changes/kbitset.diff \
	changes/kfunc.diff \
	changes/knetfile.diff \
	changes/kroundup.diff \
	changes/kstring.diff \
	changes/regidx.diff \
	changes/sam.diff \
	changes/synced_bcf_reader.diff \
	changes/tbx.diff \
	changes/thread_pool.diff \
	changes/vcf.diff \
	changes/vcf_sweep.diff \
	changes/vcfutils.diff

all: $(outputfiles) $(patchFiles)
.PHONY: all

patches: $(patchFiles)
.PHONY: patches

changes: $(changePatchFiles)
.PHONY: changes

# clean out files
clean:
	rm patches/*
	rm output/*

apply:
	cp output/* ../source/htslib/

vpath %.h include
vpath %.diff patches

# perform dstep conversion of latest
# htslib headers in sub-repo
/tmp/%.d: %.h
	$(DSTEP) -o /tmp/$(notdir $@) $<

# create patch files for patching
# converted headers to current bindings 
$(patchFiles): patches/%.diff: /tmp/%.d %.h
	diff $< ../source/htslib/$*.d > $@; [ $$? -eq 1 ]


# create patch files for viewing
# significant changes by ignoring
# changes in comment style and whitespace
$(changePatchFiles): changes/%.diff: /tmp/%.d %.h
	diff -w \
	<(sed -E 's;///.*;;g; s;//.*;;g; s;/\*\*;/*;g;' $<) \
	<(sed -E 's;///.*;;g; s;//.*;;g; s;/\*\*;/*;g;' ../source/htslib/$*.d) \
	> $@; [ $$? -eq 1 ]


# patch converted to match bindings
$(outputfiles): output/%.d: /tmp/%.d %.h %.diff
	patch $< -o $@ < patches/$*.diff