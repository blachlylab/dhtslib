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

all: $(outputfiles) $(patchFiles)
.PHONY: all

# clean out files
clean:
	rm patches/*
	rm output/*
	

patches: $(patchFiles)
.PHONY: patches

vpath %.h include
vpath %.diff patches

# perform dstep conversion of latest
# htslib headers in sub-repo
/tmp/%.d: %.h
	$(DSTEP) -o /tmp/$(notdir $@) $<

# create patch files for patching
# converted headers to current bindings 
# ignores comment changes and whitespace
$(patchFiles): patches/%.diff: %.h
	$(DSTEP) -o /tmp/$(notdir $@) $<
	diff -w /tmp/$(notdir $@) ../source/htslib/$*.d > $@; [ $$? -eq 1 ]

# patch converted to match bindings
$(outputfiles): output/%.d: %.h %.diff
	$(DSTEP) -o /tmp/$(notdir $@) $< 
	patch /tmp/$(notdir $@) -o $@ < patches/$*.diff