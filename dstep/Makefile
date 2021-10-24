
DSTEP=dstep --space-after-function-name=false

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

all: $(outputfiles)
.PHONY: all

vpath %.h ../htslib/htslib/
vpath %.d ../source/htslib/

converted/%.d: %.h
	$(DSTEP) $< -o converted/$(notdir $@)

patches/%.diff: converted/%.d %.h %.d
	diff ../source/htslib/$*.d $< | grep -E -e '[0-9]+a[0-9]+' -e '>' > $@; [ $$? -eq 1 ] 

$(outputfiles): output/%.d: %.d patches/%.diff 
	patch $< -o $@ < $(word 2,$^)



