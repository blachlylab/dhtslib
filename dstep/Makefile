
DSTEP=dstep --space-after-function-name=false

outputfiles=output/faidx.d

all: $(outputfiles)
.PHONY: all

vpath %.h include
vpath %.diff patches

$(outputfiles): output/%.d: %.h %.diff
	$(DSTEP) $< -o /tmp/$(notdir $@)
	patch /tmp/$(notdir $@) -o $@ < patches/$*.diff

output/bgzf.d: /tmp/bgzf.d patches/bgzf.diff
	patch /tmp/bgzf.d -o output/bgzf.d < patches/bgzf.diff

/tmp/bgzf.d: include/bgzf.h
	dstep include/bgzf.h -o /tmp/bgzf.d

output/cram.d: /tmp/cram.d patches/cram.diff
	$(DSTEP) -o /tmp/cram.d include/cram.h
	patch /tmp/cram.d -o output/cram.d < patches/cram.diff


