
DSTEP=dstep --space-after-function-name=false

output/bgzf.d: /tmp/bgzf.d patches/bgzf.diff
	patch /tmp/bgzf.d -o output/bgzf.d < patches/bgzf.diff

/tmp/bgzf.d: include/bgzf.h
	dstep include/bgzf.h -o /tmp/bgzf.d

output/cram.d: /tmp/cram.d patches/cram.diff
	$(DSTEP) -o /tmp/cram.d include/cram.h
	patch /tmp/cram.d -o output/cram.d < patches/cram.diff
