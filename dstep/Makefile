# Adding new rules:
#
# If patching is needed, create two rules.
#
# If no patching is needed, create a single rule.

output/bgzf.d: /tmp/bgzf.d patches/bgzf.diff
	patch /tmp/bgzf.d -o output/bgzf.d < patches/bgzf.diff

/tmp/bgzf.d: include/bgzf.h
	dstep include/bgzf.h -o /tmp/bgzf.d


