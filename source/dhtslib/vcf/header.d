module dhtslib.vcf.header;

import std.string: fromStringz, toStringz;

import htslib.vcf;

/** Wrapper around bcf_hdr_t

    Q: Why do we have VCFHeader, but not SAMHeader?
    A: Most (all)? of the SAM (bam1_t) manipulation functions do not require a ptr to the header,
    whereas almost all of the VCF (bcf1_t) manipulation functions do. Therefore, we track bcf_hdr_t*
    inside each VCFRecord; this is wrapped by VCFHeader for future convenience (for example,
    now we have @property nsamples; may move the tag reader and writer functions here?)

    In order to avoid double free()'ing an instance bcf_hdr_t,
    this wrapper will be the authoritative holder of of bcf_hdr_t ptrs,
    it shall be passed around by reference, and copies are disabled.
*/
struct VCFHeader
{
    /// Pointer to htslib BCF/VCF header struct; will be freed from VCFHeader dtor 
    bcf_hdr_t *hdr;

    // Copies have to be disabled to avoid double free()
    @disable this(this);

    ~this()
    {
        // Deallocate header
        if (this.hdr != null) bcf_hdr_destroy(this.hdr);
    }

    invariant
    {
        assert(this.hdr != null);
    }

    /// List of contigs in the header
    @property string[] sequences()
    {
        import core.stdc.stdlib : free;
        int nseqs;

        /** Creates a list of sequence names. It is up to the caller to free the list (but not the sequence names) */
        //const(char) **bcf_hdr_seqnames(const(bcf_hdr_t) *h, int *nseqs);
        const(char*)*ary = bcf_hdr_seqnames(this.hdr, &nseqs);
        if (!nseqs) return [];

        string[] ret;
        ret.reserve(nseqs);

        for(int i; i < nseqs; i++) {
            ret ~= fromStringz(ary[i]).idup;
        }

        free(cast(void*)ary);
        return ret;        
    }

    /// Number of samples in the header
    pragma(inline, true)
    @property int nsamples() { return bcf_hdr_nsamples(this.hdr); }

}