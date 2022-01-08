module dhtslib.file.iterator;

import core.stdc.stdlib;

import dhtslib.memory;
import dhtslib.file.file;
import dhtslib.util;
import htslib;

struct HtslibIterator(T, bool is_tbx = false)
if(is(T == Bam1) || is(T == Bcf1) || is(T == Kstring))
{
    HtslibFile f;
    HtsItr itr;
    T rec;
    private bool is_itr;
    private bool initialized;
    private int r;
    this(HtslibFile f)
    {
        this.f = f;
        static if(is(T == Bam1)){
            rec = Bam1(bam_init1);
        }else static if(is(T == BcfHdr)){
            rec = Bcf1(bcf_init);
        }else static if(is(T == Kstring)){
            rec = Kstring(initKstring);
            ks_initialize(rec);
        }
        this.empty;
    }

    this(HtslibFile f, HtsItr itr)
    {
        this.is_itr = true;
        this.itr = itr;
        this(f);
    }

    HtslibIterator dup()
    {
        /// init
        auto newHtsItr = cast(hts_itr_t *) malloc(hts_itr_t.sizeof);
        /// bulk copy data
        *newHtsItr = *itr;

        /// initialize and copy region list
        /// if it is not null
        if(newHtsItr.reg_list){
            /// initialize the region lists
            newHtsItr.reg_list = cast(hts_reglist_t *) malloc(itr.n_reg * hts_reglist_t.sizeof);
            newHtsItr.reg_list[0 .. newHtsItr.n_reg] = itr.reg_list[0 .. itr.n_reg];
            /// for each list
            for(auto i=0; i < newHtsItr.n_reg; i++)
            {
                /// copy all intervals
                newHtsItr.reg_list[i].intervals = cast(hts_pair_pos_t *) malloc(itr.reg_list[i].count * hts_pair_pos_t.sizeof);
                newHtsItr.reg_list[i].intervals[0 .. newHtsItr.reg_list[i].count] = itr.reg_list[i].intervals[0 .. itr.reg_list[i].count];
            }
        }

        /// initialize and copy bins list
        newHtsItr.bins.a = cast(int *) malloc(itr.bins.m * int.sizeof);
        assert(newHtsItr.bins.m >= newHtsItr.bins.n);
        newHtsItr.bins.a[0 .. newHtsItr.bins.m] = itr.bins.a[0 .. itr.bins.m];

        /// initialize and copy off list
        newHtsItr.off = cast(hts_pair64_max_t *) malloc(itr.n_off * hts_pair64_max_t.sizeof);
        newHtsItr.off[0 .. newHtsItr.n_off] = itr.off[0 .. itr.n_off];

        /// Create new HtslibIterator
        auto newItr = HtslibIterator(this.f.dup, HtsItr(newHtsItr));

        /// set private values
        newItr.r = this.r;
        newItr.initialized = this.initialized;

        /// duplicate current record
        static if(is(T == Bam1)){
            newItr.rec = Bam1(bam_dup1(this.rec));
        }else static if(is(T == Bcf1)){
            newItr.rec = Bcf1(bcf_dup(this.rec));
        }else static if(is(T == Kstring)){
            auto ks = Kstring(initKstring);
            ks_initialize(ks);
            kputs(this.rec.s, ks);
            newItr.rec = ks;
        }
        return newItr;
    }

    T front()
    {
        return rec;
    }

    void popFront()
    {
        if(!is_itr){
            static if(is(T == Bam1)){
                assert(this.f.bamHdr.initialized && this.f.bamHdr.getRef);
                this.r = sam_read1(this.f.fp, this.f.bamHdr, rec);
            }
            else static if(is(T == Bcf1)){
                assert(this.f.bcfHdr.initialized && this.f.bcfHdr.getRef);
                this.r = bcf_read(this.f.fp, this.f.bcfHdr, rec);
            }else static if(is(T == Kstring)){
                this.r = hts_getline(this.f.fp, cast(int)'\n', rec);
            }
        }else{
            if (itr.multi)
                this.r = hts_itr_multi_next(f.fp, itr, rec.getRef);
            else{
                static if(is_tbx)
                    this.r = hts_itr_next(f.fp.is_bgzf ? f.fp.fp.bgzf : null, itr, rec.getRef, f.tbx);
                else 
                    this.r = hts_itr_next(f.fp.is_bgzf ? f.fp.fp.bgzf : null, itr, rec.getRef, f.fp);
            }
        }
        
    }

    /// InputRange interface
    @property bool empty()
    {
        // TODO, itr.finished shouldn't be used
        if (!this.initialized) {
            this.popFront();
            this.initialized = true;
        }
        if(!is_itr){
            return r < 0 ? true : false;
        }else{
            assert(this.itr !is null);
            return (r < 0 || itr.finished) ? true : false;
        }
    }
}