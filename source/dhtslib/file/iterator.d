module dhtslib.file.Iterator;

import core.stdc.stdlib;

import dhtslib.memory;
import dhtslib.file.file;
import htslib.hts;

struct HtslibIterator(T)
if(is(T == Bam1) || is(T == Bcf1) || is(T == Tbx))
{
    HtslibFile f;
    HtsItr itr;
    T rec;
    private bool initialized;
    private int r;

    this(HtslibFile f, HtsItr itr)
    {
        this.f = f;
        this.itr = itr;
        this.empty;
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
        auto newItr = HtslibIterator(this.f, HtsItr(newHtsItr));
        
        return newItr;
    }

    T front()
    {
        return rec;
    }

    void popFront()
    {
        if (itr.multi)
            this.r = hts_itr_multi_next(f.fp, itr, rec.getRef);
        else
            this.r = hts_itr_next(f.fp.is_bgzf ? f.fp.fp.bgzf : null, itr, rec.getRef, f.fp);
    }

    /// InputRange interface
    @property bool empty()
    {
        // TODO, itr.finished shouldn't be used
        if (!this.initialized) {
            this.popFront();
            this.initialized = true;
        }
        assert(this.itr !is null);
        return (r < 0 || itr.finished) ? true : false;
    }
}