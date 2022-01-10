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
        }else static if(is(T == Bcf1)){
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
            assert(this.itr.initialized && this.itr.getRef);
            return (r < 0 || itr.finished) ? true : false;
        }
    }
}
debug(dhtslib_unittest) unittest
{
    import std.path:buildPath,dirName;
    import std.algorithm: count;
    import std.string: fromStringz;
    hts_log_info(__FUNCTION__, "Testing HtslibIterator SAM/BAM reading");
    
    auto fn = buildPath(dirName(dirName(dirName(dirName(__FILE__)))),"htslib","test","range.bam");
    auto f = HtslibFile(fn);
    f.loadHeader;
    auto read = f.readRecord!Bam1();
    assert(fromStringz(bam_get_qname(read)) == "HS18_09653:4:1315:19857:61712");
    
    f.seek(0);
    f.loadHeader;
    assert(f.byRecord!Bam1.count == 112);

    f.seek(0);
    f.loadHeader;
    f.loadHtsIndex;
    assert(f.query!Bam1(0,913,914).count == 1);

    f.seek(0);
    f.loadHeader;
    f.loadHtsIndex;
    assert(f.query!Bam1(0,913,934).count == 2);

    f.seek(0);
    f.loadHeader;
    f.loadHtsIndex;
    assert(f.query!Bam1(0,913,934).count == 2);
    assert(f.query!Bam1("CHROMOSOME_I:914-935").count == 2);
}

debug(dhtslib_unittest) unittest
{
    import std.algorithm: count;
    import std.path:buildPath,dirName;
    hts_log_info(__FUNCTION__, "Testing HtslibIterator VCF/BCF reading");
    auto fn = buildPath(dirName(dirName(dirName(dirName(__FILE__)))),"htslib","test","tabix","vcf_file.vcf.gz");
    auto f = HtslibFile(fn);
    f.loadHeader;
    auto read = f.readRecord!Bcf1();
    assert(read.pos == 3000149);

    import std.stdio;
    f.seek(0);
    f.loadHeader;
    assert(f.byRecord!Bcf1.count == 14);

    f.loadTabixIndex;
    
    f.seek(0);
    f.loadHeader;
    assert(f.queryTabix!Bcf1(0, 3000149, 3000151).count == 2);

    f.seek(0);
    f.loadHeader;
    assert(f.queryTabix!Bcf1("1:3000150-3000151").count == 2);
}