module dhtslib.file.iterator;

import core.stdc.stdlib;
import core.stdc.string;

import dhtslib.memory;
import dhtslib.file.file;
import dhtslib.util;
import htslib;

/** 
 * HtslibIterator is an abstraction for htslib's hts_itr_t using dhtslib.memory for reference counting.
 * HtslibIterator can be used to iterate VCF/BCF/SAM/BAM/Text files using a BAI/CSI/TBX index
 * or by simply iterating the file.
 * 
 * This struct adapts htslib's iterators into ForwardRange. It is paired with 
 * and relies on HtslibFile. 
 */
struct HtslibIterator(T)
if(is(T == Bam1) || is(T == Bcf1) || is(T == Kstring))
{
    HtslibFile f;               /// HtslibFile
    HtsItr itr;                 /// refcounted hts_itr_t
    T rec;                      /// refcounted bam1_t, bcf1_t, or kstring_t
    private Kstring line;
    private bool is_itr;        /// Using an Itr or just calling *_read functions
    private bool initialized;   /// Is the range initialized
    /// htslib read function return value
    /// -1 on eof, < -1 on err
    /// must be pointer to keep range updated when blitted
    private int * r; 

    /// ctor using only HtslibFile
    /// without an iterator
    this(HtslibFile f)
    {
        this.r = new int(0);
        this.f = f;
        static if(is(T == Bam1)){
            rec = Bam1(bam_init1);
        }else static if(is(T == Bcf1)){
            rec = Bcf1(bcf_init);
        }else static if(is(T == Kstring)){
            rec = Kstring(initKstring);
            ks_initialize(rec);
        }
        this.line = Kstring(initKstring);
        ks_initialize(this.line);
        this.empty;
    }

    /// ctor using an HtslibFile and an iterator
    this(HtslibFile f, HtsItr itr)
    {
        this.is_itr = true;
        this.itr = itr;
        this(f);
    }

    /// Duplicate the iterator
    /// must fully duplicate hts_itr_t
    HtslibIterator dup()
    {

        HtslibIterator newItr;
        // dup file
        newItr.f = this.f.dup;

        // set private values
        newItr.r = new int(*this.r);
        newItr.initialized = this.initialized;
        newItr.is_itr = this.is_itr;

        // duplicate current record
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
        auto ks2 = Kstring(initKstring);
        ks_initialize(ks2);
        kputs(ks_c_str(this.line), ks2);
        newItr.line = ks2;

        // if itr we need to deep copy it
        if(this.is_itr){
            // init
            auto newHtsItr = cast(hts_itr_t *) malloc(hts_itr_t.sizeof);

            // bulk copy data
            *newHtsItr = *itr;

            // initialize and copy region list
            // if it is not null
            if(newHtsItr.reg_list != null){

                // initialize the region lists
                newHtsItr.reg_list = cast(hts_reglist_t *) malloc(itr.n_reg * hts_reglist_t.sizeof);
                newHtsItr.reg_list[0 .. newHtsItr.n_reg] = itr.reg_list[0 .. itr.n_reg];

                // for each list
                for(auto i=0; i < newHtsItr.n_reg; i++)
                {
                    // copy all intervals
                    newHtsItr.reg_list[i].intervals = cast(hts_pair_pos_t *) malloc(itr.reg_list[i].count * hts_pair_pos_t.sizeof);
                    newHtsItr.reg_list[i].intervals[0 .. newHtsItr.reg_list[i].count] = itr.reg_list[i].intervals[0 .. itr.reg_list[i].count];

                    // copy region c string
                    // if not null
                    if(itr.reg_list[i].reg){
                        newHtsItr.reg_list[i].reg = cast(char *) malloc(strlen(itr.reg_list[i].reg) + 1);
                        memcpy(cast(void*)(newHtsItr.reg_list[i].reg),itr.reg_list[i].reg, strlen(itr.reg_list[i].reg) + 1);
                    }
                }
            }

            // initialize and copy bins list
            newHtsItr.bins.a = cast(int *) malloc(itr.bins.m * int.sizeof);
            assert(newHtsItr.bins.m >= newHtsItr.bins.n);
            newHtsItr.bins.a[0 .. newHtsItr.bins.m] = itr.bins.a[0 .. itr.bins.m];

            // initialize and copy off list
            newHtsItr.off = cast(hts_pair64_max_t *) malloc(itr.n_off * hts_pair64_max_t.sizeof);
            newHtsItr.off[0 .. newHtsItr.n_off] = itr.off[0 .. itr.n_off];
            
            // set new itr
            newItr.itr = HtsItr(newHtsItr);
        }else{
            newItr.itr = this.itr;
        }
        
        return newItr;
    }

    /// Get front of the iterator, returns  Bam1, Bcf1, or Kstring
    /// Backing bam1_t, bcf1_t, kstring_t is re-used 
    /// If you keep the result around it should be duplicated
    T front()
    {
        return rec;
    }

    /// popFront to move range forward
    /// is destructive
    void popFront()
    {
        if(!is_itr){
            static if(is(T == Bam1)){
                assert(this.f.bamHdr.initialized && this.f.bamHdr.getRef);
                *this.r = sam_read1(this.f.fp, this.f.bamHdr, rec);
            }
            else static if(is(T == Bcf1)){
                assert(this.f.bcfHdr.initialized && this.f.bcfHdr.getRef);
                *this.r = bcf_read(this.f.fp, this.f.bcfHdr, rec);
            }else static if(is(T == Kstring)){
                *this.r = hts_getline(this.f.fp, cast(int)'\n', rec);
            }
        }else{
            if (itr.multi)
                *this.r = hts_itr_multi_next(f.fp, itr, rec.getRef);
            else{
                if(f.idx != null)
                    *this.r = hts_itr_next(f.fp.is_bgzf ? f.fp.fp.bgzf : null, itr, rec.getRef, f.fp);
                else if(f.tbx != null){
                    static if(is(T == Kstring))
                        *this.r = hts_itr_next(f.fp.is_bgzf ? f.fp.fp.bgzf : null, itr, rec.getRef, f.tbx);
                    else {
                        *this.r = hts_itr_next(f.fp.is_bgzf ? f.fp.fp.bgzf : null, itr, this.line.getRef, f.tbx);
                        static if(is(T == Bam1))
                            *this.r = sam_parse1(this.line, this.f.bamHdr, this.rec);
                        else static if(is(T == Bcf1))
                            *this.r = vcf_parse(this.line, this.f.bcfHdr, this.rec);
                        ks_clear(this.line);
                    }
                    
                }
                else{
                    *r = -2;
                    hts_log_error(__FUNCTION__, "Neither tabix nor bai/csi index are loaded");
                }
                    
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
            return *r < 0 ? true : false;
        }else{
            assert(this.itr.initialized && this.itr.getRef);
            return (*r < 0 || itr.finished) ? true : false;
        }
    }
    
    /// Needed to be ForwardRange
    typeof(this) save()
    {
        return this.dup;
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

    f.seek(0);
    f.loadHeader;
    string[] regions = ["CHROMOSOME_I:900-2000","CHROMOSOME_II:900-2000"];
    assert(f.query(regions).count == 33);
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
    assert(f.query!Bcf1(0, 3000149, 3000151).count == 2);

    f.seek(0);
    f.loadHeader;
    assert(f.query!Bcf1("1:3000150-3000151").count == 2);
}

debug(dhtslib_unittest) unittest
{
    import std.algorithm: count;
    import std.path:buildPath,dirName;
    import std.string : fromStringz;
    hts_log_info(__FUNCTION__, "Testing HtslibIterator dup");

    auto fn = buildPath(dirName(dirName(dirName(dirName(__FILE__)))),"htslib","test","range.bam");
    auto f = HtslibFile(fn);
    f.loadHeader;
    f.loadHtsIndex;
    auto read = f.readRecord!Bam1();
    assert(fromStringz(bam_get_qname(read)) == "HS18_09653:4:1315:19857:61712");

    auto range1 = f.query!Bam1(0,900,2000);
    auto range2 = range1.dup;
    import std.stdio;

    assert(range1.count == 14);
    assert(range1.empty == true);
    assert(range2.count == 14);

    f.seek(0);
    f.loadHeader;

    range1 = f.query!Bam1(0,900,2000);
    range1.popFront;
    range2 = range1.dup;
    range2.popFront;

    assert(range1.count == 13);
    assert(range1.empty == true);
    assert(range2.count == 12);

    f.seek(0);
    f.loadHeader;

    range1 = f.byRecord!Bam1();
    range1.popFront;
    range2 = range1.dup;
    range2.popFront;

    assert(range1.count == 111);
    assert(range1.count == 0);
    assert(range1.empty == true);
    assert(range2.count == 110);

    f.seek(0);
    f.loadHeader;
    string[] regions = ["CHROMOSOME_I:900-2000","CHROMOSOME_II:900-2000"];
    range1 = f.query(regions);
    range2 = range1.dup;
    range2.popFront;
    assert(range1.count == 33);
    assert(range2.count == 32);

    f.seek(0);
    f.loadHeader;
    auto f2 =  f.dup;
    read = f.readRecord!Bam1();
    assert(fromStringz(bam_get_qname(read)) == "HS18_09653:4:1315:19857:61712");
    read = f2.readRecord!Bam1();
    assert(fromStringz(bam_get_qname(read)) == "HS18_09653:4:1315:19857:61712");


    fn = buildPath(dirName(dirName(dirName(dirName(__FILE__)))),"htslib","test","tabix","vcf_file.vcf");
    f = HtslibFile(fn);
    f.loadHeader;
    auto vrange1 = f.byRecord!Bcf1();
    vrange1.popFront;
    auto vrange2 = vrange1.dup;
    vrange2.popFront;

    assert(vrange1.count == 13);
    assert(vrange2.count == 12);

    fn = buildPath(dirName(dirName(dirName(dirName(__FILE__)))),"htslib","test","tabix","bed_file.bed");
    f = HtslibFile(fn);
    f.loadHeader;
    auto trange1 = f.byRecord!Kstring();
    trange1.popFront;
    auto trange2 = trange1.dup;
    trange2.popFront;

    assert(trange1.count == 14);
    assert(trange2.count == 13);

    fn = buildPath(dirName(dirName(dirName(dirName(__FILE__)))),"htslib","test","tabix","bed_file.bed.gz");
    f = HtslibFile(fn);
    f.loadHeader;
    trange1 = f.byRecord!Kstring();
    trange1.popFront;
    trange2 = trange1.dup;
    trange2.popFront;

    assert(trange1.count == 14);
    assert(trange2.count == 13);

}