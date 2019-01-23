module dhtslib.tagvalue;

import std.stdio;
import std.string:fromStringz;
import dhtslib.htslib.sam:bam_aux_get,bam1_t;

struct TagValue{
    ubyte* data;
    this(bam1_t * b,char[2] tag){
        data=bam_aux_get(b,tag);
    }
    string toString(){
        if(data !is null &&cast(char)data[0]=='Z'){
            return fromStringz(cast(char*)&data[1]).idup;
        }
        return "";
    }
}
unittest{
    import dhtslib.sam;
    import dhtslib.htslib.hts_log;
    hts_set_log_level(htsLogLevel.HTS_LOG_TRACE);
    hts_log_info(__FUNCTION__,"Testing tagvalue");
    hts_log_info(__FUNCTION__,"Loading test file");
    auto bam=SAMFile("htslib/test/range.bam",0);
    auto readrange=bam["CHROMOSOME_I",914];
    hts_log_info(__FUNCTION__,"Getting read 1");
    auto read=readrange.front();
    writeln(read.queryName);
    hts_log_info(__FUNCTION__,"Tag XA:"~TagValue(read.b,['X','A']).toString());
    assert(TagValue(read.b,['X','A']).toString()=="CHROMOSOME_X,-295167,78M1D22M,4;");
}
