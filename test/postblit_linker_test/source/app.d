module test.samreader;

import std.stdio;
import std.string;
import std.array;

import dhtslib;
import htslib;

int main()
{
    hts_set_log_level(htsLogLevel.HTS_LOG_INFO);
    auto bam = SAMFile("../../htslib/test/range.bam");
    

    auto samRecs = bam.all_records.array;
    assert(samRecs.length == 112);   // confirmed by samtools flagstat

    auto vcf = VCFReader("../../htslib/test/tabix/vcf_file.vcf");
    auto vcfRecs = vcf.array;
    assert(vcfRecs.length == 14);
    
    // auto vcf = VCFReader("../../htslib/test/tabix/vcf_file.vcf");
    // auto vcfRecs = vcf.array;
    // assert(y.length == 14);
    return 0;
}
