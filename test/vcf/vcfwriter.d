module test.vcfwriter;

import std.stdio;

import dhtslib.vcf;
import dhtslib.htslib.vcf;

int main()
{
    writeln("dhtslib ⚡ VCFWriter");

    VCFWriter w = VCFWriter("output.vcf");

    w.addHeaderLineRaw("##source=dhtslib-vcfwriterV0.4");
    w.addHeaderLineKV("phasing", "none");
    w.addHeaderLineKV("contig", "<ID=chr3,length=999999,assembly=hg19>");
    //w.addSample("SAMPLE01");

    bcf_hdr_append(w.vcfhdr.hdr, "##FILTER=<ID=triallelic,Description=\"Triallelic site\">");

    w.addInfoTag("NS", "1", "Integer", "Number of Samples With Data");
    w.addInfoTag("XFL", "1", "Float", "Floating point number(s)");
    w.addInfoTag("XF", "1", "Flag", "Bool something something mumble");

    bcf_hdr_append(w.vcfhdr.hdr, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">");

    w.writeHeader();
    auto vcfhdr = w.getHeader();

    //string[] filters = ["TRIALLELIC", "GOATS"];
    string[] filters = ["PASS", "triallelic", "nonex"];
    VCFRecord r = VCFRecord(vcfhdr, "chr3", 999, "rs321", "C", "T", 40, filters);
    r.addInfo("NS", 1);
    r.addInfo("XS", "Hello");
    r.addInfo("XFL", 2.1);
    r.addInfo("XF", true);
    w.writeRecord(r);

/+
    filters = [];
    w.addRecord("chr3", 999, "", "C,T", 40, filters);
    w.addRecord("chr3", 999999, ".", "A,T", 29.5, filters);
    w.addRecord("chr3", 1000000, ".", "A,AC", 111.222, filters);
    w.addRecord("chr3", 440, ".", "CG,C", 30, filters);

    //w.addRecord("chr1", 999, ".", "A,G", 30, filters);

    writeln(*w.rows[0]);

    bcf1_t *b = new bcf1_t;
    b.pos = 111;

    w.addRecord(b);

    w.writeFile();
+/


    return 0;
}