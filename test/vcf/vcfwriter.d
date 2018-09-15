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
    w.addSample("SAMPLE01");

    w.writeHeader();
    auto vcfhdr = w.getHeader();

    //string[] filters = ["TRIALLELIC", "GOATS"];
    string[] filters = ["PASS", "triallelic", "nonex"];
    VCFRecord r = VCFRecord(vcfhdr, "chr3", 999, "", "C", "T", 40, filters);
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