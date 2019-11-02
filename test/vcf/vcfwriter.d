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

    // These should be equivalent: raw htslib call, templated fn
    bcf_hdr_append(w.vcfhdr.hdr, "##FILTER=<ID=triallelic,Description=\"Triallelic site\">");
    w.addTag!"FILTER"("noisy", "Noisy region");

    w.addTag!"INFO"("NS", "1", "Integer", "Number of Samples With Data");
    w.addTag!"INFO"("XFL", "1", "Float", "Floating point number(s)");
    w.addTag!"INFO"("XF", "1", "Flag", "Bool something something mumble");

    // These should be equivalent: raw htslib call, templated fn with string 2nd param, templated fn with int 2nd param
    bcf_hdr_append(w.vcfhdr.hdr, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">");
    w.addTag!"FORMAT"("DDP", "1", "Integer", "De-Duplicated Read Depth (Diamond Dallas Page)");
    w.addTag!"FORMAT"("XDP", 1, "Integer", "X depth");  // pass integer as second param instead of string

    // Test vector valued tags
    w.addTag!"FORMAT"("XXX", 2, "Integer", "Array test");

    w.addSample("Sample1");
    w.addSample("Sample2");

    w.writeHeader();
    auto vcfhdr = w.getHeader();

    // for genotype/format tag value arrays
    const int x = 100;
    const int y = 200;
    const int z = 300;
    const int zz= 400;

    //string[] filters = ["TRIALLELIC", "GOATS"];
    string[] filters = ["PASS", "triallelic", "nonex"];
    VCFRecord r = new VCFRecord(vcfhdr, "chr3", 999, "rs321", "C", "T", 40, filters);
    r.addInfo("NS", 1);
    r.addInfo("XS", "Hello");
    r.addInfo("XFL", 2.1);
    r.addFormat("DP", [x,x] );
    w.writeRecord(r);

    r = new VCFRecord(vcfhdr, "chr3", 1001, "", "", "", 30, "");
    r.setAlleles("A", "G");
    r.addInfo("XF", true);
    r.addFormat("DP", [x,x] );
    w.writeRecord(r);

    r = new VCFRecord(vcfhdr, "chr3", 1002, "", "", "", 30, "PASS");
    r.setAlleles("A", "T", "TCGA");
    r.addFormat("DP", [x,x] );
    r.addFormat("XDP", [x,y]);
    w.writeRecord(r);

    r = new VCFRecord(vcfhdr, "chr3", 1003, "", "A", "G", 30, "PASS");
    r.addFormat("XXX", [x, y, z, zz]);
    w.writeRecord(r);

    r = new VCFRecord(vcfhdr, "chr3", 1004, "", "G", "A", 20, "PASS");
    r.add!"FORMAT"("DP", [x,y]);
    //r.addValue!("FORMAT", int)("DP", [x,y]);
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