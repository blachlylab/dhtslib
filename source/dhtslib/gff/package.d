/**

GFF3 Record abstraction

This module provides a readable, writeable abstraction of GFF3 records.

For reading and writing from and to bgzip-compressed/tabix-indexed files,
use dhtslib -- https://github.com/blachlylab/dhtslib

Authors: James S Blachly, MD <james.blachly@gmail.com>; Thomas Gregory <charles.gregory@osumc.edu>
License: MIT
Date: 2019-01-28
Standards: http://gmod.org/wiki/GFF3
    https://useast.ensembl.org/info/website/upload/gff3.html
    https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
    http://www.sequenceontology.org/gff3.shtml
*/
module dhtslib.gff;

public import dhtslib.gff.gffrecord;
public import dhtslib.gff.gffreader;

unittest{
    import dhtslib.coordinates;
    auto rec    = GTFRecord("chr1\tHAVANA\tgene\t11869\t14409\t.\t+\t.\tID \"ENSG00000223972.5\" ; gene_id ENSG00000223972.5 ; gene_id ENSG00000223972.5 ; gene_type transcribed_unprocessed_pseudogene ; gene_name DDX11L1 ; level 2 ; havana_gene OTTHUMG00000000961.2"); // @suppress(dscanner.style.long_line)
    auto rec_neg= GTFRecord("chr1\tHAVANA\tgene\t11869\t14409\t.\t-\t.\tID \"ENSG00000223972.5\" ; gene_id ENSG00000223972.5 ; gene_id ENSG00000223972.5 ; gene_type transcribed_unprocessed_pseudogene ; gene_name DDX11L1 ; level 2 ; havana_gene OTTHUMG00000000961.2"); // @suppress(dscanner.style.long_line)

    assert(rec.seqid=="chr1");
    assert(rec.source=="HAVANA");
    assert(rec.type=="gene");
    assert(rec.start==11_869);
    assert(rec.end==14_409);
    assert(rec.score==-1.0);
    assert(rec.strand()=='+');
    assert(rec.phase==-1);
    assert(rec["ID"] == "ENSG00000223972.5");
    assert(rec["gene_id"] == "ENSG00000223972.5");
    assert(rec["gene_type"] == "transcribed_unprocessed_pseudogene");
    assert(rec["gene_name"] == "DDX11L1");
    assert(rec["level"] == "2");
    assert(rec["havana_gene"] == "OTTHUMG00000000961.2");

    assert(rec.length == 2541);
    assert(rec.relativeStart == 1);
    assert(rec.relativeEnd == 2541);

    // Test forward and backward offsets
    assert(rec.coordinateAtOffset(2) == 11_870);
    assert(rec_neg.coordinateAtOffset(2) == 14_408);

    assert(rec.coordinateAtBegin == 11_869);
    assert(rec.coordinateAtEnd   == 14_409);

    assert(rec_neg.coordinateAtBegin == 14_409);
    assert(rec_neg.coordinateAtEnd   == 11_869);

    rec.seqid = "chr2";
    rec.source = "HAVANA1";
    rec.type = "gene1";
    rec.coordinates = OBC(11_870, 14_410);
    rec.score = 1.0;
    rec.strand = '-';
    rec.phase = 1;
    rec["ID"] = "ENSG00000223972.51";
    rec["gene_id"] = "ENSG00000223972.51";
    rec["gene_type"] = "transcribed_unprocessed_pseudogene1";
    rec["gene_name"] = "DDX11L11";
    rec["level"] = "21";
    rec["havana_gene"] = "OTTHUMG00000000961.21";

    assert(rec.seqid=="chr2");
    assert(rec.source=="HAVANA1");
    assert(rec.type=="gene1");
    assert(rec.start==11_870);
    assert(rec.end==14_410);
    assert(rec.score==1.0);
    assert(rec.strand()=='-');
    assert(rec.phase==1);
    assert(rec["ID"] == "\"ENSG00000223972.51\"");
    assert(rec["gene_id"] == "\"ENSG00000223972.51\"");
    assert(rec["gene_type"] == "\"transcribed_unprocessed_pseudogene1\"");
    assert(rec["gene_name"] == "\"DDX11L11\"");
    assert(rec["level"] == "\"21\"");
    assert(rec["havana_gene"] == "\"OTTHUMG00000000961.21\"");

    assert(rec.length == 2541);
    assert(rec.relativeStart == 1);
    assert(rec.relativeEnd == 2541);

    // Test forward and backward offsets
    assert(rec.coordinateAtOffset(2) == 14_409);

    assert(rec.coordinateAtBegin == 14_410);
    assert(rec.coordinateAtEnd   == 11_870);

    // TODO validator
    assert(rec.isValid);
}

unittest{
    import dhtslib.coordinates;
    auto rec    = GFF3Record("chr1\tHAVANA\tgene\t11869\t14409\t.\t+\t.\tID=ENSG00000223972.5;gene_id=ENSG00000223972.5;gene_id=ENSG00000223972.5;gene_type=transcribed_unprocessed_pseudogene;gene_name=DDX11L1;level=2;havana_gene=OTTHUMG00000000961.2"); // @suppress(dscanner.style.long_line)
    auto rec_neg= GFF3Record("chr1\tHAVANA\tgene\t11869\t14409\t.\t-\t.\tID=ENSG00000223972.5;gene_id=ENSG00000223972.5;gene_id=ENSG00000223972.5;gene_type=transcribed_unprocessed_pseudogene;gene_name=DDX11L1;level=2;havana_gene=OTTHUMG00000000961.2"); // @suppress(dscanner.style.long_line)

    assert(rec.seqid=="chr1");
    assert(rec.source=="HAVANA");
    assert(rec.type=="gene");
    assert(rec.start==11_869);
    assert(rec.end==14_409);
    assert(rec.score==-1.0);
    assert(rec.strand()=='+');
    assert(rec.phase==-1);
    assert(rec["ID"] == "ENSG00000223972.5");
    assert(rec["gene_id"] == "ENSG00000223972.5");
    assert(rec["gene_type"] == "transcribed_unprocessed_pseudogene");
    assert(rec["gene_name"] == "DDX11L1");
    assert(rec["level"] == "2");
    assert(rec["havana_gene"] == "OTTHUMG00000000961.2");

    assert(rec.length == 2541);
    assert(rec.relativeStart == 1);
    assert(rec.relativeEnd == 2541);

    // Test forward and backward offsets
    assert(rec.coordinateAtOffset(2) == 11_870);
    assert(rec_neg.coordinateAtOffset(2) == 14_408);

    assert(rec.coordinateAtBegin == 11_869);
    assert(rec.coordinateAtEnd   == 14_409);

    assert(rec_neg.coordinateAtBegin == 14_409);
    assert(rec_neg.coordinateAtEnd   == 11_869);

    rec.seqid = "chr2";
    rec.source = "HAVANA1";
    rec.type = "gene1";
    rec.coordinates = OBC(11_870, 14_410);
    rec.score = 1.0;
    rec.strand = '-';
    rec.phase = 1;
    rec["ID"] = "ENSG00000223972.51";
    rec["gene_id"] = "ENSG00000223972.51";
    rec["gene_type"] = "transcribed_unprocessed_pseudogene1";
    rec["gene_name"] = "DDX11L11";
    rec["level"] = "21";
    rec["havana_gene"] = "OTTHUMG00000000961.21";

    assert(rec.seqid=="chr2");
    assert(rec.source=="HAVANA1");
    assert(rec.type=="gene1");
    assert(rec.start==11_870);
    assert(rec.end==14_410);
    assert(rec.score==1.0);
    assert(rec.strand()=='-');
    assert(rec.phase==1);
    assert(rec["ID"] == "ENSG00000223972.51");
    assert(rec["gene_id"] == "ENSG00000223972.51");
    assert(rec["gene_type"] == "transcribed_unprocessed_pseudogene1");
    assert(rec["gene_name"] == "DDX11L11");
    assert(rec["level"] == "21");
    assert(rec["havana_gene"] == "OTTHUMG00000000961.21");

    assert(rec.length == 2541);
    assert(rec.relativeStart == 1);
    assert(rec.relativeEnd == 2541);

    // Test forward and backward offsets
    assert(rec.coordinateAtOffset(2) == 14_409);

    assert(rec.coordinateAtBegin == 14_410);
    assert(rec.coordinateAtEnd   == 11_870);

    // TODO validator
    assert(rec.isValid);
}

