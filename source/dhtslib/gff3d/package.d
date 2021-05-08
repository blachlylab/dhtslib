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
module dhtslib.gff3d;

public import dhtslib.gff3d.gff3record;
public import dhtslib.gff3d.gtfrecord;
