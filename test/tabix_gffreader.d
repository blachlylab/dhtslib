module test.tabix_gffreader;

import std.stdio;

import dhtslib.tabix;

int main()
{
    writeln("tabix_gffreader");

    TabixIndexedFile tf = TabixIndexedFile("/tmp/gencode.v28.basic.annotation.sorted.gff3.gz");

    //auto r = tf.region("chr1:1-14409");
    auto r = tf.region("chr1:12000-13000");

    foreach(line ; r) {
        writeln(line[0 .. 80]);
    }

    writeln(r);
/*
    writeln("Writing Range r as array:");
    writeln(r);
    writeln("Repeating again with spent r");
    writeln(r);
    writeln("Repeating again (again) with spent r");
    writeln(r);

    r = tf.region("chr1:1-14409");
    writeln("Writing Range r as rows:");
    foreach(line; r) {
        if (line.length > 80) { writeln(line[0 .. 80]); }
    }
*/
    return 0;
}
