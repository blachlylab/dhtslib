module test.tabix_gffreader;

import std.stdio;

import dhtslib.tabix;

int main()
{
    writeln("tabix_gffreader");

    TabixIndexedFile tf = TabixIndexedFile("/tmp/gencode.v28.basic.annotation.sorted.gff3.gz");

    writeln("header: ", tf.header);

    auto r = tf.region("chr1:1-14409");

    writeln("Writing Range r:");
    foreach(line; r) {
        if (line.length > 80) { writeln(line[0 .. 80]); }
    }

    return 0;
}