module test.tabix_gffreader;

import std.stdio;

import dhtslib.tabix;

int main(string[] args)
{
    writeln("tabix_gffreader");

    writeln(args[0], args[1]);
    TabixIndexedFile tf = TabixIndexedFile(args[1]);

    //auto r = tf.region("chr1:1-14409");
    auto r = tf.region("chr1:12000-12000");

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
