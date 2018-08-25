module test.bgzfreader;

import std.stdio;

import dhtslib.bgzf;

int main(string[] args)
{
    stderr.writeln("[TEST SUITE] BGZF reader");

    auto bgzf = BGZFile(args[1]);

    foreach(line; bgzf) {
        writeln(line);
    }

    stderr.writeln("[TEST SUITE] All done.");

/+
    stderr.writeln("range redux");
    foreach(line; bgzf) {
        writeln(line);
    }
+/
    return 0;
}
