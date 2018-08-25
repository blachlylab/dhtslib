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

    // Demonstrate that a new copy of the range will be empty
    stderr.writeln("[TEST SUITE] range redux");
    foreach(line; bgzf) {
        writeln(line);
    }

    return 0;
}
