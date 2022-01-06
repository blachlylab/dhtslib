module dhtslib.file.htslibfile;

import std.stdio;

import dhtslib.memory;
import htslib.hts;

struct HtslibFile
{
    HtsFile fp;
    string fn;
    File f;

}
