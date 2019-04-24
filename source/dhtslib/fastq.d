module dhtslib.fastq;
import dhtslib.bgzf;

struct FastqRecord{
    string id;
    string sequence;
    string extra;
    string qscores;
}

struct FastqFile{
    BGZFile f;
    FastqRecord rec;
    this(string fn){
        f=BGZFile(fn);
        popFront;
    }
    FastqRecord front(){
        return rec;
    }
    void popFront(){
        rec.id = f.front;
        f.popFront;
        rec.sequence = f.front;
        f.popFront;
        rec.extra = f.front;
        f.popFront;
        rec.qscores = f.front;
        f.popFront;
    }
    bool empty(){
        return f.empty;
    }

}