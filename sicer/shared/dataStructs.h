#ifndef DATASTRUCTS_H
#define DATASTRUCTS_H
#include <string>

using namespace std;

typedef struct BEDRead {
    string chrom;
    int start;
    int end;
    string name;
    int score;
    char strand;
public:
    BEDRead():
        chrom(""), start(0), end(0), name(""), score(0), strand(""[0]) {}
    BEDRead(string chr, int start, int end, string name, int scr, char strd):
        chrom(chr), start(start), end(end), name(name), score(scr), strand(strd) {}
    string toString() 
    {
        return chrom + "\t" + to_string(start) + "\t" + to_string(end) 
                + "\t" + name + "\t" + to_string(score) + "\t" + strand;
    }
} BEDRead;


#endif