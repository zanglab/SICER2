#ifndef DATAOBJECTS_h
#define DATAOBJECTS_h

#include <string>
#include <vector>
#include <map>

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

typedef struct Window {
    string chrom;
    int start;
    int end;
    int count;  // Number of reads in the window
public:
    Window():
        chrom(""), start(0), end(0), count(0) {}
    Window(string chr, int start, int end, int count):
        chrom(chr), start(start), end(end), count(count) {}
    string toString()
    {
        return chrom + "\t" + to_string(start) + "\t" 
                + to_string(end) + "\t" + to_string(count);
    }
} Window;

typedef struct Island {
    string chrom;
    int start;
    int end;
    double score;
public:
    Island():
        chrom(""), start(0), end(0), score(0.0) {}
    Island(string chr, int start, int end, double score):
        chrom(chr), start(start), end(end), score(score) {}
    string toString()
    {
        return chrom + "\t" + to_string(start) + "\t" 
                + to_string(end) + "\t" + to_string(score);
    }
} Island;

#endif