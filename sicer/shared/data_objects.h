#ifndef DATAOBJECTS_h
#define DATAOBJECTS_h

#include <string>
 #include <sstream>

using namespace std;

typedef struct BEDRead {
    string chrom;
    unsigned int start;
    unsigned int end;
    string name;
    int score;
    char strand;
public:
    BEDRead() {}
    BEDRead(string chr, unsigned int start, unsigned int end, string name, int scr, char strd):
        chrom(chr), start(start), end(end), name(name), score(scr), strand(strd) {}
    // toString is primarily for debugging
    string toString()
    {   
        ostringstream ostr;
        ostr << chrom << "\t" << start << "\t" << end << "\t" << name << "\t" << score << "\t" << strand;
        return ostr.str();
    }
} BEDRead;

typedef struct Window {
    string chrom;
    unsigned int start;
    unsigned int end;
    unsigned int count;  // Number of reads in the window
public:
    Window() {}
    Window(string chr, unsigned int start, unsigned int end, unsigned int count):
        chrom(chr), start(start), end(end), count(count) {}
    string toString()
    {
        ostringstream ostr;
        ostr << chrom << "\t" << start << "\t" << end << "\t" << count;
        return ostr.str();
    }
} Window;

typedef struct Island {
    string chrom;
    unsigned int start;
    unsigned int end;
    double score;
    // Below fields might be set optionally later
    unsigned int obs_count;
    unsigned int control_count;
    double pvalue;
    double fc;
    double alpha_stat;
public:
    Island() {}
    Island(string chr, unsigned int start, unsigned int end, double score):
        chrom(chr), start(start), end(end), score(score), obs_count(0), control_count(0),
        pvalue(0.0), fc(0.0), alpha_stat(0.0) {}
    string toString()
    {
        ostringstream ostr;
        ostr << chrom << "\t" << start << "\t" << end << "\t" << score;
        return ostr.str();
    }
} Island;

typedef struct DiffExprIsland {
    string chrom;
    unsigned int start;
    unsigned int end;
    unsigned int count_A;
    double norm_count_A;
    unsigned int count_B;
    double norm_count_B;
    double fc_A_vs_B;
    double pvalue_A_vs_B; 
    double fdr_A_vs_B;
    double fc_B_vs_A;
    double pvalue_B_vs_A;
    double fdr_B_vs_A;

public:
    DiffExprIsland() {}
    DiffExprIsland(string chr, unsigned int start, unsigned int end):
        chrom(chr), start(start), end(end), count_A(0), count_B(0) {}
    DiffExprIsland(Island island): 
        chrom(island.chrom), start(island.start), end(island.end), count_A(0), count_B(0) {}
    string toString()
    {
        ostringstream ostr;
        ostr << chrom << "\t" << start << "\t" << end;
        return ostr.str();
    }
} DiffExprIsland;


#endif