#ifndef    __GRINPUTFILE_H__
#define    __GRINPUTFILE_H__

#include <string>
#include <vector>

namespace nkchem {

using namespace std;

class GrInputFile {
public:
    enum Style { Random, Evenly, Assigned };
    const char* RandomChar = "R";
    const char* EvenlyChar = "E";

public:
    GrInputFile(const string& fn);
    ~GrInputFile(void);

    void Export(const string& fn) const;


    int nx;
    int ny;
    double rCC;
    double z;
    string resname;

    int num_OH;   int num_OH_up;   Style style_OH;   vector<int> OH_indices;
    int num_Oe;   int num_Oe_up;   Style style_Oe;   vector<int> Oe_indices;
    int num_H;    int num_H_up;    Style style_H;    vector<int> H_indices;
    int num_COOH; int num_COOH_up; Style style_COOH; vector<int> COOH_indices;

private:
    void ParseFunctionlization(const string& tag, int dup, const string& str, int& num, int& num_up, Style& sty, vector<int>& indices) const;
    string triminput(const string& str) const;
};

}

#endif //  __INPUTFILE_H__
