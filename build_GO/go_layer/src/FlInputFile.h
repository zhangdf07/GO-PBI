#ifndef    __FLINPUTFILE_H__
#define    __FLINPUTFILE_H__

#include <string>
#include <vector>

namespace nkchem {

using namespace std;

class FlInputFile {
public:
    enum Style { Random, Assigned };
    const char* RandomChar = "R";
    const char* EvenlyChar = "E";

public:
    FlInputFile(void);
    FlInputFile(const string& fn);
    ~FlInputFile(void);

    void Export(const string& fn) const;

    int nx;
    int ny;
    double rCC;
    double z;
    string resname;
    
    int num_edge_Oe;       vector<int> edge_Oe_indices;
    int num_internal_Oe;   vector<int> internal_Oe_indices;
    int num_edge_OH;       vector<int> edge_OH_indices;
    int num_internal_OH;   vector<int> internal_OH_indices;
    int num_edge_COOH;     vector<int> edge_COOH_indices;
    int num_internal_COOH; vector<int> internal_COOH_indices;    

private:
    void ParseFunctionlization(const string& tag, int dup, const string& str, int& num, vector<int>& indices) const;
    string triminput(const string& str) const;
};

}

#endif //  __FLINPUTFILE_H__
