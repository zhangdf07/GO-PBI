#ifndef    __OPTIMIZER_H__
#define    __OPTIMIZER_H__

#include <string>
#include <vector>

namespace nkchem {

using namespace std;

class Optimizer {
public:
    Optimizer(void);
    ~Optimizer(void);

    void Do(const string& fn);

private:
    // Fragment.
    int nx;
    int ny;
    double rCC;
    double z;
    string resname;
    // Decorations.
    int num_Oe;
    int num_OH;
    int num_COOH;
    // Results.
    string result_folder;
    int max_step;
    string cmd;

    // Parse the input.
    void ParseInput(const string& fn);
    string triminput(const string& str) const;
    // Mkdir.
    bool makedir(const string& str) const;

};

}

#endif //  __OPTIMIZER_H__
