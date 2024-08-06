#ifndef    __FlBUILDER_H__
#define    __FlBUILDER_H__

#include "GrFlake.h"
#include "FlInputFile.h"
#include <vector>
#include <string>

namespace nkchem {

using namespace std;

class FlBuilder {
public:
    FlBuilder(void);
    ~FlBuilder(void);

    static void Build(FlInputFile& inp, GrFlake& flake);
    static void AddGroups(const string& group_name, bool is_edge, GrFlake& flake, int num_groups, vector<int>& indices);
    static void CloseWithH(GrFlake& flake);

private:
    static void AddGroupsCore(const string& group_name, bool is_edge, GrFlake& flake, int num_groups, vector<int>& indices);
};

}

#endif //  __FlBUILDER_H__
