#ifndef    __BUILDER_H__
#define    __BUILDER_H__

#include "Graphene.h"
#include "GrInputFile.h"
#include <vector>
#include <string>

namespace nkchem {

using namespace std;

class Builder {
public:
    Builder(void);
    ~Builder(void);

    static void PrintHelp(void);
    static void AddGroups(Graphene& gr, int num, int num_up, GrInputFile::Style style, vector<int>& indices, const string& group);

private:

};

}

#endif //  __BUILDER_H__
