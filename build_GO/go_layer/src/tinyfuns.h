#ifndef    __TINYFUNS_H__
#define    __TINYFUNS_H__

#include <cstdarg>

namespace nkchem {

using namespace std;

void ErrorExit(const char* fmt, ...);
void srand(void);
double rand01(void);

}

#endif //  __TINYFUNS_H__