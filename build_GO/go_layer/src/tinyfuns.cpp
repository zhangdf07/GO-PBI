#include "tinyfuns.h"
#include <boost/random.hpp>
#include <cstdarg>
#include <cstdlib>
#include <cstdio>
#include <ctime>

using namespace nkchem;
using namespace boost;
using namespace std;

// One can define the random number generator by modification of tinyrand.
typedef lagged_fibonacci19937 TinyRand;
TinyRand rng(time(NULL));
uniform_01<TinyRand&> u01(rng);  

void nkchem::srand(void)
{
    std::srand(time(NULL));   
}

double nkchem::rand01(void)
{
    return u01();
}

void nkchem::ErrorExit(const char* fmt, ...)
{
    const int BufferSize = 2048;    
    char buffer[BufferSize];
    va_list args;
    va_start(args, fmt);
    vsprintf(buffer, fmt, args);
    va_end(args);
    printf("Error: %s\n", buffer);
    exit(1);
}
