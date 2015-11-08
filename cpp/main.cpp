#include <iostream>
#include <fstream>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include "draw.h"
#include "buddha.h"
#include "CImg-1.6.8_pre110615/CImg.h"

int main()
{
    srand (time(NULL));
    std::vector<uint64_t> its;
    its.push_back(100);
    its.push_back(500);
    its.push_back(2000);
    GenerateAndSaveHistogram(1024, 4, its, 10, "buddy.a");
    return 0;
}
