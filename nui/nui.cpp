#include "NexusParse.h"


int main(int argc, char *argv[])
{
    CNexusParse cNexusParse(argv[1], argv[2]);
    cNexusParse.ReadNexusFile();
    cNexusParse.Report();
    return 0;
}

