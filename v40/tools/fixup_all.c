#define SFA_IMPLEMENTATION
#include "../../sfa/format/sfa.h"
#include <stdio.h>
int main(int argc, char **argv) {
    for (int i = 1; i < argc; i++) {
        printf("Fixing: %s ... ", argv[i]);
        int ret = sfa_fixup_index(argv[i]);
        printf("%s\n", ret == 0 ? "OK" : "FAILED");
    }
    return 0;
}
