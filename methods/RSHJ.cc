#include "RSHJ.h"

int main() {
    RSHJ rshj(3, 1.0, "data.bin");
    rshj.preprocess();
    rshj.GSC();
    rshj.BMC();
    rshj.filter();
    return 0;
}