#include <vector>
#include <stdlib.h>
#include <cstdint>
#include <string>

std::string generateLN(const int size);

class LargeNumber {
private:
    std::vector<uint32_t> large_number;
    bool positive;

public:
    LargeNumber(std::string& stringLN);
};