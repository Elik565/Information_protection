#include <vector>
#include <stdlib.h>
#include <cstdint>
#include <string>
#define BASE 1000000000  // 9 цифр

std::string generateLN(const int size);

class LargeNumber {
    friend class LNMath;
private:
    std::vector<uint32_t> large_number;
    bool positive;

public:
    // Конструкторы
    LargeNumber() {};
    LargeNumber(std::string& stringLN);  // (парсинг строки)

    std::string toString();
};

class LNMath {
public:
    static LargeNumber sum(const LargeNumber& a, const LargeNumber& b);
    static LargeNumber sub(const LargeNumber& a, const LargeNumber& b);
    static LargeNumber mult(const LargeNumber& a, const LargeNumber& b);
    static LargeNumber div(const LargeNumber& a, const LargeNumber& b);
    static LargeNumber pow(const LargeNumber& a, const LargeNumber& b);

private:
    static int compareLN(const LargeNumber& a, const LargeNumber& b);
    static LargeNumber absSum(const LargeNumber& a, const LargeNumber& b);
    static LargeNumber absSub(const LargeNumber& a, const LargeNumber& b);
};