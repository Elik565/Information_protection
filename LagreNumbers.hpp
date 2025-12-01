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
    bool positive = true;

public:
    // Конструкторы
    LargeNumber() {};
    LargeNumber(std::string& stringLN);  // парсинг строки
    LargeNumber(int64_t LN);  // из int64_t

    std::string toString();
};

struct LucasParams {
    LargeNumber P;
    LargeNumber Q;
    LargeNumber D;
};

class LNMath {
public:
    static LargeNumber sum(const LargeNumber& a, const LargeNumber& b);
    static LargeNumber sub(const LargeNumber& a, const LargeNumber& b);
    static LargeNumber mult(const LargeNumber& a, const LargeNumber& b);
    static LargeNumber div(const LargeNumber& a, const LargeNumber& b);
    static LargeNumber pow(const LargeNumber& a, const LargeNumber& b);
    static LargeNumber sqrt(const LargeNumber&a);
    static LargeNumber gcd(const LargeNumber& a, const LargeNumber& b);  // наибольший общий делитель
    static LargeNumber lcm(const LargeNumber& a, const LargeNumber& b);  // наименьшее общее кратное

    static bool isPrimeStd(const LargeNumber& a);
    static bool sieveEratosthenes(const LargeNumber& a);
    static bool sieveAtkin(const LargeNumber& a);
    static bool LucasLehmer(const LargeNumber& a, uint64_t p);

    static bool MillerRabin(const LargeNumber&a, int iterations = 1);
    static int jacobi(LargeNumber a, const LargeNumber& n);
    static LucasParams findLucasParameters(const LargeNumber& a);
    static LargeNumber computeLucasU(const LargeNumber& n, const LargeNumber& P, 
                                const LargeNumber& Q, const LargeNumber& D);
    static bool StrongLucasTest(const LargeNumber& a);
    static bool BailliePSW(const LargeNumber& n);

private:
    static int compareLN(const LargeNumber& a, const LargeNumber& b);
    static LargeNumber absSum(const LargeNumber& a, const LargeNumber& b);
    static LargeNumber absSub(const LargeNumber& a, const LargeNumber& b);
    static LargeNumber mod(const LargeNumber& a, const LargeNumber& b);
};