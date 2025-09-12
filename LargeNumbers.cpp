#include "LagreNumbers.hpp"
#include <random>

std::string generateLN(const int size) {
    std::string LN;
    LN.reserve(size + 1);  // сразу выделяем память

    char sign = (rand() % 2 == 0) ? '-' : '+';
    LN += sign;

    LN += '1' + rand() % 9;  // предотвращаем ведущие нули

    for (int i = 0; i < size; i++) {
        char digit = rand() % 10 + '0';
        LN += digit;
    }

    return LN;
}

LargeNumber::LargeNumber(std::string& stringLN) {
    positive = stringLN[0] != '-';

    stringLN.erase(0, 1);   

    int i = stringLN.size();
    while (i > 0) {
        int begin_block = i - 9;
        if (begin_block < 0) {
            begin_block = 0;
        }

        uint32_t block = 0;
        for (int j = begin_block; j < i; ++j) {
            block = block * 10 + (stringLN[j] - '0');
        }

        large_number.push_back(block);
        i = begin_block;
    }
}