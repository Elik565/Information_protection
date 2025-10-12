#include "LagreNumbers.hpp"
#include <random>
#include <iostream>

std::string generateLN(const int size) {
    std::string LN;
    LN.reserve(size);  // сразу выделяем память

    char sign = (rand() % 1 == 0) ? '-' : '+';
    LN += sign;

    LN += '1' + rand() % 9;  // предотвращаем ведущие нули

    for (int i = 0; i < size - 1; i++) {
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

std::string LargeNumber::toString() {
    std::string str = "";

    positive ? str +=  "+" : str += "-";

    str += std::to_string(large_number.back());

    for (int i = large_number.size() - 2; i >= 0; --i) {
        std::string block_str = std::to_string(large_number[i]);

        // Добавляем ведущие нули
        str += std::string(9 - block_str.length(), '0') + block_str;
    }

    return str;
}

int LNMath::compareLN(const LargeNumber& a, const LargeNumber& b) {
    if (a.large_number.size() > b.large_number.size()) return 1;
    else if (a.large_number.size() < b.large_number.size()) return -1;
    else {
        for (int i = a.large_number.size() - 1; i >= 0; --i) {
            if (a.large_number[i] > b.large_number[i]) return 1;
            if (a.large_number[i] < b.large_number[i]) return -1;
        }
    }

    return 0;
}

LargeNumber LNMath::sum(const LargeNumber& a, const LargeNumber& b) {
    return a.positive == b.positive ? absSum(a, b) : absSub(a, b);
}

LargeNumber LNMath::absSum(const LargeNumber& a, const LargeNumber& b) {
    LargeNumber result;

    a.positive ? result.positive = true : result.positive = false;

    uint64_t carry = 0;  // для переноса
    
    // Находим максимальное количество блоков
    size_t repeat = std::max(a.large_number.size(), b.large_number.size());
    
    for (size_t i = 0; i < repeat || carry > 0; i++) {
        uint64_t blockA = (i < a.large_number.size()) ? a.large_number[i] : 0;
        uint64_t blockB = (i < b.large_number.size()) ? b.large_number[i] : 0;
        
        uint64_t sum = blockA + blockB + carry;
        
        result.large_number.push_back(sum % BASE);
        carry = sum / BASE;
    }

    return result;
}

LargeNumber LNMath::absSub(const LargeNumber& a, const LargeNumber& b) {
    LargeNumber result;

    int compare = compareLN(a, b);
    const LargeNumber* bigger;
    const LargeNumber* smaller;

    if (compare == 1) {
        bigger = &a;
        smaller = &b;
    } else if (compare == -1) {
        bigger = &b;
        smaller = &a;
    } else {
        result.large_number.push_back(0);
        return result;
    }

    int64_t carry = 0;
    for (size_t i = 0; i < bigger->large_number.size(); ++i) {
        int64_t cur = (int64_t)bigger->large_number[i] - carry;
        if (i < smaller->large_number.size()) {
            cur -= smaller->large_number[i];
        }

        if (cur < 0) {
            cur += BASE; 
            carry = 1;
        } else {
            carry = 0;
        }

        result.large_number.push_back((uint32_t)cur);
    }

    return result;
}