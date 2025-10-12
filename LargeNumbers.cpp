#include "LagreNumbers.hpp"
#include <random>
#include <iostream>

std::string generateLN(int size) {
    static std::mt19937 gen(std::random_device{}());
    std::uniform_int_distribution<int> signDist(0, 1);
    std::uniform_int_distribution<int> firstDigitDist(1, 9);
    std::uniform_int_distribution<int> digitDist(0, 9);

    std::string LN;
    LN.reserve(size + 1);  // выделяем место под знак

    LN += (signDist(gen) == 0) ? '-' : '+';
    LN += static_cast<char>('0' + firstDigitDist(gen));  // без ведущего нуля

    for (int i = 1; i < size; ++i)
        LN += static_cast<char>('0' + digitDist(gen));

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

    if (large_number[0] == 0) {
        str.erase(0);  // ноль не положительное, и не отрицательное
    }

    str += std::to_string(large_number.back());

    for (int i = large_number.size() - 2; i >= 0; --i) {
        std::string block_str = std::to_string(large_number[i]);

        // Добавляем ведущие нули
        str += std::string(9 - block_str.length(), '0') + block_str;
    }

    return str;
}

LargeNumber LNMath::sum(const LargeNumber& a, const LargeNumber& b) {
    return a.positive == b.positive ? absSum(a, b) : absSub(a, b);
}

LargeNumber LNMath::sub(const LargeNumber& a, const LargeNumber& b) {
    LargeNumber result;

    if (a.positive == b.positive) {
        int cmp = compareLN(a, b);
        if (cmp >= 0) {
            result = absSub(a, b);
            result.positive = a.positive;
        } else {
            result = absSub(b, a);
            result.positive = !a.positive; // знак меняется
        }
    } else {
        result = absSum(a, b);
        result.positive = a.positive;
    }

    // удаляем ведущие нули
    while (result.large_number.size() > 1 && result.large_number.back() == 0)
        result.large_number.pop_back();

    return result;
}

LargeNumber LNMath::mult(const LargeNumber& a, const LargeNumber& b) {
    LargeNumber result;
    a.positive == b.positive ? result.positive = true : result.positive = false;

    if (a.large_number[0] == 0 || b.large_number[0] == 0) {
        result.large_number = {0};
        return result;
    }

    size_t A = a.large_number.size();
    size_t B = b.large_number.size();
    result.large_number.assign(A + B, 0);

    for (size_t i = 0; i < A; i++) {
        uint64_t carry = 0;

        for (size_t j = 0; j < B || carry > 0; j++) {
            uint64_t blockA = a.large_number[i];
            uint64_t blockB = (j < B) ? b.large_number[j] : 0;

            uint64_t cur = result.large_number[i + j] + blockA * blockB + carry;

            result.large_number[i + j] = cur % BASE;
            carry = cur / BASE;
        }
    }

    // Убираем ведущие нули
    while (result.large_number.size() > 1 && result.large_number.back() == 0)
        result.large_number.pop_back();

    return result;
}

LargeNumber LNMath::div(const LargeNumber& a, const LargeNumber& b) {
    LargeNumber result;
    result.positive = (a.positive == b.positive);

    if (b.large_number[0] == 0) {
        throw std::runtime_error("Error: division by zero!");
    }

   if (compareLN(a, b) < 0) {
        result.large_number.push_back(0);
        result.positive = true;
        return result;
    }

    LargeNumber dividend = a;
    LargeNumber divisor = b;
    LargeNumber current;
    current.large_number.clear();

    result.large_number.assign(dividend.large_number.size(), 0);

    for (size_t i = dividend.large_number.size(); i-- > 0;) {
        // Сдвигаем блоки в current
        current.large_number.insert(current.large_number.begin(), dividend.large_number[i]);

        // Убираем ведущие нули
        while (current.large_number.size() > 1 && current.large_number.back() == 0)
            current.large_number.pop_back();

        uint32_t left = 0, right = BASE - 1, x = 0;

        // Бинарный поиск коэффициента x
        while (left <= right) {
            uint32_t mid = (left + right) / 2;

            LargeNumber temp;
            temp.large_number.push_back(mid);
            temp = mult(divisor, temp);

            if (compareLN(temp, current) <= 0) {
                x = mid;
                left = mid + 1;
            } else {
                right = mid - 1;
            }
        }

        result.large_number[i] = x;

        LargeNumber temp;
        temp.large_number.push_back(x);
        current = absSub(current, mult(divisor, temp));
    }

    // Убираем ведущие нули
    while (result.large_number.size() > 1 && result.large_number.back() == 0)
        result.large_number.pop_back();

    return result;
}

int LNMath::compareLN(const LargeNumber& a, const LargeNumber& b) {
    if (a.large_number.size() > b.large_number.size()) {
        return 1;
    } else if (a.large_number.size() < b.large_number.size()) {
        return -1;
    } else {
        for (int i = a.large_number.size() - 1; i >= 0; --i) {
            if (a.large_number[i] > b.large_number[i]) return 1;
            if (a.large_number[i] < b.large_number[i]) return -1;
        }
    }

    return 0;
}

LargeNumber LNMath::absSum(const LargeNumber& a, const LargeNumber& b) {
    LargeNumber result;
    result.positive = (a.positive == b.positive);

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
        result.positive = a.positive;
        bigger = &a;
        smaller = &b;
    } else if (compare == -1) {
        result.positive = b.positive;
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