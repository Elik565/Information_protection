#include "LagreNumbers.hpp"
#include <random>
#include <complex>
#include <cmath>
#include <iostream>
#include <algorithm>

const double PI = acos(-1);

void fft(std::vector<std::complex<double>>& a, bool invert) {
    int n = a.size();
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) swap(a[i], a[j]);
    }

    for (int len = 2; len <= n; len <<= 1) {
        double ang = 2 * PI / len * (invert ? -1 : 1);
        std::complex<double> wlen(cos(ang), sin(ang));
        for (int i = 0; i < n; i += len) {
            std::complex<double> w(1);
            for (int j = 0; j < len / 2; j++) {
                std::complex<double> u = a[i + j];
                std::complex<double> v = a[i + j + len / 2] * w;
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }

    if (invert)
        for (std::complex<double>& x : a) x /= n;
}

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

LargeNumber::LargeNumber(int64_t LN) {
    if (LN < 0) {
        positive = false;
        LN = -LN;  // берём модуль для хранения блоков
    }

    while (LN > 0) {
        large_number.push_back(static_cast<uint32_t>(LN % BASE));
        LN /= BASE;
    }

    if (large_number.empty()) {
        large_number.push_back(0);
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

    return result;
}

LargeNumber LNMath::mult(const LargeNumber& a, const LargeNumber& b) {
    std::vector<int> A, B;
    for (auto x : a.large_number) {
        for (int i = 0; i < 9; i++) {
            A.push_back(x % 10);
            x /= 10;
        }
    }
    for (auto x : b.large_number) {
        for (int i = 0; i < 9; i++) {
            B.push_back(x % 10);
            x /= 10;
        }
    }

    int n = 1;
    while (n < (int)A.size() + (int)B.size()) n <<= 1;

    std::vector<std::complex<double>> fa(A.begin(), A.end()), fb(B.begin(), B.end());
    fa.resize(n); fb.resize(n);

    fft(fa, false);
    fft(fb, false);
    for (int i = 0; i < n; i++) fa[i] *= fb[i];
        fft(fa, true);

    std::vector<uint64_t> res(n);
    for (int i = 0; i < n; i++) res[i] = (uint64_t)(fa[i].real() + 0.5);

    uint64_t carry = 0;
    for (int i = 0; i < n; i++) {
        uint64_t cur = res[i] + carry;
        res[i] = cur % 10;
        carry = cur / 10;
    }

    while (carry) {
        res.push_back(carry % 10);
        carry /= 10;
    }

    LargeNumber result;
    result.positive = (a.positive == b.positive);
    uint64_t block = 0, p = 1;
    for (size_t i = 0; i < res.size(); i++) {
        block += res[i] * p;
        p *= 10;
        if (p == BASE) {
            result.large_number.push_back(block);
            block = 0;
            p = 1;
        }
    }
    if (block > 0) result.large_number.push_back(block);

    while (result.large_number.size() > 1 && result.large_number.back() == 0)
        result.large_number.pop_back();

    return result;
}

LargeNumber LNMath::div(const LargeNumber& a, const LargeNumber& b) {
    LargeNumber result;
    result.positive = (a.positive == b.positive);

    if (b.large_number.size() == 1 && b.large_number[0] == 0)
        throw std::runtime_error("Error: division by zero!");

    // Меньшее делим на большее
    if (compareLN(a, b) < 0) {
        result.large_number.push_back(0);
        result.positive = true;
        return result;
    }

    LargeNumber dividend = a;
    LargeNumber divisor = b;
    LargeNumber current;
    result.large_number.assign(dividend.large_number.size(), 0);

    // Проходим с конца большего числа (делимого)
    for (int i = dividend.large_number.size() - 1; i >= 0; --i) {
        // Сдвигаем current влево и добавляем старший блок
        current.large_number.insert(current.large_number.begin(), dividend.large_number[i]);

        while (current.large_number.size() > 1 && current.large_number.back() == 0)
            current.large_number.pop_back();

        uint32_t left = 0, right = BASE - 1, x = 0;

        // Находим сколько раз делитель помещается в текущем остатке
        // Бинарный поиск
        while (left <= right) {
            uint32_t mid = (left + right) / 2;
            LargeNumber LNmid;
            LNmid.large_number.push_back(mid);
            LargeNumber temp = mult(divisor, LNmid);

            // Проверяем, что не вышли за границы
            if (compareLN(temp, current) <= 0) {
                x = mid;
                left = mid + 1;
            } else {
                right = mid - 1;
            }
        }

        result.large_number[i] = x;
        LargeNumber LNx;
        LNx.large_number.push_back(x);
        LargeNumber temp = mult(divisor, LNx);
        current = absSub(current, temp);
    }

    while (result.large_number.size() > 1 && result.large_number.back() == 0)
        result.large_number.pop_back();

    return result;
}

LargeNumber LNMath::pow(const LargeNumber& a, const LargeNumber& b) {
    LargeNumber result;

    if (a.large_number[0] == 1 && !b.positive) {
        result.large_number.push_back(1);
        return result;
    }

    if (!b.positive) {
        throw std::runtime_error("Error: power is negative!");
    }
    
    if (b.large_number[0] == 0) {
        result.large_number.push_back(1);
        return result;
    }

    LargeNumber base = a;
    LargeNumber pow = b;
    LargeNumber zero;
    zero.large_number.push_back(0);
    result.large_number.push_back(1);

    while (compareLN(pow, zero) > 0) {
        // Если степень нечетная, то умножаем результат на base
        if (pow.large_number[0] % 2 == 1) {
            result = mult(result, base);
        }
        base = mult(base, base);
        
        // Делим степень на 2
        LargeNumber two;
        two.large_number.push_back(2);
        pow = div(pow, two);

        // Убираем ведущие нули
        while (pow.large_number.size() > 1 && pow.large_number.back() == 0)
            pow.large_number.pop_back();
    }

    result.positive = (b.large_number[0] % 2 == 0) ? true : a.positive;
    return result;
}

LargeNumber LNMath::sqrt(const LargeNumber& a) {
    if (a.large_number.size() == 1 && a.large_number[0] == 0)
        return a;

    LargeNumber x = a; 
    LargeNumber two;
    two.large_number.push_back(2);

    LargeNumber prev;

    while (true) {
        LargeNumber nx = div(sum(x, div(a, x)), two);

        if (compareLN(nx, x) == 0 || compareLN(nx, prev) == 0) {
            break;
        }

        prev = x;
        x = nx;
    }

    LargeNumber one;
    one.large_number.push_back(1);
    while (compareLN(mult(x, x), a) > 0) {
        x = sub(x, one);
    }

    return x;
}

LargeNumber LNMath::gcd(const LargeNumber& a, const LargeNumber& b) {
    LargeNumber A = a;
    LargeNumber B = b;

    while (!(B.large_number.size() == 1 && B.large_number[0] == 0)) {
        LargeNumber temp = mod(A, B);
        
        A = B;
        B = temp;
    }

    return A;
}

LargeNumber LNMath::lcm(const LargeNumber& a, const LargeNumber& b) {
    LargeNumber gcd_val = gcd(a, b);
    LargeNumber mult_val = mult(a, b);
    LargeNumber result = div(mult_val, gcd_val);
    result.positive = true;
    return result;
}

bool LNMath::isPrimeStd(const LargeNumber& a) {
    if (!a.positive || (a.large_number.size() == 1 && (a.large_number[0] == 0 || a.large_number[0] == 1))) {
        return 0;
    }

    if (a.large_number[0] % 2 == 0 && a.large_number.size() != 1) {
        return false;
    }

    if (a.large_number[0] % 10 == 0 || a.large_number[0] % 10 == 5) {
        return false;
    }

    // Проверка на деление суммы цифр на 3 и 9
    int sum_digits = 0;
    LargeNumber copy = a;
    std::string s = copy.toString();
    for (size_t i = 1; i < s.size(); ++i) {
        sum_digits += s[i] - '0';
    }

    if (sum_digits % 3 == 0  || sum_digits % 9 == 0) {
        return (s == "+3");  // только 3 простое число
    }

    LargeNumber i;
    i.large_number.push_back(3);  // начинаем с 3
    LargeNumber one; one.large_number.push_back(1);
    LargeNumber step; step.large_number.push_back(2);  // проверяем только нечётные

    LargeNumber sqrt = LNMath::sqrt(a);

    while (compareLN(i, sqrt) <= 0) {
        LargeNumber mod_res = LNMath::mod(a, i);
        if (mod_res.large_number.size() == 1 && mod_res.large_number[0] == 0)
            return false;
        i = LNMath::sum(i, step);
    }

    return true;
}

bool LNMath::sieveEratosthenes(const LargeNumber& a) {
    if (compareLN(a, LargeNumber(2)) < 0) return false;
    if (compareLN(a, LargeNumber(2)) == 0) return true;
    if (compareLN(a, LargeNumber(3)) == 0) return true;
    if (mod(a, LargeNumber(2)).large_number[0] == 0) return false;

    LargeNumber limit = sqrt(a);  // предел проверки
    LargeNumber i(3);
    LargeNumber two(2);

    while (compareLN(i, limit) <= 0) {
        if (mod(a, i).large_number[0] == 0) return false;
        i = absSum(i, two);
    }

    return true;
}

bool LNMath::sieveAtkin(const LargeNumber& a) {
    if (compareLN(a, LargeNumber(2)) < 0) return false;
    if (compareLN(a, LargeNumber(2)) == 0) return true;
    if (compareLN(a, LargeNumber(3)) == 0) return true;
    if (mod(a, LargeNumber(2)).large_number[0] == 0) return false;
    if (mod(a, LargeNumber(3)).large_number[0] == 0) return false;

    if (mod(a, LargeNumber(2)).large_number[0] == 0) return false;

    LargeNumber limit = sqrt(a);  // предел проверки
    LargeNumber one(1), two(2), three(3), four(4);

    for (uint64_t x = 1; ; ++x) {
        LargeNumber X2 = LargeNumber(x * x);
        if (compareLN(X2, limit) > 0) break;

        for (uint64_t y = 1; ; ++y) {
            LargeNumber Y2 = LargeNumber(y * y);
            LargeNumber n;

            n = absSum(mult(four, X2), Y2);
            if (compareLN(n, limit) > 0) break;
            if (n.large_number[0] % 12 == 1 || n.large_number[0] % 12 == 5)
                if (mod(a, n).large_number[0] == 0) return false;

            n = absSum(mult(three, X2), Y2);
            if (compareLN(n, limit) > 0) break;
            if (n.large_number[0] % 12 == 7)
                if (mod(a, n).large_number[0] == 0) return false;

            if (x > y) {
                n = absSub(mult(three, X2), Y2);
                if (compareLN(n, limit) > 0) break;
                if (n.large_number[0] % 12 == 11)
                    if (mod(a, n).large_number[0] == 0) return false;
            }
        }
    }

    return true;
}

bool LNMath::LucasLehmer(const LargeNumber& a, uint64_t p) {
    if (p < 2) return false;

    // Проверяем, что p — простое
    if (!isPrimeStd(LargeNumber(p))) return false;

    LargeNumber s(4);

    for (uint64_t i = 0; i < p - 2; ++i) {
        s = absSub(mult(s, s), LargeNumber(2));
        s = mod(s, a);
    }

    return compareLN(s, LargeNumber(0)) == 0;
}

bool LNMath::MillerRabin(const LargeNumber& a, int iterations) {
    if (compareLN(a, LargeNumber(2)) < 0) return false;
    if (compareLN(a, LargeNumber(2)) == 0) return true;
    if (a.large_number[0] % 2 == 0) return false;

    LargeNumber d = absSub(a, LargeNumber(1));
    int s = 0;
    while (d.large_number[0] % 2 == 0) {
        d = div(d, LargeNumber(2));
        s++;
    }

    std::mt19937_64 gen(std::random_device{}());

    for (int iter = 0; iter < iterations; ++iter) {
        LargeNumber range = absSub(a, LargeNumber(4));
        uint64_t rnd = gen();

        uint64_t r = rnd % 1000000 + 2;
        if (compareLN(LargeNumber(r), range) > 0)
            r = 2;

        LargeNumber a_test(r);
        LargeNumber x = mod(pow(a_test, d), a);

        if (compareLN(x, LargeNumber(1)) == 0 ||
            compareLN(x, absSub(a, LargeNumber(1))) == 0)
            continue;

        bool composite = true;

        for (int r_i = 1; r_i < s; ++r_i) {
            x = mod(mult(x, x), a);
            if (compareLN(x, absSub(a, LargeNumber(1))) == 0) {
                composite = false;
                break;
            }
            if (compareLN(x, LargeNumber(1)) == 0)
                return false;
        }

        if (composite)
            return false;
    }

    return true;
}

int LNMath::jacobi(LargeNumber a, const LargeNumber& n_const) {
    if (!n_const.positive || n_const.large_number.empty() || 
        (n_const.large_number[0] % 2 == 0)) {
        return 0;
    }
    
    LargeNumber n = n_const;
    
    a = mod(a, n);

    if (compareLN(a, LargeNumber(0)) == 0) {
        return 0;
    }
    
    int result = 1;
    
    while (compareLN(a, LargeNumber(0)) != 0) {
        while (a.large_number[0] % 2 == 0) {
            uint32_t n_mod_8 = mod(n, LargeNumber(8)).large_number[0];
            
            if (n_mod_8 == 3 || n_mod_8 == 5) {
                result = -result;
            }
            
            a = div(a, LargeNumber(2));
            if (compareLN(a, LargeNumber(0)) == 0) {
                break;
            }
        }
        
        if (compareLN(a, LargeNumber(0)) == 0) {
            break;
        }
        
        // Квадратичный закон взаимности
        uint32_t a_mod_4 = mod(a, LargeNumber(4)).large_number[0];
        uint32_t n_mod_4 = mod(n, LargeNumber(4)).large_number[0];
        
        if (a_mod_4 == 3 && n_mod_4 == 3) {
            result = -result;
        }
        
        LargeNumber temp = mod(n, a);
        n = a;
        a = temp;
    }
    
    if (compareLN(n, LargeNumber(1)) == 0) {
        return result;
    }
    
    return 0;
}

LucasParams LNMath::findLucasParameters(const LargeNumber& a) {
    LucasParams params;
    if (!a.positive || a.large_number.empty() || 
        compareLN(a, LargeNumber(5)) < 0 ||
        (a.large_number[0] % 2 == 0)) {
        params.P = LargeNumber(0);
        params.Q = LargeNumber(0);
        params.D = LargeNumber(0);
        return params;
    }
    
    LargeNumber P(3), Q(1), D;
    LargeNumber two(2), four(4), five(5);
    
    int max_trials = 100;
    int trial = 0;
    
    while (trial++ < max_trials) {
        LargeNumber P_squared = mult(P, P);
        LargeNumber four_Q = mult(four, Q);
        D = sub(P_squared, four_Q);
        
        int jac = jacobi(D, a);
        
        if (jac == -1) {
            params.P = P;
            params.Q = Q;
            params.D = D;
            return params;
        }
        
        if (jac == 0) {
            params.P = LargeNumber(0);
            params.Q = LargeNumber(0);
            params.D = LargeNumber(0);
            return params;
        }

        if (compareLN(Q, LargeNumber(1)) == 0) {
            Q = LargeNumber(-1);
        } else if (compareLN(Q, LargeNumber(-1)) == 0) {
            Q = LargeNumber(2);
        } else if (compareLN(Q, LargeNumber(2)) == 0) {
            Q = LargeNumber(-2);
        } else if (compareLN(Q, LargeNumber(-2)) == 0) {
            Q = LargeNumber(3);
            P = sum(P, LargeNumber(2));
        } else {
            if (Q.positive) {
                Q.positive = false;
            } else {
                Q.positive = true;
                Q = sum(Q, LargeNumber(2));
            }
        }
    }
    
    params.P = LargeNumber(0);
    params.Q = LargeNumber(0);
    params.D = LargeNumber(0);
    return params;
}

bool LNMath::StrongLucasTest(const LargeNumber& n) {
    if (compareLN(n, LargeNumber(2)) < 0) return false;
    if (compareLN(n, LargeNumber(3)) == 0) return true;
    if (compareLN(n, LargeNumber(2)) == 0) return true;
    if (n.large_number[0] % 2 == 0) return false;
    
    LucasParams params = findLucasParameters(n);
    
    LargeNumber m = sum(n, LargeNumber(1));
    LargeNumber d = m; 
    LargeNumber s(0);
    LargeNumber two(2);
    LargeNumber zero(0);
    LargeNumber one(1);
    
    while (compareLN(d, zero) > 0 && d.large_number[0] % 2 == 0) {
        d = div(d, two);
        s = sum(s, one);
    }

    LargeNumber inverse_two = div(sum(n, one), two);
    
    std::vector<bool> bits;
    LargeNumber temp_d = d;
    
    while (compareLN(temp_d, zero) > 0) {
        LargeNumber rem = mod(temp_d, two);
        bits.push_back(compareLN(rem, one) == 0);
        temp_d = div(temp_d, two);
    }
    
    std::reverse(bits.begin(), bits.end());
    
    LargeNumber U = LargeNumber(0);
    LargeNumber V = LargeNumber(2);      
    LargeNumber U1 = LargeNumber(1);     
    LargeNumber V1 = params.P;           
    LargeNumber Qk = params.Q;           
    
    for (size_t i = 1; i < bits.size(); ++i) {
        LargeNumber U2 = mod(mult(U1, V1), n);
        LargeNumber V2 = mod(sub(mult(V1, V1), 
                               mult(two, Qk)), n);
        Qk = mod(mult(Qk, Qk), n);
        
        if (bits[i]) {
            LargeNumber term1 = mod(mult(params.P, U2), n);
            LargeNumber U_next = mod(mult(sum(term1, V2), inverse_two), n);
            
            LargeNumber term2 = mod(mult(params.D, U2), n);
            LargeNumber term3 = mod(mult(params.P, V2), n);
            LargeNumber V_next = mod(mult(sum(term2, term3), inverse_two), n);
            
            Qk = mod(mult(Qk, params.Q), n);
            
            U1 = U_next;
            V1 = V_next;
        } else {
            U1 = U2;
            V1 = V2;
        }
    }
    
    if (compareLN(mod(U1, n), zero) == 0) {
        return true;
    }
    
    LargeNumber current_V = V1;
    for (LargeNumber r(0); compareLN(r, s) < 0; r = sum(r, one)) {
        if (compareLN(mod(current_V, n), zero) == 0) {
            return true;
        }
        current_V = mod(sub(mult(current_V, current_V), two), n);
    }
    
    return false;
}


// Auxiliary functions
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

    uint32_t carry = 0;  // для переноса
    
    // Находим максимальное количество блоков
    size_t repeat = std::max(a.large_number.size(), b.large_number.size());
    
    for (size_t i = 0; i < repeat || carry > 0; i++) {
        uint32_t blockA = (i < a.large_number.size()) ? a.large_number[i] : 0;
        uint32_t blockB = (i < b.large_number.size()) ? b.large_number[i] : 0;
        
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
        result.positive = !a.positive;
        bigger = &b;
        smaller = &a;
    } else {
        result.large_number.push_back(0);
        return result;
    }

    uint32_t carry = 0;
    for (size_t i = 0; i < bigger->large_number.size(); ++i) {
        int64_t cur = (int64_t)bigger->large_number[i] - carry;
        if (i < smaller->large_number.size()) {  // проверяем, есть ли у меньшего числа еще блок
            cur -= smaller->large_number[i];  // вычитаем из блока большего числа блок меньшего
        }

        if (cur < 0) {  // если результат получился отрицательным
            cur += BASE; 
            carry = 1;  // берем у старшего блока единицу
        } else {
            carry = 0;
        }

        result.large_number.push_back(cur);
    }

    return result;
}

LargeNumber LNMath::mod(const LargeNumber& a, const LargeNumber& b) {
    LargeNumber integer_part = div(a, b);  // целая часть
    LargeNumber product = mult(integer_part, b);  // произведение
    LargeNumber remainder = sub(a, product);  // остаток

    // Убираем ведущие нули
    while (remainder.large_number.size() > 1 && remainder.large_number.back() == 0)
        remainder.large_number.pop_back();

    return remainder;
}