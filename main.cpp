#include "LagreNumbers.hpp"
#include <iostream>

int main() {
    // // Длины чисел
    // int size1 = 0;
    // int size2 = 0;

    // std::cout << "Enter size of first number: ";
    // std::cin >> size1;
    // std::cout << "Enter size of second number: ";
    // std::cin >> size2;

    // // Генерируем два числа в виде строки
    // std::string stringLN1 = generateLN(size1);
    // std::string stringLN2 = generateLN(size2);

    std::string stringLN1 = "+351270045";
    std::string stringLN2 = "+145";

    std::cout << "First number: " << stringLN1 << "\n";
    std::cout << "Second number: " << stringLN2 << "\n";

    // Парсим строки
    LargeNumber LN1(stringLN1);
    LargeNumber LN2(stringLN2);

    LargeNumber Sum = LNMath::sum(LN1, LN2);
    std::cout << "Sum = " << Sum.toString() << "\n";

    LargeNumber Diff = LNMath::sub(LN1, LN2);
    std::cout << "Diff = " << Diff.toString() << "\n";

    LargeNumber Product = LNMath::mult(LN1, LN2);
    std::cout << "Product = " << Product.toString() << "\n";

    LargeNumber Quotient = LNMath::div(LN1, LN2);
    std::cout << "Quotient = " << Quotient.toString() << "\n";

    LargeNumber Powered = LNMath::pow(LN1, LN2);
    std::cout << "Powered = " << Powered.toString() << "\n";

    LargeNumber LCM = LNMath::lcm(LN1, LN2);
    std::cout << "LCM = " << LCM.toString() << "\n";

    LargeNumber GCD = LNMath::gcd(LN1, LN2);
    std::cout << "GCD = " << GCD.toString() << "\n";

    return 0;
}