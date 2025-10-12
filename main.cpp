#include "LagreNumbers.hpp"
#include <iostream>

int main() {
    // Длины чисел
    int size1 = 0;
    int size2 = 0;

    std::cout << "Enter size of first number: ";
    std::cin >> size1;
    std::cout << "Enter size of second number: ";
    std::cin >> size2;

    // Генерируем два числа в виде строки
    std::string stringLN1 = generateLN(size1);
    std::string stringLN2 = generateLN(size2);

    
    std::cin >> stringLN1;
    std::cin >> stringLN2;

    std::cout << "First number: " << stringLN1 << "\n";
    std::cout << "Second number: " << stringLN2 << "\n";

    // Парсим строки
    LargeNumber LN1(stringLN1);
    LargeNumber LN2(stringLN2);

    LargeNumber Sum = LNMath::sum(LN1, LN2);
    std::cout << "Sum = " << Sum.toString() << "\n";

    return 0;
}