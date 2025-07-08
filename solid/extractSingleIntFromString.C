#include "extractSingleIntFromString.H"

int extractSingleIntFromString(const std::string& str) {  
    int result = -1; // 初始化结果为-1，表示没有找到数字  
    std::string numberStr;  
    bool foundDigit = false;  
  
    for (char ch : str) {  
        if (std::isdigit(ch) || (!foundDigit && ch == '-' && !numberStr.empty() == false)) {  
            // If we find a digit or (if we haven't found a digit yet and the character is a '-' sign)  
            // Note: The check for '-' ensures that '-' can only be at the start of the number.  
            numberStr += ch;  
            foundDigit = true; // Set foundDigit to true once we find the first digit (or '-' sign)  
        } else if (foundDigit && !std::isdigit(ch)) {  
            // If we find a non-digit character and we have already found a digit, stop collecting the number  
            break;  
        }  
    }  
  
    // 如果找到了有效的数字字符串，则将其转换为整数  
    if (!numberStr.empty()) {  
        // Use strtol to convert string to long, with error checking  
        char* endPtr = nullptr;  
        long number = std::strtol(numberStr.c_str(), &endPtr, 10);  
  
        // Check for conversion errors and ensure the number is within the int range  
        if (endPtr != numberStr.c_str() && *endPtr == '\0' && number >= INT_MIN && number <= INT_MAX) {  
            result = static_cast<int>(number);  
        }  
    }  
  
    return result;  
}  