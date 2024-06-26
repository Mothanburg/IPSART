#ifndef GARS_UTIL
#define GARS_UTIL

#include <fstream>
#include <iostream>
#include <string>

// 函数：将非ASCII字符替换为'0'
inline void replaceNonASCII(const std::string& inputFilePath,
                            const std::string& outputFilePath) {
  // 以二进制模式打开输入文件
  std::ifstream inputFile(inputFilePath, std::ios::binary);
  if (!inputFile.is_open()) {
    std::cerr << "Failed to open input file: " << inputFilePath << std::endl;
    return;
  }

  // 以二进制模式打开输出文件
  std::ofstream outputFile(outputFilePath, std::ios::binary);
  if (!outputFile.is_open()) {
    std::cerr << "Failed to open output file: " << outputFilePath << std::endl;
    inputFile.close();
    return;
  }

  char ch;
  while (inputFile.get(ch)) {
    // 如果字符是非ASCII字符，则替换为'0'
    if (static_cast<unsigned char>(ch) > 127) {
      ch = '0';
    }
    // 写入输出文件
    outputFile.put(ch);
  }

  inputFile.close();
  outputFile.close();
}

#endif