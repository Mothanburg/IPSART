#ifndef GARS_UTIL
#define GARS_UTIL

#include <cmath>
#include <span>
#include <type_traits>


// // 函数：将非ASCII字符替换为'0'
// inline void replaceNonASCII(const std::string &inputFilePath,
//                             const std::string &outputFilePath) {
//   // 以二进制模式打开输入文件
//   std::ifstream inputFile(inputFilePath, std::ios::binary);
//   if (!inputFile.is_open()) {
//     std::cerr << "Failed to open input file: " << inputFilePath << std::endl;
//     return;
//   }

//   // 以二进制模式打开输出文件
//   std::ofstream outputFile(outputFilePath, std::ios::binary);
//   if (!outputFile.is_open()) {
//     std::cerr << "Failed to open output file: " << outputFilePath <<
//     std::endl; inputFile.close(); return;
//   }

//   char ch;
//   while (inputFile.get(ch)) {
//     // 如果字符是非ASCII字符，则替换为'0'
//     if (static_cast<unsigned char>(ch) > 127) {
//       ch = '0';
//     }
//     // 写入输出文件
//     outputFile.put(ch);
//   }

//   inputFile.close();
//   outputFile.close();
// }

template <typename Dtype>
struct PolMat3 {
  PolMat3(long long total_len, Dtype *raw_m11, Dtype *raw_m22, Dtype *raw_m33,
          Dtype *raw_m12_r, Dtype *raw_m13_r, Dtype *raw_m23_r,
          Dtype *raw_m12_i, Dtype *raw_m13_i, Dtype *raw_m23_i)
      : m11(raw_m11, total_len),
        m22(raw_m22, total_len),
        m33(raw_m33, total_len),
        m12_r(raw_m12_r, total_len),
        m13_r(raw_m13_r, total_len),
        m23_r(raw_m23_r, total_len),
        m12_i(raw_m12_i, total_len),
        m13_i(raw_m13_i, total_len),
        m23_i(raw_m23_i, total_len) {}

  std::span<Dtype> m11, m22, m33, m12_r, m13_r, m23_r, m12_i, m13_i, m23_i;
};

template <typename Tint>
  requires std::is_integral_v<Tint>
inline Tint int_sqrt(Tint x) {
  auto res = std::sqrt(static_cast<long double>(x));
  return static_cast<Tint>(std::round(res));
}

#endif