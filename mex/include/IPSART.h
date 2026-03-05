#pragma once

#ifndef _IPSART_
#define _IPSART_

#include <string>
#include <vector>

#include <mex.hpp>

namespace ipsart {

// ============================================================================
// IPSART函数接口
// ============================================================================

/**
 * @brief Cloude-Pottier 分解
 * @param input 极化协方差矩阵元素 (C2: 4个, C3: 9个)
 * @param af MATLAB ArrayFactory
 * @return {H, a, A} - 熵、Alpha角、各向异性
 */
std::vector<matlab::data::Array>
CloudePottier(const std::vector<matlab::data::Array> &input,
              matlab::data::ArrayFactory &af);

/**
 * @brief Yamaguchi 四分量分解
 * @param input 极化协方差矩阵元素 (9个)
 * @param af MATLAB ArrayFactory
 * @return {Ps, Pd, Pv, Ph} - 表面散射、偶次散射、体散射、螺旋散射
 */
std::vector<matlab::data::Array>
Yamaguchi(const std::vector<matlab::data::Array> &input,
          matlab::data::ArrayFactory &af);

/**
 * @brief General Four component decomposition with Unitary transformation (G4U)
 * 四分量分解
 * @param input 极化协方差矩阵元素 (9个)
 * @param af MATLAB ArrayFactory
 * @return {Ps, Pd, Pv, Ph} - 表面散射、偶次散射、体散射、螺旋散射
 */
std::vector<matlab::data::Array>
G4U(const std::vector<matlab::data::Array> &input,
    matlab::data::ArrayFactory &af);

/**
 * @brief 多视处理 (Multi-looking)
 * @param input 输入图像
 * @param af MATLAB ArrayFactory
 * @return 多视处理后的图像
 */
std::vector<matlab::data::Array>
Multilook(const std::vector<matlab::data::Array> &input,
          matlab::data::ArrayFactory &af);

/**
 * @brief Refined Lee 极化滤波
 * @param input 极化矩阵元素 + look_num (C2: 5个, C3: 10个)
 * @param af MATLAB ArrayFactory
 * @return 滤波后的极化矩阵元素
 */
std::vector<matlab::data::Array>
RefinedLeeFilter(const std::vector<matlab::data::Array> &input,
                 matlab::data::ArrayFactory &af);

// ============================================================================
// 工具函数
// ============================================================================

namespace ocl {
/// 获取 OpenCL 设备名称
std::string getOpenCLDeviceName();

/// 获取 OpenCL 设备版本
std::string getOpenCLDeviceVersion();

} // namespace ocl
} // namespace ipsart

#endif