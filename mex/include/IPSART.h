#pragma once

#ifndef _IPSART_
#define _IPSART_

#include <vector>

#include <mex.hpp>

std::vector<matlab::data::Array>
CloudePottier(const std::vector<matlab::data::Array> &input,
              matlab::data::ArrayFactory &af);

std::vector<matlab::data::Array>
Yamaguchi(const std::vector<matlab::data::Array> &input,
          matlab::data::ArrayFactory &af);

std::vector<matlab::data::Array>
G4U(const std::vector<matlab::data::Array> &input,
    matlab::data::ArrayFactory &af);

std::vector<matlab::data::Array>
Multilook(const std::vector<matlab::data::Array> &input,
          matlab::data::ArrayFactory &af);

std::vector<matlab::data::Array>
RefinedLeeFilter(const std::vector<matlab::data::Array> &input,
                 matlab::data::ArrayFactory &af);

#endif