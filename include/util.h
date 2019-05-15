#pragma once

#include <algorithm>
#include <vector>
#include <sstream>
#include <string>

inline std::vector<int> string_to_vecInt(std::string s) {
    int number;
    std::vector<int> vec;
    std::stringstream ss(s);
    while (ss >> number) {
        vec.push_back(number);
    }
    return vec;
}