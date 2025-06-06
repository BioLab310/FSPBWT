//
// Created by cui on 25-6-6.
//
#ifndef TOOLS_H
#define TOOLS_H
#include <chrono>
#include<iostream>
#include <algorithm>
#include <cstring>
#include <fstream>

#include <vector>
#include <ctime>
#include <string>

#define FF 4
using namespace std;


template<typename T>
void saveStringVector(ofstream &out, const vector<T> &vec) {
    size_t size = vec.size();
    out.write((char *) &size, sizeof(size)); // Write the size of the vector
    for (const auto &str: vec) {
        size_t len = str.size(); // Get the length of the string
        out.write((char *) &len, sizeof(len)); // Write the length of the string
        out.write(str.c_str(), len); // Write the string data
    }
}

template<typename Stream, typename T>
void saveVectorWithDimension(Stream &out, const vector<T> &vec) {
    size_t size = vec.size();
    out.write((char *) &size, sizeof(size)); // Write the size of the vector
    if (size > 0) {
        size_t dim = sizeof(T); // Write the dimension of each element in the vector
        out.write((char *) &dim, sizeof(dim));
        out.write((char *) &vec[0], size * sizeof(T)); // Write the elements
    }
}

template<typename Stream, typename T>
void saveVectorVectorWithDimension(Stream &out,
                                   const vector<vector<T> > &vec) {
    size_t size = vec.size();
    out.write((char *) &size, sizeof(size)); // Write the size of the vector of vectors
    for (const auto &subVec: vec) {
        saveVectorWithDimension(out, subVec); // Save each subvector with dimension info
    }
}

template<typename Stream, typename T>
void saveVectorVectorVectorWithDimension(Stream &out,
                                         const vector<vector<vector<T> > > &vec) {
    size_t size = vec.size();
    out.write((char *) &size, sizeof(size)); // Write the size of the vector of vectors of vectors
    for (const auto &matrix: vec) {
        saveVectorVectorWithDimension(out, matrix); // Save each matrix with dimension info
    }
}

template<typename T>
void loadStringVector(ifstream &in, vector<T> &vec) {
    size_t size;
    in.read((char *) &size, sizeof(size)); // Read the size of the vector
    vec.resize(size); // Resize the vector
    for (auto &str: vec) {
        size_t len;
        in.read((char *) &len, sizeof(len)); // Read the length of the string
        str.resize(len); // Resize the string
        in.read(&str[0], len); // Read the string data
    }
}

template<typename Stream, typename T>
void loadVectorWithDimension(Stream &in, vector<T> &vec) {
    size_t size;
    in.read((char *) &size, sizeof(size)); // Read the size of the vector
    if (size > 0) {
        size_t dim;
        in.read((char *) &dim, sizeof(dim)); // Read the dimension of each element
        vec.resize(size); // Resize the vector
        in.read((char *) &vec[0], size * dim); // Read the elements
    }
}

template<typename Stream, typename T>
void loadVectorVectorWithDimension(Stream &in, vector<vector<T> > &vec) {
    size_t size;
    in.read((char *) &size, sizeof(size)); // Read the size of the vector of vectors
    vec.resize(size); // Resize the outer vector
    for (auto &subVec: vec) {
        loadVectorWithDimension(in, subVec); // Load each subvector
    }
}

template<typename Stream, typename T>
void loadVectorVectorVectorWithDimension(Stream &in,
                                         vector<vector<vector<T> > > &vec) {
    size_t size;
    in.read((char *) &size, sizeof(size)); // Read the size of the vector of vectors of vectors
    vec.resize(size); // Resize the outer vector
    for (auto &matrix: vec) {
        loadVectorVectorWithDimension(in, matrix); // Load each matrix
    }
}

template<typename T>
void print2DVector(const std::vector<std::vector<T> > &vec) {
    for (const auto &row: vec) {
        for (const auto &elem: row) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }
}


int ctz128_uint128(__uint128_t num) {
    uint64_t low_bits = num;
    uint64_t high_bits = num >> 64;

    if (low_bits != 0) {
        return __builtin_ctzll(low_bits);
    } else if (high_bits != 0) {
        return 64 + __builtin_ctzll(high_bits);
    } else {
        return 128; // 如果输入为 0，则返回 128
    }
}

int clz128_uint128(__uint128_t num) {
    uint64_t low_bits = num;
    uint64_t high_bits = num >> 64;

    if (high_bits != 0) {
        return __builtin_clzll(high_bits);
    } else if (low_bits != 0) {
        return 64 + __builtin_clzll(low_bits);
    } else {
        return 128; // 如果输入为 0，则返回 128
    }
}

// 为 unsigned __int128 类型定义 operator<<
std::ostream &operator<<(std::ostream &os, const unsigned __int128 &value) {
    // 由于 __int128 是一个非常大的整数，我们需要转换成字符串来输出
    std::string output;
    unsigned __int128 temp = value;
    do {
        output += '0' + (temp % 10);
        temp /= 10;
    } while (temp != 0);
    std::reverse(output.begin(), output.end());
    os << output;
    return os;
}

int countSetBits128(const __uint128_t &n) {
    uint64_t low_part = static_cast<uint64_t>(n); // 低位部分
    uint64_t high_part = static_cast<uint64_t>(n >> 64); // 高位部分
    return __builtin_popcountll(low_part) + __builtin_popcountll(high_part);
}

// 计算 std::vector<string> 的大小
size_t calculateSize(const std::vector<std::string> &vec) {
    size_t totalSize = vec.capacity() * sizeof(std::string); // 存储 string 对象的空间
    for (const auto &str: vec) {
        totalSize += str.capacity(); // 每个 string 对象实际存储字符串的空间
    }
    return totalSize;
}

// 计算 std::vector<int> 的大小
size_t calculateSize(const std::vector<int> &vec) {
    return vec.capacity() * sizeof(int);
}

// 计算 std::vector<std::vector<T>> 的大小
template<typename T>
size_t calculateSize(const std::vector<std::vector<T> > &vec) {
    size_t totalSize = vec.capacity() * sizeof(std::vector<T>); // 存储内层 vector 的空间
    for (const auto &innerVec: vec) {
        totalSize += innerVec.capacity() * sizeof(T); // 每个内层 vector 存储元素的空间
    }
    return totalSize;
}




#endif //TOOLS_H
