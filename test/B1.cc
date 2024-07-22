#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <unordered_map>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <sys/time.h>
#include "psm-helper.hpp"

// 生成Simhash值及其位翻转组合
std::pair<std::vector<uint64_t>, uint64_t> generate_bit_flips(const std::vector<double> &feature_vector)
{
    int bit_count = 64;
    std::vector<double> bit_weights(bit_count, 0);

    // 根据特征向量计算每个位的权重
    for (int i = 0; i < 128; ++i)
    {
        uint64_t hash = PSMHelper::MurmurHash64A(std::to_string(i));
        for (int j = 0; j < bit_count; ++j)
        {
            if (hash & (1ULL << j))
            {
                bit_weights[j] += feature_vector[i];
            }
            else
            {
                bit_weights[j] -= feature_vector[i];
            }
        }
    }

    uint64_t simhash_value = 0;
    std::pair<double, int> first_bit{0.0, -1};
    std::pair<double, int> second_bit{0.0, -1};

    // 计算Simhash值，并找到前两位最大的权重值及其位置
    for (int i = 0; i < bit_count; ++i)
    {
        double abs_weight = std::abs(bit_weights[i]);
        if (i < 32)
        {
            if (abs_weight > first_bit.first)
            {
                std::swap(first_bit, second_bit);
                first_bit.first = abs_weight;
                first_bit.second = i;
            }
            else if (abs_weight > second_bit.first)
            {
                second_bit.first = abs_weight;
                second_bit.second = i;
            }
        }
        if (bit_weights[i] >= 0)
        {
            simhash_value |= (1ULL << i);
        }
    }

    std::vector<uint64_t> bit_flip_indices;
    uint64_t upper_32_bits = (simhash_value >> 32) & 0xFFFF;

    // 生成包含单位翻转和双位翻转的组合
    bit_flip_indices.push_back(upper_32_bits ^ (1ULL << first_bit.second));
    bit_flip_indices.push_back(upper_32_bits ^ (1ULL << second_bit.second));
    bit_flip_indices.push_back((upper_32_bits ^ (1ULL << first_bit.second)) ^ (1ULL << second_bit.second));
    bit_flip_indices.push_back(upper_32_bits);

    return {bit_flip_indices, simhash_value};
}

// 计算两个64位整数的汉明距离
int calculate_hamming_distance(uint64_t x, uint64_t y)
{
    uint64_t z = x ^ y;
    int distance = 0;
    while (z)
    {
        distance += z & 1;
        z >>= 1;
    }
    return distance;
}

// 在数据库中查找近似重复的文件，并将结果写入输出文件
void find_near_duplicates(std::unordered_map<uint64_t, std::map<uint64_t, std::vector<std::string>>> &database, std::pair<std::vector<uint64_t>, uint64_t> &simhash_data, std::string &filename, int max_hamming_distance, const std::string &output_filename)
{
    uint64_t simhash_value = simhash_data.second;
    uint64_t upper_32_bits = (simhash_value >> 32) & 0xFFFF;

    if (database.size() == 0)
    {
        database[upper_32_bits][simhash_value].push_back(filename);
        return;
    }
    std::vector<std::string> duplicate_filenames;
    auto bit_flip_indices = simhash_data.first;

    // 查找位翻转组合对应的数据库条目，计算汉明距离以查找近似重复
    for (const auto &bit_flip_index : bit_flip_indices)
    {
        auto it = database.find(bit_flip_index);
        if (it != database.end())
        {
            for (const auto &entry : it->second)
            {
                if (calculate_hamming_distance(simhash_value, entry.first) <= max_hamming_distance)
                {
                    duplicate_filenames.insert(duplicate_filenames.end(), entry.second.begin(), entry.second.end());
                }
            }
        }
    }
    database[upper_32_bits][simhash_value].push_back(filename);

    // 将近似重复的文件对写入输出文件
    std::ofstream output(output_filename, std::ios::app);
    if (!output.is_open())
    {
        std::cerr << "Could not open output file " << output_filename << "\n";
        return;
    }
    for (const auto &duplicate_filename : duplicate_filenames)
    {
        output << filename << " " << duplicate_filename << "\n";
    }
    output.close();
}

// 从文件中读取图像特征数据库
void read_feature_database(std::unordered_map<uint64_t, std::map<uint64_t, std::vector<std::string>>> &database, std::string input_filename, const std::string &output_filename, int max_hamming_distance)
{
    std::ifstream input_file(input_filename);
    if (!input_file.is_open())
    {
        std::cerr << "Could not open input file " << input_filename << "\n";
        return;
    }

    int num_features;
    input_file >> num_features;
    std::string line;
    std::getline(input_file, line); // 跳过第一行
    for (int i = 0; i < num_features; ++i)
    {
        std::getline(input_file, line);
        std::istringstream iss(line);
        std::string filename;
        std::vector<double> feature_vector(128);
        iss >> filename;
        for (int j = 0; j < 128; ++j)
        {
            iss >> feature_vector[j];
        }

        auto simhash_data = generate_bit_flips(feature_vector);
        find_near_duplicates(database, simhash_data, filename, max_hamming_distance, output_filename);
    }
    input_file.close();
}

int main(int argc, char *argv[])
{
    if (argc > 50)
    {
        return 1;
    }

    std::unordered_map<uint64_t, std::map<uint64_t, std::vector<std::string>>> database;
    std::string output_filename = argv[argc - 1];
    int max_hamming_distance = 2;

    // 逐个读取输入文件并处理
    for (int i = 1; i < argc - 1; ++i)
    {
        std::string input_filename = argv[i];
        read_feature_database(database, input_filename, output_filename, max_hamming_distance);
    }

    return 0;
}
