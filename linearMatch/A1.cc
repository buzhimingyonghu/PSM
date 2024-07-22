#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <unordered_map>
#include <ctime>
#include <cmath>
#include <sys/time.h>
#include "psm-helper.hpp"

// 计算Simhash的函数
uint64_t compute_simhash(const std::vector<double> &feature)
{
    int f = 64;
    std::vector<double> t(f, 0);

    for (int i = 0; i < 128; ++i)
    {
        uint64_t hash = PSMHelper::MurmurHash64A(std::to_string(i));
        for (int j = 0; j < f; ++j)
        {
            if (hash & (1ULL << j))
            {
                t[j] += feature[i];
            }
            else
            {
                t[j] -= feature[i];
            }
        }
    }

    uint64_t simhash = 0;
    for (int i = 0; i < f; ++i)
    {
        if (t[i] >= 0)
        {
            simhash |= (1ULL << i);
        }
    }
    return simhash;
}

// 计算两个64位整数的汉明距离
int hamming_distance(uint64_t x, uint64_t y)
{
    uint64_t z = x ^ y;
    int dist = 0;
    while (z)
    {
        dist += z & 1;
        z >>= 1;
    }
    return dist;
}

// 从文件中读取图像特征数据库
std::unordered_map<uint64_t, std::vector<std::string>> read_database(const std::string &filename)
{
    std::ifstream file(filename);
    int n;
    file >> n;

    std::unordered_map<uint64_t, std::vector<std::string>> database;
    std::string line;
    std::getline(file, line); // 跳过第一行
    for (int i = 0; i < n; ++i)
    {
        std::getline(file, line);
        std::istringstream iss(line);
        std::string filename;
        std::vector<double> feature(128);
        iss >> filename;
        for (int j = 0; j < 128; ++j)
        {
            iss >> feature[j];
        }
        uint64_t simhash = compute_simhash(feature);

        database[simhash].push_back(filename);
    }
    return database;
}

// 从文件中读取查询特征
std::vector<std::vector<double>> read_queries(const std::string &filename)
{
    std::ifstream file(filename);
    int m;
    file >> m;

    std::vector<std::vector<double>> queries(m, std::vector<double>(128));
    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < 128; ++j)
        {
            file >> queries[i][j];
        }
    }
    return queries;
}
// 查找近似重复图像
std::vector<std::vector<std::string>> find_near_duplicates(
    const std::unordered_map<uint64_t, std::vector<std::string>> &database,
    const std::vector<std::vector<double>> &queries,
    int h)
{

    std::vector<std::vector<std::string>> results;

    for (const auto &query : queries)
    {
        uint64_t query_simhash = compute_simhash(query);
        std::vector<std::string> matches;
        for (auto &simhash : database)
        {
            if (hamming_distance(query_simhash, simhash.first) <= h)
            {
                matches.insert(matches.end(), simhash.second.begin(), simhash.second.end());
            }
        }
        results.push_back(matches);
    }
    return results;
}
int main(int argc, char *argv[])
{
    if (argc != 4)
    {
        exit(0);
    }

    std::string database_file = argv[1];
    std::string query_file = argv[2];
    std::string output_file = argv[3];

    auto database = read_database(database_file);

    PSMHelper::print_current_timestamp();

    auto queries = read_queries(query_file);
    int h = 2; // 汉明距离阈值
    auto results = find_near_duplicates(database, queries, h);

    PSMHelper::print_current_timestamp();

    std::ofstream output(output_file);
    for (const auto &matches : results)
    {
        output << matches.size() << "\n";
        for (const auto &filename : matches)
        {
            output << filename << "\n";
        }
    }

    return 0;
}
