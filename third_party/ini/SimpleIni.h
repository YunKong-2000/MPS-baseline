/*
 * SimpleIni - A cross-platform C++ library for reading and writing INI files
 * This is a simplified version for lightweight usage
 */

#ifndef SIMPLE_INI_H
#define SIMPLE_INI_H

#include <string>
#include <map>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>

class SimpleIni {
public:
    SimpleIni() = default;
    ~SimpleIni() = default;

    // 加载INI文件
    bool LoadFile(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            return false;
        }

        std::string line;
        std::string current_section;
        
        while (std::getline(file, line)) {
            // 移除行首行尾空白
            trim(line);
            
            // 跳过空行和注释行
            if (line.empty() || line[0] == ';' || line[0] == '#') {
                continue;
            }
            
            // 检查是否是节（section）
            if (line[0] == '[' && line.back() == ']') {
                current_section = line.substr(1, line.length() - 2);
                trim(current_section);
            } else {
                // 解析键值对
                size_t pos = line.find('=');
                if (pos != std::string::npos) {
                    std::string key = line.substr(0, pos);
                    std::string value = line.substr(pos + 1);
                    trim(key);
                    trim(value);
                    
                    // 移除值的引号（如果存在）
                    if (value.length() >= 2 && 
                        ((value[0] == '"' && value.back() == '"') ||
                         (value[0] == '\'' && value.back() == '\''))) {
                        value = value.substr(1, value.length() - 2);
                    }
                    
                    std::string full_key = current_section.empty() ? 
                                          key : current_section + "." + key;
                    data_[full_key] = value;
                }
            }
        }
        
        file.close();
        return true;
    }

    // 获取字符串值
    std::string GetValue(const std::string& section, const std::string& key, 
                        const std::string& default_value = "") const {
        std::string full_key = section.empty() ? key : section + "." + key;
        auto it = data_.find(full_key);
        return (it != data_.end()) ? it->second : default_value;
    }

    // 获取整数值
    int GetIntValue(const std::string& section, const std::string& key, 
                   int default_value = 0) const {
        std::string value = GetValue(section, key, "");
        if (value.empty()) {
            return default_value;
        }
        try {
            return std::stoi(value);
        } catch (...) {
            return default_value;
        }
    }

    // 获取浮点数值
    double GetDoubleValue(const std::string& section, const std::string& key, 
                         double default_value = 0.0) const {
        std::string value = GetValue(section, key, "");
        if (value.empty()) {
            return default_value;
        }
        try {
            return std::stod(value);
        } catch (...) {
            return default_value;
        }
    }

    // 获取布尔值
    bool GetBoolValue(const std::string& section, const std::string& key, 
                     bool default_value = false) const {
        std::string value = GetValue(section, key, "");
        if (value.empty()) {
            return default_value;
        }
        std::transform(value.begin(), value.end(), value.begin(), ::tolower);
        return (value == "true" || value == "1" || value == "yes" || value == "on");
    }

    // 检查键是否存在
    bool HasKey(const std::string& section, const std::string& key) const {
        std::string full_key = section.empty() ? key : section + "." + key;
        return data_.find(full_key) != data_.end();
    }

    // 清空所有数据
    void Clear() {
        data_.clear();
    }

private:
    std::map<std::string, std::string> data_;

    // 移除字符串首尾空白
    static void trim(std::string& str) {
        str.erase(0, str.find_first_not_of(" \t\r\n"));
        str.erase(str.find_last_not_of(" \t\r\n") + 1);
    }
};

#endif // SIMPLE_INI_H

