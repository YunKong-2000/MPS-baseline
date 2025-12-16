#pragma once
#include <string>
#include "../../third_party/ini/SimpleIni.h"

namespace mps2D {

// 配置参数基类
class ConfigParams {
public:
    ConfigParams() = default;
    virtual ~ConfigParams() = default;
    
    // 从 SimpleIni 加载参数
    virtual bool LoadFromConfig(const SimpleIni& config) = 0;
    
    // 验证参数有效性
    virtual bool Validate() const { return true; }
};

} // namespace mps2D
