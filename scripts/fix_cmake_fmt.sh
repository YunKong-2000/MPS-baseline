#!/bin/bash
# 修复CMake fmt库下载卡住的问题

echo "正在清理build目录中的fmt相关文件..."
cd "$(dirname "$0")/.."

# 清理fmt相关的下载缓存
if [ -d "build/_deps" ]; then
    echo "删除 build/_deps 目录..."
    rm -rf build/_deps
fi

# 清理CMake缓存
if [ -f "build/CMakeCache.txt" ]; then
    echo "删除 CMakeCache.txt..."
    rm -f build/CMakeCache.txt
fi

echo "清理完成！"
echo ""
echo "现在可以重新运行 cmake 配置："
echo "  cd build"
echo "  cmake .."
echo ""
echo "如果网络仍然很慢，可以考虑："
echo "1. 使用系统已安装的fmt库: sudo apt-get install libfmt-dev (Ubuntu/Debian)"
echo "2. 或者等待网络连接稳定后重试"
