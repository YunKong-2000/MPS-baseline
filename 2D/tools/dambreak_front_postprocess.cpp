#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

namespace {

constexpr int kSurfaceTypeSurface = 2;
constexpr int kSurfaceTypeSplash = 3;

struct FrontPoint {
  bool found = false;
  double x = 0.0;
  double y = 0.0;
  int particle_index = -1;
  int surface_type = -1;
};

bool ParsePointCountLine(const std::string& line, int& num_points) {
  std::istringstream iss(line);
  std::string keyword;
  std::string data_type;
  if (!(iss >> keyword >> num_points >> data_type)) {
    return false;
  }
  return keyword == "POINTS" && num_points > 0;
}

bool ParsePointLine(const std::string& line, double& x, double& y, double& z) {
  std::istringstream iss(line);
  return static_cast<bool>(iss >> x >> y >> z);
}

std::optional<int> ExtractTrailingInteger(const std::string& text) {
  int end = static_cast<int>(text.size()) - 1;
  while (end >= 0 && !std::isdigit(static_cast<unsigned char>(text[end]))) {
    --end;
  }
  if (end < 0) {
    return std::nullopt;
  }

  int begin = end;
  while (begin >= 0 && std::isdigit(static_cast<unsigned char>(text[begin]))) {
    --begin;
  }
  ++begin;

  try {
    return std::stoi(text.substr(static_cast<size_t>(begin),
                                 static_cast<size_t>(end - begin + 1)));
  } catch (...) {
    return std::nullopt;
  }
}

bool ExtractFrontFromVTK(const std::filesystem::path& vtk_path,
                         FrontPoint& front_point,
                         std::string& error_message) {
  std::ifstream vtk_file(vtk_path);
  if (!vtk_file.is_open()) {
    error_message = "无法打开文件";
    return false;
  }

  front_point = FrontPoint{};
  std::string line;
  int num_points = 0;
  bool points_found = false;

  while (std::getline(vtk_file, line)) {
    if (line.rfind("POINTS ", 0) == 0) {
      if (!ParsePointCountLine(line, num_points)) {
        error_message = "POINTS 行格式无效";
        return false;
      }
      points_found = true;
      break;
    }
  }
  if (!points_found) {
    error_message = "未找到 POINTS 段";
    return false;
  }

  std::vector<double> x_coords(num_points, 0.0);
  std::vector<double> y_coords(num_points, 0.0);
  for (int i = 0; i < num_points; ++i) {
    if (!std::getline(vtk_file, line)) {
      error_message = "POINTS 数据不足";
      return false;
    }
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    if (!ParsePointLine(line, x, y, z)) {
      error_message = "POINTS 数据格式无效";
      return false;
    }
    x_coords[i] = x;
    y_coords[i] = y;
  }

  bool surface_type_found = false;
  while (std::getline(vtk_file, line)) {
    if (line.rfind("SCALARS surface_type ", 0) == 0) {
      surface_type_found = true;
      break;
    }
  }
  if (!surface_type_found) {
    error_message = "未找到 surface_type 标量场";
    return false;
  }

  if (!std::getline(vtk_file, line) || line.rfind("LOOKUP_TABLE", 0) != 0) {
    error_message = "surface_type LOOKUP_TABLE 缺失";
    return false;
  }

  std::vector<int> surface_type_values(num_points, 0);
  for (int i = 0; i < num_points; ++i) {
    if (!std::getline(vtk_file, line)) {
      error_message = "surface_type 数据不足";
      return false;
    }
    std::istringstream iss(line);
    double raw_value = 0.0;
    if (!(iss >> raw_value)) {
      error_message = "surface_type 数据格式无效";
      return false;
    }
    surface_type_values[i] = static_cast<int>(std::lround(raw_value));
  }

  double max_x = -std::numeric_limits<double>::infinity();
  for (int i = 0; i < num_points; ++i) {
    const int st = surface_type_values[i];
    if (st != kSurfaceTypeSplash && st != kSurfaceTypeSurface) {
      continue;
    }
    if (x_coords[i] > max_x) {
      max_x = x_coords[i];
      front_point.found = true;
      front_point.x = x_coords[i];
      front_point.y = y_coords[i];
      front_point.particle_index = i;
      front_point.surface_type = st;
    }
  }

  return true;
}

void PrintUsage(const char* exe_name) {
  std::cout << "用法: " << exe_name
            << " <vtk目录> [输出csv文件] [时间步大小dt]\n"
            << "示例:\n"
            << "  " << exe_name << " ./output\n"
            << "  " << exe_name
            << " ./output ./output/dambreak_front_position.csv 0.001\n";
}

}  // namespace

int main(int argc, char* argv[]) {
  if (argc < 2) {
    PrintUsage(argv[0]);
    return 1;
  }

  const std::filesystem::path vtk_dir(argv[1]);
  if (!std::filesystem::exists(vtk_dir) || !std::filesystem::is_directory(vtk_dir)) {
    std::cerr << "错误：路径不是有效目录: " << vtk_dir << std::endl;
    return 1;
  }

  std::filesystem::path output_csv;
  if (argc >= 3) {
    output_csv = std::filesystem::path(argv[2]);
  } else {
    output_csv = vtk_dir / "dambreak_front_position.csv";
  }

  double dt = 1.0;
  if (argc >= 4) {
    try {
      dt = std::stod(argv[3]);
    } catch (...) {
      std::cerr << "错误：时间步大小 dt 不是有效数字: " << argv[3] << std::endl;
      return 1;
    }
    if (dt <= 0.0) {
      std::cerr << "错误：时间步大小 dt 必须大于0: " << argv[3] << std::endl;
      return 1;
    }
  }

  std::vector<std::filesystem::path> vtk_files;
  for (const auto& entry : std::filesystem::directory_iterator(vtk_dir)) {
    if (!entry.is_regular_file()) {
      continue;
    }
    const auto& path = entry.path();
    std::string extension = path.extension().string();
    std::transform(extension.begin(), extension.end(), extension.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (extension == ".vtk") {
      vtk_files.push_back(path);
    }
  }

  if (vtk_files.empty()) {
    std::cerr << "错误：目录中未找到 vtk 文件: " << vtk_dir << std::endl;
    return 1;
  }

  std::sort(vtk_files.begin(), vtk_files.end(),
            [](const std::filesystem::path& a, const std::filesystem::path& b) {
              return a.filename().string() < b.filename().string();
            });

  std::ofstream out_file(output_csv, std::ios::trunc);
  if (!out_file.is_open()) {
    std::cerr << "错误：无法创建输出文件: " << output_csv << std::endl;
    return 1;
  }

  out_file << std::fixed << std::setprecision(15);
  out_file << "time,front_x\n";

  int processed_count = 0;
  int skipped_count = 0;
  for (size_t order = 0; order < vtk_files.size(); ++order) {
    const auto& vtk_path = vtk_files[order];
    const std::string stem = vtk_path.stem().string();
    const std::optional<int> parsed_index = ExtractTrailingInteger(stem);
    const int step_id = parsed_index.has_value() ? *parsed_index
                                                 : static_cast<int>(order);
    const double time = static_cast<double>(step_id) * dt;

    FrontPoint front_point;
    std::string error_message;
    if (!ExtractFrontFromVTK(vtk_path, front_point, error_message)) {
      ++skipped_count;
      std::cerr << "跳过文件: " << vtk_path.filename().string()
                << "，原因: " << error_message << std::endl;
      continue;
    }

    out_file << time << ",";
    if (front_point.found) {
      out_file << front_point.x << "\n";
    } else {
      out_file << "nan\n";
    }
    ++processed_count;
  }

  std::cout << "后处理完成:\n";
  std::cout << "  输入目录: " << vtk_dir << "\n";
  std::cout << "  输出文件: " << output_csv << "\n";
  std::cout << "  处理文件数: " << processed_count << "\n";
  std::cout << "  跳过文件数: " << skipped_count << "\n";
  std::cout << "  时间步大小 dt: " << dt << "\n";

  return 0;
}
