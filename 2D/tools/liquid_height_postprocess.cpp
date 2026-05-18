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

struct Options {
  std::vector<std::filesystem::path> vtk_inputs;
  std::vector<double> x_probes;
  std::filesystem::path output_csv = "liquid_height.csv";
  double dt = 1.0;
  double x_tolerance = 1.0e-6;
};

struct FrameData {
  std::vector<double> x_coords;
  std::vector<double> y_coords;
};

bool StartsWith(const std::string& text, const std::string& prefix) {
  return text.size() >= prefix.size() && text.compare(0, prefix.size(), prefix) == 0;
}

bool IsResultVTKFile(const std::filesystem::path& file_path) {
  std::string extension = file_path.extension().string();
  std::transform(extension.begin(), extension.end(), extension.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  if (extension != ".vtk") {
    return false;
  }
  const std::string stem = file_path.stem().string();
  return StartsWith(stem, "result");
}

bool ParseDouble(const std::string& text, double& value) {
  try {
    size_t pos = 0;
    value = std::stod(text, &pos);
    return pos == text.size();
  } catch (...) {
    return false;
  }
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

std::vector<double> ParseNumbersInLine(const std::string& line) {
  std::vector<double> values;
  std::istringstream iss(line);
  double value = 0.0;
  while (iss >> value) {
    values.push_back(value);
  }
  return values;
}

void PrintUsage(const char* exe_name) {
  std::cout << "用法:\n"
            << "  " << exe_name
            << " --vtk <vtk文件或目录> [--vtk <vtk文件或目录> ...]\n"
            << "           --x <x坐标> [--x <x坐标> ...]\n"
            << "           [--x-file <x坐标文件>] [--dt <输出间隔>] [--x-tolerance <容差>]\n"
            << "           [--output <输出csv路径>]\n\n"
            << "参数说明:\n"
            << "  --vtk          输入vtk文件或目录，可重复\n"
            << "  --x            单个观测点x坐标，可重复\n"
            << "  --x-file       x坐标文件，每行一个x值（支持#注释）\n"
            << "  --dt           相邻输出帧对应时间间隔，默认1.0\n"
            << "  --x-tolerance  观测点x坐标匹配容差 |x_i - x_probe| <= tol，默认1e-6\n"
            << "  --output       输出csv路径，默认 liquid_height.csv\n\n"
            << "示例:\n"
            << "  " << exe_name
            << " --vtk ./output --x 0.5 --x 0.7 --dt 0.001 "
               "--x-tolerance 0.01 --output ./output/liquid_height.csv\n";
}

bool LoadXProbesFromFile(const std::filesystem::path& file_path,
                         std::vector<double>& x_probes,
                         std::string& error_message) {
  std::ifstream file(file_path);
  if (!file.is_open()) {
    error_message = "无法打开x坐标文件: " + file_path.string();
    return false;
  }

  std::string line;
  int line_number = 0;
  while (std::getline(file, line)) {
    ++line_number;
    if (line.empty() || line[0] == '#') {
      continue;
    }
    std::istringstream iss(line);
    double x = 0.0;
    if (!(iss >> x)) {
      error_message = "x坐标文件格式错误，行号 " + std::to_string(line_number);
      return false;
    }
    x_probes.push_back(x);
  }
  return true;
}

bool ExpandVTKInput(const std::filesystem::path& input_path,
                    std::vector<std::filesystem::path>& vtk_files,
                    std::string& error_message) {
  if (!std::filesystem::exists(input_path)) {
    error_message = "VTK输入路径不存在: " + input_path.string();
    return false;
  }

  if (std::filesystem::is_regular_file(input_path)) {
    if (!IsResultVTKFile(input_path)) {
      error_message =
          "输入文件必须是以 result 开头的 vtk 文件: " + input_path.string();
      return false;
    }
    vtk_files.push_back(input_path);
    return true;
  }

  if (std::filesystem::is_directory(input_path)) {
    std::vector<std::filesystem::path> discovered_files;
    for (const auto& entry : std::filesystem::directory_iterator(input_path)) {
      if (!entry.is_regular_file()) {
        continue;
      }
      if (IsResultVTKFile(entry.path())) {
        discovered_files.push_back(entry.path());
      }
    }
    std::sort(discovered_files.begin(), discovered_files.end(),
              [](const std::filesystem::path& a, const std::filesystem::path& b) {
                return a.filename().string() < b.filename().string();
              });
    vtk_files.insert(vtk_files.end(), discovered_files.begin(), discovered_files.end());
    return true;
  }

  error_message = "不支持的VTK输入路径类型: " + input_path.string();
  return false;
}

bool ParseOptions(int argc, char* argv[], Options& options) {
  if (argc <= 1) {
    PrintUsage(argv[0]);
    return false;
  }

  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--help" || arg == "-h") {
      PrintUsage(argv[0]);
      return false;
    }
    if (arg == "--vtk") {
      if (i + 1 >= argc) {
        std::cerr << "错误：--vtk 缺少参数\n";
        return false;
      }
      options.vtk_inputs.emplace_back(argv[++i]);
      continue;
    }
    if (arg == "--x") {
      if (i + 1 >= argc) {
        std::cerr << "错误：--x 缺少参数\n";
        return false;
      }
      double x_value = 0.0;
      if (!ParseDouble(argv[++i], x_value)) {
        std::cerr << "错误：--x 参数必须是数字\n";
        return false;
      }
      options.x_probes.push_back(x_value);
      continue;
    }
    if (arg == "--x-file") {
      if (i + 1 >= argc) {
        std::cerr << "错误：--x-file 缺少参数\n";
        return false;
      }
      std::string error_message;
      if (!LoadXProbesFromFile(argv[++i], options.x_probes, error_message)) {
        std::cerr << "错误：" << error_message << "\n";
        return false;
      }
      continue;
    }
    if (arg == "--dt") {
      if (i + 1 >= argc) {
        std::cerr << "错误：--dt 缺少参数\n";
        return false;
      }
      if (!ParseDouble(argv[++i], options.dt) || options.dt <= 0.0) {
        std::cerr << "错误：--dt 必须是正数\n";
        return false;
      }
      continue;
    }
    if (arg == "--x-tolerance") {
      if (i + 1 >= argc) {
        std::cerr << "错误：--x-tolerance 缺少参数\n";
        return false;
      }
      if (!ParseDouble(argv[++i], options.x_tolerance) || options.x_tolerance < 0.0) {
        std::cerr << "错误：--x-tolerance 必须是非负数\n";
        return false;
      }
      continue;
    }
    if (arg == "--output") {
      if (i + 1 >= argc) {
        std::cerr << "错误：--output 缺少参数\n";
        return false;
      }
      options.output_csv = argv[++i];
      continue;
    }

    std::cerr << "错误：未知参数 " << arg << "\n";
    PrintUsage(argv[0]);
    return false;
  }

  if (options.vtk_inputs.empty()) {
    std::cerr << "错误：至少需要一个VTK输入（--vtk）\n";
    return false;
  }
  if (options.x_probes.empty()) {
    std::cerr << "错误：至少需要一个观测点x坐标（--x 或 --x-file）\n";
    return false;
  }
  return true;
}

bool ReadPointCoordinates(std::ifstream& file, int num_points, FrameData& frame) {
  std::vector<double> buffer;
  buffer.reserve(static_cast<size_t>(num_points) * 3U);

  std::string line;
  while (static_cast<int>(buffer.size()) < num_points * 3 && std::getline(file, line)) {
    const std::vector<double> values = ParseNumbersInLine(line);
    buffer.insert(buffer.end(), values.begin(), values.end());
  }
  if (static_cast<int>(buffer.size()) < num_points * 3) {
    return false;
  }

  frame.x_coords.clear();
  frame.y_coords.clear();
  frame.x_coords.reserve(static_cast<size_t>(num_points));
  frame.y_coords.reserve(static_cast<size_t>(num_points));
  for (int i = 0; i < num_points; ++i) {
    frame.x_coords.push_back(buffer[static_cast<size_t>(3 * i)]);
    frame.y_coords.push_back(buffer[static_cast<size_t>(3 * i + 1)]);
  }
  return true;
}

bool ReadVTKFramePositions(const std::filesystem::path& vtk_path,
                           FrameData& frame,
                           std::string& error_message) {
  std::ifstream file(vtk_path);
  if (!file.is_open()) {
    error_message = "无法打开VTK文件";
    return false;
  }

  std::string line;
  int num_points = 0;
  bool found_points = false;

  while (std::getline(file, line)) {
    if (!StartsWith(line, "POINTS ")) {
      continue;
    }
    std::istringstream iss(line);
    std::string keyword;
    std::string vtk_type;
    if (!(iss >> keyword >> num_points >> vtk_type) || num_points <= 0) {
      error_message = "POINTS 行格式错误";
      return false;
    }
    found_points = true;
    break;
  }

  if (!found_points) {
    error_message = "未找到 POINTS 段";
    return false;
  }

  if (!ReadPointCoordinates(file, num_points, frame)) {
    error_message = "POINTS 数据不完整";
    return false;
  }
  return true;
}

double ComputeHeightAtProbeX(const FrameData& frame, double probe_x, double x_tolerance) {
  if (frame.x_coords.empty()) {
    return 0.0;
  }

  const double tolerance = x_tolerance + 1.0e-15;
  double max_y = -std::numeric_limits<double>::infinity();
  bool found = false;

  for (size_t i = 0; i < frame.x_coords.size(); ++i) {
    if (std::fabs(frame.x_coords[i] - probe_x) <= tolerance) {
      if (!found || frame.y_coords[i] > max_y) {
        max_y = frame.y_coords[i];
      }
      found = true;
    }
  }

  if (!found) {
    return 0.0;
  }
  return std::max(0.0, max_y);
}

}  // namespace

int main(int argc, char* argv[]) {
  Options options;
  if (!ParseOptions(argc, argv, options)) {
    return 1;
  }

  std::vector<std::filesystem::path> vtk_files;
  std::string error_message;
  for (const auto& input : options.vtk_inputs) {
    if (!ExpandVTKInput(input, vtk_files, error_message)) {
      std::cerr << "错误：" << error_message << "\n";
      return 1;
    }
  }

  if (vtk_files.empty()) {
    std::cerr << "错误：未发现任何VTK文件\n";
    return 1;
  }

  std::sort(vtk_files.begin(), vtk_files.end(),
            [](const std::filesystem::path& a, const std::filesystem::path& b) {
              return a.string() < b.string();
            });
  vtk_files.erase(std::unique(vtk_files.begin(), vtk_files.end()), vtk_files.end());

  std::ofstream csv_file(options.output_csv, std::ios::trunc);
  if (!csv_file.is_open()) {
    std::cerr << "错误：无法创建输出CSV文件: " << options.output_csv << "\n";
    return 1;
  }

  csv_file << std::fixed << std::setprecision(15);
  csv_file << "time";
  for (size_t i = 0; i < options.x_probes.size(); ++i) {
    csv_file << ",probe_" << i << "_height";
  }
  csv_file << "\n";

  int processed = 0;
  int skipped = 0;
  for (size_t order = 0; order < vtk_files.size(); ++order) {
    const std::filesystem::path& vtk_path = vtk_files[order];
    const std::optional<int> step_id_opt =
        ExtractTrailingInteger(vtk_path.stem().string());
    const int step_id = step_id_opt.has_value() ? *step_id_opt : static_cast<int>(order);
    const double time = static_cast<double>(step_id) * options.dt;

    FrameData frame;
    if (!ReadVTKFramePositions(vtk_path, frame, error_message)) {
      ++skipped;
      std::cerr << "跳过 " << vtk_path << "，原因: " << error_message << "\n";
      continue;
    }

    csv_file << time;
    for (double probe_x : options.x_probes) {
      const double h = ComputeHeightAtProbeX(frame, probe_x, options.x_tolerance);
      csv_file << "," << h;
    }
    csv_file << "\n";
    ++processed;
  }

  std::cout << "后处理完成\n";
  std::cout << "  输出文件: " << options.output_csv << "\n";
  std::cout << "  VTK文件总数: " << vtk_files.size() << "\n";
  std::cout << "  成功处理: " << processed << "\n";
  std::cout << "  跳过文件: " << skipped << "\n";
  std::cout << "  观测点数量: " << options.x_probes.size() << "\n";
  std::cout << "  时间间隔 dt: " << options.dt << "\n";
  std::cout << "  x匹配容差: " << options.x_tolerance << "\n";

  return (processed > 0) ? 0 : 1;
}
