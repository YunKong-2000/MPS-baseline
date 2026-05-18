#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace {

constexpr double kTinyDistance = 1.0e-12;
constexpr int kDimension2D = 2;

struct Point2D {
  double x = 0.0;
  double y = 0.0;
};

struct Options {
  enum class ProbeMethod {
    kAverageInRadius = 0,
    kNearestParticle = 1,
    kOriginalMpsTaylor = 2,
  };

  double radius = 0.02;
  int mps_probe_min_neighbors = 3;
  std::string pressure_field = "pressure";
  std::filesystem::path output_csv = "pressure_probe.csv";
  ProbeMethod method = ProbeMethod::kOriginalMpsTaylor;
  std::vector<Point2D> probe_points;
  std::vector<std::filesystem::path> vtk_inputs;
};

struct VTKFrame {
  std::vector<Point2D> positions;
  std::vector<double> pressures;
};

void PrintUsage(const char* exe_name) {
  std::cout << "用法:\n"
            << "  " << exe_name
            << " --radius <半径> --point <x> <y> [--point <x> <y> ...]\n"
            << "  " << std::string(exe_name)
            << " [--points-file <点文件>] --vtk <vtk文件或目录>"
            << " [--vtk <vtk文件或目录> ...]\n"
            << "  " << std::string(exe_name)
            << " [--pressure-field <字段名>] [--method <average|nearest|mps_taylor>]"
               " [--output <输出csv>]\n\n"
            << "参数说明:\n"
            << "  --radius          邻域半径（average/nearest/mps_taylor均使用）\n"
            << "  --mps-probe-min-neighbors  mps_taylor 方法下测压点最少邻域粒子数，默认 3\n"
            << "  --point x y       单个观测点，可重复多次\n"
            << "  --points-file     点文件，每行格式: x y（支持#注释）\n"
            << "  --vtk             输入vtk文件或包含vtk文件的目录，可重复\n"
            << "  --pressure-field  压力标量字段名，默认 pressure\n"
            << "  --method          压力测量方法：average、nearest 或 mps_taylor，默认 mps_taylor\n"
            << "                    兼容别名：lsmps（会映射为 mps_taylor）\n"
            << "  --output          输出csv路径，默认 pressure_probe.csv\n\n"
            << "示例:\n"
            << "  " << exe_name
            << " --radius 0.03 --point 0.5 0.25 --point 0.7 0.25 "
               "--vtk ./output --output ./output/probe_pressure.csv\n";
}

bool StartsWith(const std::string& text, const std::string& prefix) {
  if (text.size() < prefix.size()) {
    return false;
  }
  return text.compare(0, prefix.size(), prefix) == 0;
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

std::vector<double> ParseNumbersInLine(const std::string& line) {
  std::vector<double> values;
  std::istringstream iss(line);
  double value = 0.0;
  while (iss >> value) {
    values.push_back(value);
  }
  return values;
}

bool LoadPointsFromFile(const std::filesystem::path& file_path,
                        std::vector<Point2D>& points,
                        std::string& error_message) {
  std::ifstream file(file_path);
  if (!file.is_open()) {
    error_message = "无法打开点文件: " + file_path.string();
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
    Point2D point;
    if (!(iss >> point.x >> point.y)) {
      error_message = "点文件格式错误，行号 " + std::to_string(line_number);
      return false;
    }
    points.push_back(point);
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
    std::string extension = input_path.extension().string();
    std::transform(extension.begin(), extension.end(), extension.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (extension != ".vtk") {
      error_message = "输入文件不是vtk格式: " + input_path.string();
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
      std::string extension = entry.path().extension().string();
      std::transform(extension.begin(), extension.end(), extension.begin(),
                     [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
      if (extension == ".vtk") {
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
    if (arg == "--radius") {
      if (i + 1 >= argc) {
        std::cerr << "错误：--radius 缺少参数\n";
        return false;
      }
      if (!ParseDouble(argv[++i], options.radius) || options.radius <= 0.0) {
        std::cerr << "错误：--radius 必须是正数\n";
        return false;
      }
      continue;
    }
    if (arg == "--mps-probe-min-neighbors") {
      if (i + 1 >= argc) {
        std::cerr << "错误：--mps-probe-min-neighbors 缺少参数\n";
        return false;
      }
      double value = 0.0;
      if (!ParseDouble(argv[++i], value) || value < 1.0) {
        std::cerr << "错误：--mps-probe-min-neighbors 必须是 >=1 的整数\n";
        return false;
      }
      const int int_value = static_cast<int>(value);
      if (std::fabs(value - static_cast<double>(int_value)) > 1.0e-9) {
        std::cerr << "错误：--mps-probe-min-neighbors 必须是整数\n";
        return false;
      }
      options.mps_probe_min_neighbors = int_value;
      continue;
    }
    if (arg == "--pressure-field") {
      if (i + 1 >= argc) {
        std::cerr << "错误：--pressure-field 缺少参数\n";
        return false;
      }
      options.pressure_field = argv[++i];
      continue;
    }
    if (arg == "--method") {
      if (i + 1 >= argc) {
        std::cerr << "错误：--method 缺少参数\n";
        return false;
      }
      const std::string method = argv[++i];
      if (method == "average") {
        options.method = Options::ProbeMethod::kAverageInRadius;
      } else if (method == "nearest") {
        options.method = Options::ProbeMethod::kNearestParticle;
      } else if (method == "mps_taylor" || method == "lsmps") {
        options.method = Options::ProbeMethod::kOriginalMpsTaylor;
      } else {
        std::cerr << "错误：--method 仅支持 average、nearest 或 mps_taylor\n";
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
    if (arg == "--point") {
      if (i + 2 >= argc) {
        std::cerr << "错误：--point 需要两个参数 x y\n";
        return false;
      }
      Point2D point;
      if (!ParseDouble(argv[++i], point.x) || !ParseDouble(argv[++i], point.y)) {
        std::cerr << "错误：--point 参数必须是数字\n";
        return false;
      }
      options.probe_points.push_back(point);
      continue;
    }
    if (arg == "--points-file") {
      if (i + 1 >= argc) {
        std::cerr << "错误：--points-file 缺少参数\n";
        return false;
      }
      std::string error_message;
      if (!LoadPointsFromFile(argv[++i], options.probe_points, error_message)) {
        std::cerr << "错误：" << error_message << "\n";
        return false;
      }
      continue;
    }
    if (arg == "--vtk") {
      if (i + 1 >= argc) {
        std::cerr << "错误：--vtk 缺少参数\n";
        return false;
      }
      options.vtk_inputs.emplace_back(argv[++i]);
      continue;
    }

    std::cerr << "错误：未知参数 " << arg << "\n";
    PrintUsage(argv[0]);
    return false;
  }

  if (options.probe_points.empty()) {
    std::cerr << "错误：至少需要一个观测点（--point 或 --points-file）\n";
    return false;
  }
  if (options.vtk_inputs.empty()) {
    std::cerr << "错误：至少需要一个VTK输入（--vtk）\n";
    return false;
  }

  return true;
}

bool ReadPointCoordinates(std::ifstream& file, int num_points, std::vector<Point2D>& points) {
  points.clear();
  points.reserve(static_cast<size_t>(num_points));

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

  for (int i = 0; i < num_points; ++i) {
    points.push_back(Point2D{
        buffer[static_cast<size_t>(3 * i)],
        buffer[static_cast<size_t>(3 * i + 1)]});
  }
  return true;
}

bool ReadScalarBlock(std::ifstream& file, int value_count, std::vector<double>& values) {
  values.clear();
  values.reserve(static_cast<size_t>(value_count));

  std::string line;
  while (static_cast<int>(values.size()) < value_count && std::getline(file, line)) {
    const std::vector<double> line_values = ParseNumbersInLine(line);
    values.insert(values.end(), line_values.begin(), line_values.end());
  }

  return static_cast<int>(values.size()) >= value_count;
}

bool ReadVTKFrame(const std::filesystem::path& vtk_path,
                  const std::string& pressure_field_name,
                  VTKFrame& frame,
                  std::string& error_message) {
  std::ifstream file(vtk_path);
  if (!file.is_open()) {
    error_message = "无法打开VTK文件";
    return false;
  }

  frame = VTKFrame{};
  std::string line;
  int num_points = 0;
  bool found_points = false;

  while (std::getline(file, line)) {
    if (StartsWith(line, "POINTS ")) {
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
  }
  if (!found_points) {
    error_message = "未找到 POINTS 段";
    return false;
  }

  if (!ReadPointCoordinates(file, num_points, frame.positions)) {
    error_message = "POINTS 数据不完整";
    return false;
  }

  int point_data_count = num_points;
  bool found_point_data = false;
  while (std::getline(file, line)) {
    if (StartsWith(line, "POINT_DATA ")) {
      std::istringstream iss(line);
      std::string keyword;
      if (!(iss >> keyword >> point_data_count) || point_data_count <= 0) {
        error_message = "POINT_DATA 行格式错误";
        return false;
      }
      found_point_data = true;
      break;
    }
  }
  if (!found_point_data) {
    error_message = "未找到 POINT_DATA 段";
    return false;
  }

  bool found_pressure_field = false;
  while (std::getline(file, line)) {
    if (!StartsWith(line, "SCALARS ")) {
      continue;
    }

    std::istringstream iss(line);
    std::string token;
    std::string scalar_name;
    std::string scalar_type;
    if (!(iss >> token >> scalar_name >> scalar_type)) {
      error_message = "SCALARS 行格式错误";
      return false;
    }

    std::string lookup_line;
    if (!std::getline(file, lookup_line) || !StartsWith(lookup_line, "LOOKUP_TABLE")) {
      error_message = "SCALARS 段缺少 LOOKUP_TABLE";
      return false;
    }

    std::vector<double> scalar_values;
    if (!ReadScalarBlock(file, point_data_count, scalar_values)) {
      error_message = "SCALARS 数据不完整";
      return false;
    }
    if (scalar_name == pressure_field_name) {
      frame.pressures.assign(scalar_values.begin(),
                             scalar_values.begin() + point_data_count);
      found_pressure_field = true;
      continue;
    }
  }

  if (!found_pressure_field) {
    error_message = "未找到压力字段: " + pressure_field_name;
    return false;
  }

  if (frame.positions.size() != frame.pressures.size()) {
    error_message = "位置与压力数据数量不一致";
    return false;
  }
  return true;
}

double ComputeProbePressure(const VTKFrame& frame, const Point2D& probe, double radius) {
  const double radius_sq = radius * radius;
  double sum_pressure = 0.0;
  int count = 0;

  for (size_t i = 0; i < frame.positions.size(); ++i) {
    const double dx = frame.positions[i].x - probe.x;
    const double dy = frame.positions[i].y - probe.y;
    const double dist_sq = dx * dx + dy * dy;
    if (dist_sq <= radius_sq) {
      sum_pressure += frame.pressures[i];
      ++count;
    }
  }

  if (count == 0) {
    return 0.0;
  }
  return sum_pressure / static_cast<double>(count);
}

double ComputeProbePressureByNearest(const VTKFrame& frame,
                                     const Point2D& probe,
                                     double radius) {
  if (frame.positions.empty()) {
    return 0.0;
  }

  const double radius_sq = radius * radius;
  double best_dist_sq = std::numeric_limits<double>::max();
  size_t best_index = 0;
  bool found = false;
  for (size_t i = 0; i < frame.positions.size(); ++i) {
    const double dx = frame.positions[i].x - probe.x;
    const double dy = frame.positions[i].y - probe.y;
    const double dist_sq = dx * dx + dy * dy;
    if (dist_sq > radius_sq) {
      continue;
    }
    if (dist_sq < best_dist_sq) {
      best_dist_sq = dist_sq;
      best_index = i;
      found = true;
    }
  }

  if (!found) {
    return 0.0;
  }
  return frame.pressures[best_index];
}

double ComputeMpsWeight(double distance, double radius) {
  if (distance <= kTinyDistance || distance >= radius) {
    return 0.0;
  }
  return pow(radius / distance - 1.0, 2);
}

bool FindNearestParticleIndex(const VTKFrame& frame,
                              const Point2D& probe,
                              double radius,
                              size_t& index,
                              double& distance_sq,
                              bool only_within_radius) {
  const double radius_sq = radius * radius;
  bool found = false;
  double best_dist_sq = std::numeric_limits<double>::max();
  size_t best_index = 0;
  for (size_t k = 0; k < frame.positions.size(); ++k) {
    const double dx = frame.positions[k].x - probe.x;
    const double dy = frame.positions[k].y - probe.y;
    const double dist_sq = dx * dx + dy * dy;
    if (only_within_radius && dist_sq > radius_sq) {
      continue;
    }
    if (dist_sq < best_dist_sq) {
      best_dist_sq = dist_sq;
      best_index = k;
      found = true;
    }
  }
  if (!found) {
    return false;
  }
  index = best_index;
  distance_sq = best_dist_sq;
  return true;
}

double ComputeProbePressureByOriginalMpsTaylor(const VTKFrame& frame,
                                               const Point2D& probe,
                                               double radius,
                                               int probe_min_neighbors) {
  if (frame.positions.empty() || radius <= 0.0) {
    return 0.0;
  }

  const double radius_sq = radius * radius;
  int probe_neighbor_count = 0;
  for (size_t k = 0; k < frame.positions.size(); ++k) {
    const double dx = frame.positions[k].x - probe.x;
    const double dy = frame.positions[k].y - probe.y;
    const double dist_sq = dx * dx + dy * dy;
    if (dist_sq <= kTinyDistance * kTinyDistance || dist_sq > radius_sq) {
      continue;
    }
    ++probe_neighbor_count;
  }

  // 当测压点邻域粒子数少于3时，不寻找锚点，压力直接置零。
  if (probe_neighbor_count < probe_min_neighbors) {
    return 0.0;
  }

  size_t anchor_index = 0;
  double anchor_dist_sq = 0.0;
  if (!FindNearestParticleIndex(frame, probe, radius, anchor_index, anchor_dist_sq, true)) {
    return 0.0;
  }
  const Point2D& anchor = frame.positions[anchor_index];
  const double anchor_pressure = frame.pressures[anchor_index];

  double n0 = 0.0;
  for (size_t k = 0; k < frame.positions.size(); ++k) {
    const double dx = frame.positions[k].x - anchor.x;
    const double dy = frame.positions[k].y - anchor.y;
    const double dist_sq = dx * dx + dy * dy;
    if (dist_sq <= kTinyDistance * kTinyDistance || dist_sq > radius_sq) {
      continue;
    }
    const double dist = std::sqrt(dist_sq);
    n0 += ComputeMpsWeight(dist, radius);
  }

  if (n0 <= kTinyDistance) {
    return anchor_pressure;
  }

  double grad_x = 0.0;
  double grad_y = 0.0;
  for (size_t k = 0; k < frame.positions.size(); ++k) {
    const double dx = frame.positions[k].x - anchor.x;
    const double dy = frame.positions[k].y - anchor.y;
    const double dist_sq = dx * dx + dy * dy;
    if (dist_sq <= kTinyDistance * kTinyDistance || dist_sq > radius_sq) {
      continue;
    }
    const double dist = std::sqrt(dist_sq);
    const double w = ComputeMpsWeight(dist, radius);
    if (w <= 0.0) {
      continue;
    }

    const double factor = static_cast<double>(kDimension2D) * w *
                          (frame.pressures[k] - anchor_pressure) / (n0 * dist_sq);
    grad_x += factor * dx;
    grad_y += factor * dy;
  }

  const double rx = probe.x - anchor.x;
  const double ry = probe.y - anchor.y;
  return anchor_pressure + grad_x * rx + grad_y * ry;
}

}  // namespace

int main(int argc, char* argv[]) {
  Options options;
  if (!ParseOptions(argc, argv, options)) {
    return 1;
  }

  std::vector<std::filesystem::path> vtk_files;
  std::string error_message;
  for (const auto& input_path : options.vtk_inputs) {
    if (!ExpandVTKInput(input_path, vtk_files, error_message)) {
      std::cerr << "错误：" << error_message << std::endl;
      return 1;
    }
  }

  if (vtk_files.empty()) {
    std::cerr << "错误：未发现任何VTK文件" << std::endl;
    return 1;
  }

  std::sort(vtk_files.begin(), vtk_files.end(),
            [](const std::filesystem::path& a, const std::filesystem::path& b) {
              return a.string() < b.string();
            });
  vtk_files.erase(std::unique(vtk_files.begin(), vtk_files.end()), vtk_files.end());

  std::ofstream csv_file(options.output_csv, std::ios::trunc);
  if (!csv_file.is_open()) {
    std::cerr << "错误：无法创建输出CSV文件: " << options.output_csv << std::endl;
    return 1;
  }

  csv_file << std::fixed << std::setprecision(15);
  csv_file << "step_id";
  for (size_t i = 0; i < options.probe_points.size(); ++i) {
    csv_file << ",probe_" << i << "_pressure";
  }
  csv_file << "\n";

  int processed = 0;
  int skipped = 0;

  for (size_t step = 0; step < vtk_files.size(); ++step) {
    const auto& vtk_path = vtk_files[step];

    VTKFrame frame;
    if (!ReadVTKFrame(vtk_path, options.pressure_field, frame, error_message)) {
      ++skipped;
      std::cerr << "跳过 " << vtk_path << "，原因: " << error_message << std::endl;
      continue;
    }

    csv_file << step;
    for (const auto& probe : options.probe_points) {
      double pressure = 0.0;
      if (options.method == Options::ProbeMethod::kAverageInRadius) {
        pressure = ComputeProbePressure(frame, probe, options.radius);
      } else if (options.method == Options::ProbeMethod::kNearestParticle) {
        pressure = ComputeProbePressureByNearest(frame, probe, options.radius);
      } else {
        pressure = ComputeProbePressureByOriginalMpsTaylor(
            frame, probe, options.radius, options.mps_probe_min_neighbors);
      }
      csv_file << "," << pressure;
    }
    csv_file << "\n";
    ++processed;
  }

  std::cout << "后处理完成\n";
  std::cout << "  输出文件: " << options.output_csv << "\n";
  std::cout << "  VTK文件总数: " << vtk_files.size() << "\n";
  std::cout << "  成功处理: " << processed << "\n";
  std::cout << "  跳过文件: " << skipped << "\n";
  std::cout << "  观测点数量: " << options.probe_points.size() << "\n";
  std::cout << "  邻域半径: " << options.radius << "\n";
  std::cout << "  压力字段: " << options.pressure_field << "\n";
  std::cout << "  mps_taylor最少邻域粒子数: " << options.mps_probe_min_neighbors << "\n";
  std::cout << "  测量方法: "
            << (options.method == Options::ProbeMethod::kAverageInRadius
                    ? "average"
                    : (options.method == Options::ProbeMethod::kNearestParticle
                           ? "nearest"
                           : "mps_taylor"))
            << "\n";

  return (processed > 0) ? 0 : 1;
}

