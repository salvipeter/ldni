#include "cell-point.hh"

#include <Eigen/Dense>

using namespace Geometry;

Point3D cellPoint(const std::vector<Point3D> &points, const std::vector<Vector3D> &normals,
                  const Point3D &cell_min, const Point3D &cell_max) {
  size_t m = points.size();
  if (m >= 3) {
    // Setup QEF
    Eigen::MatrixXd A(m, 3);
    Eigen::VectorXd b(m);
    for (size_t i = 0; i < m; ++i) {
      const auto &n = normals[i];
      for (size_t j = 0; j < 3; ++j)
        A(i, j) = n[j];
      b(i) = n * points[i];
    }

    // Solve
    Eigen::Vector3d x = A.colPivHouseholderQr().solve(b);
    Point3D result(x(0), x(1), x(2));

    // Check if it is within bounds
    if (result[0] >= cell_min[0] && result[0] <= cell_max[0] &&
        result[1] >= cell_min[1] && result[1] <= cell_max[1] &&
        result[2] >= cell_min[2] && result[2] <= cell_max[2])
      return result;
  }

  // Fall back to the mass center of the crossing points
  Point3D result(0, 0, 0);
  for (const auto &p : points)
    result += p;
  return result / m;
}
