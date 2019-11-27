#include "cell-point.hh"

#include <Eigen/Dense>

using namespace Geometry;

Point3D cellPoint(const std::vector<Point3D> &points, const std::vector<Vector3D> &normals,
                  const Point3D &cell_min, const Point3D &cell_max) {
  size_t m = points.size();
  Point3D centroid(0, 0, 0);
  for (const auto &p : points)
    centroid += p;
  centroid /= m;

  // Setup QEF
  Eigen::MatrixXd A(m + 3, 3);
  Eigen::VectorXd b(m + 3);
  for (size_t i = 0; i < m; ++i) {
    const auto &n = normals[i];
    for (size_t j = 0; j < 3; ++j)
      A(i, j) = n[j];
    b(i) = n * points[i];
  }
  A(  m  , 0) = 1; A(  m  , 1) = 0; A(  m  , 2) = 0; b(  m  ) = centroid[0];
  A(m + 1, 0) = 0; A(m + 1, 1) = 1; A(m + 1, 2) = 0; b(m + 1) = centroid[1];
  A(m + 2, 0) = 0; A(m + 2, 1) = 0; A(m + 2, 2) = 1; b(m + 2) = centroid[2];

  // Solve
  Eigen::Vector3d x = A.colPivHouseholderQr().solve(b);
  Point3D result(x(0), x(1), x(2));

  // Check if it is within bounds
  if (result[0] >= cell_min[0] && result[0] <= cell_max[0] &&
      result[1] >= cell_min[1] && result[1] <= cell_max[1] &&
      result[2] >= cell_min[2] && result[2] <= cell_max[2])
    return result;

  // Fall back to the mass center of the crossing points
  return centroid;
}
