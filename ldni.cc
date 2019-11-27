#include "ldni.hh"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <optional>

#include "cell-point.hh"

using namespace Geometry;


// LDNI data generation

inline double triangleArea(double ax, double ay, double bx, double by, double cx, double cy) {
  return 0.5 * ((ax - cx) * (by - ay) - (ax - bx) * (cy - ay)); // signed area
}

// Check if the (a,b,c) triangle is hit by any rays parallel to the k0 axis,
// and update the LDNI data ray structure accordingly.
static
void checkTriangle(LDNI &ldni, size_t k0, const Point3D &a, const Point3D &b, const Point3D &c) {
  Vector3D n = ((b - a) ^ (c - a)).normalize();
  size_t k1 = (k0 + 1) % 3, k2 = (k0 + 2) % 3;
  double area = triangleArea(a[k1], a[k2], b[k1], b[k2], c[k1], c[k2]);
  if (std::abs(area) < epsilon)
    return;

  // Select the possible rays
  double minx = std::min(a[k1], std::min(b[k1], c[k1]));
  double maxx = std::max(a[k1], std::max(b[k1], c[k1]));
  double miny = std::min(a[k2], std::min(b[k2], c[k2]));
  double maxy = std::max(a[k2], std::max(b[k2], c[k2]));
  minx = std::floor((minx - ldni.bbox[0][k1]) * (ldni.res[k1] - 1) / ldni.axis[k1]);
  maxx =  std::ceil((maxx - ldni.bbox[0][k1]) * (ldni.res[k1] - 1) / ldni.axis[k1]);
  miny = std::floor((miny - ldni.bbox[0][k2]) * (ldni.res[k2] - 1) / ldni.axis[k2]);
  maxy =  std::ceil((maxy - ldni.bbox[0][k2]) * (ldni.res[k2] - 1) / ldni.axis[k2]);
  size_t mini = std::min((size_t)std::max(minx, 0.0), ldni.res[k1] - 1);
  size_t maxi = std::min((size_t)std::max(maxx, 0.0), ldni.res[k1] - 1);
  size_t minj = std::min((size_t)std::max(miny, 0.0), ldni.res[k2] - 1);
  size_t maxj = std::min((size_t)std::max(maxy, 0.0), ldni.res[k2] - 1);

  // Check the intersections
  for (size_t i = mini; i <= maxi; ++i) {
    double u = ldni.bbox[0][k1] + i * ldni.dirs[k1][k1];
    for (size_t j = minj; j <= maxj; ++j) {
      double v = ldni.bbox[0][k2] + j * ldni.dirs[k2][k2];
      // Compute barycentric coordinates
      double la = triangleArea(u, v, b[k1], b[k2], c[k1], c[k2]) / area;
      double lb = triangleArea(u, v, c[k1], c[k2], a[k1], a[k2]) / area;
      double lc = triangleArea(u, v, a[k1], a[k2], b[k1], b[k2]) / area;
      if (0 <= la && la <= 1 && 0 <= lb && lb <= 1 && 0 <= lc && lc <= 1) {
        double d = a[k0] * la + b[k0] * lb + c[k0] * lc - ldni.bbox[0][k0];
        ldni.rays[k0][i*ldni.res[k2]+j].emplace_back(d, n);
      }
    }
  }
}

LDNI mesh2ldni(const TriMesh &mesh, size_t size) {
  // Compute the bounding box
  Point3D boxmin, boxmax;
  const auto &points = mesh.points();
  boxmin = boxmax = points[0];
  for (const auto &p : points)
    for (int i = 0; i < 3; ++i) {
      boxmin[i] = std::min(boxmin[i], p[i]);
      boxmax[i] = std::max(boxmax[i], p[i]);
    }
  // Add 5%
  auto mean = (boxmin + boxmax) / 2;
  boxmin = mean + (boxmin - mean) * 1.05;
  boxmax = mean + (boxmax - mean) * 1.05;

  LDNI ldni;
  ldni.bbox[0] = boxmin; ldni.bbox[1] = boxmax;
  ldni.axis = boxmax - boxmin;

  // Compute the resolution
  double axis_delta = ldni.axis.norm() / size / std::sqrt(3);
  ldni.res[0] = std::max<size_t>((size_t)std::ceil(ldni.axis[0] / axis_delta) + 1, 2);
  ldni.res[1] = std::max<size_t>((size_t)std::ceil(ldni.axis[1] / axis_delta) + 1, 2);
  ldni.res[2] = std::max<size_t>((size_t)std::ceil(ldni.axis[2] / axis_delta) + 1, 2);

  // Ray edge vectors
  ldni.dirs[0] = Vector3D(ldni.axis[0] / (ldni.res[0] - 1), 0, 0);
  ldni.dirs[1] = Vector3D(0, ldni.axis[1] / (ldni.res[1] - 1), 0);
  ldni.dirs[2] = Vector3D(0, 0, ldni.axis[2] / (ldni.res[2] - 1));

  // Allocate memory
  ldni.rays[0].resize(ldni.res[1] * ldni.res[2]);
  ldni.rays[1].resize(ldni.res[2] * ldni.res[0]);
  ldni.rays[2].resize(ldni.res[0] * ldni.res[1]);

  // Find the ray intersections
  for (const auto &tri : mesh.triangles()) {
    for (int c = 0; c < 3; ++c)
      checkTriangle(ldni, c, points[tri[0]], points[tri[1]], points[tri[2]]);
  }
  // ... and sort them
  for (int c = 0; c < 3; ++c)
    for (auto &ray : ldni.rays[c])
      std::sort(ray.begin(), ray.end());

  return ldni;
}


// Contouring

static bool insidep(const LDNI &ldni, const std::array<size_t, 3> &index) {
  int votes = 0;
  for (int c0 = 0; c0 < 3; ++c0) {
    int c1 = (c0 + 1) % 3, c2 = (c0 + 2) % 3;
    const auto &ray = ldni.rays[c0][index[c1]*ldni.res[c2]+index[c2]];
    double distance = index[c0] * ldni.dirs[c0][c0];
    bool inside = false;
    for (const auto &dn : ray) {
      if (dn.d > distance)
        break;
      inside = !inside;
    }
    if (inside)
      ++votes;
  }
  return votes >= 2;
}

// Finds the first crossing on the edge parallel to the c0 axis,
// at the given index + 0 or 1 in the two non-axis directions (given by d1 and d2).
// Returns std::nullopt when no crossing is found.
static std::optional<std::pair<Point3D, Vector3D>>
findCrossing(const LDNI &ldni, const std::array<size_t, 3> &index, int c0, size_t d1, size_t d2) {
  static std::array<Vector3D, 3> units =
    { { { 1, 0, 0 },
        { 0, 1, 0 },
        { 0, 0, 1 } } };
  int c1 = (c0 + 1) % 3, c2 = (c0 + 2) % 3;
  const auto &ray = ldni.rays[c0][(index[c1]+d1)*ldni.res[c2]+index[c2]+d2];
  double dist_min = index[c0] * ldni.dirs[c0][c0];
  double dist_max = dist_min + ldni.dirs[c0][c0];
  for (const auto &dn : ray) {
    if (dn.d > dist_max)
      break;
    if (dn.d >= dist_min) {
      auto p = ldni.bbox[0] + units[c0] * dn.d +
        ldni.dirs[c1] * (index[c1] + d1) + ldni.dirs[c2] * (index[c2] + d2);
      return { { p, dn.n } };
    }
  }
  return {};
}

std::vector<size_t> addPoints(QuadMesh &mesh, const LDNI &ldni) {
  std::vector<size_t> cells;
  cells.reserve((ldni.res[0] - 1) * (ldni.res[1] - 1) * (ldni.res[2] - 1));

  Vector3D delta(ldni.dirs[0][0], ldni.dirs[1][1], ldni.dirs[2][2]);

  size_t point_index = 1;
  for (size_t i = 0; i < ldni.res[0] - 1; ++i) {
    for (size_t j = 0; j < ldni.res[1] - 1; ++j) {
      for (size_t k = 0; k < ldni.res[2] - 1; ++k) {

        // Check if it is an interesting cell
        bool found_inside = false, found_outside = false;
        for (size_t di = 0, vi = 0; di <= 1; ++di)
          for (size_t dj = 0; dj <= 1; ++dj)
            for (size_t dk = 0; dk <= 1; ++dk, ++vi)
              if (insidep(ldni, { i + di, j + dj, k + dk }))
                found_inside = true;
              else
                found_outside = true;
        if (!found_inside || !found_outside) {
          cells.push_back(0);
          continue;
        }

        // Compute crossing data
        std::vector<Point3D> points;
        std::vector<Vector3D> normals;
        for (int c = 0; c < 3; ++c)
          for (size_t di = 0; di <= 1; ++di)
            for (size_t dj = 0; dj <= 1; ++dj) {
              auto dn = findCrossing(ldni, { i, j, k }, c, di, dj);
              if (dn) {
                points.push_back(dn->first);
                normals.push_back(dn->second);
              }
            }
        if (points.empty()) {
          cells.push_back(0);
          continue;
        }

        Point3D origin = ldni.bbox[0] + Vector3D(delta[0] * i, delta[1] * j, delta[2] * k);
        auto surface_point = cellPoint(points, normals, origin, origin + delta);
        // auto surface_point = origin + delta / 2; // Cell centers - "blocky" mesh
        mesh.addPoint(surface_point);
        cells.push_back(point_index++);
      }
    }
  }

  return cells;
}

void addQuads(QuadMesh &mesh, const LDNI &ldni, const std::vector<size_t> &cells) {
  std::array<size_t, 3> ns = { (ldni.res[1] - 1) * (ldni.res[2] - 1), ldni.res[2] - 1, 1 };
  for (size_t c0 = 0; c0 < 3; ++c0) {
    int c1 = (c0 + 1) % 3, c2 = (c0 + 2) % 3;
    size_t ni = ns[c0], nj = ns[c1], nk = ns[c2];
    for (size_t i = 0; i < ldni.res[c0] - 1; ++i) {
      for (size_t j = 1; j < ldni.res[c1] - 1; ++j) {
        for (size_t k = 1; k < ldni.res[c2] - 1; ++k) {
          size_t index = i * ni + j * nj + k * nk;
          size_t a = cells[index], b = cells[index-nj], c = cells[index-nj-nk], d = cells[index-nk];
          if (a * b * c * d == 0)
            continue;
          static std::array<size_t, 3> index1, index2;
          index1[c0] = i;     index1[c1] = j; index1[c2] = k;
          index2[c0] = i + 1; index2[c1] = j; index2[c2] = k;
          bool inside = insidep(ldni, index1);
          if (inside == insidep(ldni, index2))
            continue;
          if (inside)
            mesh.addQuad(a, b, c, d);
          else
            mesh.addQuad(d, c, b, a);
        }
      }
    }
  }
}

QuadMesh ldni2mesh(const LDNI &ldni) {
  QuadMesh mesh;
  auto cells = addPoints(mesh, ldni);
  addQuads(mesh, ldni, cells);
  return mesh;
}


// I/O

template<typename T>
static T readType(std::istream &is) {
  static char block[sizeof(T)];
  is.read(block, sizeof(T));
  return *reinterpret_cast<T *>(block);
}

static Vector3D readVector(std::istream &is) {
  float x = readType<float>(is);
  float y = readType<float>(is);
  float z = readType<float>(is);
  return { x, y, z };
}

LDNI readLDNI(std::string filename) {
  LDNI ldni;
  std::ifstream f(filename, std::ios::binary);
  ldni.bbox[0] = readVector(f);
  ldni.bbox[1] = readVector(f);
  ldni.axis = ldni.bbox[1] - ldni.bbox[0];
  ldni.res[0] = readType<uint16_t>(f);
  ldni.res[1] = readType<uint16_t>(f);
  ldni.res[2] = readType<uint16_t>(f);
  ldni.dirs[0] = Vector3D(ldni.axis[0] / (ldni.res[0] - 1), 0, 0);
  ldni.dirs[1] = Vector3D(0, ldni.axis[1] / (ldni.res[1] - 1), 0);
  ldni.dirs[2] = Vector3D(0, 0, ldni.axis[2] / (ldni.res[2] - 1));
  for (int i = 0; i < 3; ++i) {
    size_t n = ldni.res[(i+1)%3] * ldni.res[(i+2)%3];
    ldni.rays[i].resize(n);
    for (size_t j = 0; j < n; ++j) {
      size_t m = readType<uint8_t>(f);
      ldni.rays[i][j].reserve(m);
      for (size_t k = 0; k < m; ++k) {
        Vector3D n = readVector(f);
        double d = n.norm();
        ldni.rays[i][j].emplace_back(d, n.normalize());
      }
    }
  }
  return ldni;
}

template<typename T>
static void writeType(std::ostream &os, T x) {
  os.write(reinterpret_cast<const char *>(&x), sizeof(T));
}

static void writeVector(std::ostream &os, const Vector3D &v) {
  writeType<float>(os, v[0]);
  writeType<float>(os, v[1]);
  writeType<float>(os, v[2]);
}

void writeLDNI(const LDNI &ldni, std::string filename) {
  // Limits:
  // - maximum 65535x65535x65535 resolution
  // - maximum 255 intersections on one ray
  std::ofstream f(filename, std::ios::binary);
  f.exceptions(std::ios::failbit | std::ios::badbit);
  writeVector(f, ldni.bbox[0]);
  writeVector(f, ldni.bbox[1]);
  writeType<uint16_t>(f, ldni.res[0]);
  writeType<uint16_t>(f, ldni.res[1]);
  writeType<uint16_t>(f, ldni.res[2]);
  for (int i = 0; i < 3; ++i)
    for (const auto &ray : ldni.rays[i]) {
      writeType<uint8_t>(f, ray.size());
      for (const auto &dn : ray)
        writeVector(f, dn.n * dn.d);
    }
}


// Basic quad mesh implementation (to be removed later)

void
QuadMesh::addPoint(const Point3D &p) {
  points.push_back(p);
}

void
QuadMesh::addQuad(size_t a, size_t b, size_t c, size_t d) {
  quads.push_back({a, b, c, d});
}

void
QuadMesh::writeOBJ(std::string filename) const {
  std::ofstream f(filename);
  f.exceptions(std::ios::failbit | std::ios::badbit);
  for (const auto &p : points)
    f << "v " << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
  for (const auto &t : quads)
    f << "f " << t[0] << ' ' << t[1] << ' ' << t[2] << ' ' << t[3] << std::endl;
}
