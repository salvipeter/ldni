#include "ldni.hh"

#include <fstream>

using namespace Geometry;


// Main code

LDNI mesh2ldni(const TriMesh &mesh, const std::array<size_t, 3> resolution) {
  LDNI ldni;
  return ldni;
}

QuadMesh ldni2mesh(const LDNI &ldni) {
  QuadMesh mesh;
  return mesh;
}


// I/O

template<typename T>
static T readType(std::istream &is) {
  static char block[sizeof(T)];
  is.read(block, sizeof(T));
  return *reinterpret_cast<T *>(block);
}

Vector3D readVector(std::istream &is) {
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
  ldni.res[0] = readType<uint16_t>(f);
  ldni.res[1] = readType<uint16_t>(f);
  ldni.res[2] = readType<uint16_t>(f);
  for (size_t i = 0; i < 3; ++i) {
    size_t n = ldni.res[(i+1)%3] * ldni.res[(i+2)%3];
    ldni.cells[i].resize(n);
    for (size_t j = 0; j < n; ++j) {
      size_t m = readType<uint8_t>(f);
      ldni.cells[i][j].reserve(m);
      for (size_t k = 0; k < m; ++k) {
        double d = readType<float>(f);
        Vector3D n = readVector(f);
        ldni.cells[i][j].push_back({d, n});
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
  for (size_t i = 0; i < 3; ++i)
    for (const auto &cell : ldni.cells[i]) {
      writeType<uint8_t>(f, cell.size());
      for (const auto &dn : cell) {
        writeType<float>(f, dn.first);
        writeVector(f, dn.second);
      }
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
