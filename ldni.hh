#pragma once

#include <geometry.hh>
namespace Geometry {
  struct QuadMesh {             // placeholder implementation
    using Quad = std::array<size_t, 4>;
    std::vector<Point3D> points;
    std::vector<Quad> quads;

    void addPoint(const Point3D &p);
    void addQuad(size_t a, size_t b, size_t c, size_t d);
    void writeOBJ(std::string filename) const;
  };
}

using Cell = std::vector<std::pair<double, Geometry::Vector3D>>;

struct LDNI {
  Geometry::Point3D bbox[2];
  size_t res[3];
  std::vector<Cell> cells[3];
};

LDNI mesh2ldni(const Geometry::TriMesh &mesh, const std::array<size_t, 3> resolution);

Geometry::QuadMesh ldni2mesh(const LDNI &ldni);

LDNI readLDNI(std::string filename);

void writeLDNI(const LDNI &ldni, std::string filename);
