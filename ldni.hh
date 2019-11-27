#pragma once

// Based on:
// C.L. Wang, Y-Sh. Leung, Y. Chen:
//   Solid modeling of polyhedral objects by Layered Depth-Normal Images on the GPU.
//   In: Computer-Aided Design, Vol. 42, pp. 535-544, 2010.

// Mesh building:
// T. Ju, F. Losasso, S. Schaefer, J. Warren:
//   Dual Contouring of Hermite Data.
//   In: SIGGRAPH'02, pp. 339-346, 2002.

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

struct DepthNormal {
  double d;
  Geometry::Vector3D n;
  DepthNormal(double d, const Geometry::Vector3D &n) : d(d), n(n) { }
  bool operator<(const DepthNormal &dn) const { return d < dn.d; }
};

using Ray = std::vector<DepthNormal>;

struct LDNI {
  Geometry::Point3D bbox[2];
  Geometry::Vector3D axis;     // = bbox[1] - bbox[0]
  Geometry::Vector3D dirs[3];  // = {{axis[0]/res[0],0,0},{0,axis[1]/res[1],0},{0,0,axis[2]/res[2]}}
  size_t res[3];               // # of cells in each direction
  std::vector<Ray> rays[3];
};

LDNI mesh2ldni(const Geometry::TriMesh &mesh, size_t size);

Geometry::QuadMesh ldni2mesh(const LDNI &ldni);

LDNI readLDNI(std::string filename);

void writeLDNI(const LDNI &ldni, std::string filename);
