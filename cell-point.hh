#pragma once

#include <geometry.hh>

Geometry::Point3D
cellPoint(const std::vector<Geometry::Point3D> &points,
          const std::vector<Geometry::Vector3D> &normals,
          const Geometry::Point3D &cell_min,
          const Geometry::Point3D &cell_max);
