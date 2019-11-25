#include <iostream>

#include "ldni.hh"

using namespace Geometry;

int main(int argc, char **argv) {
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] << " <model.obj> <resolution>" << std::endl;
    return 1;
  }
  size_t res = std::atoi(argv[2]);

  auto input_mesh = TriMesh::readOBJ(argv[1]);
  auto ldni = mesh2ldni(input_mesh, { res, res, res });
  writeLDNI(ldni, "/tmp/test.ldni");
  auto output_mesh = ldni2mesh(ldni);
  output_mesh.writeOBJ("/tmp/test.obj");

  // Test LDNI file
  ldni = readLDNI("/tmp/test.ldni");
  output_mesh = ldni2mesh(ldni);
  output_mesh.writeOBJ("/tmp/test2.obj");
}
