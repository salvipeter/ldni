#include <chrono>
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
  std::cout << "File loaded." << std::endl;

  std::chrono::steady_clock::time_point start, stop;
  start = std::chrono::steady_clock::now();
  auto ldni = mesh2ldni(input_mesh, res);
  stop = std::chrono::steady_clock::now();
  std::cout << "LDNI generation: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count()
            << "ms" << std::endl;
  std::cout << "Resolution: "
            << ldni.res[0] << "x" << ldni.res[1] << "x" << ldni.res[2] << std::endl;

  writeLDNI(ldni, "/tmp/test.ldni");
  std::cout << "LDNI file written (/tmp/test.ldni)." << std::endl;

  start = std::chrono::steady_clock::now();
  auto output_mesh = ldni2mesh(ldni);
  stop = std::chrono::steady_clock::now();
  std::cout << "Mesh generation: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count()
            << "ms" << std::endl;

  output_mesh.writeOBJ("/tmp/test.obj");
  std::cout << "Mesh file written (/tmp/test.obj)." << std::endl;

  // Test LDNI file
  ldni = readLDNI("/tmp/test.ldni");
  std::cout << "LDNI file read." << std::endl;
  output_mesh = ldni2mesh(ldni);
  std::cout << "Mesh generated." << std::endl;
  output_mesh.writeOBJ("/tmp/test2.obj");
  std::cout << "Mesh written (again)." << std::endl;
}
