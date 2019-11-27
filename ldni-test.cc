#include <chrono>
#include <iostream>

#include "ldni.hh"

using namespace Geometry;

int main(int argc, char **argv) {
  if (argc != 2 && argc != 3) {
    std::cerr << "Usage:" << std::endl;
    std::cerr << "a) " << argv[0] << " <model.obj> <resolution>" << std::endl;
    std::cerr << "b) " << argv[0] << " <model.ldni>" << std::endl;
    return 1;
  }

  LDNI ldni;
  std::chrono::steady_clock::time_point start, stop;

  if (argc == 2) {
    ldni = readLDNI(argv[1]);
    std::cout << "LDNI file read." << std::endl;
  }

  if (argc == 3) {
    size_t res = std::atoi(argv[2]);

    auto input_mesh = TriMesh::readOBJ(argv[1]);
    std::cout << "File loaded." << std::endl;

    start = std::chrono::steady_clock::now();
    ldni = mesh2ldni(input_mesh, res);
    stop = std::chrono::steady_clock::now();
    std::cout << "LDNI generation: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count()
              << "ms" << std::endl;
    std::cout << "Resolution: "
              << ldni.res[0] - 1 << "x" << ldni.res[1] - 1 << "x" << ldni.res[2] - 1 << std::endl;

    writeLDNI(ldni, "/tmp/test.ldni");
    std::cout << "LDNI file written (/tmp/test.ldni)." << std::endl;
  }

  start = std::chrono::steady_clock::now();
  auto output_mesh = ldni2mesh(ldni);
  stop = std::chrono::steady_clock::now();
  std::cout << "Mesh generation: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count()
            << "ms" << std::endl;

  output_mesh.writeOBJ("/tmp/test.obj");
  std::cout << "Mesh file written (/tmp/test.obj)." << std::endl;
}
