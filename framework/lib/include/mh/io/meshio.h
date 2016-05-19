#ifndef MESHIO_H
#define MESHIO_H 

#include <string>
#include <map>

#if USE_ASSIMP
#   include <assimp/mesh.h>
#   include <assimp/postprocess.h>
#   include <assimp/scene.h>
#   include <assimp/Importer.hpp>
#endif

#include "mh/3d/mesh.h"

namespace mh
{

#if USE_ASSIMP
  //  Mesh                                convertFromAssimp (const aiMesh * inputMesh, float const eps = 1e-5f);
	Mesh convertFromAssimp(const aiMesh * inputMesh);
#endif
std::vector<std::shared_ptr<Mesh> > loadMeshesFromFile (std::string path);

} // namespace mh
#endif /* MESHIO_H */
