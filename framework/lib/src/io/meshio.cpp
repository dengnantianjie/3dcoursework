#include "mh/io/meshio.h"

#include <iostream>
#include <fstream>

namespace mh
{

  std::shared_ptr<Mesh> parsePly(std::string const path)
  {
      std::ifstream f(path);
      if ( !f.is_open() )
      {
          std::cerr << "[meshio.cpp::parsePly()] could not open " << path << std::endl;
          return nullptr;
      }

      std::vector<Eigen::Vector3f> vertData; // filled with read data
      std::vector<Eigen::Vector3f> normalData; // filled empty
      std::vector<Eigen::Vector3i> faceData; // unused

      long int nVertices(-1);
      long int vertexId(0);
      bool endHeaderFound(false);
      std::string line;
      while ( std::getline(f,line) && (vertexId != nVertices) )
      {
          if ( endHeaderFound )
          {
              std::istringstream iss(line);
              Eigen::Vector3f point;
              iss >> point(0) >> point(1) >> point(2);
              vertData.push_back(point);
              ++vertexId;
              continue;
          }

          // find vertex count, then end of header
          static const std::string tag("element vertex");
          size_t index = line.find(tag);
          if ( index != std::string::npos )
          {
              std::cout << "[meshio.cpp::parsePly()] parsing vertex count from line: " << line << std::endl;
              nVertices = std::atol( line.substr(index+tag.size(),line.size()-index+tag.size()).c_str() );
              std::cout << "[meshio.cpp::parsePly()] parsed vertex count: " << nVertices << std::endl;
              vertData.reserve(nVertices);
              normalData.resize( nVertices, Eigen::Vector3f::Zero() );
          }
          else if ( line.find("end_header") != std::string::npos )
          {
              if ( nVertices == (long int)(-1) )
              {
                  std::cerr << "[meshio.cpp::parsePly()] reached end of ply header, but did not find vertex count...stopping" << std::endl;
                  return nullptr;
              }
              endHeaderFound = true;
          }
          else
              std::cout << "[meshio.cpp::parsePly()] skipping line " << line << std::endl;
      }

      f.close();

      return std::make_shared<Mesh>( vertData, normalData, faceData );
  } //...parsePly()

#   if USE_ASSIMP
    Mesh convertFromAssimp(const aiMesh * inputMesh)
    {
        std::vector<Eigen::Vector3f> vertData;
        std::vector<Eigen::Vector3f> normalData;
        std::vector<Eigen::Vector3i> faceData;

        for (size_t i = 0; i < inputMesh->mNumVertices; ++i)
        {
            vertData.push_back  (Eigen::Vector3f(inputMesh->mVertices[i].x,
                inputMesh->mVertices[i].y,
                inputMesh->mVertices[i].z));
            if (inputMesh->HasNormals())
            {
                normalData.push_back(Eigen::Vector3f(inputMesh->mNormals[i].x,
                    inputMesh->mNormals[i].y,
                    inputMesh->mNormals[i].z));
            } else {
                normalData.push_back(Eigen::Vector3f(0.0f, 0.0f, 0.0f));
            }
        }

        if (inputMesh->HasFaces())
        {
            for (size_t i = 0; i < inputMesh->mNumFaces; ++i)
            {
                const aiFace * face = &(inputMesh->mFaces[i]);
                faceData.push_back(Eigen::Vector3i(face->mIndices[0],
                    face->mIndices[1],
                    face->mIndices[2]));
            }
        }

        Mesh mesh(vertData, normalData, faceData);

        return mesh;
    } //...convertFromAssimp()
#   endif // USE_ASSIMP

    std::vector<std::shared_ptr<Mesh> > loadMeshesFromFile(std::string path)
    {
        std::vector<std::shared_ptr<Mesh> > meshes;

        // Temporary fix: assimp cannot take ply point clouds without faces
        if ( path.find("ply") == path.size() - 3 )
        {
            std::cout << "[meshio.cpp::loadMeshesFromFile()] loading " << path << std::endl;
            std::shared_ptr<Mesh> mesh = parsePly(path);
            if ( mesh )
                meshes.emplace_back( mesh );
            else
                std::cerr << "Could not parse ply file " << path << std::endl;
        }
        else
        {
#       if USE_ASSIMP
            Assimp::Importer importer;

            const aiScene *scene = importer.ReadFile(path.c_str(), aiProcess_JoinIdenticalVertices
            );

            MH_ASSERT(scene);

            meshes.reserve(scene->mNumMeshes);
            for (size_t i = 0; i < scene->mNumMeshes; ++i)
            {
                auto mesh = std::make_shared<Mesh>(convertFromAssimp(scene->mMeshes[i]));
                meshes.push_back(mesh);
            }
#       else
            std::cerr << "Enable Assimp library to load other types of meshes" << std::endl;
#       endif // USE_ASSIMP
        }

        std::cout << "[meshio.cpp::loadMeshesFromFile()] \"meshes\" has " << meshes.size() << " mesh(es)" << std::endl;
        for ( size_t meshId = 0; meshId != meshes.size(); ++meshId )
        {
            Mesh const& mesh = *meshes.at(meshId);
            std::cout << "[meshio.cpp::loadMeshesFromFile()]\t Mesh" << meshId << " has " << mesh.getVertices().size() << " vertices" << std::endl;
        }
        return meshes;
    } //...loadMeshesFromFile()
} // namespace mh