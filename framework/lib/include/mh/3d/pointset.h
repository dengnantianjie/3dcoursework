#ifndef POINTSET_H
#define POINTSET_H 

#include "mh/base/defs.h"
#include "mh/base/imports.h"

#include <memory>
#include <vector>

#include <GL/glew.h>

#include "eigen3/Eigen/Geometry"

#include "mh/3d/vertex.h"

namespace mh
{

class Pointset
{
public:
    static const int POSITION_LOCATION = 0;
    static const int NORMAL_LOCATION   = 1;
    static const int COLOR_LOCATION    = 2;

public:
    Pointset (const std::vector<Eigen::Vector3f> & vertData,
              const std::vector<Eigen::Vector3f> & normalData);

    const std::vector<std::shared_ptr<Vertex> > & getVertices() const { return m_verts; }

protected:

private:
    // index
    std::size_t m_idx;

    // raw data - not accessible directly through public interface
    // std::vecto
    std::vector<std::shared_ptr<Vertex> > m_verts;

}; // class Pointset

} // namespace mh
#endif /* POINTSET_H */
