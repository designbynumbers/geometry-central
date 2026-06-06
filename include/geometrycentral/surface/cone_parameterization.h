#pragma once

#include "geometrycentral/surface/edge_length_geometry.h"
#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/utilities/vector2.h"

#include <memory>
#include <vector>

namespace geometrycentral {
namespace surface {

// Result of cone parameterization: a flattened cut disk.
struct ConeParameterizationResult {
  // The cut mesh (a topological disk) and its intrinsic geometry. The cones and
  // the original boundary together form `cutMesh`'s single boundary loop.
  std::unique_ptr<ManifoldSurfaceMesh> cutMesh;
  std::unique_ptr<EdgeLengthGeometry> cutGeometry;

  // UV coordinates of the flattening, per vertex of cutMesh.
  VertexData<Vector2> uvs;

  // The interior cone vertices that were placed (vertices of the *input* mesh).
  std::vector<Vertex> cones;

  // For each boundary edge of cutMesh, the index of the original (intrinsic)
  // edge it came from. The two sides of a cut slit share the same value; ordinary
  // boundary edges get distinct values. The seamless flatten guarantees both
  // sides of each shared value have equal length in the UV layout.
  EdgeData<int> boundarySeamId;
};

// Conformally flatten a surface with automatically-placed cone singularities,
// using a seamless Boundary First Flattening (Sawhney & Crane, ACM TOG 2017):
//
//   1. place `nCones` interior cones minimizing area distortion (cone_placement),
//   2. cut geodesic slits connecting the cones to the boundary (cone_cut),
//   3. flatten the cut disk with the cones' curvature prescribed, laying out the
//      boundary so that the two sides of every slit have EQUAL length (the
//      seamless "closeLengths" step) -- so the cut can be re-joined without
//      stretch mismatch.
//
// Requires a surface WITH boundary (e.g. a disk patch): the uncut scale-factor
// solve is anchored at the original boundary. Depends only on intrinsic geometry.
ConeParameterizationResult parameterizeBFFwithCones(ManifoldSurfaceMesh& mesh,
                                                    IntrinsicGeometryInterface& geo, size_t nCones);

} // namespace surface
} // namespace geometrycentral
