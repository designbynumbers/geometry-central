#pragma once

#include "geometrycentral/numerical/linear_algebra_utilities.h"
#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"

#include <vector>

namespace geometrycentral {
namespace surface {

// Result of automatic cone-singularity placement.
struct ConePlacementResult {
  // The chosen interior cone vertices, in vertex iteration order.
  std::vector<Vertex> cones;

  // Target angle (prescribed curvature) at each vertex, normalized so the total
  // prescribed curvature satisfies the discrete Gauss-Bonnet theorem
  // (sum == 2*pi*eulerCharacteristic). Entries are nonzero only at interior cone
  // vertices and at boundary vertices (which carry their exterior angle); all
  // other (ordinary interior) vertices are 0.
  VertexData<double> coneAngles;
};

// Automatic cone placement via the greedy CETM strategy of
//   Sawhney & Crane, "Boundary First Flattening", ACM TOG 2017.
//
// Selects up to `nCones` interior cone singularities that minimize area
// distortion: starting from an initial anchor set (all boundary vertices for a
// surface with boundary, or the single extreme-curvature vertex for a closed
// surface), it repeatedly solves for the conformal scale-factor field and adds
// the vertex of largest |scale factor| as the next cone. The resulting cone
// angles are normalized to satisfy discrete Gauss-Bonnet.
//
// `nCones` is the total number of interior cones requested (the closed-surface
// anchor counts toward it). The mesh may be closed or have boundary; it need not
// be a disk. Depends only on the intrinsic geometry.
ConePlacementResult computeConePlacement(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geo,
                                         size_t nCones);

} // namespace surface
} // namespace geometrycentral
