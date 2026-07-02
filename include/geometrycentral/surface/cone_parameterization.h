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

  // The interior cone vertices. From `flattenWithGivenConesAndCut`, these are
  // exactly the `cones` argument (vertices of that call's `mesh`). From
  // `parameterizeBFFwithCones`, these are reported on the *input* mesh (not the
  // internally-built intrinsic mesh the cut actually lives on).
  std::vector<Vertex> cones;

  // For each boundary edge of cutMesh, the index of the original (intrinsic)
  // edge it came from. The two sides of a cut slit share the same value; ordinary
  // boundary edges get distinct values. The seamless flatten guarantees both
  // sides of each shared value have equal length in the UV layout.
  EdgeData<int> boundarySeamId;

  // For each halfedge of cutMesh, the halfedge of the pre-cut mesh it was copied
  // from -- the map surgery::cutAlongEdges produces internally. Set for every
  // interior halfedge; the new exterior halfedges introduced along the cut (i.e.
  // cutMesh's boundary-loop halfedges) are Halfedge(). Lets a caller stitch
  // per-corner data between cutMesh and the pre-cut mesh per-HALFEDGE, which
  // cutAlongEdges' face-cycle-based copy does not otherwise preserve a
  // consistent starting point for.
  //
  // From `flattenWithGivenConesAndCut`, these halfedges belong to that call's
  // `mesh` argument -- directly usable by the caller, who owns it. From
  // `parameterizeBFFwithCones`, they belong to the internally-built intrinsic
  // mesh the cut actually lives on (NOT the function's `mesh` parameter, and
  // not otherwise exposed) -- so this field is only meaningful when calling
  // `flattenWithGivenConesAndCut` directly.
  HalfedgeData<Halfedge> parentHalfedge;
};

// Flatten a surface with a caller-chosen cone set and cut, using seamless
// Boundary First Flattening (Sawhney & Crane, ACM TOG 2017): steps 2-3 of
// `parameterizeBFFwithCones` below, factored out as a primitive for callers who
// place cones and/or route cuts themselves (e.g. on their own intrinsic
// triangulation) instead of using `computeConePlacement`/`computeConeCut`.
//
// `cones` are interior vertices of `mesh`; `coneAngles` prescribes their target
// curvature (only entries at `cones` are read -- pass the corresponding
// `ConePlacementResult::coneAngles` directly, or build your own). `cutEdges` is
// any edge set of `mesh` that, cut, yields a single topological disk with every
// vertex in `cones` on the resulting boundary loop -- it need not come from
// `computeConeCut`. That disk-and-cones-on-boundary precondition is asserted
// (throws std::runtime_error) rather than silently producing garbage.
//
// `nCutEdges == 0` with an empty `cones` is plain BFF (mesh is already a disk).
// Multiple disjoint cuts (several cones, several separate trees to the boundary)
// are fine. Boundary scale factors are held at zero. Deterministic: no RNG,
// index-order iteration only. Requires a surface WITH boundary: the uncut
// scale-factor solve is anchored at the original boundary. Depends only on
// intrinsic geometry.
ConeParameterizationResult flattenWithGivenConesAndCut(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geo,
                                                        const std::vector<Vertex>& cones,
                                                        const VertexData<double>& coneAngles,
                                                        const EdgeData<char>& cutEdges);

// Conformally flatten a surface with automatically-placed cone singularities:
//
//   1. place `nCones` interior cones minimizing area distortion (cone_placement),
//   2. cut geodesic slits connecting the cones to the boundary (cone_cut),
//   3. flatten with the cones' curvature prescribed (flattenWithGivenConesAndCut).
//
// Requires a surface WITH boundary (e.g. a disk patch). Depends only on
// intrinsic geometry.
ConeParameterizationResult parameterizeBFFwithCones(ManifoldSurfaceMesh& mesh,
                                                    IntrinsicGeometryInterface& geo, size_t nCones);

} // namespace surface
} // namespace geometrycentral
