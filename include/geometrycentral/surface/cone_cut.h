#pragma once

#include "geometrycentral/surface/flip_geodesics.h"
#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"

#include <memory>
#include <vector>

namespace geometrycentral {
namespace surface {

// Result of cone cutting: a geodesic cut, defined on an intrinsic triangulation.
struct ConeCutResult {
  // FlipEdgeNetwork owning the intrinsic triangulation the cut is defined on.
  // The cut edges below are keyed on `network->mesh` (the intrinsic mesh), and
  // edge lengths come from `network->tri`. Keeping this alive keeps the cut and
  // its geometry valid.
  std::unique_ptr<FlipEdgeNetwork> network;

  // Geodesic cut edges (keyed on network->mesh). Cutting the intrinsic mesh
  // along these (e.g. with surgery::cutAlongEdges) yields a topological disk
  // passing through the cones.
  EdgeData<bool> cutEdges;
};

// Compute geodesic cone-slit cuts which reduce the surface to a single
// topological disk passing through the given cone vertices.
//
// The cut TOPOLOGY is an approximate minimum-weight (by edge length) Steiner
// tree connecting all boundary loops and cone vertices, built with gc's
// MarkedDisjointSets (the cutter of Sawhney & Crane, "Boundary First
// Flattening", ACM TOG 2017). The cut GEOMETRY is then straightened to
// geodesics by a FlipEdgeNetwork on an internally-built intrinsic
// triangulation, so the slits are smooth rather than jagged mesh-edge polylines.
//
// Boundary edges are never cut. Cones are vertices of `mesh`. Cutting along the
// result (e.g. with surgery::cutAlongEdges) yields a single topological disk.
// For a closed surface at least two cones are required (the tree connecting them
// must have at least one edge).
ConeCutResult computeConeCut(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geo,
                             const std::vector<Vertex>& cones);

} // namespace surface
} // namespace geometrycentral
