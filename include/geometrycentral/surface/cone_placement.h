#pragma once

#include "geometrycentral/numerical/linear_algebra_utilities.h"
#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"

#include <vector>

namespace geometrycentral {
namespace surface {

// Result of automatic cone-singularity placement.
struct ConePlacementResult {
  // The chosen interior cone vertices (any caller-supplied `initialCones` plus
  // any added by the greedy loop), in vertex iteration order.
  std::vector<Vertex> cones;

  // Target angle (prescribed curvature) at each vertex, normalized so the total
  // prescribed curvature satisfies the discrete Gauss-Bonnet theorem
  // (sum == 2*pi*eulerCharacteristic). Entries are nonzero only at interior cone
  // vertices and at boundary vertices (which carry their exterior angle); all
  // other (ordinary interior) vertices are 0.
  VertexData<double> coneAngles;

  // The conformal scale-factor field from the final solve (the same field the
  // greedy loop uses internally to pick the next cone). Zero on the anchor set
  // (boundary vertices, any closed-surface auto-seed, and all cones). Useful for
  // gate metrics computed on the field itself (e.g. area-weighted percentiles)
  // rather than just the vertex the greedy step happened to pick.
  VertexData<double> uField;
};

struct ConePlacementOptions {
  // Cones already chosen (kept fixed, added to the anchor set before the greedy
  // loop runs). Empty = current behavior. Must be interior vertices.
  std::vector<Vertex> initialCones;

  // Number of NEW cones the greedy loop may add, on top of `initialCones` (and,
  // for a closed surface, the automatic extreme-curvature seed). 0 = no greedy
  // step: just solve the scale-factor field and prescribe/normalize angles for
  // `initialCones` as they stand -- the "angle refresh" needed when re-flattening
  // a patch whose cone set the caller manages itself.
  size_t maxNewCones = 1;

  // Optional early stop: if > 0, stop the greedy loop as soon as the largest
  // |u| among free (non-anchor) vertices drops below this threshold. Checked
  // after every solve, including the initial one -- so a patch already under
  // the threshold adds no cones even if `maxNewCones > 0`.
  double stopMaxU = 0.0;
};

// Automatic cone placement via the greedy CETM strategy of
//   Sawhney & Crane, "Boundary First Flattening", ACM TOG 2017.
//
// Starting from an initial anchor set (all boundary vertices for a surface with
// boundary, or the single extreme-curvature vertex for a closed surface, plus
// `opts.initialCones`), repeatedly solves for the conformal scale-factor field
// and adds the vertex of largest |scale factor| as the next cone, up to
// `opts.maxNewCones` additions (or until `opts.stopMaxU` is satisfied). The
// resulting cone angles are normalized to satisfy discrete Gauss-Bonnet.
//
// The mesh may be closed or have boundary; it need not be a disk. Depends only
// on the intrinsic geometry.
ConePlacementResult computeConePlacement(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geo,
                                         const ConePlacementOptions& opts);

// Convenience overload: place `nCones` total interior cones (the closed-surface
// auto-seed, if any, counts toward it) with no caller-supplied initial set.
// Equivalent to `computeConePlacement(mesh, geo, {{}, nCones - autoSeedCount, 0})`.
ConePlacementResult computeConePlacement(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geo,
                                         size_t nCones);

} // namespace surface
} // namespace geometrycentral
