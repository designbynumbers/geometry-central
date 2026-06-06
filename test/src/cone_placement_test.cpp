#include "geometrycentral/surface/cone_placement.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "load_test_meshes.h"

#include "gtest/gtest.h"

#include <cmath>
#include <set>

using namespace geometrycentral;
using namespace geometrycentral::surface;

class ConePlacementSuite : public MeshAssetSuite {};

// Helper: total prescribed curvature.
static double angleSum(const VertexData<double>& coneAngles, ManifoldSurfaceMesh& mesh) {
  double s = 0.0;
  for (Vertex v : mesh.vertices()) s += coneAngles[v];
  return s;
}

// On a closed genus-0 surface the requested number of interior cones are placed,
// they are distinct interior vertices, and the prescribed angles satisfy the
// discrete Gauss-Bonnet theorem (sum == 2*pi*chi == 4*pi).
TEST_F(ConePlacementSuite, ClosedGaussBonnet) {
  for (const MeshAsset& a : {getAsset("spot.ply", true), getAsset("sphere_small.ply", true)}) {
    a.printThyName();
    ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
    VertexPositionGeometry& geom = *a.geometry;
    ASSERT_EQ(mesh.nBoundaryLoops(), 0);

    size_t nCones = 5;
    ConePlacementResult res = computeConePlacement(mesh, geom, nCones);

    EXPECT_EQ(res.cones.size(), nCones);

    std::set<size_t> seen;
    for (Vertex v : res.cones) {
      EXPECT_FALSE(v.isBoundary());
      seen.insert(v.getIndex());
    }
    EXPECT_EQ(seen.size(), nCones); // distinct

    double target = 2.0 * M_PI * mesh.eulerCharacteristic();
    EXPECT_NEAR(angleSum(res.coneAngles, mesh), target, 1e-6);
  }
}

// On surfaces with boundary, cones are interior vertices, the requested count is
// placed, and Gauss-Bonnet holds (for non-zero Euler characteristic).
TEST_F(ConePlacementSuite, BoundaryInteriorCones) {
  for (const MeshAsset& a : boundaryMeshes()) {
    if (!a.isTriangular || a.manifoldMesh == nullptr || a.geometry == nullptr) continue;
    a.printThyName();
    ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
    VertexPositionGeometry& geom = *a.geometry;
    if (mesh.nInteriorVertices() < 3) continue;

    size_t nCones = 3;
    ConePlacementResult res = computeConePlacement(mesh, geom, nCones);

    EXPECT_EQ(res.cones.size(), nCones);
    for (Vertex v : res.cones) EXPECT_FALSE(v.isBoundary());

    if (mesh.eulerCharacteristic() != 0) {
      double target = 2.0 * M_PI * mesh.eulerCharacteristic();
      EXPECT_NEAR(angleSum(res.coneAngles, mesh), target, 1e-6);
    }
  }
}

// The greedy placement is deterministic: identical inputs give identical cones.
TEST_F(ConePlacementSuite, Deterministic) {
  MeshAsset a = getAsset("spot.ply", true);
  ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
  VertexPositionGeometry& geom = *a.geometry;

  ConePlacementResult r1 = computeConePlacement(mesh, geom, 6);
  ConePlacementResult r2 = computeConePlacement(mesh, geom, 6);

  ASSERT_EQ(r1.cones.size(), r2.cones.size());
  for (size_t i = 0; i < r1.cones.size(); i++) {
    EXPECT_EQ(r1.cones[i].getIndex(), r2.cones[i].getIndex());
  }
}
