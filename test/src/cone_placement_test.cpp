#include "geometrycentral/surface/cone_placement.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "load_test_meshes.h"

#include "gtest/gtest.h"

#include <cmath>
#include <memory>
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

// A triangulated, curved square patch: a topological disk with interior
// curvature (gc ships no clean disk asset).
static std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>>
makeBumpyDisk(int n) {
  std::vector<Vector3> pos;
  for (int i = 0; i <= n; i++)
    for (int j = 0; j <= n; j++) {
      double x = (double)i / n, y = (double)j / n;
      double z = 0.4 * std::exp(-8.0 * ((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5)));
      pos.push_back(Vector3{x, y, z});
    }
  auto idx = [&](int i, int j) { return (size_t)(i * (n + 1) + j); };
  std::vector<std::vector<size_t>> polys;
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
      polys.push_back({idx(i, j), idx(i + 1, j), idx(i + 1, j + 1)});
      polys.push_back({idx(i, j), idx(i + 1, j + 1), idx(i, j + 1)});
    }
  return makeManifoldSurfaceMeshAndGeometry(polys, pos);
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

// ConePlacementOptions::maxNewCones = 0 with a nonempty initialCones is an
// "angle refresh": it must return exactly the coneAngles (and uField) that the
// tail of an equivalent greedy run arrives at, since both solves anchor the
// same final cone set and the linear solve depends only on that set.
TEST(ConePlacementOptionsSuite, AngleRefreshMatchesGreedyTail) {
  std::unique_ptr<ManifoldSurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geom;
  std::tie(mesh, geom) = makeBumpyDisk(8);

  // Greedy run: place 2 cones one at a time, remembering the final set.
  ConePlacementOptions o1;
  o1.maxNewCones = 2;
  ConePlacementResult greedy = computeConePlacement(*mesh, *geom, o1);
  ASSERT_EQ(greedy.cones.size(), 2u);

  // Angle refresh: hand that exact final set back in as initialCones, ask for
  // no new cones.
  ConePlacementOptions o2;
  o2.initialCones = greedy.cones;
  o2.maxNewCones = 0;
  ConePlacementResult refresh = computeConePlacement(*mesh, *geom, o2);

  ASSERT_EQ(refresh.cones.size(), greedy.cones.size());
  for (size_t i = 0; i < refresh.cones.size(); i++)
    EXPECT_EQ(refresh.cones[i].getIndex(), greedy.cones[i].getIndex());

  for (Vertex v : mesh->vertices()) {
    EXPECT_NEAR(refresh.coneAngles[v], greedy.coneAngles[v], 1e-12);
    EXPECT_NEAR(refresh.uField[v], greedy.uField[v], 1e-12);
  }

  // uField is zero on the anchor set (boundary + cones).
  for (Vertex v : mesh->vertices())
    if (v.isBoundary()) EXPECT_EQ(refresh.uField[v], 0.0);
  for (Vertex c : refresh.cones) EXPECT_EQ(refresh.uField[c], 0.0);
}

// ConePlacementOptions::stopMaxU stops the greedy loop as soon as the gate is
// satisfied, even if maxNewCones budget remains -- including adding zero cones
// when the initial (pre-loop) field is already under threshold.
TEST(ConePlacementOptionsSuite, StopMaxUEarlyStop) {
  std::unique_ptr<ManifoldSurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geom;
  std::tie(mesh, geom) = makeBumpyDisk(8);

  // A huge threshold is satisfied immediately: no cones added even though
  // maxNewCones budget is nonzero.
  ConePlacementOptions oNone;
  oNone.maxNewCones = 5;
  oNone.stopMaxU = 1e6;
  ConePlacementResult none = computeConePlacement(*mesh, *geom, oNone);
  EXPECT_EQ(none.cones.size(), 0u);

  // An unreachable threshold (0, i.e. disabled by convention) exhausts the
  // full budget.
  ConePlacementOptions oFull;
  oFull.maxNewCones = 5;
  oFull.stopMaxU = 0.0;
  ConePlacementResult full = computeConePlacement(*mesh, *geom, oFull);
  EXPECT_EQ(full.cones.size(), 5u);

  // A threshold just under the pre-loop max|u| must stop after (at most) the
  // first addition -- adding the peak vertex as a new anchor can only reduce
  // the remaining max|u|, so the gate is satisfied well before the 5-cone
  // budget is exhausted.
  double maxU0 = 0.0;
  for (Vertex v : mesh->vertices())
    if (!v.isBoundary()) maxU0 = std::max(maxU0, std::abs(none.uField[v]));
  ConePlacementOptions oPartial;
  oPartial.maxNewCones = 5;
  oPartial.stopMaxU = 0.95 * maxU0;
  ConePlacementResult partial = computeConePlacement(*mesh, *geom, oPartial);
  EXPECT_LT(partial.cones.size(), 5u);
}
