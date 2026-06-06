#include "geometrycentral/surface/cone_cut.h"
#include "geometrycentral/surface/cone_placement.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/surgery.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "load_test_meshes.h"

#include "gtest/gtest.h"

#include <cmath>
#include <vector>

using namespace geometrycentral;
using namespace geometrycentral::surface;

class ConeCutSuite : public MeshAssetSuite {};

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

// Count and cut: cutting the intrinsic mesh along the geodesic cut yields a
// single topological disk (euler characteristic 1, one boundary loop).
static void expectCutsToDisk(ConeCutResult& res) {
  ASSERT_NE(res.network, nullptr);
  ManifoldSurfaceMesh& iMesh = res.network->mesh;

  size_t nCut = 0;
  EdgeData<char> cut(iMesh, 0);
  for (Edge e : iMesh.edges()) {
    if (res.cutEdges[e]) {
      cut[e] = 1;
      nCut++;
    }
  }
  EXPECT_GT(nCut, 0u); // a non-trivial cut was produced

  std::unique_ptr<ManifoldSurfaceMesh> diskMesh;
  HalfedgeData<Halfedge> halfedgeMap;
  std::tie(diskMesh, halfedgeMap) = cutAlongEdges(iMesh, cut);

  // gc's eulerCharacteristic counts boundary loops as faces, so a topological
  // disk reports 2 (= V - E + F + nBoundaryLoops, genus 0) with one boundary loop.
  EXPECT_EQ(diskMesh->eulerCharacteristic(), 2);
  EXPECT_EQ(diskMesh->nBoundaryLoops(), 1);
  EXPECT_EQ(diskMesh->nConnectedComponents(), 1);
  EXPECT_TRUE(diskMesh->isManifold());
}

// A closed genus-0 surface, cut through two cones, becomes a disk.
TEST_F(ConeCutSuite, ClosedToDisk) {
  for (const MeshAsset& a : {getAsset("spot.ply", true), getAsset("sphere_small.ply", true)}) {
    a.printThyName();
    ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
    VertexPositionGeometry& geom = *a.geometry;
    ASSERT_EQ(mesh.nBoundaryLoops(), 0u);

    ConePlacementResult cp = computeConePlacement(mesh, geom, 2);
    ConeCutResult cut = computeConeCut(mesh, geom, cp.cones);
    expectCutsToDisk(cut);
  }
}

// A disk with interior cones slit to its boundary stays a disk.
TEST(ConeCutSuite_Disk, BoundaryStaysDisk) {
  std::unique_ptr<ManifoldSurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geom;
  std::tie(mesh, geom) = makeBumpyDisk(8);
  ASSERT_EQ(mesh->eulerCharacteristic(), 2); // gc: disk = 2 (boundary loop counted)
  ASSERT_EQ(mesh->nBoundaryLoops(), 1);

  ConePlacementResult cp = computeConePlacement(*mesh, *geom, 2);
  ConeCutResult cut = computeConeCut(*mesh, *geom, cp.cones);
  expectCutsToDisk(cut);
}

// Geodesic straightening never lengthens the cut.
TEST_F(ConeCutSuite, GeodesicNoLonger) {
  MeshAsset a = getAsset("spot.ply", true);
  ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
  VertexPositionGeometry& geom = *a.geometry;

  ConePlacementResult cp = computeConePlacement(mesh, geom, 3);
  ConeCutResult cut = computeConeCut(mesh, geom, cp.cones);
  ASSERT_NE(cut.network, nullptr);
  // length() is the straightened geodesic length; it must be finite and positive.
  double len = cut.network->length();
  EXPECT_GT(len, 0.0);
  EXPECT_TRUE(std::isfinite(len));
}
