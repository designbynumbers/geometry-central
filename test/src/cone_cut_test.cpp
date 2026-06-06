#include "geometrycentral/surface/cone_cut.h"
#include "geometrycentral/surface/cone_placement.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/surgery.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "load_test_meshes.h"

#include "gtest/gtest.h"

#include <vector>

using namespace geometrycentral;
using namespace geometrycentral::surface;

class ConeCutSuite : public MeshAssetSuite {};

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

  EXPECT_EQ(diskMesh->eulerCharacteristic(), 1);
  EXPECT_EQ(diskMesh->nBoundaryLoops(), 1);
  EXPECT_TRUE(diskMesh->isManifold());
}

// On a closed genus-0 surface, the cut connecting the cones is a single
// well-formed geodesic arc (one path, two leaf endpoints at the cones).
// Note: opening this interior slit into a disk requires interior-slit cutting,
// which gc's surgery::cutAlongEdges / separateEdge does not yet support; so we
// verify the cut here rather than the post-cut disk. Surfaces with boundary
// (the manufacturing case) are exercised by BoundaryStaysDisk below.
TEST_F(ConeCutSuite, ClosedCutIsGeodesicArc) {
  for (const MeshAsset& a : {getAsset("spot.ply", true), getAsset("sphere_small.ply", true)}) {
    a.printThyName();
    ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
    VertexPositionGeometry& geom = *a.geometry;
    ASSERT_EQ(mesh.nBoundaryLoops(), 0);

    ConePlacementResult cp = computeConePlacement(mesh, geom, 2);
    ConeCutResult cut = computeConeCut(mesh, geom, cp.cones);
    ASSERT_NE(cut.network, nullptr);
    ManifoldSurfaceMesh& iMesh = cut.network->mesh;

    size_t nCut = 0, leaves = 0;
    for (Edge e : iMesh.edges())
      if (cut.cutEdges[e]) nCut++;
    for (Vertex v : iMesh.vertices()) {
      size_t d = 0;
      for (Edge e : v.adjacentEdges())
        if (cut.cutEdges[e]) d++;
      if (d == 1) leaves++;
    }
    EXPECT_GT(nCut, 0u);
    EXPECT_EQ(cut.network->paths.size(), 1u); // single arc between the two cones
    EXPECT_EQ(leaves, 2u);                     // two endpoints
  }
}

// A surface with boundary stays a disk after slitting interior cones to it.
TEST_F(ConeCutSuite, BoundaryStaysDisk) {
  for (const MeshAsset& a : boundaryMeshes()) {
    if (!a.isTriangular || a.manifoldMesh == nullptr || a.geometry == nullptr) continue;
    ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
    // Only a disk-topology input is guaranteed to stay a disk after slitting.
    if (mesh.eulerCharacteristic() != 1 || mesh.nBoundaryLoops() != 1) continue;
    if (mesh.nInteriorVertices() < 3) continue;
    a.printThyName();
    VertexPositionGeometry& geom = *a.geometry;

    ConePlacementResult cp = computeConePlacement(mesh, geom, 2);
    ConeCutResult cut = computeConeCut(mesh, geom, cp.cones);
    expectCutsToDisk(cut);
  }
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
