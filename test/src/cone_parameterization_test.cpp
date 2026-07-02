#include "geometrycentral/surface/cone_cut.h"
#include "geometrycentral/surface/cone_parameterization.h"
#include "geometrycentral/surface/cone_placement.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "gtest/gtest.h"

#include <cmath>
#include <map>
#include <memory>
#include <set>
#include <vector>

using namespace geometrycentral;
using namespace geometrycentral::surface;

// Build a triangulated, curved square patch -- a topological disk with genuine
// interior Gaussian curvature (a bump), so cone placement is meaningful.
static std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>>
makeBumpyDisk(int n) {
  std::vector<Vector3> pos;
  pos.reserve((n + 1) * (n + 1));
  for (int i = 0; i <= n; i++) {
    for (int j = 0; j <= n; j++) {
      double x = (double)i / n, y = (double)j / n;
      double z = 0.4 * std::exp(-8.0 * ((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5)));
      pos.push_back(Vector3{x, y, z});
    }
  }
  auto idx = [&](int i, int j) { return (size_t)(i * (n + 1) + j); };
  std::vector<std::vector<size_t>> polys;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      polys.push_back({idx(i, j), idx(i + 1, j), idx(i + 1, j + 1)});
      polys.push_back({idx(i, j), idx(i + 1, j + 1), idx(i, j + 1)});
    }
  }
  return makeManifoldSurfaceMeshAndGeometry(polys, pos);
}

static double uvEdgeLength(const VertexData<Vector2>& uvs, Edge e) {
  return (uvs[e.halfedge().tailVertex()] - uvs[e.halfedge().tipVertex()]).norm();
}

static Edge findEdge(Vertex a, Vertex b) {
  for (Halfedge he : a.outgoingHalfedges())
    if (he.tipVertex() == b) return he.edge();
  throw std::runtime_error("findEdge: no edge between the given vertices");
}

// The seamless guarantee: the two sides of every slit have EQUAL length in the
// UV layout, so the cut can be re-sewn without stretch mismatch. Also: a disk
// with finite UVs.
TEST(ConeParamSuite, SeamlessSlitLengths) {
  std::unique_ptr<ManifoldSurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geom;
  std::tie(mesh, geom) = makeBumpyDisk(8);
  ASSERT_EQ(mesh->eulerCharacteristic(), 2); // gc: disk = 2 (boundary loop counted)
  ASSERT_EQ(mesh->nBoundaryLoops(), 1);

  ConeParameterizationResult res = parameterizeBFFwithCones(*mesh, *geom, 2);
  ASSERT_NE(res.cutMesh, nullptr);
  ManifoldSurfaceMesh& cm = *res.cutMesh;

  EXPECT_EQ(cm.eulerCharacteristic(), 2);
  EXPECT_EQ(cm.nBoundaryLoops(), 1);
  for (Vertex v : cm.vertices()) {
    EXPECT_TRUE(std::isfinite(res.uvs[v].x));
    EXPECT_TRUE(std::isfinite(res.uvs[v].y));
  }

  std::map<int, std::vector<Edge>> byId;
  for (Edge e : cm.edges()) {
    if (e.isBoundary() && res.boundarySeamId[e] >= 0) {
      byId[res.boundarySeamId[e]].push_back(e);
    }
  }
  size_t pairs = 0;
  double maxRelDiff = 0.0;
  for (auto& kv : byId) {
    if (kv.second.size() != 2) continue; // a slit contributes two boundary edges
    pairs++;
    double la = uvEdgeLength(res.uvs, kv.second[0]);
    double lb = uvEdgeLength(res.uvs, kv.second[1]);
    double rel = std::abs(la - lb) / std::max(1e-12, 0.5 * (la + lb));
    maxRelDiff = std::max(maxRelDiff, rel);
  }
  EXPECT_GT(pairs, 0u); // cones were actually slit
  EXPECT_LT(maxRelDiff, 1e-6) << "slit sides differ in length (not seamless)";
}

// Determinism.
TEST(ConeParamSuite, Deterministic) {
  std::unique_ptr<ManifoldSurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geom;
  std::tie(mesh, geom) = makeBumpyDisk(8);

  ConeParameterizationResult r1 = parameterizeBFFwithCones(*mesh, *geom, 2);
  ConeParameterizationResult r2 = parameterizeBFFwithCones(*mesh, *geom, 2);
  ASSERT_EQ(r1.cones.size(), r2.cones.size());
  for (size_t i = 0; i < r1.cones.size(); i++)
    EXPECT_EQ(r1.cones[i].getIndex(), r2.cones[i].getIndex());
}

// flattenWithGivenConesAndCut accepts a cut that did NOT come from
// computeConeCut: hand-route a single interior cone to the boundary along a
// straight chain of plain mesh edges (no geodesic straightening at all), and
// confirm the seamless guarantee -- the two UV-space slit legs, paired via
// boundarySeamId, have equal length -- still holds on this external cut.
TEST(ConeParamSuite, ExternalCutEqualSlitLegs) {
  const int n = 8;
  std::unique_ptr<ManifoldSurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geom;
  std::tie(mesh, geom) = makeBumpyDisk(n);
  auto idx = [&](int i, int j) { return (size_t)(i * (n + 1) + j); };

  Vertex cone = mesh->vertex(idx(n / 2, n / 2)); // interior
  ASSERT_FALSE(cone.isBoundary());

  ConePlacementOptions po;
  po.initialCones = {cone};
  po.maxNewCones = 0;
  ConePlacementResult cp = computeConePlacement(*mesh, *geom, po);

  EdgeData<char> cutEdges(*mesh, 0);
  for (int j = n / 2; j < n; j++) {
    cutEdges[findEdge(mesh->vertex(idx(n / 2, j)), mesh->vertex(idx(n / 2, j + 1)))] = 1;
  }

  ConeParameterizationResult res =
      flattenWithGivenConesAndCut(*mesh, *geom, {cone}, cp.coneAngles, cutEdges);
  ASSERT_NE(res.cutMesh, nullptr);
  EXPECT_EQ(res.cutMesh->nBoundaryLoops(), 1);

  std::map<int, std::vector<Edge>> byId;
  for (Edge e : res.cutMesh->edges())
    if (e.isBoundary() && res.boundarySeamId[e] >= 0) byId[res.boundarySeamId[e]].push_back(e);

  size_t pairs = 0;
  for (auto& kv : byId) {
    if (kv.second.size() != 2) continue;
    pairs++;
    double la = uvEdgeLength(res.uvs, kv.second[0]);
    double lb = uvEdgeLength(res.uvs, kv.second[1]);
    EXPECT_NEAR(la, lb, 1e-9) << "slit sides differ in length (not seamless) on an external cut";
  }
  EXPECT_GT(pairs, 0u);
}

// Multiple disjoint slits (several cones, several separate cut trees that
// share no vertex) must work: two interior cones, each hand-routed to a
// DIFFERENT, far-apart point on the boundary along its own straight chain of
// plain mesh edges. Confirms the disk-and-cones-on-boundary precondition
// holds for a genuinely disconnected cut-edge graph, and that both slits get
// their seamless equal-length-legs pairing.
TEST(ConeParamSuite, MultipleDisjointSlitsCut) {
  const int n = 12;
  std::unique_ptr<ManifoldSurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geom;
  std::tie(mesh, geom) = makeBumpyDisk(n);
  auto idx = [&](int i, int j) { return (size_t)(i * (n + 1) + j); };

  Vertex coneA = mesh->vertex(idx(3, 6));
  Vertex coneB = mesh->vertex(idx(9, 6));
  ASSERT_FALSE(coneA.isBoundary());
  ASSERT_FALSE(coneB.isBoundary());

  ConePlacementOptions po;
  po.initialCones = {coneA, coneB};
  po.maxNewCones = 0;
  ConePlacementResult cp = computeConePlacement(*mesh, *geom, po);

  EdgeData<char> cutEdges(*mesh, 0);
  for (int j = 6; j > 0; j--) // coneA -> boundary at j=0
    cutEdges[findEdge(mesh->vertex(idx(3, j)), mesh->vertex(idx(3, j - 1)))] = 1;
  for (int j = 6; j < n; j++) // coneB -> boundary at j=n
    cutEdges[findEdge(mesh->vertex(idx(9, j)), mesh->vertex(idx(9, j + 1)))] = 1;

  ConeParameterizationResult res =
      flattenWithGivenConesAndCut(*mesh, *geom, {coneA, coneB}, cp.coneAngles, cutEdges);
  ASSERT_NE(res.cutMesh, nullptr);
  EXPECT_EQ(res.cutMesh->nBoundaryLoops(), 1);

  std::map<int, std::vector<Edge>> byId;
  for (Edge e : res.cutMesh->edges())
    if (e.isBoundary() && res.boundarySeamId[e] >= 0) byId[res.boundarySeamId[e]].push_back(e);

  size_t pairs = 0;
  for (auto& kv : byId) {
    if (kv.second.size() != 2) continue;
    pairs++;
    double la = uvEdgeLength(res.uvs, kv.second[0]);
    double lb = uvEdgeLength(res.uvs, kv.second[1]);
    EXPECT_NEAR(la, lb, 1e-9) << "slit sides differ in length (not seamless) on a disjoint external cut";
  }
  // Each cut EDGE (not each cone) gets its own seamId pair: 6 edges per slit,
  // 2 disjoint slits.
  EXPECT_EQ(pairs, 12u);
}

// parameterizeBFFwithCones(mesh, geo, n) must be exactly equivalent to its
// documented decomposition: placement + computeConeCut + a manual transfer to
// the intrinsic mesh + flattenWithGivenConesAndCut. A regression guard for the
// wrapper -- if a future edit makes it diverge from that decomposition (wrong
// argument, stale field, etc.), this catches it via non-identical UVs.
TEST(ConeParamSuite, DelegationEquivalence) {
  std::unique_ptr<ManifoldSurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geom;
  std::tie(mesh, geom) = makeBumpyDisk(8);

  ConeParameterizationResult viaWrapper = parameterizeBFFwithCones(*mesh, *geom, 2);

  ConePlacementResult cp = computeConePlacement(*mesh, *geom, 2);
  ConeCutResult cc = computeConeCut(*mesh, *geom, cp.cones);
  ManifoldSurfaceMesh& iMesh = cc.network->mesh;
  IntrinsicGeometryInterface& iGeo = *cc.network->tri;

  std::vector<Vertex> iCones;
  VertexData<double> coneAngleI(iMesh, 0.0);
  for (Vertex c : cp.cones) {
    Vertex ic = iMesh.vertex(c.getIndex());
    iCones.push_back(ic);
    coneAngleI[ic] = cp.coneAngles[c];
  }
  EdgeData<char> cutChar(iMesh, 0);
  for (Edge e : iMesh.edges())
    if (cc.cutEdges[e]) cutChar[e] = 1;

  ConeParameterizationResult viaManual = flattenWithGivenConesAndCut(iMesh, iGeo, iCones, coneAngleI, cutChar);

  ASSERT_EQ(viaWrapper.cutMesh->nVertices(), viaManual.cutMesh->nVertices());
  for (size_t i = 0; i < viaWrapper.cutMesh->nVertices(); i++) {
    Vector2 a = viaWrapper.uvs[viaWrapper.cutMesh->vertex(i)];
    Vector2 b = viaManual.uvs[viaManual.cutMesh->vertex(i)];
    EXPECT_EQ(a.x, b.x);
    EXPECT_EQ(a.y, b.y);
  }
}

// Every interior halfedge of cutMesh maps (via parentHalfedge) to an input
// halfedge with the same intrinsic length; the two cut copies of a slit edge
// map to the SAME input edge (that's what makes them a pair at all).
// Exercises flattenWithGivenConesAndCut directly (not the parameterizeBFFwithCones
// wrapper): parentHalfedge maps back to THAT call's `mesh` argument, which is
// only directly usable by a caller who owns it -- exactly pendel's situation,
// passing its own intrinsic spine. (Via the wrapper, parentHalfedge maps back
// to the wrapper's internally-built, unexposed intrinsic mesh instead -- see
// the field's doc comment.)
TEST(ConeParamSuite, ParentHalfedgeRoundTrip) {
  std::unique_ptr<ManifoldSurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geom;
  std::tie(mesh, geom) = makeBumpyDisk(8);

  // Place cones and cut on `*mesh` directly, the same mesh flattenWithGivenConesAndCut
  // is called on below -- so parentHalfedge maps back to a mesh this test owns
  // and can check lengths against.
  ConePlacementResult cp = computeConePlacement(*mesh, *geom, 2);
  ConeCutResult cc = computeConeCut(*mesh, *geom, cp.cones);
  ManifoldSurfaceMesh& iMesh = cc.network->mesh;
  IntrinsicGeometryInterface& iGeo = *cc.network->tri;
  iGeo.requireEdgeLengths();

  std::vector<Vertex> iCones;
  VertexData<double> coneAngleI(iMesh, 0.0);
  for (Vertex c : cp.cones) {
    Vertex ic = iMesh.vertex(c.getIndex());
    iCones.push_back(ic);
    coneAngleI[ic] = cp.coneAngles[c];
  }
  EdgeData<char> cutChar(iMesh, 0);
  for (Edge e : iMesh.edges())
    if (cc.cutEdges[e]) cutChar[e] = 1;

  ConeParameterizationResult res = flattenWithGivenConesAndCut(iMesh, iGeo, iCones, coneAngleI, cutChar);
  ManifoldSurfaceMesh& cm = *res.cutMesh;

  for (Halfedge he : cm.halfedges()) {
    Halfedge parent = res.parentHalfedge[he];
    if (he.isInterior()) {
      ASSERT_NE(parent, Halfedge()) << "interior cutMesh halfedge has no parent";
      EXPECT_NEAR(res.cutGeometry->edgeLengths[he.edge()], iGeo.edgeLengths[parent.edge()], 1e-9);
    }
  }

  // Slit pairs: both boundary edges sharing a seam id must trace back to the
  // same iMesh edge.
  std::map<int, std::set<size_t>> origEdgeOfSeam;
  for (Edge ce : cm.edges()) {
    if (!ce.isBoundary() || res.boundarySeamId[ce] < 0) continue;
    Halfedge cutHe = ce.halfedge();
    Halfedge parent = cutHe.isInterior() ? res.parentHalfedge[cutHe] : res.parentHalfedge[cutHe.twin()];
    ASSERT_NE(parent, Halfedge());
    origEdgeOfSeam[res.boundarySeamId[ce]].insert(parent.edge().getIndex());
  }
  for (auto& kv : origEdgeOfSeam) EXPECT_EQ(kv.second.size(), 1u) << "seam id " << kv.first;
}
