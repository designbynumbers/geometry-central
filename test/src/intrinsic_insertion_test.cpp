// Tests for arbitrary vertex insertion/deletion on
// IntegerCoordinatesIntrinsicTriangulation: targeted regression tests plus
// randomized stress tests which validate both combinatorial consistency
// (common subdivision construction) and geometric consistency (common
// subdivision edge lengths computed three independent ways).
#include "geometrycentral/surface/integer_coordinates_intrinsic_triangulation.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "load_test_meshes.h"

#include "gtest/gtest.h"

#include <random>
#include <sstream>

using namespace geometrycentral;
using namespace geometrycentral::surface;

class IntrinsicInsertionSuite : public MeshAssetSuite {};

namespace {

// Validate the common subdivision, combinatorially and geometrically.
// Returns an error string ("" if ok).
// NOTE: invalidates element handles (compresses the mesh)
std::string checkCS(IntegerCoordinatesIntrinsicTriangulation& tri, VertexPositionGeometry& origGeometry) {
  try {
    CommonSubdivision& cs = tri.getCommonSubdivision();
    cs.constructMesh();

    // Lengths from extrinsic vertex positions
    const VertexData<Vector3>& posCS = cs.interpolateAcrossA(origGeometry.vertexPositions);
    VertexPositionGeometry csGeo(*cs.mesh, posCS);
    csGeo.requireEdgeLengths();
    EdgeData<double> lengthsFromPosA = csGeo.edgeLengths;
    csGeo.unrequireEdgeLengths();

    // Lengths from input edge lengths
    origGeometry.requireEdgeLengths();
    EdgeData<double> lengthsFromLenA = cs.interpolateEdgeLengthsA(origGeometry.edgeLengths);

    // Lengths from intrinsic edge lengths
    EdgeData<double> lengthsFromLenB = cs.interpolateEdgeLengthsB(tri.edgeLengths);

    double scale = lengthsFromPosA.toVector().norm();
    double errAB = (lengthsFromPosA.toVector() - lengthsFromLenB.toVector()).norm() / scale;
    double errAA = (lengthsFromPosA.toVector() - lengthsFromLenA.toVector()).norm() / scale;

    // The CS mesh is rebuilt each time, so drop it to keep memory in check
    if (errAA > 1e-5 || errAB > 1e-5) {
      std::ostringstream err;
      err << "CS geometry mismatch: relative errAA = " << errAA << ", errAB = " << errAB;
      return err.str();
    }
  } catch (std::exception& e) {
    return std::string("CS validation failed: ") + e.what();
  }
  return "";
}

// Sample the i'th *live* face/edge/vertex (mesh may be uncompressed)
template <typename SetType>
auto sampleNth(const SetType& set, size_t n, size_t count) -> decltype(*set.begin()) {
  size_t i = n % count;
  for (auto elem : set) {
    if (i == 0) return elem;
    i--;
  }
  throw std::runtime_error("sampleNth ran off the end");
}

void probeMesh(const MeshAsset& a, unsigned seed, int nOps, int checkEvery = 50) {
  a.printThyName();
  ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
  VertexPositionGeometry& origGeometry = *a.geometry;

  IntegerCoordinatesIntrinsicTriangulation tri(mesh, origGeometry);
  tri.flipToDelaunay(); // make sure edges have crossings

  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> unit(0.05, 0.95);
  std::uniform_int_distribution<int> opDist(0, 9);
  std::uniform_int_distribution<size_t> idxDist(0, 1u << 30);

  ManifoldSurfaceMesh& im = *tri.intrinsicMesh;

  int opCount = 0;
  for (int iOp = 0; iOp < nOps; iOp++) {
    int op = opDist(rng);
    std::ostringstream desc;
    try {
      if (op < 4) {
        // Face split at random interior point
        Face f = sampleNth(im.faces(), idxDist(rng), im.nFaces());
        double b0 = unit(rng), b1 = unit(rng), b2 = unit(rng);
        Vector3 bary = Vector3{b0, b1, b2} / (b0 + b1 + b2);
        desc << "splitFace(" << f << ", " << bary << ")";
        tri.splitFace(f, bary);
      } else if (op < 7) {
        // Edge split at random point
        Edge e = sampleNth(im.edges(), idxDist(rng), im.nEdges());
        double t = unit(rng);
        desc << "splitEdge(" << e << " [n=" << tri.normalCoordinates[e] << ", bdy=" << e.isBoundary() << "], " << t
             << ")";
        tri.splitEdge(e, t);
      } else if (op < 9) {
        // Remove a random inserted vertex (vertex whose location is not a
        // vertex of the input mesh)
        Vertex v;
        size_t start = idxDist(rng);
        size_t nV = im.nVertices();
        for (size_t probe = 0; probe < nV; probe++) {
          Vertex cand = sampleNth(im.vertices(), start + probe, nV);
          if (tri.vertexLocations[cand].type != SurfacePointType::Vertex) {
            v = cand;
            break;
          }
        }
        if (v == Vertex()) continue; // nothing to remove
        desc << "removeInsertedVertex(" << v << " [at " << tri.vertexLocations[v] << ", bdy=" << v.isBoundary()
             << "])";
        tri.removeInsertedVertex(v);
      } else {
        // insertVertex at a random edge point (tests vertex-type dispatch)
        Edge e = sampleNth(im.edges(), idxDist(rng), im.nEdges());
        double t = unit(rng);
        desc << "insertVertex(EdgePoint " << e << " [n=" << tri.normalCoordinates[e] << ", bdy=" << e.isBoundary()
             << "], " << t << ")";
        tri.insertVertex(SurfacePoint(e, t));
      }
      opCount++;
    } catch (std::exception& err) {
      ADD_FAILURE() << "op " << iOp << ": " << desc.str() << " threw: " << err.what();
      return;
    }

    // Periodically validate
    if (iOp % checkEvery == checkEvery - 1) {
      std::string csErr = checkCS(tri, origGeometry);
      if (!csErr.empty()) {
        ADD_FAILURE() << "after op " << iOp << " (" << desc.str() << "): " << csErr;
        return;
      }
    }
  }
  std::cout << "  completed " << opCount << " ops" << std::endl;
  std::string csErr = checkCS(tri, origGeometry);
  EXPECT_EQ(csErr, "");
}

} // namespace

// Regression test: splitting a crossed edge exactly at one of its
// input-edge crossing parameters must record the new vertex as a point ON
// that input edge (Edge-typed location), not as a face point with a
// degenerate (~zero) barycentric component. Such degenerate face points
// break consumers which dispatch on the location type (curve tracing,
// identifyInputCurveRange, common subdivision construction); found by a
// downstream project whose insertion parameters come from the common
// subdivision itself.
TEST_F(IntrinsicInsertionSuite, SplitEdgeAtCrossingParameter) {
  auto a = getAsset("spot.ply", true);
  ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
  VertexPositionGeometry& origGeometry = *a.geometry;

  IntegerCoordinatesIntrinsicTriangulation tri(mesh, origGeometry);
  tri.delaunayRefine(25.0);
  tri.intrinsicMesh->compress();

  // Read each crossed edge's first transverse crossing parameter off the
  // common subdivision
  struct Target {
    Edge e;
    double t;
  };
  std::vector<Target> targets;
  {
    CommonSubdivision& cs = tri.getCommonSubdivision();
    for (Edge e : tri.intrinsicMesh->edges()) {
      if (tri.normalCoordinates[e] <= 0) continue;
      for (CommonSubdivisionPoint* p : cs.pointsAlongB[e]) {
        if (p->intersectionType != CSIntersectionType::EDGE_TRANSVERSE) continue;
        if (p->posB.type != SurfacePointType::Edge) continue;
        targets.push_back({e, p->posB.tEdge});
        break;
      }
      if (targets.size() >= 25) break;
    }
  }
  ASSERT_GT(targets.size(), 0u);

  // Insert exactly at the crossing parameters; every resulting location must
  // name the element it lies on (in particular, no face point may have a
  // numerically-zero barycentric component)
  for (const Target& tg : targets) {
    if (tg.e.isDead()) continue;
    Vertex v = tri.insertVertex(SurfacePoint(tg.e, tg.t));
    ASSERT_NE(v, Vertex());
    SurfacePoint loc = tri.vertexLocations[v];
    if (loc.type == SurfacePointType::Face) {
      double minBary = std::min(loc.faceCoords.x, std::min(loc.faceCoords.y, loc.faceCoords.z));
      EXPECT_GT(minBary, 1e-9) << "degenerate face-point location " << loc;
    } else if (loc.type == SurfacePointType::Edge) {
      EXPECT_GE(loc.tEdge, 0.);
      EXPECT_LE(loc.tEdge, 1.);
    }
  }

  EXPECT_EQ(checkCS(tri, origGeometry), "");
}

// Every vertex inserted during a (possibly seam-constrained) delaunayRefine
// must be recorded with a location that names the element it lies on: no
// Face-typed location may carry a numerically-zero barycentric component.
// Uses the insertion callbacks so every vertex is checked at creation time.
TEST_F(IntrinsicInsertionSuite, RefineRecordsNoDegenerateLocations) {
  for (const MeshAsset& a : {getAsset("spot.ply", true), getAsset("lego.ply", true)}) {
    a.printThyName();
    ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
    VertexPositionGeometry& origGeometry = *a.geometry;

    IntegerCoordinatesIntrinsicTriangulation tri(mesh, origGeometry);

    // Mark some edges to simulate constrained (seam) refinement: a few
    // edge paths walked from arbitrary vertices
    EdgeData<bool> marked(*tri.intrinsicMesh, false);
    int nMarked = 0;
    for (size_t iV = 0; iV < tri.intrinsicMesh->nVertices(); iV += 37) {
      Halfedge he = tri.intrinsicMesh->vertex(iV).halfedge();
      for (int step = 0; step < 8; step++) {
        if (he.edge().isBoundary()) break;
        marked[he.edge()] = true;
        nMarked++;
        he = he.next().next().twin().next(); // wander
      }
    }
    ASSERT_GT(nMarked, 0);
    tri.setMarkedEdges(marked);

    // Check every inserted vertex's location at creation time
    int nChecked = 0, nDegenerate = 0;
    auto checkLoc = [&](Vertex v) {
      nChecked++;
      SurfacePoint loc = tri.vertexLocations[v];
      if (loc.type == SurfacePointType::Face) {
        double minBary = std::min(loc.faceCoords.x, std::min(loc.faceCoords.y, loc.faceCoords.z));
        if (minBary < 1e-12) {
          nDegenerate++;
          ADD_FAILURE() << "degenerate location recorded at creation: " << loc;
        }
      }
    };
    tri.faceInsertionCallbackList.push_back([&](Face f, Vertex v) { checkLoc(v); });
    tri.edgeSplitCallbackList.push_back([&](Edge e, Halfedge he1, Halfedge he2) { checkLoc(he1.vertex()); });

    tri.delaunayRefine(25.);

    std::cout << "  checked " << nChecked << " insertions (" << nDegenerate << " degenerate)" << std::endl;
    EXPECT_GT(nChecked, 0);
    EXPECT_EQ(nDegenerate, 0);
    EXPECT_EQ(checkCS(tri, origGeometry), "");
  }
}

// Unit trigger from the round-2 downstream report: a face split whose
// barycentric coordinate has a ~1e-18 component (i.e. a point numerically
// on an intrinsic edge of a crossed face) must still record a
// properly-typed location.
TEST_F(IntrinsicInsertionSuite, SplitFaceAtNearEdgeBary) {
  auto a = getAsset("spot.ply", true);
  ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
  VertexPositionGeometry& origGeometry = *a.geometry;

  IntegerCoordinatesIntrinsicTriangulation tri(mesh, origGeometry);
  tri.delaunayRefine(25.0);
  tri.intrinsicMesh->compress();

  int nTested = 0;
  for (size_t iF = 0; iF < tri.intrinsicMesh->nFaces() && nTested < 15; iF++) {
    Face f = tri.intrinsicMesh->face(iF);
    bool crossed = false;
    for (Edge e : f.adjacentEdges()) crossed = crossed || (tri.normalCoordinates[e] > 0);
    if (!crossed) continue;

    Vector3 bary{1.04494e-18, 0.41, 1. - 0.41 - 1.04494e-18};
    Vertex v = tri.splitFace(f, bary);
    ASSERT_NE(v, Vertex());
    SurfacePoint loc = tri.vertexLocations[v];
    if (loc.type == SurfacePointType::Face) {
      double minBary = std::min(loc.faceCoords.x, std::min(loc.faceCoords.y, loc.faceCoords.z));
      EXPECT_GT(minBary, 1e-12) << "degenerate location " << loc;
    }
    nTested++;
  }
  EXPECT_GT(nTested, 0);
  EXPECT_EQ(checkCS(tri, origGeometry), "");
}

// Regression test: splitting an interior edge which has crossings (n > 0)
// must produce valid (in particular, nonnegative) normal coordinates on the
// new cross edges. An unsigned-arithmetic bug used to wrap these to -1,
// marking the cross edges as "running along an input edge" and corrupting
// the correspondence.
TEST_F(IntrinsicInsertionSuite, SplitEdgeWithCrossings) {
  auto a = getAsset("fox.ply", true);
  ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
  VertexPositionGeometry& origGeometry = *a.geometry;

  IntegerCoordinatesIntrinsicTriangulation tri(mesh, origGeometry);
  tri.flipToDelaunay();
  ManifoldSurfaceMesh& im = *tri.intrinsicMesh;

  int nSplit = 0;
  for (size_t iE = 0; iE < im.nEdges() && nSplit < 25; iE++) {
    Edge e = im.edge(iE);
    if (tri.normalCoordinates[e] <= 0) continue;

    Vertex v = tri.splitEdge(e, 0.5);
    nSplit++;
    ASSERT_NE(v, Vertex());
    for (Edge ve : v.adjacentEdges()) {
      // a freshly inserted vertex is not on any input edge or vertex, so
      // none of its edges may be shared with the input mesh
      EXPECT_GE(tri.normalCoordinates[ve], 0);
    }
  }
  EXPECT_GT(nSplit, 0);
  EXPECT_EQ(checkCS(tri, origGeometry), "");
}

// Regression test: splitting a shared edge (n = -1) leaves inserted vertices
// along the input edge; faces incident on such vertices have curves
// emanating from a corner (half-integer corner coordinates), and splitting
// those faces or their edges must still work.
TEST_F(IntrinsicInsertionSuite, SplitAdjacentToSharedSplit) {
  auto a = getAsset("fox.ply", true);
  ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
  VertexPositionGeometry& origGeometry = *a.geometry;

  IntegerCoordinatesIntrinsicTriangulation tri(mesh, origGeometry);
  tri.flipToDelaunay();
  ManifoldSurfaceMesh& im = *tri.intrinsicMesh;

  int nTested = 0;
  for (size_t iE = 0; iE < im.nEdges() && nTested < 10; iE++) {
    Edge e = im.edge(iE);
    if (tri.normalCoordinates[e] != -1) continue;

    // Split the shared edge, leaving an inserted vertex on the input edge
    Vertex v = tri.splitEdge(e, 0.75);
    ASSERT_NE(v, Vertex());

    // Split an edge of v which has crossings, and split a neighboring face;
    // these exercise the emanating-curve cases
    for (Edge ve : v.adjacentEdges()) {
      if (tri.normalCoordinates[ve] > 0) {
        tri.splitEdge(ve, 0.5);
        break;
      }
    }
    for (Face f : v.adjacentFaces()) {
      tri.splitFace(f, Vector3{0.3, 0.3, 0.4});
      break;
    }
    nTested++;
  }
  EXPECT_GT(nTested, 0);
  EXPECT_EQ(checkCS(tri, origGeometry), "");
}

// Insert vertices on boundary edges, then remove them again
TEST_F(IntrinsicInsertionSuite, BoundaryInsertDelete) {
  auto a = getAsset("lego.ply", true);
  ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
  VertexPositionGeometry& origGeometry = *a.geometry;

  IntegerCoordinatesIntrinsicTriangulation tri(mesh, origGeometry);
  tri.flipToDelaunay();
  ManifoldSurfaceMesh& im = *tri.intrinsicMesh;

  // Find boundary edges and split them
  std::vector<Vertex> insertedBoundaryVertices;
  for (size_t iE = 0; iE < im.nEdges(); iE++) {
    Edge e = im.edge(iE);
    if (!e.isBoundary()) continue;
    Vertex v = tri.splitEdge(e, 0.4);
    ASSERT_NE(v, Vertex());
    EXPECT_TRUE(v.isBoundary());
    EXPECT_EQ(tri.vertexLocations[v].type, SurfacePointType::Edge);
    insertedBoundaryVertices.push_back(v);
    if (insertedBoundaryVertices.size() >= 20) break;
  }
  ASSERT_GT(insertedBoundaryVertices.size(), 0u);

  // Remove them again (in reverse order, just to mix things up)
  int nRemoved = 0;
  for (auto it = insertedBoundaryVertices.rbegin(); it != insertedBoundaryVertices.rend(); ++it) {
    Face f = tri.removeInsertedVertex(*it);
    if (f != Face()) nRemoved++;
  }
  EXPECT_EQ(nRemoved, (int)insertedBoundaryVertices.size());
  EXPECT_EQ(checkCS(tri, origGeometry), "");
}

// Regression test: tracing an input edge through an inserted vertex must
// scan the vertex's whole fan of corners for the curve's continuation.
// traceNextCurve used to throw on the first corner with no crossing, which
// fails whenever the shared sub-edges left by a shared-edge split are
// later flipped away (e.g. by flipToDelaunay), leaving an Edge-typed vertex
// with only transverse continuations.
TEST_F(IntrinsicInsertionSuite, TraceThroughInsertedVertexAfterFlips) {
  auto a = getAsset("fox.ply", true);
  ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
  VertexPositionGeometry& origGeometry = *a.geometry;

  IntegerCoordinatesIntrinsicTriangulation tri(mesh, origGeometry);
  tri.flipToDelaunay();
  ManifoldSurfaceMesh& im = *tri.intrinsicMesh;

  // Split shared edges, then flip the two shared sub-edges away so that the
  // input edge continues through the inserted vertex transversally on both
  // sides. (Splitting the surrounding faces first adds edges at the inserted
  // vertex, which makes the shared sub-edges geometrically flippable.)
  int nIsolated = 0;
  size_t nOrigEdges = im.nEdges();
  for (size_t iE = 0; iE < nOrigEdges; iE++) {
    Edge e = im.edge(iE);
    if (tri.normalCoordinates[e] != -1 || e.isBoundary()) continue;

    Vertex v = tri.splitEdge(e, 0.5);
    ASSERT_NE(v, Vertex());

    std::vector<Face> neighborhood;
    for (Face f : v.adjacentFaces()) neighborhood.push_back(f);
    for (Face f : neighborhood) tri.splitFace(f, Vector3{1. / 3., 1. / 3., 1. / 3.});

    std::vector<Edge> sharedSubEdges;
    for (Edge ve : v.adjacentEdges()) {
      if (tri.normalCoordinates[ve] < 0) sharedSubEdges.push_back(ve);
    }
    for (Edge se : sharedSubEdges) tri.flipEdgeIfPossible(se);

    bool anyShared = false;
    for (Edge ve : v.adjacentEdges()) anyShared = anyShared || (tri.normalCoordinates[ve] < 0);
    if (!anyShared) {
      nIsolated++;
      // Directly trace the input edge through the isolated vertex, in both
      // directions; the trace must continue through it (multiple components)
      Edge inputE = tri.vertexLocations[v].edge;
      NormalCoordinatesCompoundCurve cc = tri.traceInputHalfedge(inputE.halfedge());
      EXPECT_GE(cc.components.size(), 2u);
      NormalCoordinatesCompoundCurve ccRev = tri.traceInputHalfedge(inputE.halfedge().twin());
      EXPECT_GE(ccRev.components.size(), 2u);
    }

    if (nIsolated >= 5) break;
  }
  EXPECT_GT(nIsolated, 0); // at least one vertex is in the regression configuration
  EXPECT_EQ(checkCS(tri, origGeometry), "");
}

TEST_F(IntrinsicInsertionSuite, ClosedMeshProbe) {
  probeMesh(getAsset("fox.ply", true), 7, 500);
}

TEST_F(IntrinsicInsertionSuite, ClosedMeshProbe2) {
  probeMesh(getAsset("sphere_small.ply", true), 11, 500);
}

TEST_F(IntrinsicInsertionSuite, ClosedMeshProbe3) {
  probeMesh(getAsset("bob_small.ply", true), 17, 500);
}

TEST_F(IntrinsicInsertionSuite, TinyMeshProbe) {
  probeMesh(getAsset("tet.obj", true), 19, 300);
}

TEST_F(IntrinsicInsertionSuite, BoundaryMeshProbe) {
  probeMesh(getAsset("lego.ply", true), 13, 500);
}

TEST_F(IntrinsicInsertionSuite, BoundaryMeshProbe2) {
  probeMesh(getAsset("cat_head.obj", true), 23, 500);
}

TEST_F(IntrinsicInsertionSuite, ClosedMeshProbe4) {
  probeMesh(getAsset("spot.ply", true), 29, 300);
}

TEST_F(IntrinsicInsertionSuite, ManySeedsFox) {
  for (unsigned seed = 100; seed < 110; seed++) {
    probeMesh(getAsset("fox.ply", true), seed, 200, 100);
  }
}

TEST_F(IntrinsicInsertionSuite, ManySeedsLego) {
  for (unsigned seed = 200; seed < 210; seed++) {
    probeMesh(getAsset("lego.ply", true), seed, 200, 100);
  }
}
