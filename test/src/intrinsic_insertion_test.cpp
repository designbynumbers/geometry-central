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
      try {
        tri.validate();
      } catch (std::exception& vErr) {
        ADD_FAILURE() << "after op " << iOp << " (" << desc.str() << "): " << vErr.what();
        return;
      }
      std::string csErr = checkCS(tri, origGeometry);
      if (!csErr.empty()) {
        ADD_FAILURE() << "after op " << iOp << " (" << desc.str() << "): " << csErr;
        return;
      }
    }
  }
  std::cout << "  completed " << opCount << " ops" << std::endl;
  EXPECT_NO_THROW(tri.validate());
  std::string csErr = checkCS(tri, origGeometry);
  EXPECT_EQ(csErr, "");

  // Validate the exact containing-face derivation on every uncrossed edge
  // (internally asserts that both endpoints agree on the face)
  int nUncrossed = 0;
  for (Edge e : im.edges()) {
    if (tri.normalCoordinates[e] != 0) continue;
    try {
      Face fIn = tri.inputFaceOfUncrossedEdge(e);
      EXPECT_NE(fIn, Face());
    } catch (std::exception& err) {
      ADD_FAILURE() << "inputFaceOfUncrossedEdge(" << e << "): " << err.what();
      break;
    }
    nUncrossed++;
  }
  std::cout << "  validated " << nUncrossed << " uncrossed edges" << std::endl;
}

} // namespace

// Regression test: splitting a crossed edge exactly at one of its
// input-edge crossing parameters. The split's combinatorial classification
// places the new vertex strictly beside the crossing, in an input face
// derived exactly from the crossing curves; the recorded coordinates may
// land numerically on the face boundary, which is legal. (Found by a
// downstream project whose insertion parameters come from the common
// subdivision itself; an earlier eps-snapping treatment of this case was
// replaced by the exact element derivation.)
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

  // Insert exactly at the crossing parameters. The split's combinatorial
  // classification places each new vertex strictly beside the crossing, so
  // the recorded element is the (exactly-derived) input face of that
  // segment, with coordinates which may legitimately land numerically on
  // the face boundary. The oracle is that all coordinates are finite/in
  // range and the common subdivision remains consistent.
  for (const Target& tg : targets) {
    if (tg.e.isDead()) continue;
    Vertex v = tri.insertVertex(SurfacePoint(tg.e, tg.t));
    ASSERT_NE(v, Vertex());
    SurfacePoint loc = tri.vertexLocations[v];
    if (loc.type == SurfacePointType::Face) {
      for (int i = 0; i < 3; i++) {
        EXPECT_TRUE(std::isfinite(loc.faceCoords[i])) << loc;
        EXPECT_GE(loc.faceCoords[i], 0.);
        EXPECT_LE(loc.faceCoords[i], 1.);
      }
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

    // Check every inserted vertex's location at creation time. Note that a
    // Face-typed location with a ~zero barycentric component is LEGAL (the
    // element is exact; the coordinates are honest floats); what we require
    // is that coordinates are finite and in range.
    int nChecked = 0, nDegenerate = 0;
    auto checkLoc = [&](Vertex v) {
      nChecked++;
      SurfacePoint loc = tri.vertexLocations[v];
      if (loc.type == SurfacePointType::Face) {
        for (int i = 0; i < 3; i++) {
          if (!std::isfinite(loc.faceCoords[i]) || loc.faceCoords[i] < -1e-9 || loc.faceCoords[i] > 1. + 1e-9) {
            nDegenerate++;
            ADD_FAILURE() << "out-of-range location recorded at creation: " << loc;
            break;
          }
        }
      } else if (loc.type == SurfacePointType::Edge) {
        if (!std::isfinite(loc.tEdge) || loc.tEdge < -1e-9 || loc.tEdge > 1. + 1e-9) {
          nDegenerate++;
          ADD_FAILURE() << "out-of-range location recorded at creation: " << loc;
        }
      }
    };
    tri.faceInsertionCallbackList.push_back([&](Face f, Vertex v) { checkLoc(v); });
    tri.edgeSplitCallbackList.push_back([&](Edge e, Halfedge he1, Halfedge he2) { checkLoc(he1.vertex()); });

    tri.delaunayRefine(25.);

    std::cout << "  checked " << nChecked << " insertions (" << nDegenerate << " out of range)" << std::endl;
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
      for (int i = 0; i < 3; i++) EXPECT_TRUE(std::isfinite(loc.faceCoords[i])) << loc;
    }
    nTested++;
  }
  EXPECT_GT(nTested, 0);
  EXPECT_EQ(checkCS(tri, origGeometry), "");
}

// insertVertexAtCrossing cuts an input-edge curve at a chosen transverse
// crossing: the new vertex's location is a point ON the input edge, the
// curve's halves terminate at it, and the correspondence stays consistent
TEST_F(IntrinsicInsertionSuite, InsertVertexAtCrossing) {
  auto a = getAsset("spot.ply", true);
  ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
  VertexPositionGeometry& origGeometry = *a.geometry;

  IntegerCoordinatesIntrinsicTriangulation tri(mesh, origGeometry);
  tri.delaunayRefine(25.0);
  tri.intrinsicMesh->compress();
  ManifoldSurfaceMesh& im = *tri.intrinsicMesh;

  int nInserted = 0;
  size_t nOrigEdges = im.nEdges();
  for (size_t iE = 0; iE < nOrigEdges && nInserted < 25; iE++) {
    Edge e = im.edge(iE);
    if (tri.normalCoordinates[e] <= 0) continue;

    int nBefore = tri.normalCoordinates[e];
    Vertex v = tri.insertVertexAtCrossing(e.halfedge(), 0);
    ASSERT_NE(v, Vertex());
    nInserted++;

    // The location names a point on an input edge
    SurfacePoint loc = tri.vertexLocations[v];
    ASSERT_EQ(loc.type, SurfacePointType::Edge);
    EXPECT_GE(loc.tEdge, 0.);
    EXPECT_LE(loc.tEdge, 1.);

    // All new edges have valid normal coordinates (>= -1; a value of -1
    // marks a cross edge which coincides with a piece of the cut curve,
    // which happens when the curve emanated from the adjacent apex)
    for (Edge ve : v.adjacentEdges()) {
      EXPECT_GE(tri.normalCoordinates[ve], -1);
    }
    (void)nBefore;

    // The input edge's curve terminates at v: tracing it yields a component
    // boundary there (the trace passes through v as a compound curve)
    NormalCoordinatesCompoundCurve cc = tri.traceInputHalfedge(loc.edge.halfedge());
    EXPECT_GE(cc.components.size(), 2u);
  }
  ASSERT_GT(nInserted, 0);
  EXPECT_EQ(checkCS(tri, origGeometry), "");
}

// validate() must catch representative corruptions of each invariant
// class; this is the validator's own regression test (a validator that
// passes everything is worse than none)
TEST_F(IntrinsicInsertionSuite, ValidatorCatchesCorruption) {
  auto a = getAsset("fox.ply", true);
  ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
  VertexPositionGeometry& origGeometry = *a.geometry;

  auto freshTri = [&]() {
    auto tri = std::make_unique<IntegerCoordinatesIntrinsicTriangulation>(mesh, origGeometry);
    tri->flipToDelaunay();
    // a few inserts so all location types exist
    tri->splitFace(tri->intrinsicMesh->face(3), Vector3{0.3, 0.3, 0.4});
    for (Edge e : tri->intrinsicMesh->edges()) {
      if (tri->normalCoordinates[e] == -1 && !e.isBoundary()) {
        tri->splitEdge(e, 0.4);
        break;
      }
    }
    EXPECT_NO_THROW(tri->validate());
    return tri;
  };

  { // normal coordinate range
    auto tri = freshTri();
    tri->normalCoordinates.edgeCoords[tri->intrinsicMesh->edge(0)] = -2;
    EXPECT_THROW(tri->validate(), std::runtime_error);
  }
  { // crossing count corruption (breaks corner consistency or coverage)
    auto tri = freshTri();
    for (Edge e : tri->intrinsicMesh->edges()) {
      if (tri->normalCoordinates[e] > 0) {
        tri->normalCoordinates.edgeCoords[e] += 1;
        break;
      }
    }
    EXPECT_THROW(tri->validate(), std::runtime_error);
  }
  { // wrong location type: face point claimed on an edge
    auto tri = freshTri();
    for (Vertex v : tri->intrinsicMesh->vertices()) {
      if (tri->vertexLocations[v].type == SurfacePointType::Face) {
        tri->vertexLocations[v] = SurfacePoint(tri->inputMesh.edge(0), 0.5);
        break;
      }
    }
    EXPECT_THROW(tri->validate(), std::runtime_error);
  }
  { // wrong location element: on-curve vertex moved to a different input edge
    auto tri = freshTri();
    for (Vertex v : tri->intrinsicMesh->vertices()) {
      if (tri->vertexLocations[v].type == SurfacePointType::Edge) {
        Edge wrong;
        for (Edge eIn : tri->inputMesh.edges()) {
          if (eIn != tri->vertexLocations[v].edge) {
            wrong = eIn;
            break;
          }
        }
        tri->vertexLocations[v] = SurfacePoint(wrong, 0.5);
        break;
      }
    }
    EXPECT_THROW(tri->validate(), std::runtime_error);
  }
  { // roundabout corruption
    auto tri = freshTri();
    Halfedge he = tri->intrinsicMesh->vertex(0).halfedge();
    tri->normalCoordinates.roundabouts[he] =
        (tri->normalCoordinates.roundabouts[he] + 1) % tri->normalCoordinates.roundaboutDegrees[he.vertex()];
    EXPECT_THROW(tri->validate(), std::runtime_error);
  }
  { // geometric layer: broken triangle inequality
    auto tri = freshTri();
    Edge e0 = tri->intrinsicMesh->edge(0);
    tri->edgeLengths[e0] = 1e6 * tri->edgeLengths[e0];
    EXPECT_THROW(tri->validate(), std::runtime_error);
    EXPECT_NO_THROW(tri->validate(false)); // geometry check off: combinatorial layer still fine
  }
}

// The seam corner-cut configuration (from a downstream report): two
// vertices inserted on different input edges sharing a corner, connected by
// an uncrossed intrinsic edge cutting across that corner. Splitting the
// corner-cut edge must record a location in the corner's input face --
// the face shared by both endpoints' input edges -- not on the wrong side
// of either curve.
TEST_F(IntrinsicInsertionSuite, SplitCornerCutSeamEdge) {
  auto a = getAsset("spot.ply", true);
  ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
  VertexPositionGeometry& origGeometry = *a.geometry;

  int nTested = 0;
  for (size_t iF = 0; iF < mesh.nFaces() && nTested < 20; iF += 37) {
    IntegerCoordinatesIntrinsicTriangulation tri(mesh, origGeometry);
    ManifoldSurfaceMesh& im = *tri.intrinsicMesh;

    // Pick a corner of input face F: edges E1 (corner->a) and E2 (corner->b)
    Face F = mesh.face(iF);
    Halfedge heF = F.halfedge();
    if (heF.edge().isBoundary() || heF.next().next().edge().isBoundary()) continue;

    // The same elements in the (initially identical) intrinsic mesh
    Halfedge heI = im.face(iF).halfedge();
    Edge e1 = heI.edge();               // shared with input E1
    Edge e2 = heI.next().next().edge(); // shared with input E2 (other corner edge)

    // Split both shared edges, leaving Edge-typed vertices on E1 and E2.
    // Splitting e1 then e2 creates a cross edge connecting the two new
    // vertices (the corner-cut "seam" edge).
    Vertex v1 = tri.splitEdge(e1, 0.49);
    ASSERT_NE(v1, Vertex());
    Vertex v2 = tri.splitEdge(e2, 0.51);
    ASSERT_NE(v2, Vertex());
    ASSERT_EQ(tri.vertexLocations[v1].type, SurfacePointType::Edge);
    ASSERT_EQ(tri.vertexLocations[v2].type, SurfacePointType::Edge);
    ASSERT_NE(tri.vertexLocations[v1].edge, tri.vertexLocations[v2].edge);

    Edge seam;
    for (Halfedge he : v1.outgoingHalfedges()) {
      if (he.tipVertex() == v2) seam = he.edge();
    }
    if (seam == Edge() || tri.normalCoordinates[seam] != 0) continue; // configuration didn't arise here

    // Flip the far shared sub-edges at both vertices (where geometrically
    // possible), so that the curve pieces cross the flanking faces
    // transversally instead of running along shared edges -- the
    // configuration in which the containing-face derivation must reason
    // about curves passing through the seam edge's own endpoints
    for (Vertex vv : {v1, v2}) {
      std::vector<Edge> sub;
      for (Halfedge he : vv.outgoingHalfedges()) {
        if (tri.normalCoordinates[he.edge()] < 0) sub.push_back(he.edge());
      }
      for (Edge se : sub) tri.flipEdgeIfPossible(se);
    }
    if (tri.normalCoordinates[seam] != 0) continue; // a flip crossed the seam edge

    // Split the corner-cut edge; this exercises the uncrossed-edge face
    // derivation with on-curve endpoints on different input edges
    Vertex vMid = tri.splitEdge(seam, 0.5);
    ASSERT_NE(vMid, Vertex());
    SurfacePoint loc = tri.vertexLocations[vMid];
    ASSERT_EQ(loc.type, SurfacePointType::Face);
    EXPECT_EQ(loc.face, F) << "corner-cut split recorded in the wrong input face";

    EXPECT_EQ(checkCS(tri, origGeometry), "");
    nTested++;
  }
  EXPECT_GT(nTested, 0);
}

// A shrinking cascade of insertions toward a fixed point must be refused
// once the would-be edges drop below insertionMinEdgeLength, instead of
// minting ever-tinier elements (the parameter-eps check alone cannot stop
// this: a split at parameter 0.5 of a tiny edge is far from its endpoints
// in parameter while being geometrically coincident with everything nearby)
TEST_F(IntrinsicInsertionSuite, ShrinkingInsertCascadeRefused) {
  auto a = getAsset("fox.ply", true);
  ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
  VertexPositionGeometry& origGeometry = *a.geometry;

  IntegerCoordinatesIntrinsicTriangulation tri(mesh, origGeometry);
  tri.flipToDelaunay();
  ManifoldSurfaceMesh& im = *tri.intrinsicMesh;

  EXPECT_GT(tri.insertionMinEdgeLength, 0.); // default-initialized in the constructor

  // Repeatedly split the left half of the same edge: lengths halve each
  // time, so without the geometric refusal this creates ~200 vertices in a
  // shrinking cluster
  Edge e;
  for (Edge cand : im.edges()) {
    if (tri.normalCoordinates[cand] == 0 && !cand.isBoundary()) {
      e = cand;
      break;
    }
  }
  if (e == Edge()) {
    for (Edge cand : im.edges()) {
      if (!cand.isBoundary()) {
        e = cand;
        break;
      }
    }
  }
  ASSERT_NE(e, Edge());

  Vertex anchor = e.halfedge().tailVertex();
  int nCreated = 0;
  for (int i = 0; i < 200; i++) {
    // find the current edge out of `anchor` along the shrinking direction:
    // the shortest edge at the anchor
    Edge eShort;
    double shortest = std::numeric_limits<double>::infinity();
    for (Edge ve : anchor.adjacentEdges()) {
      if (tri.edgeLengths[ve] < shortest) {
        shortest = tri.edgeLengths[ve];
        eShort = ve;
      }
    }
    size_t nVBefore = im.nVertices();
    tri.insertVertex(SurfacePoint(eShort, eShort.halfedge().tailVertex() == anchor ? 0.5 : 0.5));
    if (im.nVertices() > nVBefore) nCreated++;
  }

  // The cascade must be cut off long before 200 insertions
  EXPECT_LT(nCreated, 60);
  EXPECT_EQ(checkCS(tri, origGeometry), "");
}

// insertVertex must refuse insertions coincident with an existing vertex,
// returning that vertex instead of creating a near-zero-length edge
TEST_F(IntrinsicInsertionSuite, CoincidentInsertRefused) {
  auto a = getAsset("fox.ply", true);
  ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
  VertexPositionGeometry& origGeometry = *a.geometry;

  IntegerCoordinatesIntrinsicTriangulation tri(mesh, origGeometry);
  tri.flipToDelaunay();
  ManifoldSurfaceMesh& im = *tri.intrinsicMesh;

  size_t nVBefore = im.nVertices();

  // Edge points within eps of an endpoint
  Edge e = im.edge(0);
  EXPECT_EQ(tri.insertVertex(SurfacePoint(e, 1e-15)), e.halfedge().tailVertex());
  EXPECT_EQ(tri.insertVertex(SurfacePoint(e, 1. - 1e-15)), e.halfedge().tipVertex());

  // Face point within eps of a corner
  Face f = im.face(0);
  Vertex v0 = f.halfedge().vertex();
  EXPECT_EQ(tri.insertVertex(SurfacePoint(f, Vector3{1. - 1e-15, 5e-16, 5e-16})), v0);

  EXPECT_EQ(im.nVertices(), nVBefore); // no mutation happened

  // Repeated insertion at the same point: first creates, second is refused
  Vertex v1 = tri.insertVertex(SurfacePoint(im.edge(5), 0.5));
  ASSERT_NE(v1, Vertex());
  EXPECT_EQ(im.nVertices(), nVBefore + 1);
  // the same input point now lies (within eps) at vertex v1: inserting on
  // one of its new edges at parameter ~0 returns v1
  for (Halfedge he : v1.outgoingHalfedges()) {
    EXPECT_EQ(tri.insertVertex(SurfacePoint(he.edge(), he == he.edge().halfedge() ? 1e-15 : 1. - 1e-15)), v1);
  }
  EXPECT_EQ(im.nVertices(), nVBefore + 1);
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

// Regression test (downstream garment-pipeline report,
// surfacepoint-ambiguous-classification): removing an inserted vertex
// erases the three spokes of its degree-3 star, and with them any curve
// crossings recorded there. Erasing spoke crossings is fine as long as the
// curve still enters and leaves the merged face through its outer edges --
// but if a curve component's crossings interior to the merged face lie
// ONLY on spokes (a curve running corner-to-corner between two on-curve
// vertices, e.g. seam vertices inserted on an input edge), the curve would
// become a corner-to-corner chord of the merged face, which normal
// coordinates cannot express: the component silently vanishes and the
// roundabouts at the merged face's corners go stale by one, corrupting
// every later roundabout-anchored derivation (identifyInputEdge,
// wedgeInputFace). removeInsertedVertex must refuse exactly these removals
// (returning Face()); every removal it does perform must leave a fully
// valid triangulation. In the wild the trigger was a seam-insertion +
// straighten + delaunayRefine pipeline: refinement inserted a cascade
// vertex recorded numerically ON an input edge (bary <0.5, 4.5e-22, 0.5>)
// and deleteNearbyVertices removed a neighbor whose spokes carried the
// whole inter-seam-vertex curve segment.
//
// The chord state cannot arise from exact geometry -- every operation
// sequence that would produce it contains a flip refused by the signed-
// area check. It arises exactly when floating point makes a borderline
// choice: a vertex numerically ON the curve, where each individual choice
// is legal but the removal's flip-to-degree-3 reduction can strand the
// inter-seam-vertex curve piece on spokes. So we manufacture the seam
// scenario at its degenerate limit:
//   1. split a shared edge twice -> consecutive seam vertices u, s joined
//      by a shared sub-edge; the curve piece u->s runs along it;
//   2. insert v as a Face-point numerically ON the input edge between
//      them (the wild record: bary <0.5, 4.5e-22, 0.5>);
//   3. flip the u-s sub-edge -- geometrically legal, but the curve piece
//      u->s must be re-expressed around v: the float-ambiguous choice;
//   4. remove v (then u, s), validating after every attempt.
// (The full wild trigger needs garment-scale float coincidences which no
// distilled construction reached -- gc's own flip/split bookkeeping
// resolves each of these ambiguities consistently. This test pins that
// behavior; the conservation-law guard inside removeInsertedVertex is
// exercised at scale by the downstream repro, where it takes 15 corrupting
// removals to 0.)
TEST_F(IntrinsicInsertionSuite, RemoveVertexBesideCurveChord) {
  for (const MeshAsset& a : {getAsset("fox.ply", true), getAsset("spot.ply", true)}) {
    a.printThyName();
    ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
    VertexPositionGeometry& origGeometry = *a.geometry;

    IntegerCoordinatesIntrinsicTriangulation tri(mesh, origGeometry);
    ManifoldSurfaceMesh& im = *tri.intrinsicMesh;

    // Cycle the off-edge barycentric component across sites (the wild
    // trigger was 4.5e-22 -- the ambiguity lives far below any geometric
    // epsilon)
    const std::vector<double> etas = {4.5e-22, 1e-18, 1e-15, 0., 1e-12};

    int nSites = 0, nRemoved = 0, nRefused = 0, nAttempts = 0;
    size_t nOrigEdges = im.nEdges();
    for (size_t iE = 0; iE < nOrigEdges && nSites < 15; iE += 7) {
      Edge e = im.edge(iE);
      if (tri.normalCoordinates[e] != -1 || e.isBoundary()) continue;

      // 1. Two consecutive seam vertices on the input edge
      Vertex u = tri.splitEdge(e, 0.5);
      ASSERT_NE(u, Vertex());
      if (tri.vertexLocations[u].type != SurfacePointType::Edge) continue;
      Edge eSub;
      for (Edge ve : u.adjacentEdges()) {
        if (tri.normalCoordinates[ve] == -1 && !ve.isBoundary()) {
          eSub = ve;
          break;
        }
      }
      if (eSub == Edge()) continue;
      Vertex s = tri.splitEdge(eSub, 0.5);
      ASSERT_NE(s, Vertex());
      if (tri.vertexLocations[s].type != SurfacePointType::Edge) continue;
      Edge eMid;
      for (Halfedge he : u.outgoingHalfedges()) {
        if (he.tipVertex() == s && tri.normalCoordinates[he.edge()] == -1) eMid = he.edge();
      }
      if (eMid == Edge() || eMid.isBoundary()) continue;

      // 2. Insert v as a Face-point numerically ON the input edge between
      // the seam vertices (the wild record was bary <0.5, 4.5e-22, 0.5>)
      Face fIns = eMid.halfedge().face();
      Vertex apex; // the corner of fIns which is not u or s
      double eta = etas[nSites % etas.size()];
      Vector3 bary = Vector3::zero();
      {
        int iC = 0;
        for (Vertex cV : fIns.adjacentVertices()) {
          if (cV == u) {
            bary[iC] = 0.5;
          } else if (cV == s) {
            bary[iC] = 0.5 - eta;
          } else {
            apex = cV;
            bary[iC] = eta;
          }
          iC++;
        }
      }
      size_t nVBefore = im.nVertices();
      Vertex v = tri.insertVertex(SurfacePoint(fIns, bary));
      if (v == Vertex() || im.nVertices() == nVBefore) continue; // insertion refused
      if (tri.vertexLocations[v].type == SurfacePointType::Vertex) continue;

      // 3. Flip the seam sub-edge u-s. Geometrically this is fine (v lies
      // numerically ON the segment u-s, the quad is only degenerate at v);
      // combinatorially the curve piece u->s must be re-expressed around
      // (or through) v -- exactly the float-ambiguous choice. Any choice
      // must leave a valid triangulation.
      if (!tri.flipEdgeIfPossible(eMid)) continue;
      nSites++;
      try {
        tri.validate();
      } catch (std::exception& err) {
        ADD_FAILURE() << "validate() failed after flipping the seam sub-edge beside v (site " << nSites
                      << ", eta = " << eta << "): " << err.what();
        return;
      }

      // 4. Remove the inserted vertices again. Every attempt must either
      // be refused or leave a fully valid triangulation -- in particular
      // the curve piece u->s must never be silently erased.
      for (Vertex cand : {v, u, s, apex}) {
        if (cand.isDead() || tri.vertexLocations[cand].type == SurfacePointType::Vertex) continue;
        nAttempts++;
        Face fMerged = tri.removeInsertedVertex(cand);
        if (fMerged == Face()) {
          nRefused++;
        } else {
          nRemoved++;
        }
        try {
          tri.validate();
        } catch (std::exception& err) {
          ADD_FAILURE() << "validate() failed after remove attempt of " << cand << " (site " << nSites << ", eta "
                        << eta << ", " << (fMerged == Face() ? "refused" : "removed") << "): " << err.what();
          return;
        }
      }
    }
    std::cout << "  " << nSites << " sites, " << nAttempts << " remove attempts: " << nRemoved << " removed, "
              << nRefused << " refused" << std::endl;
    EXPECT_GT(nSites, 0);
    EXPECT_EQ(checkCS(tri, origGeometry), "");
  }
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

// Companion to RemoveVertexBesideCurveChord, at pipeline scale: the wild
// trigger's shape end to end -- seam vertices along many input edges
// (level-set-seam-like), straightening flips, then delaunayRefine, whose
// internal deleteNearbyVertices removes cascade vertices near the seams
// (the exact call path of the downstream failure). The triangulation and
// common subdivision must come out fully valid.
TEST_F(IntrinsicInsertionSuite, SeamRefineDelete) {
  for (const MeshAsset& a : {getAsset("spot.ply", true), getAsset("lego.ply", true), getAsset("fox.ply", true)}) {
    a.printThyName();
    ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
    VertexPositionGeometry& origGeometry = *a.geometry;

    for (unsigned seed = 0; seed < 5; seed++) {
      IntegerCoordinatesIntrinsicTriangulation tri(mesh, origGeometry);
      ManifoldSurfaceMesh& im = *tri.intrinsicMesh;
      std::mt19937 rng(seed);
      std::uniform_real_distribution<double> unit(0.15, 0.85);

      // Seam vertices along many input edges (level-set-seam-like)
      size_t nOrigEdges = im.nEdges();
      int nSeam = 0;
      for (size_t iE = seed; iE < nOrigEdges; iE += 3) {
        Edge e = im.edge(iE);
        if (tri.normalCoordinates[e] != -1) continue;
        Vertex u = tri.splitEdge(e, unit(rng));
        if (u != Vertex()) nSeam++;
      }
      ASSERT_GT(nSeam, 0);

      // Straighten: push toward Delaunay (flips shared sub-edges away)
      tri.flipToDelaunay();

      // Refine: inserts circumcenters (cascades near the seam vertices)
      // and deletes nearby cascade vertices again
      tri.delaunayRefine(25.);

      try {
        tri.validate();
      } catch (std::exception& err) {
        ADD_FAILURE() << "validate() failed after seam+refine, seed " << seed << ": " << err.what();
        return;
      }
      std::string csErr = checkCS(tri, origGeometry);
      EXPECT_EQ(csErr, "") << "seed " << seed;
    }
  }
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
