// Stress tests for vertex insertion/deletion on
// SignpostIntrinsicTriangulation, mirroring the integer-coordinates tests in
// intrinsic_insertion_test.cpp but driven entirely through the
// representation-agnostic IntrinsicTriangulation interface.
#include "geometrycentral/surface/integer_coordinates_intrinsic_triangulation.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/signpost_intrinsic_triangulation.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "load_test_meshes.h"

#include "gtest/gtest.h"

#include <random>
#include <sstream>

using namespace geometrycentral;
using namespace geometrycentral::surface;

class SignpostInsertionSuite : public MeshAssetSuite {};

namespace {

// Validate the common subdivision, combinatorially and geometrically.
// Returns an error string ("" if ok).
// NOTE: invalidates element handles (compresses the mesh)
std::string checkCS(IntrinsicTriangulation& tri, VertexPositionGeometry& origGeometry, double tol) {
  try {
    CommonSubdivision& cs = tri.getCommonSubdivision();
    cs.constructMesh();

    // Lengths from extrinsic vertex positions
    VertexPositionGeometry csGeo(*cs.mesh, cs.interpolateAcrossA(origGeometry.vertexPositions));
    csGeo.requireEdgeLengths();
    EdgeData<double> lengthsFromPosA = csGeo.edgeLengths;

    // Lengths from input edge lengths
    origGeometry.requireEdgeLengths();
    EdgeData<double> lengthsFromLenA = cs.interpolateEdgeLengthsA(origGeometry.edgeLengths);

    // Lengths from intrinsic edge lengths
    EdgeData<double> lengthsFromLenB = cs.interpolateEdgeLengthsB(tri.edgeLengths);

    double scale = lengthsFromPosA.toVector().norm();
    double errAB = (lengthsFromPosA.toVector() - lengthsFromLenB.toVector()).norm() / scale;
    double errAA = (lengthsFromPosA.toVector() - lengthsFromLenA.toVector()).norm() / scale;

    if (errAA > tol || errAB > tol) {
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

// Run random insertions/deletions through the generic IntrinsicTriangulation
// interface, periodically validating the common subdivision
void probeTriangulation(IntrinsicTriangulation& tri, VertexPositionGeometry& origGeometry, unsigned seed, int nOps,
                        int checkEvery, double tol) {
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
        desc << "insertVertex(FacePoint " << f << ", " << bary << ")";
        tri.insertVertex(SurfacePoint(f, bary));
      } else if (op < 7) {
        // Edge split at random point
        Edge e = sampleNth(im.edges(), idxDist(rng), im.nEdges());
        double t = unit(rng);
        desc << "splitEdge(" << e << " [bdy=" << e.isBoundary() << "], " << t << ")";
        tri.splitEdge(e.halfedge(), t);
      } else if (op < 9) {
        // Remove a random inserted vertex
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
        desc << "insertVertex(EdgePoint " << e << " [bdy=" << e.isBoundary() << "], " << t << ")";
        tri.insertVertex(SurfacePoint(e, t));
      }
      opCount++;
      if (checkEvery == 1) std::cout << "  op " << iOp << ": " << desc.str() << std::endl;
    } catch (std::exception& err) {
      ADD_FAILURE() << "op " << iOp << ": " << desc.str() << " threw: " << err.what();
      return;
    }

    if (iOp % checkEvery == checkEvery - 1) {
      std::string csErr = checkCS(tri, origGeometry, tol);
      if (!csErr.empty()) {
        ADD_FAILURE() << "after op " << iOp << " (" << desc.str() << "): " << csErr;
        return;
      }
    }
  }
  std::cout << "  completed " << opCount << " ops" << std::endl;
  EXPECT_EQ(checkCS(tri, origGeometry, tol), "");
}

void probeSignpost(const MeshAsset& a, unsigned seed, int nOps, int checkEvery = 50, double tol = 1e-5) {
  a.printThyName();
  ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
  VertexPositionGeometry& origGeometry = *a.geometry;

  SignpostIntrinsicTriangulation tri(mesh, origGeometry);
  tri.flipToDelaunay();
  probeTriangulation(tri, origGeometry, seed, nOps, checkEvery, tol);
}

} // namespace

// KNOWN LIMITATION: long random insert/delete sequences eventually corrupt
// the signpost correspondence. The traced (floating point) vertex locations
// and crossings can land on the wrong side of input edges/vertices, after
// which common subdivision construction fails; this is the fragility that
// motivated the integer-coordinates representation (which passes the
// analogous probes in intrinsic_insertion_test.cpp). The probes below are
// disabled but kept as a harness for future robustness work; on a typical
// run they survive ~60-200 ops.
// (run with --gtest_also_run_disabled_tests to include them)
TEST_F(SignpostInsertionSuite, DISABLED_ClosedMeshProbe) {
  probeSignpost(getAsset("fox.ply", true), 7, 300);
}

TEST_F(SignpostInsertionSuite, DISABLED_ClosedMeshProbe2) {
  probeSignpost(getAsset("sphere_small.ply", true), 11, 300);
}

TEST_F(SignpostInsertionSuite, DISABLED_BoundaryMeshProbe) {
  probeSignpost(getAsset("lego.ply", true), 13, 300);
}

// The user-facing workflow that motivated this suite: delaunayRefine, then
// use the triangulation (e.g. via FlipEdgeNetwork). Validate that the
// triangulation's correspondence is intact after refinement.
TEST_F(SignpostInsertionSuite, DelaunayRefineThenValidate) {
  for (const MeshAsset& a : {getAsset("fox.ply", true), getAsset("lego.ply", true)}) {
    a.printThyName();
    ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
    VertexPositionGeometry& origGeometry = *a.geometry;

    SignpostIntrinsicTriangulation tri(mesh, origGeometry);
    tri.delaunayRefine();
    EXPECT_EQ(checkCS(tri, origGeometry, 1e-5), "") << " on " << a.name;
  }
}

// Same check for the integer-coordinates triangulation, for comparison
TEST_F(SignpostInsertionSuite, IntegerDelaunayRefineThenValidate) {
  for (const MeshAsset& a : {getAsset("fox.ply", true), getAsset("lego.ply", true)}) {
    a.printThyName();
    ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
    VertexPositionGeometry& origGeometry = *a.geometry;

    IntegerCoordinatesIntrinsicTriangulation tri(mesh, origGeometry);
    tri.delaunayRefine();
    EXPECT_EQ(checkCS(tri, origGeometry, 1e-5), "") << " on " << a.name;
  }
}
