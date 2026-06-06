#include "geometrycentral/surface/cone_parameterization.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "gtest/gtest.h"

#include <cmath>
#include <map>
#include <memory>
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
