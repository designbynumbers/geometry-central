#include "geometrycentral/surface/cone_cut.h"

#include "geometrycentral/utilities/disjoint_sets.h"

#include <queue>
#include <tuple>

namespace geometrycentral {
namespace surface {

namespace {

// Approximate minimum-weight Steiner tree (by edge length) connecting all cones
// and boundary loops, returned as a set of original-mesh edges. This mirrors the
// BFF reference Cutter using gc's MarkedDisjointSets: grow a Prim/Dijkstra
// forest seeded at every marked component (boundary loops + cones), and whenever
// an edge bridges two distinct marked components, carve the cut along the two
// parent-pointer chains plus the bridge edge.
EdgeData<bool> steinerTreeCut(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geo,
                              const std::vector<Vertex>& cones) {
  geo.requireEdgeLengths();

  size_t nV = mesh.nVertices();
  EdgeData<bool> inPath(mesh, false);
  MarkedDisjointSets ds(nV);
  VertexData<Halfedge> parentHe(mesh, Halfedge()); // parentHe[v]: tail v, tip its predecessor
  VertexData<size_t> vIdx = mesh.getVertexIndices();

  // Entry: (cumulative weight, candidate vertex index, halfedge index). The
  // halfedge goes v2 (in tree) -> v1 (candidate): he.vertex()==v2, tipVertex==v1.
  using Entry = std::tuple<double, size_t, size_t>;
  std::priority_queue<Entry, std::vector<Entry>, std::greater<Entry>> pq;

  auto pushNeighbors = [&](Vertex v, double baseWeight) {
    for (Halfedge he : v.outgoingHalfedges()) {
      Edge e = he.edge();
      if (e.isBoundary()) continue; // boundary edges are never cut
      pq.push(Entry(baseWeight + geo.edgeLengths[e], vIdx[he.tipVertex()], he.getIndex()));
    }
  };

  // Seed: merge + mark each boundary loop.
  for (BoundaryLoop b : mesh.boundaryLoops()) {
    Vertex prev;
    bool first = true;
    for (Halfedge he : b.adjacentHalfedges()) {
      Vertex v = he.vertex();
      if (first)
        first = false;
      else
        ds.merge(vIdx[prev], vIdx[v]);
      prev = v;
      pushNeighbors(v, 0.0);
    }
    if (!first) ds.mark(vIdx[prev]);
  }

  // Seed: mark each cone.
  for (Vertex c : cones) {
    pushNeighbors(c, 0.0);
    ds.mark(vIdx[c]);
  }

  // Grow the forest.
  while (!pq.empty()) {
    Entry top = pq.top();
    pq.pop();
    double weight = std::get<0>(top);
    Halfedge he = mesh.halfedge(std::get<2>(top));
    Vertex v2 = he.vertex();
    Vertex v1 = he.tipVertex();
    size_t i1 = vIdx[v1], i2 = vIdx[v2];

    if (ds.find(i1) == ds.find(i2)) continue;

    // Bridging two marked components: carve v1's chain, the bridge, v2's chain.
    if (ds.isMarked(i1) && ds.isMarked(i2)) {
      for (Halfedge h = parentHe[v1]; h != Halfedge(); h = parentHe[h.tipVertex()]) inPath[h.edge()] = true;
      inPath[he.edge()] = true;
      for (Halfedge h = parentHe[v2]; h != Halfedge(); h = parentHe[h.tipVertex()]) inPath[h.edge()] = true;
    }

    parentHe[v1] = he.twin(); // from v1 to v2
    ds.merge(i1, i2);
    pushNeighbors(v1, weight);
  }

  geo.unrequireEdgeLengths();
  return inPath;
}

} // namespace

ConeCutResult computeConeCut(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geo,
                             const std::vector<Vertex>& cones) {

  // 1. Cut topology: combinatorial Steiner tree over cones + boundary.
  EdgeData<bool> seed = steinerTreeCut(mesh, geo, cones);

  // 2. Anchor the cut endpoints (cones + boundary) so straightening can't slide
  //    a slit off a cone or off the boundary.
  VertexData<bool> markedVerts(mesh, false);
  for (Vertex v : mesh.vertices())
    if (v.isBoundary()) markedVerts[v] = true;
  for (Vertex c : cones)
    markedVerts[c] = true;

  // 3. Straighten to geodesics on an internally-built intrinsic triangulation.
  std::unique_ptr<FlipEdgeNetwork> net =
      FlipEdgeNetwork::constructFromEdgeSet(mesh, geo, seed, markedVerts);
  if (net && !net->paths.empty()) {
    net->iterativeShorten();
  }

  // 4. Collect the straightened cut edges (on the intrinsic mesh).
  ConeCutResult result;
  if (net) {
    EdgeData<bool> cutEdges(net->mesh, false);
    for (auto& path : net->paths) {
      for (Halfedge he : path->getHalfedgeList()) {
        cutEdges[he.edge()] = true;
      }
    }
    result.cutEdges = std::move(cutEdges);
  } else {
    result.cutEdges = EdgeData<bool>(mesh, false);
  }
  result.network = std::move(net);
  return result;
}

} // namespace surface
} // namespace geometrycentral
