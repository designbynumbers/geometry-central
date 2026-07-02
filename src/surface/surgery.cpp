#include "geometrycentral/surface/surgery.h"

#include <queue>

namespace geometrycentral {
namespace surface {

std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, HalfedgeData<Halfedge>> cutAlongEdges(ManifoldSurfaceMesh& origMesh,
                                                                                const EdgeData<char>& origCut) {

  // Create a copy of the input mesh
  std::unique_ptr<ManifoldSurfaceMesh> mesh = origMesh.copy();

  // Initialize parent references
  // The built-in dynamic container updates will automatically keep this container in sync as we modify the mesh
  HalfedgeData<Halfedge> parentHalfedges(*mesh, Halfedge());
  for (size_t i = 0; i < mesh->nHalfedges(); i++) {
    parentHalfedges[i] = origMesh.halfedge(i);
  }

  // Transfer the cut data to the new mesh.
  EdgeData<char> cut = origCut.reinterpretTo(*mesh);
  cut.setDefault(false);

  // Cut along all cut edges
  // TODO use this instead of tree trick below
  // TODO modifying while iterating
  /*
  for(Edge e : mesh->edges()) {
    if(cut[e]) {
      mesh->separateEdge(e); // all the hard works happens here
    }
  }
  */


  // TODO Right now separateEdge() can only handle cutting along a tree, so walk along tree
  EdgeData<char> considered(*mesh, false);

  // A valid cut may be several disjoint trees (e.g. several cone-to-boundary
  // slits that share no vertex), so repeat until every cut edge has been
  // processed. Within each tree, prefer to start from a leaf edge that
  // already touches the mesh boundary and walk inward: every subsequent
  // separateEdge() call then has exactly one already-boundary endpoint (the
  // just-extended loop) and one still-interior endpoint, which is always a
  // supported case (separateEdge's Case 2). Starting from the opposite end
  // (e.g. an interior cone) instead opens a fresh, separate mini boundary
  // loop first, which must eventually MERGE with the original boundary loop
  // to finish the tree -- unimplemented (separateEdge's Case 3). A tree with
  // no boundary touch point at all (a closed-surface cone-to-cone bridge) has
  // no such preference to make; any leaf works.
  while (true) {
    Edge firstEdge;
    bool anyUnconsidered = false;
    for (Edge e : mesh->edges()) {
      if (!cut[e] || considered[e]) continue;
      anyUnconsidered = true;

      int topCount = 0;
      for (Edge en : e.halfedge().vertex().adjacentEdges()) {
        if (cut[en]) topCount++;
      }
      int botCount = 0;
      for (Edge en : e.halfedge().twin().vertex().adjacentEdges()) {
        if (cut[en]) botCount++;
      }
      bool boundaryLeaf =
          (topCount == 1 && e.halfedge().vertex().isBoundary()) ||
          (botCount == 1 && e.halfedge().twin().vertex().isBoundary());
      if (boundaryLeaf) {
        firstEdge = e;
        break; // prefer a boundary-touching leaf outright
      }
      if (firstEdge == Edge() && (topCount == 1 || botCount == 1)) {
        firstEdge = e; // remember the first non-boundary leaf as a fallback
      }
    }
    if (!anyUnconsidered) break; // every cut edge has been processed
    if (firstEdge == Edge())
      throw std::runtime_error("could not find leaf edge. must cut along simple disk tree. see note");

    std::queue<Edge> queue;
    queue.emplace(firstEdge);
    considered[firstEdge] = true;

    // process until queue is empty
    while (!queue.empty()) {
      Edge e = queue.front();
      queue.pop();

      // Cache to restore after separate
      Halfedge oldParent = parentHalfedges[e.halfedge()];
      Halfedge oldParentT = parentHalfedges[e.halfedge().twin()];

      Halfedge newHe, newHeOpp;
      std::tie(newHe, newHeOpp) = mesh->separateEdge(e); // all the hard works happens here

      // Keep the parent map accurate
      parentHalfedges[newHe] = oldParent;
      parentHalfedges[newHeOpp] = oldParentT;

      // Add new neighbors for processing
      for (Edge en : e.halfedge().vertex().adjacentEdges()) {
        if (cut[en] && !considered[en]) {
          considered[en] = true;
          queue.emplace(en);
        }
      }
      for (Edge en : e.halfedge().twin().vertex().adjacentEdges()) {
        if (cut[en] && !considered[en]) {
          considered[en] = true;
          queue.emplace(en);
        }
      }
    }
  }

  // Clear out any stale boundary parents
  for (Halfedge he : mesh->exteriorHalfedges()) {
    parentHalfedges[he] = Halfedge();
  }

  return std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, HalfedgeData<Halfedge>>{std::move(mesh), parentHalfedges};
}


} // namespace surface
} // namespace geometrycentral

