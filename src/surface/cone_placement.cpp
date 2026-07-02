#include "geometrycentral/surface/cone_placement.h"

#include "geometrycentral/numerical/linear_solvers.h"

#include <cmath>
#include <limits>

namespace geometrycentral {
namespace surface {

namespace {

// Solve for the interior conformal scale factors anchored at the cones/boundary:
//   A_nn u_n = -K_n,   u = 0 on the anchor set (cones + boundary),
// where `n` is the set of free (non-anchor) vertices marked true in `isN`.
// Returns a full-length vector (0 outside `n`). If there is no anchor (every
// vertex free), the system is singular, so we return all zeros -- the caller's
// greedy step then deterministically seeds the first cone.
Vector<double> solveScaleFactors(const SparseMatrix<double>& A, const Vector<double>& K,
                                 const Vector<bool>& isN) {
  size_t nV = static_cast<size_t>(A.rows());

  size_t nFree = 0;
  for (size_t i = 0; i < nV; i++)
    if (isN(i)) nFree++;

  Vector<double> u = Vector<double>::Zero(nV);
  if (nFree == 0 || nFree == nV) return u; // nothing to solve, or no anchor

  BlockDecompositionResult<double> decomp = blockDecomposeSquare(A, isN);
  SparseMatrix<double> Ann = decomp.AA;
  shiftDiagonal(Ann, 1e-12); // guard against round-off indefiniteness

  Vector<double> Kn, Kanchor;
  decomposeVector(decomp, K, Kn, Kanchor);

  PositiveDefiniteSolver<double> solver(Ann);
  Vector<double> un = solver.solve(Vector<double>(-Kn));

  Vector<double> uAnchor = Vector<double>::Zero(Kanchor.rows());
  return reassembleVector(decomp, un, uAnchor);
}

} // namespace

ConePlacementResult computeConePlacement(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geo,
                                         const ConePlacementOptions& opts) {

  geo.requireCotanLaplacian();
  geo.requireVertexAngleSums();

  SparseMatrix<double> A = geo.cotanLaplacian;
  VertexData<size_t> vIdx = mesh.getVertexIndices();
  size_t nV = mesh.nVertices();

  // Discrete Gaussian curvature / boundary exterior angle:
  //   interior: 2*pi - angleSum,   boundary: pi - angleSum.
  Vector<double> K(nV);
  for (Vertex v : mesh.vertices()) {
    double base = v.isBoundary() ? M_PI : 2.0 * M_PI;
    K(vIdx[v]) = base - geo.vertexAngleSums[v];
  }

  // Initialize the anchor (cone) set.
  VertexData<char> isCone(mesh, 0);
  bool hasBoundary = mesh.nBoundaryLoops() > 0;
  if (hasBoundary) {
    // All boundary vertices anchor the flattening (they are not counted as cones).
    for (Vertex v : mesh.vertices())
      if (v.isBoundary()) isCone[v] = 1;
  } else {
    // Closed surface: seed with the single extreme-curvature vertex (largest for
    // chi > 0, smallest for chi < 0). For chi == 0 (torus) no seed is placed.
    int chi = mesh.eulerCharacteristic();
    if (chi != 0) {
      Vertex best;
      double bestK = (chi > 0) ? -std::numeric_limits<double>::infinity()
                               : std::numeric_limits<double>::infinity();
      for (Vertex v : mesh.vertices()) {
        double kv = K(vIdx[v]);
        if ((chi > 0 && kv > bestK) || (chi < 0 && kv < bestK)) {
          bestK = kv;
          best = v;
        }
      }
      isCone[best] = 1;
    }
  }

  // Caller-supplied cones join the anchor set immediately (before any greedy
  // step), on top of the boundary/auto-seed anchors above.
  for (Vertex v : opts.initialCones)
    if (!v.isBoundary()) isCone[v] = 1;

  auto buildIsN = [&]() {
    Vector<bool> isN(nV);
    for (Vertex v : mesh.vertices())
      isN(vIdx[v]) = (!v.isBoundary() && !isCone[v]);
    return isN;
  };

  // Greedy: add the vertex of largest |scale factor|, up to opts.maxNewCones
  // times (or until opts.stopMaxU is satisfied).
  Vector<double> u = solveScaleFactors(A, K, buildIsN());
  for (size_t it = 0; it < opts.maxNewCones; it++) {
    Vertex best;
    double maxAbs = -1.0;
    for (Vertex v : mesh.vertices()) {
      if (v.isBoundary() || isCone[v]) continue;
      double a = std::abs(u(vIdx[v]));
      if (a > maxAbs) { // strict: ties resolve to the lowest vertex index
        maxAbs = a;
        best = v;
      }
    }
    if (maxAbs < 0.0) break; // no remaining candidate vertex
    if (opts.stopMaxU > 0.0 && maxAbs < opts.stopMaxU) break; // gate already satisfied
    isCone[best] = 1;
    u = solveScaleFactors(A, K, buildIsN());
  }

  // Prescribe cone angles (CETM): C = K + A*u at cone/boundary vertices.
  Vector<double> Au = A * u;
  VertexData<double> coneAngles(mesh, 0.0);
  double total = 0.0;
  for (Vertex v : mesh.vertices()) {
    size_t i = vIdx[v];
    if (v.isBoundary() || isCone[v]) {
      coneAngles[v] = K(i) + Au(i);
      total += coneAngles[v];
    }
  }

  // Normalize so the total prescribed curvature satisfies discrete Gauss-Bonnet.
  double chi = static_cast<double>(mesh.eulerCharacteristic());
  if (std::abs(total) > 1e-8) {
    double f = (2.0 * M_PI * chi) / total;
    for (Vertex v : mesh.vertices())
      coneAngles[v] *= f;
  }

  std::vector<Vertex> cones;
  for (Vertex v : mesh.vertices())
    if (!v.isBoundary() && isCone[v]) cones.push_back(v);

  VertexData<double> uField(mesh, 0.0);
  for (Vertex v : mesh.vertices()) uField[v] = u(vIdx[v]);

  geo.unrequireCotanLaplacian();
  geo.unrequireVertexAngleSums();

  return {cones, coneAngles, uField};
}

ConePlacementResult computeConePlacement(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geo,
                                         size_t nCones) {
  // Reproduce the pre-ConePlacementOptions semantics exactly: nCones counted the
  // closed-surface auto-seed (if any) toward its total, so the greedy loop only
  // added nCones - autoSeedCount new cones.
  bool hasBoundary = mesh.nBoundaryLoops() > 0;
  size_t autoSeedCount = (!hasBoundary && mesh.eulerCharacteristic() != 0) ? 1 : 0;

  ConePlacementOptions opts;
  opts.maxNewCones = (nCones > autoSeedCount) ? (nCones - autoSeedCount) : 0;
  return computeConePlacement(mesh, geo, opts);
}

} // namespace surface
} // namespace geometrycentral
