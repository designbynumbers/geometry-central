#include "geometrycentral/surface/cone_parameterization.h"

#include "geometrycentral/numerical/linear_algebra_utilities.h"
#include "geometrycentral/numerical/linear_solvers.h"
#include "geometrycentral/surface/cone_cut.h"
#include "geometrycentral/surface/cone_placement.h"
#include "geometrycentral/surface/surgery.h"

#include <cmath>
#include <unordered_map>
#include <vector>

namespace geometrycentral {
namespace surface {

namespace {

// Solve A_nn x_n = rhs_n on the interior (free) set, returning a full-length
// vector with the prescribed values `bdy` on the anchor (boundary) set.
Vector<double> dirichletSolve(const SparseMatrix<double>& A, const Vector<bool>& isInterior,
                              const Vector<double>& rhsFull, const Vector<double>& bdyFull) {
  BlockDecompositionResult<double> decomp = blockDecomposeSquare(A, isInterior);
  SparseMatrix<double> Aii = decomp.AA;
  shiftDiagonal(Aii, 1e-12);

  Vector<double> rhsI, rhsB, bdyI, bdyB;
  decomposeVector(decomp, rhsFull, rhsI, rhsB);
  decomposeVector(decomp, bdyFull, bdyI, bdyB);

  // A_ii x_i = rhs_i - A_ib * bdy_b
  Vector<double> b = rhsI - decomp.AB * bdyB;
  PositiveDefiniteSolver<double> solver(Aii);
  Vector<double> xI = solver.solve(b);

  return reassembleVector(decomp, xI, bdyB);
}

} // namespace

ConeParameterizationResult flattenWithGivenConesAndCut(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geo,
                                                        const std::vector<Vertex>& cones,
                                                        const VertexData<double>& coneAngles,
                                                        const EdgeData<char>& cutEdges) {

  geo.requireCotanLaplacian();
  geo.requireVertexAngleSums();
  geo.requireEdgeLengths();

  size_t niV = mesh.nVertices();
  VertexData<size_t> iVIdx = mesh.getVertexIndices();

  // === Scale-factor field on the UNCUT surface, with the cone curvature
  //     subtracted from the source so curvature concentrates at the cones.
  //     Boundary scale factors are held at zero. ===
  Vector<double> Kminus(niV); // -(K - C)
  for (Vertex v : mesh.vertices()) {
    double base = v.isBoundary() ? M_PI : 2.0 * M_PI;
    double K = base - geo.vertexAngleSums[v];
    Kminus(iVIdx[v]) = -(K - coneAngles[v]);
  }
  Vector<bool> isInteriorI(niV);
  for (Vertex v : mesh.vertices()) isInteriorI(iVIdx[v]) = !v.isBoundary();
  Vector<double> zeroBdy = Vector<double>::Zero(niV);
  Vector<double> aUncut = dirichletSolve(geo.cotanLaplacian, isInteriorI, Kminus, zeroBdy);

  // === Cut the mesh into a disk. ===
  size_t nCutEdges = 0;
  for (Edge e : mesh.edges())
    if (cutEdges[e]) nCutEdges++;

  std::unique_ptr<ManifoldSurfaceMesh> cutMesh;
  HalfedgeData<Halfedge> parentHe;
  if (nCutEdges > 0) {
    std::tie(cutMesh, parentHe) = cutAlongEdges(mesh, cutEdges);
  } else {
    cutMesh = mesh.copy();
    parentHe = HalfedgeData<Halfedge>(*cutMesh);
    for (size_t i = 0; i < cutMesh->nHalfedges(); i++) parentHe[i] = mesh.halfedge(i);
  }

  if (cutMesh->nBoundaryLoops() != 1) {
    throw std::runtime_error("flattenWithGivenConesAndCut: cutting `mesh` along `cutEdges` must yield a single "
                             "topological disk (got " +
                             std::to_string(cutMesh->nBoundaryLoops()) + " boundary loops)");
  }

  // For a cut-mesh element, the corresponding pre-cut-mesh halfedge (an
  // interior halfedge, whose parent is set).
  auto parentInterior = [&](Halfedge cutHe) -> Halfedge {
    if (cutHe.isInterior() && parentHe[cutHe] != Halfedge()) return parentHe[cutHe];
    return parentHe[cutHe.twin()];
  };

  // === Transfer scale factors, edge lengths and slit pairing onto the cut
  //     disk. The two sides of a slit share an origin vertex/edge, so they
  //     receive identical scale factors and a shared seam id. Also verifies the
  //     disk-and-cones-on-boundary precondition: every copy of a cone vertex
  //     must be on the cut boundary (else the cone's prescribed angle deficit
  //     never reaches the boundary layout, silently producing a wrong shape). ===
  VertexData<char> isConeOrig(mesh, 0);
  for (Vertex c : cones) isConeOrig[c] = 1;

  size_t ncV = cutMesh->nVertices();
  VertexData<size_t> cVIdx = cutMesh->getVertexIndices();

  Vector<double> aCut(ncV);
  for (Vertex cv : cutMesh->vertices()) {
    Halfedge origHe;
    for (Halfedge he : cv.outgoingHalfedges()) {
      Halfedge p = parentInterior(he);
      if (p != Halfedge()) {
        origHe = p;
        break;
      }
    }
    if (!cv.isBoundary() && isConeOrig[origHe.vertex()]) {
      throw std::runtime_error("flattenWithGivenConesAndCut: a cone vertex has an interior copy after cutting -- "
                               "`cutEdges` does not fully separate every cone onto the cut boundary");
    }
    aCut(cVIdx[cv]) = aUncut(iVIdx[origHe.vertex()]);
  }

  EdgeData<double> cutLen(*cutMesh, 0.0);
  EdgeData<int> seamId(*cutMesh, -1);
  for (Edge ce : cutMesh->edges()) {
    Halfedge p = parentInterior(ce.halfedge());
    cutLen[ce] = geo.edgeLengths[p.edge()];
    if (ce.isBoundary()) seamId[ce] = static_cast<int>(p.edge().getIndex());
  }

  std::unique_ptr<EdgeLengthGeometry> cutGeom(new EdgeLengthGeometry(*cutMesh, cutLen));
  cutGeom->requireCotanLaplacian();
  cutGeom->requireVertexAngleSums();

  // === Boundary curvatures of the cut disk (the BFF "k - du/dn"). ===
  Vector<bool> isInteriorC(ncV);
  for (Vertex v : cutMesh->vertices()) isInteriorC(cVIdx[v]) = !v.isBoundary();
  BlockDecompositionResult<double> decompC = blockDecomposeSquare(cutGeom->cotanLaplacian, isInteriorC);

  Vector<double> aCutI, g; // interior scale factors, boundary scale factors
  decomposeVector(decompC, aCut, aCutI, g);

  // du/dn = -(A_ib^T a_i + A_bb g)
  Vector<double> dudn = -(decompC.AB.transpose() * aCutI + decompC.BB * g);

  // boundary geodesic curvature k = pi - angleSum, in block (B) order
  Vector<double> kBdy(g.rows());
  for (Vertex v : cutMesh->vertices()) {
    if (!v.isBoundary()) continue;
    size_t bi = decompC.newInds(cVIdx[v]); // index within the boundary block
    kBdy(bi) = M_PI - cutGeom->vertexAngleSums[v];
  }
  Vector<double> ktildeBdy = kBdy - dudn;

  // === Lay out the boundary, with seamless (paired) closure so the two sides
  //     of each slit get equal length. ===
  // Traverse the single boundary loop in order.
  std::vector<Vertex> loopV;
  std::vector<Edge> loopE; // loopE[k] connects loopV[k] -> loopV[k+1]
  for (BoundaryLoop bl : cutMesh->boundaryLoops()) {
    for (Halfedge he : bl.adjacentHalfedges()) {
      loopV.push_back(he.vertex());
      loopE.push_back(he.edge());
    }
    break; // a disk has exactly one boundary loop
  }
  size_t m = loopV.size();

  auto bIdxOf = [&](Vertex v) { return decompC.newInds(cVIdx[v]); };

  // target lengths and tangents per boundary edge
  Vector<double> lstar(m), Tx(m), Ty(m);
  double phi = 0.0;
  for (size_t k = 0; k < m; k++) {
    Vertex va = loopV[k];
    Vertex vb = loopV[(k + 1) % m];
    double ua = g(bIdxOf(va));
    double ub = g(bIdxOf(vb));
    double l = cutLen[loopE[k]];
    lstar(k) = std::exp(0.5 * (ua + ub)) * l;
    phi += ktildeBdy(bIdxOf(va)); // turn at the edge's start vertex
    Tx(k) = std::cos(phi);
    Ty(k) = std::sin(phi);
  }

  // Group edges by seam id (slit pairs share an id); ordinary boundary edges are
  // singletons. Accumulate a per-group length, inverse-mass and tangent.
  std::unordered_map<int, int> groupOf;
  std::vector<double> gL, gNinv, gTx, gTy;
  std::vector<int> edgeGroup(m);
  for (size_t k = 0; k < m; k++) {
    int sid = seamId[loopE[k]];
    int gi;
    auto it = groupOf.find(sid);
    if (sid >= 0 && it != groupOf.end()) {
      gi = it->second;
    } else {
      gi = static_cast<int>(gL.size());
      if (sid >= 0) groupOf[sid] = gi;
      gL.push_back(0.0);
      gNinv.push_back(0.0);
      gTx.push_back(0.0);
      gTy.push_back(0.0);
    }
    edgeGroup[k] = gi;
    gL[gi] = lstar(k); // paired members share equal lstar (seamless scale factors)
    gNinv[gi] += 1.0 / cutLen[loopE[k]];
    gTx[gi] += Tx(k);
    gTy[gi] += Ty(k);
  }
  size_t eN = gL.size();

  // Closure projection in the reduced (per-group) system:
  //   L <- L - Ninv T^T (T Ninv T^T)^{-1} (T L)
  for (size_t i = 0; i < eN; i++) gNinv[i] = 1.0 / gNinv[i];
  double TL0 = 0.0, TL1 = 0.0;       // T * L
  double m00 = 0.0, m01 = 0.0, m11 = 0.0; // T Ninv T^T (symmetric 2x2)
  for (size_t i = 0; i < eN; i++) {
    TL0 += gTx[i] * gL[i];
    TL1 += gTy[i] * gL[i];
    m00 += gTx[i] * gNinv[i] * gTx[i];
    m01 += gTx[i] * gNinv[i] * gTy[i];
    m11 += gTy[i] * gNinv[i] * gTy[i];
  }
  double det = m00 * m11 - m01 * m01;
  if (std::abs(det) > 1e-12) {
    // solve [m] s = TL
    double s0 = (m11 * TL0 - m01 * TL1) / det;
    double s1 = (-m01 * TL0 + m00 * TL1) / det;
    for (size_t i = 0; i < eN; i++) {
      gL[i] -= gNinv[i] * (gTx[i] * s0 + gTy[i] * s1);
    }
  }

  // Accumulate boundary positions from the seamless edge lengths.
  std::vector<double> posX(m), posY(m);
  double re = 0.0, im = 0.0;
  for (size_t k = 0; k < m; k++) {
    posX[k] = -re; // minus sign: boundary is traversed clockwise
    posY[k] = im;
    double L = gL[edgeGroup[k]];
    re += L * Tx(k);
    im += L * Ty(k);
  }

  // === Harmonically extend the boundary positions into the interior. ===
  Vector<double> bdyX = Vector<double>::Zero(ncV), bdyY = Vector<double>::Zero(ncV);
  for (size_t k = 0; k < m; k++) {
    bdyX(cVIdx[loopV[k]]) = posX[k];
    bdyY(cVIdx[loopV[k]]) = posY[k];
  }
  Vector<double> zeroFull = Vector<double>::Zero(ncV);
  Vector<double> X = dirichletSolve(cutGeom->cotanLaplacian, isInteriorC, zeroFull, bdyX);
  Vector<double> Y = dirichletSolve(cutGeom->cotanLaplacian, isInteriorC, zeroFull, bdyY);

  ConeParameterizationResult result;
  result.uvs = VertexData<Vector2>(*cutMesh);
  for (Vertex v : cutMesh->vertices()) {
    result.uvs[v] = Vector2{X(cVIdx[v]), Y(cVIdx[v])};
  }
  result.cones = cones;
  result.boundarySeamId = seamId;
  result.parentHalfedge = std::move(parentHe);
  result.cutMesh = std::move(cutMesh);
  result.cutGeometry = std::move(cutGeom);

  geo.unrequireCotanLaplacian();
  geo.unrequireVertexAngleSums();
  geo.unrequireEdgeLengths();

  return result;
}

ConeParameterizationResult parameterizeBFFwithCones(ManifoldSurfaceMesh& mesh,
                                                    IntrinsicGeometryInterface& geo, size_t nCones) {

  // === 1. Place cones and route geodesic slits to the boundary. ===
  ConePlacementResult cp = computeConePlacement(mesh, geo, nCones);
  ConeCutResult cc = computeConeCut(mesh, geo, cp.cones);

  ManifoldSurfaceMesh& iMesh = cc.network->mesh; // intrinsic mesh the cut lives on
  IntrinsicGeometryInterface& iGeo = *cc.network->tri;

  // Cones and their prescribed angles, transferred to the intrinsic mesh by
  // vertex index (constructFromEdgeSet's triangulation shares mesh's vertex set).
  std::vector<Vertex> iCones;
  iCones.reserve(cp.cones.size());
  VertexData<double> coneAngleI(iMesh, 0.0);
  for (Vertex c : cp.cones) {
    Vertex ic = iMesh.vertex(c.getIndex());
    iCones.push_back(ic);
    coneAngleI[ic] = cp.coneAngles[c];
  }

  EdgeData<char> cutChar(iMesh, 0);
  for (Edge e : iMesh.edges())
    if (cc.cutEdges[e]) cutChar[e] = 1;

  // === 2-3. Flatten the cut disk with the cones' curvature prescribed. ===
  ConeParameterizationResult result = flattenWithGivenConesAndCut(iMesh, iGeo, iCones, coneAngleI, cutChar);
  result.cones = cp.cones; // report on the INPUT mesh, not the intrinsic one
  return result;
}

} // namespace surface
} // namespace geometrycentral
