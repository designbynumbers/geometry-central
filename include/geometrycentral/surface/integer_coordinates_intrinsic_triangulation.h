#pragma once

#include "geometrycentral/surface/barycentric_coordinate_helpers.h"
#include "geometrycentral/surface/common_subdivision.h"
#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/intrinsic_mollification.h"
#include "geometrycentral/surface/intrinsic_triangulation.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/normal_coordinates.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/surface/trace_geodesic.h"
#include "geometrycentral/utilities/elementary_geometry.h"

#include <deque>

namespace geometrycentral {
namespace surface {

class IntegerCoordinatesIntrinsicTriangulation : public IntrinsicTriangulation {

public:
  // Construct an intrinsic triangulation which sits atop this input mesh.
  // Initially, the input triangulation will just be a copy of the input mesh.
  IntegerCoordinatesIntrinsicTriangulation(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& inputGeom,
                                           double mollifyEPS = 1e-5);

  // ======================================================
  //                   Core Members
  // ======================================================

  // The actual normal coordinates (and roundabouts) encoding the triangulation. These normal coordinates are defined on
  // top of the _intrinsic_ mesh---for each intrinsic edge, the encode how many original edges cross it.
  NormalCoordinates normalCoordinates;

  // insertVertex() refuses to insert a vertex which coincides with an
  // existing one to within this tolerance (in barycentric parameter), and
  // returns the existing vertex instead. Near-coincident vertices create
  // near-zero-length intrinsic edges, whose degenerate geometry makes
  // subsequent floating point classifications unreliable. Set to 0 to
  // disable. (This refuses an operation; it never alters recorded data.)
  double insertionCoincidenceEPS = 1e-12;

  // ======================================================
  // ======== Queries & Accessors
  // ======================================================

  // See intrinsic_triangulation.h for docs.

  SurfacePoint equivalentPointOnIntrinsic(const SurfacePoint& pointOnInput) override;

  SurfacePoint equivalentPointOnInput(const SurfacePoint& pointOnIntrinsic) override;

  EdgeData<std::vector<SurfacePoint>> traceAllIntrinsicEdgesAlongInput() override;
  std::vector<SurfacePoint> traceIntrinsicHalfedgeAlongInput(Halfedge intrinsicHe) override;

  EdgeData<std::vector<SurfacePoint>> traceAllInputEdgesAlongIntrinsic() override;
  std::vector<SurfacePoint> traceInputHalfedgeAlongIntrinsic(Halfedge inputHe) override;

  bool checkEdgeOriginal(Edge e) const override;

  // ======================================================
  // ======== Low-Level Mutators
  // ======================================================

  // See intrinsic_triangulation.h for docs.

  bool flipEdgeIfNotDelaunay(Edge e) override;

  bool flipEdgeIfPossible(Edge e) override;

  Vertex insertVertex(SurfacePoint newPositionOnIntrinsic) override;

  Face removeInsertedVertex(Vertex v) override;

  Halfedge splitEdge(Halfedge he, double tSplit) override;

  // Check if an edge can be flipped geometrically, as defined by the (relative) signed areas of the resulting
  // triangles; positive values mean flippable.
  double checkFlip(Edge e);

  // Insert circumcenter or split segment
  Vertex insertCircumcenterOrSplitSegment(Face f, bool verbose = false);

  // Insert a vertex exactly AT the crossingIndex'th transverse crossing of
  // an input-edge curve along he (crossings indexed from he's tail, as in
  // computeEdgeSplitData). Unlike splitEdge -- whose combinatorial
  // classification always places the new vertex strictly *between*
  // crossings -- this cuts the curve at the new vertex: the curve's two
  // halves terminate there, and the vertex's location is recorded as a
  // point ON the corresponding input edge. This is the principled way to
  // obtain an intrinsic vertex lying on an input edge at a chosen crossing.
  Vertex insertVertexAtCrossing(Halfedge he, int crossingIndex, bool verbose = false);

  Vertex splitFace(Face f, Vector3 bary, bool verbose = false);
  Vertex splitEdge(Edge e, double bary, bool verbose = false);
  Halfedge splitInteriorEdge(Halfedge he, double bary, bool verbose = false);
  Halfedge splitBoundaryEdge(Halfedge he, double bary, bool verbose = false);

  // Move a vertex `v` in direction `vec`, represented as a vector in the
  // vertex's tangent space.
  Vertex moveVertex(Vertex v, Vector2 vec);

  // Assumes intrinsicEdgeLengths is up to date
  void updateCornerAngle(Corner c);

  // Assumes cornerAngles, vertexAngleSums exist and are up to date
  void updateHalfedgeVectorsInVertex(Vertex v);

  // ======================================================
  //                Low-Level Queries
  // ======================================================

  // Takes in a halfedge of the intrinsic mesh whose edge's normal coordinate
  // is negative (meaning that it lies along an edge of the input mesh) and
  // returns the halfedge in the input mesh pointing in the same direction
  // e.vertex() must live in both meshes
  Halfedge getSharedInputEdge(Halfedge e) const;

  // Takes in an intrinsic point, represented as an intrinsic face and barycentric coordinate,
  // and computes the corresponding point on the input mesh, as well as the normal coordinates of
  // the edges connecting this new point to f's vertices.
  std::pair<SurfacePoint, std::array<int, 3>> computeFaceSplitData(Face f, Vector3 bary, bool verbose = false);

  // Computes the corresponding point on the input mesh, as well as which segment the point belongs to.
  // The segment is specified by an integer 0 <= segmentIndex <= normalCoordinates[he.edge()]
  std::pair<SurfacePoint, size_t> computeEdgeSplitData(Halfedge he, double tBary);

  // Compute the number of vertices in the common subdivision
  // i.e. intrinsicMesh->nVertices() plus the sum of the normal coordinates
  size_t nSubdividedVertices() const;

  // HACK: represents arcs parallel to a mesh edge with a single pair {-n,
  // he} where n is the number of arcs parallel to he.edge() Trace an edge
  // of the input mesh over the intrinsic triangulation
  NormalCoordinatesCompoundCurve traceInputEdge(Edge e, bool verbose = false) const;
  NormalCoordinatesCompoundCurve traceInputHalfedge(Halfedge inputHe, bool verbose = false) const;

  std::pair<bool, NormalCoordinatesCurve> traceNextCurve(const NormalCoordinatesCurve& oldCurve,
                                                         bool verbose = false) const;

  // Inverse to traceInputEdge
  Halfedge identifyInputEdge(const NormalCoordinatesCurve& path, bool verbose = false) const;

  // Identify the input edge along which a transverse curve component lies,
  // together with the range of the input edge spanned by the component.
  // Returns (inputHe, tStart, tEnd): the component runs along inputHe from
  // tStart to tEnd in inputHe's parameterization. In particular, a point at
  // parameter t along the component lies at tStart + t * (tEnd - tStart)
  // along inputHe.
  // Unlike identifyInputEdge, this handles components which begin or end at
  // inserted vertices lying on the input edge (e.g. vertices added by
  // splitting an edge which runs along an input edge); such components span
  // only part of the input edge.
  std::tuple<Halfedge, double, double> identifyInputCurveRange(const NormalCoordinatesCurve& path) const;

  // Identify shared halfedge, throw exception if halfedge is not shared
  // (i.e. edgeCoords[he.edge()] must be negative)
  Halfedge identifyInputEdge(Halfedge he) const;

  std::array<Vector2, 3> vertexCoordinatesInFace(Face face) const;

  void setFixedEdges(const EdgeData<bool>& fixedEdges);

  // If f is entirely contained in some face of the input mesh, return that
  // face Otherwise return Face()
  Face getParentFace(Face f) const;

private:
  // Implementation details

  // Construct the common subdivision for the current triangulation.
  void constructCommonSubdivision() override;

  // === Exact input-element derivation
  // The element of the input mesh on which a point lies is always determined
  // combinatorially (from normal coordinates, roundabouts, and the
  // already-exact elements of existing vertex locations), never by
  // inspecting floating point coordinate values.

  // Given a transverse curve component and which side of it a region lies
  // on, return the input face containing that region. The component must
  // have been produced by tracing through the halfedge bounding the region
  // (topologicalTrace*(he, ...)), so the trace crosses he positively: he
  // then points to the curve's left, and an input halfedge has its face on
  // its left, making this a purely combinatorial lookup.
  // `regionOnTailSide` says whether the region contains he's tail (true)
  // or its tip (false).
  Face inputFaceBesideCurve(const NormalCoordinatesCurve& curve, bool regionOnTailSide) const;

  // For a shared intrinsic halfedge (normal coordinate < 0), return the
  // input halfedge pointing in the same direction (resolved via roundabouts
  // when both endpoints are original vertices, and via the endpoints'
  // recorded positions along the input edge otherwise).
  Halfedge inputHalfedgeAlongShared(Halfedge he) const;

  // Return the input face containing a maximal connected region of
  // uncrossed intrinsic faces (every intrinsic edge interior to the region
  // has normal coordinate <= 0 and is not shared). Searches the region for
  // an exact anchor: a Face-typed vertex location, a shared boundary edge
  // (resolved via roundabouts), or a crossed boundary edge (resolved via
  // inputFaceBesideCurve).
  Face inputFaceOfUncrossedRegion(Face f) const;

  // The input face containing the interior of an uncrossed (normal
  // coordinate == 0) intrinsic edge.
  Face inputFaceOfUncrossedEdge(Edge e) const;

  // Find the most-preferred flippable edge incident on v (highest checkFlip
  // score, preferring loop edges, skipping fixed edges). Returns Edge() if
  // no incident edge can be flipped.
  Edge bestFlippableEdgeAround(Vertex v);

  // Remove an inserted vertex lying on the boundary (i.e. a vertex created
  // by splitting a boundary edge). Returns the new face, or Face() if the
  // vertex could not be removed.
  Face removeInsertedBoundaryVertex(Vertex v);
};

FaceData<Vector2> interpolateTangentVectorsB(const IntegerCoordinatesIntrinsicTriangulation& tri,
                                             const CommonSubdivision& cs, const FaceData<Vector2>& dataB);


} // namespace surface
} // namespace geometrycentral

#include "integer_coordinates_intrinsic_triangulation.ipp"
