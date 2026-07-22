#pragma once

#include "geometrycentral/surface/common_subdivision.h"
#include "geometrycentral/surface/edge_length_geometry.h"
#include "geometrycentral/surface/embedded_geometry_interface.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/utilities/elementary_geometry.h"


namespace geometrycentral {
namespace surface {

// A general abstract class for representing intrinsic triangulations.
// In particular, this class encodes not just an intrinsic triangulation, but furthermore an intrinsic triangulation
// which sits on top of some original domain. This motivates many additional operations which involve the correspondence
// with the original mesh.
//
//
// Several different underlying intrinsic triangulation datastructures support this paradigm:
// - SignpostIntrinsicTriangulation
// - IntegerCoordinatesIntrinsicTriangulation
// - ExplicitOverlayIntrinsicTriangulation (TODO implement one day)
// - EdgeLengthIntrinsicTriangulation
//
// See the SIGGRAPH 2021 Course "Geometry Processing with Intrinsic Triangulations" by Nicholas Sharp, Mark Gillespie,
// and Keenan Crane for an introduction to these techniques.

// Default maxInsertions for delaunayRefine(): an automatic budget of
// 10 * nFaces + 10000, computed at call time. A divergent refinement (numerics
// or unrefinable input) on an unlimited budget is an infinite loop, so the
// default is finite; pass INVALID_IND explicitly for the old unlimited
// behavior. Hitting the budget is reported via reachedInsertionBudget.
const size_t AUTO_INSERTION_BUDGET = INVALID_IND - 1;

// Result of delaunayRefine(). The refinement loop is guaranteed to return in
// finite time; this struct reports how it terminated and what, if anything, was
// left unrefined, so callers can decide how to proceed (e.g. accept the partial
// result, relax thresholds, or repair the input near the reported faces).
struct DelaunayRefinementResult {

  // The loop ran to natural completion (all work queues drained). If false, one
  // of the flags below explains the early exit.
  bool completed = false;

  // Terminated early because maxInsertions was reached.
  bool reachedInsertionBudget = false;

  // Terminated early because the internal edge-flip budget was exhausted. This
  // indicates numerical inconsistency in the floating-point Delaunay test (a
  // flip cycle); the mesh may not be fully Delaunay.
  bool flipBudgetExhausted = false;

  // Terminated early because a confirmed insert/delete stall was detected: a
  // full window of insertions produced no net growth in vertex count AND no
  // reduction in the number of criterion-violating faces (see
  // refinementStallWindow). The mesh is left valid and Delaunay; the faces
  // that could not be refined are reported in unrefinedFaces.
  bool stallDetected = false;

  // Statistics.
  size_t nFlips = 0;
  size_t nInsertions = 0;        // vertices successfully inserted
  size_t nDeletions = 0;         // previously-inserted vertices removed (Chew's algorithm)
  size_t nRefusedInsertions = 0; // insertions refused by robustness guards: the new vertex
                                 // would have minted a degenerate or below-floor edge,
                                 // snapped onto an existing vertex, or failed to trace

  // Faces which still violate the refinement criterion on exit. Handles are
  // valid until the next mutation/compression of the intrinsic mesh.
  std::vector<Face> unrefinedFaces;

  // All goals met: completed normally with no face violating the criterion.
  bool success() const { return completed && unrefinedFaces.empty(); }
};

class IntrinsicTriangulation : public EdgeLengthGeometry {

public:
  // Construct an intrinsic triangulation which sits atop this input mesh. Initially, the input triangulation will
  // just be a copy of the input mesh.
  IntrinsicTriangulation(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& inputGeom);
  virtual ~IntrinsicTriangulation();

  // ======================================================
  // ======== Core Members
  // ======================================================

  // The underlying surface on which the intrinsic triangulation has been constructed
  ManifoldSurfaceMesh& inputMesh;
  IntrinsicGeometryInterface& inputGeom;

  // The connectivity of the intrinsic triangulation
  // note that somewhat confusingly, there is a .mesh reference which points to this same mesh,
  // inherited from the geometry interface
  std::unique_ptr<ManifoldSurfaceMesh> intrinsicMesh;

  // The geometry of the intrinsic triangulation is defined by the member
  // EdgeData<double> edgeLengths, which we inherit from the IntrinsicGeometryInterface/EdgeLengthGeometry class

  // Vertex locations for the intrinsic triangulation
  // (i.e., for each vertex in the intrinsic triangulation, what is its location on the input surface)
  VertexData<SurfacePoint> vertexLocations;

  // NOTE: To enable use to make efficient use of the surface tracers, this class always automatically updates the
  // halfedgeVectorsInVertex and halfedgeVectorsInFace geometry members. Could remove this requirement if we change the
  // way the tracer works.

  // Marked edges, which cannot be removed.
  // (set to an array which holds true if an edge is fixed, and should not be flipped)
  // A callback is automatically registered which will update this array as edge splits are performed, so if a marked
  // edge is split the two resulting edges will be marked.
  // Note that if no marked edges have been set, this array will be uninitialized; the helpers isFixed() etc account for
  // this possiblity.
  EdgeData<bool> markedEdges;
  void setMarkedEdges(const EdgeData<bool>& markedEdges);
  void clearMarkedEdges();
  // Is this a marked or boundary edge?
  bool isFixed(Edge e) const;
  bool isOnFixedEdge(Vertex v) const; // boundary vertex or on fixed edge

  // Parameters
  double triangleTestEPS = 1e-6; // used for numerical checks in mesh operations

  // === Robustness / termination controls for delaunayRefine().
  //
  // Insertion length floor (packing guard): refinement refuses to insert a
  // vertex whose shortest incident edge would be below
  //   refinementMinRelativeLength * min(shortest edge at call time, circumradiusThresh).
  // In exact arithmetic Chew's algorithm never shrinks the shortest edge (for
  // angle bounds <= 30 degrees), so a generous relative floor only triggers
  // when floating-point error or unrefinable input (e.g. tiny cone angles or
  // sharp fixed-edge wedges) has broken the termination argument; the offending
  // faces are reported in DelaunayRefinementResult::unrefinedFaces instead of
  // looping forever. Set to 0 to disable the guard.
  double refinementMinRelativeLength = 1e-3;

  // Stall guard: refinement evaluates progress over tumbling windows of this
  // many insertions, and declares a stall (stops, reports stallDetected + the
  // unrefined faces) only when consecutive windows show no progress on ALL
  // THREE of: (1) net growth of the vertices inserted by this call, (2)
  // deletion of vertices that predate this call (each drains a finite pool, so
  // it is progress by itself -- restructuring refines on constrained meshes
  // are legitimately deletion-dominant with little or no net growth), and (3)
  // the number of criterion-violating faces (checked by an O(n) sweep run at
  // most once per window, only when (1) and (2) already show no progress).
  // Set to 0 to disable.
  size_t refinementStallWindow = 100;


  // ======================================================
  // ======== Queries & Accessors
  // ======================================================

  // Trace out the edges of the intrinsic triangulation along the surface of the input mesh.
  // Each path is ordered along edge.halfedge(), and includes both the start and end points
  virtual EdgeData<std::vector<SurfacePoint>> traceAllIntrinsicEdgesAlongInput();
  // Trace just one halfedge, ordered in the direction the halfedge points.
  virtual std::vector<SurfacePoint> traceIntrinsicHalfedgeAlongInput(Halfedge intrinsicHe) = 0;

  // Trace out the edges of the intrinsic triangulation along the surface of the input mesh.
  // Each path is ordered along edge.halfedge(), and includes both the start and end points
  virtual EdgeData<std::vector<SurfacePoint>> traceAllInputEdgesAlongIntrinsic();
  // Trace just one halfedge, ordered in the direction the halfedge points.
  virtual std::vector<SurfacePoint> traceInputHalfedgeAlongIntrinsic(Halfedge inputHe) = 0;

  // Get a reference to the common subdivision. May construct it from scratch if this is the first time
  // it is needed. The intrinsic triangulation manages the lifetime of the subdivision---it will be deallocated if (a)
  // this object is deleted, or (b) the triangulation is mutated, invalidating the common subdivision. Be sure to copy
  // it if you want to retain it through those operations.
  // WARNING: constructing the subdivision compresses the intrinsic mesh, which INVALIDATES all
  // intrinsic element handles (Vertex/Edge/Face/...) held by the caller. Re-acquire any handles
  // after calling this (sample them afterwards, or store indices/positions instead), exactly as
  // after an explicit mesh.compress().
  CommonSubdivision& getCommonSubdivision();

  // Given a point on the input triangulation, returns the corresponding point on the intrinsic triangulation
  virtual SurfacePoint equivalentPointOnIntrinsic(const SurfacePoint& pointOnInput) = 0;

  // Given a point on the intrinsic triangulation, returns the corresponding point on the input triangulation
  virtual SurfacePoint equivalentPointOnInput(const SurfacePoint& pointOnIntrinsic) = 0;

  // Given data defined on the vertices of the input triangulation, samples it to the vertices of the intrinsic
  // triangulation
  template <typename T>
  VertexData<T> sampleFromInput(const VertexData<T>& dataOnInput);

  // Given data defined on the vertices of the intrinsic triangulation, restrict it to the vertices of the input
  // triangulation
  template <typename T>
  VertexData<T> restrictToInput(const VertexData<T>& dataOnIntrinsic);

  // Returns true if the intrinsic triangulation (or edge) satisifies the intrinsic Delaunay criterion
  bool isDelaunay();
  bool isDelaunay(Edge e);

  // Returns the smallest angle in the intrinsic triangulation, in degrees
  double minAngleDegrees() const;

  // Minimum angle (degrees) between two fixed (marked or boundary) edges that
  // meet at a vertex -- the smallest angular sector spanned between consecutive
  // fixed edges around any vertex. This is an UPPER BOUND on the min angle any
  // delaunayRefine can achieve: a wedge bounded by two fixed edges cannot be
  // widened, so no refinement exceeds it. With only boundary edges fixed (e.g.
  // Stage A) it reports the sharpest boundary corner. Returns -1 if no vertex
  // has two incident fixed edges. Invariant: minFixedAngleDegrees >= minAngle.
  double minFixedAngleDegrees() const;

  // Returns the smallest angle in the intrinsic triangulation, in degrees, among faces whose vertices
  // all have angle sum at least minAngleSum, and which are not contained in an original face incident
  // on a vertex with small angle sum
  double minAngleDegreesAtValidFaces(double minAngleSum) const;

  // If f is entirely contained in some face of the input mesh, return that
  // face Otherwise return Face()
  Face getParentFace(Face f) const;

  // Check if edge is shared with input mesh
  virtual bool checkEdgeOriginal(Edge e) const = 0;

  // Immediate computation of corner angle
  double getCornerAngle(Corner c) const;

  // ======================================================
  // ======== High-Level Mutators
  // ======================================================
  //
  // Call once to build a useful triangulation

  // Flips edges in the intrinsic triangulation until is satisfies the intrinsic Delaunay criterion
  void flipToDelaunay();

  // Perform intrinsic Delaunay refinement the intrinsic triangulation until it simultaneously:
  //   - satisfies the intrinsic Delaunay criterion
  //   - has no angles smaller than `angleThreshDegrees` (values > 30 degrees have no
  //     termination guarantee even in exact arithmetic)
  //   - has no triangles larger than `circumradiusThresh`
  // Terminates no matter what after maxInsertions insertions. The default is an automatic
  // finite budget (see AUTO_INSERTION_BUDGET); pass INVALID_IND for unlimited.
  //
  // Always returns in finite time (given a finite insertion budget): when the goals
  // cannot be met (numerical trouble, unrefinable input corners, budget), the robustness
  // guards (see refinementMinRelativeLength / refinementStallWindow above) stop the loop
  // and the returned DelaunayRefinementResult reports what happened and which
  // faces remain unrefined, rather than iterating forever.
  DelaunayRefinementResult delaunayRefine(double angleThreshDegrees = 25.,
                                          double circumradiusThresh = std::numeric_limits<double>::infinity(),
                                          size_t maxInsertions = AUTO_INSERTION_BUDGET);


  // General version of intrinsic Delaunay refinement, taking a function which will be called
  // to determine if a triangle should be refined.
  // Returns when all triangles pass this function, when maxInsertions is exceeded, or when
  // the robustness guards conclude no further progress can be made (see the returned
  // DelaunayRefinementResult).
  // lengthFloor is the absolute insertion length floor (packing guard); if negative it is
  // computed as refinementMinRelativeLength * (shortest edge at call time).
  DelaunayRefinementResult delaunayRefine(const std::function<bool(Face)>& shouldRefine,
                                          size_t maxInsertions = AUTO_INSERTION_BUDGET, double lengthFloor = -1.);


  // ======================================================
  // ======== Low-Level Mutators
  // ======================================================
  //
  // Basic operations to modify the intrinsic triangulation
  // NOTE: The individual operations to not call refreshQuantities(), so you should call it if you want quantities
  // updated.

  // If the edge is not Delaunay, flip it. Returns true if flipped.
  virtual bool flipEdgeIfNotDelaunay(Edge e) = 0;

  // If the edge can be flipped, flip it (must be combinatorially flippable and inside a convex quad). Returns true if
  // flipped.
  virtual bool flipEdgeIfPossible(Edge e) = 0;

  // Insert a new vertex in to the intrinsic triangulation
  virtual Vertex insertVertex(SurfacePoint newPositionOnIntrinsic) = 0;

  // Insert the circumcenter of a face in to the triangulation. Returns the newly created intrinsic vertex.
  Vertex insertCircumcenter(Face f);

  // Insert the barycenter of a face in to the triangulation. Returns the newly created intrinsic vertex.
  Vertex insertBarycenter(Face f);

  // Remove an (inserted) vertex from the triangulation.
  // Note: if something goes terribly (numerically?) wrong, will exit without removing the vertex.
  virtual Face removeInsertedVertex(Vertex v) = 0;

  // Split an edge
  virtual Halfedge splitEdge(Halfedge he, double tSplit) = 0;


  // ==== Misc
  // Recover t-values after tracing
  // Note that really we ought to just report these back from the tracing routine itself, which computes them
  // internally. We don't have a nice API for passing that data around, so this lazily recovers it
  std::vector<double> recoverTraceTValues(const std::vector<SurfacePoint>& edgeTrace);

  // ======================================================
  // ======== Callbacks
  // ======================================================
  //
  // Get called whenever mesh mutations occur. Register a callback by inserting it in to this list.
  //

  // edge E if flipped
  std::list<std::function<void(Edge)>> edgeFlipCallbackList;

  // old face F is split by new vertex V
  std::list<std::function<void(Face, Vertex)>> faceInsertionCallbackList;

  // old edge E is split to halfedge HE1,HE2 both with he.vertex() as split vertex
  std::list<std::function<void(Edge, Halfedge, Halfedge)>> edgeSplitCallbackList;


  // ======================================================
  // ======== Geometry Immediates
  // ======================================================

  // Computed on the intrinsic triangulation
  double shortestEdge(Face f) const;

  // Isometrically lay out the vertices around a halfedge in 2D coordinates
  // he points from vertex 2 to 0; others are numbered CCW
  std::array<Vector2, 4> layoutDiamond(Halfedge he);
  std::array<Vector2, 3> vertexCoordinatesInTriangle(Face face);

  // Repopulate the member halfedgeVectorInFace
  void updateFaceBasis(Face f); // assumes edgeLengths are valid


protected:
  // Absolute insertion length floor currently in force. Set (and restored to -1) by
  // delaunayRefine(); consulted by insertCircumcenter(), which refuses insertions that
  // would mint an edge shorter than this. <= 0 means no floor is active.
  double activeRefinementLengthFloor = -1.;

  // The current common subdivision. This member will only be populated if the subdivision is valid.
  // Implementations must sure to call triangulationChanged() below after any modifications, which handles clearing out
  // the common subdivision.
  std::unique_ptr<CommonSubdivision> commonSubdivision;
  virtual void constructCommonSubdivision() = 0; // ensure commonSubdivision is populated

  // Must be called any time the intrinsic triangulation is modified.
  void triangulationChanged();

  // Callback helpers
  void invokeEdgeFlipCallbacks(Edge e);
  void invokeFaceInsertionCallbacks(Face f, Vertex v);
  void invokeEdgeSplitCallbacks(Edge e, Halfedge he1, Halfedge he2);
};


} // namespace surface
} // namespace geometrycentral


#include "geometrycentral/surface/intrinsic_triangulation.ipp"
