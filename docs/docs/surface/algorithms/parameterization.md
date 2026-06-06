This section describes the algorithms in geometry-central for surface parameterization, which compute maps from surfaces meshes to the plane.

Note that these procedures all depend on the _intrinsic_ geometry of a surface (via the `IntrinsicGeometryInterface`). Therefore, these routines can be run on abstract geometric domains as well as traditional surfaces in 3D.

## Boundary First Flattening
![parameterized face](/media/bff.png){: style="max-width: 15em; display: block; margin-left: auto; margin-right: auto;"}
This algorithm is described in the paper [Boundary First Flattening](http://www.cs.cmu.edu/~kmcrane/Projects/BoundaryFirstFlattening/paper.pdf). It computes a _conformal_ parameterization of a surface mesh, allowing the user to specify target angles or scale factors along the boundary of the mesh. The input mesh must be a topological disk.

`#include "geometrycentral/surface/boundary_first_flattening.h"`

### Single Parameterizations

A one-off utility function is provided to compute single parameterizations. Repeated parameterizations of the same mesh can be computed more efficiently using the utility class `BFF` below.

Example:
```cpp
#include "geometrycentral/surface/boundary_first_flattening.h"
#include "geometrycentral/surface/meshio.h"

// Load a mesh
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;
std::tie(mesh, geometry) = readManifoldSurfaceMesh(filename);

VertexData<Vector2> parameterization = parameterizeBFF(*mesh, *geometry);
```

??? func "`#!cpp VertexData<Vector2> parameterizeBFF(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom)`"
    Conformally parameterize the input mesh. Picks boundary conditions to minimize area distortion (i.e. sets the conformal scale factor to 0 along the boundary).

??? func "`#!cpp VertexData<Vector2> parameterizeBFFfromScaleFactors(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom, const VertexData<double>& boundaryScaleFactors)`"
    Conformally parameterize the input mesh, setting the scale factors to the given values along the boundary.
    Although `boundaryScaleFactors` is a `VertexData` object which stores values at all vertices, only the values at boundary vertices are used by the algorithm. All other values are ignored.
    
??? func "`#!cpp VertexData<Vector2> parameterizeBFFfromExteriorAngles(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom, const VertexData<double>& exteriorAngles)`"
    Conformally parameterize the input mesh so that the boundary vertices of the parameterized mesh have the given exterior angles. The exterior angles must sum up to $2\pi$ along the boundary.
    
    Although `exteriorAngles` is a `VertexData` object which stores values at all vertices, only the values at boundary vertices are used by the algorithm. All other values are ignored.

### Repeated Parameterization

The stateful class `BFF` does precomputation when constructed to efficiently compute many parameterizations of the same mesh.

Example:
```cpp
#include "geometrycentral/surface/boundary_first_flattening.h"
#include "geometrycentral/surface/meshio.h"

// Load a mesh
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;
std::tie(mesh, geometry) = readManifoldSurfaceMesh(filename);

VertexData<double> boundaryScaleFactors = /* target scale factors */;
VertexData<double> exteriorAngles = /* target exteriorAngles */;

BFF bff(*mesh, *geometry);
VertexData<Vector2> parameterization1 = bff.flattenFromScaleFactors(boundaryScaleFactors);
VertexData<Vector2> parameterization2 = bff.flattenFromExteriorAngles(exteriorAngles);
```

??? func "`#!cpp BFF::BFF(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom)`"

    Create a new solver for boundary first flattening. Most precomputation is done immediately when the object is constructed, although some additional precomputation may be done lazily later on.

??? func "`#!cpp VertexData<Vector2> BFF::flatten()`"

    Compute a conformal parameterization which minimizes area distortion (i.e. sets the scale factor to 0 along the boundary).
    
??? func "`#!cpp VertexData<Vector2> BFF::flattenFromScaleFactors(const VertexData<double>& boundaryScaleFactors)`"

    Compute a conformal parameterization with the given scale factor along the boundary.
    Although `boundaryScaleFactors` is a `VertexData` object which stores values at all vertices, only the values at boundary vertices are used by the algorithm. All other values are ignored.
    
??? func "`#!cpp VertexData<Vector2> BFF::flattenFromExteriorAngles(const VertexData<double>& exteriorAngles)`"

    Compute a conformal parameterization with the given exterior angles along the boundary. The exterior angles must sum up to $2\pi$ along the boundary.
    
    Although `exteriorAngles` is a `VertexData` object which stores values at all vertices, only the values at boundary vertices are used by the algorithm. All other values are ignored.
    
### Citation

If this algorithm contributes to academic work, please cite the following paper:

```bib
@article{Sawhney:2017:BFF,
author = {Sawhney, Rohan and Crane, Keenan},
title = {Boundary First Flattening},
journal = {ACM Transactions on Graphics (TOG)},
volume = {37},
number = {1},
month = dec,
year = {2017},
issn = {0730-0301},
pages = {5:1--5:14},
articleno = {5},
numpages = {14},
url = {http://doi.acm.org/10.1145/3132705},
doi = {10.1145/3132705},
acmid = {3132705},
publisher = {ACM},
address = {New York, NY, USA}
}
```

## Cone Parameterization

Surfaces with significant Gaussian curvature cannot be flattened to the plane with low area distortion. _Cone parameterization_ addresses this by introducing a small number of _cone singularities_, which concentrate the curvature at isolated points, and cutting the surface open along geodesic slits connecting the cones (and the boundary) so it becomes a topological disk. The disk is then flattened with [Boundary First Flattening](#boundary-first-flattening), prescribing the cone curvatures and laying out the boundary so that the **two sides of every cut have equal length** -- a _seamless_ map, whose cuts can be re-joined without stretch mismatch. This follows the cone strategy of the [Boundary First Flattening](http://www.cs.cmu.edu/~kmcrane/Projects/BoundaryFirstFlattening/paper.pdf) paper.

These routines depend only on the _intrinsic_ geometry (via the `IntrinsicGeometryInterface`).

### Cone Placement

`#include "geometrycentral/surface/cone_placement.h"`

??? func "`#!cpp ConePlacementResult computeConePlacement(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geo, size_t nCones)`"
    Automatically select up to `nCones` interior cone singularities which minimize area distortion, using the greedy "CETM" strategy: starting from an anchor set (all boundary vertices for a surface with boundary, or the extreme-curvature vertex for a closed surface), repeatedly solve for the conformal scale-factor field and add the vertex of largest scale factor as the next cone. The mesh may be closed or have boundary.

    Returns a `ConePlacementResult` with:

    - `#!cpp std::vector<Vertex> cones` -- the chosen interior cone vertices.
    - `#!cpp VertexData<double> coneAngles` -- the prescribed target angle at each cone (and boundary) vertex, normalized so the total prescribed curvature satisfies the discrete Gauss-Bonnet theorem (sum $= 2\pi\chi$); zero at ordinary interior vertices.

### Cone Cutting

`#include "geometrycentral/surface/cone_cut.h"`

??? func "`#!cpp ConeCutResult computeConeCut(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geo, const std::vector<Vertex>& cones)`"
    Compute geodesic slits which reduce the surface to a single topological disk passing through the given cones. The cut _topology_ is an approximate minimum-weight (by edge length) Steiner tree connecting all boundary loops and cone vertices; the cut _geometry_ is then straightened to geodesics by a [`FlipEdgeNetwork`](flip_geodesics.md) on an internally-built intrinsic triangulation, so the slits are smooth rather than jagged mesh-edge polylines. Boundary edges are never cut. For a closed surface at least two cones are required.

    Returns a `ConeCutResult` with:

    - `#!cpp std::unique_ptr<FlipEdgeNetwork> network` -- owns the intrinsic triangulation the cut is defined on (`network->mesh`), with edge lengths from `network->tri`.
    - `#!cpp EdgeData<bool> cutEdges` -- the geodesic cut edges (keyed on `network->mesh`). Cutting along these with `surgery::cutAlongEdges` yields a disk.

### Seamless Cone Parameterization

`#include "geometrycentral/surface/cone_parameterization.h"`

??? func "`#!cpp ConeParameterizationResult parameterizeBFFwithCones(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geo, size_t nCones)`"
    Conformally flatten a surface with `nCones` automatically-placed cone singularities, producing a seamless map of the resulting cut disk. Composes `computeConePlacement` and `computeConeCut`, then flattens the cut disk with the cone curvatures prescribed and the two sides of each slit constrained to equal length. Requires a surface _with boundary_ (the scale-factor solve is anchored at the original boundary).

    Returns a `ConeParameterizationResult` with:

    - `#!cpp std::unique_ptr<ManifoldSurfaceMesh> cutMesh` and `#!cpp std::unique_ptr<EdgeLengthGeometry> cutGeometry` -- the cut disk and its intrinsic geometry.
    - `#!cpp VertexData<Vector2> uvs` -- the flattening, per vertex of `cutMesh`.
    - `#!cpp std::vector<Vertex> cones` -- the interior cone vertices that were placed (vertices of the input `mesh`).
    - `#!cpp EdgeData<int> boundarySeamId` -- for each boundary edge of `cutMesh`, the id of the original edge it came from; the two sides of a slit share an id, and the flattening guarantees both sides have equal length.

Example:
```cpp
#include "geometrycentral/surface/cone_parameterization.h"
#include "geometrycentral/surface/meshio.h"

// Load a mesh with boundary (e.g. a disk patch)
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;
std::tie(mesh, geometry) = readManifoldSurfaceMesh(filename);

// Flatten with 4 cones, seamlessly
ConeParameterizationResult result = parameterizeBFFwithCones(*mesh, *geometry, 4);
// result.uvs holds the UV coordinates on result.cutMesh
```

These routines use the same citation as Boundary First Flattening above.
