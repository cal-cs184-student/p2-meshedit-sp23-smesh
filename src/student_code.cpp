#include "student_code.h"
#include "mutablePriorityQueue.h"

using namespace std;

namespace CGL
{

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (class member).
   *
   * @param points A vector of points in 2D
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector2D> BezierCurve::evaluateStep(std::vector<Vector2D> const &points)
  { 
    // TODO Part 1.
//    intermediate_points = []
//    for i in range(len(points)-1):
//        p_i_new = (1-t) * points[i] + t * points[i+1];
//        intermediate_points.append(p_i_new);

    std::vector<Vector2D> intermediate_points;
    for (int i = 0; i < points.size()-1; i++) {
        Vector2D p_i_new = (1 - BezierCurve::t) * points[i] + BezierCurve::t * points[i+1];
        intermediate_points.push_back(p_i_new);
    }

//    return std::vector<Vector2D>();
    return intermediate_points;
  }

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (function parameter).
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector3D> BezierPatch::evaluateStep(std::vector<Vector3D> const &points, double t) const
  {
    // TODO Part 2.
    std::vector<Vector3D> intermediate_points;

    for (int i = 0; i < points.size()-1; i++) {
      Vector3D p_i_new = (1 - t) * points[i] + t * points[i+1];
      intermediate_points.push_back(p_i_new);
    }

    return intermediate_points;
  }

  /**
   * Fully evaluates de Casteljau's algorithm for a vector of points at scalar parameter t
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate1D(std::vector<Vector3D> const &points, double t) const
  {
    // TODO Part 2.

    std::vector<Vector3D> curr_intermediate_points = evaluateStep(points, t);
    while (curr_intermediate_points.size() > 1) {
        curr_intermediate_points = evaluateStep(curr_intermediate_points, t);
    }

//    std::cout << curr_intermediate_points[0] << endl;
    return curr_intermediate_points[0];
  }


  /**
   * Evaluates the Bezier patch at parameter (u, v)
   *
   * @param u         Scalar interpolation parameter
   * @param v         Scalar interpolation parameter (along the other axis)
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate(double u, double v) const 
  {  
    // TODO Part 2.
    std::vector<Vector3D> final_control_points;
    for (int i = 0; i < BezierPatch::controlPoints.size(); i++) {
        Vector3D p_i = evaluate1D(BezierPatch::controlPoints[i], u);
        final_control_points.push_back(p_i);
    }

    Vector3D final_point = evaluate1D(final_control_points, v);
    return final_point;
  }


  Vector3D Vertex::normal( void ) const
  {
    // TODO Part 3.
    // Returns an approximate unit normal at this vertex, computed by
    // taking the area-weighted average of the normals of neighboring
    // triangles, then normalizing.
    HalfedgeCIter h = halfedge();
    Vector3D wnormal(0,0,0);

    do {
        // edge case: do no include boundary surfaces
        if (!h->face()->isBoundary()) {
            // get the vertices of triangle and area
            Vector3D v_0 = position;
            Vector3D v_1 = h->next()->vertex()->position;
            Vector3D v_2 = h->next()->next()->vertex()->position;

            double area = cross(v_0, v_1).norm() / 2; // norm of the cross product gives the area of the parallelogram
            wnormal += h->face()->normal() * area;
        }
        h = h->twin()->next();
    } while (h != halfedge());

//    return Vector3D();
    return wnormal.unit();
  }


  EdgeIter HalfedgeMesh::flipEdge( EdgeIter e0 )
  {
    // TODO Part 4.
    // This method should flip the given edge and return an iterator to the flipped edge.
    return EdgeIter();
  }

  VertexIter HalfedgeMesh::splitEdge( EdgeIter e0 )
  {
    // TODO Part 5.
    // This method should split the given edge and return an iterator to the newly inserted vertex.
    // The halfedge of this vertex should point along the edge that was split, rather than the new edges.
    return VertexIter();
  }



  void MeshResampler::upsample( HalfedgeMesh& mesh )
  {
    // TODO Part 6.
    // This routine should increase the number of triangles in the mesh using Loop subdivision.
    // One possible solution is to break up the method as listed below.

    // 1. Compute new positions for all the vertices in the input mesh, using the Loop subdivision rule,
    // and store them in Vertex::newPosition. At this point, we also want to mark each vertex as being
    // a vertex of the original mesh.
    
    // 2. Compute the updated vertex positions associated with edges, and store it in Edge::newPosition.
    
    // 3. Split every edge in the mesh, in any order. For future reference, we're also going to store some
    // information about which subdivide edges come from splitting an edge in the original mesh, and which edges
    // are new, by setting the flat Edge::isNew. Note that in this loop, we only want to iterate over edges of
    // the original mesh---otherwise, we'll end up splitting edges that we just split (and the loop will never end!)
    
    // 4. Flip any new edge that connects an old and new vertex.

    // 5. Copy the new vertex positions into final Vertex::position.

  }
}
