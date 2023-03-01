#include "student_code.h"
#include "mutablePriorityQueue.h"
#include "algorithm"

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
    if (!e0->isBoundary()) {
        // do inner and outer boundaries is the edge is not a boundary
//        list<HalfedgeIter> inner_halfedges;
//        HalfedgeIter h = e0->halfedge();
//        HalfedgeIter start = e0->halfedge();
//
//        inner_halfedges.push_back(h);
//        do {
//            inner_halfedges.push_back(h->next());
//            h = h->next();
//        } while (h != start);
//
//        // get all the outer edges
//        list<HalfedgeIter> outer_halfedges;
//        outer_halfedges.push_back(h->twin());
//        do {
//            outer_halfedges.push_back(h->twin());
//            h = h->next();
//        } while (h != start);

        // inner edges
        HalfedgeIter h_in = e0->halfedge();
        HalfedgeIter h_in_1 = h_in->next();
        HalfedgeIter h_in_2 = h_in_1->next();
        HalfedgeIter h_in_3 = h_in->twin();
        HalfedgeIter h_in_4 = h_in_3->next();
        HalfedgeIter h_in_5 = h_in_4->next();
        //outer edges
        HalfedgeIter h_out_6 = h_in_1->twin();
        HalfedgeIter h_out_7 = h_in_2->twin();
        HalfedgeIter h_out_8 = h_in_4->twin();
        HalfedgeIter h_out_9 = h_in_5->twin();
        // vertices
        VertexIter b = h_in->vertex(); //v0
        VertexIter a = h_in_2->vertex(); //v2
        VertexIter c = h_in_3->vertex(); //v1
        VertexIter d = h_in_5->vertex(); //v3
        // faces
        FaceIter f1 = h_in->face();
        FaceIter f2 = h_in_3->face();
        // edges
        EdgeIter e1 = h_in_1->edge();
        EdgeIter e2 = h_in_2->edge();
        EdgeIter e3 = h_in_4->edge();
        EdgeIter e4 = h_in_5->edge();

        // recreate the mesh
        // next, twin, vertex, edge, face
        h_in->setNeighbors(h_in_1,h_in_3,a , e0, f1);
        h_in_1->setNeighbors(h_in_2,h_out_9,d,e4,f1);
        h_in_2->setNeighbors(h_in,h_out_6,c,e1,f1);
        h_in_3->setNeighbors(h_in_4,h_in,d, e0,f2);
        h_in_4->setNeighbors(h_in_5,h_out_7,a,e2,f2);
        h_in_5->setNeighbors(h_in_3,h_out_8,b,e3,f2);

        h_out_6->setNeighbors(h_out_6->next(),h_in_2,a ,e1,h_out_6->face());
        h_out_7->setNeighbors(h_out_7->next(),h_in_4,b,e2,h_out_7->face());
        h_out_8->setNeighbors(h_out_8->next(),h_in_5, d, e3, h_out_8->face());
        h_out_9->setNeighbors(h_out_9->next(),h_in_1,c,e4,h_out_9->face());

        // Vertices
        b->halfedge() = h_in_5;
        c->halfedge() = h_in_2;
        a->halfedge() = h_in_4;
        d->halfedge() = h_in_1;

        // Edges
        e0->halfedge() = h_in;
        e1->halfedge() = h_in_2;
        e2->halfedge() = h_in_4;
        e3->halfedge() = h_in_5;
        e4->halfedge() = h_in_1;

        // Faces
        f1->halfedge() = h_in;
        f2->halfedge() = h_in_3;

    }

//    return EdgeIter();
    return e0;
  }

  VertexIter HalfedgeMesh::splitEdge( EdgeIter e0 )
  {
    // TODO Part 5.
    // This method should split the given edge and return an iterator to the newly inserted vertex.
    // The halfedge of this vertex should point along the edge that was split, rather than the new edges.
      if (!e0->isBoundary()) {
          // all halfedges
          HalfedgeIter h_in = e0->halfedge();
          HalfedgeIter h_in_1 = h_in->next();
          HalfedgeIter h_in_2 = h_in_1->next();
          HalfedgeIter h_in_3 = h_in->twin();
          HalfedgeIter h_in_4 = h_in_3->next();
          HalfedgeIter h_in_5 = h_in_4->next();

          HalfedgeIter h_out_6 = h_in_1->twin();
          HalfedgeIter h_out_7 = h_in_2->twin();
          HalfedgeIter h_out_8 = h_in_4->twin();
          HalfedgeIter h_out_9 = h_in_5->twin();

          // all edges
          EdgeIter e1 = h_in_1->edge();
          EdgeIter e2 = h_in_2->edge();
          EdgeIter e3 = h_in_4->edge();
          EdgeIter e4 = h_in_5->edge();

          // all faces
          FaceIter f1 = h_in->face();
          FaceIter f2 = h_in_3->face();

          // all vertices
          VertexIter b = h_in->vertex();
          VertexIter c = h_in_3->vertex();
          VertexIter a = h_in_2->vertex();
          VertexIter d = h_in_5->vertex();


          // *************** new elements ***************
          // new inner edges
          HalfedgeIter h_in_10 = newHalfedge();
          HalfedgeIter h_in_11 = newHalfedge();
          HalfedgeIter h_in_12 = newHalfedge();
          HalfedgeIter h_in_13 = newHalfedge();
          HalfedgeIter h_in_14 = newHalfedge();
          HalfedgeIter h_in_15 = newHalfedge();

          // new edges
          EdgeIter e5 = newEdge();
          EdgeIter e6 = newEdge();
          EdgeIter e7 = newEdge();

          // new faces
          FaceIter f3 = newFace();
          FaceIter f4 = newFace();

          // new vertex
          VertexIter m = newVertex();


          // *************** updates ***************
          // update halfedge connections
          h_in->setNeighbors(h_in_1, h_in_3, m, e0, f1);
          h_in_1->setNeighbors(h_in_2, h_out_6, c, e1, f1);
          h_in_2->setNeighbors(h_in, h_in_11, a, e5, f1);
          h_in_3->setNeighbors(h_in_4, h_in, c, e0, f2);
          h_in_4->setNeighbors(h_in_5, h_in_15, m, e7, f2);
          h_in_5->setNeighbors(h_in_3, h_out_9, d, e4, f2);
          h_out_6->setNeighbors(h_out_6->next(), h_in_1, a, e1, h_out_6->face());
          h_out_7->setNeighbors(h_out_7->next(), h_in_12, b, e2, h_out_7->face());
          h_out_8->setNeighbors(h_out_8->next(), h_in_14, d, e3, h_out_8->face());
          h_out_9->setNeighbors(h_out_9->next(), h_in_5, c, e4, h_out_9->face());
          h_in_10->setNeighbors(h_in_11, h_in_13, b, e6, f3);
          h_in_11->setNeighbors(h_in_12, h_in_2, m, e5, f3);
          h_in_12->setNeighbors(h_in_10, h_out_7, a, e2, f3);
          h_in_13->setNeighbors(h_in_14, h_in_10, m, e6, f4);
          h_in_14->setNeighbors(h_in_15, h_out_8, b, e3, f4);
          h_in_15->setNeighbors(h_in_13, h_in_4, d, e7, f4);

          // updates for new vertex
          m->isNew = 1;
          m->position = 0.5 * (b->position + c->position);

          // updates vertex halfedges
          m->halfedge() = h_in;
          b->halfedge() = h_in_10;
          c->halfedge() = h_in_1;
          a->halfedge() = h_in_12;
          d->halfedge() = h_in_5;


          // updates edges
          e0->halfedge() = h_in;
          e1->halfedge() = h_in_1;
          e2->halfedge() = h_in_12;
          e3->halfedge() = h_in_14;
          e4->halfedge() = h_in_5;
          e5->halfedge() = h_in_2;
          e6->halfedge() = h_in_10;
          e7->halfedge() = h_in_4;

          // updates edges status that touches the vertices
          e0->isNew = 0;
          e5->isNew = 1;
          e6->isNew = 0;
          e7->isNew = 1;

          // add 2 new faces
          f1->halfedge() = h_in;
          f2->halfedge() = h_in_3;
          f3->halfedge() = h_in_10;
          f4->halfedge() = h_in_13;

          return m;

      }

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

      for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
          HalfedgeIter h0 = e->halfedge();
          HalfedgeIter h1 = h0->next();
          HalfedgeIter h2 = h1->next();
          HalfedgeIter h3 = h0->twin();
          HalfedgeIter h4 = h3->next();
          HalfedgeIter h5 = h4->next();
          VertexIter v0 = h0->vertex();
          VertexIter v1 = h3->vertex();
          VertexIter v2 = h2->vertex();
          VertexIter v3 = h5->vertex();

          e->newPosition =  (v0->position + v1->position) * 3.0 / 8.0 +  (v2->position + v3->position) * 1.0 / 8.0;
          e->isNew = 0;
      }

      // 2. Compute the updated vertex positions associated with edges, and store it in Edge::newPosition.
      for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
          HalfedgeIter h = v->halfedge();
          Vector3D og_pos_sum(0, 0, 0);
          do {
              og_pos_sum += h->twin()->vertex()->position;
              h = h->twin()->next();
          } while (h != v->halfedge());

          float n = (float) v->degree();
          float u = (n == 3.0) ? (3.0 / 16.0) : (3.0 / (8.0 * n));

          v->newPosition = (1.0 - n * u) * v->position + u * og_pos_sum;
          v->isNew = 0;

      }

      // 3. Split every edge in the mesh, in any order. For future reference, we're also going to store some
      // information about which subdivide edges come from splitting an edge in the original mesh, and which edges
      // are new, by setting the flat Edge::isNew. Note that in this loop, we only want to iterate over edges of
      // the original mesh---otherwise, we'll end up splitting edges that we just split (and the loop will never end!)
      for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
          VertexIter v1 = e->halfedge()->vertex();
          VertexIter v2 = e->halfedge()->twin()->vertex();

          if(!(v1->isNew || v2->isNew)) {
              VertexIter v = mesh.splitEdge(e);
              v->newPosition = e->newPosition;
          }
      }

      // 4. Flip any new edge that connects an old and new vertex.
      for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
          VertexIter v1 = e->halfedge()->vertex();
          VertexIter v2 = e->halfedge()->twin()->vertex();

          if(e->isNew && (v1->isNew == 1 || v2->isNew == 1)) {
              mesh.flipEdge(e);
          }
      }


      // 5. Copy the new vertex positions into final Vertex::position.
      for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
          v->isNew = 0;
          v->position = v->newPosition;
      }

  }
}
