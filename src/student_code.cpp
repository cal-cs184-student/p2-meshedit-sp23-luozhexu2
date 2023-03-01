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
    std::vector<Vector2D> BezierCurve::evaluateStep(std::vector<Vector2D> const& points)
    {
        // TODO Part 1.
        int number = points.size();
        std::vector<Vector2D> newpoints; // the number of points in the iterator
        if (number == 1) {
            return points;
        }
        else {
            for (int i = 0; i <= number - 2; i++) {
                Vector2D lerp = (1.0 - t) * points[i] + t * points[i + 1];
                newpoints.push_back(lerp);
            }
            return newpoints;
        }
    }

    /**
     * Evaluates one step of the de Casteljau's algorithm using the given points and
     * the scalar parameter t (function parameter).
     *
     * @param points    A vector of points in 3D
     * @param t         Scalar interpolation parameter
     * @return A vector containing intermediate points or the final interpolated vector
     */
    std::vector<Vector3D> BezierPatch::evaluateStep(std::vector<Vector3D> const& points, double t) const
    {
        // TODO Part 2.
        int number = points.size();
        std::vector<Vector3D> newpoints;
        if (number == 1) {
            return points;
        }
        else {
            for (int i = 0; i <= number - 2; i++) {
                Vector3D lerp = (1.0 - t) * points[i] + t * points[i + 1];
                newpoints.push_back(lerp);
            }
            return newpoints;
        }
    }

    /**
     * Fully evaluates de Casteljau's algorithm for a vector of points at scalar parameter t
     *
     * @param points    A vector of points in 3D
     * @param t         Scalar interpolation parameter
     * @return Final interpolated vector
     */
    Vector3D BezierPatch::evaluate1D(std::vector<Vector3D> const& points, double t) const
    {
        // TODO Part 2.
        if (points.size() == 1) {
            return points[0];
        }
        else {
            return evaluate1D(evaluateStep(points, t), t);
        }
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
        int rows = controlPoints.size();
        std::vector<Vector3D> curvepoints;
        // for each row of the points, get the curve point
        for (int i = 0; i < rows; i++) {
            curvepoints.push_back(evaluate1D(controlPoints[i], u));
        }
        // do BezierCurve evaluation again for the curve points
        return evaluate1D(curvepoints, v);

    }

    Vector3D Vertex::normal(void) const
    {
        // TODO Part 3.
        // Returns an approximate unit normal at this vertex, computed by
        // taking the area-weighted average of the normals of neighboring
        // triangles, then normalizing.

        Vector3D result(0, 0, 0);
        HalfedgeCIter h = halfedge(); // initialization

        do {
            Vector3D v0 = position;
            Vector3D v1 = h->next()->vertex()->position;
            Vector3D v2 = h->next()->next()->vertex()->position;
            // two edges
            Vector3D s1 = v1 - v0;
            Vector3D s2 = v2 - v0;
            // compute the normal vector using cross product
            Vector3D n = cross(v1 - v0, v2 - v0);
            // compute the area
            double area = (cross(s1, s2).norm()) / 2.0;
            // use area as the weight to compute the result
            result += area * n;
            h = h->twin()->next();

        } while (h != halfedge());

        return result.unit();

    }

    EdgeIter HalfedgeMesh::flipEdge(EdgeIter e0)
    {
        // TODO Part 4.
        // This method should flip the given edge and return an iterator to the flipped edge.

          // check if the edge reaches the boundary
        if (e0->isBoundary())
            return e0;

        // the figure referenced by the code are shown below:
        /* C        C
          /|\      / \
         A | D -> A---D
          \|/      \ /
           B        B  */

           // initialize the halfedges
        HalfedgeIter h0 = e0->halfedge(); // initially BC
        HalfedgeIter h1 = h0->twin(); // initially CB
        HalfedgeIter CA = h0->next();
        HalfedgeIter AB = CA->next();
        HalfedgeIter BD = h1->next();
        HalfedgeIter DC = BD->next();

        // initialize the vertices
        VertexIter B = h0->vertex();
        VertexIter C = CA->vertex();
        VertexIter A = AB->vertex();
        VertexIter D = DC->vertex();

        // initialize the faces
        FaceIter CAD = h0->face();
        FaceIter ABD = h1->face();

        // update the halfedges
        h0->setNeighbors(DC, h1, A, e0, CAD);
        h1->setNeighbors(AB, h0, D, e0, ABD);
        CA->setNeighbors(h0, CA->twin(), C, CA->edge(), CAD);
        AB->setNeighbors(BD, AB->twin(), A, AB->edge(), ABD);
        BD->setNeighbors(h1, BD->twin(), B, BD->edge(), ABD);
        DC->setNeighbors(CA, DC->twin(), D, DC->edge(), CAD);

        // update the vertices
        C->halfedge() = CA;
        A->halfedge() = AB;
        B->halfedge() = BD;
        D->halfedge() = DC;

        // update the faces
        CAD->halfedge() = h0;
        ABD->halfedge() = h1;

        return e0;
    }

    VertexIter HalfedgeMesh::splitEdge(EdgeIter e0)
    {
        // TODO Part 5.
        // This method should split the given edge and return an iterator to the newly inserted vertex.
        // The halfedge of this vertex should point along the edge that was split, rather than the new edges.

          // check if the edge reaches the boundary
        if (e0->isBoundary())
            return VertexIter();

        // the figure referenced by the code are shown below:
        /* C        C
          /|\      /|\
         A | D -> A-M-D
          \|/      \|/
           B        B  */

           // initialize the old halfedges
        HalfedgeIter h0 = e0->halfedge(); // initially BC
        HalfedgeIter h1 = h0->twin(); // initially CB
        HalfedgeIter CA = h0->next();
        HalfedgeIter AB = CA->next();
        HalfedgeIter BD = h1->next();
        HalfedgeIter DC = BD->next();

        // create new halfedges, h0->BM, h1->MB
        HalfedgeIter MC = newHalfedge();
        HalfedgeIter CM = newHalfedge();
        HalfedgeIter AM = newHalfedge();
        HalfedgeIter MA = newHalfedge();
        HalfedgeIter MD = newHalfedge();
        HalfedgeIter DM = newHalfedge();

        // initialize the vertices
        VertexIter B = h0->vertex();
        VertexIter C = CA->vertex();
        VertexIter A = AB->vertex();
        VertexIter D = DC->vertex();
        VertexIter M = newVertex();
        M->isNew = true;

        // initialize the edges
        EdgeIter C_A = CA->edge();
        EdgeIter A_B = AB->edge();
        EdgeIter B_D = BD->edge();
        EdgeIter D_C = DC->edge();
        EdgeIter M_C = newEdge();
        EdgeIter A_M = newEdge();
        EdgeIter M_D = newEdge();
        A_M->isNew = true;
        M_D->isNew = true;

        // initialize the faces
        FaceIter ABM = h0->face();
        FaceIter MBD = h1->face();
        FaceIter AMC = newFace();
        FaceIter CMD = newFace();

        // update the M's position
        M->position = 0.5 * (B->position + C->position);

        // update the halfedges
        h0->setNeighbors(MA, h1, B, e0, ABM);
        h1->setNeighbors(BD, h0, M, e0, MBD);
        MA->setNeighbors(AB, AM, M, A_M, ABM);
        AM->setNeighbors(MC, MA, A, A_M, ABM);
        MD->setNeighbors(DC, DM, M, M_D, CMD);
        DM->setNeighbors(h1, MD, D, M_D, MBD);
        MC->setNeighbors(CA, CM, M, M_C, AMC);
        CM->setNeighbors(MD, MC, C, M_C, CMD);

        AB->setNeighbors(h0, AB->twin(), A, A_B, ABM);
        BD->setNeighbors(DM, BD->twin(), B, B_D, MBD);
        DC->setNeighbors(CM, DC->twin(), D, D_C, CMD);
        CA->setNeighbors(AM, CA->twin(), C, C_A, AMC);

        // update the vertices
        B->halfedge() = BD;
        C->halfedge() = CA;
        A->halfedge() = AB;
        D->halfedge() = DC;
        M->halfedge() = MC;

        // update the edges
        e0->halfedge() = h0;
        C_A->halfedge() = CA;
        A_B->halfedge() = AB;
        B_D->halfedge() = BD;
        D_C->halfedge() = DC;
        A_M->halfedge() = AM;
        M_C->halfedge() = MC;
        M_D->halfedge() = MD;

        // update the faces
        ABM->halfedge() = h0;
        MBD->halfedge() = h1;
        CMD->halfedge() = CM;
        AMC->halfedge() = MC;
        return M;

    }



    void MeshResampler::upsample(HalfedgeMesh& mesh)
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

          // 1.Compute new positions for all the old vertices
        for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
            v->isNew = false; // set it to be false
            HalfedgeIter h = v->halfedge();
            Vector3D original_position = v->position;
            Vector3D original_neighbor_position_sum(0, 0, 0);
            int n = 0;
            float u;

            // find all the neighbors of the vertex
            do {
                n++;
                original_neighbor_position_sum += h->twin()->vertex()->position;
                h = h->twin()->next();
            } while (h != v->halfedge());

            // use the formula to compute the final result
            if (n == 3) {
                u = 3.0 / 16.0;
            }
            else {
                u = 3.0 / (8.0 * n);
            }

            v->newPosition = (1 - n * u) * original_position + u * original_neighbor_position_sum;

        }

        // 2.Compute the updated vertex positions associated with edges
        for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
            // Set isNew to false for all old edges
            e->isNew = false;

            /* C
              /|\
             A | D
              \|/
               B  */

            HalfedgeIter h0 = e->halfedge();
            HalfedgeIter h1 = h0->twin();
            HalfedgeIter CA = h0->next();
            HalfedgeIter AB = CA->next();
            HalfedgeIter BD = h1->next();
            HalfedgeIter DC = BD->next();

            // find the position of the four vertices
            Vector3D pB = h0->vertex()->position;
            Vector3D pC = CA->vertex()->position;
            Vector3D pA = AB->vertex()->position;
            Vector3D pD = DC->vertex()->position;

            // use the formula to compute the final result
            e->newPosition = 3.0 / 8.0 * (pB + pC) + 1.0 / 8.0 * (pA + pD);

        }

        //3.Split every edge in the mesh
        for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
            /* C
              /|\
             A | D
              \|/
               B  */
            HalfedgeIter h0 = e->halfedge();
            HalfedgeIter h1 = h0->twin();
            VertexIter B = h0->vertex();
            VertexIter C = h1->vertex();

            // check whether both vertices are new
            if (B->isNew == false && C->isNew == false) {
                VertexIter new_v = mesh.splitEdge(e);
                new_v->newPosition = e->newPosition;
            }
        }

        //4. Flip any new edge that connects an old and new vertex
        for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
            /* C
              /|\
             A | D
              \|/
               B  */
            HalfedgeIter h0 = e->halfedge();
            HalfedgeIter h1 = h0->twin();
            VertexIter B = h0->vertex();
            VertexIter C = h1->vertex();

            // check whether one vertex is new, the other is old
            if ((B->isNew ^ C->isNew) && e->isNew == true) {
                mesh.flipEdge(e);
            }
        }

        //5. Copy the new vertex positions into final Vertex::position
        for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
            v->position = v->newPosition;
        }
    }
}
