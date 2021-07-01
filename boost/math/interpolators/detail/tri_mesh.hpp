// Copyright Robert J. Renka 1990-1996
// Copyright Matt Borland, 2021
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Source: R. J. Renka (1996). Algorithm 751: TRIPACK: a constrained two-dimensional Delaunay triangulation package. 
//         ACM Transactions on Mathematical Software. 22, 1-8.

#include <cstddef>
#include <cmath>
#include <stdexcept>
#include <iterator>
#include <limits>
#include <vector>
#include <algorithm>

namespace boost { namespace math { namespace interpolators { namespace detail {

template <typename ForwardIterator, typename Real = ForwardIterator::value_type>
class tri_mesh
{
private:
    // Nodal indexes in counterclockwise order of the verticies of a triangle
    struct nodal_index
    {
        Real i1;
        Real i2;
        Real i3;
    };

    // Triangle 
    struct triangle
    {
        Real x1;
        Real y1;
        Real x2;
        Real y2;
        Real x3;
        Real y3;

        // Circumcenter
        Real xc;
        Real yc;

        Real circumradius;
        Real area;
        Real aspect_ratio;

        triangle(Real x1_, Real y1_, Real x2_, Real y2_, Real x3_, Real y3_) : x1 {x1_}, y1 {y1_}, 
                                                                               x2 {x2_}, y2 {y2_}, 
                                                                               x3 {x3_}, y3 {y3_} {}
    }

    const ForwardIterator x_begin_;
    const ForwardIterator x_end_;
    const ForwardIterator y_begin_;
    const ForwardIterator y_end_;

    std::vector<std::size_t> nodes_;
    std::size_t node_count_;
    std::size_t boundary_node_count_;
    std::size_t arc_count_;
    std::size_t triangle_count_;

    // Set of nodal indexes. In order to distinguish between interior and boundary nodes, the last neighbor of each
    // boundary node is represented by the negative of its index
    const std::vector<long> list;

    // Set of indexes to list in one-to-one correspondence with the elements of list
    const std::vector<std::size_t> lptr;

    // Set of indexes to adjacency lists
    const std::vector<std::size_t> lend;

    // Index to the first empty location in list and lptr
    const std::size_t lnew;

    // Determines whether node N0 is to the left or right of the line through N1-N2 as viewed by an observer at N1
    // facing N2
    inline bool left(Real x1, Real y1, Real x2, Real y2, Real x0, Real y0) const noexcept;

    // True iff C is forward of A->B
    // iff <A->B, A->C> >= 0
    inline bool forward(Real xa, Real ya, Real xb, Real yb, Real xc, Real yc) const noexcept;

    // Locates a point P relative to a triangulation. If P is contained in a triangle, the vertexes are returned.
    // Otherwise, the indexes of the right most and leftmost visible boundary nodes are returned
    nodal_index find_triangle(std::size_t index, Real x, Real y) const;

    // Given a set of triangulation nodes update with a new node at position K
    // Returns the location of the new node
    std::size_t add_node(std::size_t index, Real x, Real y, std::size_t search_index, std::size_t node_location);

    void build_nodes();

    void circum(triangle& t);

public:
    tri_mesh(ForwardIterator x_begin, ForwardIterator x_end, ForwardIterator y_begin, ForwardIterator y_end);
};

template <typename ForwardIterator, typename Real>
tri_mesh<ForwardIterator, Real>::tri_mesh(ForwardIterator x_begin, ForwardIterator x_end, ForwardIterator y_begin, ForwardIterator y_end)
{   
    // Validate inputs
    node_count_ = std::distance(x_begin, x_end);
    if(node_count_ != std::distance(y_begin, y_end))
    {
        throw std::logic_error("X and Y must be the same length.\n");
    }
    else if(node_count_ < 3)
    {
        throw std::domain_error("X and Y must have at least three nodes for meshing.\n");
    }

    x_begin_ = x_begin;
    x_end_ = x_end;
    y_begin_ = y_begin;
    y_end_ = y_end;

    // Calculate the tolerance for calculations (10*epsilon)
    const Real tol = std::numeric_limits<Real>::epsilon();
    
    //
    // Store the first triangle
    //
    list.resize(node_count_);
    lptr.resize(node_count_);
    lend.resize(node_count_);

    if(!left(*x_begin, *y_begin, *std::next(x_begin, 1), *std::next(y_begin, 1), *std::next(x_begin, 2), *std::next(y_begin, 2)))
    {
        // The initial triangle is 1, 3, 2
        list[0] = 3;
        lptr[0] = 2;
        list[1] = -2;
        lptr[1] = 1;
        lend[0] = 2;

        list[2] = 1;
        lptr[2] = 4;
        list[3] = -3;
        lptr[3] = 3;
        lend[1] = 4;

        list[4] = 2;
        lptr[4] = 6;
        list[5] = -1;
        lptr[5] = 5;
        lend[2] = 6;
    } 
    else if(!left(*std::next(x_begin, 1), *std::next(y_begin, 1), 
                  *x_begin, *y_begin, 
                  *std::next(x_begin, 2), *std::next(y_begin, 2)))
    {
        // The initial triangle is 1, 2, 3
        list[0] = 2;
        lptr[0] = 2;
        list[1] = -3;
        lptr[1] = 1;
        lend[0] = 2;

        list[2] = 3;
        lptr[2] = 4;
        list[3] = -1;
        lptr[3] = 3;
        lend[1] = 4;

        list[4] = 1;
        lptr[4] = 6;
        list[5] = -2;
        lptr[5] = 5;
        lend[2] = 6;
    }
    else
    {
        // The first three nodes are co-linear
        throw std::logic_error("The first three nodes must not be co-linear.\n");
    }

    // Initialize lnew and add the remaining nodes
}

template <typename ForwardIterator, typename Real>
inline bool tri_mesh<ForwardIterator, Real>::left(Real x1, Real y1, Real x2, Real y2, Real x0, Real y0) const noexcept
{
    // Components of vector N1->N2
    const Real dx1 {x2 - x1};
    const Real dy1 {y2 - y1};
    
    // Components of vector N1->N0
    const Real dx2 {x0 - x1};
    const Real dy2 {y0 - y1};

    return dx1*dy2 >= dx2*dy1;
}

template <typename ForwardIterator, typename Real>
inline bool tri_mesh<ForwardIterator, Real>::forward(Real xa, Real ya, Real xb, Real yb, Real xc, Real yc) const noexcept
{
    return ((xb-xa) * (xc-xa) + (yb-ya) * (yc-ya)) >= 0;
}

template <typename ForwardIterator, typename Real>
tri_mesh<ForwardIterator, Real>::nodal_index tri_mesh<ForwardIterator, Real>::find_triangle(std::size_t index, Real x, Real y) const
{
    const Real xp = x;
    const Real yp = y;

    // Set n1 = nf and nl to the first and last neighbors of n0.
    auto n0 = index;
    auto lp = lend[index];
    auto nl = list[lp];
    lp = lptr[lp];
    auto nf = list[lp];
    auto n1 = nf;
    auto np = nl;
    auto npp = nf;
}

template <typename ForwardIterator, typename Real>
std::size_t tri_mesh<ForwardIterator, Real>::add_node(std::size_t index, Real x, Real y, std::size_t search_index, std::size_t node_location)
{
    // Call triangle find
}

// Given a sequence of points in the plane, this function computes the signed area bounded by the closed polygonal
// curve which passes through the points in the specified order.
template <typename RandomAccessContainer, typename Real = RandomAccessContainer::value_type>
Real area(const RandomAccessContainer& x, const RandomAccessContainer& y, const RandomAccessContainer& nodes)
{
    std::size_t node_1 = 0;
    std::size_t node_2 = nodes.back();
    Real partial_area = Real(0);

    if(nodal_boundaries > 3)
    {
        for(std::size_t i = 0; i < nodes.size(); ++i)
        {
            node_1 = node_2;
            node_2 = nodes_[i];

            partial_area += (x[node_2] - x[node_1]) * (y[node_1] + y[node_2]);
        }
    }

    // A contains twice the negative signed area of the region
    return -partial_area/2;
}

template <typename ForwardIterator, typename Real>
void tri_mesh<ForwardIterator, Real>::build_nodes()
{
    std::size_t start_node = 0;
    std::size_t lp = lend[start_node];

    while(list[start_node] > 0)
    {
        ++start_node;
    }

    nodes_.push_back(start_node);
    std::size_t k = 1;
    std::size_t current_node = start_node;

    // Traverse the boundary in counter-clockwise order
    lp = lend[current_node];
    current_node = list[lp];

    while(start_node != current_node)
    {
        ++k;
        nodes_.push_back(current_node)
        lp = lend[current_node];
        current_node = list[lp];
    }

    boundary_node_count_ = k;
    triangle_count_ = 2*node_count_ - boundary_node_count_ - 2;
    arc_count_ = triangle_count_ + node_count_ - 1;
}

template <typename ForwardIterator, typename Real>
void tri_mesh<ForwardIterator, Real>::circum(triangle& t)
{
    using std::abs;
    using std::sqrt;
    using std::fpclassify;

    std::array<Real, 3> u {t.x3 - t.x2, t.x1 - t.x3, t.x2 - t.x1};
    std::array<Real, 3> v {t.y3 - t.y2, t.y1 - t.y3, t.y2 - t.y1};

    t.area = (u[0] * v[1] - u[1] * v[0]) / 2;

    if(fpclassify(t.area) == FP_ZERO)
    {
        t.aspect_ratio = 0;
        return;
    }

    std::array<Real, 3> squared_distance {t.x1 * t.x1 + t.y1 * t.y1,
                                          t.x2 * t.x2 + t.y2 * t.y2,
                                          t.x3 * t.x3 + t.y3 * t.y3};

    // Compute factors of t.xc and t.yc
    Real fx = 0;
    Real fy = 0;

    for(std::size_t i = 0; i < u.size(); ++i)
    {
        fx += squared_distance[i] * v[i];
        fy += squared_distance[i] * u[i];
    }

    t.xc = fx / (4 * t.area);
    t.yc = fy / (4 * t.area);

    t.circumradius = sqrt((t.xc - t.x1)*(t.xc - t.x1) + (t.yc - t.y1)*(t.yc - t.y1));

    // Compute the squared edge lengths and aspect ratio
    std::array<Real, 3> squared_vector {u[0] * u[0] + v[0] * v[0],
                                        u[1] * u[1] + v[1] * v[1],
                                        u[2] * u[2] + v[2] * v[2]};

    t.aspect_ratio = 2 * abs(t.area) / (sqrt(squared_vector[0]) + sqrt(squared_vector[1]) + 
                                        sqrt(squared_vector[2]) * t.circumradius);
}

}}}} // Namespaces
