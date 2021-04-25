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

    const ForwardIterator x_begin_;
    const ForwardIterator x_end_;
    const ForwardIterator y_begin_;
    const ForwardIterator y_end_;

    const std::size_t node_count_;

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

}}}} // Namespaces
