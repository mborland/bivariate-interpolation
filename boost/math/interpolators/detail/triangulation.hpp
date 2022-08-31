//  (C) Copyright Robert J. Renka 1990.
//  (C) Copyright Matt Borland 2022.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  Derived from Algorithm 751: TRIPACK: a constrained two-dimensional Delaunay triangulation package
//  https://dl.acm.org/doi/pdf/10.1145/225545.225546

#ifndef BOOST_MATH_INTERPOLATORS_DETAIL_TRIANGULATION_HPP
#define BOOST_MATH_INTERPOLATORS_DETAIL_TRIANGULATION_HPP

#include <iterator>
#include <type_traits>
#include <boost/math/tools/promotion.hpp>

namespace boost { namespace math { namespace interpolators {

// Returns the signed area bounded by a polygonal curve, such as a constraint curve
template <typename RAIter, typename RAIter2>
auto polygonal_area (RAIter x_begin, RAIter x_end, RAIter y_begin, RAIter y_end, RAIter2 nodes_begin, RAIter2 nodes_end)
{
    using Real = typename boost::math::tools::promote_arg<typename std::iterator_traits<RAIter>::value_type>::type;

    const auto node_length = std::distance(nodes_begin, nodes_end);

    if (node_length <= 3)
    {
        return static_cast<Real>(0);
    }

    Real area = 0;
    auto node_2 = *nodes_end;
    while (nodes_begin != nodes_end)
    {
        auto node_1 = node_2;
        node_2 = nodes_begin;

        area += (*std::next(x_begin, node_2) - *std::next(x_begin, node_1)) *
                (*std::next(y_begin, node_1) + *std::next(y_begin, node_2));

        ++nodes_begin;
    }

    return -area / 2;
}

template <typename RAContainer, typename RAContainer2>
inline auto polygonal_area (RAContainer x, RAContainer y, RAContainer2 nodes)
{
    return polygonal_area(std::cbegin(x), std::cend(x),
                          std::cbegin(y), std::cend(y),
                          std::cbegin(nodes), std::cend(nodes));
}

}}} // Namespaces

#endif // BOOST_MATH_INTERPOLATORS_DETAIL_TRIANGULATION_HPP
