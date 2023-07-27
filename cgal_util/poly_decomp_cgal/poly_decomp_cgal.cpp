#include <cxxopts.hpp>
#include <boost/mpl/long.hpp>
#include <boost/mpl/list/aux_/tag.hpp>
#include <boost/mpl/aux_/config/msvc.hpp>
#include <boost/mpl/aux_/config/workaround.hpp>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/partition_2.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/property_map.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <cassert>
#include <list>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Partition_traits_2<K, CGAL::Pointer_property_map<K::Point_2>::type > Partition_traits_2;
typedef Partition_traits_2::Point_2                         Point_2;
typedef Partition_traits_2::Polygon_2                       Polygon_2;  // a polygon of indices
typedef std::list<Polygon_2>                                Polygon_list;
/*
	  v4     v2
	  | \   /|
	  |  \ / |
	  |  v3  |
	  |      |
	  v0-----v1
 */
int main(int argc, char **argv)
{
	cxxopts::Options options(argv[0], " - example command line options");
	options.add_options()
		("verts", "A list of ints", cxxopts::value<std::vector<double>>())
		;
	auto cmd = options.parse(argc, argv);
	std::vector<double> verts = cmd["verts"].as<std::vector<double>>();

	std::vector<K::Point_2> points;
	for (int i = 0; i < verts.size(); i=i+2) {
		points.push_back(K::Point_2(verts[i], verts[i + 1]));
	}

	Partition_traits_2 traits(CGAL::make_property_map(points));
	Polygon_2 polygon;
	for (int i = 0; i < points.size(); i++) {
		polygon.push_back(i);
	}
	
	Polygon_list partition_polys;
	CGAL::optimal_convex_partition_2(polygon.vertices_begin(),
		polygon.vertices_end(),
		std::back_inserter(partition_polys),
		traits);

	std::cout << partition_polys.size() << ' ';
	for (const Polygon_2& poly : partition_polys) {
		std::cout << poly.size() << ' ';
		for (Point_2 p : poly.container()) {
			std::cout << points[p] <<  ' ';
		}
	}

	assert(CGAL::partition_is_valid_2(polygon.vertices_begin(),
		polygon.vertices_end(),
		partition_polys.begin(),
		partition_polys.end(),
		traits));
	return 0;
}