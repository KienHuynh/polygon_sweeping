#include <cxxopts.hpp>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Triangular_expansion_visibility_2.h>
#include <CGAL/Simple_polygon_visibility_2.h>
#include <CGAL/Arrangement_2.h>
#include "arr_print.h"
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Quotient.h>
#include <istream>
#include <vector>
typedef CGAL::Exact_predicates_exact_constructions_kernel               Kernel;
typedef Kernel::Point_2                                                 Point_2;
typedef Kernel::Segment_2                                               Segment_2;
typedef CGAL::Arr_segment_traits_2<Kernel>                              Traits_2;
typedef CGAL::Arrangement_2<Traits_2>                                   Arrangement_2;
typedef Arrangement_2::Face_handle                                      Face_handle;
typedef Arrangement_2::Edge_const_iterator                              Edge_const_iterator;
typedef Arrangement_2::Halfedge_const_handle                    Halfedge_const_handle;
typedef Arrangement_2::Ccb_halfedge_circulator                          Ccb_halfedge_circulator;

// Define the used visibility class
typedef CGAL::Triangular_expansion_visibility_2<Arrangement_2>  TEV;

union myUnion {
	double dValue;
	uint64_t iValue;
};


void interior_vis(std::vector<Segment_2> segments, std::vector<double> q_) {
	Arrangement_2 env;
	CGAL::insert_non_intersecting_curves(env, segments.begin(), segments.end());
	Point_2 q(q_[0], q_[1]);
	// find the face of the query point
	// (usually you may know that by other means)
	Arrangement_2::Face_const_handle * face;
	CGAL::Arr_naive_point_location<Arrangement_2> pl(env);
	CGAL::Arr_point_location_result<Arrangement_2>::Type obj = pl.locate(q);
	// The query point locates in the interior of a face
	face = boost::get<Arrangement_2::Face_const_handle>(&obj);
	
	// compute regularized visibility area
	// Define visibiliy object type that computes regularized visibility area
	typedef CGAL::Simple_polygon_visibility_2<Arrangement_2, CGAL::Tag_true> RSPV;
	Arrangement_2 regular_output;
	RSPV regular_visibility(env);
	regular_visibility.compute_visibility(q, *face, regular_output);

	std::vector<std::vector<double>> output;
	std::vector<std::vector<std::vector<double>>> edges;

	for (Edge_const_iterator eit = regular_output.edges_begin(); eit != regular_output.edges_end(); ++eit) {
		//std::cout << "[" << eit->source()->point() << " -> " << eit->target()->point() << "]" << std::endl;
		std::vector<std::vector<double>> edge;
		double source_x = CGAL::to_double(eit->source()->point().x());
		double source_y = CGAL::to_double(eit->source()->point().y());
		double target_x = CGAL::to_double(eit->target()->point().x());
		double target_y = CGAL::to_double(eit->target()->point().y());
		edges.push_back(std::vector<std::vector<double>> { {source_x, source_y}, { target_x, target_y }});
	}

	// The edges returned above are not ordered
	// The block below finds the ordering of the vertices
	output.push_back(edges[0][0]);
	output.push_back(edges[0][1]);
	edges.erase(edges.begin());
	while (edges.size() > 1) {
		for (int i = 0; i < edges.size(); i++) {
			auto e = edges[i];
			if (e[0][0] == output.back()[0] && e[0][1] == output.back()[1]) {
				output.push_back(e[1]);
				edges.erase(edges.begin() + i);
				break;
			}
			if (e[1][0] == output.back()[0] && e[1][1] == output.back()[1]) {
				output.push_back(e[0]);
				edges.erase(edges.begin() + i);
				break;
			}
		}
	}

	for (auto o : output) {
		myUnion num0, num1;
		num0.dValue = o[0];
		num1.dValue = o[1];
		std::cout << num0.iValue << " " << num1.iValue << " ";
	}
}


void boundary_vis(std::vector<Segment_2> segments, std::vector<double> q_, std::vector<double> pre_q_, bool print_hex) {
	Point_2 q(q_[0], q_[1]);
	Point_2 pre_q(pre_q_[0], pre_q_[1]);

	Arrangement_2 env;
	CGAL::insert_non_intersecting_curves(env, segments.begin(), segments.end());

	// compute regularized vis area using half edge  
	//Find the halfedge whose target is the query point.
	Point_2 query_point = q;
	Halfedge_const_handle he = env.halfedges_begin();
	while (he->source()->point() != pre_q || he->target()->point() != q)
		he++;

	//visibility query
	Arrangement_2 output_arr;
	TEV tev(env);
	Face_handle fh = tev.compute_visibility(query_point, he, output_arr);
	//print out the visibility region.
	Arrangement_2::Ccb_halfedge_circulator curr = fh->outer_ccb();
	//print_arrangement(env);
	//print_arrangement(output_arr);
	std::vector<Point_2> output;
	output.push_back(curr->target()->point());
	while (++curr != fh->outer_ccb()) {
		output.push_back(curr->target()->point());
	}

	if (print_hex) {
		std::cout << std::hex;
		for (auto o : output) {
			myUnion num0, num1;
			num0.dValue = CGAL::to_double(o.x());
			num1.dValue = CGAL::to_double(o.y());
			std::cout << num0.iValue << " " << num1.iValue << " ";
		}
	}
	else {
		for (auto o : output) {
			std::cout << o.x() << " " << o.y() << " ";
		}
	}
}



int main(int argc, char **argv) {
	cxxopts::Options options(argv[0], " - example command line options");
	options.add_options()
		("verts", "A list of floats", cxxopts::value<std::vector<double>>())
		("query", "Query point", cxxopts::value<std::vector<double>>())
		("pre_query", "Point on the polygon BEFORE the query point (CCW)", cxxopts::value<std::vector<double>>()->default_value("0,0"))
		("boundary", "Specify if computing visibility for point on boundary", cxxopts::value<bool>()->default_value("false"))
		("print_hex", "Specify if print the numbers in hex", cxxopts::value<bool>()->default_value("true"));
	auto cmd = options.parse(argc, argv);
	std::vector<double> verts = cmd["verts"].as<std::vector<double>>();
	std::vector<double> q_ = cmd["query"].as<std::vector<double>>();
	std::vector<double> pre_q_ = cmd["pre_query"].as<std::vector<double>>();
	bool boundary = cmd["boundary"].as<bool>();
	bool print_hex = cmd["print_hex"].as<bool>();

	//create environment
	std::vector<Point_2> points;
	for (int i = 0; i < verts.size(); i += 2) {
		points.push_back(Point_2(verts[i], verts[i + 1]));
	}

	std::vector<Segment_2> segments;
	int s = points.size();
	for (int i = 0; i < points.size() - 1; i++) {
		segments.push_back(Segment_2(points[i], points[i + 1]));
	}
	segments.push_back(Segment_2(points.back(), points.front()));

	if (boundary) {
		boundary_vis(segments, q_, pre_q_, print_hex);
	}
	else {
		interior_vis(segments, q_);
	}

	return 0;
}