#include <cxxopts.hpp>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/property_map.h>
#include <iostream>


struct FaceInfo2
{
	FaceInfo2() {}
	int nesting_level;
	bool in_domain() {
		return nesting_level % 2 == 1;
	}
};

typedef CGAL::Exact_predicates_inexact_constructions_kernel			K;
typedef CGAL::Triangulation_vertex_base_2<K>						Vb;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2, K>		Fbb;
typedef CGAL::Constrained_triangulation_face_base_2<K, Fbb>			Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>				TDS;
typedef CGAL::Exact_predicates_tag									Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag>	CDT;
typedef CDT::Point													Point;
typedef CGAL::Polygon_2<K>											Polygon_2;
typedef CDT::Face_handle											Face_handle;
typedef CDT::Vertex_handle											Vertex_handle;


void
mark_domains(CDT& ct,
	Face_handle start,
	int index,
	std::list<CDT::Edge>& border)
{
	if (start->info().nesting_level != -1) {
		return;
	}
	std::list<Face_handle> queue;
	queue.push_back(start);
	while (!queue.empty()) {
		Face_handle fh = queue.front();
		queue.pop_front();
		if (fh->info().nesting_level == -1) {
			fh->info().nesting_level = index;
			for (int i = 0; i < 3; i++) {
				CDT::Edge e(fh, i);
				Face_handle n = fh->neighbor(i);
				if (n->info().nesting_level == -1) {
					if (ct.is_constrained(e)) border.push_back(e);
					else queue.push_back(n);
				}
			}
		}
	}
}


//explore set of facets connected with non constrained edges,
//and attribute to each such set a nesting level.
//We start from facets incident to the infinite vertex, with a nesting
//level of 0. Then we recursively consider the non-explored facets incident
//to constrained edges bounding the former set and increase the nesting level by 1.
//Facets in the domain are those with an odd nesting level.
void
mark_domains(CDT& cdt)
{
	for (CDT::Face_handle f : cdt.all_face_handles()) {
		f->info().nesting_level = -1;
	}
	std::list<CDT::Edge> border;
	mark_domains(cdt, cdt.infinite_face(), 0, border);
	while (!border.empty()) {
		CDT::Edge e = border.front();
		border.pop_front();
		Face_handle n = e.first->neighbor(e.second);
		if (n->info().nesting_level == -1) {
			mark_domains(cdt, n, e.first->info().nesting_level + 1, border);
		}
	}
}


int main(int argc, char **argv)
{
	cxxopts::Options options(argv[0], " - example command line options");
	options.add_options()
		("verts", "A list of ints", cxxopts::value<std::vector<int>>())
		;
	auto cmd = options.parse(argc, argv);
	std::vector<int> verts = cmd["verts"].as<std::vector<int>>();

	CDT cdt;

	int i0 = 1;
	while (true) {
		int n = verts[i0-1];

		Polygon_2 poly;
		for (int i = i0; i < n*2 + i0; i = i + 2) {
			poly.push_back(K::Point_2(verts[i], verts[i + 1]));
		}
		cdt.insert_constraint(poly.vertices_begin(), poly.vertices_end(), true);

		i0 += (n*2 + 1);
		if (i0 >= verts.size()) break;
	}

	


	//construct two non-intersecting nested polygons
	//Polygon_2 polygon1;
	//polygon1.push_back(Point(0, 0));
	//polygon1.push_back(Point(10, 0));
	//polygon1.push_back(Point(10, 10));
	//polygon1.push_back(Point(0, 10));
	//Polygon_2 polygon2;
	//polygon2.push_back(Point(1, 1));
	//polygon2.push_back(Point(2, 1));
	//polygon2.push_back(Point(2, 2));
	//polygon2.push_back(Point(1, 2));
	//Polygon_2 polygon3;
	//polygon2.push_back(Point(4, 4));
	//polygon2.push_back(Point(5, 4));
	//polygon2.push_back(Point(5, 5));
	//polygon2.push_back(Point(4, 5));

	////Insert the polygons into a constrained triangulation
	//CDT cdt;
	//cdt.insert_constraint(polygon1.vertices_begin(), polygon1.vertices_end(), true);

	//cdt.insert_constraint(polygon2.vertices_begin(), polygon2.vertices_end(), true);
	//cdt.insert_constraint(polygon3.vertices_begin(), polygon3.vertices_end(), true);

	//Mark facets that are inside the domain bounded by the polygon
	mark_domains(cdt);
	int count = 0;
	for (Face_handle f : cdt.finite_face_handles())
	{
		if (f->info().in_domain()) {
			for (int i = 0; i < 3; i++) {
				Vertex_handle vh = f->vertex(i);
				Point p = vh->point();
				std::cout << p.x() << " " << p.y() << " ";
			}
			std::cout << ',';
			++count;
		}
	}
	return 0;
}