#include <iostream>
#include <vector>
#include<CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include<CGAL/Polygon_2.h>
#include<CGAL/create_offset_polygons_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polygon_2<K>   CGALPolygon;

int main()
{
  /*const double x[28] = { 3.4641,3.4641,3.4641,3.4641,3.4641,3.4641,3.4641,3.4641,3.4641,3.4641,3.4641,3.4641,3.4641,-6.52835e-17,-6.33865e-17,-5.78057e-17,-4.88654e-17,-3.70853e-17,-2.31499e-17,-7.86906e-18,7.86906e-18,2.31499e-17,3.70853e-17,4.88654e-17,5.78057e-17,6.33865e-17,6.52835e-17,3.4641, };
    const double y[28] = { 16.9915,15.4955,13.0989,9.94113,6.20559,2.10939,-2.10939,-6.20559,-9.94113,-13.0989,-15.4955,-16.9915,-17.5,-19.5,-18.9334,-17.2664,-14.596,-11.0773,-6.9148,-2.35047,2.35047,6.9148,11.0773,14.596,17.2664,18.9334,19.5,17.5, };*/
  double x[28] = { 3.4641015529632568,3.4641015529632568,3.4641015529632568,3.4641015529632568,3.4641015529632568,3.4641015529632568,3.4641015529632568,3.4641015529632568,3.4641015529632568,3.4641015529632568,3.4641015529632568,3.4641015529632568,3.4641015529632568,-6.5283512761760263e-17,-6.3386490615069335e-17,-5.7805677252904074e-17,-4.8865411228468343e-17,-3.7085263204542273e-17,-2.314985300804045e-17,-7.8690580132517397e-18,7.8690580132513576e-18,2.3149853008040067e-17,3.7085263204541891e-17,4.8865411228467961e-17,5.7805677252903704e-17,6.3386490615068965e-17,6.5283512761759894e-17,3.4641015529632568, };
  double y[28] = { 16.991481781005859,15.495480537414551,13.09893798828125,9.9411334991455078,6.2055854797363281,2.1093919277191162,-2.1093919277191162,-6.2055854797363281,-9.9411334991455078,-13.09893798828125,-15.495480537414551,-16.991481781005859,-17.5,-19.5,-18.933364868164063,-17.266391754150391,-14.595959663391113,-11.077262878417969,-6.9147953987121582,-2.3504652976989746,2.3504652976989746,6.9147953987121582,11.077262878417969,14.595959663391113,17.266391754150391,18.933364868164063,19.5,17.5, };

  CGALPolygon cgalPolygon;
  typedef boost::shared_ptr<CGALPolygon> PolygonPtr;
  for (int i = 0; i < 28; ++i) {
    auto point = K::Point_2(x[i], y[i]);
    cgalPolygon.push_back(K::Point_2(x[i], y[i]));
  }
  if (cgalPolygon.is_clockwise_oriented()) {
    cgalPolygon.reverse_orientation();
  }
  std::vector<PolygonPtr> exteriorSkeleton = CGAL::create_exterior_skeleton_and_offset_polygons_2(1e-5, cgalPolygon);
  std::cout << "Uncreachable point";

  return 0;
}
