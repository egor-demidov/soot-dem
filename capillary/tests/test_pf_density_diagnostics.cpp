#include "../src/pf_density_diagnostics.h"
#include <sstream>
#include <iostream>
using namespace pf_dem;
void require(bool x,const char*m){if(!x)throw std::runtime_error(m);}
int main(){
 Parameters p;p.dx=.1;p.rho0=1000;
 Sphere a;a.radius=.2;a.mass=1;
 std::vector<Sample> samples={{{.25,0,0},{1,2,3},800,.001,false},{{.3,0,0},{},1100,.001,false},{{},{},1000,.001,true}};
 System s(p,samples,2,{a},{{2,0,{.2,0,0}}});
 Stage rates;rates.fluid=pf_liquid::empty_rates(3);rates.fluid.density[0]=-1600;
 auto original=s.state;auto [lo,hi]=density_extreme_indices(s);require(lo==0&&hi==1,"Extreme indices");
 auto probe=density_probe(s,rates,lo);
 require(std::abs(probe.rho_ratio-.8)<1e-15,"Signed density ratio");
 require(probe.pressure<0,"Depletion has negative EOS pressure");
 require(probe.continuity_divergence==2,"Continuity divergence sign");
 require(std::abs(probe.surface_gap-.05)<1e-15&&probe.nearest_body==0,"Surface distance");
 require(probe.liquid_neighbors==1&&probe.boundary_neighbors==1,"Support neighbor counts");
 require(s.state[0].rho==original[0].rho&&s.state[0].v==original[0].v,"Probe is read-only");
 std::ostringstream out;write_density_header(out);write_density_probe(out,s,1e-15,"min",probe);
 require(out.str().find("rho_over_rho0")!=std::string::npos,"CSV header");
 require(out.str().find(",min,0,")!=std::string::npos,"CSV extreme label");
 std::cout<<"Density diagnostic tests passed\n";
}
