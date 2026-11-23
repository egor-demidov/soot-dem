#include "../src/pf_snapshot.h"
#include <iostream>
using namespace pf_snapshot;
void require(bool v,const char *s){if(!v) throw std::runtime_error(s);}
int main(){
    std::stringstream idx("step,time,file\n0,0,pf_snapshot_0.csv\n23377,7.0002326908245482,pf_snapshot_23377.csv\n");
    const auto m=find_metadata(idx,"pf_snapshot_23377.csv");
    require(m.step==23377&&m.time==7.0002326908245482,"Snapshot index time");
    std::vector<pf_liquid::Sample> a={{{0,0,0},{},1,.1,false},{{0,0,-.2},{},1,.1,true}};
    std::stringstream csv("x,y,z,rho,vx,vy,vz\n0.1,0.2,0.3,1.01,0.4,0.5,0.6\n");
    auto b=load(csv,a,1);
    require(b[0].rho==1.01 && b[0].v[2]==.6 && b[1].solid && b[1].x==a[1].x,"CSV fields and wall retention");
    require(std::abs(fluid_momentum(b,1,false)[0]-.04)<1e-15,"Momentum rebase");
    std::stringstream bad("x,y,z,rho,vx,vy,vz\n0,0,0,1,0,0,0\n0,0,0,1,0,0,0\n");
    bool rejected=false; try{(void)load(bad,a,1);}catch(const std::runtime_error&){rejected=true;}
    require(rejected,"Reject extra samples");
    std::cout<<"PF legacy snapshot metadata, state recovery and rejection passed\n";
}
