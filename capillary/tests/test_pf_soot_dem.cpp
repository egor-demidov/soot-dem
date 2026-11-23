#include "../src/pf_dem.h"
#include "../src/pf_soot_dem.h"
#include <sstream>
#include <iostream>
using namespace pf_dem;
void require(bool ok,const char *message){if(!ok)throw std::runtime_error(message);}
double norm(pf_rigid::Vector x){return std::sqrt(pf_liquid::dot(x,x));}
System fixture(bool neck=true,bool fixed=false){
    Parameters p;p.dx=.1;p.sfs=.003;
    Sphere a;a.radius=.24;a.mass=1;a.center={-.25,0,0};a.velocity={.001,.02,0};a.angular_momentum={.0001,.0002,0};
    Sphere b=a;b.center={.25,0,0};b.velocity={-.001,-.01,.001};b.angular_momentum={0,0,.0001};
    b.fixed=fixed;if(fixed){b.velocity={};b.angular_momentum={};}
    Sphere c=a;c.center={0,.6,0};c.velocity={};c.angular_momentum={};
    std::vector<Sample>s={{{0,0,.25},{.001,.002,.003},1.001,.001,false},{{0,.25,.2},{},.999,.001,false}};
    std::vector<Site> sites;
    for(size_t j=0;j<3;++j)for(auto r:{Vector{-.19,0,0},Vector{.19,0,0}}){sites.push_back({s.size(),j,r});s.push_back({{}, {},1,.001,true});}
    ContactModel cm;cm.stiffness=cm.damping=cm.friction=0;
    System out(p,s,2,{a,b,c},sites,cm);
    SootSettings settings;settings.enabled=true;settings.neck_mode=2;if(neck)settings.edges={{0,1}};
    out.mechanics=std::make_shared<SootDemLaw>(settings,out.bodies);return out;
}
std::string history(const MechanicalLaw &law){std::ostringstream out;out<<std::setprecision(17);law.save(out);return out.str();}
std::vector<double> neck_force_row(const SootDemLaw &law,const std::vector<Sphere> &b){
    std::ostringstream out;out<<std::setprecision(17);law.write_neck_forces(out,b,12,.1);
    std::string row=out.str();std::replace(row.begin(),row.end(),',',' ');
    std::istringstream in(row);std::vector<double> values;
    for(double value;in>>value;)values.push_back(value);
    return values;
}
int main(){
    // Frozen equilibration must not age contact springs or break loaded necks.
    auto eq=fixture(true,true);auto history_before=history(*eq.mechanics);
    auto start_pose=eq.bodies[0];eq.begin_equilibration(3e-5);
    eq.advance(1e-5);
    require(history(*eq.mechanics)==history_before,"Equilibration advanced mechanical history");
    require(eq.bodies[0].center==start_pose.center&&eq.bodies[0].orientation==start_pose.orientation,"Equilibration moved soot particle");
    eq.save("pf_soot_equilibration.tmp");auto eq_resume=System::load("pf_soot_equilibration.tmp",SootDemLaw::restore);std::remove("pf_soot_equilibration.tmp");
    require(eq_resume.equilibration_active&&history(*eq_resume.mechanics)==history_before,"Equilibration restart changed mechanical history");
    eq_resume.advance(eq_resume.equilibration_end-eq_resume.time);
    require(!eq_resume.equilibration_active&&history(*eq_resume.mechanics)==history_before,"Release advanced mechanical history");
    require(eq_resume.bodies[0].velocity==start_pose.velocity&&eq_resume.bodies[1].fixed,"Soot release lost initial state/fixity");
    eq_resume.advance(1e-5);
    require(history(*eq_resume.mechanics)!=history_before,"Released soot mechanics stayed frozen");
    require(norm(eq_resume.momentum_error())<1e-13&&norm(eq_resume.momentum_error(true))<1e-13,"Soot equilibration support balance");
    auto s=fixture(),dry=fixture(false);auto r=s.rates(),d=dry.rates();
    require(s.mechanics->neck_count()==1,"Explicit neck graph missing");
    require(r.fluid.force==d.fluid.force&&r.fluid.density==d.fluid.density,"Neck acquired fluid volume or quadrature");
    require(r.mechanical[0].force[0]>1.9,"Tensile neck did not pull particles together");
    // Compare canonical and reciprocal loads to the actual native functor.
    auto law=std::dynamic_pointer_cast<SootDemLaw>(s.mechanics);
    SootDemLaw::Model direct(law->model());
    std::vector<Eigen::Vector3d>x,v,w,z(3,Eigen::Vector3d::Zero());
    for(const auto &b:s.bodies){x.emplace_back(b.center[0],b.center[1],b.center[2]);v.emplace_back(b.velocity[0],b.velocity[1],b.velocity[2]);auto q=b.omega();w.emplace_back(q[0],q[1],q[2]);}
    auto ij=direct(0,1,x,v,z,w,0),ji=direct(1,0,x,v,z,w,0);
    for(int k=0;k<3;++k){
        require(std::abs(r.mechanical[0].force[k]-ij.first[k]*s.bodies[0].mass)<1e-13,"Native force conversion");
        require(std::abs(r.mechanical[1].force[k]-ji.first[k]*s.bodies[1].mass)<1e-13,"Reciprocal force differs from native functor");
        require(std::abs(r.mechanical[1].torque[k]-ji.second[k]*s.bodies[1].inertia())<1e-13,"Reciprocal torque differs from native functor");
    }
    // Repeated force queries and abandoned midpoint trials must not age springs.
    auto before=history(*law);for(int i=0;i<4;++i)(void)s.rates();
    require(history(*law)==before,"Force query changed committed history");
    auto force_row=neck_force_row(*law,s.bodies);
    require(force_row.size()==8&&force_row[0]==12&&force_row[1]==.1&&force_row[2]==0&&force_row[3]==1,
            "Neck forces CSV identifiers or field count");
    require(std::abs(force_row[4]-2.)<1e-13&&force_row[5]==0&&force_row[6]==0&&force_row[7]==0,
            "Initial neck force CSV must match tensile spring");
    require(history(*law)==before,"Neck force output changed mechanical history");
    law->begin_step(s.bodies,0,.001);(void)law->evaluate(s.bodies,.0005);law->rollback();
    require(history(*law)==before,"Rejected trial changed springs/graph");
    for(int i=0;i<15;++i)s.advance(1e-5);
    require(history(*law)!=before,"Accepted steps did not advance springs");
    require(norm(s.momentum_error())<1e-13&&norm(s.momentum_error(true))<1e-13,"Neck/fluid momentum conservation");
    // Canonical sliding history advances dt once, not once per force query.
    auto slide=fixture();for(auto &b:slide.bodies){b.velocity={};b.angular_momentum={};}
    slide.bodies[0].velocity={0,.1,0};
    auto sl=std::dynamic_pointer_cast<SootDemLaw>(slide.mechanics);
    sl->begin_step(slide.bodies,0,.002);
    for(int k=0;k<5;++k)(void)sl->evaluate(slide.bodies,.001);
    sl->commit(slide.bodies,slide.bodies,.001,.002);
    auto xi=std::get<0>(sl->model().pf_sinter_model().pf_springs()[1]);
    require(std::abs(xi[1]-.0002)<1e-15,"Sliding spring advanced more than one accepted dt");
    force_row=neck_force_row(*sl,slide.bodies);
    require(std::abs(force_row[5]-20*xi.norm())<1e-15,"Sliding force CSV mismatch");
    const auto &rolling=std::get<1>(sl->model().pf_sinter_model().pf_springs()[1]);
    require(std::abs(force_row[6]-5*rolling.norm())<1e-15,"Rolling force CSV mismatch");
    // Failure after midpoint preparation rolls back trial history and motion.
    auto failed=fixture();auto initial=failed.rates();auto saved_history=history(*failed.mechanics);
    failed.evaluator=[](const auto &state,const auto&,const auto&){auto result=empty_rates(state.size());result.density[0]=-1e12;return result;};
    bool caught=false;try{failed.advance(.001,&initial);}catch(const std::runtime_error&){caught=true;}
    require(caught&&failed.time==0&&history(*failed.mechanics)==saved_history,"Rejected fluid step committed mechanical state");
    // Check torsional elastic memory: it must remain axial during force queries.
    auto spin=fixture();for(auto &b:spin.bodies){b.velocity={};b.angular_momentum={};}
    spin.bodies[0].angular_momentum={.2*spin.bodies[0].inertia(),0,0};
    auto spinlaw=std::dynamic_pointer_cast<SootDemLaw>(spin.mechanics);
    spinlaw->begin_step(spin.bodies,0,.01);spinlaw->commit(spin.bodies,spin.bodies,.005,.01);
    auto torque=spinlaw->evaluate(spin.bodies,.01);
    require(std::abs(torque[0].torque[0]+5*.24*.24*.2*.01)<1e-13,"Axial torsion spring lost its history");
    force_row=neck_force_row(*spinlaw,spin.bodies);
    require(std::abs(force_row[7]-5*.24*.2*.01)<1e-13,"Torsional force CSV mismatch");
    auto fixed=fixture(true,true);for(int i=0;i<10;++i)fixed.advance(1e-5);
    require(norm(fixed.momentum_error())<1e-13&&norm(fixed.momentum_error(true))<1e-13,"Fixed neck support balance");
    // Exact history + liquid + quaternion restart, not only endpoint import.
    auto whole=fixture(),split=fixture();
    for(int k=0;k<12;++k)whole.advance(1e-5);
    for(int k=0;k<5;++k)split.advance(1e-5);
    split.save("pf_soot_restart.tmp");auto restored=System::load("pf_soot_restart.tmp",SootDemLaw::restore);std::remove("pf_soot_restart.tmp");
    for(int k=5;k<12;++k)restored.advance(1e-5);
    for(size_t i=0;i<whole.state.size();++i)require(whole.state[i].x==restored.state[i].x&&whole.state[i].v==restored.state[i].v&&whole.state[i].rho==restored.state[i].rho,"Restart fluid state differs");
    for(size_t j=0;j<whole.bodies.size();++j)require(whole.bodies[j].center==restored.bodies[j].center&&whole.bodies[j].velocity==restored.bodies[j].velocity&&whole.bodies[j].orientation==restored.bodies[j].orientation&&whole.bodies[j].angular_momentum==restored.bodies[j].angular_momentum,"Restart particle state differs");
    require(history(*whole.mechanics)==history(*restored.mechanics),"Restart spring/graph histories differ");
    // Breaking occurs only when the accepted step commits; dry contact follows.
    auto breaking=fixture();SootSettings weak;weak.neck_mode=2;weak.edges={{0,1}};weak.v[26]=.01;
    breaking.mechanics=std::make_shared<SootDemLaw>(weak,breaking.bodies);
    (void)breaking.rates();require(breaking.mechanics->neck_count()==1,"Trial force query broke a neck");
    breaking.advance(1e-5);require(breaking.mechanics->neck_count()==0,"Accepted endpoint did not break tensile neck");
    breaking.save("pf_soot_broken.tmp");auto broken=System::load("pf_soot_broken.tmp",SootDemLaw::restore);std::remove("pf_soot_broken.tmp");
    require(broken.mechanics->neck_count()==0,"Restart healed a broken neck");
    std::ostringstream broken_rows;
    std::dynamic_pointer_cast<SootDemLaw>(broken.mechanics)->write_neck_forces(broken_rows,broken.bodies,broken.step,broken.time);
    require(broken_rows.str().empty(),"Broken neck appeared in force CSV");
    std::cout<<"Native soot-dem forces/torques, invisible necks, transactional history, axial torsion, breakage, conservation and exact restart passed\n";
}
