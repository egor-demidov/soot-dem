#include "../src/pf_dem_initial.h"
#include "../src/pf_dem_add_liquid.h"
#include <iostream>
using namespace pf_dem;
void require(bool c,const char *m){if(!c)throw std::runtime_error(m);}
double norm(Vector a){return std::sqrt(dot(a,a));}
System fixture(bool fixed=false){
    Parameters p;p.dx=.1;p.mu=.01;p.sff=.001;p.sfs=.003;
    Sphere a;a.center={-.18,0,0};a.radius=.2;a.mass=.8;a.velocity={.01,.02,0};a.angular_momentum={.001,0,.002};
    Sphere b=a;b.center={.18,0,0};b.velocity={-.02,0,.01};b.fixed=fixed;if(fixed){b.velocity={};b.angular_momentum={};}
    Sphere c=a;c.center={0,.45,.05};c.radius=.15;c.velocity={0,-.01,0};
    std::vector<Sample>s={{{0,0,.24},{.01,-.02,.003},1.001,.001,false},
                           {{0,.24,.22},{-.01,.01,0},.999,.001,false}};
    std::vector<Site> sites;
    for(size_t j=0;j<3;++j)for(auto local:{Vector{-.13,0,.05},Vector{.13,0,.05}}){
        sites.push_back({s.size(),j,local});s.push_back({{}, {},1,.001,true});}
    ContactModel cm;cm.stiffness=25;cm.damping=.02;cm.friction=.3;cm.slip_speed=.1;
    return System(p,std::move(s),2,{a,b,c},sites,cm);
}
int main(){
    // An off-center level of a single radial kernel is a convex sphere.
    LiquidSurfaceField spherical_field;
    spherical_field.accumulate({.3,.4,0},.5,1.);
    require(std::abs(spherical_field.mean_curvature(.2)-10.)<1e-12,
            "Liquid mean curvature must be positive and equal to 1/R on a sphere");
    // Equilibration clamps rigid motion, retains fluid evolution and survives
    // a mid-hold restart. Restored initial momentum belongs to the support.
    auto eq=fixture(true),eq_split=fixture(true);
    const auto initial_bodies=eq.bodies;const auto initial_fluid=eq.state[0];
    const double hold=3e-5;
    eq.begin_equilibration(hold);eq_split.begin_equilibration(hold);
    for(const auto &b:eq.bodies)require(b.fixed&&b.velocity==Vector{}&&b.angular_momentum==Vector{},"Equilibration did not freeze rigid motion");
    require(norm(eq.momentum_error())<1e-14&&norm(eq.momentum_error(true))<1e-14,"Initial clamp impulse not accounted");
    eq.advance(1e-5);eq_split.advance(1e-5);
    eq_split.save("equilibration.tmp");auto eq_restored=System::load("equilibration.tmp");std::remove("equilibration.tmp");
    require(eq_restored.equilibration_active&&eq_restored.equilibration_end==hold,"Restart lost equilibration clock");
    for(auto *run:{&eq,&eq_restored}){
        for(size_t j=0;j<run->bodies.size();++j)require(run->bodies[j].center==initial_bodies[j].center&&run->bodies[j].orientation==initial_bodies[j].orientation,"Held particle pose changed");
        require(run->state[0].x!=initial_fluid.x&&run->state[0].v!=initial_fluid.v,"Fluid did not evolve during equilibration");
        bool cross_refused=false;try{run->advance(hold);}catch(const std::runtime_error&){cross_refused=true;}
        require(cross_refused,"Step across release boundary accepted");
        run->advance(hold-run->time);
        require(!run->equilibration_active&&run->release_states.empty(),"Particles not released at end of equilibration");
        for(size_t j=0;j<run->bodies.size();++j){const auto &a=run->bodies[j],&b=initial_bodies[j];
            require(a.center==b.center&&a.orientation==b.orientation&&a.velocity==b.velocity&&a.angular_momentum==b.angular_momentum&&a.fixed==b.fixed,"Release did not restore original motion/fixity");}
        require(norm(run->momentum_error())<1e-13&&norm(run->momentum_error(true))<1e-13,"Equilibration support impulse balance");
        run->advance(1e-5);
        require(run->bodies[0].center!=initial_bodies[0].center,"DEM did not resume after equilibration");
        require(run->bodies[1].center==initial_bodies[1].center&&run->bodies[1].fixed,"Permanent fixed body released");
    }
    for(size_t i=0;i<eq.state.size();++i)require(eq.state[i].x==eq_restored.state[i].x&&eq.state[i].rho==eq_restored.state[i].rho&&eq.state[i].v==eq_restored.state[i].v,"Equilibration restart differs from uninterrupted run");
    eq_restored.save("equilibration_released.tmp");auto after_release=System::load("equilibration_released.tmp");std::remove("equilibration_released.tmp");
    require(!after_release.equilibration_active&&after_release.equilibration_end==hold&&!after_release.bodies[0].fixed,"Post-release restart refroze particles");
    // Common-point friction must conserve total angular momentum with overlap.
    auto s=fixture();auto c=contacts(s.bodies,s.contact_model,s.parameters.dx);
    Vector f{},tau{};for(size_t j=0;j<s.bodies.size();++j){f=add(f,c.loads[j].force);
        tau=add(tau,add(cross(s.bodies[j].center,c.loads[j].force),c.loads[j].torque));}
    require(norm(f)<1e-15&&norm(tau)<1e-15,"Frictional contact force/torque reciprocity");
    require(c.power<=0&&c.energy>0,"Contact elasticity and dissipative power");
    auto original_mass=s.state[0].mass+s.state[1].mass;
    for(int k=0;k<30;++k)s.advance(1e-5);
    require(norm(s.momentum_error())<1e-13,"Three-body coupled linear conservation");
    require(norm(s.momentum_error(true))<1e-13,"Three-body coupled angular conservation");
    require(original_mass==s.state[0].mass+s.state[1].mass,"No liquid mass restoration/source");
    require(s.bodies[0].orientation!=Quaternion{1,0,0,0},"Quaternion evolves under spin");
    // Fixed particle load is assigned to support impulse, not lost.
    auto fixed=fixture(true);const auto x=fixed.bodies[1].center;
    for(int k=0;k<30;++k)fixed.advance(1e-5);
    require(fixed.bodies[1].center==x,"Fixed body moved");
    require(norm(fixed.supportP)>0&&norm(fixed.momentum_error())<1e-13&&norm(fixed.momentum_error(true))<1e-13,"Fixed support impulse ledger");
    // Actual N-body split trajectory, including friction/contact, must restart.
    auto continuous=fixture(),split=fixture();
    for(int k=0;k<12;++k)continuous.advance(1e-5);
    for(int k=0;k<5;++k)split.advance(1e-5);
    split.save("pf_dem_test_checkpoint.tmp");auto restored=System::load("pf_dem_test_checkpoint.tmp");
    std::remove("pf_dem_test_checkpoint.tmp");
    for(int k=5;k<12;++k)restored.advance(1e-5);
    require(continuous.time==restored.time&&continuous.step==restored.step,"Restart clock mismatch");
    for(size_t i=0;i<continuous.state.size();++i)
        require(continuous.state[i].x==restored.state[i].x&&continuous.state[i].v==restored.state[i].v&&continuous.state[i].rho==restored.state[i].rho,"Restart liquid/site mismatch");
    for(size_t j=0;j<continuous.bodies.size();++j)
        require(continuous.bodies[j].center==restored.bodies[j].center&&continuous.bodies[j].velocity==restored.bodies[j].velocity&&continuous.bodies[j].angular_momentum==restored.bodies[j].angular_momentum&&continuous.bodies[j].orientation==restored.bodies[j].orientation,"Restart DEM mismatch");
    // Pure external DEM force callback receives synchronized stage state and
    // is counted as an applied impulse rather than a conservation defect.
    auto callback=fixture();std::vector<double>times;
    callback.extra_loads=[&](const auto &b,double t){times.push_back(t);std::vector<Load>loads(b.size());loads[2].force={.1,.2,.3};loads[2].torque={.01,0,0};return loads;};
    callback.advance(1e-5);
    require(times.size()==2&&times[0]==0&&times[1]==.5e-5,"Callback did not receive both stages");
    require(norm(callback.momentum_error())<1e-13&&norm(callback.momentum_error(true))<1e-13,"External impulse accounting");
    bool refused=false;try{callback.save("unused.tmp");}catch(const std::runtime_error&){refused=true;}
    require(refused,"Unserialized external law must not silently checkpoint");
    // Interior overlap quadrature has nonzero partition weight, not a hole.
    Parameters p;p.dx=.1;p.sff=0;p.mu=0;
    Sphere a;a.radius=.2;Sphere b=a;b.center={.2,0,0};
    std::vector<Sample> st={{{.12,.1,.25},{},1,.001,false},{{.12,0,0},{},1,.001,true},{{.3,0,0},{},1,.001,true}};
    System buried(p,st,1,{a,b},{{1,0,{.12,0,0}},{2,1,{.1,0,0}}});
    auto got=buried.rates();auto weighted=buried.state;weighted[1].quadrature_weight=.5;
    auto expected=evaluate(weighted,{{0,1},{0,2}},p);
    require(got.fluid.force==expected.force,"Union-volume partition weight differs");
    // Coincident volume grids in an overlap sum to ONE copy, not zero/two.
    Sphere ca;ca.center={-.1,0,0};ca.radius=.3;ca.mass=1;Sphere cb=ca;cb.center={.1,0,0};
    std::vector<Sample> ov={{{0,.31,0},{.01,0,0},1.01,.001,false},{{},{},1,.001,true},{{},{},1,.001,true}};
    System union_test(p,ov,1,{ca,cb},{{1,0,{.1,0,0}},{2,1,{-.1,0,0}}});
    auto union_rate=union_test.rates();
    auto single=evaluate(union_test.state,{{0,1}},p);
    require(norm(subtract(union_rate.fluid.force[0],single.force[0]))<1e-18,"Overlap did not count volume exactly once");
    require(norm(union_rate.fluid.force[0])>0,"Overlap region deleted");
    Vector reaction{};for(const auto &l:union_rate.hydrodynamic)reaction=add(reaction,l.force);
    require(norm(add(reaction,union_rate.fluid.force[0]))<1e-18,"Weighted force reciprocity");
    Vector torque=cross(union_test.state[0].x,union_rate.fluid.force[0]);
    for(size_t j=0;j<2;++j)torque=add(torque,add(cross(union_test.bodies[j].center,union_rate.hydrodynamic[j].force),union_rate.hydrodynamic[j].torque));
    require(norm(torque)<1e-18,"Weighted angular reciprocity");
    // Every pair-rate term uses the same quadrature factor.
    Parameters wp=p;wp.mu=.02;wp.sfs=.003;
    auto probe=union_test.state;probe[0].v={.01,.02,0};
    auto full=evaluate(probe,{{0,1}},wp);probe[1].quadrature_weight=.5;
    auto half=evaluate(probe,{{0,1}},wp);
    for(size_t i=0;i<probe.size();++i){
        require(norm(subtract(half.force[i],scale(full.force[i],.5)))<1e-18,"Weighted pair force");
        require(std::abs(half.density[i]-.5*full.density[i])<1e-18,"Weighted continuity");
        require(std::abs(half.viscous_row[i]-.5*full.viscous_row[i])<1e-18,"Weighted viscous row");
        require(std::abs(half.conservative_virial[i]-.5*full.conservative_virial[i])<1e-18,"Weighted virial");}
    require(std::abs(half.pair_energy-.5*full.pair_energy)<1e-18,"Weighted pair energy");
    require(std::abs(half.viscous_power-.5*full.viscous_power)<1e-18,"Weighted dissipation");
    union_test.save("overlap_legacy.tmp");
    {std::ifstream f("overlap_legacy.tmp");std::ostringstream text;text<<f.rdbuf();auto data=text.str();data.replace(0,7,"PFDEM 2");std::istringstream lines(data);std::string line;std::ofstream out("overlap_legacy.tmp");size_t row=0;while(std::getline(lines,line)){if(row==4)line.erase(line.find_last_of(' '));if(line.rfind("EQUILIBRATION ",0)!=0)out<<line<<'\n';++row;}}
    bool rejected_old=false;try{auto old=System::load("overlap_legacy.tmp");}catch(const std::runtime_error&){rejected_old=true;}
    require(rejected_old,"Legacy overlap checkpoint accepted different physics");
    auto legacy=System::load("overlap_legacy.tmp",{},true);require(legacy.legacy_overlap_mask,"Explicit legacy continuation ignored");
    auto legacy_rate=legacy.rates();require(norm(legacy_rate.fluid.force[0])==0,"Legacy masking not preserved");
    legacy.save("overlap_legacy_saved.tmp");auto again=System::load("overlap_legacy_saved.tmp");
    require(again.legacy_overlap_mask,"Legacy choice lost on checkpoint round trip");
    std::remove("overlap_legacy.tmp");std::remove("overlap_legacy_saved.tmp");
    // Initial liquid is the liquid geometry minus the entire solid union.
    InitialGeometry g;g.bodies={a,b};g.droplets={{{0,0,0},.45}};g.lower={-.6,-.6,-.6};g.upper={.6,.6,.6};
    auto generated=initialize(p,g);
    for(size_t i=0;i<generated.nf;++i)require(!g.in_solid(generated.state[i].x),"Liquid initialized in solid");
    for(const auto &site:generated.sites)require(generated.state[site.particle].mass==p.rho0*std::pow(p.dx,3),"Boundary quadrature mass changed");
    InitialGeometry sessile;
    sessile.wall=FlatWall{.13,.8};
    sessile.droplets={{{.025,-.031,.13},.4}};
    p.sfs=.003;
    auto drop=initialize(p,sessile);
    require(drop.bodies.empty()&&drop.sites.empty()&&drop.state.size()>drop.nf,
            "Sessile calibration requires a static boundary without DEM bodies");
    for(size_t i=0;i<drop.state.size();++i)
        require((drop.state[i].x[2]<sessile.wall->z)==(i>=drop.nf),
                "Wall plane must separate solid and liquid samples");
    const auto wall_site=drop.state[drop.nf];
    drop.advance(1e-5);
    require(drop.state[drop.nf].x==wall_site.x&&drop.state[drop.nf].v==Vector{},
            "Static wall moved");
    require(norm(drop.supportP)>0&&norm(drop.momentum_error())<1e-13&&norm(drop.momentum_error(true))<1e-13,
            "Flat wall support impulse accounting");
    drop.save("pf_dem_wall.tmp");
    auto sessile_restart=System::load("pf_dem_wall.tmp");
    std::remove("pf_dem_wall.tmp");
    drop.advance(1e-5);sessile_restart.advance(1e-5);
    require(drop.state[0].x==sessile_restart.state[0].x&&drop.supportP==sessile_restart.supportP,
            "Flat wall restart differs");
    std::cout<<"N-body reciprocal coupling, friction torque, fixed supports, exact restart, callback, and solid-union initialization passed\n";
}
