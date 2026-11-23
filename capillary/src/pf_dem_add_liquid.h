#pragma once
#include "pf_dem.h"

#include <cstdint>
#include <numeric>
#include <random>
#include <tuple>
#include <unordered_set>
#include <utility>

namespace pf_dem {

enum class AdditionMode {
    liquid_vapor,
    uniform_liquid_vapor,
    uniform_solid_vapor,
    curvature_based
};

inline AdditionMode parse_addition_mode(const std::string &value) {
    if(value=="liquid-vapor")return AdditionMode::liquid_vapor;
    if(value=="uniform-liquid-vapor")return AdditionMode::uniform_liquid_vapor;
    if(value=="uniform-solid-vapor")return AdditionMode::uniform_solid_vapor;
    if(value=="curvature-based")return AdditionMode::curvature_based;
    throw std::invalid_argument("Invalid liquid addition mode: "+value);
}

struct KelvinAddition {
    double saturation_ratio;
    double kelvin_length;
    std::uint64_t seed=0;
};

inline bool supports_multiple_layers(AdditionMode mode) {
    return mode==AdditionMode::liquid_vapor || mode==AdditionMode::curvature_based;
}

struct LiquidAddition {
    size_t layers=0;
    size_t particles=0;
    size_t eligible_sites=0;
    double volume_before=0.;
    double added_volume=0.;
    double added_mass=0.;
    double mass_scale=1.;
    size_t resolved_solid_vapor_sites=0;
    size_t resolved_liquid_vapor_sites=0;
    size_t eligible_solid_vapor_sites=0;
    size_t added_solid_vapor_particles=0;
};

// Derivatives use coordinates scaled by the kernel support.
struct LiquidSurfaceField {
    double value=0.;
    Vector gradient{};
    std::array<Vector,3> hessian{};

    void accumulate(const Vector &d,double q,double normalization) {
        const double t=1-q;
        value+=normalization*t*t*t*t*(1+4*q);
        const double isotropic=-20*normalization*t*t*t;
        gradient=add(gradient,scale(d,isotropic));
        for(size_t i=0;i<3;++i)hessian[i][i]+=isotropic;
        if(q==0)return;
        const auto unit=scale(d,1/q);
        const double radial=60*normalization*q*t*t;
        for(size_t i=0;i<3;++i)for(size_t j=0;j<3;++j)
            hessian[i][j]+=radial*unit[i]*unit[j];
    }

    double mean_curvature(double support) const {
        const double g=std::sqrt(dot(gradient,gradient));
        const auto n=scale(gradient,1/g);
        double tangent_trace=0.;
        for(size_t i=0;i<3;++i)tangent_trace+=hessian[i][i]-n[i]*dot(hessian[i],n);
        return -tangent_trace/(2*g*support);
    }
};

// Select one layer against the checkpoint and all previously planned layers.
inline LiquidAddition plan_liquid_layer(const System &system,double target,size_t max_samples,
                                         AdditionMode mode,const std::optional<KelvinAddition> &kelvin,
                                         std::mt19937_64 &rng,std::vector<Sample> &staged) {
    const size_t capacity=max_samples-system.state.size();
    if(staged.size()>=capacity)
        throw std::runtime_error("Liquid addition would exceed max_samples");
    const auto liquid_sample=[&](size_t i)->const Sample& {
        return i<system.nf?system.state[i]:staged[i-system.nf];
    };

    using pf_neighbors::Cell;
    using pf_neighbors::CellHash;
    using pf_neighbors::cell_of;
    const double dx=system.parameters.dx;
    const double H=system.parameters.support();
    const double min_spacing=.85*dx, reach=1.45*dx;
    const double solid_clearance=.5*dx, solid_shell=1.5*dx;
    const bool curvature=mode==AdditionMode::curvature_based;
    const bool use_liquid=mode!=AdditionMode::uniform_solid_vapor;
    const bool use_solid=!use_liquid || curvature;
    const bool spread=mode==AdditionMode::uniform_liquid_vapor || mode==AdditionMode::uniform_solid_vapor;
    const double log_nu=kelvin?std::log(kelvin->saturation_ratio):0.;
    constexpr double pi=3.14159265358979323846;
    const double kernel_factor=21/(2*pi*H*H*H);
    const double nominal_mass=system.parameters.rho0*dx*dx*dx;

    using Grid=std::unordered_map<Cell,std::vector<size_t>,CellHash>;
    Grid liquid_cells,kernel_cells,body_cells;
    std::unordered_set<Cell,CellHash> seeds;
    double body_cell_width=0.;
    for(const auto &body:system.bodies)
        body_cell_width=std::max(body_cell_width,body.radius+(use_solid?solid_shell:solid_clearance));
    for(size_t j=0;j<system.bodies.size();++j)
        body_cells[cell_of(system.bodies[j].center,body_cell_width)].push_back(j);

    // A search radius no larger than the cell width fits in these 27 cells.
    auto visit_neighbors=[](const Grid &grid,Cell c,auto &&visit){
        for(int i=-1;i<=1;++i)for(int j=-1;j<=1;++j)for(int k=-1;k<=1;++k){
            const auto found=grid.find({c.x+i,c.y+j,c.z+k});
            if(found==grid.end())continue;
            for(size_t id:found->second)visit(id);
        }
    };

    auto in_solid=[&](const Vector &x){
        bool inside=false;
        visit_neighbors(body_cells,cell_of(x,body_cell_width),[&](size_t j){
            const auto &body=system.bodies[j];
            const auto d=subtract(x,body.center);
            const double r=body.radius-1e-8*dx;
            if(r>0 && dot(d,d)<r*r)inside=true;
        });
        return inside;
    };

    auto surface_field=[&](const Vector &x){
        LiquidSurfaceField field;
        visit_neighbors(kernel_cells,cell_of(x,H),[&](size_t i){
            const auto &a=liquid_sample(i);
            const auto d=scale(subtract(x,a.x),1/H);
            const double q=std::sqrt(dot(d,d));
            if(q<1)field.accumulate(d,q,(a.mass/a.rho)*kernel_factor);
        });
        return field;
    };
    auto liquid_curvature=[&](const Vector &x)->std::optional<double>{
        auto point=x;
        // Project the vapor-side site to the nearby 0.5 isosurface.
        for(int iteration=0;iteration<12;++iteration){
            const auto field=surface_field(point);
            const double g=std::sqrt(dot(field.gradient,field.gradient));
            if(g<=1e-8)break;
            const double residual=.5-field.value;
            if(std::abs(residual)<=1e-5){
                if(!in_solid(point))return field.mean_curvature(H);
                break;
            }
            const double distance=std::clamp(H*residual/g,-.5*dx,.5*dx);
            point=add(point,scale(field.gradient,distance/g));
            const auto shift=subtract(point,x);
            if(dot(shift,shift)>reach*reach)break;
        }
        return std::nullopt;
    };

    LiquidAddition report;
    const size_t liquid_count=system.nf+staged.size();
    for(size_t i=0;i<liquid_count;++i){
        const auto &a=liquid_sample(i);
        if(staged.empty())report.volume_before+=a.mass/a.rho;
        const Cell c=cell_of(a.x,dx);
        liquid_cells[cell_of(a.x,reach)].push_back(i);
        kernel_cells[cell_of(a.x,H)].push_back(i);
        if(use_liquid)
            for(int ix=-1;ix<=1;++ix)for(int iy=-1;iy<=1;++iy)for(int iz=-1;iz<=1;++iz)
                seeds.insert({c.x+ix,c.y+iy,c.z+iz});
    }

    auto lattice_position=[dx](pf_neighbors::Index i){
        return (static_cast<double>(i)+.5)*dx;
    };
    if(use_solid){
        // Enumerate shell intervals along z to avoid scanning solid interiors.
        for(const auto &body:system.bodies){
            const double inner=body.radius+solid_clearance,outer=body.radius+solid_shell;
            const Vector extent{outer,outer,outer};
            const auto lo=cell_of(subtract(body.center,extent),dx);
            const auto hi=cell_of(add(body.center,extent),dx);
            for(auto i=lo.x;i<=hi.x;++i)for(auto j=lo.y;j<=hi.y;++j){
                const double x=lattice_position(i)-body.center[0];
                const double y=lattice_position(j)-body.center[1];
                const double radial2=x*x+y*y;
                if(radial2>outer*outer)continue;
                const double near_z=std::sqrt(std::max(0.,inner*inner-radial2));
                const double far_z=std::sqrt(outer*outer-radial2);
                auto insert_interval=[&](double a,double b){
                    const auto first=cell_of({0,0,body.center[2]+a},dx).z;
                    const auto last=cell_of({0,0,body.center[2]+b},dx).z;
                    for(auto k=first;k<=last;++k)seeds.insert({i,j,k});
                };
                insert_interval(-far_z,-near_z);
                insert_interval(near_z,far_z);
            }
        }
    }

    struct Candidate {
        Cell cell;
        Sample sample;
        double surface_score;
        bool solid_vapor=false;
    };
    std::vector<Candidate> candidates;
    auto add_candidate=[&](Candidate candidate,std::optional<double> kappa){
        if(curvature){
            if(!kappa || !std::isfinite(*kappa))return;
            if(candidate.solid_vapor)++report.resolved_solid_vapor_sites;
            else ++report.resolved_liquid_vapor_sites;
            const double exponent=(.5*kelvin->kelvin_length)*(*kappa);
            if(exponent>=log_nu)return;
            // Log of nu - exp(exponent), including close-to-zero weights.
            candidate.surface_score=log_nu+std::log(-std::expm1(exponent-log_nu));
            if(candidate.solid_vapor)++report.eligible_solid_vapor_sites;
        }
        candidates.push_back(candidate);
    };
    for(const auto &c:seeds){
        const Vector x{lattice_position(c.x),lattice_position(c.y),lattice_position(c.z)};
        double nearest2=std::numeric_limits<double>::infinity();
        visit_neighbors(liquid_cells,cell_of(x,reach),[&](size_t i){
            const auto d=subtract(x,liquid_sample(i).x);
            nearest2=std::min(nearest2,dot(d,d));
        });
        const bool dry=nearest2>reach*reach;
        if(dry?!use_solid:(!use_liquid || nearest2<min_spacing*min_spacing))continue;
        bool clear=true;
        size_t nearest_body=system.bodies.size();
        double solid_gap=solid_shell;
        visit_neighbors(body_cells,cell_of(x,body_cell_width),[&](size_t j){
            const auto &body=system.bodies[j];
            const auto d=subtract(x,body.center);
            const double exclusion=body.radius+solid_clearance;
            const double r2=dot(d,d);
            if(r2<exclusion*exclusion)clear=false;
            if(dry){
                const double gap=std::sqrt(r2)-body.radius;
                if(gap<solid_gap || (gap==solid_gap && j<nearest_body)){
                    solid_gap=gap;
                    nearest_body=j;
                }
            }
        });
        if(!clear || (dry && nearest_body==system.bodies.size()))continue;

        double occupancy=0.,density=0.;
        Vector velocity{};
        visit_neighbors(kernel_cells,cell_of(x,H),[&](size_t i){
            const auto &a=liquid_sample(i);
            const auto d=subtract(x,a.x);
            const double q=std::sqrt(dot(d,d))/H;
            if(q>=1)return;
            const double t=1-q;
            const double weight=(a.mass/a.rho)*kernel_factor*t*t*t*t*(1+4*q);
            occupancy+=weight;
            density+=weight*a.rho;
            velocity=add(velocity,scale(a.v,weight));
        });
        if(dry && occupancy<=.1){
            const auto &body=system.bodies[nearest_body];
            const auto offset=subtract(x,body.center);
            if(curvature){
                const auto point=add(body.center,scale(offset,body.radius/(body.radius+solid_gap)));
                if(in_solid(point) || surface_field(point).value>.1)continue;
            }
            const Vector v=add(body.velocity,cross(body.omega(),offset));
            add_candidate({c,{x,v,system.parameters.rho0,nominal_mass,false},
                           -std::abs(solid_gap-dx)/dx,true},1/body.radius);
        } else if(!dry && occupancy>0 && occupancy<=.5){
            add_candidate({c,{x,scale(velocity,1/occupancy),density/occupancy,
                              nominal_mass,false},occupancy},
                          curvature?liquid_curvature(x):std::nullopt);
        }
    }
    if(curvature && candidates.empty()){
        if(report.resolved_solid_vapor_sites+report.resolved_liquid_vapor_sites==0)
            throw std::runtime_error("Curvature-based addition found no eligible surface sites with resolved curvature");
        throw std::runtime_error("Curvature-based addition has non-positive Kelvin weights at every eligible surface site");
    }

    // Rank by proximity to the target surface. Hash breaks score ties
    // without systematically coating one side of the box first.
    std::sort(candidates.begin(),candidates.end(),[](const Candidate &a,const Candidate &b){
        if(a.surface_score!=b.surface_score)return a.surface_score>b.surface_score;
        const auto ha=CellHash{}(a.cell),hb=CellHash{}(b.cell);
        if(ha!=hb)return ha<hb;
        return std::tie(a.cell.x,a.cell.y,a.cell.z)<std::tie(b.cell.x,b.cell.y,b.cell.z);
    });
    const size_t count_limit=std::min(candidates.size(),capacity-staged.size());
    report.eligible_sites=candidates.size();

    // Stratify by spatial median. Interleaving the two halves at each split
    // distributes any prefix across the whole candidate set in O(N log N).
    std::vector<size_t> order(candidates.size());
    std::iota(order.begin(),order.end(),size_t{0});
    if(curvature){
        // Independent exponential clocks give successive weighted draws
        // without replacement. Rank their logarithms to avoid overflow.
        for(auto &candidate:candidates){
            const double u=(static_cast<double>(rng()>>12)+.5)*0x1p-52;
            candidate.surface_score=std::log(-std::log(u))-candidate.surface_score;
        }
        std::sort(order.begin(),order.end(),[&](size_t a,size_t b){
            return std::tie(candidates[a].surface_score,a)<std::tie(candidates[b].surface_score,b);
        });
    }
    if(spread && order.size()>1){
        auto indices=order;
        auto split=[&](auto &&self,size_t first,size_t end,size_t output,size_t stride)->void {
            if(end-first==1){
                order[output]=indices[first];
                return;
            }
            auto lo=candidates[indices[first]].sample.x,hi=lo;
            for(size_t k=first+1;k<end;++k)for(int d=0;d<3;++d){
                const double x=candidates[indices[k]].sample.x[d];
                lo[d]=std::min(lo[d],x);
                hi[d]=std::max(hi[d],x);
            }
            const auto extent=subtract(hi,lo);
            const size_t axis=static_cast<size_t>(std::max_element(extent.begin(),extent.end())-extent.begin());
            const size_t mid=first+(end-first+1)/2;
            // Vary the split direction to avoid favoring one coordinate sign.
            const bool reverse=CellHash{}({static_cast<pf_neighbors::Index>(output),
                                          static_cast<pf_neighbors::Index>(stride),
                                          static_cast<pf_neighbors::Index>(axis)})%2;
            auto key=[&](size_t i){
                const double x=candidates[i].sample.x[axis];
                return std::pair<double,size_t>{reverse?-x:x,i};
            };
            std::nth_element(indices.begin()+first,indices.begin()+mid,indices.begin()+end,
                             [&](size_t a,size_t b){return key(a)<key(b);});
            self(self,first,mid,output,2*stride);
            self(self,mid,end,output+stride,2*stride);
        };
        split(split,0,indices.size(),0,1);
    }

    double sum=0.,nominal_volume=0.,best_scale_error=std::numeric_limits<double>::infinity();
    for(size_t i=0;i<count_limit;++i){
        sum+=nominal_mass/candidates[order[i]].sample.rho;
        const double error=std::abs(target/sum-1.);
        if(error<best_scale_error){
            nominal_volume=sum;
            report.particles=i+1;
            best_scale_error=error;
        }
        if(sum>=target)break;
    }

    // If this layer cannot hold the remainder, use its available sites at
    // nominal mass and rebuild candidates on the new surface.
    const bool continue_layers=supports_multiple_layers(mode) && count_limit>0
                               && sum<target && best_scale_error>.2
                               && count_limit<capacity-staged.size();
    if(continue_layers){
        report.particles=count_limit;
        nominal_volume=sum;
    } else if(report.particles==0 || best_scale_error>.2){
        std::ostringstream message;
        message<<"Cannot add requested volume in the available surface layer with added masses "
               <<"within 20% of nominal. Eligible sites: "<<count_limit
               <<", requested volume: "<<target<<", tested nominal volume: "<<sum
               <<". Change the dose or dx.";
        throw std::runtime_error(message.str());
    }
    report.mass_scale=continue_layers?1.:target/nominal_volume;
    for(size_t i=0;i<report.particles;++i){
        auto sample=candidates[order[i]].sample;
        sample.mass*=report.mass_scale;
        report.added_mass+=sample.mass;
        report.added_volume+=sample.mass/sample.rho;
        if(curvature && candidates[order[i]].solid_vapor)++report.added_solid_vapor_particles;
        staged.push_back(sample);
    }
    return report;
}

inline LiquidAddition add_surface_liquid(System &system,double target,size_t max_samples,
                                        AdditionMode mode=AdditionMode::liquid_vapor,
                                        std::optional<KelvinAddition> kelvin=std::nullopt) {
    if(!(target>0 && std::isfinite(target)) || system.bodies.empty() ||
       system.equilibration_active)
        throw std::invalid_argument("Liquid addition needs positive volume and an unfrozen DEM checkpoint");
    max_samples=std::min(max_samples,size_t{10000000}); // checkpoint reader limit
    if(system.state.size()>=max_samples)
        throw std::runtime_error("Liquid addition would exceed max_samples");
    const bool curvature=mode==AdditionMode::curvature_based;
    if(curvature!=kelvin.has_value() ||
       (kelvin && !(std::isfinite(kelvin->saturation_ratio) && kelvin->saturation_ratio>0 &&
                    std::isfinite(kelvin->kelvin_length) && kelvin->kelvin_length>=0)))
        throw std::invalid_argument("Curvature-based addition requires a positive saturation ratio and a non-negative Kelvin length");

    std::mt19937_64 rng(kelvin?kelvin->seed:0);
    std::vector<Sample> staged;
    LiquidAddition total;
    const double tolerance=1e-12*target;
    do {
        const double remaining=target-total.added_volume;
        LiquidAddition layer;
        try {
            layer=plan_liquid_layer(system,remaining,max_samples,mode,kelvin,rng,staged);
        } catch(const std::runtime_error &error) {
            if(total.layers==0)throw;
            std::ostringstream message;
            message<<(curvature?"Curvature-based":"Liquid-vapor")
                   <<" addition stopped at layer "<<total.layers+1
                   <<" with "<<remaining<<" m^3 unfilled: "<<error.what();
            throw std::runtime_error(message.str());
        }
        if(total.layers==0)total.volume_before=layer.volume_before;
        ++total.layers;
        total.particles+=layer.particles;
        total.eligible_sites+=layer.eligible_sites;
        total.added_volume+=layer.added_volume;
        total.added_mass+=layer.added_mass;
        total.mass_scale=layer.mass_scale;
        total.resolved_solid_vapor_sites+=layer.resolved_solid_vapor_sites;
        total.resolved_liquid_vapor_sites+=layer.resolved_liquid_vapor_sites;
        total.eligible_solid_vapor_sites+=layer.eligible_solid_vapor_sites;
        total.added_solid_vapor_particles+=layer.added_solid_vapor_particles;
    } while(supports_multiple_layers(mode) && target-total.added_volume>tolerance);

    // Commit once, so a blocked later layer cannot leave a partial dose.
    system.append_liquid(staged);
    return total;
}

} // namespace pf_dem
