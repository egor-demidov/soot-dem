#pragma once
#include "pf_dem_initial.h"
#include "pf_soot_settings.h"
#include <filesystem>
namespace pf_dem {
struct RunInput {
    Parameters p;InitialGeometry geometry;ContactModel contact;
    std::string signature;
    SootSettings soot;
    double end_time=.5,dt_scale=.75;
    long long max_steps=1000000,record_every=100,dump_every=500,checkpoint_every=500;
    int threads=1;
};
inline RunInput read_input(const std::string &path) {
    std::ifstream f(path);std::string magic;int version;RunInput a;size_t n;
    if(!(f>>magic>>version)||magic!="PFDEM_INPUT"||(version<1||version>3))throw std::runtime_error("Invalid prepared DEM input");
    f>>a.signature;
    if(a.signature.size()!=64||a.signature.find_first_not_of("0123456789abcdef")!=std::string::npos)throw std::runtime_error("Invalid input signature");
    f>>a.p.dx>>a.p.rho0>>a.p.sound>>a.p.mu>>a.p.sff>>a.p.sfs;
    f>>a.contact.stiffness>>a.contact.damping>>a.contact.friction>>a.contact.slip_speed;
    f>>a.end_time>>a.dt_scale>>a.max_steps>>a.record_every>>a.dump_every>>a.checkpoint_every>>a.threads;
    for(double &v:a.geometry.lower)f>>v;
    for(double &v:a.geometry.upper)f>>v;
    f>>a.geometry.max_samples;
    if(!(f>>n)||(version==3?n!=0:n==0)||n>100000)throw std::runtime_error("Invalid prepared DEM body count");
    a.geometry.bodies.resize(n);
    for(auto &b:a.geometry.bodies){for(double &v:b.center)f>>v;
        f>>b.radius>>b.mass;int fixed=-1;f>>fixed;
        if(fixed!=0&&fixed!=1)throw std::runtime_error("Invalid particle fixity");
        b.fixed=fixed;
        for(double &v:b.velocity)f>>v;
        Vector omega{};for(double &v:omega)f>>v;
        b.angular_momentum=scale(omega,b.inertia());}
    int mode=-1;f>>mode>>a.geometry.alpha>>a.geometry.threshold>>a.geometry.rm;
    if(mode!=0&&mode!=1)throw std::runtime_error("Invalid liquid mode");
    a.geometry.levelset=mode;
    if(!(f>>n)||n>100000)throw std::runtime_error("Invalid droplet count");
    a.geometry.droplets.resize(n);
    for(auto &d:a.geometry.droplets){for(double &v:d.center)f>>v;f>>d.radius;}
    if(version==2)a.soot.read(f,a.geometry.bodies.size());
    if(version==3) {
        int wall=-1;
        if(!(f>>wall) || (wall!=0 && wall!=1))throw std::runtime_error("Invalid calibration wall flag");
        if(wall) {
            a.geometry.wall.emplace();
            f>>a.geometry.wall->z>>a.geometry.wall->half_width;
        }
    }
    if(!f)throw std::runtime_error("Truncated prepared DEM input");
    std::string extra;
    if(f>>extra)throw std::runtime_error("Trailing prepared DEM input");
    a.p.validate();a.contact.validate();a.geometry.validate();
    if(!(a.end_time>0&&std::isfinite(a.end_time))||!(a.dt_scale>0&&a.dt_scale<=1)||
       a.max_steps<1||a.record_every<1||a.dump_every<0||a.checkpoint_every<0||a.threads<1)
        throw std::runtime_error("Invalid DEM runtime options");
    return a;
}
inline void write_vtp(const System &s,const std::string &name) {
    // Liquid point cloud. Density, velocity and pressure can be processed with
    // the same ParaView point-cloud reconstruction as the native SPH output.
    std::ofstream f(name);f<<std::setprecision(17);
    f<<"<?xml version=\"1.0\"?><VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\"><PolyData><Piece NumberOfPoints=\""<<s.nf<<"\" NumberOfVerts=\""<<s.nf<<"\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\"><Points><DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for(size_t i=0;i<s.nf;++i)for(double v:s.state[i].x)f<<v<<' ';
    f<<"</DataArray></Points><PointData>";
    for(const auto &field:{std::string("Density"),std::string("Pressure"),std::string("Mass")}){
        f<<"<DataArray type=\"Float64\" Name=\""<<field<<"\" format=\"ascii\">";
        for(size_t i=0;i<s.nf;++i){
            const auto &sample=s.state[i];
            double value=sample.mass;
            if(field=="Density")value=sample.rho;
            else if(field=="Pressure")value=pressure(sample.rho,s.parameters);
            f<<value<<' ';
        }
        f<<"</DataArray>";}
    f<<"<DataArray type=\"Float64\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"ascii\">";
    for(size_t i=0;i<s.nf;++i)for(double v:s.state[i].v)f<<v<<' ';
    f<<"</DataArray></PointData><Verts><DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">";
    for(size_t i=0;i<s.nf;++i)f<<i<<' ';
    f<<"</DataArray><DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">";
    for(size_t i=1;i<=s.nf;++i)f<<i<<' ';
    f<<"</DataArray></Verts></Piece></PolyData></VTKFile>\n";
    if(!f)throw std::runtime_error("Failed liquid VTP output");
}
inline void write_wall_vtp(const FlatWall &wall,const Droplet &drop,double dx) {
    // Draw the solid face over the full extent of the sampled plate.
    const double h=std::ceil(wall.half_width/dx)*dx;
    const auto &c=drop.center;
    std::ofstream f("wall.vtp");
    f<<std::setprecision(17)
      <<"<?xml version=\"1.0\"?><VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">"
      <<"<PolyData><Piece NumberOfPoints=\"4\" NumberOfVerts=\"0\" NumberOfLines=\"0\" "
      <<"NumberOfStrips=\"0\" NumberOfPolys=\"1\"><Points><DataArray type=\"Float64\" "
      <<"NumberOfComponents=\"3\" format=\"ascii\">";
    for(auto corner:{Vector{-h,-h,0},Vector{h,-h,0},Vector{h,h,0},Vector{-h,h,0}})
        f<<c[0]+corner[0]<<' '<<c[1]+corner[1]<<' '<<wall.z<<' ';
    f<<"</DataArray></Points><Polys><DataArray type=\"Int64\" Name=\"connectivity\" "
      <<"format=\"ascii\">0 1 2 3</DataArray><DataArray type=\"Int64\" Name=\"offsets\" "
      <<"format=\"ascii\">4</DataArray></Polys></Piece></PolyData></VTKFile>\n";
    if(!f)throw std::runtime_error("Failed flat wall VTP output");
}
inline void write_sphere_vtp(const System &s,const std::string &name) {
    constexpr int lat=12,lon=24;std::vector<Vector> local;
    local.push_back({0,0,1});
    constexpr double pi=3.14159265358979323846;
    for(int a=1;a<lat;++a)for(int b=0;b<lon;++b){double t=pi*a/lat,p=2*pi*b/lon;local.push_back({std::sin(t)*std::cos(p),std::sin(t)*std::sin(p),std::cos(t)});}
    size_t south=local.size();local.push_back({0,0,-1});
    std::vector<std::array<size_t,3>>tri;
    for(int j=0;j<lon;++j)tri.push_back({0,1+size_t(j),1+size_t((j+1)%lon)});
    for(int a=0;a<lat-2;++a)for(int j=0;j<lon;++j){size_t i=1+a*lon+j,k=1+a*lon+(j+1)%lon;
        tri.push_back({i,i+lon,k+lon});tri.push_back({i,k+lon,k});}
    for(int j=0;j<lon;++j)tri.push_back({south,south-lon+size_t((j+1)%lon),south-lon+size_t(j)});
    std::ofstream f(name);f<<std::setprecision(17);
    f<<"<?xml version=\"1.0\"?><VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\"><PolyData><Piece NumberOfPoints=\""<<local.size()*s.bodies.size()<<"\" NumberOfPolys=\""<<tri.size()*s.bodies.size()<<"\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\"><Points><DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">";
    for(const auto &b:s.bodies)for(auto p:local)for(double v:add(b.center,rotate(b.orientation,scale(p,b.radius))))f<<v<<' ';
    f<<"</DataArray></Points><PointData><DataArray type=\"Int32\" Name=\"ParticleId\" format=\"ascii\">";
    for(size_t j=0;j<s.bodies.size();++j)for(size_t i=0;i<local.size();++i)f<<j<<' ';
    f<<"</DataArray><DataArray type=\"Float64\" Name=\"MaterialZ\" format=\"ascii\">";
    for(size_t j=0;j<s.bodies.size();++j)for(auto p:local)f<<p[2]<<' ';
    f<<"</DataArray></PointData><Polys><DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">";
    for(size_t j=0;j<s.bodies.size();++j)for(auto t:tri)for(size_t v:t)f<<j*local.size()+v<<' ';
    f<<"</DataArray><DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">";
    for(size_t j=1;j<=tri.size()*s.bodies.size();++j)f<<3*j<<' ';
    f<<"</DataArray></Polys></Piece></PolyData></VTKFile>\n";
    if(!f)throw std::runtime_error("Failed particle VTP output");
}
// Diagnostic lines only: no radius, mass, quadrature sites or fluid exclusion.
inline void write_necks_vtp(const System &s,const std::string &name) {
    auto edges=s.mechanics?s.mechanics->necks():std::vector<std::pair<size_t,size_t>>{};
    std::ofstream f(name);f<<std::setprecision(17);
    f<<"<?xml version=\"1.0\"?><VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\"><PolyData><Piece NumberOfPoints=\""<<s.bodies.size()<<"\" NumberOfLines=\""<<edges.size()<<"\" NumberOfPolys=\"0\" NumberOfVerts=\"0\" NumberOfStrips=\"0\"><Points><DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">";
    for(const auto &b:s.bodies)for(double x:b.center)f<<x<<' ';
    f<<"</DataArray></Points><Lines><DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">";
    for(auto e:edges)f<<e.first<<' '<<e.second<<' ';
    f<<"</DataArray><DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">";
    for(size_t i=1;i<=edges.size();++i)f<<2*i<<' ';
    f<<"</DataArray></Lines></Piece></PolyData></VTKFile>\n";
    if(!f)throw std::runtime_error("Cannot write neck visualization");
}
} // namespace pf_dem
