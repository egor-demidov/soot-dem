#include "../src/viscous_stability.h"
#include <cassert>

// Unequal-mass pair: m1*b1=m2*b2 guarantees a reciprocal viscous force.
double pair_error(int refinement) {
    double u1=0., u2=1.;
    const double m1=1., m2=.001, conductance=.1;
    const double b1=conductance/m1, b2=conductance/m2;
    const double end=.1;
    const double bound=viscous_step_limit(std::max(b1,b2))/refinement;
    double time=0.;
    while (time<end) {
        const double dt=std::min(bound,end-time);
        const double old1=u1, old2=u2;
        const double energy=.5*(m1*u1*u1+m2*u2*u2);
        u1+=dt*b1*(old2-old1);
        u2+=dt*b2*(old1-old2);
        assert(u1>=0 && u1<=1 && u2>=0 && u2<=1);
        assert(std::abs(m1*u1+m2*u2-m2)<1e-14);
        assert(.5*(m1*u1*u1+m2*u2*u2)<=energy+1e-16);
        time+=dt;
    }
    return std::abs((u2-u1)-std::exp(-(b1+b2)*end));
}
int main() {
    assert(std::abs(harmonic_viscosity(.05,.00005)-.0001/1.001)<1e-15);
    const double inner=viscous_rate_numerator(.05,-3.,.1,.2,2.);
    const double wall=viscous_rate_numerator(.05,-3.,.1,.2,4.);
    assert(std::abs(wall-2*inner)<1e-15);
    assert(std::abs(viscous_row_rate(inner,.001)-1000*inner)<1e-12);
    assert(pair_error(4)<pair_error(2));
    assert(pair_error(2)<pair_error(1));
    assert(std::isinf(viscous_step_limit(0)));
    bool rejected=false;
    try { viscous_rate_numerator(.05,1.,.1,.2,2.); }
    catch (const std::runtime_error &) { rejected=true; }
    assert(rejected);
    // A frozen explicit force over .04 overshoots for an isolated rate=100.
    assert(1.-.04*100.<-1.);
    assert(viscous_step_limit(100)*100<=.5);
}
