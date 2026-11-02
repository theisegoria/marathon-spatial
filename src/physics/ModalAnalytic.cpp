// Closed-form modal models: free bar, simply supported plate, block estimate.
#include "marathon_spatial/physics/Modal.hpp"
#include <algorithm>
#include <cmath>
#include <functional>
#include <numbers>
#include <stdexcept>
#include <string>

namespace mspatial::physics {
namespace {
constexpr double Pi=std::numbers::pi;
// Roots of cos(b) cosh(b) = 1, the free-free beam (Rao, Mechanical Vibrations, table 8.4).
constexpr std::array<double,8> FreeBeam{4.730040745,7.853204624,10.995607838,14.137165491,17.278759657,20.420352245,23.561944902,26.703537556};
// Free-free beam shape at x in [0, 1], unit mean square.
double freeBeamShape(double b,double x){
    const double s=(std::cosh(b)-std::cos(b))/(std::sinh(b)-std::sin(b));
    return std::cosh(b*x)+std::cos(b*x)-s*(std::sinh(b*x)+std::sin(b*x));
}
struct Pending {double omega;std::vector<float> gains;};
void needPositive(double v,const char* what){if(!std::isfinite(v)||!(v>0))throw std::invalid_argument(std::string("analytic modal model: ")+what+" must be finite and positive");}
float clean(double v){return float(v)+0.f;} // no negative zero: the model's encoding has one form per value
}

ModalBody modalBodyForBox(Vec3 size,float wall,const ModalElastic& material){
    std::array<float,3> e{size.x,size.y,size.z};std::sort(e.begin(),e.end(),std::greater<>());
    ModalBody b;b.material=material;b.size={e[0],e[1],e[2]};
    if(wall>0&&wall<e[2]){b.shape=ModalShape::Plate;b.size.z=wall;}
    else if(e[2]<.2f*e[1])b.shape=ModalShape::Plate;
    else if(e[1]<.2f*e[0])b.shape=ModalShape::Bar;
    else b.shape=ModalShape::Block;
    return b;
}
ModalModel analyticModalModel(const ModalBody& body){
    const auto& t=body.material;const double E=t.youngModulusPa,rho=t.densityKgPerM3,nu=t.poissonRatio,eta=t.lossFactor;
    needPositive(E,"Young modulus");needPositive(rho,"density");needPositive(body.size.x,"size x");needPositive(body.size.y,"size y");needPositive(body.size.z,"size z");
    if(!std::isfinite(nu)||nu<0||nu>.495)throw std::invalid_argument("analytic modal model: Poisson ratio outside [0, 0.495]");
    if(!std::isfinite(eta)||eta<0)throw std::invalid_argument("analytic modal model: loss factor must be finite and non-negative");
    needPositive(body.maximumFrequencyHz,"maximum frequency");
    if(!(body.minimumFrequencyHz>=1)||body.minimumFrequencyHz>=body.maximumFrequencyHz)throw std::invalid_argument("analytic modal model: band invalid");
    if(!body.maximumModes||body.maximumModes>ModalAnalysisSettings::MaximumModes||body.samplesPerAxis<2||body.samplesPerAxis>12)throw std::invalid_argument("analytic modal model: mode or sample count outside bounds");
    if(!std::isfinite(body.airDampingPerSecond)||body.airDampingPerSecond<0)throw std::invalid_argument("analytic modal model: air damping invalid");
    const double a=body.size.x,b=body.size.y,c=body.size.z;const double mass=rho*a*b*c;
    const std::uint32_t N=body.samplesPerAxis;
    ModalModel m;
    auto sample=[&](Vec3 p,Vec3 n){m.samples.push_back({Vec3{clean(p.x),clean(p.y),clean(p.z)},n,std::uint32_t(m.samples.size())});};
    std::vector<Pending> modes;
    const double lowest=2*Pi*body.minimumFrequencyHz,highest=2*Pi*body.maximumFrequencyHz;
    const auto admit=[&](double omega){return omega>=lowest&&omega<=highest;};
    const double unit=1/std::sqrt(mass);
    if(body.shape==ModalShape::Plate){
        for(std::uint32_t i=0;i<N;++i)for(std::uint32_t j=0;j<N;++j)
            sample({float((i+.5)/N*a-a/2),float((j+.5)/N*b-b/2),float(c/2)},{0,0,1});
        const double stiffness=std::sqrt(E*c*c/(12*(1-nu*nu)*rho));
        for(std::uint32_t p=1;p<=16;++p)for(std::uint32_t q=1;q<=16;++q){
            const double omega=Pi*Pi*(p*p/(a*a)+q*q/(b*b))*stiffness;
            if(!admit(omega))continue;
            Pending out{omega,{}};out.gains.reserve(N*N);
            for(std::uint32_t i=0;i<N;++i)for(std::uint32_t j=0;j<N;++j)
                out.gains.push_back(clean(2*unit*std::sin(p*Pi*(i+.5)/N)*std::sin(q*Pi*(j+.5)/N)));
            modes.push_back(std::move(out));
        }
    }else if(body.shape==ModalShape::Bar){
        const std::uint32_t along=2*N;
        for(std::uint32_t i=0;i<along;++i)sample({float((i+.5)/along*a-a/2),0,float(c/2)},{0,0,1}); // bending through the thickness
        for(std::uint32_t i=0;i<along;++i)sample({float((i+.5)/along*a-a/2),float(b/2),0},{0,1,0}); // bending through the width
        for(std::uint32_t plane=0;plane<2;++plane){
            const double depth=plane?b:c;
            for(const double root:FreeBeam){
                const double omega=root*root*depth/(a*a)*std::sqrt(E/(12*rho));
                if(!admit(omega))continue;
                Pending out{omega,std::vector<float>(2*along,0.f)};
                for(std::uint32_t i=0;i<along;++i)out.gains[plane*along+i]=clean(unit*freeBeamShape(root,(i+.5)/along));
                modes.push_back(std::move(out));
            }
        }
    }else{
        const std::array<Vec3,3> normals{Vec3{1,0,0},Vec3{0,1,0},Vec3{0,0,1}};
        for(std::uint32_t f=0;f<3;++f)for(std::uint32_t i=0;i<N;++i)for(std::uint32_t j=0;j<N;++j){
            const double u=(i+.5)/N,v=(j+.5)/N;
            sample(f==0?Vec3{float(a/2),float(u*b-b/2),float(v*c-c/2)}:f==1?Vec3{float(u*a-a/2),float(b/2),float(v*c-c/2)}:Vec3{float(u*a-a/2),float(v*b-b/2),float(c/2)},normals[f]);}
        const double speed=std::sqrt(E/(2*(1+nu)*rho)); // shear: the free modes of a solid body are shear-dominated
        for(std::uint32_t l=0;l<=6;++l)for(std::uint32_t p=0;p<=6;++p)for(std::uint32_t q=0;q<=6;++q){
            if(l+p+q==0)continue;
            const double omega=Pi*speed*std::sqrt(l*l/(a*a)+p*p/(b*b)+q*q/(c*c));
            if(!admit(omega))continue;
            const double scale=unit*std::pow(std::numbers::sqrt2,double((l>0)+(p>0)+(q>0)));
            Pending out{omega,{}};out.gains.reserve(m.samples.size());
            for(const auto& s:m.samples)
                out.gains.push_back(clean(scale*std::cos(l*Pi*(s.position.x/a+.5))*std::cos(p*Pi*(s.position.y/b+.5))*std::cos(q*Pi*(s.position.z/c+.5))));
            modes.push_back(std::move(out));
        }
    }
    std::stable_sort(modes.begin(),modes.end(),[](const Pending& l,const Pending& r){return l.omega<r.omega;});
    std::vector<Pending> kept;
    for(auto& mode:modes){
        const double xi=eta/2+body.airDampingPerSecond/(2*mode.omega);
        if(xi>=.95)continue;
        double energy=0;for(float g:mode.gains)energy+=double(g)*g;
        m.modes.push_back({mode.omega,mode.omega/(2*Pi),xi,std::sqrt(energy/double(mode.gains.size()))});
        kept.push_back(std::move(mode));
        if(kept.size()==body.maximumModes)break;
    }
    if(kept.empty())throw std::invalid_argument("analytic modal model: no mode rings inside the band");
    m.gains.reserve(kept.size()*m.samples.size());
    for(std::size_t s=0;s<m.samples.size();++s)for(const auto& mode:kept)m.gains.push_back(mode.gains[s]);
    m.material={rho,E,nu,body.airDampingPerSecond,0.};
    m.latticeSpacingMeters=std::min({a,b,c});m.minimumFrequencyHz=body.minimumFrequencyHz;m.maximumFrequencyHz=body.maximumFrequencyHz;
    m.latticeNodes=std::max<std::uint32_t>(4,std::uint32_t(m.samples.size()));m.latticeTetrahedra=5;m.latticeMassKg=mass;
    m.boundsMinimum={clean(-a/2),clean(-b/2),clean(-c/2)};m.boundsMaximum={clean(a/2),clean(b/2),clean(c/2)};
    m.identity=modalModelIdentity(m);
    return m;
}
}
