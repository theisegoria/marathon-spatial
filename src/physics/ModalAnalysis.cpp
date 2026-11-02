// Cook-time modal analysis: voxel lattice -> 5 P1 tetrahedra per cell -> dense
// lumped-mass FEM -> Accelerate dsyevr -> audible elastic modes with surface
// gains and a canonical, SHA-256-named encoding.
#define ACCELERATE_NEW_LAPACK 1
#include <Accelerate/Accelerate.h>
#include <CommonCrypto/CommonDigest.h>
#include "marathon_spatial/physics/Modal.hpp"
#include <algorithm>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>

namespace mspatial::physics {
namespace {
using D3=std::array<double,3>;
D3 operator-(D3 a,D3 b){return {a[0]-b[0],a[1]-b[1],a[2]-b[2]};}
D3 operator+(D3 a,D3 b){return {a[0]+b[0],a[1]+b[1],a[2]+b[2]};}
D3 operator*(D3 a,double s){return {a[0]*s,a[1]*s,a[2]*s};}
double dot3(D3 a,D3 b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
D3 cross3(D3 a,D3 b){return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};}
double norm3(D3 a){return std::sqrt(dot3(a,a));}
bool finite(Vec3 v){return std::isfinite(v.x)&&std::isfinite(v.y)&&std::isfinite(v.z);}
std::array<unsigned char,32> sha256(std::span<const unsigned char> bytes){
    std::array<unsigned char,32> digest{};CC_SHA256_CTX context;CC_SHA256_Init(&context);
    for(std::size_t at=0;at<bytes.size();){ // CC_LONG is 32 bits
        const auto part=static_cast<CC_LONG>(std::min<std::size_t>(bytes.size()-at,0x40000000));
        CC_SHA256_Update(&context,bytes.data()+at,part);at+=part;
    }
    CC_SHA256_Final(digest.data(),&context);
    return digest;
}
D3 toD3(Vec3 v){return {double(v.x),double(v.y),double(v.z)};}
float canonical(float v){return v==0.f?0.f:v;}   // -0 -> +0: one encoding per value
double canonical(double v){return v==0.?0.:v;}
Vec3 canonical(Vec3 v){return {canonical(v.x),canonical(v.y),canonical(v.z)};}
[[noreturn]] void reject(const std::string& what){throw std::invalid_argument("modal analysis: "+what);}
std::string fmt(double v){char b[64];std::snprintf(b,sizeof b,"%.6g",v);return b;}

constexpr std::uint32_t MaximumMeshVertices=65536, MaximumMeshTriangles=65536;
constexpr std::uint64_t MaximumGridCells=1u<<20;           // inside-test grid incl. empty cells
constexpr double MaximumWindingEvaluations=1e8;            // grid cells x triangles
constexpr double RigidRelativeTolerance=1e-10;             // |lambda| <= tol * ||A||inf (roundoff ~1e-15 observed)
constexpr double RigidSubspaceTolerance=1e-6;

void checkMaterial(const ModalMaterial& m,const std::string& who){
    const auto bad=[&](const char* what){throw std::invalid_argument(who+": "+what);};
    if(!std::isfinite(m.densityKgPerM3)||!(m.densityKgPerM3>0)||m.densityKgPerM3>1e5)bad("density must be finite in (0, 1e5] kg/m^3");
    if(!std::isfinite(m.youngModulusPa)||!(m.youngModulusPa>0)||m.youngModulusPa>1e13)bad("Young modulus must be finite in (0, 1e13] Pa");
    if(!std::isfinite(m.poissonRatio)||m.poissonRatio<0||m.poissonRatio>.495)bad("Poisson ratio must be in [0, 0.495] (P1 tetrahedra lock towards 0.5)");
    if(!std::isfinite(m.rayleighAlphaPerSecond)||m.rayleighAlphaPerSecond<0||m.rayleighAlphaPerSecond>1e5)bad("Rayleigh alpha must be finite in [0, 1e5] 1/s");
    if(!std::isfinite(m.rayleighBetaSeconds)||m.rayleighBetaSeconds<0||m.rayleighBetaSeconds>1e-2)bad("Rayleigh beta must be finite in [0, 1e-2] s");
}
void checkBand(double lo,double hi,const std::string& who){
    if(!std::isfinite(lo)||!std::isfinite(hi)||lo<1||!(hi>lo)||hi>100000)
        throw std::invalid_argument(who+": frequency band must satisfy 1 <= minimum < maximum <= 100000 Hz");
}

struct MeshData{std::vector<D3> v;std::vector<std::array<std::uint32_t,3>> t;D3 lo{},hi{};double volume{},diagonal{};};
MeshData admitMesh(const ModalMesh& mesh){
    MeshData d;
    if(mesh.vertices.size()<4||mesh.vertices.size()>MaximumMeshVertices)reject("mesh needs 4.."+std::to_string(MaximumMeshVertices)+" vertices");
    if(mesh.indices.size()%3)reject("index count is not a multiple of 3");
    const auto triangles=mesh.indices.size()/3;
    if(triangles<4||triangles>MaximumMeshTriangles)reject("mesh needs 4.."+std::to_string(MaximumMeshTriangles)+" triangles");
    d.v.reserve(mesh.vertices.size());
    for(const auto& p:mesh.vertices){if(!finite(p))reject("non-finite vertex");d.v.push_back(toD3(p));}
    d.lo=d.hi=d.v[0];
    for(const auto& p:d.v)for(int i=0;i<3;++i){d.lo[i]=std::min(d.lo[i],p[i]);d.hi[i]=std::max(d.hi[i],p[i]);}
    d.diagonal=norm3(d.hi-d.lo);
    if(!(d.diagonal>0))reject("mesh has zero extent");
    std::vector<std::pair<std::uint32_t,std::uint32_t>> edges;edges.reserve(triangles*3);
    for(std::size_t t=0;t<triangles;++t){
        std::array<std::uint32_t,3> tri{mesh.indices[3*t],mesh.indices[3*t+1],mesh.indices[3*t+2]};
        for(auto i:tri)if(i>=d.v.size())reject("triangle "+std::to_string(t)+" index out of range");
        if(tri[0]==tri[1]||tri[1]==tri[2]||tri[0]==tri[2])reject("triangle "+std::to_string(t)+" repeats a vertex");
        const D3 n=cross3(d.v[tri[1]]-d.v[tri[0]],d.v[tri[2]]-d.v[tri[0]]);
        if(!(norm3(n)>1e-12*d.diagonal*d.diagonal))reject("triangle "+std::to_string(t)+" is degenerate");
        d.volume+=dot3(d.v[tri[0]],cross3(d.v[tri[1]],d.v[tri[2]]))/6;
        for(int e=0;e<3;++e)edges.push_back({tri[e],tri[(e+1)%3]});
        d.t.push_back(tri);
    }
    std::sort(edges.begin(),edges.end());
    for(std::size_t i=1;i<edges.size();++i)if(edges[i]==edges[i-1])reject("directed edge repeated: mesh is non-manifold or inconsistently wound");
    for(const auto& e:edges)if(!std::binary_search(edges.begin(),edges.end(),std::pair{e.second,e.first}))reject("open mesh: edge without an opposite half-edge");
    if(!(d.volume>0))reject("mesh encloses non-positive volume (inward winding)");
    return d;
}
std::array<unsigned char,32> meshDigest(const ModalMesh& mesh){
    std::vector<unsigned char> b{'M','M','S','H'};
    const auto u32=[&](std::uint32_t v){for(int i=0;i<4;++i)b.push_back(static_cast<unsigned char>(v>>(8*i)));};
    u32(static_cast<std::uint32_t>(mesh.vertices.size()));u32(static_cast<std::uint32_t>(mesh.indices.size()));
    for(const auto& p:mesh.vertices){u32(std::bit_cast<std::uint32_t>(p.x));u32(std::bit_cast<std::uint32_t>(p.y));u32(std::bit_cast<std::uint32_t>(p.z));}
    for(auto i:mesh.indices)u32(i);
    return sha256(b);
}
// Generalized winding number (Van Oosterom-Strackee solid angles): robust for
// query points on triangle edges or vertex rays, unlike ray parity.
double winding(const MeshData& m,D3 p){
    double sum=0;
    for(const auto& t:m.t){
        const D3 a=m.v[t[0]]-p,b=m.v[t[1]]-p,c=m.v[t[2]]-p;
        const double la=norm3(a),lb=norm3(b),lc=norm3(c);
        sum+=2*std::atan2(dot3(a,cross3(b,c)),la*lb*lc+dot3(a,b)*lc+dot3(a,c)*lb+dot3(b,c)*la);
    }
    return sum/(4*std::numbers::pi);
}
// Closest point on triangle (Ericson, Real-Time Collision Detection 5.1.5).
D3 closestOnTriangle(D3 p,D3 a,D3 b,D3 c){
    const D3 ab=b-a,ac=c-a,ap=p-a;const double d1=dot3(ab,ap),d2=dot3(ac,ap);
    if(d1<=0&&d2<=0)return a;
    const D3 bp=p-b;const double d3=dot3(ab,bp),d4=dot3(ac,bp);
    if(d3>=0&&d4<=d3)return b;
    const double vc=d1*d4-d3*d2;
    if(vc<=0&&d1>=0&&d3<=0)return a+ab*(d1/(d1-d3));
    const D3 cp=p-c;const double d5=dot3(ab,cp),d6=dot3(ac,cp);
    if(d6>=0&&d5<=d6)return c;
    const double vb=d5*d2-d1*d6;
    if(vb<=0&&d2>=0&&d6<=0)return a+ac*(d2/(d2-d6));
    const double va=d3*d6-d5*d4;
    if(va<=0&&(d4-d3)>=0&&(d5-d6)>=0)return b+(c-b)*((d4-d3)/((d4-d3)+(d5-d6)));
    const double den=1/(va+vb+vc);
    return a+ab*(vb*den)+ac*(vc*den);
}

// Five-tetrahedron split of a cube, corners v = x | y<<1 | z<<2. Alternating
// with cell parity makes every shared face diagonal agree (conforming mesh).
// Measured on a free-free bar it is far less stiff in bending than the 6-tet
// Kuhn split (6.5% versus 17% above Euler-Bernoulli at 3 cells per thickness).
constexpr std::array<std::array<std::array<std::uint8_t,4>,5>,2> CubeSplit{{
    {{{0,1,2,4},{1,2,3,7},{1,4,5,7},{2,4,6,7},{1,2,4,7}}},
    {{{0,1,3,5},{0,2,3,6},{0,4,5,6},{3,5,6,7},{0,3,5,6}}}}};
constexpr std::array<std::array<std::uint8_t,4>,6> FaceCorners{{{0,2,4,6},{1,3,5,7},{0,1,4,5},{2,3,6,7},{0,1,2,3},{4,5,6,7}}};

struct Lattice{
    std::array<std::uint32_t,3> n{};D3 origin{};double h{};
    std::vector<std::uint8_t> inside;               // per grid cell
    std::vector<std::uint32_t> cells;               // included grid cell ids, ascending
    std::vector<std::int32_t> nodeOfGridNode;
    std::vector<D3> nodes;
    std::vector<std::array<std::uint32_t,4>> tets;  // positively oriented
    std::vector<double> nodeMass;
    std::vector<std::uint8_t> surface;              // per node
    double volume{};
    [[nodiscard]] std::uint64_t cellId(std::uint32_t i,std::uint32_t j,std::uint32_t k) const {return (std::uint64_t(i)*n[1]+j)*n[2]+k;}
    [[nodiscard]] std::uint64_t gridNodeId(std::uint32_t i,std::uint32_t j,std::uint32_t k) const {return (std::uint64_t(i)*(n[1]+1)+j)*(n[2]+1)+k;}
};

class SingleThreadedLapack{ // per-thread Accelerate setting: one reduction order, one result
public:
    SingleThreadedLapack():previous_(BLASGetThreading()){BLASSetThreading(BLAS_THREADING_SINGLE_THREADED);}
    ~SingleThreadedLapack(){BLASSetThreading(previous_);}
    SingleThreadedLapack(const SingleThreadedLapack&)=delete;SingleThreadedLapack& operator=(const SingleThreadedLapack&)=delete;
private:
    BLAS_THREADING previous_;
};

// Mass-scaled stiffness A = M^-1/2 K M^-1/2 (column-major, full), so A y = lambda y
// with u = M^-1/2 y mass-normalized (u^T M u = 1). Lumped mass is diagonal, which
// makes the scaling exact and lets dsyevr compute only the lowest eigenpairs.
void assembleScaledStiffness(const Lattice& l,const ModalMaterial& m,std::vector<double>& a){
    const std::size_t n=l.nodes.size()*3;
    a.assign(n*n,0.);
    const double lambda=m.youngModulusPa*m.poissonRatio/((1+m.poissonRatio)*(1-2*m.poissonRatio));
    const double mu=m.youngModulusPa/(2*(1+m.poissonRatio));
    std::vector<double> s(l.nodes.size());
    for(std::size_t i=0;i<s.size();++i)s[i]=1/std::sqrt(l.nodeMass[i]);
    for(const auto& t:l.tets){
        const D3 e1=l.nodes[t[1]]-l.nodes[t[0]],e2=l.nodes[t[2]]-l.nodes[t[0]],e3=l.nodes[t[3]]-l.nodes[t[0]];
        const double det=dot3(e1,cross3(e2,e3)),volume=det/6;
        std::array<D3,4> g{};g[1]=cross3(e2,e3)*(1/det);g[2]=cross3(e3,e1)*(1/det);g[3]=cross3(e1,e2)*(1/det);
        g[0]=(g[1]+g[2]+g[3])*-1.;
        for(int p=0;p<4;++p)for(int q=0;q<4;++q){
            const double gg=dot3(g[p],g[q]),scale=volume*s[t[p]]*s[t[q]];
            for(int i=0;i<3;++i)for(int j=0;j<3;++j){
                const double k=lambda*g[p][i]*g[q][j]+mu*g[q][i]*g[p][j]+(i==j?mu*gg:0.);
                a[(std::size_t(t[p])*3+i)+(std::size_t(t[q])*3+j)*n]+=k*scale;
            }
        }
    }
}
struct Eigen{std::vector<double> values,vectors;std::uint32_t count{};};
Eigen solveLowest(std::vector<double>& a,std::uint32_t n,std::uint32_t upper){
    const char jobz='V',range='I',uplo='L';
    const __LAPACK_int N=static_cast<__LAPACK_int>(n),lda=N,il=1,iu=static_cast<__LAPACK_int>(upper),ldz=N;
    const double vl=0,vu=0,abstol=0;
    __LAPACK_int m=0,info=0,lwork=-1,liwork=-1,iworkQuery=0;double workQuery=0;
    Eigen e;e.values.assign(n,0.);e.vectors.assign(std::size_t(n)*upper,0.);
    std::vector<__LAPACK_int> isuppz(2*std::size_t(upper));
    SingleThreadedLapack single;
    dsyevr_(&jobz,&range,&uplo,&N,a.data(),&lda,&vl,&vu,&il,&iu,&abstol,&m,e.values.data(),e.vectors.data(),&ldz,isuppz.data(),&workQuery,&lwork,&iworkQuery,&liwork,&info);
    if(info!=0)reject("dsyevr workspace query failed (info "+std::to_string(info)+")");
    lwork=static_cast<__LAPACK_int>(workQuery);liwork=iworkQuery;
    std::vector<double> work(static_cast<std::size_t>(std::max<__LAPACK_int>(lwork,1)));
    std::vector<__LAPACK_int> iwork(static_cast<std::size_t>(std::max<__LAPACK_int>(liwork,1)));
    dsyevr_(&jobz,&range,&uplo,&N,a.data(),&lda,&vl,&vu,&il,&iu,&abstol,&m,e.values.data(),e.vectors.data(),&ldz,isuppz.data(),work.data(),&lwork,iwork.data(),&liwork,&info);
    if(info!=0)reject("dsyevr failed (info "+std::to_string(info)+")");
    if(m!=iu)reject("dsyevr returned "+std::to_string(m)+" of "+std::to_string(iu)+" eigenpairs");
    e.count=static_cast<std::uint32_t>(m);
    return e;
}
}

void validateModalMaterial(const ModalMaterial& material,std::string_view who){checkMaterial(material,std::string(who));}
ModalMaterial modalMaterialFromLossFactor(double densityKgPerM3,double youngModulusPa,double poissonRatio,double lossFactor){
    if(!std::isfinite(lossFactor)||lossFactor<0)throw std::invalid_argument("modal material: invalid loss factor");
    const double low=2*std::numbers::pi*200,high=2*std::numbers::pi*5000;
    ModalMaterial m{densityKgPerM3,youngModulusPa,poissonRatio,lossFactor*low*high/(low+high),lossFactor/(low+high)};
    checkMaterial(m,"modal material");
    return m;
}

ModalAnalysisResult analyzeModal(const ModalMesh& mesh,const ModalMaterial& material,const ModalAnalysisSettings& s){
    using Clock=std::chrono::steady_clock;
    const auto started=Clock::now();
    checkMaterial(material,"modal analysis");
    checkBand(s.minimumFrequencyHz,s.maximumFrequencyHz,"modal analysis");
    if(!std::isfinite(s.latticeSpacingMeters)||!(s.latticeSpacingMeters>0))reject("lattice spacing must be finite and positive");
    if(s.maximumModes<1||s.maximumModes>ModalAnalysisSettings::MaximumModes)reject("maximum modes must be in 1.."+std::to_string(ModalAnalysisSettings::MaximumModes));
    if(s.maximumSurfaceSamples<1||s.maximumSurfaceSamples>ModalAnalysisSettings::MaximumSurfaceSamples)reject("maximum surface samples must be in 1.."+std::to_string(ModalAnalysisSettings::MaximumSurfaceSamples));
    if(s.maximumDegreesOfFreedom<24||s.maximumDegreesOfFreedom>ModalAnalysisSettings::MaximumDegreesOfFreedom)reject("maximum degrees of freedom must be in 24.."+std::to_string(ModalAnalysisSettings::MaximumDegreesOfFreedom));
    const MeshData md=admitMesh(mesh);
    ModalAnalysisReport report;report.meshVolume=md.volume;

    // Lattice: cubic cells centred on the mesh bounds; a cell is solid when its
    // centre has winding number > 1/2.
    Lattice l;l.h=s.latticeSpacingMeters;
    for(int i=0;i<3;++i){
        const double cells=std::ceil((md.hi[i]-md.lo[i])/l.h-1e-3); // float vertices: a 0.1% overhang is not a new cell
        if(!(cells<=double(MaximumGridCells)))reject("lattice spacing is too fine for the mesh bounds; coarsen latticeSpacingMeters");
        l.n[i]=static_cast<std::uint32_t>(std::max(1.,cells));
        l.origin[i]=(md.lo[i]+md.hi[i])/2-l.n[i]*l.h/2;
    }
    report.gridCells=std::uint64_t(l.n[0])*l.n[1]*l.n[2];
    if(report.gridCells>MaximumGridCells)reject("lattice grid has "+std::to_string(report.gridCells)+" cells above "+std::to_string(MaximumGridCells)+"; coarsen latticeSpacingMeters");
    if(double(report.gridCells)*double(md.t.size())>MaximumWindingEvaluations)reject("inside test would need more than 1e8 solid-angle evaluations; simplify the mesh or coarsen latticeSpacingMeters");
    l.inside.assign(report.gridCells,0);
    for(std::uint32_t i=0;i<l.n[0];++i)for(std::uint32_t j=0;j<l.n[1];++j)for(std::uint32_t k=0;k<l.n[2];++k){
        const D3 c{l.origin[0]+(i+.5)*l.h,l.origin[1]+(j+.5)*l.h,l.origin[2]+(k+.5)*l.h};
        if(winding(md,c)>.5){l.inside[l.cellId(i,j,k)]=1;l.cells.push_back(static_cast<std::uint32_t>(l.cellId(i,j,k)));}
    }
    if(l.cells.empty())reject("no lattice cell centre lies inside the mesh; refine latticeSpacingMeters");
    const auto decode=[&](std::uint64_t id){return std::array<std::uint32_t,3>{static_cast<std::uint32_t>(id/(std::uint64_t(l.n[1])*l.n[2])),static_cast<std::uint32_t>(id/l.n[2]%l.n[1]),static_cast<std::uint32_t>(id%l.n[2])};};
    constexpr std::array<std::array<int,3>,6> Neighbour{{{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}}};
    const auto neighbourInside=[&](std::array<std::uint32_t,3> c,int f){
        std::array<std::int64_t,3> q{};
        for(int a=0;a<3;++a){q[a]=std::int64_t(c[a])+Neighbour[f][a];if(q[a]<0||q[a]>=std::int64_t(l.n[a]))return false;}
        return l.inside[l.cellId(std::uint32_t(q[0]),std::uint32_t(q[1]),std::uint32_t(q[2]))]!=0;
    };
    { // Face-connected: cells meeting only along an edge or a vertex would hinge.
        std::vector<std::uint8_t> seen(report.gridCells,0);std::uint32_t components=0;std::vector<std::uint64_t> queue;
        for(auto start:l.cells){
            if(seen[start])continue;
            ++components;queue.assign(1,start);seen[start]=1;
            while(!queue.empty()){
                const auto id=queue.back();queue.pop_back();const auto c=decode(id);
                for(int f=0;f<6;++f)if(neighbourInside(c,f)){
                    const auto nid=l.cellId(std::uint32_t(std::int64_t(c[0])+Neighbour[f][0]),std::uint32_t(std::int64_t(c[1])+Neighbour[f][1]),std::uint32_t(std::int64_t(c[2])+Neighbour[f][2]));
                    if(!seen[nid]){seen[nid]=1;queue.push_back(nid);}
                }
            }
        }
        if(components!=1)reject("lattice has "+std::to_string(components)+" face-connected components (disconnected mesh or features thinner than the spacing); refine latticeSpacingMeters or analyze parts separately");
    }
    const std::uint64_t gridNodes=std::uint64_t(l.n[0]+1)*(l.n[1]+1)*(l.n[2]+1);
    l.nodeOfGridNode.assign(gridNodes,-1);
    const auto corner=[&](std::array<std::uint32_t,3> c,int v){return l.gridNodeId(c[0]+(v&1),c[1]+((v>>1)&1),c[2]+((v>>2)&1));};
    for(auto id:l.cells){const auto c=decode(id);for(int v=0;v<8;++v)l.nodeOfGridNode[corner(c,v)]=0;}
    std::uint32_t nodes=0;
    for(std::uint32_t i=0;i<=l.n[0];++i)for(std::uint32_t j=0;j<=l.n[1];++j)for(std::uint32_t k=0;k<=l.n[2];++k){
        auto& slot=l.nodeOfGridNode[l.gridNodeId(i,j,k)];
        if(slot<0)continue;
        slot=static_cast<std::int32_t>(nodes++);
        l.nodes.push_back({l.origin[0]+i*l.h,l.origin[1]+j*l.h,l.origin[2]+k*l.h});
    }
    const std::uint64_t dofs=std::uint64_t(nodes)*3;
    if(dofs>s.maximumDegreesOfFreedom)
        reject("lattice has "+std::to_string(dofs)+" degrees of freedom above the "+std::to_string(s.maximumDegreesOfFreedom)+
               " dense-solve bound; coarsen latticeSpacingMeters to about "+fmt(l.h*std::cbrt(double(dofs)/s.maximumDegreesOfFreedom)*1.05)+" m");
    report.latticeCells=static_cast<std::uint32_t>(l.cells.size());report.latticeNodes=nodes;report.degreesOfFreedom=static_cast<std::uint32_t>(dofs);
    l.nodeMass.assign(nodes,0.);l.surface.assign(nodes,0);
    for(auto id:l.cells){
        const auto c=decode(id);
        std::array<std::uint32_t,8> v{};for(int q=0;q<8;++q)v[q]=static_cast<std::uint32_t>(l.nodeOfGridNode[corner(c,q)]);
        for(const auto& t:CubeSplit[(c[0]+c[1]+c[2])&1]){
            std::array<std::uint32_t,4> tet{v[t[0]],v[t[1]],v[t[2]],v[t[3]]};
            const double det=dot3(l.nodes[tet[1]]-l.nodes[tet[0]],cross3(l.nodes[tet[2]]-l.nodes[tet[0]],l.nodes[tet[3]]-l.nodes[tet[0]]));
            if(det<0)std::swap(tet[1],tet[2]);
            const double volume=std::fabs(det)/6;
            l.volume+=volume;
            for(auto p:tet)l.nodeMass[p]+=material.densityKgPerM3*volume/4; // row-sum lumping (exact for P1)
            l.tets.push_back(tet);
        }
        for(int f=0;f<6;++f)if(!neighbourInside(c,f))for(auto q:FaceCorners[f])l.surface[v[q]]=1;
    }
    report.latticeTetrahedra=static_cast<std::uint32_t>(l.tets.size());report.latticeVolume=l.volume;
    double totalMass=0;for(auto m:l.nodeMass)totalMass+=m;
    report.latticeMassKg=totalMass;

    // Analytic rigid subspace in mass-scaled coordinates (3 translations, 3
    // rotations about the centre of mass), orthonormalized twice.
    const std::uint32_t n=static_cast<std::uint32_t>(dofs);
    std::vector<std::vector<double>> rigid(6,std::vector<double>(n,0.));
    D3 centre{};for(std::uint32_t p=0;p<nodes;++p)centre=centre+l.nodes[p]*(l.nodeMass[p]/totalMass);
    for(std::uint32_t p=0;p<nodes;++p){
        const double root=std::sqrt(l.nodeMass[p]);const D3 r=l.nodes[p]-centre;
        for(int d=0;d<3;++d){
            rigid[d][3*p+d]=root;
            D3 axis{};axis[d]=1;const D3 rot=cross3(axis,r);
            for(int i=0;i<3;++i)rigid[3+d][3*p+i]=root*rot[i];
        }
    }
    for(int pass=0;pass<2;++pass)for(int a=0;a<6;++a){
        for(int b=0;b<a;++b){double d=0;for(std::uint32_t i=0;i<n;++i)d+=rigid[a][i]*rigid[b][i];for(std::uint32_t i=0;i<n;++i)rigid[a][i]-=d*rigid[b][i];}
        double len=0;for(auto x:rigid[a])len+=x*x;len=std::sqrt(len);
        for(auto& x:rigid[a])x/=len;
    }

    std::vector<double> a;
    const auto assembled=Clock::now();
    report.latticeSeconds=std::chrono::duration<double>(assembled-started).count();
    std::vector<ModalMode> modes;std::vector<std::uint32_t> modeColumn;Eigen eigen;
    std::uint32_t upper=std::min<std::uint32_t>(n,6+s.maximumModes);
    for(int pass=0;pass<2;++pass){
        assembleScaledStiffness(l,material,a);
        double normInf=0;
        for(std::uint32_t r=0;r<n;++r){double row=0;for(std::uint32_t c=0;c<n;++c)row+=std::fabs(a[r+std::size_t(c)*n]);normInf=std::max(normInf,row);}
        eigen=solveLowest(a,n,upper);
        ++report.solverPasses;report.eigenpairsComputed=eigen.count;
        report.rigidTolerance=RigidRelativeTolerance*normInf;
        if(eigen.count<7)reject("lattice has fewer than 7 eigenpairs; refine latticeSpacingMeters");
        for(int r=0;r<6;++r){
            report.rigidEigenvalues[r]=eigen.values[r];
            if(!(std::fabs(eigen.values[r])<=report.rigidTolerance))reject("rigid-mode check: eigenvalue "+std::to_string(r)+" = "+fmt(eigen.values[r])+" is not ~0 (tolerance "+fmt(report.rigidTolerance)+")");
        }
        report.firstElasticEigenvalue=eigen.values[6];
        if(!(eigen.values[6]>report.rigidTolerance))reject("rigid-mode check: more than six near-zero modes (a hinge or mechanism in the lattice); refine latticeSpacingMeters");
        report.rigidSubspaceResidual=0;
        for(int r=0;r<6;++r){
            const double* z=eigen.vectors.data()+std::size_t(r)*n;
            std::vector<double> residual(z,z+n);
            for(int b=0;b<6;++b){double d=0;for(std::uint32_t i=0;i<n;++i)d+=rigid[b][i]*z[i];for(std::uint32_t i=0;i<n;++i)residual[i]-=d*rigid[b][i];}
            double len=0;for(auto x:residual)len+=x*x;
            report.rigidSubspaceResidual=std::max(report.rigidSubspaceResidual,std::sqrt(len));
        }
        if(!(report.rigidSubspaceResidual<=RigidSubspaceTolerance))reject("rigid-mode check: the six zero modes leave the rigid-body subspace (residual "+fmt(report.rigidSubspaceResidual)+")");
        modes.clear();modeColumn.clear();report.modesBelowBand=0;report.overdampedModes=0;
        bool reachedTop=false;
        for(std::uint32_t r=6;r<eigen.count&&modes.size()<s.maximumModes;++r){
            const double omega=std::sqrt(eigen.values[r]),hz=omega/(2*std::numbers::pi);
            if(hz<s.minimumFrequencyHz){++report.modesBelowBand;continue;}
            if(hz>s.maximumFrequencyHz){reachedTop=true;break;}
            const double xi=(material.rayleighAlphaPerSecond/omega+material.rayleighBetaSeconds*omega)/2;
            if(!(xi<1)){++report.overdampedModes;continue;}
            modes.push_back({omega,hz,xi,0});modeColumn.push_back(r);
        }
        // The first pass computes 6 + maximumModes pairs; modes below the band
        // or overdamped ones displace audible ones, so one wider pass recovers them.
        const std::uint32_t wider=std::min<std::uint32_t>(n,6+s.maximumModes+report.modesBelowBand+report.overdampedModes);
        if(reachedTop||modes.size()>=s.maximumModes||wider<=upper)break;
        upper=wider;
    }
    a.clear();a.shrink_to_fit();
    report.solveSeconds=std::chrono::duration<double>(Clock::now()-assembled).count();
    if(modes.empty())reject("no undamped elastic mode inside "+fmt(s.minimumFrequencyHz)+".."+fmt(s.maximumFrequencyHz)+" Hz");

    // Surface samples: farthest-point subset of lattice surface nodes.
    std::vector<std::uint32_t> surface;
    for(std::uint32_t p=0;p<nodes;++p)if(l.surface[p])surface.push_back(p);
    report.surfaceNodes=static_cast<std::uint32_t>(surface.size());
    D3 surfaceCentre{};for(auto p:surface)surfaceCentre=surfaceCentre+l.nodes[p]*(1./double(surface.size()));
    const std::size_t sampleCount=std::min<std::size_t>(surface.size(),s.maximumSurfaceSamples);
    std::vector<std::uint32_t> chosen;chosen.reserve(sampleCount);
    std::vector<double> nearest(surface.size(),std::numeric_limits<double>::infinity());
    std::size_t pick=0;double best=-1;
    for(std::size_t i=0;i<surface.size();++i){const double d=norm3(l.nodes[surface[i]]-surfaceCentre);if(d>best){best=d;pick=i;}}
    while(chosen.size()<sampleCount){
        const D3 reference=l.nodes[surface[pick]];
        chosen.push_back(surface[pick]);
        best=-1;
        for(std::size_t i=0;i<surface.size();++i){
            nearest[i]=std::min(nearest[i],norm3(l.nodes[surface[i]]-reference));
            if(nearest[i]>best){best=nearest[i];pick=i;}
        }
        if(!(best>0))break;
    }

    ModalModel model;
    model.material={canonical(material.densityKgPerM3),canonical(material.youngModulusPa),canonical(material.poissonRatio),canonical(material.rayleighAlphaPerSecond),canonical(material.rayleighBetaSeconds)};
    model.latticeSpacingMeters=s.latticeSpacingMeters;model.minimumFrequencyHz=s.minimumFrequencyHz;model.maximumFrequencyHz=s.maximumFrequencyHz;
    model.latticeNodes=nodes;model.latticeTetrahedra=report.latticeTetrahedra;model.latticeMassKg=totalMass;
    model.boundsMinimum=canonical(Vec3{float(md.lo[0]),float(md.lo[1]),float(md.lo[2])});
    model.boundsMaximum=canonical(Vec3{float(md.hi[0]),float(md.hi[1]),float(md.hi[2])});
    model.sourceMeshSHA256=meshDigest(mesh);
    for(auto node:chosen){
        double dmin=std::numeric_limits<double>::infinity();std::size_t first=0;D3 closest{};
        std::vector<double> distance(md.t.size());
        for(std::size_t t=0;t<md.t.size();++t){
            const D3 q=closestOnTriangle(l.nodes[node],md.v[md.t[t][0]],md.v[md.t[t][1]],md.v[md.t[t][2]]);
            distance[t]=norm3(q-l.nodes[node]);
            if(distance[t]<dmin){dmin=distance[t];first=t;closest=q;}
        }
        D3 normal{};
        for(std::size_t t=0;t<md.t.size();++t)if(distance[t]<=dmin+1e-7*md.diagonal){
            const auto& tri=md.t[t];normal=normal+cross3(md.v[tri[1]]-md.v[tri[0]],md.v[tri[2]]-md.v[tri[0]]); // area-weighted
        }
        if(!(norm3(normal)>1e-12*md.diagonal*md.diagonal)){const auto& tri=md.t[first];normal=cross3(md.v[tri[1]]-md.v[tri[0]],md.v[tri[2]]-md.v[tri[0]]);}
        normal=normal*(1/norm3(normal));
        Vec3 nf{float(normal[0]),float(normal[1]),float(normal[2])};
        const float fl=std::sqrt(nf.x*nf.x+nf.y*nf.y+nf.z*nf.z);nf={nf.x/fl,nf.y/fl,nf.z/fl};
        model.samples.push_back({canonical(Vec3{float(closest[0]),float(closest[1]),float(closest[2])}),canonical(nf),node});
    }
    const std::size_t S=model.samples.size(),M=modes.size();
    model.gains.assign(S*M,0.f);
    for(std::size_t k=0;k<M;++k){
        const double* z=eigen.vectors.data()+std::size_t(modeColumn[k])*n;
        // Mass-normalized shape u = M^-1/2 y; sign fixed so the largest |u| is positive.
        std::size_t arg=0;double peak=-1;
        for(std::uint32_t i=0;i<n;++i){const double u=std::fabs(z[i])/std::sqrt(l.nodeMass[i/3]);if(u>peak){peak=u;arg=i;}}
        const double sign=z[arg]<0?-1.:1.;
        double sumSquares=0;
        for(std::size_t si=0;si<S;++si){
            const auto& sample=model.samples[si];const std::uint32_t p=sample.node;
            const double inv=sign/std::sqrt(l.nodeMass[p]);
            const double g=inv*(z[3*p]*double(sample.normal.x)+z[3*p+1]*double(sample.normal.y)+z[3*p+2]*double(sample.normal.z));
            sumSquares+=g*g;
            model.gains[si*M+k]=canonical(static_cast<float>(g));
        }
        modes[k].surfaceRms=canonical(std::sqrt(sumSquares/double(S)));
    }
    model.modes=std::move(modes);
    model.identity=modalModelIdentity(model);
    return {std::move(model),report};
}

std::uint32_t nearestModalSample(const ModalModel& model,Vec3 p){
    if(model.samples.empty())throw std::invalid_argument("modal model has no surface samples");
    if(!finite(p))throw std::invalid_argument("modal contact position is not finite");
    std::uint32_t best=0;double bestDistance=std::numeric_limits<double>::infinity();
    for(std::uint32_t i=0;i<model.samples.size();++i){
        const D3 d=toD3(model.samples[i].position)-toD3(p);const double dd=dot3(d,d);
        if(dd<bestDistance){bestDistance=dd;best=i;}
    }
    return best;
}

// ---- canonical encoding ---------------------------------------------------
namespace {
constexpr std::array<unsigned char,4> Magic{'M','M','D','L'};
struct Writer{
    std::vector<unsigned char> b;
    void u32(std::uint32_t v){for(int i=0;i<4;++i)b.push_back(static_cast<unsigned char>(v>>(8*i)));}
    void u64(std::uint64_t v){for(int i=0;i<8;++i)b.push_back(static_cast<unsigned char>(v>>(8*i)));}
    void f32(float v){u32(std::bit_cast<std::uint32_t>(v));}
    void f64(double v){u64(std::bit_cast<std::uint64_t>(v));}
    void vec(Vec3 v){f32(v.x);f32(v.y);f32(v.z);}
};
struct Reader{
    std::span<const unsigned char> b;std::size_t at{};
    void need(std::size_t n){if(b.size()-at<n)throw std::invalid_argument("modal model decode: truncated");}
    std::uint32_t u32(){need(4);std::uint32_t v=0;for(int i=0;i<4;++i)v|=std::uint32_t(b[at+i])<<(8*i);at+=4;return v;}
    std::uint64_t u64(){need(8);std::uint64_t v=0;for(int i=0;i<8;++i)v|=std::uint64_t(b[at+i])<<(8*i);at+=8;return v;}
    float f32(){return std::bit_cast<float>(u32());}
    double f64(){return std::bit_cast<double>(u64());}
    Vec3 vec(){const float x=f32(),y=f32(),z=f32();return {x,y,z};}
};
void writePayload(Writer& w,const ModalModel& m){
    for(auto c:Magic)w.b.push_back(c);
    w.u32(ModalModel::FormatVersion);w.u32(0);
    w.f64(m.material.densityKgPerM3);w.f64(m.material.youngModulusPa);w.f64(m.material.poissonRatio);w.f64(m.material.rayleighAlphaPerSecond);w.f64(m.material.rayleighBetaSeconds);
    w.f64(m.latticeSpacingMeters);w.f64(m.minimumFrequencyHz);w.f64(m.maximumFrequencyHz);
    w.u32(m.latticeNodes);w.u32(m.latticeTetrahedra);w.f64(m.latticeMassKg);
    w.vec(m.boundsMinimum);w.vec(m.boundsMaximum);
    for(auto c:m.sourceMeshSHA256)w.b.push_back(c);
    w.u32(static_cast<std::uint32_t>(m.modes.size()));w.u32(static_cast<std::uint32_t>(m.samples.size()));
    for(const auto& mode:m.modes){w.f64(mode.angularFrequency);w.f64(mode.dampingRatio);w.f64(mode.surfaceRms);}
    for(const auto& sample:m.samples){w.vec(sample.position);w.vec(sample.normal);w.u32(sample.node);}
    for(auto g:m.gains)w.f32(g);
}
[[noreturn]] void invalidModel(const std::string& what){throw std::invalid_argument("modal model: "+what);}
bool canonicalFinite(double v){return std::isfinite(v)&&!(v==0&&std::signbit(v));}
bool canonicalFinite(Vec3 v){return canonicalFinite(double(v.x))&&canonicalFinite(double(v.y))&&canonicalFinite(double(v.z));}
}

void validateModalModel(const ModalModel& m){
    const auto& t=m.material;
    for(double v:{t.densityKgPerM3,t.youngModulusPa,t.poissonRatio,t.rayleighAlphaPerSecond,t.rayleighBetaSeconds,m.latticeSpacingMeters,m.minimumFrequencyHz,m.maximumFrequencyHz,m.latticeMassKg})
        if(!canonicalFinite(v))invalidModel("non-finite or negative-zero scalar");
    checkMaterial(t,"modal model");
    checkBand(m.minimumFrequencyHz,m.maximumFrequencyHz,"modal model");
    if(!(m.latticeSpacingMeters>0))invalidModel("lattice spacing must be positive");
    if(m.latticeNodes<4||std::uint64_t(m.latticeNodes)*3>ModalAnalysisSettings::MaximumDegreesOfFreedom)invalidModel("lattice node count outside 4..1200");
    if(m.latticeTetrahedra<5||m.latticeTetrahedra%5)invalidModel("lattice tetrahedron count is not a positive multiple of 5");
    if(!(m.latticeMassKg>0))invalidModel("lattice mass must be positive");
    if(!canonicalFinite(m.boundsMinimum)||!canonicalFinite(m.boundsMaximum)||m.boundsMinimum.x>m.boundsMaximum.x||m.boundsMinimum.y>m.boundsMaximum.y||m.boundsMinimum.z>m.boundsMaximum.z)invalidModel("invalid bounds");
    if(m.modes.empty()||m.modes.size()>ModalAnalysisSettings::MaximumModes)invalidModel("mode count outside 1..256");
    if(m.samples.empty()||m.samples.size()>ModalAnalysisSettings::MaximumSurfaceSamples)invalidModel("surface sample count outside 1..1024");
    double previous=0;
    for(const auto& mode:m.modes){
        if(!canonicalFinite(mode.angularFrequency)||!canonicalFinite(mode.dampingRatio)||!canonicalFinite(mode.surfaceRms))invalidModel("non-finite mode");
        if(mode.frequencyHz!=mode.angularFrequency/(2*std::numbers::pi))invalidModel("mode frequency is not angularFrequency / 2 pi");
        if(mode.frequencyHz<m.minimumFrequencyHz||mode.frequencyHz>m.maximumFrequencyHz)invalidModel("mode outside the model band");
        if(mode.angularFrequency<previous)invalidModel("modes are not in ascending frequency");
        if(mode.dampingRatio<0||!(mode.dampingRatio<1))invalidModel("damping ratio outside [0, 1)");
        if(mode.surfaceRms<0)invalidModel("negative surface RMS");
        previous=mode.angularFrequency;
    }
    std::vector<std::uint32_t> nodes;
    for(const auto& sample:m.samples){
        if(!canonicalFinite(sample.position)||!canonicalFinite(sample.normal))invalidModel("non-finite or negative-zero sample");
        const double len=std::sqrt(double(sample.normal.x)*sample.normal.x+double(sample.normal.y)*sample.normal.y+double(sample.normal.z)*sample.normal.z);
        if(std::fabs(len-1)>1e-5)invalidModel("sample normal is not unit length");
        if(sample.node>=m.latticeNodes)invalidModel("sample node out of range");
        nodes.push_back(sample.node);
    }
    std::sort(nodes.begin(),nodes.end());
    if(std::adjacent_find(nodes.begin(),nodes.end())!=nodes.end())invalidModel("duplicate sample node");
    if(m.gains.size()!=m.samples.size()*m.modes.size())invalidModel("gain table size is not samples x modes");
    for(auto g:m.gains)if(!canonicalFinite(double(g)))invalidModel("non-finite or negative-zero gain");
}
std::array<unsigned char,32> modalModelIdentity(const ModalModel& m){
    validateModalModel(m);
    Writer w;writePayload(w,m);
    return sha256(w.b);
}
std::vector<unsigned char> encodeModalModel(const ModalModel& m){
    validateModalModel(m);
    Writer w;writePayload(w,m);
    const auto digest=sha256(w.b);
    w.b.insert(w.b.end(),digest.begin(),digest.end());
    return std::move(w.b);
}
ModalModel decodeModalModel(std::span<const unsigned char> bytes){
    if(bytes.size()<4+8+32)throw std::invalid_argument("modal model decode: truncated");
    const auto payload=bytes.first(bytes.size()-32);
    const auto digest=sha256(payload);
    if(!std::equal(digest.begin(),digest.end(),bytes.end()-32))throw std::invalid_argument("modal model decode: SHA-256 mismatch");
    Reader r{payload};
    r.need(4);
    if(!std::equal(Magic.begin(),Magic.end(),payload.begin()))throw std::invalid_argument("modal model decode: bad magic");
    r.at=4;
    if(r.u32()!=ModalModel::FormatVersion)throw std::invalid_argument("modal model decode: unsupported version");
    if(r.u32()!=0)throw std::invalid_argument("modal model decode: nonzero reserved flags");
    ModalModel m;
    m.material.densityKgPerM3=r.f64();m.material.youngModulusPa=r.f64();m.material.poissonRatio=r.f64();m.material.rayleighAlphaPerSecond=r.f64();m.material.rayleighBetaSeconds=r.f64();
    m.latticeSpacingMeters=r.f64();m.minimumFrequencyHz=r.f64();m.maximumFrequencyHz=r.f64();
    m.latticeNodes=r.u32();m.latticeTetrahedra=r.u32();m.latticeMassKg=r.f64();
    m.boundsMinimum=r.vec();m.boundsMaximum=r.vec();
    r.need(32);std::copy_n(payload.begin()+static_cast<std::ptrdiff_t>(r.at),32,m.sourceMeshSHA256.begin());r.at+=32;
    const std::uint32_t modeCount=r.u32(),sampleCount=r.u32();
    if(modeCount<1||modeCount>ModalAnalysisSettings::MaximumModes||sampleCount<1||sampleCount>ModalAnalysisSettings::MaximumSurfaceSamples)
        throw std::invalid_argument("modal model decode: mode or sample count out of bounds");
    const std::size_t expected=std::size_t(modeCount)*24+std::size_t(sampleCount)*28+std::size_t(modeCount)*sampleCount*4;
    if(payload.size()-r.at!=expected)throw std::invalid_argument("modal model decode: payload length does not match the counts");
    m.modes.resize(modeCount);
    for(auto& mode:m.modes){mode.angularFrequency=r.f64();mode.frequencyHz=mode.angularFrequency/(2*std::numbers::pi);mode.dampingRatio=r.f64();mode.surfaceRms=r.f64();}
    m.samples.resize(sampleCount);
    for(auto& sample:m.samples){sample.position=r.vec();sample.normal=r.vec();sample.node=r.u32();}
    m.gains.resize(std::size_t(modeCount)*sampleCount);
    for(auto& g:m.gains)g=r.f32();
    try{validateModalModel(m);}catch(const std::invalid_argument& e){throw std::invalid_argument(std::string("modal model decode: ")+e.what());}
    std::copy(digest.begin(),digest.end(),m.identity.begin());
    return m;
}
}
