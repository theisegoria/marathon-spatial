// Marathon Spatial physics audio. Everything renders offline into memory, and no
// test opens an audio device. Printed numbers are measurements; they say nothing
// about how anything sounds.
#define ACCELERATE_NEW_LAPACK 1
#include <Accelerate/Accelerate.h>
#include <AudioToolbox/AudioToolbox.h>
#include "marathon_spatial/physics/AudioAsset.hpp"
#include "marathon_spatial/physics/Contact.hpp"
#include "marathon_spatial/physics/Fracture.hpp"
#include "marathon_spatial/physics/Materials.hpp"
#include "marathon_spatial/physics/Modal.hpp"
#include "marathon_spatial/physics/SampleLayers.hpp"
#include "marathon_spatial/physics/Stream.hpp"
#include <chrono>
#include <thread>
#include <algorithm>
#include <atomic>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <new>
#include <numbers>
#include <numeric>
#include <unistd.h>
#include <vector>

using namespace mspatial::physics;
namespace fs=std::filesystem;
constexpr double Pi=std::numbers::pi;
constexpr double Fs=48000;

namespace {
int failures=0;
void check(bool ok,const char* what){if(!ok){++failures;std::printf("FAIL: %s\n",what);}}
void near(double a,double b,double tolerance,const char* what){if(!(std::fabs(a-b)<=tolerance)){++failures;std::printf("FAIL: %s: %.12g vs %.12g (tolerance %.3g)\n",what,a,b,tolerance);}}
template<class F> bool throws(F&& f){try{f();}catch(const std::exception&){return true;}return false;}
std::atomic<std::uint64_t> allocations{0};
std::atomic<bool> counting{false};
}
void* operator new(std::size_t n){if(counting.load())allocations.fetch_add(1);if(void* p=std::malloc(n?n:1))return p;throw std::bad_alloc();}
void operator delete(void* p) noexcept {std::free(p);}
void operator delete(void* p,std::size_t) noexcept {std::free(p);}

namespace {
struct Mesh {
    std::vector<Vec3> v;std::vector<std::uint32_t> i;
    [[nodiscard]] ModalMesh view() const {return {v,i};}
};
Mesh box(Vec3 s){
    Mesh m;
    for(int k=0;k<8;++k)m.v.push_back({(k&1)?s.x:0,(k&2)?s.y:0,(k&4)?s.z:0});
    m.i={0,2,3,0,3,1, 4,5,7,4,7,6, 0,1,5,0,5,4, 2,6,7,2,7,3, 0,4,6,0,6,2, 1,3,7,1,7,5};
    return m;
}
std::vector<float> renderBus(ModalSynth& s,std::size_t frames){
    std::vector<float> y(frames);std::size_t at=0;
    while(at<frames){const auto n=std::min<std::size_t>(frames-at,s.settings().maximumBlockFrames);s.render(std::span(y).subspan(at,n),{},n);at+=n;}
    return y;
}
std::vector<float> renderContact(ContactAudio& a,std::size_t frames,std::size_t block=480){
    std::vector<float> y(frames);std::size_t at=0;
    while(at<frames){const auto n=std::min(frames-at,block);a.render(std::span(y).subspan(at,n),{},n);at+=n;}
    return y;
}
double rms(std::span<const float> y){double s=0;for(float v:y)s+=double(v)*v;return std::sqrt(s/double(y.size()));}
// Hann power spectrum, N = 2^log2n, from the start of y.
std::vector<double> powerSpectrum(std::span<const float> y,int log2n){
    const std::size_t N=std::size_t(1)<<log2n;
    std::vector<float> re(N,0.f),im(N,0.f);
    for(std::size_t i=0;i<std::min(N,y.size());++i)re[i]=y[i]*float(.5-.5*std::cos(2*Pi*double(i)/double(N)));
    auto setup=vDSP_create_fftsetup(log2n,FFT_RADIX2);
    DSPSplitComplex z{re.data(),im.data()};
    vDSP_fft_zip(setup,&z,1,log2n,FFT_FORWARD);
    vDSP_destroy_fftsetup(setup);
    std::vector<double> p(N/2);
    for(std::size_t k=0;k<N/2;++k)p[k]=double(re[k])*re[k]+double(im[k])*im[k];
    return p;
}
double centroidHz(std::span<const float> y,int log2n=15){
    const auto p=powerSpectrum(y,log2n);double top=0,bottom=0;
    for(std::size_t k=1;k<p.size();++k){top+=p[k]*double(k);bottom+=p[k];}
    return top/bottom*Fs/double(std::size_t(1)<<log2n);
}
double peakHz(std::span<const float> y,int log2n,double loHz,double hiHz){
    const auto p=powerSpectrum(y,log2n);const double bin=Fs/double(std::size_t(1)<<log2n);
    std::size_t best=std::size_t(loHz/bin);
    for(std::size_t k=best;k<p.size()&&double(k)*bin<=hiHz;++k)if(p[k]>p[best])best=k;
    return double(best)*bin;
}
const ModalElastic Steel{200e9,7850,.29,6e-4};

constexpr const char* GoldenBarIdentity="805d7252e13c57c776b81c0486feb8054f84eda9b6dd44b4cf66fa2e5b5a8c0d";
// ---- Modal frequencies against analytic bars and plates ---------------------
void analyticClosedForms(){
    ModalBody bar;bar.shape=ModalShape::Bar;bar.size={.5f,.03f,.01f};bar.material=Steel;bar.maximumModes=64;
    const auto b=analyticModalModel(bar);
    std::vector<double> expected;
    for(double root:{4.730040745,7.853204624,10.995607838,14.137165491,17.278759657,20.420352245,23.561944902,26.703537556})
        for(double depth:{.01,.03}){
            const double hz=root*root*depth/(.25)*std::sqrt(200e9/(12*7850.))/(2*Pi);
            if(hz>=20&&hz<=20000)expected.push_back(hz);
        }
    std::sort(expected.begin(),expected.end());
    check(b.modes.size()==expected.size(),"bar: one mode per admitted root and plane");
    for(std::size_t k=0;k<std::min(expected.size(),b.modes.size());++k)near(b.modes[k].frequencyHz,expected[k],1e-5*expected[k],"bar mode against Euler-Bernoulli (float-size tolerance)");
    ModalBody plate;plate.shape=ModalShape::Plate;plate.size={.6f,.4f,.003f};plate.material=Steel;plate.maximumModes=256;
    const auto p=analyticModalModel(plate);
    const double D=200e9*.003*.003*.003/(12*(1-.29*.29));
    std::vector<double> exact;
    for(int m=1;m<=16;++m)for(int n=1;n<=16;++n){
        const double hz=Pi/2*(m*m/(.6*.6)+n*n/(.4*.4))*std::sqrt(D/(7850*.003));
        if(hz<=20000)exact.push_back(hz);
    }
    std::sort(exact.begin(),exact.end());
    exact.resize(std::min<std::size_t>(exact.size(),p.modes.size()));
    check(!exact.empty(),"plate has modes");
    for(std::size_t k=0;k<exact.size();++k)near(p.modes[k].frequencyHz,exact[k],1e-5*exact[k],"plate mode against Kirchhoff simply supported");
    // Mass-normalized shape: sum over samples of phi^2 approximates 1 / mass in the mean.
    double square=0;for(std::size_t s=0;s<p.samples.size();++s)square+=double(p.gain(s,0))*p.gain(s,0);
    near(square/double(p.samples.size())*p.latticeMassKg,1,.02,"plate mode 1 is mass-normalized over the samples");
    std::printf("analytic: bar %zu modes %.1f..%.1f Hz, plate %zu modes %.1f..%.1f Hz\n",b.modes.size(),b.modes.front().frequencyHz,b.modes.back().frequencyHz,p.modes.size(),p.modes.front().frequencyHz,p.modes.back().frequencyHz);
}
void femBar(ModalAnalysisResult& out){
    const double L=.3,b=.05,t=.025,h=t/3;
    const auto mesh=box({float(L),float(b),float(t)});
    const ModalMaterial panel{1800,3e9,.25,5,2e-8};
    ModalAnalysisSettings s;s.latticeSpacingMeters=h;s.maximumModes=64;
    out=analyzeModal(mesh.view(),panel,s);
    near(out.report.latticeMassKg,1800*L*b*t,1e-9,"FEM bar mass");
    check(out.report.rigidSubspaceResidual<1e-6,"FEM bar: six rigid modes span the rigid subspace");
    double worst=0;
    for(double depth:{t,b})for(double root:{4.730040745,7.853204624}){
        const double hz=root*root*depth/(L*L)*std::sqrt(3e9/(12*1800.))/(2*Pi);
        double nearest=1e300;
        for(const auto& m:out.model.modes)if(std::fabs(m.frequencyHz/hz-1)<std::fabs(nearest))nearest=m.frequencyHz/hz-1;
        std::printf("FEM bar: depth %.3f root %.3f Euler-Bernoulli %.1f Hz, nearest FEM %+.2f%%\n",depth,root,hz,100*nearest);
        worst=std::max(worst,std::fabs(nearest));
    }
    check(std::fabs(out.model.modes[0].frequencyHz/(4.730040745*4.730040745*t/(L*L)*std::sqrt(3e9/(12*1800.))/(2*Pi))-1)<=.10,"FEM bar first mode within 10% of Euler-Bernoulli");
    check(worst<=.20,"FEM bar first four bending modes within 20%");
}
// Free square Kirchhoff plate. No closed form exists, so the reference is a
// converged Rayleigh-Ritz solution in Legendre products, independent of the FEM.
std::vector<double> ritzFreePlateHz(double a,double thickness,double E,double rho,double nu,int N){
    const int Q=64;std::vector<double> nodes(Q),weights(Q);
    for(int i=0;i<Q;++i){
        double x=std::cos(Pi*(i+.75)/(Q+.5));
        for(int it=0;it<100;++it){
            double p0=1,p1=x;for(int n=2;n<=Q;++n){const double p=((2*n-1)*x*p1-(n-1)*p0)/n;p0=p1;p1=p;}
            const double dp=Q*(x*p1-p0)/(x*x-1),dx=p1/dp;x-=dx;
            if(std::fabs(dx)<1e-15){weights[i]=2/((1-x*x)*dp*dp);break;}
        }
        nodes[i]=x;
    }
    std::vector<std::vector<double>> P(N,std::vector<double>(Q)),D1=P,D2=P;
    for(int q=0;q<Q;++q){
        const double x=nodes[q];
        P[0][q]=1;D1[0][q]=0;D2[0][q]=0;P[1][q]=x;D1[1][q]=1;D2[1][q]=0;
        for(int n=1;n+1<N;++n){
            P[n+1][q]=((2*n+1)*x*P[n][q]-n*P[n-1][q])/(n+1);
            D1[n+1][q]=((2*n+1)*(P[n][q]+x*D1[n][q])-n*D1[n-1][q])/(n+1);
            D2[n+1][q]=((2*n+1)*(2*D1[n][q]+x*D2[n][q])-n*D2[n-1][q])/(n+1);
        }
    }
    const auto integrate=[&](const auto& f,const auto& g){std::vector<double> e(std::size_t(N)*N);for(int i=0;i<N;++i)for(int j=0;j<N;++j){double s=0;for(int q=0;q<Q;++q)s+=weights[q]*f[i][q]*g[j][q];e[std::size_t(i)*N+j]=s;}return e;};
    const auto E00=integrate(P,P),E11=integrate(D1,D1),E22=integrate(D2,D2),E20=integrate(D2,P),E02=integrate(P,D2);
    const std::size_t n=std::size_t(N)*N;std::vector<double> K(n*n),M(n*n);
    const double D=E*thickness*thickness*thickness/(12*(1-nu*nu)),stiff=D*std::pow(2/a,2),mass=rho*thickness*std::pow(a/2,2);
    const auto at=[&](const std::vector<double>& e,int i,int j){return e[std::size_t(i)*N+j];};
    for(int m=0;m<N;++m)for(int r=0;r<N;++r)for(int p=0;p<N;++p)for(int q=0;q<N;++q){
        const std::size_t row=std::size_t(m)*N+r,col=std::size_t(p)*N+q;
        K[row*n+col]=stiff*(at(E22,m,p)*at(E00,r,q)+at(E00,m,p)*at(E22,r,q)+nu*(at(E20,m,p)*at(E02,r,q)+at(E02,m,p)*at(E20,r,q))+2*(1-nu)*at(E11,m,p)*at(E11,r,q));
        M[row*n+col]=mass*at(E00,m,p)*at(E00,r,q);
    }
    __LAPACK_int itype=1,order=__LAPACK_int(n),lda=order,ldb=order,lwork=-1,info=0;double query=0;std::vector<double> w(n);
    dsygv_(&itype,"N","L",&order,K.data(),&lda,M.data(),&ldb,w.data(),&query,&lwork,&info);
    lwork=__LAPACK_int(query);std::vector<double> work(std::size_t(lwork),0.);
    dsygv_(&itype,"N","L",&order,K.data(),&lda,M.data(),&ldb,w.data(),work.data(),&lwork,&info);
    check(info==0,"Ritz plate reference solved");
    std::vector<double> hz;
    for(std::size_t i=3;i<std::min<std::size_t>(n,11);++i)hz.push_back(std::sqrt(std::max(w[i],0.))/(2*Pi)); // three rigid modes first
    return hz;
}
void femPlate(){
    const double a=.18,t=.02;const double E=3e9,rho=1800,nu=.25;
    const auto mesh=box({float(a),float(a),float(t)});
    const ModalMaterial m{rho,E,nu,0,1e-9};
    ModalAnalysisSettings s;s.latticeSpacingMeters=.01;s.maximumModes=12;s.maximumSurfaceSamples=32;
    const auto r=analyzeModal(mesh.view(),m,s);
    const auto coarse=ritzFreePlateHz(a,t,E,rho,nu,10),fine=ritzFreePlateHz(a,t,E,rho,nu,12);
    double converged=0;for(std::size_t k=0;k<coarse.size();++k)converged=std::max(converged,std::fabs(coarse[k]/fine[k]-1));
    check(converged<1e-3,"Ritz reference converged to 0.1%");
    double worst=0;
    for(std::size_t k=0;k<5;++k){
        const double e=r.model.modes[k].frequencyHz/fine[k]-1;worst=std::max(worst,std::fabs(e));
        std::printf("FEM free plate mode %zu: %.1f Hz, Kirchhoff Ritz %.1f Hz (%+.2f%%)\n",k+1,r.model.modes[k].frequencyHz,fine[k],100*e);
    }
    // Measured +11..21 %: at the DOF bound this plate has two P1 cells through its thickness, which is stiff in twist.
    check(worst<=.22,"FEM free plate first five modes within 22% of Kirchhoff (P1 elements, two cells through the thickness)");
    for(std::size_t k=0;k<5;++k)check(r.model.modes[k].frequencyHz>fine[k],"P1 elements err on the stiff side");
    near(r.model.modes[3].frequencyHz,r.model.modes[4].frequencyHz,1e-6*r.model.modes[3].frequencyHz,"the square plate's degenerate pair stays degenerate");
    near(fine[1]/fine[0],coarse[1]/coarse[0],1e-3,"Ritz ratio stable");
    std::printf("FEM plate: %u DOF, solve %.2f s; Ritz convergence %.2e\n",r.report.degreesOfFreedom,r.report.solveSeconds,converged);
}
// ---- Energy decay and exact impulse response --------------------------------
void decay(const ModalModel& source){
    ModalModel m=source;
    for(double xi:{.002,.01}){
        for(auto& mode:m.modes)mode.dampingRatio=xi;
        m.identity=modalModelIdentity(m);
        ModalSynthSettings s;s.cullAmplitude=0;s.maximumModesPerVoice=1;s.maximumBlockFrames=4096;
        ModalSynth synth(s);const auto h=synth.registerModel(std::make_shared<const ModalModel>(m));
        const auto r=synth.excite({.instance=1,.model=h,.sample=3u,.impulseNs=.01,.contactSeconds=5e-5});
        check(r.modesExcited==1,"one strongest mode");
        const double omega=[&]{double best=0,hz=0;for(std::size_t k=0;k<m.modes.size();++k){const double a=std::fabs(double(m.gain(3,k))*halfSinePulseSpectrum(m.modes[k].angularFrequency,5e-5)*m.modes[k].surfaceRms);if(a>best){best=a;hz=m.modes[k].angularFrequency;}}return hz;}();
        const auto e0=synth.voice(0).energy;
        const std::size_t frames=48000;const auto y=renderBus(synth,frames);
        const std::size_t window=std::size_t(std::ceil(20*Fs*2*Pi/omega));
        // log RMS over windows against time: the slope is -xi omega.
        std::vector<double> t,l;
        // Stop 50 dB down: below that a fast mode reaches the bottom of float.
        const double span=6/(xi*omega);
        for(std::size_t at=0;at+window<frames&&double(at)/Fs<span;at+=window){t.push_back(double(at)/Fs);l.push_back(std::log(rms(std::span(y).subspan(at,window))));}
        const std::size_t n=t.size();
        check(n>=4,"enough envelope windows");
        double tm=0,lm=0;for(std::size_t i=0;i<n;++i){tm+=t[i];lm+=l[i];}tm/=double(n);lm/=double(n);
        double top=0,bottom=0;for(std::size_t i=0;i<n;++i){top+=(t[i]-tm)*(l[i]-lm);bottom+=(t[i]-tm)*(t[i]-tm);}
        const double slope=top/bottom,expected=-xi*omega;
        std::printf("decay: xi %.3f at %.0f Hz: envelope slope %.3f /s, expected %.3f /s (%+.3f%%)\n",xi,omega/(2*Pi),slope,expected,100*(slope/expected-1));
        near(slope,expected,.01*std::fabs(expected),"envelope slope matches xi omega within 1%");
        near(synth.voice(0).energy/e0,std::exp(2*expected*double(frames)/Fs),.02*std::exp(2*expected*double(frames)/Fs),"voice energy falls as exp(-2 xi omega t)");
    }
    // Rendered output against the closed form, every mode, every sample.
    ModalSynthSettings s;s.cullAmplitude=0;s.maximumModesPerVoice=256;s.maximumBlockFrames=4096;
    ModalSynth synth(s);const auto model=std::make_shared<const ModalModel>(source);const auto h=synth.registerModel(model);
    const double J=.01,tau=2e-4;
    (void)synth.excite({.instance=1,.model=h,.sample=5u,.impulseNs=J,.contactSeconds=tau});
    const auto y=renderBus(synth,4096);
    double worst=0,peak=0;
    for(std::size_t n=0;n<y.size();++n){
        double sum=0;
        for(std::size_t k=0;k<source.modes.size();++k){
            const auto& m=source.modes[k];const double wd=m.angularFrequency*std::sqrt(1-m.dampingRatio*m.dampingRatio);
            if(wd/(2*Pi)>.45*Fs)continue;
            sum+=m.surfaceRms*m.angularFrequency*double(source.gain(5,k))*J*halfSinePulseSpectrum(m.angularFrequency,tau)/wd*std::exp(-m.dampingRatio*m.angularFrequency*double(n)/Fs)*std::sin(wd*double(n)/Fs);
        }
        worst=std::max(worst,std::fabs(sum-double(y[n])));peak=std::max(peak,std::fabs(sum));
    }
    std::printf("closed form: worst |rendered - exact| %.3g, peak %.3g\n",worst,peak);
    check(worst<=2e-6*peak,"rendered samples equal the closed-form modal sum");
}
void encoding(const ModalModel& model){
    {   // The format is frozen: this closed-form model has this name. No LAPACK is involved, so the bytes are exact.
        ModalBody b;b.shape=ModalShape::Bar;b.size={.5f,.03f,.01f};b.material=Steel;b.maximumModes=8;
        const auto id=analyticModalModel(b).identity;
        char hex[65]{};for(int i=0;i<32;++i)std::snprintf(hex+2*i,3,"%02x",id[std::size_t(i)]);
        std::printf("golden identity %s\n",hex);
        check(std::strcmp(hex,GoldenBarIdentity)==0,"canonical encoding unchanged");
    }
    const auto bytes=encodeModalModel(model);
    const auto back=decodeModalModel(bytes);
    check(back.identity==model.identity&&encodeModalModel(back)==bytes,"encoding round trip");
    auto bad=bytes;bad[bad.size()/2]^=1;
    check(throws([&]{(void)decodeModalModel(bad);}),"tampered bytes rejected");
    const auto mesh=box({.2f,.1f,.05f});ModalAnalysisSettings s;s.latticeSpacingMeters=.025;s.maximumModes=24;s.maximumSurfaceSamples=32;
    const auto m=modalMaterial(*findMaterialPreset("steel"));
    check(encodeModalModel(analyzeModal(mesh.view(),m,s).model)==encodeModalModel(analyzeModal(mesh.view(),m,s).model),"same cook, same bytes");
    check(throws([&]{Mesh open=mesh;open.i.resize(open.i.size()-3);(void)analyzeModal(open.view(),m,s);}),"open mesh rejected");
}
std::shared_ptr<const ModalModel> panelModel(){
    ModalBody b;b.shape=ModalShape::Plate;b.size={.8f,.5f,.004f};b.material=Steel;b.maximumModes=48;b.airDampingPerSecond=2;
    return std::make_shared<const ModalModel>(analyticModalModel(b));
}
// ---- Contact models ---------------------------------------------------------
std::vector<float> force(const ContactForceParameters& p,const SurfaceProfile* surface,std::size_t frames,std::uint32_t seed=7){
    ContactForce f(Fs,seed);f.set(p,surface);std::vector<float> y(frames);f.generate(y);return y;
}
void contactForces(){
    const SurfaceProfile surface(.001,3);
    {   // Profile statistics
        double square=0,lag=0,previous=surface.read(0,0);const int n=60000;const double dx=surface.spacingM();
        for(int i=1;i<n;++i){const double h=surface.read(i*dx,0);square+=h*h;lag+=h*previous;previous=h;}
        near(std::sqrt(square/n),1,.1,"surface profile unit RMS");
        near(lag/square,std::exp(-2*Pi/32),.02,"surface profile exponential correlation");
        double wide=0;for(int i=0;i<n;++i){const double h=surface.read(i*dx,.016);wide+=h*h;}
        check(wide<.5*square,"reading through a wide patch filters the profile");
    }
    ContactForceParameters slide;slide.speedMPerS=.5;slide.normalForceN=20;slide.dynamicFriction=.4;slide.fluctuation=.25;
    const auto slow=force(slide,&surface,1<<17);
    near(rms(slow),.25*.4*20,.1*.25*.4*20,"scrape force RMS is fluctuation x mu x N (10%)");
    {auto p=slide;p.normalForceN=40;near(rms(force(p,&surface,1<<17))/rms(slow),2,1e-3,"scrape force is linear in normal force");}
    auto quick=slide;quick.speedMPerS=2;
    const auto fast=force(quick,&surface,1<<17);
    const double c1=centroidHz(slow),c2=centroidHz(fast);
    std::printf("scrape: RMS %.3f N; centroid %.0f Hz at .5 m/s, %.0f Hz at 2 m/s\n",rms(slow),c1,c2);
    check(c2>1.5*c1,"scrape spectrum moves up with slip speed");
    {   // The same spot gives the same force, whatever the speed.
        auto a=slide;a.speedMPerS=.25;a.smoothingSeconds=0;auto b=slide;b.speedMPerS=.5;b.smoothingSeconds=0;
        const auto x=force(a,&surface,8000),y=force(b,&surface,4000);
        double xy=0,xx=0,yy=0;
        for(std::size_t i=0;i+1<y.size();++i){xy+=double(x[2*i+1])*y[i];xx+=double(x[2*i+1])*x[2*i+1];yy+=double(y[i])*y[i];}
        const double correlation=xy/std::sqrt(xx*yy);
        std::printf("scrape position lock: correlation %.5f\n",correlation);
        check(correlation>.98,"force follows position on the surface");
    }
    {   // Stick-slip period 2 (mu_s - mu_d) N / (k v).
        ContactForceParameters p;p.speedMPerS=.02;p.normalForceN=10;p.staticFriction=.6;p.dynamicFriction=.4;p.fluctuation=0;p.stickStiffnessNPerM=2e4;p.smoothingSeconds=0;
        const auto y=force(p,&surface,1<<16);
        const double period=2*(.6-.4)*10/(2e4*.02),hz=1/period;
        const double measured=peakHz(y,16,hz*.5,hz*1.5);
        std::printf("stick-slip: expected %.2f Hz, measured %.2f Hz\n",hz,measured);
        near(measured,hz,.03*hz+Fs/65536,"stick-slip fundamental");
        auto faster=p;faster.speedMPerS=2;
        check(rms(force(faster,&surface,1<<14))<.05*rms(y),"stick-slip stops above the Stribeck speed");
    }
    {   // Rolling
        ContactForceParameters r;r.kind=ContactKind::Roll;r.speedMPerS=.5;r.normalForceN=5;r.fluctuation=.1;r.radiusM=.02;r.effectiveModulusPa=effectiveModulusPa(Steel,Steel);
        const auto y=force(r,&surface,1<<17);
        auto heavy=r;heavy.normalForceN=5000;
        const double scrape=centroidHz(slow),roll=centroidHz(y),rollHeavy=centroidHz(force(heavy,&surface,1<<17));
        std::printf("rolling: centroid %.0f Hz (5 N), %.0f Hz (5 kN); scrape at the same speed %.0f Hz\n",roll,rollHeavy,scrape);
        check(rollHeavy<.8*roll,"a larger Hertz patch filters more of the surface");
        check(roll<=scrape,"rolling is no brighter than scraping over one surface");
        auto flat=r;flat.outOfRound=.5;flat.flats=3;flat.fluctuation=.01;flat.smoothingSeconds=0;
        const double hz=3*.5/(2*Pi*.02);
        const double measured=peakHz(force(flat,&surface,1<<17),17,hz*.5,hz*1.5);
        std::printf("rolling flats: expected %.2f Hz, measured %.2f Hz\n",hz,measured);
        near(measured,hz,2*Fs/131072,"out-of-round rolling beats at flats x turn rate");
        auto hard=r;hard.fluctuation=4;
        const auto opened=force(hard,&surface,1<<15);
        check(*std::min_element(opened.begin(),opened.end())>=-5.0001f,"contact force never pulls");
    }
}
ContactAudioSettings audioSettings(double cull=1e-6){
    ContactAudioSettings s;s.synth.maximumVoices=16;s.synth.maximumModesPerVoice=48;s.synth.maximumBlockFrames=1024;s.synth.cullAmplitude=cull;return s;
}
void contactVoices(){
    const auto model=panelModel();const auto surface=std::make_shared<const SurfaceProfile>(.001,9);
    const auto steel=*findMaterialPreset("steel");
    const auto play=[&](bool first,bool second,std::size_t block,bool hit){
        ContactAudio a(audioSettings(0));const auto h=a.registerModel(model);const auto s=a.registerSurface(surface);
        ContactUpdate u;u.model=h;u.surface=s;u.force.normalForceN=20;u.force.speedMPerS=.7;u.force.fluctuation=steel.texture.slideFluctuation;
        if(first){u.contact=11;u.sample=3u;a.update(u);}
        if(second){u.contact=12;u.sample=9u;u.force.kind=ContactKind::Roll;u.force.fluctuation=.05;a.update(u);}
        if(hit)a.impact({.instance=99,.model=h,.sample=1u,.impulseNs=.05,.contactSeconds=1e-4},.0312);
        return renderContact(a,24000,block);
    };
    const auto both=play(true,true,480,true),one=play(true,false,480,false),two=play(false,true,480,false),strike=play(false,false,480,true);
    double worst=0,peak=0;
    for(std::size_t i=0;i<both.size();++i){worst=std::max(worst,std::fabs(double(both[i])-one[i]-two[i]-strike[i]));peak=std::max(peak,std::fabs(double(both[i])));}
    std::printf("superposition: worst residual %.3g, peak %.3g\n",worst,peak);
    check(peak>0&&worst<=1e-5*peak,"two contacts on one body and an impact are linear");
    check(play(true,true,480,true)==both,"contact audio repeats bit for bit");
    check(play(true,true,64,true)==both&&play(true,true,1024,true)==both&&play(true,true,333,true)==both,"bits do not depend on the render block size");
    {   // Culling on: still independent of the block size.
        const auto run=[&](std::size_t block){ContactAudio a(audioSettings());const auto h=a.registerModel(model);a.impact({.instance=9,.model=h,.sample=1u,.impulseNs=.05,.contactSeconds=1e-4},.01);return renderContact(a,96000,block);};
        check(run(480)==run(77),"culling does not follow the block size");
    }
    {   // Sample accuracy.
        ContactAudio a(audioSettings(0));const auto h=a.registerModel(model);
        a.impact({.instance=1,.model=h,.sample=1u,.impulseNs=.05,.contactSeconds=1e-4},1000/Fs);
        const auto y=renderContact(a,4096,256);
        check(std::all_of(y.begin(),y.begin()+1000,[](float v){return v==0.f;}),"silence before the hit's frame");
        check(std::fabs(y[1000])<1e-9f&&std::fabs(y[1001])>1e-5f,"a hit starts on its frame (the sine is 0 at frame 1000)");
    }
    {   // End: the force ramps out, the slot frees, the body rings down.
        ContactAudio a(audioSettings());const auto h=a.registerModel(model);const auto s=a.registerSurface(surface);
        ContactUpdate u;u.contact=5;u.model=h;u.surface=s;u.sample=3u;u.force.normalForceN=30;u.force.speedMPerS=1;
        a.update(u);
        const auto held=renderContact(a,24000);
        check(a.activeContacts()==1,"contact holds its slot");
        a.end(5,.01);
        const auto tail=renderContact(a,48000*6);
        check(a.activeContacts()==0,"contact frees after its ramp");
        check(rms(std::span(tail).last(4800))<.05*rms(std::span(held).last(4800)),"the body rings down after the contact ends");
        check(a.force(5)==nullptr,"released contact has no generator");
    }
    {   // Slots run out cleanly.
        ContactAudioSettings st=audioSettings();st.maximumContacts=2;
        ContactAudio a(st);const auto h=a.registerModel(model);const auto s=a.registerSurface(surface);
        ContactUpdate u;u.model=h;u.surface=s;u.force.normalForceN=30;u.force.speedMPerS=1;u.sample=0u;
        u.contact=1;check(a.update(u),"contact 1");u.contact=2;check(a.update(u),"contact 2");u.contact=3;check(!a.update(u),"third contact refused");
        check(a.stats().contactsRefused==1,"refusal counted");
    }
    {   // Impact solution.
        const auto r=solveImpact({.normalSpeedMPerS=2,.massAKg=.5,.massBKg=std::numeric_limits<double>::infinity(),.radiusM=.02,.restitution=.6,.a=Steel,.b=Steel});
        near(r.impulseNs,1.6*.5*2,1e-12,"J = (1 + e) m v");
        const double Estar=1/(2*(1-.29*.29)/200e9);
        near(r.contactSeconds,2.87*std::pow(.5*.5/(.02*Estar*Estar*2),.2),1e-12,"Hertz time");
        near(hertzPatchRadiusM(10,.02,Estar),std::cbrt(3*10*.02/(4*Estar)),1e-15,"Hertz patch");
    }
    {   // No allocation after setup: updates, hits and render.
        ContactAudio a(audioSettings());const auto h=a.registerModel(model);const auto s=a.registerSurface(surface);
        const auto layerSample=std::make_shared<MonoSample>();layerSample->sampleRateHz=Fs;layerSample->frames.assign(4800,0.f);
        for(std::size_t i=0;i<4800;++i)layerSample->frames[i]=float(std::sin(.2*double(i)));
        std::vector<SampleLayer> layers(1);layers[0].kind=LayerKind::Granular;layers[0].sample=layerSample;
        LayerInstrument instrument(std::move(layers),{});instrument.trigger();
        std::vector<float> bus(512);std::vector<float> voices(std::size_t(16)*512);std::vector<float> mono(512);
        ContactUpdate u;u.contact=4;u.model=h;u.surface=s;u.localPosition={0,0,.002f};u.force.normalForceN=30;u.force.speedMPerS=1;
        counting=true;allocations=0;
        for(int i=0;i<400;++i){
            u.force.speedMPerS=1+.5*std::sin(i*.1);u.force.kind=i%2?ContactKind::Slide:ContactKind::Roll;
            a.update(u);
            if(i%7==0)a.impact({.instance=std::uint64_t(50+i%3),.model=h,.localPosition={0,0,.002f},.impulseNs=.01,.contactSeconds=1e-4},.001*(i%5));
            a.render(bus,voices,512);
            instrument.setParameter(ParameterSpeed,i*.01);instrument.render(mono);
            if(i==300)a.end(4);
        }
        counting=false;
        check(allocations==0,"update, impact, end and render do not allocate");
        std::printf("allocations during 400 audio blocks: %llu\n",(unsigned long long)allocations.load());
    }
}
void fracture(){
    ContactAudio a(audioSettings());
    const auto& glass=*findMaterialPreset("glass");const auto& concrete=*findMaterialPreset("concrete");
    const auto bank=makeDebrisBank(a,glass,1e-4,.2,6);
    check(bank.entries.size()==6,"bank size");
    FractureEvent e;e.id=77;e.parentMassKg=2;e.energyJ=20;e.fragments=300;e.heightM=1.2;e.floor=concrete.elastic;e.seed=4;
    const auto plan=planFracture(e,bank);
    near(std::accumulate(plan.fragmentMassKg.begin(),plan.fragmentMassKg.end(),0.),2,1e-9,"fragments conserve mass");
    double kinetic=0;for(std::size_t i=0;i<plan.fragmentMassKg.size();++i)kinetic+=.5*plan.fragmentMassKg[i]*plan.fragmentSpeedMPerS[i]*plan.fragmentSpeedMPerS[i];
    near(kinetic,.4*20,1e-9,"fragment kinetic energy is the released share");
    check(std::is_sorted(plan.hits.begin(),plan.hits.end(),[](const DebrisHit& x,const DebrisHit& y){return x.timeSeconds<y.timeSeconds;}),"hits ascend in time");
    {   // N(>m) ~ m^-b: slope of the rank-size curve, away from the truncated ends.
        auto m=plan.fragmentMassKg;std::sort(m.begin(),m.end(),std::greater<>());
        const double b=std::log(250./25.)/std::log(m[24]/m[249]);
        std::printf("fracture: %zu hits over %.2f s, mass exponent measured %.3f (set %.2f)\n",plan.hits.size(),plan.hits.back().timeSeconds,b,e.massExponent);
        near(b,e.massExponent,.25,"fragment size distribution follows the power law");
    }
    {   // Bounces of one fragment lose speed by e, and the flight time between them is 2 v / g.
        std::vector<DebrisHit> one;for(const auto& h:plan.hits)if(h.fragment==0)one.push_back(h);
        check(one.size()>=2,"fragment bounces");
        for(std::size_t i=1;i<one.size();++i){
            near(one[i].speedMPerS,one[i-1].speedMPerS*glass.restitution,1e-12,"restitution per bounce");
            near(one[i].timeSeconds-one[i-1].timeSeconds,2*one[i].speedMPerS/e.gravityMPerS2,1e-12,"ballistic flight between bounces");
        }
    }
    const auto again=planFracture(e,bank);
    check(again.hits.size()==plan.hits.size()&&std::equal(again.hits.begin(),again.hits.end(),plan.hits.begin(),[](const DebrisHit& x,const DebrisHit& y){return x.timeSeconds==y.timeSeconds&&x.impulseNs==y.impulseNs&&x.sample==y.sample;}),"fracture plan is deterministic");
    check(scheduleFracture(a,e,plan,bank)==plan.hits.size(),"plan fits the queue");
    std::size_t peakVoices=0;std::vector<float> bus(512);double peak=0;std::vector<float> all;
    for(int i=0;i<300;++i){a.render(bus,{},512);peakVoices=std::max<std::size_t>(peakVoices,a.synth().activeVoices());for(float v:bus)peak=std::max(peak,double(std::fabs(v)));all.insert(all.end(),bus.begin(),bus.end());}
    std::printf("fracture render: played %llu, refused %llu, peak voices %zu, peak %.3g\n",(unsigned long long)a.stats().impactsPlayed,(unsigned long long)a.stats().impactsRefused,peakVoices,peak);
    check(peakVoices<=bank.entries.size(),"a whole cluster stays within one voice per bank entry");
    check(a.stats().impactsPlayed>.5*double(plan.hits.size()),"most hits sound");
    check(std::isfinite(peak)&&peak>0,"fracture renders");
    FractureEvent bad=e;bad.fragments=0;check(throws([&]{(void)planFracture(bad,bank);}),"no fragments rejected");
}
void presets(){
    check(materialPresets().size()>=16,"preset table");
    for(const auto& p:materialPresets()){
        const auto m=modalMaterial(p);
        validateModalMaterial(m,p.name);
        const auto model=analyticModalModel(modalBodyForBox({.3f,.2f,float(p.panelThicknessM)},0,p.elastic));
        check(!model.modes.empty(),"preset has modes");
    }
    check(findMaterialPreset("glass")&&!findMaterialPreset("unobtainium"),"preset lookup");
    const auto& s=*findMaterialPreset("steel");const auto m=modalMaterial(s);
    for(double hz:{200.,5000.}){const double w=2*Pi*hz;near((m.rayleighAlphaPerSecond/w+m.rayleighBetaSeconds*w)/2,s.elastic.lossFactor/2,1e-12,"loss factor anchors at 200 Hz and 5 kHz");}
}
// ---- Bring your own audio ---------------------------------------------------
template<class T> void put(std::ofstream& o,T v,bool big=false){unsigned char b[sizeof(T)];std::memcpy(b,&v,sizeof(T));if(big)std::reverse(b,b+sizeof(T));o.write(reinterpret_cast<const char*>(b),sizeof(T));}
std::vector<std::int16_t> tone(double hz,double rate,std::size_t frames){std::vector<std::int16_t> y(frames);for(std::size_t i=0;i<frames;++i)y[i]=std::int16_t(std::lround(20000*std::sin(2*Pi*hz*double(i)/rate)));return y;}
void writeWav(const fs::path& p,const std::vector<std::int16_t>& y,std::uint32_t rate,std::uint32_t start,std::uint32_t end){
    std::ofstream o(p,std::ios::binary);
    const std::uint32_t data=std::uint32_t(y.size()*2),smpl=36+24;
    o.write("RIFF",4);put<std::uint32_t>(o,4+(8+16)+(8+smpl)+(8+data));o.write("WAVE",4);
    o.write("fmt ",4);put<std::uint32_t>(o,16);put<std::uint16_t>(o,1);put<std::uint16_t>(o,1);put<std::uint32_t>(o,rate);put<std::uint32_t>(o,rate*2);put<std::uint16_t>(o,2);put<std::uint16_t>(o,16);
    o.write("smpl",4);put<std::uint32_t>(o,smpl);
    for(int i=0;i<7;++i)put<std::uint32_t>(o,0);
    put<std::uint32_t>(o,1);put<std::uint32_t>(o,0); // one loop, no sampler data
    put<std::uint32_t>(o,0);put<std::uint32_t>(o,0);put<std::uint32_t>(o,start);put<std::uint32_t>(o,end);put<std::uint32_t>(o,0);put<std::uint32_t>(o,0);
    o.write("data",4);put<std::uint32_t>(o,data);
    for(auto v:y)put(o,v);
}
void writeAiff(const fs::path& p,const std::vector<std::int16_t>& y,std::uint32_t start,std::uint32_t end){
    std::ofstream o(p,std::ios::binary);
    const std::uint32_t data=8+std::uint32_t(y.size()*2);
    const std::uint32_t mark=2+2*(2+4+6);        // two markers with five-letter names ("begin", " end ")
    const std::uint32_t total=4+(8+18)+(8+mark)+(8+20)+(8+data);
    o.write("FORM",4);put<std::uint32_t>(o,total,true);o.write("AIFF",4);
    o.write("COMM",4);put<std::uint32_t>(o,18,true);put<std::int16_t>(o,1,true);put<std::uint32_t>(o,std::uint32_t(y.size()),true);put<std::int16_t>(o,16,true);
    const unsigned char rate[10]={0x40,0x0E,0xBB,0x80,0,0,0,0,0,0}; // 48000 as an 80-bit extended
    o.write(reinterpret_cast<const char*>(rate),10);
    o.write("MARK",4);put<std::uint32_t>(o,mark,true);put<std::uint16_t>(o,2,true);
    put<std::int16_t>(o,1,true);put<std::uint32_t>(o,start,true);o.put(5);o.write("begin",5);
    put<std::int16_t>(o,2,true);put<std::uint32_t>(o,end,true);o.put(5);o.write("final",5);
    o.write("INST",4);put<std::uint32_t>(o,20,true);
    for(int i=0;i<6;++i)o.put(i==0?60:(i==3||i==5?127:0));
    put<std::int16_t>(o,0,true);put<std::int16_t>(o,1,true);put<std::int16_t>(o,1,true);put<std::int16_t>(o,2,true);put<std::int16_t>(o,0,true);put<std::int16_t>(o,0,true);put<std::int16_t>(o,0,true);
    o.write("SSND",4);put<std::uint32_t>(o,data,true);put<std::uint32_t>(o,0,true);put<std::uint32_t>(o,0,true);
    for(auto v:y)put(o,v,true);
}
// CAF goes through Apple's writer, so the region parser is tested against the OS and not our own bytes.
bool writeCaf(const fs::path& p,const std::vector<std::int16_t>& y,double rate,std::uint32_t start,std::uint32_t end){
    AudioStreamBasicDescription f{};f.mSampleRate=rate;f.mFormatID=kAudioFormatLinearPCM;f.mFormatFlags=kAudioFormatFlagIsSignedInteger|kAudioFormatFlagIsPacked;f.mChannelsPerFrame=1;f.mBitsPerChannel=16;f.mFramesPerPacket=1;f.mBytesPerFrame=f.mBytesPerPacket=2;
    const auto text=p.string();auto url=CFURLCreateFromFileSystemRepresentation(nullptr,reinterpret_cast<const UInt8*>(text.data()),CFIndex(text.size()),false);
    AudioFileID file{};const auto status=AudioFileCreateWithURL(url,kAudioFileCAFType,&f,kAudioFileFlags_EraseFile,&file);CFRelease(url);
    if(status!=noErr){std::printf("CAF create status %d\n",int(status));return false;}
    UInt32 bytes=UInt32(y.size()*2);
    bool ok=AudioFileWriteBytes(file,false,0,&bytes,y.data())==noErr;
    std::vector<char> storage(sizeof(AudioFileRegionList)+sizeof(AudioFileMarker)*2);
    auto* list=reinterpret_cast<AudioFileRegionList*>(storage.data());
    list->mSMPTE_TimeType=0;list->mNumberRegions=1;
    auto& region=list->mRegions[0];
    region.mRegionID=1;region.mName=CFSTR("loop");region.mFlags=kAudioFileRegionFlag_LoopEnable|kAudioFileRegionFlag_PlayForward;region.mNumberMarkers=2;
    region.mMarkers[0]={};region.mMarkers[0].mFramePosition=start;region.mMarkers[0].mMarkerID=1;region.mMarkers[0].mName=CFSTR("a");
    region.mMarkers[1]={};region.mMarkers[1].mFramePosition=end;region.mMarkers[1].mMarkerID=2;region.mMarkers[1].mName=CFSTR("b");
    const auto set=AudioFileSetProperty(file,kAudioFilePropertyRegionList,UInt32(storage.size()),list);
    if(set!=noErr){std::printf("CAF region status %d\n",int(set));ok=false;}
    AudioFileClose(file);
    return ok;
}
void assets(const fs::path& dir){
    const auto y=tone(440,48000,4800);
    const auto compare=[&](const AudioAsset& a,const char* what){
        check(a.channels==1&&a.frames()==y.size(),what);
        double worst=0;for(std::size_t i=0;i<std::min<std::size_t>(y.size(),a.samples.size());++i)worst=std::max(worst,std::fabs(double(a.samples[i])-y[i]/32768.));
        check(worst<=1e-6,what);
    };
    const auto wav=dir/"t.wav",aiff=dir/"t.aiff",caf=dir/"t.caf",wav44=dir/"r.wav";
    writeWav(wav,y,48000,1000,3999);
    auto a=loadAudioFile(wav);compare(a,"WAV samples");
    check(a.container=="wav"&&a.loop&&a.loop->startFrame==1000&&a.loop->endFrame==4000,"WAV smpl loop (inclusive end)");
    writeAiff(aiff,y,1200,3600);
    a=loadAudioFile(aiff);compare(a,"AIFF samples");
    check(a.container=="aiff"&&a.loop&&a.loop->startFrame==1200&&a.loop->endFrame==3600,"AIFF INST and MARK loop");
    check(writeCaf(caf,y,48000,500,4500),"CAF written by AudioToolbox");
    a=loadAudioFile(caf);compare(a,"CAF samples");
    check(a.container=="caf"&&a.loop&&a.loop->startFrame==500&&a.loop->endFrame==4500,"CAF region loop");
    writeWav(wav44,tone(1000,44100,44100),44100,4410,39689);
    const auto r=loadAudioFile(wav44,{.targetSampleRateHz=48000});
    check(r.sampleRateHz==48000&&std::llabs(std::int64_t(r.frames())-48000)<=4,"resampled length");
    near(peakHz(std::span<const float>(r.samples).subspan(4096),15,500,2000),1000,3,"resampling keeps the pitch");
    check(r.loop&&r.loop->startFrame==4800&&r.loop->endFrame==43200,"loop points follow the resample");
    check(throws([&]{(void)loadAudioFile(dir/"missing.wav");}),"a missing file throws");
    {std::ofstream junk(dir/"junk.wav",std::ios::binary);junk<<"not audio at all, just some bytes";}
    check(throws([&]{(void)loadAudioFile(dir/"junk.wav");}),"junk throws");
}
std::shared_ptr<const MonoSample> sineSample(double hz,std::size_t frames,double rate=Fs,double amplitude=.5){
    auto s=std::make_shared<MonoSample>();s->sampleRateHz=rate;s->frames.resize(frames);
    for(std::size_t i=0;i<frames;++i)s->frames[i]=float(amplitude*std::sin(2*Pi*hz*double(i)/rate));
    return s;
}
std::vector<float> renderLayers(LayerInstrument& l,std::size_t frames,std::size_t block=480){
    std::vector<float> y(frames);for(std::size_t at=0;at<frames;at+=block){const auto n=std::min(block,frames-at);l.render(std::span(y).subspan(at,n));}return y;
}
void layers(){
    near(evaluateCurve({0,10,0,1,2,false},5),.25,1e-12,"power curve");
    near(evaluateCurve({0,1,100,10000,1,true},.5),1000,1e-9,"log curve");
    near(evaluateCurve({0,1,1,3},7),3,1e-12,"curves clamp");
    {   // One-shot length and pitch
        std::vector<SampleLayer> v(1);v[0].sample=sineSample(1000,4800);v[0].pitch=2;
        LayerInstrument l(std::move(v),{});check(!l.finished(),"not started is not finished");l.trigger();
        const auto y=renderLayers(l,4800);
        check(l.finished(),"one-shot ends");
        check(rms(std::span(y).first(2300))>.3&&rms(std::span(y).last(2300))<1e-6,"pitch 2 plays in half the time");
        near(peakHz(std::span(y).first(2300),12,500,4000),2000,Fs/4096,"pitch 2 doubles the frequency");
        std::vector<SampleLayer> w(1);w[0].sample=sineSample(1000,4800,24000);
        LayerInstrument half(std::move(w),{});half.trigger();
        const auto z=renderLayers(half,9700);
        near(peakHz(z,12,100,4000),1000,Fs/4096,"a 24 kHz sample keeps its pitch at 48 kHz");
        check(rms(std::span(z).subspan(9000,500))>.2&&half.finished(),"and lasts its real duration");
    }
    {   // Loop seam: a loop that does not close on itself.
        auto s=std::make_shared<MonoSample>(*sineSample(437,9600));s->loop=LoopRegion{1000,8000};
        std::vector<SampleLayer> plain(1);plain[0].kind=LayerKind::Loop;plain[0].sample=s;plain[0].crossfadeSeconds=0;
        auto fade=plain;fade[0].crossfadeSeconds=.02;
        LayerInstrument a(std::move(plain),{}),b(std::move(fade),{});a.trigger();b.trigger();
        const auto hard=renderLayers(a,48000),soft=renderLayers(b,48000);
        const auto jump=[](const std::vector<float>& y){double worst=0;for(std::size_t i=1;i<y.size();++i)worst=std::max(worst,std::fabs(double(y[i])-y[i-1]));return worst;};
        std::printf("loop: largest step %.4f hard, %.4f crossfaded (a pure tone steps at most %.4f)\n",jump(hard),jump(soft),.5*2*Pi*437/Fs);
        check(!b.finished(),"a loop keeps going");
        check(jump(hard)>3*jump(soft),"the crossfade removes the seam");
        check(jump(soft)<1.5*.5*2*Pi*437/Fs,"crossfaded loop has no step above the tone's own");
        near(rms(std::span(soft).first(9600)),rms(std::span(soft).last(9600)),.1*rms(soft),"loop level holds");
    }
    {   // Granular
        std::vector<SampleLayer> v(1);v[0].kind=LayerKind::Granular;v[0].sample=sineSample(1000,48000);
        v[0].grainDensityPerSecond=200;v[0].grainSeconds=.02;v[0].grainPitchJitter=0;v[0].grainPositionJitter=.4;
        auto copy=v;
        LayerInstrument l(std::move(v),{.seed=5});l.trigger();
        const auto y=renderLayers(l,96000);
        const double started=double(l.stats().grainsStarted);
        std::printf("granular: %.0f grains in 2 s (expected 400), RMS %.3f, peak %.0f Hz\n",started,rms(y),peakHz(y,15,100,5000));
        near(started,400,3*std::sqrt(400.)+1,"grain arrivals are Poisson at the set density");
        near(peakHz(y,15,100,5000),1000,2*Fs/32768,"grains keep the source pitch");
        LayerInstrument same(std::move(copy),{.seed=5});same.trigger();
        check(renderLayers(same,96000,37)==y,"granular bits do not depend on the block size");
        l.release(.1);
        (void)renderLayers(l,9600);
        check(l.finished(),"release finishes the cloud");
    }
    {   // Mapping: speed to gain and to grain density, a user control to pitch.
        std::vector<SampleLayer> v(2);
        v[0].kind=LayerKind::Loop;v[0].sample=sineSample(500,48000);
        v[0].mappings.push_back({ParameterSpeed,LayerTarget::Gain,{0,4,0,1},.01});
        v[0].mappings.push_back({ParameterFirstUser,LayerTarget::Pitch,{0,1,1,2},.01});
        v[1].kind=LayerKind::Granular;v[1].sample=sineSample(3000,48000);v[1].gain=0;
        v[1].mappings.push_back({ParameterForce,LayerTarget::GrainDensity,{0,100,0,500},.01});
        LayerInstrument l(std::move(v),{});l.trigger();
        l.setParameter(ParameterSpeed,1);(void)renderLayers(l,9600);const auto quiet=renderLayers(l,9600);
        l.setParameter(ParameterSpeed,3);(void)renderLayers(l,9600);const auto loud=renderLayers(l,9600);
        near(rms(loud)/rms(quiet),3,.06,"speed maps to gain");
        l.setParameter(ParameterFirstUser,1);(void)renderLayers(l,9600);
        near(peakHz(renderLayers(l,8192),13,200,3000),1000,2*Fs/8192,"a user control maps to pitch");
        check(l.stats().grainsStarted==0,"zero force, no grains");
        l.setParameter(ParameterForce,100);(void)renderLayers(l,48000);
        near(double(l.stats().grainsStarted),500,3*std::sqrt(500.)+40,"force maps to grain density");
    }
    {   // A recording drives modes: its energy reaches the body's resonances.
        ContactAudio a(audioSettings(0));const auto model=panelModel();const auto h=a.registerModel(model);
        const auto mode=model->modes[12];
        std::vector<SampleLayer> v(1);v[0].kind=LayerKind::Loop;v[0].sample=sineSample(mode.frequencyHz,48000);
        LayerInstrument source(std::move(v),{});source.trigger();
        ContactUpdate u;u.contact=3;u.model=h;u.sample=6u;u.force.kind=ContactKind::Source;u.force.source=&source;u.force.forceScaleN=40;
        check(a.update(u),"source contact starts");
        const auto y=renderContact(a,1<<16);
        check(rms(std::span(y).last(9600))>0,"source contact sounds");
        near(peakHz(y,16,20,8000),mode.frequencyHz,3*Fs/65536,"a recorded tone at a mode rings that mode");
    }
}
void stream(){
    const auto model=panelModel();const auto surface=std::make_shared<const SurfaceProfile>(.001,9);
    // The same facts through the stream give the same bus as a direct call, and each voice its own share.
    ContactAudio direct(audioSettings(0));
    {
        const auto h=direct.registerModel(model);
        direct.impact({.instance=8,.model=h,.sample=1u,.impulseNs=.05,.contactSeconds=1e-4});
        direct.impact({.instance=9,.model=h,.sample=2u,.impulseNs=.03,.contactSeconds=1e-4});
    }
    const auto expected=renderContact(direct,4800,480);
    ContactAudio audio(audioSettings(0));
    const auto h=audio.registerModel(model);const auto s=audio.registerSurface(surface);
    PhysicsStream live(audio,64);
    check(live.postImpact({.instance=8,.model=h,.sample=1u,.impulseNs=.05,.contactSeconds=1e-4}),"post impact");
    check(live.postImpact({.instance=9,.model=h,.sample=2u,.impulseNs=.03,.contactSeconds=1e-4}),"post impact");
    check(!live.postImpact({.instance=0,.model=h,.sample=1u,.impulseNs=.05}),"zero key refused on the control side");
    check(!live.postImpact({.instance=1,.model=99,.sample=1u,.impulseNs=.05}),"unknown model refused");
    ContactUpdate bad;bad.contact=1;bad.model=h;bad.sample=0u;
    check(!live.post(bad),"a slide without a surface is refused on the control side");
    check(live.droppedCommands()==3,"refusals are counted");
    std::vector<float> a(480),b(480),sum(4800,0.f);
    const std::uint32_t slots=audio.synth().settings().maximumVoices;
    std::vector<float> rest(480);
    for(std::size_t at=0;at<4800;at+=480){
        check(blockSourceCallback(&live.voice(0),a.data(),480)==480,"voice pull");
        check(blockSourceCallback(&live.voice(1),b.data(),480)==480,"voice pull");
        for(std::uint32_t v=2;v<slots;++v)(void)blockSourceCallback(&live.voice(v),rest.data(),480);
        for(std::size_t i=0;i<480;++i)sum[at+i]=a[i]+b[i];
    }
    check(live.instanceIn(0)==8&&live.instanceIn(1)==9,"voices publish their instances");
    double worst=0;for(std::size_t i=0;i<4800;++i)worst=std::max(worst,std::fabs(double(sum[i])-expected[i]));
    near(worst,0,1e-6,"voices through the stream sum to the direct bus");
    // Pulls of different sizes on each voice, as propagation lead and rate changes make them, lose and repeat nothing.
    {
        ContactAudio uneven(audioSettings(0));
        const auto m=uneven.registerModel(model);
        PhysicsStream voices(uneven,64);
        check(voices.postImpact({.instance=8,.model=m,.sample=1u,.impulseNs=.05,.contactSeconds=1e-4}),"post impact");
        check(voices.postImpact({.instance=9,.model=m,.sample=2u,.impulseNs=.03,.contactSeconds=1e-4}),"post impact");
        std::vector<float> first(4800),second(4800);
        const std::size_t sizes[]={37,480,611,5,352,480,123};
        for(std::size_t a=0,b=0,i=0;a<4800||b<4800;++i){
            const std::size_t n=std::min<std::size_t>(480,4800-a),m=std::min<std::size_t>(sizes[i%7],4800-b);
            (void)blockSourceCallback(&voices.voice(0),first.data()+a,n);a+=n;
            (void)blockSourceCallback(&voices.voice(1),second.data()+b,m);b+=m;
        }
        double residual=0;for(std::size_t i=0;i<4800;++i)residual=std::max(residual,std::fabs(double(first[i])+second[i]-expected[i]));
        near(residual,0,1e-6,"uneven pulls on each voice still sum to the direct bus");
    }
    // A block clock runs the physics with nobody pulling, and a late reader starts close to now, not at old audio.
    {
        ContactAudio clocked(audioSettings(0));
        const auto m=clocked.registerModel(model);
        PhysicsStream voices(clocked,64);
        check(voices.postImpact({.instance=8,.model=m,.sample=1u,.impulseNs=.05,.contactSeconds=1e-4}),"post impact");
        check(voices.instanceIn(0)==0,"nothing has run yet");
        blockCallback(&voices,128,0);
        check(voices.instanceIn(0)==8,"the block clock starts an impact with no voice pulling");
        std::vector<float> early(128);
        (void)blockSourceCallback(&voices.voice(0),early.data(),128);
        check(rms(early)>0,"the voice reads what the clock rendered");
        (void)blockSourceCallback(&voices.voice(0),early.data(),128);
        check(rms(early)==0,"a pull never advances a clocked synth");
        for(int i=0;i<4000;++i)blockCallback(&voices,128,0);
        std::vector<float> late(128);
        (void)blockSourceCallback(&voices.voice(0),late.data(),128);
        // Four maximum blocks of history at 1024 frames is 85 ms at 48 kHz, so a reader is never later than that.
        check(rms(late)<rms(std::span(expected).first(480)),"a late reader hears the decayed sound, not the old attack");
    }
    // Two threads. The control side floods the queue, and the render side never allocates.
    std::atomic<bool> running{true};
    std::thread control([&]{
        ContactUpdate u;u.contact=4;u.model=h;u.surface=s;u.sample=3u;u.force.normalForceN=20;u.force.speedMPerS=1;
        for(int i=0;running.load();++i){u.force.speedMPerS=1+.5*std::sin(i*.01);(void)live.post(u);if(i%64==0)(void)live.postImpact({.instance=20,.model=h,.sample=1u,.impulseNs=.01,.contactSeconds=1e-4});std::this_thread::yield();}
    });
    std::vector<float> mono(480);
    counting=true;allocations=0;
    for(int i=0;i<2000;++i)live.renderBus(mono);
    counting=false;running=false;control.join();
    // operator new is global, so the control thread is counted too, and it does not allocate either.
    check(allocations==0,"a live stream does not allocate on either side");
    check(rms(mono)>0&&std::isfinite(rms(mono)),"the stream produced sound");
    std::printf("stream: bus residual %.3g; flood dropped %llu commands\n",worst,(unsigned long long)live.droppedCommands());
    // A finished source ends through the callback.
    std::vector<SampleLayer> v(1);v[0].sample=sineSample(1000,480);
    LayerInstrument once(std::move(v),{});once.trigger();
    std::vector<float> y(480);
    check(blockSourceCallback(&once,y.data(),480)==480,"source plays");
    check(blockSourceCallback(&once,y.data(),480)==0,"a finished source ends");
}
}

// Models and surfaces registered while a render thread runs the stream, as a host builds a body's model on first sight; and drives
// posted through the stream. A full model table refuses cleanly.
void streamGrows(){
    const auto model=panelModel();
    ContactAudioSettings settings=audioSettings(0);
    settings.synth.maximumModels=8;
    ContactAudio audio(settings);
    const auto first=audio.registerModel(model);
    PhysicsStream live(audio,256);
    std::atomic<bool> running{true};
    std::atomic<double> heard{0};
    std::thread render([&]{
        std::vector<float> mono(256);
        const std::uint32_t slots=audio.synth().settings().maximumVoices;
        while(running.load(std::memory_order_acquire)){
            blockCallback(&live,256,0);
            double sum=0;
            for(std::uint32_t v=0;v<slots;++v){(void)blockSourceCallback(&live.voice(v),mono.data(),256);for(float x:mono)sum+=double(x)*x;}
            heard.store(heard.load(std::memory_order_relaxed)+sum,std::memory_order_relaxed);
        }
    });
    std::vector<ModalModelHandle> handles{first};
    for(int i=1;i<8;++i){
        handles.push_back(audio.registerModel(model));
        check(live.postImpact({.instance=std::uint64_t(100+i),.model=handles.back(),.sample=1u,.impulseNs=.05,.contactSeconds=1e-4}),"a model registered while streaming takes a hit");
        std::this_thread::sleep_for(std::chrono::milliseconds(2));
    }
    bool full=false;
    try{(void)audio.registerModel(model);}catch(const std::invalid_argument&){full=true;}
    check(full,"a full model table refuses the next model");
    ModalDrive d;d.instance=200;d.model=handles[3];d.sample=0u;d.forceRmsN=.2;d.cutoffHz=800;
    check(live.postDrive(d),"a drive posts through the stream");
    ModalDrive bad=d;bad.instance=0;
    check(!live.postDrive(bad),"a drive with no key is refused");
    std::this_thread::sleep_for(std::chrono::milliseconds(30));
    running.store(false,std::memory_order_release);
    render.join();
    check(heard.load()>0,"hits on late models and the drive sound");
    check(audio.synth().modelCount()==8,"eight models published");
}

int main(){
    analyticClosedForms();
    ModalAnalysisResult bar;
    femBar(bar);
    femPlate();
    decay(bar.model);
    encoding(bar.model);
    contactForces();
    contactVoices();
    fracture();
    presets();
    const auto dir=fs::temp_directory_path()/("mspatial-physics-"+std::to_string(getpid()));
    fs::create_directories(dir);
    assets(dir);
    fs::remove_all(dir);
    layers();
    stream();
    streamGrows();
    std::printf(failures?"marathon spatial physics: %d FAILED\n":"marathon spatial physics: all passed\n",failures);
    return failures?1:0;
}
