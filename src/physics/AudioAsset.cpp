#include "marathon_spatial/physics/AudioAsset.hpp"
#include <AudioToolbox/AudioToolbox.h>
#include <CoreFoundation/CoreFoundation.h>
#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <memory>
#include <stdexcept>

namespace mspatial::physics {
namespace {
[[noreturn]] void fail(const std::filesystem::path& path,const std::string& what,OSStatus status=noErr){
    throw std::runtime_error("audio file "+path.string()+": "+what+(status!=noErr?" (status "+std::to_string(status)+")":""));
}
struct File {
    std::FILE* f;
    explicit File(const std::filesystem::path& p):f(std::fopen(p.c_str(),"rb")){}
    ~File(){if(f)std::fclose(f);}
    File(const File&)=delete;File& operator=(const File&)=delete;
    bool read(void* to,std::size_t n){return std::fread(to,1,n,f)==n;}
    bool skip(std::int64_t n){return fseeko(f,n,SEEK_CUR)==0;}
    std::vector<unsigned char> bytes(std::uint64_t n){
        std::vector<unsigned char> b;
        if(n>(1u<<24))return b; // no metadata chunk is this large
        b.resize(n);
        if(!read(b.data(),n))b.clear();
        return b;
    }
};
std::uint32_t le32(const unsigned char* p){return std::uint32_t(p[0])|std::uint32_t(p[1])<<8|std::uint32_t(p[2])<<16|std::uint32_t(p[3])<<24;}
std::uint16_t be16(const unsigned char* p){return std::uint16_t(p[0]<<8|p[1]);}
std::uint32_t be32(const unsigned char* p){return std::uint32_t(p[0])<<24|std::uint32_t(p[1])<<16|std::uint32_t(p[2])<<8|p[3];}
std::uint64_t be64(const unsigned char* p){return std::uint64_t(be32(p))<<32|be32(p+4);}

std::optional<LoopRegion> wavLoop(File& f){
    for(;;){
        unsigned char h[8];
        if(!f.read(h,8))return std::nullopt;
        const std::uint32_t size=le32(h+4);
        if(std::memcmp(h,"smpl",4)!=0){if(!f.skip(std::int64_t(size)+(size&1)))return std::nullopt;continue;}
        const auto b=f.bytes(size);
        if(b.size()<36+24||le32(b.data()+28)==0)return std::nullopt;
        // A WAV loop end is the last frame played, so it is inclusive.
        const std::uint32_t start=le32(b.data()+36+8),end=le32(b.data()+36+12);
        if(end<start)return std::nullopt;
        return LoopRegion{start,std::uint64_t(end)+1};
    }
}
std::optional<LoopRegion> aiffLoop(File& f){
    std::vector<std::pair<std::int16_t,std::uint32_t>> markers;
    std::int16_t begin=0,end=0;bool looped=false,marks=false;
    for(;;){
        unsigned char h[8];
        if(!f.read(h,8))break;
        const std::uint32_t size=be32(h+4);
        if(std::memcmp(h,"MARK",4)==0){
            const auto b=f.bytes(size);
            if(b.size()<2)return std::nullopt;
            std::size_t at=2;
            for(std::uint16_t i=0,n=be16(b.data());i<n;++i){
                if(at+7>b.size())return std::nullopt;
                markers.push_back({std::int16_t(be16(b.data()+at)),be32(b.data()+at+2)});
                const std::size_t text=1+b[at+6]; // count byte and text, padded to even
                at+=6+text+(text&1);
            }
            marks=true;
            if(size&1)f.skip(1);
        }else if(std::memcmp(h,"INST",4)==0){
            const auto b=f.bytes(size);
            if(b.size()<14)return std::nullopt;
            looped=be16(b.data()+8)!=0;begin=std::int16_t(be16(b.data()+10));end=std::int16_t(be16(b.data()+12));
            if(size&1)f.skip(1);
        }else if(!f.skip(std::int64_t(size)+(size&1)))break;
    }
    if(!looped||!marks)return std::nullopt;
    std::optional<std::uint32_t> a,b;
    for(const auto& m:markers){if(m.first==begin)a=m.second;if(m.first==end)b=m.second;}
    if(!a||!b||*b<=*a)return std::nullopt;
    return LoopRegion{*a,*b};
}
std::optional<LoopRegion> cafLoop(File& f){
    unsigned char skip[4];
    if(!f.read(skip,4))return std::nullopt; // version and flags
    for(;;){
        unsigned char h[12];
        if(!f.read(h,12))return std::nullopt;
        const auto size=static_cast<std::int64_t>(be64(h+4));
        if(size<0)return std::nullopt;      // a data chunk of unknown size runs to the end
        if(std::memcmp(h,"regn",4)!=0){if(!f.skip(size))return std::nullopt;continue;}
        const auto b=f.bytes(std::uint64_t(size));
        if(b.size()<8)return std::nullopt;
        constexpr std::size_t MarkerBytes=28; // sizeof(CAFMarker), packed
        std::size_t at=8;
        for(std::uint32_t r=0,n=be32(b.data()+4);r<n;++r){
            if(at+12>b.size())return std::nullopt;
            const std::uint32_t flags=be32(b.data()+at+4),count=be32(b.data()+at+8);
            at+=12;
            if(at+std::size_t(count)*MarkerBytes>b.size())return std::nullopt;
            if((flags&kCAFRegionFlag_LoopEnable)&&count>=2){
                double lo=1e300,hi=-1;
                for(std::uint32_t k=0;k<count;++k){
                    const double position=std::bit_cast<double>(be64(b.data()+at+std::size_t(k)*MarkerBytes+4));
                    lo=std::min(lo,position);hi=std::max(hi,position);
                }
                if(hi>lo&&lo>=0)return LoopRegion{std::uint64_t(lo),std::uint64_t(hi)};
            }
            at+=std::size_t(count)*MarkerBytes;
        }
        return std::nullopt;
    }
}
std::string containerOf(const char id[4],const char form[4]){
    if(std::memcmp(id,"RIFF",4)==0&&std::memcmp(form,"WAVE",4)==0)return "wav";
    if(std::memcmp(id,"FORM",4)==0&&(std::memcmp(form,"AIFF",4)==0||std::memcmp(form,"AIFC",4)==0))return "aiff";
    if(std::memcmp(id,"caff",4)==0)return "caf";
    return {};
}
}

std::optional<LoopRegion> readLoopRegion(const std::filesystem::path& path){
    File f(path);
    if(!f.f)fail(path,"cannot open");
    char h[12];
    if(!f.read(h,12))return std::nullopt;
    const auto kind=containerOf(h,h+8);
    if(kind=="wav")return wavLoop(f);
    if(kind=="aiff")return aiffLoop(f);
    if(kind=="caf"){fseeko(f.f,4,SEEK_SET);return cafLoop(f);}
    return std::nullopt;
}

AudioAsset loadAudioFile(const std::filesystem::path& path,const AudioLoadOptions& options){
    if(!std::isfinite(options.targetSampleRateHz)||options.targetSampleRateHz<0)fail(path,"invalid target sample rate");
    AudioAsset asset;
    {
        File f(path);
        if(!f.f)fail(path,"cannot open");
        char h[12]{};
        if(f.read(h,12)){asset.container=containerOf(h,h+8);if(asset.container.empty())asset.container.assign(h,4);}
    }
    const std::string native=path.string();
    auto url=CFURLCreateFromFileSystemRepresentation(kCFAllocatorDefault,reinterpret_cast<const UInt8*>(native.data()),CFIndex(native.size()),false);
    if(!url)fail(path,"cannot make a URL");
    ExtAudioFileRef file{};
    const OSStatus opened=ExtAudioFileOpenURL(url,&file);
    CFRelease(url);
    if(opened!=noErr)fail(path,"cannot decode",opened);
    const std::unique_ptr<std::remove_pointer_t<ExtAudioFileRef>,decltype(&ExtAudioFileDispose)> guard(file,&ExtAudioFileDispose);
    AudioStreamBasicDescription source{};UInt32 size=sizeof source;
    if(const auto s=ExtAudioFileGetProperty(file,kExtAudioFileProperty_FileDataFormat,&size,&source);s!=noErr)fail(path,"no format",s);
    SInt64 fileFrames=0;size=sizeof fileFrames;
    if(const auto s=ExtAudioFileGetProperty(file,kExtAudioFileProperty_FileLengthFrames,&size,&fileFrames);s!=noErr)fail(path,"no length",s);
    if(source.mChannelsPerFrame<1||source.mChannelsPerFrame>64||!(source.mSampleRate>0))fail(path,"unsupported format");
    asset.channels=source.mChannelsPerFrame;
    asset.sampleRateHz=options.targetSampleRateHz>0?options.targetSampleRateHz:source.mSampleRate;
    const double ratio=asset.sampleRateHz/source.mSampleRate;
    const auto expected=static_cast<std::uint64_t>(std::ceil(double(std::max<SInt64>(fileFrames,0))*ratio))+16;
    if(expected>options.maximumFrames)fail(path,"longer than the frame limit");
    AudioStreamBasicDescription client{};
    client.mSampleRate=asset.sampleRateHz;client.mFormatID=kAudioFormatLinearPCM;
    client.mFormatFlags=kAudioFormatFlagIsFloat|kAudioFormatFlagIsPacked|kAudioFormatFlagsNativeEndian;
    client.mChannelsPerFrame=asset.channels;client.mBitsPerChannel=32;client.mFramesPerPacket=1;
    client.mBytesPerFrame=client.mBytesPerPacket=4*asset.channels;
    if(const auto s=ExtAudioFileSetProperty(file,kExtAudioFileProperty_ClientDataFormat,sizeof client,&client);s!=noErr)fail(path,"cannot convert to float",s);
    asset.samples.resize(expected*asset.channels);
    std::uint64_t done=0;
    while(done<expected){
        AudioBufferList list{};list.mNumberBuffers=1;
        UInt32 frames=static_cast<UInt32>(std::min<std::uint64_t>(expected-done,1u<<16));
        list.mBuffers[0]={asset.channels,frames*client.mBytesPerFrame,asset.samples.data()+done*asset.channels};
        if(const auto s=ExtAudioFileRead(file,&frames,&list);s!=noErr)fail(path,"read failed",s);
        if(frames==0)break;
        done+=frames;
    }
    asset.samples.resize(done*asset.channels);
    if(auto loop=readLoopRegion(path)){
        const std::uint64_t start=static_cast<std::uint64_t>(std::llround(double(loop->startFrame)*ratio));
        const std::uint64_t end=std::min(done,static_cast<std::uint64_t>(std::llround(double(loop->endFrame)*ratio)));
        if(end>start)asset.loop=LoopRegion{start,end};
    }
    return asset;
}

MonoSample toMono(const AudioAsset& a){
    MonoSample m;m.sampleRateHz=a.sampleRateHz;m.loop=a.loop;
    const auto frames=a.frames();
    m.frames.resize(frames);
    for(std::uint64_t i=0;i<frames;++i){
        double sum=0;
        for(std::uint32_t c=0;c<a.channels;++c)sum+=a.samples[i*a.channels+c];
        m.frames[i]=float(sum/a.channels);
    }
    return m;
}
}
