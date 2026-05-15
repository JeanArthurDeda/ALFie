#include <stdio.h>
#include <numeric>
#include <ranges>
#include <random>
#include <chrono>
#include <execution>
#include "rapidjson/document.h"
#include <functional>
#include <ATen/ATen.h>
#include <torch/script.h>
#define __STDC_LIB_EXT1__
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

enum Mode
{
    ExportCVS,
    TestModel
};

Mode constexpr mode = TestModel;

#pragma region Types

using f32 = float;
#define f32_min     (-FLT_MAX)
#define f32_max     (FLT_MAX)
#define f32_epsilon (FLT_EPISILON)

using f64 = double;
#define f64_min     (-DBL_MAX)
#define f64_max     (DBL_MAX)
#define f64_epsilon (DBL_EPISILON)

using i8 = char;
#define i8_min (-0x80)
#define i8_max ( 0x7f)

using u8 = unsigned char;
#define u8_min (0x0)
#define u8_max (0xff)

using i16 = short;
#define i16_min (-0x8000)
#define i16_max ( 0x7fff)

using u16 = unsigned short;
#define u16_min (0x0)
#define u16_max (0xffff)

using i32 = int;
#define i32_min (-80000000)
#define i32_max (-7fffffff)

using u32 = unsigned int;
#define u32_min (0x0)
#define u32_max (0xffffffff)

using i64 = long long;
#define i64_min (0x8000000000000000)
#define i64_max (0x7fffffffffffffff)

typedef unsigned long long u64;
#define u64_min (0x0)
#define u64_max (0xffffffffffffffff)

typedef size_t s64;
#define s64_min (0x0)
#define s64_max (0xffffffffffffffff)

#define ARRAY_SIZE(what) (sizeof (what) / sizeof (what[0]))

inline bool getSign(f32 f) { return signbit(f); }

template <class T> T min(T const a, T const b) { return a < b ? a : b; }
template <class T> T max(T const a, T const b) { return a > b ? a : b; }
template <class T> T clamp(T const v, T const minValue, T const maxValue) { return min(maxValue, max(v, minValue)); }

class Vec3 {
public:

    enum Axis {
        X = 0,
        Y,
        Z,

        AxisCount
    };

    enum SignedAxis {
        PosX = 0,
        PosY,
        PosZ,
        NegX,
        NegY,
        NegZ,

        SignedAxisCount
    };

    enum Octant {
        XPosYPosZPos = 0,
        XNegYPosZPos,
        XPosYNegZPos,
        XNegYNegZPos,
        XPosYPosZNeg,
        XNegYPosZNeg,
        XPosYNegZNeg,
        XNegYNegZNeg,

        OctantCount
    };

    Vec3() = default;

    Vec3(f32 const x, f32 const y, f32 const z) {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    Vec3(u32 const color) {
        x = static_cast<float>((color & 0xff0000) >> 16) / 255.0f;
        y = static_cast<float>((color & 0x00ff00) >> 8) / 255.0f;
        z = static_cast<float>((color & 0x0000ff) >> 0) / 255.0f;
    }

    bool operator == (Vec3 const& other) const {
        return x == other.x && y == other.y && z == other.z;
    }

    bool operator != (Vec3 const& other) const {
        return x != other.x || y != other.y || z != other.z;
    }

    Vec3 const& operator = (Vec3 const& other) {
        x = other.x;
        y = other.y;
        z = other.z;
        return *this;
    }

    Vec3 const operator - () const {
        Vec3 ret(-x, -y, -z);
        return ret;
    }

    Vec3 const& operator += (Vec3 const& other) {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }

    Vec3 const operator + (Vec3 const& other) const {
        Vec3 s(*this);
        s += other;
        return s;
    }

    Vec3 const& operator -= (Vec3 const& other) {
        x -= other.x;
        y -= other.y;
        z -= other.z;
        return *this;
    }

    Vec3 const operator - (Vec3 const& other) const {
        Vec3 s(*this);
        s -= other;
        return s;
    }

    Vec3 const& operator *= (f32 value) {
        x *= value;
        y *= value;
        z *= value;
        return *this;
    }

    Vec3 const& operator *= (Vec3 const value) {
        x *= value.x;
        y *= value.y;
        z *= value.z;
        return *this;
    }

    Vec3 const& operator /= (f32 value) {
        f32 const oov = 1.0f / value;
        x *= oov;
        y *= oov;
        z *= oov;
        return *this;
    }

    Vec3 operator / (f32 value) const {
        Vec3 v = *this;
        return v /= value;
    }

    f32 operator | (Vec3 const& other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    Vec3 operator * (f32 const scale) const {
        return Vec3(x * scale, y * scale, z * scale);
    }

    Vec3 operator *(Vec3 const& scale) const {
        return Vec3(x * scale.x, y * scale.y, z * scale.z);
    }

    Vec3 operator % (Vec3 const& other) const {
        return Vec3(y * other.z - z * other.y,
            z * other.x - x * other.z,
            x * other.y - y * other.x);
    }

    f32 length() const {
        return sqrtf(*this | *this);
    }

    f32 lengthSquared() const {
        return *this | *this;
    }

    Vec3 const& normalize() {
        f32 const ool = 1.0f / length();
        x *= ool;
        y *= ool;
        z *= ool;
        return *this;
    }

    Vec3 normalized() const {
        Vec3 ret = *this;
        return ret.normalize();
    }

    f32 operator[] (u32 index) const {
        return (&x)[index];
    }

    f32& operator[] (u32 index) {
        return (&x)[index];
    }

    f32 const* get() const {
        return &x;
    }

    f32* get() {
        return &x;
    }

    Vec3 saturated()
    {
        Vec3 r;
        r.x = clamp(x, 0.0f, 1.0f);
        r.y = clamp(y, 0.0f, 1.0f);
        r.z = clamp(z, 0.0f, 1.0f);
        return r;
    }

    u32 color(u8 const alpha) const {
        u32 const a = static_cast<u32>(alpha) << 24;
        u32 const r = static_cast<u32>(x * 255.0f);
        u32 const g = static_cast<u32>(y * 255.0f) << 8;
        u32 const b = static_cast<u32>(z * 255.0f) << 16;
        return a | r | g | b;
    }
    Axis getAxis() const {
        u32 axis = 0;
        for (u32 i = 1; i < 3; ++i) {
            if (fabs(operator[](i)) > fabs(operator[](axis))) {
                axis = i;
            }
        }
        return (Axis)axis;
    }

    u32 getNumDimensions() const {
        u32 num = 0;
        for (u32 i = 0; i < 3; ++i)
            if (operator[](i))
                num++;
        return num;
    }

    SignedAxis getSignedAxis() const {
        u32 axis = 0;
        for (u32 i = 1; i < 3; ++i) {
            if (fabs(operator[](i)) > fabs(operator[](axis))) {
                axis = i;
            }
        }
        if (getSign(operator[](axis))) {
            axis += 3;
        }
        return (SignedAxis)axis;
    }

    Octant getOctant() const {
        u32 octant = 0;
        for (u32 i = 0; i < 3; ++i)
            octant += getSign(operator[](i)) << i;
        return (Octant)octant;
    }

    void read(rapidjson::Value::ConstArray const & array)
    {
        x = array[0].GetFloat();
        y = array[1].GetFloat();
        z = array[2].GetFloat();
    }

    static Vec3 const& getAxis(Axis const axis) { return worldAxis[axis]; }
    static Vec3 const& getSignedAxis(SignedAxis const axis) { return worldAxis[axis]; }

public:
    f32 x, y, z;

    static Vec3 min;
    static Vec3 max;
    static Vec3 one;
    static Vec3 zero;
    static Vec3 xAxis;
    static Vec3 yAxis;
    static Vec3 zAxis;
    static Vec3 worldAxis[SignedAxisCount];
};

Vec3 Vec3::min(f32_min, f32_min, f32_min);
Vec3 Vec3::max(f32_max, f32_max, f32_max);
Vec3 Vec3::zero(0.0f, 0.0f, 0.0f);
Vec3 Vec3::one(1.0f, 1.0f, 1.0f);
Vec3 Vec3::xAxis(1.0f, 0.0f, 0.0f);
Vec3 Vec3::yAxis(0.0f, 1.0f, 0.0f);
Vec3 Vec3::zAxis(0.0f, 0.0f, 1.0f);
Vec3 Vec3::worldAxis[SignedAxisCount] = { Vec3::xAxis, Vec3::yAxis, Vec3::zAxis, -Vec3::xAxis, -Vec3::yAxis, -Vec3::zAxis };

class RandomGenerator
{
public:
    RandomGenerator() : random(static_cast<int>(time(nullptr))) { }
    f32 uniform() const { return static_cast<f32>(random()) / static_cast<f32>(random.max()); }
    f32 uniformRange(f32 const min, f32 const max) const { return min + uniform() * (max - min); }
    u32 uniformRange(i32 const min, i32 const max) const { return min + (max - min == 0 ? 0 : (random() % (max - min))); }
    void setSeed(u32 const seed) { random.seed(seed); }
protected:
    mutable std::mt19937 random;
} random;

constexpr f32 pi = 3.14159265358979323846f;

#pragma endregion

struct tDiskPoint { f32 x; f32 y; };
tDiskPoint sampleDisk(f32 const x, f32 const y, f32 const radius, f32 const f = 1.0f / 2.0f)
{
    tDiskPoint ret;
    auto const r = pow(random.uniform(), f) * radius;
    auto const a = random.uniformRange(0.0f, 2.0f * pi);
    ret.x = x + cos(a) * r;
    ret.y = y + sin(a) * r;
    return ret;
}

template<class T> T* alloc(u32 const size) { return reinterpret_cast<T*>(malloc(size * sizeof(T))); }

char* readJsonFile(std::string const name)
{
    FILE* f = nullptr;
    fopen_s(&f, name.c_str(), "rb");
    if (!f) return 0x0;
    fseek(f, 0, SEEK_END);
    auto const size = ftell(f);
    fseek(f, 0, SEEK_SET);
    auto data = alloc<char>(size + 1);
    memset(data, 0, size + 1);
    fread(data, size, 1, f);
    return data;
}

rapidjson::Document readJson(std::string const name)
{
    auto data = readJsonFile(name);
    rapidjson::Document d;
    d.Parse(data);
    free(data);
    return d;
}

struct tRenderSetup
{
    i32 w;
    i32 h;
    f32 cam_fov_rad;
    Vec3 cam_pos;
    Vec3 cam_forward;
    Vec3 cam_up;
    Vec3 cam_right;
    f32 ar;
    f32 hty;
    f32 htx;

    void read(std::string const json_file_name)
    {
        printf("reading json %s\n", json_file_name.c_str());
        auto const d = readJson(json_file_name);

        w = d["w"].GetInt();
        h = d["h"].GetInt();
        cam_fov_rad = d["cam_fov_rad"].GetFloat();
        cam_pos.read(d["cam_pos"].GetArray());
        cam_forward.read(d["cam_forward"].GetArray());
        cam_up.read(d["cam_up"].GetArray());
        cam_right.read(d["cam_right"].GetArray());
        ar = d["ar"].GetFloat();
        htx = d["htx"].GetFloat();
        hty = d["hty"].GetFloat();
    }
};

struct tGbuffer
{
    Vec3 p = Vec3::zero;
    Vec3 n = Vec3::zero;
    Vec3 k_d = Vec3::zero;
    Vec3 k_s = Vec3::zero;
    f32 k_r = 0.0f;
    Vec3 k_e = Vec3::zero;
    f32 k_es = 0.0f;

    void read(rapidjson::Value const &value)
    {
        p.read(value["p"].GetArray());
        n.read(value["n"].GetArray());
        k_d.read(value["k_d"].GetArray());
        k_s.read(value["k_s"].GetArray());
        k_r = value["k_r"].GetFloat();
        k_e.read(value["k_e"].GetArray());
        k_es = value["k_es"].GetFloat();
    }

    static std::vector<tGbuffer> read(tRenderSetup const& setup, std::string const json_file_name, std::string const cache_file_name)
    {
        std::vector<tGbuffer> gbuffer(setup.w * setup.h);

        FILE* f = nullptr;
        fopen_s(&f, cache_file_name.c_str(), "rb");
        if (!f)
        {
            printf("reading json %s\n", json_file_name.c_str());
            auto const d = readJson(json_file_name);
            for (int i = 0; i < setup.w * setup.h; ++i)
                gbuffer[i].read(d[i]);

            printf("writing bin cache %s\n", cache_file_name.c_str());
            fopen_s(&f, cache_file_name.c_str(), "wb");
            if (!f) return std::vector<tGbuffer>();
            fwrite(gbuffer.data(), sizeof(tGbuffer), setup.w * setup.h, f);
            fclose(f);
            return gbuffer;
        }
        printf("reading bin cache %s\n", cache_file_name.c_str());
        fread(gbuffer.data(), sizeof(tGbuffer), setup.w * setup.h, f);
        fclose(f);
        return gbuffer;
    }
};

struct tReservoir
{
    struct tSample
    {
        Vec3 wi = Vec3::zero;
        f32 pdf = 0.0f;
        f32 mis_w = 0.0f;
        f32 cos_theta = 0.0f;

        Vec3 l_pos = Vec3::zero;
        Vec3 l_nor = Vec3::zero;
        Vec3 li = Vec3::zero;

        void read(rapidjson::Value const& value)
        {
            wi.read(value["wi"].GetArray());
            pdf = value["pdf"].GetFloat();
            mis_w = value["mis_w"].GetFloat();
            cos_theta = value["cos_theta"].GetFloat();

            l_pos.read(value["l_data"]["pos"].GetArray());
            l_nor.read(value["l_data"]["nor"].GetArray());
            li.read(value["l_data"]["li"].GetArray());
        }
    };

    tSample s;
    f32 w_sum = 0.0f;
    Vec3 c_sum = Vec3::zero;
    f32 m = 0.0f;

    tReservoir const &operator += (tReservoir const& other)
    {
        w_sum += other.w_sum;
        c_sum += other.c_sum;
        m += other.m;
        if (random.uniform() * w_sum <= other.w_sum)
            s = other.s;
        return *this;
    }

    void read(rapidjson::Value const& value)
    {
        if (value.IsNull()) return;

        s.read(value["s"]);
        w_sum = value["w_sum"].GetFloat();
        c_sum.read(value["c_sun"].GetArray());
        m = value["m"].GetFloat();
    }

    static std::vector<tReservoir> read(tRenderSetup const& setup, std::string const json_file_name, std::string const cache_file_name)
    {
        std::vector<tReservoir> reservoirs(setup.w * setup.h);

        FILE* f = nullptr;
        fopen_s(&f, cache_file_name.c_str(), "rb");
        if (!f)
        {
            printf("reading json %s\n", json_file_name.c_str());
            auto const d = readJson(json_file_name);
            for (auto i = 0; i < setup.w * setup.h; ++i)
                reservoirs[i].read(d[i]);

            printf("writing bin cache %s\n", cache_file_name.c_str());
            fopen_s(&f, cache_file_name.c_str(), "wb");
            if (!f) return std::vector<tReservoir>();
            fwrite(reservoirs.data(), sizeof(tReservoir), setup.w * setup.h, f);
            fclose(f);
            return reservoirs;
        }
        printf("reading bin cache %s\n", cache_file_name.c_str());
        fread(reservoirs.data(), sizeof(tReservoir), setup.w * setup.h, f);
        fclose(f);
        return reservoirs;
    }
};

void spatial(tRenderSetup const& setup, std::vector<tReservoir> const &src, std::vector<tGbuffer> const &gbuffer, std::vector<tReservoir> & dst, i32 const radius, i32 const num,
    std::function <void (i32 const dst_ofs, tReservoir &dst_r, tGbuffer const &dst_g, i32 const src_ofs, tReservoir const &src_r, tGbuffer const& src_g)> join)
{
    assert(src.size() == dst.size());
    assert(src.size() == gbuffer.size());
    assert(src.size() == setup.w * setup.h);

    auto start = std::chrono::high_resolution_clock::now();
    auto const w = setup.w;
    auto const h = setup.h;
    for (auto i = 0; i < h; ++i)
    {
        for (auto j = 0; j < w; ++j)
        {
            auto const dst_ofs = i * w + j;
            auto const& dst_g = gbuffer[dst_ofs];
            auto& dst_r = dst[dst_ofs];

            if (dst_g.k_es != 0.0f) continue;

            for (int c = 0; c < num; ++c)
            {
                auto const s = sampleDisk(static_cast<f32>(j), static_cast<f32>(i), static_cast<f32>(radius));
                auto const ii = static_cast<i32>(s.y);
                auto const jj = static_cast<i32>(s.x);
                if (ii < 0 || ii >= h || jj < 0 || jj >= w) continue;
                auto const src_ofs = ii * w + jj;
                auto const& src_g = gbuffer[src_ofs];
                auto src_r = src[src_ofs];

                if (src_g.k_es != 0.0f || src_r.m == 0.0f) continue;

                join(dst_ofs, dst_r, dst_g, src_ofs, src_r, src_g);
            }

            dst[i * w + j] = dst_r;
        }
        auto now = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(now - start);
        if (duration.count() > 5000)
        {
            start = now;
            auto const r = static_cast<f32>(i) / static_cast<f32>(h-1);
            printf("\tspatial join %f%%\n", r * 100.0f);
            //render(setup, dst_reservoirs, gbuffer, "intermediate.png");
        }
    }
}

void render(std::vector<tReservoir> const &src, std::vector<tGbuffer> const &gbuffer, u32 * const output)
{
    assert(src.size() == gbuffer.size());

    for (auto i = 0; i < src.size(); ++i)
    {
        auto const& r = src[i];
        auto const& g = gbuffer[i];
        if (g.k_es != 0.0f)
        {
            output[i] = (g.k_e * g.k_es).saturated().color(255);
            continue;
        }
        if (r.m == 0.0f)
        {
            output[i] = 0x0ff000000;
            continue;
        }
        output[i] = (r.c_sum / r.m).saturated().color(255);
    }
}

void render(tRenderSetup const& setup, std::vector<tReservoir> const&src, std::vector<tGbuffer> const &gbuffer, std::string const png_file_name)
{
    assert(src.size() == gbuffer.size());
    assert(src.size() == setup.w * setup.h);

    auto data = alloc<u32>(static_cast<u32>(src.size()));
    if (!data) return;
    memset(data, 0, src.size() * sizeof(u32));
    render(src, gbuffer, data);
    stbi_write_png(png_file_name.c_str(), setup.w, setup.h, 4, data, setup.w * 4);
    free(data);
}

void render(tRenderSetup const& setup, std::vector<tReservoir> const&src, std::vector<tGbuffer> const&gbuffer, std::function<u32(tRenderSetup const& setup, tReservoir const& r, tGbuffer const& g)> renderFunction, std::string const pngFileName)
{
    auto data = alloc<u32>(setup.w * setup.h);
    if (!data) return;
    memset(data, 0, setup.w * setup.h * sizeof(u32));
    for (auto i = 0; i < setup.w * setup.h; ++i)
    {
        auto const& r = src[i];
        auto const& g = gbuffer[i];
        data[i] = renderFunction(setup, r, g);
    }
    stbi_write_png(pngFileName.c_str(), setup.w, setup.h, 4, data, setup.w * 4);
    free(data);
}

float compute_loss(std::vector<tGbuffer> const& gbuffer, u32 const * const d1, u32 const * const d2)
{
    auto loss = 0.0f;
    for (auto i = 0; i < gbuffer.size(); ++i)
    {
        auto const& g = gbuffer[i];
        if (g.k_es != 0.0f) continue;
        auto const c1 = Vec3(d1[i]);
        auto const c2 = Vec3(d2[i]);
        loss += (c2 - c1) | (c2 - c1);
    }
    loss /= static_cast<f32>(gbuffer.size());
    return loss;
}

struct tJoin
{
    tJoin() = default;
    tJoin(tRenderSetup const& s, tReservoir& dst_r, tGbuffer const& dst_g, i32 const src_ofs, tReservoir const& src_r, tGbuffer const& src_g)
    {
        auto const dst_wo = (s.cam_pos - dst_g.p).normalized();
        auto const src_wo = (s.cam_pos - src_g.p).normalized();
        auto const dst_h = (dst_r.s.wi + dst_wo).normalized();
        auto const src_h = (src_r.s.wi + src_wo).normalized();

        d = (dst_g.p - src_g.p).length() / 3.0f;
        n_cos = dst_g.n | src_g.n;
        h_cos = dst_h | src_h;
        d_r = dst_g.k_r;
        d_k_d = dst_g.k_d;
        s_r = src_g.k_r;
        s_k_d = src_g.k_d;

        this->src_ofs = src_ofs;
    }
    f32 d = 0.0f;
    f32 n_cos = 0.0f;
    f32 h_cos = 0.0f;
    f32 d_r = 0.0f;
    f32 s_r = 0.0f;
    Vec3 d_k_d = Vec3::zero;
    Vec3 s_k_d = Vec3::zero;

    i32 src_ofs = 0;

    f32 loss = 0.0f;

    bool operator == (tJoin const& other) const
    {
        return  d == other.d &&
            n_cos == other.n_cos &&
            h_cos == other.h_cos &&
            d_r == other.d_r &&
            s_r == other.s_r &&
            d_k_d == other.d_k_d &&
            s_k_d == other.s_k_d;
    }

    static void writeMaxAbsLoss(std::string const file_name, f32 const max_abs_loss)
    {
        FILE* f = nullptr;
        fopen_s(&f, file_name.c_str(), "wt");
        if (!f) return;
        fprintf(f, "max_abs_loss: %f\n", max_abs_loss);
        fclose(f);
    }

    static void writeCSVHeader(std::string const csv_file_name)
    {
        FILE* f = nullptr;
        fopen_s(&f, csv_file_name.c_str(), "wt");
        if (!f) return;
        fprintf(f, "d,n_cos,h_cos,d_r,d_k_d_r,d_k_d_g,d_k_d_b,s_r,s_k_d_r,s_k_d_g,s_k_d_b,loss\n");
        fflush(f);
        fclose(f);
    }

    void writeCSV(std::string const csv_file_name) const
    {
        FILE* f = nullptr;
        fopen_s(&f, csv_file_name.c_str(), "at");
        if (!f) return;
        fprintf(f, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", d, n_cos, h_cos, d_r, d_k_d.x, d_k_d.y, d_k_d.z, s_r, s_k_d.x, s_k_d.y, s_k_d.z, loss);
        fflush(f);
        fclose(f);
    }
};

class Model
{
public:
    Model(std::string const file_name)
    {
        printf("loading model %s\n\n", file_name.c_str());
        try {
            model = torch::jit::load(file_name);
        }
        catch (const c10::Error& e) {
            std::cerr << "Error loading the model: " << e.msg() << std::endl;
        }
    }

    float evaluate(f32 const d, f32 const n_cos, f32 const h_cos, f32 const d_r, f32 const s_r, Vec3 const &d_k_d, Vec3 const &s_k_d)
    {
        std::vector<float> input_data = { d, n_cos, h_cos, d_r, d_k_d.x, d_k_d.y, d_k_d.z, s_r, s_k_d.x, s_k_d.y, s_k_d.z };
        torch::Tensor input = torch::tensor(input_data).reshape({ 1, 11 }).to(torch::kFloat32);
        torch::Tensor output = model.forward({ input }).toTensor();
        return output.item<float>();
    }
protected:
    torch::jit::script::Module model;

};


int main()
{
    random.setSeed(0xCafeBabe);

    // load reference
    printf("load reference image\n");
    i32 reference_w = 0;
    i32 reference_h = 0;
    i32 reference_b = 0;
    auto const *reference = reinterpret_cast<u32*>(stbi_load("../output/cycles.png", &reference_w, &reference_h, &reference_b, 4));

    // load data
    tRenderSetup setup;
    setup.read("../output/config.json");
    auto const gbuffer = tGbuffer::read(setup, "../output/gbuffer.json", "./gbuffer.bin");
    auto const src_reservoirs = tReservoir::read(setup, "../output/reservoirs.json", "./reservoirs.bin");
    render(setup, src_reservoirs, gbuffer, "init_shadow.png");


    if (mode == ExportCVS)
    {
        // =======================================================================
        // Export Join CSV
        // =======================================================================

        auto constexpr join_radius = 10;
        auto constexpr join_num = 10;
        // do spatial and record joins
        printf("record spatial joins\n");
        auto reserved = std::vector<tJoin>();
        reserved.reserve(join_num);
        auto dst_reservoirs = src_reservoirs;
        std::vector<std::vector<tJoin>> all_joins(setup.w * setup.h, reserved);
        auto num_joins = 0;
        spatial(setup, src_reservoirs, gbuffer, dst_reservoirs, join_radius, join_num,
            [&setup, &all_joins, &num_joins](i32 const dst_ofs, tReservoir& dst_r, tGbuffer const& dst_g, i32 const src_ofs, tReservoir const& src_r, tGbuffer const& src_g)
            {
                all_joins[dst_ofs].push_back(tJoin(setup, dst_r, dst_g, src_ofs, src_r, src_g));
                dst_r += src_r;
                num_joins++;
            });
        auto with_all = alloc<u32>(setup.w * setup.h);
        render(dst_reservoirs, gbuffer, with_all);
        auto const with_all_loss = compute_loss(gbuffer, with_all, reference);
        stbi_write_png("spatial.png", setup.w, setup.h, 4, with_all, setup.w * 4);

        printf("computing excluded joins loss\n");
        auto constexpr csv_file = "./samples.csv";
        tJoin::writeCSVHeader(csv_file);

        for (auto dst_ofs = 0; dst_ofs < setup.w * setup.h; ++dst_ofs)
        {
            auto& joins = all_joins[dst_ofs];
            auto const reference_color = Vec3(reference[dst_ofs]);
            auto const with_color = with_all[dst_ofs];
            auto const oneOverWH = 1.0f / static_cast<f32>(setup.w * setup.h);
            for (auto& excluded : joins)
            {
                auto dst = src_reservoirs[dst_ofs];
                for (auto const& j : joins)
                {
                    if (excluded == j) continue;
                    dst += src_reservoirs[j.src_ofs];
                }

                auto const without_color = (dst.c_sum / dst.m).saturated();

                auto const p_without_loss = ((reference_color - without_color) | (reference_color - without_color));
                auto const p_with_loss = (reference_color - with_color) | (reference_color - with_color);
                excluded.loss = p_without_loss - p_with_loss;
            }
        }

        printf("computing min/max excluded joins loss\n");
        auto max_loss = f32_min;
        auto min_loss = f32_max;
        for (auto ofs = 0; ofs < setup.w * setup.h; ++ofs)
        {
            auto const& joins = all_joins[ofs];
            for (auto const& join : joins)
            {
                min_loss = min(min_loss, join.loss);
                max_loss = max(max_loss, join.loss);
            }
        }
        printf("min_loss %f max_loss %f\n", min_loss, max_loss);

        printf("convert to tanh loss\n");
        auto const max_abs_loss = max(fabs(min_loss), fabs(max_loss));
        printf("tanh max_abs_loss %f\n", max_abs_loss);
        tJoin::writeMaxAbsLoss("samples_max_abs_loss.txt", max_abs_loss);
        for (auto ofs = 0; ofs < setup.w * setup.h; ++ofs)
        {
            auto& joins = all_joins[ofs];
            for (auto& join : joins)
                join.loss /= max_abs_loss;
        }

        printf("write csv\n");
        auto start = std::chrono::high_resolution_clock::now();
        auto num_done = 0;
        tJoin::writeCSVHeader("samples.csv");
        for (auto ofs = 0; ofs < setup.w * setup.h; ++ofs)
        {
            auto const& joins = all_joins[ofs];
            for (auto const& join : joins)
            {
                join.writeCSV("samples.csv");
                num_done++;
            }

            auto now = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(now - start);
            if (duration.count() > 5000)
            {
                start = now;
                auto const r = static_cast<f32>(num_done) / static_cast<f32>(num_joins);
                printf("write csv %f%%\n", r * 100.0f);
                //render(setup, dst_reservoirs, gbuffer, "intermediate.png");
            }
        }
        
        return 0x0;
    }

    // =======================================================================
    // Test model
    // =======================================================================

    auto model = Model("join_predictor.pt");
    auto constexpr join_radius = 10;
    auto constexpr join_num = 10;

    // Analitic
    printf("all reservoirs join.\n");
    auto dst_reservoirs = src_reservoirs;
    spatial(setup, src_reservoirs, gbuffer, dst_reservoirs, join_radius, join_num,
        [&s = setup](i32 const dst_ofs, tReservoir& dst_r, tGbuffer const& dst_g, i32 const src_ofs, tReservoir const& src_r, tGbuffer const& src_g)
        {
            dst_r += src_r;
        });
    render(setup, dst_reservoirs, gbuffer, "join_all.png");

    // Analitic
    printf("analytic reservoirs join.\n");
    dst_reservoirs = src_reservoirs;
    spatial(setup, src_reservoirs, gbuffer, dst_reservoirs, join_radius, join_num,
        [&s = setup](i32 const dst_ofs, tReservoir& dst_r, tGbuffer const& dst_g, i32 const src_ofs, tReservoir const& src_r, tGbuffer const& src_g)
        {
            auto const dst_wo = (s.cam_pos - dst_g.p).normalized();
            auto const src_wo = (s.cam_pos - src_g.p).normalized();
            auto const dst_h = (dst_r.s.wi + dst_wo).normalized();
            auto const src_h = (src_r.s.wi + src_wo).normalized();

            //if (max(0.0f, 1.0f - abs(dst_g.k_r - src_g.k_r) * 2.0f) < 0.2f) return;
            if ((dst_g.p - src_g.p).length() > 0.1f) return;
            if ((dst_g.n | src_g.n) < 0.9f) return;
            if ((dst_h | src_h) < 0.83f) return;
            
            
            dst_r += src_r;
        });
    render(setup, dst_reservoirs, gbuffer, "join_analytic.png");

    // neural
    printf("neural reservoirs join.\n");
    dst_reservoirs = src_reservoirs;
    spatial(setup, src_reservoirs, gbuffer, dst_reservoirs, join_radius, join_num,
        [&s = setup, &model](i32 const dst_ofs, tReservoir& dst_r, tGbuffer const& dst_g, i32 const src_ofs, tReservoir const& src_r, tGbuffer const& src_g)
        {
            auto const dst_wo = (s.cam_pos - dst_g.p).normalized();
            auto const src_wo = (s.cam_pos - src_g.p).normalized();
            auto const dst_h = (dst_r.s.wi + dst_wo).normalized();
            auto const src_h = (src_r.s.wi + src_wo).normalized();

            auto const d = (dst_g.p - src_g.p).length() / 3.0f;
            auto const n_cos = dst_g.n | src_g.n;
            auto const h_cos = dst_h | src_h;
            auto const d_r = dst_g.k_r;
            auto const d_k_d = dst_g.k_d;
            auto const s_r = src_g.k_r;
            auto const s_k_d = src_g.k_d;

            auto r = model.evaluate(d, n_cos, h_cos, d_r, s_r, d_k_d, s_k_d);
            if (r < 0.05f) return;

            dst_r += src_r;
        });
    render(setup, dst_reservoirs, gbuffer, "join_neural.png");

    return 0x0;
}