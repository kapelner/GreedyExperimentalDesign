#include <Rcpp.h>
#include <vector>
#include <string>
#include <iostream>
#include <numeric>
#include <cstring>

#ifdef GED_USE_WGPU
#include <webgpu/webgpu.h>
#include <webgpu/wgpu.h>

namespace ged_wgpu_backend {

static void handle_adapter_request(WGPURequestAdapterStatus status, WGPUAdapter adapter, WGPUStringView message, void* userdata, void* userdata2) {
    if (status == WGPURequestAdapterStatus_Success) *(WGPUAdapter*)userdata = adapter;
}
static void handle_device_request(WGPURequestDeviceStatus status, WGPUDevice device, WGPUStringView message, void* userdata, void* userdata2) {
    if (status == WGPURequestDeviceStatus_Success) *(WGPUDevice*)userdata = device;
}
static void handle_buffer_map(WGPUMapAsyncStatus status, WGPUStringView message, void* userdata, void* userdata2) {
    *(bool*)userdata = (status == WGPUMapAsyncStatus_Success);
}

struct GpuContext {
    WGPUInstance instance = nullptr;
    WGPUAdapter  adapter  = nullptr;
    WGPUDevice   device   = nullptr;
    WGPUQueue    queue    = nullptr;
};
static GpuContext g_ctx;

static GpuContext& get_context() {
    if (!g_ctx.instance) {
        WGPUInstanceDescriptor idesc = {}; idesc.nextInChain = NULL;
        g_ctx.instance = wgpuCreateInstance(&idesc);
        WGPURequestAdapterCallbackInfo acb = {.callback = handle_adapter_request, .userdata1 = &g_ctx.adapter};
        wgpuInstanceRequestAdapter(g_ctx.instance, NULL, acb);
        WGPURequestDeviceCallbackInfo dcb = {.callback = handle_device_request, .userdata1 = &g_ctx.device};
        wgpuAdapterRequestDevice(g_ctx.adapter, NULL, dcb);
        g_ctx.queue = wgpuDeviceGetQueue(g_ctx.device);
    }
    return g_ctx;
}

void release_context() {
    if (g_ctx.instance) {
        wgpuQueueRelease(g_ctx.queue);
        wgpuDeviceRelease(g_ctx.device);
        wgpuAdapterRelease(g_ctx.adapter);
        wgpuInstanceRelease(g_ctx.instance);
        g_ctx = GpuContext{};
    }
}

const char* OBJECTIVE_SHADER = R"(
struct Params { n: u32, p: u32, nT: u32, numIT: u32, numIC: u32 };
@group(0) @binding(0) var<uniform> params: Params;
@group(0) @binding(1) var<storage, read> X: array<f32>;
@group(0) @binding(2) var<storage, read> Sinv: array<f32>;
@group(0) @binding(3) var<storage, read> i_Ts: array<u32>;
@group(0) @binding(4) var<storage, read> i_Cs: array<u32>;
@group(0) @binding(5) var<storage, read> avg_T: array<f32>;
@group(0) @binding(6) var<storage, read> avg_C: array<f32>;
@group(0) @binding(7) var<storage, read_write> results: array<f32>;
var<workgroup> sSinv: array<f32, 1024>;
var<workgroup> sAvgT: array<f32, 32>;
var<workgroup> sAvgC: array<f32, 32>;
@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) global_id: vec3<u32>,
        @builtin(local_invocation_id) lid: vec3<u32>) {
    let tid = lid.x; let p = params.p;
    for (var ti = tid; ti < p * p; ti += 64u) { sSinv[ti] = Sinv[ti]; }
    for (var ti = tid; ti < p; ti += 64u) { sAvgT[ti] = avg_T[ti]; sAvgC[ti] = avg_C[ti]; }
    workgroupBarrier();
    let idx = global_id.x;
    if (idx >= params.numIT * params.numIC) { return; }
    let t_idx = idx / params.numIC; let c_idx = idx % params.numIC;
    let i_T = i_Ts[t_idx]; let i_C = i_Cs[c_idx];
    let nT = f32(params.nT); let nC = f32(params.n - params.nT);
    var val = 0.0f;
    for (var i = 0u; i < p; i++) {
        let x_iT_i = X[i * params.n + i_T]; let x_iC_i = X[i * params.n + i_C];
        let diff_i = (sAvgT[i] - x_iT_i/nT + x_iC_i/nT) - (sAvgC[i] - x_iC_i/nC + x_iT_i/nC);
        var temp = 0.0f;
        for (var j = 0u; j < p; j++) {
            let x_iT_j = X[j * params.n + i_T]; let x_iC_j = X[j * params.n + i_C];
            let diff_j = (sAvgT[j] - x_iT_j/nT + x_iC_j/nT) - (sAvgC[j] - x_iC_j/nC + x_iT_j/nC);
            temp += sSinv[i * p + j] * diff_j;
        }
        val += diff_i * temp;
    }
    results[idx] = val;
}
)";

const char* BATCH_OBJECTIVE_SHADER = R"(
struct Params { n: u32, r: u32 };
@group(0) @binding(0) var<uniform> p: Params;
@group(0) @binding(1) var<storage, read> K: array<f32>;
@group(0) @binding(2) var<storage, read> W: array<f32>;
@group(0) @binding(3) var<storage, read_write> out_res: array<f32>;
var<workgroup> partial: array<f32, 256>;
var<workgroup> sw: array<f32, 256>;
@compute @workgroup_size(256)
fn main(@builtin(workgroup_id) wid: vec3<u32>, @builtin(local_invocation_id) lid: vec3<u32>) {
    let di = wid.x; if (di >= p.r) { return; }
    let tid = lid.x; let n = p.n; let r = p.r; var s = 0.0;
    for (var jt = 0u; jt < n; jt += 256u) {
        let j = jt + tid;
        sw[tid] = select(0.0, 2.0 * W[j * r + di] - 1.0, j < n);
        workgroupBarrier();
        let tile_end = min(256u, n - jt);
        for (var i = tid; i < n; i += 256u) {
            let wi = 2.0 * W[i * r + di] - 1.0;
            for (var k = 0u; k < tile_end; k++) { s += wi * K[(jt + k) * n + i] * sw[k]; }
        }
        workgroupBarrier();
    }
    partial[tid] = s; workgroupBarrier();
    for (var stride = 128u; stride > 0u; stride >>= 1u) {
        if (tid < stride) { partial[tid] += partial[tid + stride]; }
        workgroupBarrier();
    }
    if (tid == 0u) { out_res[di] = partial[0]; }
}
)";

const char* MULTIPLE_KERNEL_SHADER = R"(
struct Params { n: u32, m: u32, r: u32, s: f32 };
@group(0) @binding(0) var<uniform> p: Params;
@group(0) @binding(1) var<storage, read> Kgrams: array<f32>;
@group(0) @binding(2) var<storage, read> weights: array<f32>;
@group(0) @binding(3) var<storage, read> initial_objs: array<f32>;
@group(0) @binding(4) var<storage, read> max_reds: array<f32>;
@group(0) @binding(5) var<storage, read> W: array<f32>;
@group(0) @binding(6) var<storage, read_write> out_res: array<f32>;
var<workgroup> partial: array<f32, 256>;
var<workgroup> sw: array<f32, 256>;
const LOG10_E: f32 = 0.4342944819;
@compute @workgroup_size(256)
fn main(@builtin(workgroup_id) wid: vec3<u32>, @builtin(local_invocation_id) lid: vec3<u32>) {
    let di = wid.x; if (di >= p.r) { return; }
    let tid = lid.x; let n = p.n; let r = p.r; var aggregate = 0.0;
    for (var km = 0u; km < p.m; km++) {
        let offset = km * n * n; var thread_sum = 0.0;
        for (var jt = 0u; jt < n; jt += 256u) {
            let j = jt + tid;
            sw[tid] = select(0.0, W[j * r + di], j < n);
            workgroupBarrier();
            let tile_end = min(256u, n - jt);
            for (var i = tid; i < n; i += 256u) {
                let wi = W[i * r + di];
                for (var k = 0u; k < tile_end; k++) { thread_sum += wi * Kgrams[offset + (jt + k) * n + i] * sw[k]; }
            }
            workgroupBarrier();
        }
        partial[tid] = thread_sum; workgroupBarrier();
        for (var stride = 128u; stride > 0u; stride >>= 1u) {
            if (tid < stride) { partial[tid] += partial[tid + stride]; }
            workgroupBarrier();
        }
        if (tid == 0u) {
            var val = partial[0];
            if (p.s != 0.0) { val = log(max(initial_objs[km] / val, 1e-10)) * LOG10_E / max_reds[km]; }
            aggregate += weights[km] * val;
        }
        workgroupBarrier();
    }
    if (tid == 0u) { out_res[di] = aggregate; }
}
)";

const char* DISTANCE_TILE_SHADER = R"(
struct Params { n: u32, p: u32 };
@group(0) @binding(0) var<uniform> params: Params;
@group(0) @binding(1) var<storage, read> X: array<f32>;
@group(0) @binding(2) var<storage, read_write> D: array<f32>;
var<workgroup> sXi: array<f32, 256>;
var<workgroup> sXj: array<f32, 256>;
@compute @workgroup_size(16, 16)
fn main(@builtin(workgroup_id) wid: vec3<u32>, @builtin(local_invocation_id) lid: vec3<u32>) {
    let tx = lid.x; let ty = lid.y;
    let i = wid.x * 16u + tx; let j = wid.y * 16u + ty;
    let n = params.n; let p = params.p;
    var d = 0.0;
    for (var kt = 0u; kt < p; kt += 16u) {
        let feat = kt + ty;
        let pi = wid.x * 16u + tx; let pj = wid.y * 16u + tx;
        sXi[ty * 16u + tx] = select(0.0, X[feat * n + pi], feat < p && pi < n);
        sXj[ty * 16u + tx] = select(0.0, X[feat * n + pj], feat < p && pj < n);
        workgroupBarrier();
        if (i < n && j < n && i < j) {
            let tile_end = min(16u, p - kt);
            for (var f = 0u; f < tile_end; f++) {
                let diff = sXi[f * 16u + tx] - sXj[f * 16u + ty]; d += diff * diff;
            }
        }
        workgroupBarrier();
    }
    if (i < n && j < n && i <= j) {
        if (i == j) { D[i * n + i] = 0.0; } else { D[j * n + i] = d; D[i * n + j] = d; }
    }
}
)";

Rcpp::NumericMatrix distance_matrix(Rcpp::NumericMatrix X, int device_id) {
    const uint32_t n = X.nrow(), p = X.ncol();
    auto& ctx = get_context(); WGPUInstance& instance = ctx.instance; WGPUDevice& device = ctx.device; WGPUQueue& queue = ctx.queue;
    WGPUShaderSourceWGSL wgsl = {}; wgsl.chain.sType = WGPUSType_ShaderSourceWGSL; wgsl.code.data = DISTANCE_TILE_SHADER; wgsl.code.length = strlen(DISTANCE_TILE_SHADER);
    WGPUShaderModuleDescriptor smd = {}; smd.nextInChain = (WGPUChainedStruct*)&wgsl;
    WGPUShaderModule shader = wgpuDeviceCreateShaderModule(device, &smd);
    uint64_t xSize = (uint64_t)n * p * 4, dSize = (uint64_t)n * n * 4;
    WGPUBufferDescriptor bd_x = {.usage=WGPUBufferUsage_Storage|WGPUBufferUsage_CopyDst, .size=xSize}; WGPUBuffer xB = wgpuDeviceCreateBuffer(device, &bd_x);
    std::vector<float> x_f(n*p); for(size_t i=0; i<n*p; i++) x_f[i]=(float)REAL(X)[i]; wgpuQueueWriteBuffer(queue, xB, 0, x_f.data(), xSize);
    WGPUBufferDescriptor bd_d = {.usage=WGPUBufferUsage_Storage|WGPUBufferUsage_CopySrc, .size=dSize}; WGPUBuffer dB = wgpuDeviceCreateBuffer(device, &bd_d);
    WGPUBufferDescriptor bd_rb = {.usage=WGPUBufferUsage_MapRead|WGPUBufferUsage_CopyDst, .size=dSize}; WGPUBuffer rbB = wgpuDeviceCreateBuffer(device, &bd_rb);
    WGPUBufferDescriptor bd_p = {.usage=WGPUBufferUsage_Uniform|WGPUBufferUsage_CopyDst, .size=8}; WGPUBuffer pB = wgpuDeviceCreateBuffer(device, &bd_p);
    uint32_t pd[]={n,p}; wgpuQueueWriteBuffer(queue, pB, 0, pd, 8);
    WGPUBindGroupLayoutEntry be[3] = {}; for(int i=0; i<3; i++){ be[i].binding=i; be[i].visibility=WGPUShaderStage_Compute; be[i].buffer.type=(i==0?WGPUBufferBindingType_Uniform:(i==1?WGPUBufferBindingType_ReadOnlyStorage:WGPUBufferBindingType_Storage)); }
    WGPUBindGroupLayoutDescriptor bgld = {.entryCount=3, .entries=be}; WGPUBindGroupLayout bgl = wgpuDeviceCreateBindGroupLayout(device, &bgld);
    WGPUBindGroupEntry bge[3] = {}; for(int i=0; i<3; i++){ bge[i].binding=i; bge[i].size=(i==0?8:(i==1?xSize:dSize)); }
    bge[0].buffer=pB; bge[1].buffer=xB; bge[2].buffer=dB;
    WGPUBindGroupDescriptor bgd = {.layout=bgl, .entryCount=3, .entries=bge}; WGPUBindGroup bg = wgpuDeviceCreateBindGroup(device, &bgd);
    WGPUPipelineLayoutDescriptor pld = {.bindGroupLayoutCount=1, .bindGroupLayouts=&bgl}; WGPUPipelineLayout pll = wgpuDeviceCreatePipelineLayout(device, &pld);
    WGPUComputePipelineDescriptor cpd = {.layout=pll, .compute={.module=shader, .entryPoint={.data="main",.length=4}}};
    WGPUComputePipeline pipe = wgpuDeviceCreateComputePipeline(device, &cpd);
    WGPUCommandEncoder enc = wgpuDeviceCreateCommandEncoder(device, NULL); WGPUComputePassEncoder pass = wgpuCommandEncoderBeginComputePass(enc, NULL);
    wgpuComputePassEncoderSetPipeline(pass, pipe); wgpuComputePassEncoderSetBindGroup(pass, 0, bg, 0, NULL);
    wgpuComputePassEncoderDispatchWorkgroups(pass, (n+15)/16, (n+15)/16, 1); wgpuComputePassEncoderEnd(pass);
    wgpuCommandEncoderCopyBufferToBuffer(enc, dB, 0, rbB, 0, dSize); WGPUCommandBuffer cmds = wgpuCommandEncoderFinish(enc, NULL); wgpuQueueSubmit(queue, 1, &cmds);
    bool mapped=false; WGPUBufferMapCallbackInfo mcb = {.callback=handle_buffer_map, .userdata1=&mapped};
    wgpuBufferMapAsync(rbB, WGPUMapMode_Read, 0, dSize, mcb); while(!mapped) wgpuInstanceProcessEvents(instance);
    const float* ptr = (const float*)wgpuBufferGetConstMappedRange(rbB, 0, dSize);
    Rcpp::NumericMatrix res(n, n); for(uint32_t i=0; i<n*n; i++) res[i]=(double)ptr[i];
    wgpuBufferUnmap(rbB); wgpuBufferRelease(xB); wgpuBufferRelease(dB); wgpuBufferRelease(rbB); wgpuBufferRelease(pB);
    wgpuBindGroupRelease(bg); wgpuBindGroupLayoutRelease(bgl); wgpuComputePipelineRelease(pipe); wgpuShaderModuleRelease(shader);
    return res;
}

const char* KERNEL_TILE_SHADER = R"(
struct Params { n: u32, p: u32, kernel_type: u32, poly_s: u32, gamma: f32 };
@group(0) @binding(0) var<uniform> params: Params;
@group(0) @binding(1) var<storage, read> X: array<f32>;
@group(0) @binding(2) var<storage, read_write> K: array<f32>;
var<workgroup> sXi: array<f32, 256>;
var<workgroup> sXj: array<f32, 256>;
@compute @workgroup_size(16, 16)
fn main(@builtin(workgroup_id) wid: vec3<u32>, @builtin(local_invocation_id) lid: vec3<u32>) {
    let tx = lid.x; let ty = lid.y;
    let i = wid.x * 16u + tx; let j = wid.y * 16u + ty;
    let n = params.n; let p = params.p;
    var acc = 0.0;
    for (var kt = 0u; kt < p; kt += 16u) {
        let feat = kt + ty;
        let pi = wid.x * 16u + tx; let pj = wid.y * 16u + tx;
        sXi[ty * 16u + tx] = select(0.0, X[feat * n + pi], feat < p && pi < n);
        sXj[ty * 16u + tx] = select(0.0, X[feat * n + pj], feat < p && pj < n);
        workgroupBarrier();
        if (i < n && j < n && i <= j) {
            let tile_end = min(16u, p - kt);
            if (params.kernel_type == 1u || params.kernel_type == 2u) {
                for (var f = 0u; f < tile_end; f++) { acc += sXi[f * 16u + tx] * sXj[f * 16u + ty]; }
            } else {
                for (var f = 0u; f < tile_end; f++) { let d = sXi[f * 16u + tx] - sXj[f * 16u + ty]; acc += d * d; }
            }
        }
        workgroupBarrier();
    }
    if (i < n && j < n && i <= j) {
        var val = 0.0;
        if (params.kernel_type == 1u) { 
            var base = 1.0 + acc / f32(params.poly_s);
            val = 1.0;
            for (var k = 0u; k < params.poly_s; k++) { val *= base; }
        }
        else if (params.kernel_type == 2u) { val = exp(params.gamma * acc); }
        else if (params.kernel_type == 3u) { val = exp(-params.gamma * sqrt(acc)); }
        else if (params.kernel_type == 4u) { val = 1.0 / sqrt(params.gamma * acc + 1.0); }
        else { val = exp(-params.gamma * acc); }
        K[i * n + j] = val;
        if (i != j) { K[j * n + i] = val; }
    }
}
)";
// kernel_type: 1=poly  2=exponential  3=laplacian  4=inv_mult_quad  5=gaussian

Rcpp::NumericMatrix kernel_matrix(Rcpp::NumericMatrix X, const std::string& kernel, int poly_s, double gamma, int device_id) {
    const uint32_t n = X.nrow(), p = X.ncol();
    uint32_t kt = 5; // gaussian default
    if (kernel == "poly") kt = 1;
    else if (kernel == "exponential") kt = 2;
    else if (kernel == "laplacian") kt = 3;
    else if (kernel == "inv_mult_quad") kt = 4;
    auto& ctx = get_context(); WGPUInstance& instance = ctx.instance; WGPUDevice& device = ctx.device; WGPUQueue& queue = ctx.queue;
    WGPUShaderSourceWGSL wgsl = {}; wgsl.chain.sType = WGPUSType_ShaderSourceWGSL; wgsl.code.data = KERNEL_TILE_SHADER; wgsl.code.length = strlen(KERNEL_TILE_SHADER);
    WGPUShaderModuleDescriptor smd = {}; smd.nextInChain = (WGPUChainedStruct*)&wgsl;
    WGPUShaderModule shader = wgpuDeviceCreateShaderModule(device, &smd);
    uint64_t xSize = (uint64_t)n * p * 4, kSize = (uint64_t)n * n * 4;
    WGPUBufferDescriptor bd_x = {.usage=WGPUBufferUsage_Storage|WGPUBufferUsage_CopyDst, .size=xSize}; WGPUBuffer xB = wgpuDeviceCreateBuffer(device, &bd_x);
    std::vector<float> x_f(n*p); for(size_t i=0; i<n*p; i++) x_f[i]=(float)REAL(X)[i]; wgpuQueueWriteBuffer(queue, xB, 0, x_f.data(), xSize);
    WGPUBufferDescriptor bd_k = {.usage=WGPUBufferUsage_Storage|WGPUBufferUsage_CopySrc, .size=kSize}; WGPUBuffer kB = wgpuDeviceCreateBuffer(device, &bd_k);
    WGPUBufferDescriptor bd_rb = {.usage=WGPUBufferUsage_MapRead|WGPUBufferUsage_CopyDst, .size=kSize}; WGPUBuffer rbB = wgpuDeviceCreateBuffer(device, &bd_rb);
    struct { uint32_t n, p, kt, ps; float g; } pd = {n, p, kt, (uint32_t)poly_s, (float)gamma};
    WGPUBufferDescriptor bd_p = {.usage=WGPUBufferUsage_Uniform|WGPUBufferUsage_CopyDst, .size=20}; WGPUBuffer pB = wgpuDeviceCreateBuffer(device, &bd_p);
    wgpuQueueWriteBuffer(queue, pB, 0, &pd, 20);
    WGPUBindGroupLayoutEntry be[3] = {}; for(int i=0; i<3; i++){ be[i].binding=i; be[i].visibility=WGPUShaderStage_Compute; be[i].buffer.type=(i==0?WGPUBufferBindingType_Uniform:(i==1?WGPUBufferBindingType_ReadOnlyStorage:WGPUBufferBindingType_Storage)); }
    WGPUBindGroupLayoutDescriptor bgld = {.entryCount=3, .entries=be}; WGPUBindGroupLayout bgl = wgpuDeviceCreateBindGroupLayout(device, &bgld);
    WGPUBindGroupEntry bge[3] = {}; for(int i=0; i<3; i++){ bge[i].binding=i; bge[i].size=(i==0?20:(i==1?xSize:kSize)); }
    bge[0].buffer=pB; bge[1].buffer=xB; bge[2].buffer=kB;
    WGPUBindGroupDescriptor bgd = {.layout=bgl, .entryCount=3, .entries=bge}; WGPUBindGroup bg = wgpuDeviceCreateBindGroup(device, &bgd);
    WGPUPipelineLayoutDescriptor pld = {.bindGroupLayoutCount=1, .bindGroupLayouts=&bgl}; WGPUPipelineLayout pll = wgpuDeviceCreatePipelineLayout(device, &pld);
    WGPUComputePipelineDescriptor cpd = {.layout=pll, .compute={.module=shader, .entryPoint={.data="main",.length=4}}};
    WGPUComputePipeline pipe = wgpuDeviceCreateComputePipeline(device, &cpd);
    WGPUCommandEncoder enc = wgpuDeviceCreateCommandEncoder(device, NULL); WGPUComputePassEncoder pass = wgpuCommandEncoderBeginComputePass(enc, NULL);
    wgpuComputePassEncoderSetPipeline(pass, pipe); wgpuComputePassEncoderSetBindGroup(pass, 0, bg, 0, NULL);
    wgpuComputePassEncoderDispatchWorkgroups(pass, (n+15)/16, (n+15)/16, 1); wgpuComputePassEncoderEnd(pass);
    wgpuCommandEncoderCopyBufferToBuffer(enc, kB, 0, rbB, 0, kSize); WGPUCommandBuffer cmds = wgpuCommandEncoderFinish(enc, NULL); wgpuQueueSubmit(queue, 1, &cmds);
    bool mapped=false; WGPUBufferMapCallbackInfo mcb = {.callback=handle_buffer_map, .userdata1=&mapped};
    wgpuBufferMapAsync(rbB, WGPUMapMode_Read, 0, kSize, mcb); while(!mapped) wgpuInstanceProcessEvents(instance);
    const float* ptr = (const float*)wgpuBufferGetConstMappedRange(rbB, 0, kSize);
    Rcpp::NumericMatrix res(n, n); for(uint32_t i=0; i<n*n; i++) res[i]=(double)ptr[i];
    wgpuBufferUnmap(rbB); wgpuBufferRelease(xB); wgpuBufferRelease(kB); wgpuBufferRelease(rbB); wgpuBufferRelease(pB);
    wgpuBindGroupRelease(bg); wgpuBindGroupLayoutRelease(bgl); wgpuComputePipelineRelease(pipe); wgpuShaderModuleRelease(shader);
    return res;
}
Rcpp::NumericVector objective_vals(Rcpp::NumericMatrix W, Rcpp::NumericMatrix Kgram, int device_id) {
    const uint32_t r = W.nrow(), n = W.ncol();
    auto& ctx = get_context(); WGPUInstance& instance = ctx.instance; WGPUDevice& device = ctx.device; WGPUQueue& queue = ctx.queue;
    WGPUShaderSourceWGSL wgsl = {}; wgsl.chain.sType = WGPUSType_ShaderSourceWGSL; wgsl.code.data = BATCH_OBJECTIVE_SHADER; wgsl.code.length = strlen(BATCH_OBJECTIVE_SHADER);
    WGPUShaderModuleDescriptor smd = {}; smd.nextInChain = (WGPUChainedStruct*)&wgsl;
    WGPUShaderModule shader = wgpuDeviceCreateShaderModule(device, &smd);
    uint64_t kS=(uint64_t)n*n*4, wS=(uint64_t)r*n*4, rS=(uint64_t)r*4;
    WGPUBufferDescriptor bd_k = {.usage=WGPUBufferUsage_Storage|WGPUBufferUsage_CopyDst,.size=kS}; WGPUBuffer kB=wgpuDeviceCreateBuffer(device,&bd_k);
    std::vector<float> k_f(n*n); for(size_t i=0; i<n*n; i++) k_f[i]=(float)Kgram[i]; wgpuQueueWriteBuffer(queue, kB, 0, k_f.data(), kS);
    WGPUBufferDescriptor bd_w = {.usage=WGPUBufferUsage_Storage|WGPUBufferUsage_CopyDst,.size=wS}; WGPUBuffer wB=wgpuDeviceCreateBuffer(device,&bd_w);
    std::vector<float> w_f(r*n); for(size_t i=0; i<r*n; i++) w_f[i]=(float)REAL(W)[i]; wgpuQueueWriteBuffer(queue, wB, 0, w_f.data(), wS);
    WGPUBufferDescriptor bd_re = {.usage=WGPUBufferUsage_Storage|WGPUBufferUsage_CopySrc,.size=rS}; WGPUBuffer reB=wgpuDeviceCreateBuffer(device,&bd_re);
    WGPUBufferDescriptor bd_rb = {.usage=WGPUBufferUsage_MapRead|WGPUBufferUsage_CopyDst,.size=rS}; WGPUBuffer rbB=wgpuDeviceCreateBuffer(device,&bd_rb);
    WGPUBufferDescriptor bd_p = {.usage=WGPUBufferUsage_Uniform|WGPUBufferUsage_CopyDst,.size=8}; WGPUBuffer pB=wgpuDeviceCreateBuffer(device,&bd_p);
    uint32_t pd[]={n,r}; wgpuQueueWriteBuffer(queue, pB, 0, pd, 8);
    WGPUBindGroupLayoutEntry be[4]={}; for(int i=0; i<4; i++){ be[i].binding=i; be[i].visibility=WGPUShaderStage_Compute; be[i].buffer.type=(i==0?WGPUBufferBindingType_Uniform:(i==3?WGPUBufferBindingType_Storage:WGPUBufferBindingType_ReadOnlyStorage)); }
    WGPUBindGroupLayoutDescriptor bgld = {.entryCount=4, .entries=be}; WGPUBindGroupLayout bgl=wgpuDeviceCreateBindGroupLayout(device, &bgld);
    WGPUBindGroupEntry bge[4]={}; for(int i=0; i<4; i++){ bge[i].binding=i; bge[i].size=(i==0?8:(i==1?kS:(i==2?wS:rS))); }
    bge[0].buffer=pB; bge[1].buffer=kB; bge[2].buffer=wB; bge[3].buffer=reB;
    WGPUBindGroupDescriptor bgd={.layout=bgl,.entryCount=4,.entries=bge}; WGPUBindGroup bg=wgpuDeviceCreateBindGroup(device,&bgd);
    WGPUPipelineLayoutDescriptor pld = {.bindGroupLayoutCount=1, .bindGroupLayouts=&bgl}; WGPUPipelineLayout pll = wgpuDeviceCreatePipelineLayout(device, &pld);
    WGPUComputePipelineDescriptor cpd = {.layout=pll, .compute={.module=shader, .entryPoint={.data="main",.length=4}}};
    WGPUComputePipeline pipe = wgpuDeviceCreateComputePipeline(device, &cpd);
    WGPUCommandEncoder enc=wgpuDeviceCreateCommandEncoder(device,NULL); WGPUComputePassEncoder pass=wgpuCommandEncoderBeginComputePass(enc,NULL);
    wgpuComputePassEncoderSetPipeline(pass,pipe); wgpuComputePassEncoderSetBindGroup(pass,0,bg,0,NULL);
    wgpuComputePassEncoderDispatchWorkgroups(pass,r,1,1); wgpuComputePassEncoderEnd(pass);
    wgpuCommandEncoderCopyBufferToBuffer(enc,reB,0,rbB,0,rS); WGPUCommandBuffer cmds=wgpuCommandEncoderFinish(enc,NULL); wgpuQueueSubmit(queue,1,&cmds);
    bool mapped=false; WGPUBufferMapCallbackInfo mcb={.callback=handle_buffer_map,.userdata1=&mapped};
    wgpuBufferMapAsync(rbB,WGPUMapMode_Read,0,rS,mcb); while(!mapped) wgpuInstanceProcessEvents(instance);
    const float* ptr=(const float*)wgpuBufferGetConstMappedRange(rbB,0,rS);
    Rcpp::NumericVector rrres(r); for(uint32_t i=0; i<r; i++) rrres[i]=(double)ptr[i];
    wgpuBufferUnmap(rbB); wgpuBufferRelease(kB); wgpuBufferRelease(wB); wgpuBufferRelease(reB); wgpuBufferRelease(rbB); wgpuBufferRelease(pB);
    wgpuBindGroupRelease(bg); wgpuBindGroupLayoutRelease(bgl); wgpuComputePipelineRelease(pipe); wgpuShaderModuleRelease(shader);
    return rrres;
}

Rcpp::NumericVector multiple_kernel_objective_vals(Rcpp::NumericMatrix W, Rcpp::List Kgrams, Rcpp::NumericVector weights, Rcpp::NumericVector initial_objs, Rcpp::NumericVector running_sums, Rcpp::NumericVector max_reds, double maximum_gain_scaling, int device_id) {
    const uint32_t r = W.nrow(), n = W.ncol(), m = Kgrams.size();
    auto& ctx = get_context(); WGPUInstance& instance = ctx.instance; WGPUDevice& device = ctx.device; WGPUQueue& queue = ctx.queue;
    WGPUShaderSourceWGSL wgsl = {}; wgsl.chain.sType = WGPUSType_ShaderSourceWGSL; wgsl.code.data = MULTIPLE_KERNEL_SHADER; wgsl.code.length = strlen(MULTIPLE_KERNEL_SHADER);
    WGPUShaderModuleDescriptor smd = {}; smd.nextInChain = (WGPUChainedStruct*)&wgsl;
    WGPUShaderModule shader = wgpuDeviceCreateShaderModule(device, &smd);
    std::vector<float> k_f(m*n*n); for(uint32_t k=0; k<m; k++){Rcpp::NumericMatrix Kg=Kgrams[k];for(uint32_t i=0; i<n*n; i++) k_f[k*n*n+i]=(float)Kg[i];}
    std::vector<float> w_f(r*n); for(size_t i=0; i<r*n; i++) w_f[i]=(float)REAL(W)[i];
    std::vector<float> wts_f(m); for(uint32_t k=0; k<m; k++) wts_f[k]=(float)weights[k];
    std::vector<float> in_f(m); for(uint32_t k=0; k<m; k++) in_f[k]=(float)initial_objs[k];
    std::vector<float> mr_f(m); for(uint32_t k=0; k<m; k++) mr_f[k]=(float)max_reds[k];
    uint64_t kS=k_f.size()*4, wS=w_f.size()*4, mS=m*4, rS=r*4;
    WGPUBufferDescriptor bd_k={.usage=WGPUBufferUsage_Storage|WGPUBufferUsage_CopyDst,.size=kS}; WGPUBuffer kB=wgpuDeviceCreateBuffer(device,&bd_k); wgpuQueueWriteBuffer(queue,kB,0,k_f.data(),kS);
    WGPUBufferDescriptor bd_w={.usage=WGPUBufferUsage_Storage|WGPUBufferUsage_CopyDst,.size=wS}; WGPUBuffer wB=wgpuDeviceCreateBuffer(device,&bd_w); wgpuQueueWriteBuffer(queue,wB,0,w_f.data(),wS);
    WGPUBufferDescriptor bd_m={.usage=WGPUBufferUsage_Storage|WGPUBufferUsage_CopyDst,.size=mS};
    WGPUBuffer wtB=wgpuDeviceCreateBuffer(device,&bd_m); wgpuQueueWriteBuffer(queue,wtB,0,wts_f.data(),mS);
    WGPUBuffer inB=wgpuDeviceCreateBuffer(device,&bd_m); wgpuQueueWriteBuffer(queue,inB,0,in_f.data(),mS);
    WGPUBuffer maB=wgpuDeviceCreateBuffer(device,&bd_m); wgpuQueueWriteBuffer(queue,maB,0,mr_f.data(),mS);
    WGPUBufferDescriptor bd_r={.usage=WGPUBufferUsage_Storage|WGPUBufferUsage_CopySrc,.size=rS}; WGPUBuffer reB=wgpuDeviceCreateBuffer(device,&bd_r);
    WGPUBufferDescriptor bd_rb={.usage=WGPUBufferUsage_MapRead|WGPUBufferUsage_CopyDst,.size=rS}; WGPUBuffer rbB=wgpuDeviceCreateBuffer(device,&bd_rb);
    WGPUBufferDescriptor bd_p={.usage=WGPUBufferUsage_Uniform|WGPUBufferUsage_CopyDst,.size=16}; WGPUBuffer pB=wgpuDeviceCreateBuffer(device,&bd_p);
    struct { uint32_t n, m, r; float s; } pd = {n, m, r, (float)maximum_gain_scaling}; wgpuQueueWriteBuffer(queue, pB, 0, &pd, 16);
    WGPUBindGroupLayoutEntry be[7]={}; for(int i=0; i<7; i++){ be[i].binding=i; be[i].visibility=WGPUShaderStage_Compute; be[i].buffer.type=(i==0?WGPUBufferBindingType_Uniform:(i==6?WGPUBufferBindingType_Storage:WGPUBufferBindingType_ReadOnlyStorage)); }
    WGPUBindGroupLayoutDescriptor bgld = {.entryCount=7,.entries=be}; WGPUBindGroupLayout bgl=wgpuDeviceCreateBindGroupLayout(device, &bgld);
    WGPUBindGroupEntry ge[7]={}; for(int i=0; i<7; i++){ ge[i].binding=i; ge[i].size=(i==0?16:(i==6?rS:(i==1?kS:(i==5?wS:mS)))); }
    ge[0].buffer=pB; ge[1].buffer=kB; ge[2].buffer=wtB; ge[3].buffer=inB; ge[4].buffer=maB; ge[5].buffer=wB; ge[6].buffer=reB;
    WGPUBindGroupDescriptor bgd={.layout=bgl,.entryCount=7,.entries=ge}; WGPUBindGroup bg=wgpuDeviceCreateBindGroup(device,&bgd);
    WGPUPipelineLayoutDescriptor pld={.bindGroupLayoutCount=1,.bindGroupLayouts=&bgl}; WGPUPipelineLayout pll=wgpuDeviceCreatePipelineLayout(device,&pld);
    WGPUComputePipelineDescriptor cpd={.layout=pll,.compute={.module=shader,.entryPoint={.data="main",.length=4}}};
    WGPUComputePipeline pipe=wgpuDeviceCreateComputePipeline(device,&cpd);
    WGPUCommandEncoder enc=wgpuDeviceCreateCommandEncoder(device,NULL); WGPUComputePassEncoder pass=wgpuCommandEncoderBeginComputePass(enc,NULL);
    wgpuComputePassEncoderSetPipeline(pass,pipe); wgpuComputePassEncoderSetBindGroup(pass,0,bg,0,NULL);
    wgpuComputePassEncoderDispatchWorkgroups(pass,r,1,1); wgpuComputePassEncoderEnd(pass);
    wgpuCommandEncoderCopyBufferToBuffer(enc,reB,0,rbB,0,rS); WGPUCommandBuffer cmds=wgpuCommandEncoderFinish(enc,NULL); wgpuQueueSubmit(queue,1,&cmds);
    bool mapped=false; WGPUBufferMapCallbackInfo mcb={.callback=handle_buffer_map,.userdata1=&mapped};
    wgpuBufferMapAsync(rbB,WGPUMapMode_Read,0,rS,mcb); while(!mapped) wgpuInstanceProcessEvents(instance);
    const float* ptr=(const float*)wgpuBufferGetConstMappedRange(rbB,0,rS);
    Rcpp::NumericVector res(r); for(uint32_t i=0; i<r; i++) res[i]=(double)ptr[i];
    wgpuBufferUnmap(rbB); wgpuBufferRelease(kB); wgpuBufferRelease(wB); wgpuBufferRelease(wtB); wgpuBufferRelease(inB); wgpuBufferRelease(maB); wgpuBufferRelease(reB); wgpuBufferRelease(rbB); wgpuBufferRelease(pB);
    wgpuBindGroupRelease(bg); wgpuBindGroupLayoutRelease(bgl); wgpuComputePipelineRelease(pipe); wgpuShaderModuleRelease(shader);
    return res;
}

const char* ARGMIN_SHADER = R"(
struct Params { n: u32 };
@group(0) @binding(0) var<uniform> p: Params;
@group(0) @binding(1) var<storage, read> vals: array<f32>;
@group(0) @binding(2) var<storage, read_write> out_idx: array<u32>;
var<workgroup> ws_val: array<f32, 256>;
var<workgroup> ws_idx: array<u32, 256>;
@compute @workgroup_size(256)
fn main(@builtin(local_invocation_id) lid: vec3<u32>) {
    let tid = lid.x; var my_val = 1e30f; var my_idx = 0u;
    for (var i = tid; i < p.n; i += 256u) {
        let v = vals[i]; if (v < my_val) { my_val = v; my_idx = i; }
    }
    ws_val[tid] = my_val; ws_idx[tid] = my_idx; workgroupBarrier();
    for (var stride = 128u; stride > 0u; stride >>= 1u) {
        if (tid < stride) {
            if (ws_val[tid + stride] < ws_val[tid]) { ws_val[tid] = ws_val[tid + stride]; ws_idx[tid] = ws_idx[tid + stride]; }
        }
        workgroupBarrier();
    }
    if (tid == 0u) { out_idx[0] = ws_idx[0]; }
}
)";

Rcpp::IntegerVector full_greedy_search(Rcpp::NumericMatrix X, Rcpp::NumericMatrix Sinv, Rcpp::IntegerVector start_indicT, int max_iters, int device_id) {
    const uint32_t n = X.nrow(), p = X.ncol();
    const uint32_t nT = (uint32_t)std::accumulate(start_indicT.begin(), start_indicT.end(), 0);
    const uint32_t nC = n - nT;
    auto& ctx = get_context(); WGPUInstance& instance = ctx.instance; WGPUDevice& device = ctx.device; WGPUQueue& queue = ctx.queue;

    std::vector<float> x_f(n*p); for(size_t i=0; i<n*p; i++) x_f[i]=(float)REAL(X)[i];
    std::vector<float> s_f(p*p); for(size_t i=0; i<p*p; i++) s_f[i]=(float)REAL(Sinv)[i];
    std::vector<uint32_t> ci(n); for(uint32_t i=0; i<n; i++) ci[i]=(uint32_t)start_indicT[i];
    std::vector<uint32_t> i_Ts, i_Cs; for(uint32_t i=0; i<n; i++) { if(ci[i]) i_Ts.push_back(i); else i_Cs.push_back(i); }
    uint32_t numIT=i_Ts.size(), numIC=i_Cs.size(); uint64_t resS=(uint64_t)numIT*numIC*4;

    // initial avg (incremental updates after each swap avoid O(n*p) recompute per iteration)
    std::vector<float> at(p,0), ac(p,0);
    for(uint32_t i=0; i<n; i++){ if(ci[i]) for(uint32_t j=0; j<p; j++) at[j]+=x_f[j*n+i]; else for(uint32_t j=0; j<p; j++) ac[j]+=x_f[j*n+i]; }
    for(uint32_t j=0; j<p; j++){ at[j]/=nT; ac[j]/=nC; }

    auto mkBuf = [&](WGPUBufferUsage usage, uint64_t size) -> WGPUBuffer {
        WGPUBufferDescriptor bd = {.usage=usage, .size=size}; return wgpuDeviceCreateBuffer(device, &bd);
    };
    WGPUBuffer xB  = mkBuf(WGPUBufferUsage_Storage|WGPUBufferUsage_CopyDst, (uint64_t)n*p*4);
    WGPUBuffer sB  = mkBuf(WGPUBufferUsage_Storage|WGPUBufferUsage_CopyDst, (uint64_t)p*p*4);
    WGPUBuffer itB = mkBuf(WGPUBufferUsage_Storage|WGPUBufferUsage_CopyDst, (uint64_t)numIT*4);
    WGPUBuffer icB = mkBuf(WGPUBufferUsage_Storage|WGPUBufferUsage_CopyDst, (uint64_t)numIC*4);
    WGPUBuffer atB = mkBuf(WGPUBufferUsage_Storage|WGPUBufferUsage_CopyDst, (uint64_t)p*4);
    WGPUBuffer acB = mkBuf(WGPUBufferUsage_Storage|WGPUBufferUsage_CopyDst, (uint64_t)p*4);
    WGPUBuffer rB  = mkBuf(WGPUBufferUsage_Storage|WGPUBufferUsage_CopySrc, resS);
    WGPUBuffer prB = mkBuf(WGPUBufferUsage_Uniform|WGPUBufferUsage_CopyDst, 32);
    WGPUBuffer apB = mkBuf(WGPUBufferUsage_Uniform|WGPUBufferUsage_CopyDst, 16);  // argmin params
    WGPUBuffer arB = mkBuf(WGPUBufferUsage_Storage|WGPUBufferUsage_CopySrc, 4);   // argmin result index
    WGPUBuffer arbB= mkBuf(WGPUBufferUsage_MapRead|WGPUBufferUsage_CopyDst, 4);   // readback (4 bytes vs resS)

    wgpuQueueWriteBuffer(queue,xB,0,x_f.data(),(uint64_t)n*p*4);
    wgpuQueueWriteBuffer(queue,sB,0,s_f.data(),(uint64_t)p*p*4);
    wgpuQueueWriteBuffer(queue,itB,0,i_Ts.data(),(uint64_t)numIT*4);
    wgpuQueueWriteBuffer(queue,icB,0,i_Cs.data(),(uint64_t)numIC*4);
    wgpuQueueWriteBuffer(queue,atB,0,at.data(),(uint64_t)p*4);
    wgpuQueueWriteBuffer(queue,acB,0,ac.data(),(uint64_t)p*4);
    uint32_t n_pairs=numIT*numIC; wgpuQueueWriteBuffer(queue,apB,0,&n_pairs,4);

    auto mkShader = [&](const char* src) -> WGPUShaderModule {
        WGPUShaderSourceWGSL wgsl={}; wgsl.chain.sType=WGPUSType_ShaderSourceWGSL; wgsl.code.data=src; wgsl.code.length=strlen(src);
        WGPUShaderModuleDescriptor smd={}; smd.nextInChain=(WGPUChainedStruct*)&wgsl; return wgpuDeviceCreateShaderModule(device,&smd);
    };
    WGPUShaderModule obj_shader  = mkShader(OBJECTIVE_SHADER);
    WGPUShaderModule am_shader   = mkShader(ARGMIN_SHADER);

    // objective pipeline (8 bindings)
    WGPUBindGroupLayoutEntry obj_be[8]={}; for(int i=0; i<8; i++){ obj_be[i].binding=i; obj_be[i].visibility=WGPUShaderStage_Compute; obj_be[i].buffer.type=(i==0?WGPUBufferBindingType_Uniform:(i==7?WGPUBufferBindingType_Storage:WGPUBufferBindingType_ReadOnlyStorage)); }
    WGPUBindGroupLayoutDescriptor obj_bgld={.entryCount=8,.entries=obj_be}; WGPUBindGroupLayout obj_bgl=wgpuDeviceCreateBindGroupLayout(device,&obj_bgld);
    WGPUBindGroupEntry obj_bge[8]={}; for(int i=0; i<8; i++){ obj_bge[i].binding=i; obj_bge[i].size=(i==0?32:(i==1?(uint64_t)n*p*4:(i==2?(uint64_t)p*p*4:(i==3?(uint64_t)numIT*4:(i==4?(uint64_t)numIC*4:(i==5?(uint64_t)p*4:(i==6?(uint64_t)p*4:resS))))))); }
    obj_bge[0].buffer=prB; obj_bge[1].buffer=xB; obj_bge[2].buffer=sB; obj_bge[3].buffer=itB; obj_bge[4].buffer=icB; obj_bge[5].buffer=atB; obj_bge[6].buffer=acB; obj_bge[7].buffer=rB;
    WGPUBindGroupDescriptor obj_bgd={.layout=obj_bgl,.entryCount=8,.entries=obj_bge}; WGPUBindGroup obj_bg=wgpuDeviceCreateBindGroup(device,&obj_bgd);
    WGPUPipelineLayoutDescriptor obj_pld={.bindGroupLayoutCount=1,.bindGroupLayouts=&obj_bgl}; WGPUPipelineLayout obj_pll=wgpuDeviceCreatePipelineLayout(device,&obj_pld);
    WGPUComputePipelineDescriptor obj_cpd={.layout=obj_pll,.compute={.module=obj_shader,.entryPoint={.data="main",.length=4}}};
    WGPUComputePipeline obj_pipe=wgpuDeviceCreateComputePipeline(device,&obj_cpd);

    // argmin pipeline (3 bindings: uniform params, read-only results, write-only index output)
    WGPUBindGroupLayoutEntry am_be[3]={}; am_be[0].binding=0; am_be[0].visibility=WGPUShaderStage_Compute; am_be[0].buffer.type=WGPUBufferBindingType_Uniform; am_be[1].binding=1; am_be[1].visibility=WGPUShaderStage_Compute; am_be[1].buffer.type=WGPUBufferBindingType_ReadOnlyStorage; am_be[2].binding=2; am_be[2].visibility=WGPUShaderStage_Compute; am_be[2].buffer.type=WGPUBufferBindingType_Storage;
    WGPUBindGroupLayoutDescriptor am_bgld={.entryCount=3,.entries=am_be}; WGPUBindGroupLayout am_bgl=wgpuDeviceCreateBindGroupLayout(device,&am_bgld);
    WGPUBindGroupEntry am_bge[3]={}; am_bge[0].binding=0; am_bge[0].buffer=apB; am_bge[0].size=16; am_bge[1].binding=1; am_bge[1].buffer=rB; am_bge[1].size=resS; am_bge[2].binding=2; am_bge[2].buffer=arB; am_bge[2].size=4;
    WGPUBindGroupDescriptor am_bgd={.layout=am_bgl,.entryCount=3,.entries=am_bge}; WGPUBindGroup am_bg=wgpuDeviceCreateBindGroup(device,&am_bgd);
    WGPUPipelineLayoutDescriptor am_pld={.bindGroupLayoutCount=1,.bindGroupLayouts=&am_bgl}; WGPUPipelineLayout am_pll=wgpuDeviceCreatePipelineLayout(device,&am_pld);
    WGPUComputePipelineDescriptor am_cpd={.layout=am_pll,.compute={.module=am_shader,.entryPoint={.data="main",.length=4}}};
    WGPUComputePipeline am_pipe=wgpuDeviceCreateComputePipeline(device,&am_cpd);

    for(int iter=0; iter<max_iters; iter++) {
        uint32_t pd[]={n,p,nT,numIT,numIC}; wgpuQueueWriteBuffer(queue,prB,0,pd,20);
        wgpuQueueWriteBuffer(queue,atB,0,at.data(),(uint64_t)p*4);
        wgpuQueueWriteBuffer(queue,acB,0,ac.data(),(uint64_t)p*4);
        WGPUCommandEncoder enc=wgpuDeviceCreateCommandEncoder(device,NULL); WGPUComputePassEncoder pass=wgpuCommandEncoderBeginComputePass(enc,NULL);
        wgpuComputePassEncoderSetPipeline(pass,obj_pipe); wgpuComputePassEncoderSetBindGroup(pass,0,obj_bg,0,NULL);
        wgpuComputePassEncoderDispatchWorkgroups(pass,(numIT*numIC+63)/64,1,1);
        wgpuComputePassEncoderSetPipeline(pass,am_pipe); wgpuComputePassEncoderSetBindGroup(pass,0,am_bg,0,NULL);
        wgpuComputePassEncoderDispatchWorkgroups(pass,1,1,1);
        wgpuComputePassEncoderEnd(pass);
        wgpuCommandEncoderCopyBufferToBuffer(enc,arB,0,arbB,0,4); WGPUCommandBuffer cmds=wgpuCommandEncoderFinish(enc,NULL); wgpuQueueSubmit(queue,1,&cmds);
        bool mapped=false; WGPUBufferMapCallbackInfo mcb={.callback=handle_buffer_map,.userdata1=&mapped};
        wgpuBufferMapAsync(arbB,WGPUMapMode_Read,0,4,mcb); while(!mapped) wgpuInstanceProcessEvents(instance);
        const uint32_t* ptr=(const uint32_t*)wgpuBufferGetConstMappedRange(arbB,0,4);
        uint32_t best_flat=ptr[0];
        wgpuBufferUnmap(arbB);
        uint32_t bT=i_Ts[best_flat/numIC], bC=i_Cs[best_flat%numIC]; ci[bT]=0; ci[bC]=1;
        i_Ts.clear(); i_Cs.clear(); for(uint32_t i=0; i<n; i++) if(ci[i]) i_Ts.push_back(i); else i_Cs.push_back(i);
        // incremental avg update: O(p) instead of O(n*p)
        for(uint32_t j=0; j<p; j++){
            at[j]=(at[j]*nT - x_f[j*n+bT] + x_f[j*n+bC])/nT;
            ac[j]=(ac[j]*nC - x_f[j*n+bC] + x_f[j*n+bT])/nC;
        }
        wgpuQueueWriteBuffer(queue,itB,0,i_Ts.data(),(uint64_t)numIT*4);
        wgpuQueueWriteBuffer(queue,icB,0,i_Cs.data(),(uint64_t)numIC*4);
    }
    Rcpp::IntegerVector res(n); for(uint32_t i=0; i<n; i++) res[i]=(int)ci[i];
    wgpuBufferRelease(xB); wgpuBufferRelease(sB); wgpuBufferRelease(itB); wgpuBufferRelease(icB); wgpuBufferRelease(atB); wgpuBufferRelease(acB); wgpuBufferRelease(rB); wgpuBufferRelease(prB); wgpuBufferRelease(apB); wgpuBufferRelease(arB); wgpuBufferRelease(arbB);
    wgpuBindGroupRelease(obj_bg); wgpuBindGroupLayoutRelease(obj_bgl); wgpuComputePipelineRelease(obj_pipe); wgpuShaderModuleRelease(obj_shader);
    wgpuBindGroupRelease(am_bg); wgpuBindGroupLayoutRelease(am_bgl); wgpuComputePipelineRelease(am_pipe); wgpuShaderModuleRelease(am_shader);
    wgpuPipelineLayoutRelease(obj_pll); wgpuPipelineLayoutRelease(am_pll);
    return res;
}

} // namespace ged_wgpu_backend
#endif
