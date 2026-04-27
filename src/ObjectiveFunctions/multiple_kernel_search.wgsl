struct Params {
    n: u32,
    m: u32,
    num_i_Ts: u32,
    num_i_Cs: u32,
    maximum_gain_scaling: f32,
};

@group(0) @binding(0) var<uniform> params: Params;
@group(0) @binding(1) var<storage, read> Kgrams: array<f32>;        // m * n * n (packed)
@group(0) @binding(2) var<storage, read> weights: array<f32>;       // m
@group(0) @binding(3) var<storage, read> initial_objs: array<f32>;  // m
@group(0) @binding(4) var<storage, read> running_sums: array<f32>;  // m
@group(0) @binding(5) var<storage, read> max_reds: array<f32>;      // m
@group(0) @binding(6) var<storage, read> w: array<f32>;             // n (current design in {1, -1})
@group(0) @binding(7) var<storage, read> i_Ts: array<u32>;          // num_i_Ts
@group(0) @binding(8) var<storage, read> i_Cs: array<u32>;          // num_i_Cs
@group(0) @binding(9) var<storage, read_write> results: array<f32>; // num_i_Ts * num_i_Cs

const LOG10_E: f32 = 0.4342944819; // log10(e)

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) global_id: vec3<u32>) {
    let idx = global_id.x;
    let num_swaps = params.num_i_Ts * params.num_i_Cs;
    if (idx >= num_swaps) { return; }

    let t_idx = idx / params.num_i_Cs;
    let c_idx = idx % params.num_i_Cs;
    
    let i_T = i_Ts[t_idx];
    let i_C = i_Cs[c_idx];

    var aggregate_obj = 0.0;

    for (var k = 0u; k < params.m; k = k + 1u) {
        let n = params.n;
        let k_offset = k * n * n;
        
        // Compute qcd for this kernel: sum_{l != t, c} w[l] * (K[t][l] - K[c][l])
        var qcd = 0.0;
        for (var l = 0u; l < n; l = l + 1u) {
            if (l == i_T || l == i_C) { continue; }
            let k_tl = Kgrams[k_offset + i_T * n + l];
            let k_cl = Kgrams[k_offset + i_C * n + l];
            qcd = qcd + (w[l] * (k_tl - k_cl));
        }

        let current_k = running_sums[k] - 4.0 * qcd;
        
        // pct_off_from_best(k) = 1 - log10(initial_k / current_k) / (max_red_k * scaling)
        // log10(a/b) = log(a/b) * LOG10_E
        let ratio = initial_objs[k] / current_k;
        let log10_ratio = log(max(ratio, 1e-10)) * LOG10_E;
        
        let pct_red = 1.0 - log10_ratio / (max_reds[k] * params.maximum_gain_scaling);
        aggregate_obj = aggregate_obj + (weights[k] * pct_red);
    }
    
    results[idx] = aggregate_obj;
}
