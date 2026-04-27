struct Params {
    n: u32,
    p: u32,
    nT: u32,
    num_i_Ts: u32,
    num_i_Cs: u32,
};

@group(0) @binding(0) var<uniform> params: Params;
@group(0) @binding(1) var<storage, read> X: array<f32>; // n x p, column-major? R is column-major.
@group(0) @binding(2) var<storage, read> Sinv: array<f32>; // p x p
@group(0) @binding(3) var<storage, read> i_Ts: array<u32>;
@group(0) @binding(4) var<storage, read> i_Cs: array<u32>;
@group(0) @binding(5) var<storage, read> current_avg_T: array<f32>; // p
@group(0) @binding(6) var<storage, read> current_avg_C: array<f32>; // p
@group(0) @binding(7) var<storage, read_write> results: array<f32>; // size = num_i_Ts * num_i_Cs

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) global_id: vec3<u32>) {
    let idx = global_id.x;
    let num_swaps = params.num_i_Ts * params.num_i_Cs;
    if (idx >= num_swaps) { return; }

    let t_idx = idx / params.num_i_Cs;
    let c_idx = idx % params.num_i_Cs;
    
    let i_T = i_Ts[t_idx];
    let i_C = i_Cs[c_idx];

    // Compute updated averages
    // avg_T_new = avg_T_old - X[i_T]/nT + X[i_C]/nT
    // avg_C_new = avg_C_old - X[i_C]/nC + X[i_T]/nC
    let nC = params.n - params.nT;

    // Compute Mahal: diff^T * Sinv * diff
    var val = 0.0;
    for (var i = 0u; i < params.p; i = i + 1u) {
        let x_i_T_i = X[i_T * params.p + i];
        let x_i_C_i = X[i_C * params.p + i];
        let avg_T_i = current_avg_T[i] - (x_i_T_i / f32(params.nT)) + (x_i_C_i / f32(params.nT));
        let avg_C_i = current_avg_C[i] - (x_i_C_i / f32(nC)) + (x_i_T_i / f32(nC));
        let diff_i = avg_T_i - avg_C_i;

        var temp = 0.0;
        for (var j = 0u; j < params.p; j = j + 1u) {
            let x_i_T_j = X[i_T * params.p + j];
            let x_i_C_j = X[i_C * params.p + j];
            let avg_T_j = current_avg_T[j] - (x_i_T_j / f32(params.nT)) + (x_i_C_j / f32(params.nT));
            let avg_C_j = current_avg_C[j] - (x_i_C_j / f32(nC)) + (x_i_T_j / f32(nC));
            let diff_j = avg_T_j - avg_C_j;
            temp = temp + Sinv[i * params.p + j] * diff_j;
        }
        val = val + diff_i * temp;
    }
    
    results[idx] = val;
}
