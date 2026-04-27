struct Params {
    n: u32,
    r: u32,
};

@group(0) @binding(0) var<uniform> params: Params;
@group(0) @binding(1) var<storage, read> Kgram: array<f32>; // n x n, flat
@group(0) @binding(2) var<storage, read> W: array<f32>;     // r x n, row-major (+1 or -1)
@group(0) @binding(3) var<storage, read_write> result: array<f32>; // length r

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) global_id: vec3<u32>) {
    let design_idx = global_id.x;
    let n = params.n;
    if (design_idx >= params.r) {
        return;
    }

    var sum = 0.0;

    for (var i = 0u; i < n; i = i + 1u) {
        let wi = W[design_idx * n + i];
        var inner_sum = 0.0;
        for (var j = 0u; j < n; j = j + 1u) {
            inner_sum = inner_sum + Kgram[i * n + j] * W[design_idx * n + j];
        }
        sum = sum + wi * inner_sum;
    }

    result[design_idx] = sum;
}
