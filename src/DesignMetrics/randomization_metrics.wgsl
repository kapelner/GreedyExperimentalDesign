struct Params {
    n: u32,
    r: u32,
};

@group(0) @binding(0) var<uniform> params: Params;
@group(0) @binding(1) var<storage, read> designs: array<u32>; // r x n, row-major 0/1 assignments
@group(0) @binding(2) var<storage, read_write> results: array<f32>; // n x n matrix, flat

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) global_id: vec3<u32>) {
    let pair_idx = global_id.x;
    let n = params.n;
    let r = params.r;
    let num_pairs = n * (n - 1u) / 2u;

    if (pair_idx >= num_pairs) { return; }

    // Map pair_idx to (i1, i2)
    var i1 = 0u;
    var i2 = 0u;
    var count = 0u;
    var found = false;
    for (var i = 0u; i < n - 1u; i = i + 1u) {
        for (var j = i + 1u; j < n; j = j + 1u) {
            if (count == pair_idx) {
                i1 = i;
                i2 = j;
                found = true;
                break;
            }
            count = count + 1u;
        }
        if (found) { break; }
    }

    var num_same_group = 0u;
    for (var j = 0u; j < r; j = j + 1u) {
        let val1 = designs[j * n + i1];
        let val2 = designs[j * n + i2];
        if (val1 == val2) {
            num_same_group = num_same_group + 1u;
        }
    }

    results[i1 * n + i2] = f32(num_same_group) / f32(r);
}
