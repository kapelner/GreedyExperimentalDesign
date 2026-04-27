struct Params {
    n: u32,
    p: u32,
    n_over_two: u32,
    num_designs: u32,
};

@group(0) @binding(0) var<uniform> params: Params;
@group(0) @binding(1) var<storage, read> X: array<f32>; // n x p, flat
@group(0) @binding(2) var<storage, read> Sinv: array<f32>; // p x p, flat
@group(0) @binding(3) var<storage, read_write> results: array<f32>; // results, length = num_combinations

fn choose(n: u32, k_in: u32) -> u32 {
    if (k_in > n) {
        return 0u;
    }
    var k = k_in;
    if (k > n - k) {
        k = n - k;
    }
    var result = 1u;
    for (var i = 1u; i <= k; i = i + 1u) {
        result = (result * (n - k + i)) / i;
    }
    return result;
}

fn combination_mask(index: u32, n: u32, k: u32) -> u32 {
    var remaining_index = index;
    var remaining_k = k;
    var mask = 0u;

    for (var pos = 0u; pos < n; pos = pos + 1u) {
        if (remaining_k == 0u) {
            break;
        }
        let with_pos = choose(n - pos - 1u, remaining_k - 1u);
        if (remaining_index < with_pos) {
            mask = mask | (1u << pos);
            remaining_k = remaining_k - 1u;
        } else {
            remaining_index = remaining_index - with_pos;
        }
    }

    return mask;
}

fn is_treatment(mask: u32, i: u32) -> bool {
    return ((mask >> i) & 1u) == 1u;
}

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) global_id: vec3<u32>) {
    let design_idx = global_id.x;
    let n = params.n;
    let p = params.p;

    if (design_idx >= params.num_designs) {
        return;
    }
    if (n > 32u || params.n_over_two == 0u || params.n_over_two >= n) {
        results[design_idx] = 3.402823e38;
        return;
    }

    let mask = combination_mask(design_idx, n, params.n_over_two);
    let nT = f32(params.n_over_two);
    let nC = f32(n - params.n_over_two);

    var val = 0.0;
    for (var a = 0u; a < p; a = a + 1u) {
        var sum_T_a = 0.0;
        var sum_C_a = 0.0;
        for (var row_a = 0u; row_a < n; row_a = row_a + 1u) {
            let x = X[row_a * p + a];
            if (is_treatment(mask, row_a)) {
                sum_T_a = sum_T_a + x;
            } else {
                sum_C_a = sum_C_a + x;
            }
        }
        let diff_a = (sum_T_a / nT) - (sum_C_a / nC);

        var inner = 0.0;
        for (var b = 0u; b < p; b = b + 1u) {
            var sum_T_b = 0.0;
            var sum_C_b = 0.0;
            for (var row_b = 0u; row_b < n; row_b = row_b + 1u) {
                let x = X[row_b * p + b];
                if (is_treatment(mask, row_b)) {
                    sum_T_b = sum_T_b + x;
                } else {
                    sum_C_b = sum_C_b + x;
                }
            }
            let diff_b = (sum_T_b / nT) - (sum_C_b / nC);
            inner = inner + Sinv[a * p + b] * diff_b;
        }
        val = val + diff_a * inner;
    }

    results[design_idx] = val;
}
