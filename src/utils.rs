/// compute the dot product between a vector and a matrix
pub fn dot(v: [f64; 3], m: [[f64; 3]; 3]) -> [f64; 3] {
    let mut out = [0f64; 3];
    for i in 0..3 {
        out[i] = v[0] * m[0][i] + v[1] * m[1][i] + v[2] * m[2][i]
    }
    return out;
}

/// compute the dot product between two vectors
pub fn vdot(a: [f64; 3], b: [f64; 3]) -> f64 {
    let mut out = 0f64;
    for i in 0..3 {
        out += a[i] * b[i]
    }
    return out;
}

/// compute the norm of a vector
pub fn norm(a: [f64; 3]) -> f64 {
    a.iter().map(|a| a.powi(2)).sum::<f64>().powf(0.5)
}

/// compute M.T * M
pub fn transpose_square(m: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    [[vdot([m[0][0], m[1][0], m[2][0]], [m[0][0], m[1][0], m[2][0]]),
      vdot([m[0][0], m[1][0], m[2][0]], [m[0][1], m[1][1], m[2][1]]),
      vdot([m[0][0], m[1][0], m[2][0]], [m[0][2], m[1][2], m[2][2]])],
     [vdot([m[0][1], m[1][1], m[2][1]], [m[0][0], m[1][0], m[2][0]]),
      vdot([m[0][1], m[1][1], m[2][1]], [m[0][1], m[1][1], m[2][1]]),
      vdot([m[0][1], m[1][1], m[2][1]], [m[0][2], m[1][2], m[2][2]])],
     [vdot([m[0][2], m[1][2], m[2][2]], [m[0][0], m[1][0], m[2][0]]),
      vdot([m[0][2], m[1][2], m[2][2]], [m[0][1], m[1][1], m[2][1]]),
      vdot([m[0][2], m[1][2], m[2][2]], [m[0][2], m[1][2], m[2][2]])]]
}
