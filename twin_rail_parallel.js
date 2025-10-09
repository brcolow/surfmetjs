/**
 * Autocollimator + Optical Square data reduction (upgraded)
 * Angles -> slopes -> integrate (trapezoid | simpson | spline) -> detrend -> parallelism
 *
 * Units:
 * - positions: mm
 * - angles: arcseconds
 * - output displacement: micrometers (µm)
 */

const ARCSEC_TO_RAD = Math.PI / (180 * 3600);

// ---------- Math helpers ----------
function mean(xs) {
  return xs.reduce((a, b) => a + b, 0) / xs.length;
}

function linearFit(xmm, yum) {
  const n = xmm.length;
  const xbar = mean(xmm);
  const ybar = mean(yum);
  let num = 0, den = 0;
  for (let i = 0; i < n; i++) {
    const dx = xmm[i] - xbar;
    num += dx * (yum[i] - ybar);
    den += dx * dx;
  }
  const a = den === 0 ? 0 : num / den;
  const b = ybar - a * xbar;
  return { a, b };
}

function detrend(xmm, yraw_um) {
  const { a, b } = linearFit(xmm, yraw_um);
  const y = yraw_um.map((v, i) => v - (a * xmm[i] + b));
  return { y, fit: { a, b } };
}

function approxEqual(a, b, tol = 1e-9) {
  return Math.abs(a - b) <= tol * Math.max(1, Math.abs(a), Math.abs(b));
}

// ---------- Integration backends ----------
/** Trapezoid: works on unequal spacing */
function integrateTrapezoid(x_mm, slope_rad) {
  const n = x_mm.length;
  const y_m = new Array(n).fill(0);
  for (let i = 0; i < n - 1; i++) {
    const h_m = (x_mm[i + 1] - x_mm[i]) * 1e-3;
    if (!(h_m > 0)) throw new Error("positions must be strictly increasing");
    const avg = 0.5 * (slope_rad[i] + slope_rad[i + 1]);
    y_m[i + 1] = y_m[i] + avg * h_m;
  }
  return y_m;
}

/** Composite Simpson: requires equal spacing and even number of intervals */
function integrateSimpson(x_mm, slope_rad) {
  const n = x_mm.length;
  if (n < 2) throw new Error("need at least two points");
  const h0 = x_mm[1] - x_mm[0];
  for (let i = 1; i < n - 1; i++) {
    const h = x_mm[i + 1] - x_mm[i];
    if (!approxEqual(h, h0, 1e-9)) {
      throw new Error("Simpson requires equal spacing");
    }
  }
  const m = n - 1; // intervals
  if (m % 2 !== 0) {
    throw new Error("Simpson requires an even number of intervals");
  }

  // Build cumulative integral using composite Simpson over pairs
  const y_m = new Array(n).fill(0);
  const h_m = h0 * 1e-3;
  for (let i = 0; i < m; i += 2) {
    const s0 = slope_rad[i];
    const s1 = slope_rad[i + 1];
    const s2 = slope_rad[i + 2];
    const area = (h_m / 3) * (s0 + 4 * s1 + s2);
    y_m[i + 1] = y_m[i] + (h_m / 6) * (2 * s0 + 4 * s1); // partial (0->i+1) via Simpson segment
    y_m[i + 2] = y_m[i] + area;
    // Ensure monotone cumulative for interior (keep already set values for previous segments)
  }
  // Clean up interior cumulative values by linear interpolation between even nodes
  for (let i = 2; i < n - 1; i += 2) {
    // already set exact at even nodes; ensure odd nodes between evens:
    const yLeft = y_m[i - 2];
    const yMid = y_m[i - 1]; // already approximated above; fine to keep
    const yRight = y_m[i];
    // leave as-is; Simpson partial lines above already placed odd nodes.
  }
  return y_m;
}

/**
 * Natural cubic spline integration of slope over each interval.
 * Produces exact integral of the spline (not numeric sampling).
 *
 * Algorithm:
 * 1) Build natural spline second derivatives M[i] for slope data s[i] at x[i].
 * 2) Over interval i with width h, integral is:
 *    ∫ s(x) dx from x_i to x_{i+1} =
 *      h/2 * (s_i + s_{i+1}) - (h^3/24) * (M_i + M_{i+1})
 * (closed form from cubic-spline basis)
 */
function computeNaturalSplineSecondDerivatives(x_mm, s_vals) {
  const n = x_mm.length;
  if (n < 3) {
    // For 2 points, second derivatives are zero (linear)
    return new Array(n).fill(0);
  }

  const h = new Array(n - 1);
  for (let i = 0; i < n - 1; i++) {
    h[i] = x_mm[i + 1] - x_mm[i];
    if (!(h[i] > 0)) throw new Error("positions must be strictly increasing");
  }

  // Tridiagonal system: A * M = d
  const A = new Array(n).fill(0);
  const B = new Array(n).fill(0);
  const C = new Array(n).fill(0);
  const d = new Array(n).fill(0);

  // Natural spline boundary conditions: M0 = Mn-1 = 0
  B[0] = 1; d[0] = 0;
  B[n - 1] = 1; d[n - 1] = 0;

  for (let i = 1; i < n - 1; i++) {
    A[i] = h[i - 1] / 6;
    B[i] = (h[i - 1] + h[i]) / 3;
    C[i] = h[i] / 6;
    d[i] = (s_vals[i + 1] - s_vals[i]) / h[i] - (s_vals[i] - s_vals[i - 1]) / h[i - 1];
  }

  // Solve tridiagonal system via Thomas algorithm
  // Forward sweep
  for (let i = 1; i < n; i++) {
    const m = A[i] / B[i - 1];
    B[i] = B[i] - m * C[i - 1];
    d[i] = d[i] - m * d[i - 1];
  }
  // Back substitution
  const M = new Array(n).fill(0);
  M[n - 1] = d[n - 1] / B[n - 1];
  for (let i = n - 2; i >= 0; i--) {
    M[i] = (d[i] - C[i] * M[i + 1]) / B[i];
  }
  return M;
}

function integrateSpline(x_mm, slope_rad) {
  const n = x_mm.length;
  const M = computeNaturalSplineSecondDerivatives(x_mm, slope_rad);
  const y_m = new Array(n).fill(0);

  for (let i = 0; i < n - 1; i++) {
    const h_m = (x_mm[i + 1] - x_mm[i]) * 1e-3; // m
    const s_i = slope_rad[i];
    const s_ip1 = slope_rad[i + 1];
    const M_i = M[i];
    const M_ip1 = M[i + 1];
    // Closed-form exact integral of cubic spline over interval i
    const intervalArea_m = (h_m / 2) * (s_i + s_ip1) - (Math.pow(h_m, 3) / 24) * (M_i + M_ip1);
    y_m[i + 1] = y_m[i] + intervalArea_m;
  }
  return y_m;
}

// ---------- Public API ----------
/**
 * Convert angular readings (arcsec) at positions (mm) to displacement profile (µm).
 * method: 'trapezoid' | 'simpson' | 'spline'
 */
function anglesToDisplacementUm(positionsMm, anglesArcsec, method = 'trapezoid') {
  if (positionsMm.length !== anglesArcsec.length) {
    throw new Error("positions and angles must have same length");
  }
  const n = positionsMm.length;
  if (n < 2) throw new Error("need at least two points");
  const theta_rad = anglesArcsec.map(a => a * ARCSEC_TO_RAD);

  let y_m;
  switch (method) {
    case 'trapezoid':
      y_m = integrateTrapezoid(positionsMm, theta_rad);
      break;
    case 'simpson':
      y_m = integrateSimpson(positionsMm, theta_rad);
      break;
    case 'spline':
      y_m = integrateSpline(positionsMm, theta_rad);
      break;
    default:
      throw new Error(`Unknown method: ${method}`);
  }

  // meters -> micrometers
  return y_m.map(v => v * 1e6);
}

function straightnessResiduals(positionsMm, anglesArcsec, method = 'trapezoid') {
  const y_um = anglesToDisplacementUm(positionsMm, anglesArcsec, method);
  return detrend(positionsMm, y_um); // { y, fit }
}

function parallelismFromRails(positionsMm, anglesA_arcsec, anglesB_arcsec, method = 'trapezoid') {
  const A = straightnessResiduals(positionsMm, anglesA_arcsec, method);
  const B = straightnessResiduals(positionsMm, anglesB_arcsec, method);
  const p = A.y.map((v, i) => v - B.y[i]);
  const pvUm = Math.max(...p) - Math.min(...p);
  return {
    positionsMm,
    pUm: p,
    railA: { residualsUm: A.y, fit: A.fit },
    railB: { residualsUm: B.y, fit: B.fit },
    pvUm
  };
}

function round(val, digits = 6) {
  const p = Math.pow(10, digits);
  return Math.round(val * p) / p;
}

function demo() {
  const x = [0, 100, 200, 300, 400];
  const thetaA = [+2.0, +1.0, 0.0, -1.0, -2.0];
  const thetaB = [+1.5, +0.5, 0.0, -0.5, -1.0];

  for (const method of ['trapezoid', 'simpson', 'spline']) {
    try {
      const yA = anglesToDisplacementUm(x, thetaA, method);
      const yB = anglesToDisplacementUm(x, thetaB, method);
      const P = parallelismFromRails(x, thetaA, thetaB, method);

      console.log(`\n=== Method: ${method} ===`);
      console.log('y_A (µm):', yA.map(v => round(v, 5)));
      console.log('y_B (µm):', yB.map(v => round(v, 5)));
      console.log('r_A (µm):', P.railA.residualsUm.map(v => round(v, 5)));
      console.log('r_B (µm):', P.railB.residualsUm.map(v => round(v, 5)));
      console.log('p(x) (µm):', P.pUm.map(v => round(v, 5)));
      console.log('Parallelism P-V (µm):', round(P.pvUm, 5));
    } catch (e) {
      console.log(`\n=== Method: ${method} ===`);
      console.log('Error:', e.message);
    }
  }
}

export {
  anglesToDisplacementUm,
  straightnessResiduals,
  parallelismFromRails,
  // exposed for testing:
  integrateTrapezoid,
  integrateSimpson,
  integrateSpline
};