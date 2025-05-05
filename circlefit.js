import { inv, norm, matrix, multiply, pinv, subtract, transpose } from 'mathjs';

/*
 * FitCircle Java class for fitting circles to 2D coordinate data
 *
 * Copyright 2009 2010 Michael Doube
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * Kåsa circle fit.
 * 
 * Adapted from Java https://github.com/mdoube/BoneJ/blob/1c5a90873f5d2dda0294dabc9698f002670aa816/src/org/doube/geometry/FitCircle.java#L53
 * 
 * @param {*} points 2d array containing x,y points (double[n][2])
 * @returns circle object centered at (a, b) with radius r
 */
function kasaFit(points) {
  const nPoints = points.length;
  if (nPoints < 3) {
    throw new Error("Too few points: must have at least 3 but was: " + nPoints);
  }
  const z = twoDArray(nPoints, 1);
  const xy1 = twoDArray(nPoints, 3);
  for (let n = 0; n < nPoints; n++) {
    z[n][0] = points[n][0] * points[n][0] + points[n][1] * points[n][1];
    xy1[n][0] = points[n][0];
    xy1[n][1] = points[n][1];
    xy1[n][2] = 1;
  }

  const XY1 = matrix(xy1);
  const Z = matrix(z);
  const P = XY1.size()[0] === XY1.size()[1] ? multiply(inv(XY1), Z) : multiply(pinv(XY1), Z);

  const p0 = P.get([0, 0]);
  const p1 = P.get([1, 0]);
  const p2 = P.get([2, 0]);
  return { "a": p0 / 2, "b": p1 / 2, "r": Math.sqrt((p0 * p0 + p1 * p1) / 4 + p2) };
}

/**
 * Fits a circle to a set of 2D points using the Levenberg–Marquardt algorithm.
 *
 * Adapted from Nikolai Chernov's original MATLAB script for geometric circle fitting
 * using a damped Gauss–Newton method (Levenberg–Marquardt), via Michael Doube's Java port in BoneJ:
 * https://github.com/mdoube/BoneJ/blob/1c5a90873f5d2dda0294dabc9698f002670aa816/src/org/doube/geometry/FitCircle.java
 *
 * Original publication:
 * Al-Sharadqah & Chernov (2009). "Error analysis for circle fitting algorithms."
 * Electronic Journal of Statistics, 3, 886–911. https://doi.org/10.1214/09-EJS419
 *
 * This function performs geometric total least squares (TLS) fitting — minimizing the
 * orthogonal (Euclidean) distance from each point to the estimated circle.
 * The algorithm uses a Levenberg–Marquardt iteration to optimize center (a, b) and radius r.
 *
 * @param {Array<Array<number>>} points - An array of [x, y] coordinate pairs.
 * @param {number} [lambdaIni=1] - Initial damping parameter (higher = more gradient descent-like).
 * @returns {{a: number, b: number, r: number}} Best-fit circle parameters (center and radius).
 * 
 * @see https://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm
 * @see https://en.wikipedia.org/wiki/Total_least_squares
 * @see https://en.wikipedia.org/wiki/Circle_fitting
 * @see https://github.com/mdoube/BoneJ/blob/1c5a90873f5d2dda0294dabc9698f002670aa816/src/org/doube/geometry/FitCircle.java#L500
 */
function circleFitOrthogonal(points, lambdaIni = 1) {
  const nPoints = points.length;
  if (nPoints < 3) {
    throw new Error("Too few points: must have at least 3 but was: " + nPoints);
  }

  // Initial guess using Kåsa method (algebraic OLS fit to x² + y² + Dx + Ey + F = 0).
  // This is a fast linear approximation, but biased and not geometric.
  // Used only to seed the nonlinear TLS refinement.
  const guess = kasaFit(points);
  let x = guess.a
  let y = guess.b;
  let r = guess.r;
  const par = [[x, y, r]];
  let Par = matrix(par);      // Current parameters
  let ParTemp = matrix(par);  // Trial parameters
  const epsilon = 1e-6;
  let progress = epsilon;
  const iterMax = 50;
  let lambda_sqrt = Math.sqrt(lambdaIni);

  let f = 0; // Initial cost
  let j = twoDArray(nPoints + 3, 3); // Jacobian
  let g = twoDArray(nPoints + 3, 1); // Residuals

  // Compute initial residuals and Jacobian.
  // This matches the logic from Chernov’s `CurrentIteration` function in MATLAB.
  for (let i = 0; i < nPoints; i++) {
    const dX = points[i][0] - x;
    const dY = points[i][1] - y;
    const d = Math.sqrt(dX * dX + dY * dY);
    j[i][0] = -dX / d;
    j[i][1] = -dY / d;
    j[i][2] = -1;
    g[i][0] = d - r;
    f += (d - r) * (d - r);
  }

  const J = matrix(j);
  const G = matrix(g);

  let fTemp = 0;
  const jTemp = twoDArray(nPoints + 3, 3);
  const gTemp = twoDArray(nPoints + 3, 1);

  for (let iter = 0; iter < iterMax; iter++) {
    let safety = 0;

    // Inner loop: adjust lambda to ensure improvement.
    // Matches the "secondary loop" in Chernov’s original MATLAB (while (1)).
    while (safety < 100) {
      safety++;

      // Reset fTemp before each inner safety iteration.
      // In the original MATLAB, fTemp is recomputed from scratch each time.
      // During porting, this reset was initially omitted, causing accumulation bugs.
      fTemp = 0;

      // Augment Jacobian and residuals with damping terms for LM step.
      J.set([nPoints, 0], lambda_sqrt);
      J.set([nPoints + 1, 1], lambda_sqrt);
      J.set([nPoints + 2, 2], lambda_sqrt);
      G.set([nPoints, 0], 0);
      G.set([nPoints + 1, 0], 0);
      G.set([nPoints + 2, 0], 0);

      // Solve for update step ΔPar using pseudo-inverse (robust to non-square J).
      // This is equivalent to solving (JᵀJ + λI) Δθ = Jᵀg.
      const DelPar = J.size()[0] === J.size()[1] ? multiply(inv(J), G) : multiply(pinv(J), G);
      progress = norm(DelPar, 'fro') / (norm(Par, 'fro') + epsilon);
      if (progress < epsilon) {
        break;
      }

      // Trial update
      ParTemp = subtract(Par, transpose(DelPar));
      x = ParTemp.get([0, 0]);
      y = ParTemp.get([0, 1]);
      r = ParTemp.get([0, 2]);

      // Recompute residuals and Jacobian for trial step.
      for (let i = 0; i < nPoints; i++) {
        const dX = points[i][0] - x;
        const dY = points[i][1] - y;
        const d = Math.sqrt(dX * dX + dY * dY);
        jTemp[i][0] = -dX / d;
        jTemp[i][1] = -dY / d;
        jTemp[i][2] = -1;
        gTemp[i][0] = d - r;
        fTemp += (d - r) * (d - r);
      }

      // Accept update if it improves cost and radius remains valid.
      if (fTemp < f && r > 0) {
        lambda_sqrt /= 2;
        break;
      }

      // Otherwise, increase lambda and try again.
      lambda_sqrt *= 2;
      continue;
    }
    if (progress < epsilon) {
      break;
    }

    // Accept the new parameters and update state
    Par = ParTemp;
    j = jTemp;
    g = gTemp;
    f = fTemp;
  }
  return {
    a: Par.get([0, 0]),
    b: Par.get([0, 1]),
    r: Par.get([0, 2]),
  };
}

function twoDArray(rows, cols) {
  let x = new Array(rows);

  for (let i = 0; i < rows; i++) {
    x[i] = new Array(cols).fill(0);
  }
  return x;
}

export { circleFitOrthogonal, kasaFit }
