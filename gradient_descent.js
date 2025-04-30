import {
  computeConvexHull,
  distanceSquared,
  getInitialSolution
} from './utils.js';

/**
 * Performs gradient descent to optimize the circle parameters.
 * 
 * @param {Array<Array<number>>} points - The list of points to evaluate.
 * @param {string} circleType - The type of circle ("MIC", "MCC", or "MZC").
 * @param {Object} initialEstimate - Initial estimate of the circle parameters.
 * @param {number} [learningRate=1e-11] - The step size for gradient updates.
 * @param {number} [maxIterations=1000] - The maximum number of iterations.
 * @returns {Object} The best solution found as an Object: {center: [x, y], radius: r}
 */
function gradientDescent(points, circleType, initialEstimate, learningRate = 0.00000000001, maxIterations = 1000) {
  let params;
  if (circleType === "MZC") {
    const initialSolution = getInitialSolution(points, circleType, initialEstimate);
    params = [
      (initialSolution.outerCircle.center[0] + initialSolution.innerCircle.center[0]) / 2,
      (initialSolution.outerCircle.center[1] + initialSolution.innerCircle.center[1]) / 2,
      initialSolution.radius
    ];
  } else {
    const initialSolution = getInitialSolution(points, circleType, initialEstimate);
    params = [initialSolution.center[0], initialSolution.center[1], initialSolution.radius];
  }

  let previousObjective = Infinity;
  let previousParams = params;
  for (let i = 0; i < maxIterations; i++) {
    const grad = gradient2(params, points, circleType);
    params = params.map((p, i) => p - learningRate * grad[i]);
    const currentObjective = getObjectiveFunction(circleType, params, points);
    // console.log(`Iteration ${i}: Objective value: ${currentObjective}`);
    if (currentObjective > previousObjective) {
      console.log("Objective function increased, stopping at iteration: " + i);
      params = previousParams;
      break;
    }
    previousParams = params;
    previousObjective = currentObjective;
  }

  if (circleType === "MZC") {
    return computeMZCResult([params[0], params[1]], points);
  } else if (circleType === "MIC") {
    return computeMICResult([params[0], params[1]], params[2], points);
  } else if (circleType === "MCC") {
    return computeMCCResult([params[0], params[1]], params[2], points);
  }
}

/**
 * Performs a gradient descent with adaptive learning rate (backtracking line search).
 * 
 * @param {Array<Array<number>>} points - The list of points to evaluate.
 * @param {string} circleType - The type of circle ("MIC", "MCC", or "MZC").
 * @param {Object} initialEstimate - Initial estimate of the circle parameters.
 * @param {number} [gamma=200] - The approximation parameter to simulate max-min functions with log-sum-exp (default: 200)
 * @param {number} [learningRateDecay=0.75] - The decay factor for the learning rate.
 * @param {number} [slopeDamper=0.25] - The damping factor for the slope at the current point.
 * @param {number} [maxIterations=1000] - The maximum number of iterations.
 * @param {number} [maxIterationsLineSearch=100] - The maximum number of (inner) iterations during line searchs.
 * @param {number} [epsilon=1e-8] - Stop if the gradient norm is below epsilon.
 * @returns {Object} The best solution found as an Object: {center: [x, y], radius: r}
 */
function adaptiveGradientDescent(points, circleType, initialEstimate, gamma = 200, learningRateDecay = 0.75, slopeDamper = 0.25, maxIterations = 1000, maxIterationsLineSearch = 100, epsilon = 1e-8) {
  let params;
  if (circleType === "MZC") {
    const initialSolution = getInitialSolution(points, circleType, initialEstimate);
    params = [
      (initialSolution.outerCircle.center[0] + initialSolution.innerCircle.center[0]) / 2,
      (initialSolution.outerCircle.center[1] + initialSolution.outerCircle.center[1]) / 2
    ];
  } else {
    const initialSolution = getInitialSolution(points, circleType, initialEstimate);
    params = [initialSolution.center[0], initialSolution.center[1], initialSolution.radius];
  }

  for (let i = 0; i < maxIterations; i++) {
    const grad = gradient(params, points, circleType, 100, gamma);
    if (grad[0] ** 2 + grad[1] ** 2 < epsilon) {
      console.log("Stopping gradient descent due to small gradient.");
      break;
    }
    const dampedSlope = slopeDamper * (grad[0] ** 2 + grad[1] ** 2);
    let t = 1.0; // Line search parameter
    const objective = getObjectiveFunction(circleType, params, points, gamma);
    for (let j = 0; j < maxIterationsLineSearch; j++) {
      let paramsT = params.map((p, i) => p - t * grad[i]);
      let objectiveT = getObjectiveFunction(circleType, paramsT, points, gamma);
      if (objectiveT < objective - t * dampedSlope) {
        params = paramsT;
        break; // Accept backtracked params - t * grad
      } else {
        t *= learningRateDecay; // scale down the steplength
      }
    }
  }

  printGradientDescentSummary(params, points, circleType, gamma);

  if (circleType === "MZC") {
    return computeMZCResult([params[0], params[1]], points);
  } else if (circleType === "MIC") {
    return computeMICResult([params[0], params[1]], params[2], points);
  } else if (circleType === "MCC") {
    return computeMCCResult([params[0], params[1]], params[2], points);
  }
}

function armijoWolfeGradientDescent(points, circleType, initialEstimate, gamma = 200, maxIterations = 1000, epsilon = 1e-16) {
  let params;
  if (circleType === "MZC") {
    const initialSolution = getInitialSolution(points, circleType, initialEstimate);
    params = [
      (initialSolution.outerCircle.center[0] + initialSolution.innerCircle.center[0]) / 2,
      (initialSolution.outerCircle.center[1] + initialSolution.outerCircle.center[1]) / 2
    ];
  } else {
    const initialSolution = getInitialSolution(points, circleType, initialEstimate);
    params = [initialSolution.center[0], initialSolution.center[1], initialSolution.radius];
  }

  for (let i = 0; i < maxIterations; i++) {
    const grad = gradient(params, points, circleType, 100, gamma);

    const gradNormSquared = grad.reduce((sum, g) => sum + g * g, 0);
    if (gradNormSquared < epsilon) {
      console.log(`Converged: small gradient norm at iteration ${i}`);
      break;
    }

    params = lineSearchArmijoWolfe(params, grad, points, circleType);
  }

  printGradientDescentSummary(params, points, circleType, gamma);

  if (circleType === "MZC") {
    return computeMZCResult([params[0], params[1]], points);
  } else if (circleType === "MIC") {
    return computeMICResult([params[0], params[1]], params[2], points);
  } else if (circleType === "MCC") {
    return computeMCCResult([params[0], params[1]], params[2], points);
  }
}

/**
 * Prints final optimization diagnostics after gradient descent.
 *
 * @param {Array<number>} params - Final optimized parameters [x, y, r].
 * @param {Array<Array<number>>} points - List of 2D points evaluated.
 * @param {string} circleType - Circle type ("MIC", "MCC", "MZC").
 * @param {number} gamma - Smoothing parameter (only needed for MZC).
 */
function printGradientDescentSummary(params, points, circleType, gamma = 200) {
  const finalObjective = getObjectiveFunction(circleType, params, points, gamma);
  const grad = gradient(params, points, circleType, 100, gamma);
  const gradNorm = Math.sqrt(grad.reduce((sum, g) => sum + g * g, 0));

  console.log("=== Gradient Descent Summary ===");
  console.log(`Final objective value: ${finalObjective.toExponential(6)}`);
  console.log(`Final gradient norm: ${gradNorm.toExponential(6)}`);
  console.log(`Final center: (${params[0]}, ${params[1]})`);
  if (params.length === 3) {
    console.log(`Final radius: ${params[2]}`);
  }
  console.log("================================");
}

/**
 * Performs Armijo-Wolfe line search to find an optimal step size.
 *
 * @param {Array<number>} params - Current circle parameters [x, y, r].
 * @param {Array<number>} grad - Current gradient [dx, dy, dr].
 * @param {Array<Array<number>>} points - The list of points.
 * @param {string} circleType - The circle type ("MIC", "MCC", or "MZC").
 * @param {number} [initialStepSize=1.0] - Initial step size guess.
 * @param {number} [c1=1e-4] - Armijo sufficient decrease parameter.
 * @param {number} [c2=0.9] - Wolfe curvature parameter.
 * @param {number} [maxLineSearchIterations=100] - Maximum number of attempts to calculate new params.
 * @returns {Array<number>} New parameters after a good step.
 */
function lineSearchArmijoWolfe(params, grad, points, circleType, initialStepSize = 1.0, c1 = 1e-4, c2 = 0.9, maxLineSearchIterations = 500) {
  let stepSize = initialStepSize;
  const currentObjective = getObjectiveFunction(circleType, params, points);
  const gradDot = grad.reduce((sum, g) => sum + g * g, 0); // ||grad||^2

  let attempts = 0;
  while (true) {
    const newParams = params.map((p, i) => p - stepSize * grad[i]);
    const newObjective = getObjectiveFunction(circleType, newParams, points);
    const newGrad = gradient(newParams, points, circleType);

    const gradInnerProduct = grad.reduce((sum, g, i) => sum + g * newGrad[i], 0);

    // Armijo condition (sufficient decrease)
    if (newObjective > currentObjective - c1 * stepSize * gradDot) {
      stepSize *= 0.5; // Too small decrease → reduce step size
    }
    // Wolfe condition (sufficient curvature)
    else if (gradInnerProduct < c2 * gradDot) {
      stepSize *= 2.0; // Step was too short → increase step size
    }
    else {
      return newParams; // Both conditions satisfied: accept step
    }

    attempts++;

    if (stepSize < 1e-15 || attempts > maxLineSearchIterations) {
      console.warn("Line search failed: step size too small or too many iterations, returning original params.");
      return params; // Abort line search
    }
  }
}

/**
 * Objective function for Minimum Inscribed Circle (MIC).
 * Attempts to minimize the radius while penalizing circles that don't enclose all points.
 * 
 * @param {Array<number>} params - The circle parameters [x, y, r].
 * @param {Array<Array<number>>} points - The list of points to evaluate.
 * @returns {number} The objective value for the MIC.
 */
function objectiveFunctionMIC(params, points) {
  const [x, y, r] = params;
  let penalty = 0;
  for (const point of points) {
    const dist = distanceSquared(point, [x, y]);
    if (dist < r * r) {
      penalty += 100;
    }
  }
  return -r + penalty;
}

/**
 * Objective function for Maximum Circumscribed Circle (MCC).
 * Attempts to maximize the radius while penalizing circles that exclude any points.
 * 
 * @param {Array<number>} params - The circle parameters [x, y, r].
 * @param {Array<Array<number>>} points - The list of points to evaluate.
 * @returns {number} The objective value for the MCC.
 */
function objectiveFunctionMCC(params, points) {
  const [x, y, r] = params;
  let penalty = 0;

  for (const point of points) {
    const dist = Math.sqrt(distanceSquared(point, [x, y]));
    const violation = dist - r;

    if (violation > 0) {
      // Smooth quadratic penalty: gentle push if outside
      penalty += 0.5 * violation * violation;
    }
  }

  // Main goal: minimize r while covering all points
  return -r + penalty;
}
/**
 * Objective function for Maximum Zone Circle (MZC).
 * Computes the zone width (difference between maximum and minimum distance to center)
 * and adds a penalty if the proposed circle radius is too small to enclose all points.
 * 
 * @param {Array<number>} params - The circle parameters [x, y, r] (center coordinates and radius).
 * @param {Array<Array<number>>} points - The list of points to evaluate.
 * @returns {number} The penalized objective value for the MZC.
 */
function objectiveFunctionMZC(params, points) {
  const [x, y, r] = params;

  let rMin = Infinity;
  let rMax = -Infinity;

  for (const point of points) {
    const dist = Math.sqrt(distanceSquared(point, [x, y]));
    rMin = Math.min(rMin, dist);
    rMax = Math.max(rMax, dist);
  }

  const zoneWidth = rMax - rMin;
  let penalty = 0;
  if (rMin > r) {
    penalty = (rMin - r) * 100;
  }
  return zoneWidth + penalty;
}

/**
 * Smoothed objective function for Maximum Zone Circle (MZC) using log-sum-exp approximation.
 * Approximates the maximum and minimum radii differences to provide a smoother optimization landscape.
 * 
 * @param {Array<number>} params - The circle center parameters [x, y] (no radius parameter).
 * @param {Array<Array<number>>} points - The list of points to evaluate.
 * @param {number} [gamma=200] - The smoothing parameter; higher values approximate the true max/min closer.
 * @returns {number} The approximated objective value for the MZC.
 */
function objectiveFunctionMZCSmooth(params, points, gamma = 200) {
  const [x, y] = params;

  const radii = points.map(point => Math.sqrt(distanceSquared(point, [x, y])));
  let rMax = Math.max(...radii);
  let rMin = Math.min(...radii);

  let sumExpPos = radii.reduce((sum, r) => Math.exp(gamma * (r - rMax)) + sum, 0);
  let sumExpNeg = radii.reduce((sum, r) => Math.exp(gamma * (rMin - r)) + sum, 0);

  return (Math.log(sumExpPos) - Math.log(sumExpNeg)) / gamma + rMax - rMin;
}

/**
 * Computes the gradient of the objective function for a given circle type (MIC, MCC, or MZC),
 * using smooth approximations to ensure continuous behavior near constraint boundaries.
 *
 * For MIC and MCC, a tanh smoothing is used instead of sharp penalties to avoid zero gradients.
 * For MZC, log-sum-exp smoothing is used to approximate the maximum and minimum radii.
 * 
 * @param {Array<number>} params - The current circle parameters [x, y, r].
 * @param {Array<Array<number>>} points - The list of points to evaluate.
 * @param {string} circleType - The type of circle ("MIC", "MCC", or "MZC").
 * @param {number} [penalty=100] - The penalty factor for constraint violations (default: 100).
 * @param {number} [gamma=200] - The smoothing parameter controlling sharpness (default: 200).
 * @returns {Array<number>} The computed gradient [dx, dy, dr] for center (x, y) and radius r.
 */
function gradient(params, points, circleType, penalty = 100, gamma = 200) {
  const [x, y, r] = params;
  let dx = 0, dy = 0, dr = 0;

  if (circleType === "MZC") {
    // Smooth MZC using log-sum-exp approximation
    const radii = [];
    let rMin = Infinity;
    let rMax = -Infinity;

    for (const [px, py] of points) {
      const ri = Math.sqrt(distanceSquared([px, py], [x, y])) + 1e-9;
      radii.push(ri);
      if (ri < rMin) rMin = ri;
      if (ri > rMax) rMax = ri;
    }

    const expPos = radii.map(ri => Math.exp(gamma * (ri - rMax)));
    const expNeg = radii.map(ri => Math.exp(gamma * (rMin - ri)));

    const sumExpPos = expPos.reduce((sum, val) => sum + val, 0);
    const sumExpNeg = expNeg.reduce((sum, val) => sum + val, 0);

    const weights = radii.map((ri, i) => (expPos[i] / sumExpPos - expNeg[i] / sumExpNeg) / ri);

    dx = -weights.reduce((sum, w, i) => sum + w * (points[i][0] - x), 0);
    dy = -weights.reduce((sum, w, i) => sum + w * (points[i][1] - y), 0);
    // dr remains zero for MZC
  } else {
    // MIC and MCC handled together
    const sign = (circleType === "MIC") ? +1 : -1;

    for (const [px, py] of points) {
      const dist = Math.sqrt(distanceSquared([px, py], [x, y])) + 1e-9;
      const violation = (circleType === "MIC") ? (r - dist) : (dist - r);
      const smoothPenalty = penalty * Math.tanh(gamma * violation);

      dx += (x - px) * smoothPenalty / dist;
      dy += (y - py) * smoothPenalty / dist;
      dr += sign * smoothPenalty;
    }
  }

  return [dx, dy, dr];
}

/**
 * Computes the gradient of the objective function for a given circle type.
 * The gradient is used to iteratively adjust the circle parameters to optimize the objective.
 * 
 * @param {Array<number>} params - The current circle parameters [x, y, r].
 * @param {Array<Array<number>>} points - The list of points to evaluate.
 * @param {string} circleType - The type of circle ("MIC", "MCC", or "MZC").
 * @param {number} penalty - The penalty factor for invalid configurations (default: 100).
 * @param {number} [gamma=200] - The approximation parameter to simulate max-min functions with log-sum-exp (default: 200). Remark: This parameter is only needed for the mzc gradient. Maybe refactor code?
 * @returns {Array<number>} The gradient [dx, dy, dr] for x, y, and r respectively.
 */
function gradient2(params, points, circleType, penalty = 100, gamma = 200) {
  const [x, y, r] = params;

  let dx = 0, dy = 0, dr = 0;

  if (circleType === "MZC") {
    const radii = [];

    let rMin = Infinity;
    let rMax = -Infinity;

    for (const [px, py] of points) {
      const ri = Math.sqrt(distanceSquared([px, py], [x, y]));
      radii.push(ri);
      if (ri < rMin) rMin = ri;
      if (ri > rMax) rMax = ri;
    }

    const expPos = radii.map(r => Math.exp(gamma * (r - rMax)));
    const expNeg = radii.map(r => Math.exp(gamma * (rMin - r)));
    const sumExpPos = expPos.reduce((sum, ep) => ep + sum, 0);
    const sumExpNeg = expNeg.reduce((sum, en) => en + sum, 0);
    const weights = radii.map((r, i) => (expPos[i] / sumExpPos - expNeg[i] / sumExpNeg) / r);
    dx = -weights.reduce((sum, w, i) => sum + w * (points[i][0] - x), 0);
    dy = -weights.reduce((sum, w, i) => sum + w * (points[i][1] - y), 0);

  } else {
    for (const point of points) {
      const [px, py] = point;
      const dist = Math.sqrt(distanceSquared(point, [x, y]));

      if (circleType === "MIC" && dist < r) {
        const factor = penalty * (r - dist) / (dist + 1e-9);
        dx += (x - px) * factor;
        dy += (y - py) * factor;
        dr -= factor;
      } else if (circleType === "MCC" && dist > r) {
        const factor = penalty * (dist - r) / (dist + 1e-9);
        dx += (x - px) * factor;
        dy += (y - py) * factor;
        dr += factor;
      }
    }
  }
  return [dx, dy, dr];
}

/**
 * Evaluates and returns the objective value based on the circle type.
 *
 * @param {string} circleType - The type of circle ("MIC", "MCC", or "MZC").
 * @param {Array<number>} params - The initial parameters (e.g., [x, y]) for the objective function.
 * @param {Array<[number, number]>} points - The set of 2D points to evaluate against.
 * @param {number} [gamma=200] - (Only for MZC) Approximation parameter for simulating max-min behavior using log-sum-exp.
 * @returns {number} The computed objective value.
 *
 */
const getObjectiveFunction = (circleType, params, points, gamma = 200) => {
  switch (circleType) {
    case "MIC":
      return objectiveFunctionMIC(params, points);
    case "MCC":
      return objectiveFunctionMCC(params, points);
    case "MZC":
      return objectiveFunctionMZCSmooth(params, points, gamma);
    default:
      throw new Error("Invalid circleType");
  }
};

/**
 * Computes the final result for the Minimum Zone Circle (MZC) after optimization.
 * 
 * The MZC minimizes the radial separation (zone width) between the largest and smallest 
 * distances from the center to the profile points. 
 * 
 * After optimization, this function computes:
 * - Outer radius (distance to farthest point)
 * - Inner radius (distance to nearest point)
 * - Zone width (outer radius - inner radius)
 * - Roundness (half of the zone width, per ISO 1101 definition)
 * 
 * @param {Array<number>} center - The optimized center of the MZC [x, y].
 * @param {Array<Array<number>>} points - The set of 2D points defining the profile.
 * @returns {Object} An object containing:
 *   - {Array<number>} center - The center [x, y] of the MZC.
 *   - {number} outerRadius - The maximum distance from center to a point (outer circle).
 *   - {number} innerRadius - The minimum distance from center to a point (inner circle).
 *   - {number} zoneWidth - The radial separation between outer and inner circles.
 *   - {number} roundness - Half of the zone width (standard roundness definition).
 */
function computeMZCResult(center, points) {
  const radii = points.map(p => Math.hypot(p[0] - center[0], p[1] - center[1]));
  const outerRadius = Math.max(...radii);
  const innerRadius = Math.min(...radii);
  const zoneWidth = outerRadius - innerRadius;
  const roundness = zoneWidth / 2;

  return {
    center: center,
    outerRadius: outerRadius,
    innerRadius: innerRadius,
    zoneWidth: zoneWidth,
    roundness: roundness,
    circleType: "MZC"
  };
}

/**
 * Computes the final result for the Maximum Inscribed Circle (MIC) after optimization.
 * 
 * The MIC is the largest circle fully contained within the point cloud.
 * After optimization, this function computes the maximum outward deviation 
 * from the MIC (RONt_MIC), which quantifies the roundness error.
 * 
 * @param {Array<number>} center - The optimized center of the MIC [x, y].
 * @param {number} radius - The optimized radius of the MIC.
 * @param {Array<Array<number>>} points - The set of 2D points defining the profile.
 * @returns {Object} An object containing:
 *   - {Array<number>} center - The center [x, y] of the MIC.
 *   - {number} radius - The radius of the MIC.
 *   - {number} roundness - The maximum outward deviation from the MIC (RONt_MIC).
 */
function computeMICResult(center, radius, points) {
  const deviations = points.map(([px, py]) => Math.hypot(px - center[0], py - center[1]) - radius);
  const maxDeviation = Math.max(...deviations);
  return {
    center: center,
    radius: radius,
    roundness: maxDeviation,
    circleType: "MIC"
  };
}

/**
 * Computes the final result for the Minimum Circumscribed Circle (MCC) after optimization.
 * 
 * The MCC is the smallest circle that fully encloses the point cloud.
 * After optimization, this function computes the maximum inward deviation 
 * from the MCC (RONt_MCC), which quantifies the roundness error.
 * 
 * @param {Array<number>} center - The optimized center of the MCC [x, y].
 * @param {number} radius - The optimized radius of the MCC.
 * @param {Array<Array<number>>} points - The set of 2D points defining the profile.
 * @returns {Object} An object containing:
 *   - {Array<number>} center - The center [x, y] of the MCC.
 *   - {number} radius - The radius of the MCC.
 *   - {number} roundness - The maximum inward deviation from the MCC (RONt_MCC).
 */
function computeMCCResult(center, radius, points) {
  const deviations = points.map(([px, py]) => radius - Math.hypot(px - center[0], py - center[1]));
  const maxDeviation = Math.max(...deviations);
  return {
    center: center,
    radius: radius,
    roundness: maxDeviation,
    circleType: "MCC"
  };
}

/**
 * Nicely prints the result of a circle optimization (MIC, MCC, or MZC).
 * Automatically detects the circle type based on the 'circleType' field.
 *
 * @param {Object} result - The result object.
 */
function printCircleResult(result) {
  const typeMap = {
    MIC: "Maximum Inscribed Circle (MIC)",
    MCC: "Minimum Circumscribed Circle (MCC)",
    MZC: "Minimum Zone Circle (MZC)"
  };

  const label = typeMap[result.circleType] || "Unknown Circle Type";
  console.log(`--- ${label} ---`);
  console.log(`Center: (${result.center[0].toFixed(10)}, ${result.center[1].toFixed(10)})`);

  if (result.circleType === "MZC") {
    console.log(`Outer Radius: ${result.outerRadius.toFixed(10)}`);
    console.log(`Inner Radius: ${result.innerRadius.toFixed(10)}`);
    console.log(`Zone Width: ${result.zoneWidth.toFixed(10)}`);
    console.log(`Roundness: ${result.roundness.toFixed(10)}`);
  } else if (result.circleType === "MIC" || result.circleType === "MCC") {
    console.log(`Radius: ${result.radius.toFixed(10)}`);
    console.log(`Roundness: ${result.roundness.toFixed(10)}`);
  } else {
    console.log("Unknown result format: " + result.radius);
  }
}

// Fix over-zealous exports just for testing.
export { adaptiveGradientDescent, gradient, gradientDescent, getObjectiveFunction, printCircleResult, objectiveFunctionMZCSmooth }
