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
  const convexHull = computeConvexHull(points);
  const initialSolution = getInitialSolution(points, circleType, initialEstimate, convexHull);

  let params = [initialSolution.center[0], initialSolution.center[1], initialSolution.radius];
  let previousObjective = Infinity;
  let previousParams = params;
  for (let i = 0; i < maxIterations; i++) {
    const grad = gradient(params, points, circleType);
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

  return { center: [params[0], params[1]], radius: params[2] };
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
  const initialSolution = getInitialSolution(points, circleType, initialEstimate, null);

  let params;
  if (circleType === "MZC") {
    params = [
      (initialSolution.outerCircle.center[0] + initialSolution.innerCircle.center[0]) / 2,
      (initialSolution.outerCircle.center[1] + initialSolution.innerCircle.center[1]) / 2,
      null // We will not use the variable r
    ];
  } else {
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
      // let objectiveT = objectiveFunctionMZC(paramsT, points, gamma);
      let objectiveT = getObjectiveFunction(circleType, paramsT, points, gamma);
      if (objectiveT < objective - t * dampedSlope) {
        params = paramsT;
        break; // Accept backtracked params - t * grad
      } else {
        t *= learningRateDecay; // scale down the steplength
      }
    }
  }
  // Here again in MZC, radius will not be used
  return { center: [params[0], params[1]], radius: params[2] };
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
    const dist = distanceSquared(point, [x, y]);
    if (dist > r * r) {
      penalty += 100;
    }
  }
  return r + penalty;
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
function gradient(params, points, circleType, penalty = 100, gamma = 200) {
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

// Fix over-zealous exports just for testing.
export { adaptiveGradientDescent, gradientDescent, objectiveFunctionMZCSmooth, gradient }
