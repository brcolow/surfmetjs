import { computeConvexHull, distanceSquared, getInitialSolution } from './utils.js';

const gamma = 50; // approximation parameter to simulate max-min functions with log-sum-exp (softmax) 

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

  let params = [ initialSolution.center[0], initialSolution.center[1], initialSolution.radius ];
  let previousObjective = Infinity;
  let previousParams = params;
  for (let i = 0; i < maxIterations; i++) {
    const grad = gradient(params, points, circleType);
    params = params.map((p, i) => p - learningRate * grad[i]);
    const currentObjective = getObjectiveFunction(circleType)(params, points);
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
 * Performs a gradient descent with adaptive learning rate (Backtracking line search)
 * Remark: The code here is only for MZC
 * 
 * @param {Array<Array<number>>} points - The list of points to evaluate.
 * @param {string} circleType - The type of circle ("MIC", "MCC", or "MZC").
 * @param {Object} initialEstimate - Initial estimate of the circle parameters.
 * @param {number} [learningRateDecay=0.75] - The decay factor for the learning rate.
 * @param {number} [slopeDamper=0.25] - The damping factor for the slope at the current point.
 * @param {number} [maxIterations=1000] - The maximum number of iterations.
 * @param {number} [epsilon=1e-8] - Stop if the gradient norm is below epsilon.
 * @returns {Object} The best solution found as an Object: {center: [x, y], radius: r}
 */
function adaptiveGradientDescent(points, circleType, initialEstimate, learningRateDecay = 0.75, slopeDamper=0.25, maxIterations = 1000, epsilon = 1e-8) {
  const convexHull = computeConvexHull(points);
  const initialSolution = getInitialSolution(points, circleType, initialEstimate, convexHull);

  let params = [
    (initialSolution.outerCircle.center[0] + initialSolution.innerCircle.center[0]) / 2,
    (initialSolution.outerCircle.center[1] + initialSolution.innerCircle.center[1]) / 2,
    null // We will not use the variable r
  ];

  for (let i = 0; i < maxIterations; i++) {
    const grad = gradient(params, points, circleType);
    if (grad[0] ** 2 + grad[1] ** 2 < epsilon) {
      console.log("Stopping gradient descent due to small gradient.");
      break;
    }
    const dampedSlope = slopeDamper * (grad[0] ** 2 + grad[1] ** 2);
    let t = 1.0; // Line search parameter
    const objective = getObjectiveFunction(circleType)(params, points);
    while (true) {
      let paramsT = params.map((p, i) => p - grad[i] );
      let objectiveT = getObjectiveFunction(circleType)(paramsT, points);
      if (objectiveT < objective - t * dampedSlope) {
        params = paramsT;
        break; // Accept backtracked params + t * grad
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
 * Attempts to minimize the zone width (difference between the maximum and minimum radius)
 * while penalizing invalid inner circles.
 * 
 * @param {Array<number>} params - The circle parameters [x, y, r].
 * @param {Array<Array<number>>} points - The list of points to evaluate.
 * @returns {number} The objective value for the MZC.
 */
function objectiveFunctionMZC(params, points) {
  const [x, y, r] = params;
  // r will not play any role and will be ignored in the following computations

  const radii = points.map(point => Math.sqrt(distanceSquared(point, [x, y])));
  let rMax = Math.max(radii);
  let rMin = Math.min(radii);

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
 * @returns {Array<number>} The gradient [dx, dy, dr] for x, y, and r respectively.
 */
function gradient(params, points, circleType, penalty = 100) {
  const [x, y, r] = params;
  // Again, r and dr will not be used
  let dx = 0, dy = 0, dr = 0;

  if (circleType === "MZC") {

    const radii = points.map(point => Math.sqrt(distanceSquared(point, [x, y])));
    let rMax = Math.max(radii);
    let rMin = Math.min(radii);

    let expPos = radii.map(r => Math.exp(gamma * (r - rMax)));
    let expNeg = radii.map(r => Math.exp(gamma * (rMin - r)));

    let sumExpPos = expPos.reduce((sum, ep) => ep + sum, 0);
    let sumExpNeg = expNeg.reduce((sum, en) => en + sum, 0);

    let weights = radii.map((r, i) => (expPos[i] / sumExpPos - expNeg[i] / sumExpNeg) / r);

    dx = - weights.reduce((sum, w, i) => sum + w * (points[i][0] - x))
    dy = - weights.reduce((sum, w, i) => sum + w * (points[i][1] - y))

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
 * Selects the appropriate objective function based on the circle type.
 * 
 * @param {string} circleType - The type of circle ("MIC", "MCC", or "MZC").
 * @returns {Function} The objective function for the specified circle type.
 */
const getObjectiveFunction = (circleType) => {
  switch (circleType) {
      case "MIC":
          return objectiveFunctionMIC;
      case "MCC":
          return objectiveFunctionMCC;
      case "MZC":
          return objectiveFunctionMZC;
      default:
          // Handle invalid circleType
          throw new Error("Invalid circleType");
  }
};

export { adaptiveGradientDescent, gradientDescent }
