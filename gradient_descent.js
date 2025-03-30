import { computeConvexHull, distanceSquared, getInitialSolution } from './utils.js';

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
  let dx = 0, dy = 0, dr = 0;

  if (circleType === "MZC") {
    let rMax = -Infinity;
    let rMin = Infinity;
    let rMaxPoint = null;
    let rMinPoint = null;

    for (const point of points) {
      const dist = Math.sqrt(distanceSquared(point, [x, y]));
      if (dist > rMax) {
        rMax = dist;
        rMaxPoint = point;
      }
      if (dist < rMin) {
        rMin = dist;
        rMinPoint = point;
      }
    }

    if (rMaxPoint) {
      const [px, py] = rMaxPoint;
      const dist = Math.sqrt(distanceSquared(rMaxPoint, [x, y]));
      const factor = 1 / (dist + 1e-9); // Avoid division by zero
      dx += (x - px) * factor;
      dy += (y - py) * factor;
      dr += 1;
    }

    if (rMinPoint) {
      const [px, py] = rMinPoint;
      const dist = Math.sqrt(distanceSquared(rMinPoint, [x, y]));
      const factor = 1 / (dist + 1e-9); // Avoid division by zero
      dx -= (x - px) * factor;
      dy -= (y - py) * factor;
      dr -= 1;
    }

    if (rMin > r) {
      dr += penalty * (rMin - r);
    }
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

export { gradientDescent }
