import { getInitialSolution, computeConvexHull } from './simulated_annealing.js';

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

function gradient(params, points, penalty = 100) {
  const [x, y, r] = params;
  let dx = 0, dy = 0, dr = -1;
  for (const point of points) {
    const [px, py] = point;
    const dist = distanceSquared(point, [x, y]);
    if (dist < r * r) {
      const factor = penalty * (r - dist) / dist;
      dx += (x - px) * factor;
      dy += (y - py) * factor;
      dr += factor;
    }
  }
  return [dx, dy, dr];
}

function findMIC(points, leastSquaresCircle, learningRate = 0.00000000001, maxIterations = 1000) {
  const convexHull = computeConvexHull(points);
  const initialSolution = getInitialSolution(points, "MIC", leastSquaresCircle, convexHull);

  let params = [ initialSolution.center[0], initialSolution.center[1], initialSolution.radius ];
  let previousObjective = Infinity;
  let previousParams = params;
  for (let i = 0; i < maxIterations; i++) {
    const grad = gradient(params, points);
    params = params.map((p, i) => p - learningRate * grad[i]);
    const currentObjective = objectiveFunctionMIC(params, points);
   // console.log(`Iteration ${i}: Objective value: ${currentObjective}`);
    if (currentObjective > previousObjective) {
      console.log("Objective function increased, stopping at iteration: " + i);
      params = previousParams;
      break;
    }
    previousParams = params;
    previousObjective = currentObjective;
 }

  return params;
}

export { findMIC }
