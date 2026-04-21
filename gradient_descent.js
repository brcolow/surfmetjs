import {
  distanceToSegmentSquared,
  isPointInPolygon
} from './utils.js';
import {
  refineRoundnessCenter,
  roundnessResultFromCenter
} from './roundness_refinement.js';

const DEFAULT_ADAPTIVE_OPTIONS = {
  gamma: 200,
  initialStepSize: 1,
  learningRateDecay: 0.75,
  slopeDamper: 0.25,
  maxIterations: 1000,
  maxIterationsLineSearch: 100,
  epsilon: 1e-8,
  verbose: false
};

const DISTANCE_EPSILON = 1e-12;
const OUTSIDE_MIC_PENALTY = 1e6;

function getInitialCenter(points, initialEstimate) {
  if (Array.isArray(initialEstimate)) {
    return [initialEstimate[0], initialEstimate[1]];
  }

  if (initialEstimate?.center) {
    return [initialEstimate.center[0], initialEstimate.center[1]];
  }

  if (Number.isFinite(initialEstimate?.a) && Number.isFinite(initialEstimate?.b)) {
    return [initialEstimate.a, initialEstimate.b];
  }

  return meanPoint(points);
}

function getInitialCenters(points, initialEstimate) {
  return uniqueCenters([
    getInitialCenter(points, initialEstimate),
    meanPoint(points),
    boundingBoxCenter(points)
  ]);
}

function uniqueCenters(centers) {
  const seen = new Set();
  const result = [];

  for (const center of centers) {
    if (!center || !Number.isFinite(center[0]) || !Number.isFinite(center[1])) {
      continue;
    }

    const key = `${center[0]},${center[1]}`;
    if (!seen.has(key)) {
      seen.add(key);
      result.push(center);
    }
  }

  return result;
}

function meanPoint(points) {
  const sum = points.reduce(
    ([sumX, sumY], [x, y]) => [sumX + x, sumY + y],
    [0, 0]
  );

  return [sum[0] / points.length, sum[1] / points.length];
}

function normalizeAdaptiveOptions(optionsOrGamma, learningRateDecay, slopeDamper) {
  if (typeof optionsOrGamma === 'number') {
    return {
      ...DEFAULT_ADAPTIVE_OPTIONS,
      gamma: optionsOrGamma,
      learningRateDecay: learningRateDecay ?? DEFAULT_ADAPTIVE_OPTIONS.learningRateDecay,
      slopeDamper: slopeDamper ?? DEFAULT_ADAPTIVE_OPTIONS.slopeDamper
    };
  }

  return {
    ...DEFAULT_ADAPTIVE_OPTIONS,
    ...optionsOrGamma
  };
}

/**
 * Performs gradient descent to optimize the circle center. Radius is derived
 * from the optimized center so MIC/MCC/MZC results are always geometrically
 * consistent with the point set.
 *
 * @param {Array<Array<number>>} points - The list of points to evaluate.
 * @param {string} circleType - The type of circle ("MIC", "MCC", or "MZC").
 * @param {Object|Array<number>} initialEstimate - Initial center estimate.
 * @param {number} [learningRate=1] - Initial step size for line search.
 * @param {number} [maxIterations=1000] - The maximum number of iterations.
 * @returns {Object} The best solution found.
 */
function gradientDescent(points, circleType, initialEstimate, learningRate = 1, maxIterations = 1000) {
  return adaptiveGradientDescent(points, circleType, initialEstimate, {
    initialStepSize: learningRate,
    maxIterations
  });
}

/**
 * Performs gradient descent with backtracking line search.
 *
 * @param {Array<Array<number>>} points - The list of points to evaluate.
 * @param {string} circleType - The type of circle ("MIC", "MCC", or "MZC").
 * @param {Object|Array<number>} initialEstimate - Initial center estimate.
 * @param {Object|number} [optionsOrGamma] - Options object, or legacy gamma value.
 * @param {number} [legacyLearningRateDecay] - Legacy decay value.
 * @param {number} [legacySlopeDamper] - Legacy slope damper value.
 * @returns {Object} The best solution found.
 */
function adaptiveGradientDescent(
  points,
  circleType,
  initialEstimate,
  optionsOrGamma = {},
  legacyLearningRateDecay,
  legacySlopeDamper
) {
  assertPoints(points);

  const {
    gamma,
    initialStepSize,
    learningRateDecay,
    slopeDamper,
    maxIterations,
    maxIterationsLineSearch,
    epsilon,
    verbose
  } = normalizeAdaptiveOptions(optionsOrGamma, legacyLearningRateDecay, legacySlopeDamper);

  let bestCenter = null;
  let bestObjective = Infinity;

  const initialCenters = getInitialCenters(points, initialEstimate)
    .filter(center => circleType !== "MIC" || isPointInPolygon(center, points));
  const centers = initialCenters.length > 0 ? initialCenters : [meanPoint(points)];

  for (const initialCenter of centers) {
    const center = optimizeCenter(
      initialCenter,
      points,
      circleType,
      gamma,
      initialStepSize,
      learningRateDecay,
      slopeDamper,
      maxIterations,
      maxIterationsLineSearch,
      epsilon
    );
    const objective = getObjectiveFunction(circleType, center, points, gamma);

    if (objective < bestObjective) {
      bestCenter = center;
      bestObjective = objective;
    }
  }

  const refinedCenter = refineRoundnessCenter(points, circleType, bestCenter);

  if (verbose) {
    printGradientDescentSummary(refinedCenter, points, circleType, gamma);
  }

  return roundnessResultFromCenter(points, circleType, refinedCenter);
}

function optimizeCenter(
  initialCenter,
  points,
  circleType,
  gamma,
  initialStepSize,
  learningRateDecay,
  slopeDamper,
  maxIterations,
  maxIterationsLineSearch,
  epsilon
) {
  let center = [...initialCenter];
  let objective = getObjectiveFunction(circleType, center, points, gamma);

  for (let iteration = 0; iteration < maxIterations; iteration++) {
    const [dx, dy] = gradient(center, points, circleType, 100, gamma);
    const gradNormSquared = dx * dx + dy * dy;

    if (gradNormSquared <= epsilon) {
      break;
    }

    let stepSize = initialStepSize;
    let accepted = false;
    let bestCandidate = center;
    let bestCandidateObjective = objective;

    for (let lineStep = 0; lineStep < maxIterationsLineSearch; lineStep++) {
      const candidate = [
        center[0] - stepSize * dx,
        center[1] - stepSize * dy
      ];
      const candidateObjective = getObjectiveFunction(circleType, candidate, points, gamma);

      if (candidateObjective < bestCandidateObjective) {
        bestCandidate = candidate;
        bestCandidateObjective = candidateObjective;
      }

      if (candidateObjective <= objective - slopeDamper * stepSize * gradNormSquared) {
        center = candidate;
        objective = candidateObjective;
        accepted = true;
        break;
      }

      stepSize *= learningRateDecay;
    }

    if (!accepted && bestCandidateObjective < objective) {
      center = bestCandidate;
      objective = bestCandidateObjective;
      accepted = true;
    }

    if (!accepted) {
      break;
    }
  }

  return center;
}

function assertPoints(points) {
  if (!Array.isArray(points) || points.length === 0) {
    throw new Error("At least one point is required.");
  }
}

function boundingBoxCenter(points) {
  const xs = points.map(([x]) => x);
  const ys = points.map(([, y]) => y);
  return [
    (Math.min(...xs) + Math.max(...xs)) / 2,
    (Math.min(...ys) + Math.max(...ys)) / 2
  ];
}

/**
 * Prints final optimization diagnostics after gradient descent.
 *
 * @param {Array<number>} center - Final optimized center [x, y].
 * @param {Array<Array<number>>} points - List of 2D points evaluated.
 * @param {string} circleType - Circle type ("MIC", "MCC", "MZC").
 * @param {number} gamma - Smoothing parameter.
 */
function printGradientDescentSummary(center, points, circleType, gamma = 200) {
  const finalObjective = getObjectiveFunction(circleType, center, points, gamma);
  const grad = gradient(center, points, circleType, 100, gamma);
  const gradNorm = Math.hypot(grad[0], grad[1]);

  console.log("=== Gradient Descent Summary ===");
  console.log(`Final objective value: ${finalObjective.toExponential(6)}`);
  console.log(`Final gradient norm: ${gradNorm.toExponential(6)}`);
  console.log(`Final center: (${center[0]}, ${center[1]})`);
  console.log("================================");
}

function boundaryDistancesFromCenter(center, points) {
  return points.map((point, index) => {
    const previous = points[(index + points.length - 1) % points.length];
    return Math.sqrt(distanceToSegmentSquared(center, previous, point));
  });
}

function boundaryClearance(center, points) {
  return Math.min(...boundaryDistancesFromCenter(center, points));
}

function smoothRadiusStats(center, points, gamma = 200) {
  const [x, y] = center;
  const radii = points.map(([px, py]) => ({
    dx: x - px,
    dy: y - py,
    radius: Math.max(Math.hypot(px - x, py - y), DISTANCE_EPSILON)
  }));

  const rawRadii = radii.map(({ radius }) => radius);
  const rMax = Math.max(...rawRadii);
  const rMin = Math.min(...rawRadii);
  const expMax = rawRadii.map(radius => Math.exp(gamma * (radius - rMax)));
  const expMin = rawRadii.map(radius => Math.exp(gamma * (rMin - radius)));
  const sumExpMax = expMax.reduce((sum, value) => sum + value, 0);
  const sumExpMin = expMin.reduce((sum, value) => sum + value, 0);
  const maxWeights = expMax.map(value => value / sumExpMax);
  const minWeights = expMin.map(value => value / sumExpMin);

  return {
    radii,
    smoothMax: rMax + Math.log(sumExpMax) / gamma,
    smoothMin: rMin - Math.log(sumExpMin) / gamma,
    smoothRange: rMax - rMin + (Math.log(sumExpMax) - Math.log(sumExpMin)) / gamma,
    maxWeights,
    minWeights
  };
}

function objectiveFunctionMIC(params, points, gamma = 200) {
  if (points.length < 3 || !isPointInPolygon(params, points)) {
    return OUTSIDE_MIC_PENALTY + boundaryClearance(params, points);
  }

  const distances = boundaryDistancesFromCenter(params, points);
  const dMin = Math.min(...distances);
  const sumExpMin = distances.reduce(
    (sum, distance) => Math.exp(gamma * (dMin - distance)) + sum,
    0
  );
  const smoothMin = dMin - Math.log(sumExpMin) / gamma;
  return -smoothMin;
}

function objectiveFunctionMCC(params, points, gamma = 200) {
  return smoothRadiusStats(params, points, gamma).smoothMax;
}

/**
 * Smoothed objective function for Minimum Zone Circle (MZC) using log-sum-exp.
 *
 * @param {Array<number>} params - The circle center parameters [x, y].
 * @param {Array<Array<number>>} points - The list of points to evaluate.
 * @param {number} [gamma=200] - Higher values approximate true max/min closer.
 * @returns {number} The approximated MZC zone width.
 */
function objectiveFunctionMZCSmooth(params, points, gamma = 200) {
  const [x, y] = params;
  const radii = points.map(point => Math.sqrt((point[0] - x) ** 2 + (point[1] - y) ** 2));
  const rMax = Math.max(...radii);
  const rMin = Math.min(...radii);
  const sumExpPos = radii.reduce((sum, radius) => Math.exp(gamma * (radius - rMax)) + sum, 0);
  const sumExpNeg = radii.reduce((sum, radius) => Math.exp(gamma * (rMin - radius)) + sum, 0);

  return (Math.log(sumExpPos) - Math.log(sumExpNeg)) / gamma + rMax - rMin;
}

/**
 * Computes the gradient of the smoothed objective for MIC, MCC, or MZC.
 *
 * @param {Array<number>} params - The current circle center [x, y].
 * @param {Array<Array<number>>} points - The list of points to evaluate.
 * @param {string} circleType - The type of circle ("MIC", "MCC", or "MZC").
 * @param {number} [_penalty=100] - Kept for backwards-compatible call sites.
 * @param {number} [gamma=200] - Smoothing parameter.
 * @returns {Array<number>} The gradient [dx, dy, 0].
 */
function gradient(params, points, circleType, _penalty = 100, gamma = 200) {
  if (circleType === "MIC" || circleType === "MZC") {
    return finiteDifferenceGradient(
      center => getObjectiveFunction(circleType, center, points, gamma),
      params
    );
  }

  const { radii, maxWeights, minWeights } = smoothRadiusStats(params, points, gamma);
  let dx = 0;
  let dy = 0;

  for (let i = 0; i < radii.length; i++) {
    const derivativeX = radii[i].dx / radii[i].radius;
    const derivativeY = radii[i].dy / radii[i].radius;
    let weight;

    switch (circleType) {
      case "MIC":
        weight = -minWeights[i];
        break;
      case "MCC":
        weight = maxWeights[i];
        break;
      default:
        throw new Error("Invalid circleType");
    }

    dx += weight * derivativeX;
    dy += weight * derivativeY;
  }

  return [dx, dy, 0];
}

function finiteDifferenceGradient(objective, center) {
  const stepX = 1e-6 * Math.max(1, Math.abs(center[0]));
  const stepY = 1e-6 * Math.max(1, Math.abs(center[1]));
  const dx = (
    objective([center[0] + stepX, center[1]]) -
    objective([center[0] - stepX, center[1]])
  ) / (2 * stepX);
  const dy = (
    objective([center[0], center[1] + stepY]) -
    objective([center[0], center[1] - stepY])
  ) / (2 * stepY);

  return [dx, dy, 0];
}

/**
 * Evaluates and returns the objective value based on the circle type.
 *
 * @param {string} circleType - The type of circle ("MIC", "MCC", or "MZC").
 * @param {Array<number>} params - Center parameters [x, y]; a radius value is ignored.
 * @param {Array<[number, number]>} points - The set of 2D points.
 * @param {number} [gamma=200] - Smoothing parameter.
 * @returns {number} The computed objective value.
 */
const getObjectiveFunction = (circleType, params, points, gamma = 200) => {
  switch (circleType) {
    case "MIC":
      return objectiveFunctionMIC(params, points, gamma);
    case "MCC":
      return objectiveFunctionMCC(params, points, gamma);
    case "MZC":
      return objectiveFunctionMZCSmooth(params, points, gamma);
    default:
      throw new Error("Invalid circleType");
  }
};

/**
 * Nicely prints the result of a circle optimization (MIC, MCC, or MZC).
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

export {
  adaptiveGradientDescent,
  gradient,
  gradientDescent,
  getObjectiveFunction,
  objectiveFunctionMZCSmooth,
  printCircleResult
};
