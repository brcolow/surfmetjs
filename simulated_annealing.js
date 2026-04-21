import {
  computeConvexHull,
  isPointInPolygon
} from './utils.js';
import {
  refineRoundnessCenter,
  roundnessResultFromCenter
} from './roundness_refinement.js';

/**
 * Simulated annealing optimizer with equilibrium steps.
 *
 * The annealed state is the candidate center. MIC uses clearance to the
 * reconstructed closed profile boundary; MCC/MZC use extracted profile points.
 *
 * @param {Array<Array<number>>} points - The input point set.
 * @param {string} circleType - "MIC", "MCC", or "MZC".
 * @param {Object|Array<number>} initialEstimate - Initial estimate.
 * @param {number} stepSize - Initial perturbation scale.
 * @param {Object} [options] - Optional configuration parameters.
 * @param {number} [options.initialTemperature=1000] - Starting temperature.
 * @param {Object} [options.coolingType={ type: 'logarithmic', rate: 0.1 }] - Cooling schedule object.
 * @param {number} [options.maxIterations=1000] - Number of temperature drops.
 * @param {number} [options.equilibriumSteps=20] - Number of neighbors tried per temperature.
 * @param {number} [options.maxNeighborIterations=5000] - Max tries to generate a valid neighbor.
 * @param {Function} [options.random=Math.random] - Random source returning values in [0, 1).
 * @returns {Object} The best solution found.
 */
function simulatedAnnealing(
  points,
  circleType,
  initialEstimate,
  stepSize,
  {
    initialTemperature = 1000,
    coolingType = { type: 'logarithmic', rate: 0.1 },
    maxIterations = 1000,
    equilibriumSteps = 20,
    maxNeighborIterations = 5000,
    random = Math.random
  } = {}
) {
  assertAnnealingInputs(points, circleType);
  assertRandom(random);

  const convexHull = computeConvexHull(points.map(point => [...point]));
  const baseStepSize = getBaseStepSize(points, stepSize);
  let currentSolution = getInitialSolution(points, circleType, initialEstimate, convexHull);
  let currentEnergy = calculateEnergy(currentSolution, circleType);
  let bestSolution = currentSolution;
  let bestEnergy = currentEnergy;

  for (let i = 0; i < maxIterations; i++) {
    const temperature = temperatureSchedule(initialTemperature, i, coolingType);
    const temperatureRatio = initialTemperature > 0 ? temperature / initialTemperature : 1;
    const adaptiveStepSize = baseStepSize * Math.max(temperatureRatio, 1e-6);

    for (let j = 0; j < equilibriumSteps; j++) {
      const neighbor = generateNeighbor(
        points,
        convexHull,
        circleType,
        currentSolution,
        maxNeighborIterations,
        adaptiveStepSize,
        random
      );

      if (!neighbor) {
        continue;
      }

      const newEnergy = calculateEnergy(neighbor, circleType);
      const deltaEnergy = newEnergy - currentEnergy;

      if (shouldAcceptNeighbor(deltaEnergy, temperature, random)) {
        currentSolution = neighbor;
        currentEnergy = newEnergy;

        if (newEnergy < bestEnergy) {
          bestSolution = neighbor;
          bestEnergy = newEnergy;
        }
      }
    }
  }

  const refinedCenter = refineRoundnessCenter(
    points,
    circleType,
    centerFromSolution(bestSolution)
  );
  return roundnessResultFromCenter(points, circleType, refinedCenter);
}

function assertAnnealingInputs(points, circleType) {
  if (!Array.isArray(points) || points.length < 3) {
    throw new Error("At least three points are required.");
  }

  if (!["MIC", "MCC", "MZC"].includes(circleType)) {
    throw new Error("Invalid circleType");
  }
}

function assertRandom(random) {
  if (typeof random !== 'function') {
    throw new Error("random must be a function.");
  }
}

function getBaseStepSize(points, stepSize) {
  if (Number.isFinite(stepSize) && stepSize > 0) {
    return stepSize;
  }

  return boundingBoxDiagonal(points) / 100;
}

function getInitialSolution(points, circleType, initialEstimate, convexHull) {
  const centerDomain = circleType === "MIC" ? points : convexHull;
  const candidates = getInitialCenters(points, initialEstimate)
    .filter(center => isAllowedCenter(center, circleType, points, centerDomain));

  if (candidates.length === 0) {
    candidates.push(meanPoint(points));
  }

  return candidates
    .map(center => solutionFromCenter(points, circleType, center))
    .reduce((best, candidate) => {
      return calculateEnergy(candidate, circleType) < calculateEnergy(best, circleType)
        ? candidate
        : best;
    });
}

function getInitialCenters(points, initialEstimate) {
  const centers = [];

  if (Array.isArray(initialEstimate)) {
    centers.push([initialEstimate[0], initialEstimate[1]]);
  }

  if (initialEstimate?.center) {
    centers.push([initialEstimate.center[0], initialEstimate.center[1]]);
  }

  if (Number.isFinite(initialEstimate?.a) && Number.isFinite(initialEstimate?.b)) {
    centers.push([initialEstimate.a, initialEstimate.b]);
  }

  centers.push(meanPoint(points));
  centers.push(boundingBoxCenter(points));

  return uniqueCenters(centers);
}

function uniqueCenters(centers) {
  const seen = new Set();
  const result = [];

  for (const center of centers) {
    if (!Number.isFinite(center[0]) || !Number.isFinite(center[1])) {
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

function boundingBoxCenter(points) {
  const xs = points.map(([x]) => x);
  const ys = points.map(([, y]) => y);
  return [
    (Math.min(...xs) + Math.max(...xs)) / 2,
    (Math.min(...ys) + Math.max(...ys)) / 2
  ];
}

function boundingBoxDiagonal(points) {
  const xs = points.map(([x]) => x);
  const ys = points.map(([, y]) => y);
  return Math.hypot(
    Math.max(...xs) - Math.min(...xs),
    Math.max(...ys) - Math.min(...ys)
  );
}

function centerFromSolution(solution) {
  return solution.center || solution.outerCircle.center;
}

function solutionFromCenter(points, circleType, center) {
  return roundnessResultFromCenter(points, circleType, center);
}

function shouldAcceptNeighbor(deltaEnergy, temperature, random) {
  if (deltaEnergy < 0) {
    return true;
  }

  if (temperature <= 0 || !Number.isFinite(deltaEnergy)) {
    return false;
  }

  return random() < Math.exp(-deltaEnergy / temperature);
}

/**
 * Computes the temperature based on the current iteration and cooling schedule.
 *
 * @param {number} initialTemperature - The starting temperature.
 * @param {number} i - The current iteration index.
 * @param {Object} coolingType - The cooling schedule, with type and rate properties.
 * @returns {number} - The temperature for the current iteration.
 */
function temperatureSchedule(initialTemperature, i, coolingType) {
  switch (coolingType.type) {
    case 'exponential':
      if (coolingType.rate <= 0 || coolingType.rate >= 1) {
        throw new Error("The cooling rate for exponential cooling must be strictly between 0 and 1.");
      }
      return initialTemperature * coolingType.rate ** i;
    case 'linear':
      if (coolingType.rate <= 0) {
        throw new Error("The cooling rate for linear cooling must be strictly greater than 0.");
      }
      return Math.max(0, initialTemperature - coolingType.rate * i);
    case 'logarithmic':
      return initialTemperature / Math.log(i + 2);
    default:
      throw new Error("Unknown cooling schedule type: " + coolingType.type);
  }
}

/**
 * Generates a valid neighbor solution for MIC, MCC, or MZC by perturbing the
 * current center.
 *
 * @param {Array<Array<number>>} points - 2D point cloud.
 * @param {Array<Array<number>>} convexHull - Convex hull for center constraint.
 * @param {string} circleType - "MIC", "MCC", or "MZC".
 * @param {Object} currentSolution - The current solution object.
 * @param {number} maxNeighborIterations - Max attempts to find a valid neighbor.
 * @param {number} stepSize - Perturbation scale.
 * @returns {Object|null} - Valid neighbor or null if none found.
 */
function generateNeighbor(
  points,
  convexHull,
  circleType,
  currentSolution,
  maxNeighborIterations,
  stepSize,
  random
) {
  const currentCenter = centerFromSolution(currentSolution);

  for (let attempt = 0; attempt < maxNeighborIterations; attempt++) {
    const newCenter = [
      currentCenter[0] + randomBetween(random, -stepSize, stepSize),
      currentCenter[1] + randomBetween(random, -stepSize, stepSize)
    ];

    if (!isAllowedCenter(newCenter, circleType, points, convexHull)) {
      continue;
    }

    return solutionFromCenter(points, circleType, newCenter);
  }

  return null;
}

function randomBetween(random, min, max) {
  return random() * (max - min) + min;
}

function isAllowedCenter(center, circleType, points, convexHull) {
  const domain = circleType === "MIC" ? points : convexHull;
  return isPointInPolygon(center, domain);
}

/**
 * Calculates the energy (cost) of a given circle solution. Lower is better.
 *
 * @param {Object} solution - The circle solution to evaluate.
 * @param {string} circleType - The type of circle ("MIC", "MCC", or "MZC").
 * @returns {number} - The energy score.
 */
function calculateEnergy(solution, circleType) {
  switch (circleType) {
    case "MIC":
      return -solution.radius;
    case "MCC":
      return solution.radius;
    case "MZC":
      return solution.zoneWidth;
    default:
      throw new Error("Invalid circleType");
  }
}

export { simulatedAnnealing };
