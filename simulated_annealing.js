import { areConcentric, computeConvexHull, distanceSquared, findFarthestPoint, findNearestPoint, getInitialSolution, getRandomBetween, isPointInPolygon } from './utils.js';

/**
 * Simulated annealing optimizer with equilibrium steps and (optional) final gradient descent refinement.
 *
 * @param {Array<Array<number>>} points - The input point set.
 * @param {string} circleType - "MIC", "MCC", or "MZC".
 * @param {Object} initialEstimate - Initial estimate of circle parameters.
 * @param {number} [initialTemperature=1000] - Starting temperature.
 * @param {Object} [coolingType={ type: 'logarithmic', rate: 0.1 }] - Cooling schedule object.
 * @param {number} [maxIterations=1000] - Number of temperature drops.
 * @param {number} [equilibriumSteps=20] - Number of neighbors tried per temperature.
 * @param {number} [maxNeighborIterations=5000] - Max tries to generate valid neighbor.
 * @param {number} stepSize - Perturbation scale.
 * @returns {Object} The best solution found as an Object: {center: [x, y], radius: r}
 */
function simulatedAnnealing(points, circleType, initialEstimate, initialTemperature = 1000, coolingType = { type: 'logarithmic', rate: 0.1 }, maxIterations = 1000, equilibriumSteps = 20, maxNeighborIterations = 5000, stepSize) {
  const initialSolution = getInitialSolution(points, circleType, initialEstimate);

  let currentSolution = initialSolution;
  let bestSolution = initialSolution;
  const convexHull = computeConvexHull(points);
  let currentEnergy = calculateEnergy(currentSolution, circleType, points, convexHull);
  let bestEnergy = currentEnergy;

  const minTemperature = 1e-5;
  let didRefine = false;

  for (let i = 0; i < maxIterations; i++) {
    const temperature = temperatureSchedule(initialTemperature, i, coolingType);
    const adaptiveStepSize = stepSize * temperature;

    for (let j = 0; j < equilibriumSteps; j++) {
      const neighbor = generateNeighbor(points, convexHull, circleType, currentSolution, maxNeighborIterations, adaptiveStepSize);
      if (!neighbor) {
        continue;
      }

      const newEnergy = calculateEnergy(neighbor, circleType, points, convexHull);
      const deltaEnergy = newEnergy - currentEnergy;

      const accept = deltaEnergy < 0 || Math.exp(-deltaEnergy / Math.max(temperature, minTemperature)) > Math.random();

      if (accept) {
        currentSolution = neighbor;
        currentEnergy = newEnergy;

        if (newEnergy < bestEnergy) {
          bestSolution = neighbor;
          bestEnergy = newEnergy;
        }
      }
    }

    if (i + 1 === maxIterations) {
      console.log("Reached max iterations for SA: " + i);
    }
  }

  // console.log("Best solution found at iteration = " + bestIteration + ", energy = " + bestEnergy + ", temperature = " + bestTemperature);
  return bestSolution;
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
        throw new Error("The cooling rate for linear cooling must strictly greater than 0.");
      }
      return initialTemperature * coolingType.rate ** i;
    case 'logarithmic':
      return initialTemperature / Math.log(i + 2);
  }
}

/**
 * Generates a valid neighbor solution for MIC, MCC, or MZC by perturbing the current solution.
 *
 * @param {Array<Array<number>>} points - 2D point cloud.
 * @param {Array<Array<number>>} convexHull - Convex hull for center constraint.
 * @param {string} circleType - "MIC", "MCC", or "MZC".
 * @param {Object} currentSolution - The current solution object.
 * @param {number} maxNeighborIterations - Max attempts to find a valid neighbor.
 * @param {number} stepSize - Perturbation scale.
 * @returns {Object|null} - Valid neighbor or null if none found.
 */
function generateNeighbor(points, convexHull, circleType, currentSolution, maxNeighborIterations, stepSize) {
  for (let attempt = 0; attempt < maxNeighborIterations; attempt++) {
    if (circleType === "MZC") {
      // Clone concentric center and radii
      let newCenter = [...currentSolution.outerCircle.center];
      let innerRadius = currentSolution.innerCircle.radius;
      let outerRadius = currentSolution.outerCircle.radius;

      // Perturb center
      if (Math.random() < 0.6) {
        newCenter = [
          newCenter[0] + getRandomBetween(-stepSize, stepSize),
          newCenter[1] + getRandomBetween(-stepSize, stepSize)
        ];
        if (!isPointInPolygon(newCenter, convexHull)) continue;
      }

      // Perturb inner and/or outer radii independently
      if (Math.random() < 0.6) {
        innerRadius += getRandomBetween(-stepSize, stepSize);
        innerRadius = Math.max(0, innerRadius);
      }
      if (Math.random() < 0.6) {
        outerRadius += getRandomBetween(-stepSize, stepSize);
        outerRadius = Math.max(innerRadius + 1e-6, outerRadius); // Ensure outer > inner
      }

      const neighbor = {
        innerCircle: { center: newCenter, radius: innerRadius },
        outerCircle: { center: newCenter, radius: outerRadius }
      };

      const valid = points.every(p => {
        const d = Math.sqrt(distanceSquared(p, newCenter));
        return d >= innerRadius && d <= outerRadius;
      });

      if (valid) return neighbor;
    } else {
      // MIC or MCC
      const changeCenter = Math.random() < 0.6;
      const changeRadius = Math.random() < 0.6;

      let newCenter = [...currentSolution.center];
      let newRadius = currentSolution.radius;

      if (changeCenter) {
        const dx = getRandomBetween(-stepSize, stepSize);
        const dy = getRandomBetween(-stepSize, stepSize);
        newCenter = [newCenter[0] + dx, newCenter[1] + dy];
        if (!isPointInPolygon(newCenter, convexHull)) continue;
      }

      if (changeRadius) {
        const dr = getRandomBetween(0, stepSize);
        if (circleType === "MIC") {
          newRadius += dr;
        } else if (circleType === "MCC") {
          newRadius = Math.max(0, newRadius - dr);
        }
      }

      const neighbor = {
        center: newCenter,
        radius: newRadius
      };

      const valid = points.every(p => {
        const d2 = distanceSquared(p, newCenter);
        if (circleType === "MIC") {
          return d2 >= newRadius * newRadius;
        } else {
          // MCC
          return d2 <= newRadius * newRadius;
        }
      });

      if (valid) {
        return neighbor;
      }
    }
  }

  return null;
}


/**
 * Calculates the energy (cost) of a given circle solution based on geometric criteria.
 * Lower energy means a better solution.
 *
 * For MIC: maximize the minimum distance to any point (i.e., maximize inscribed radius).
 * For MCC: minimize the maximum distance to any point (i.e., minimize circumscribed radius).
 * For MZC: minimize the radial difference between the outer and inner concentric circles,
 *          ensuring all points lie within the band.
 *
 * @param {Object} solution - The circle solution to evaluate.
 * @param {string} circleType - The type of circle ("MIC", "MCC", or "MZC").
 * @param {Array<Array<number>>} points - The set of 2D points.
 * @param {Array<Array<number>>} convexHull - The convex hull of the point set (for MZC validation).
 * @returns {number} - The energy score; lower values are better.
 */
function calculateEnergy(solution, circleType, points, convexHull) {
  switch (circleType) {
    case "MIC": {
      const minDist = Math.min(...points.map(p => Math.sqrt(distanceSquared(p, solution.center))));
      return -minDist; // Maximize min distance
    }

    case "MCC": {
      const maxDist = Math.max(...points.map(p => Math.sqrt(distanceSquared(p, solution.center))));
      return maxDist; // Minimize max distance
    }

    case "MZC": {
      const inner = solution.innerCircle;
      const outer = solution.outerCircle;

      // Penalty if points are outside the band or not concentric
      const outerViolations = points.some(p => Math.sqrt(distanceSquared(p, outer.center)) > outer.radius);
      const innerViolations = points.some(p => Math.sqrt(distanceSquared(p, inner.center)) < inner.radius);
      const notConcentric = !areConcentric(inner, outer);

      if (outerViolations || innerViolations || notConcentric) {
        return Infinity;
      }

      const bandWidth = outer.radius - inner.radius;
      return bandWidth; // Minimize width of the band
    }
  }
}

export { simulatedAnnealing };
