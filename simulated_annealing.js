import { areConcentric, computeConvexHull, distanceSquared, findFarthestPoint, findNearestPoint, getInitialSolution, getRandomBetween, isPointInPolygon } from './utils.js';
import { gradientDescent } from './gradient_descent.js';

/**
 * Simulated annealing optimizer with equilibrium steps and (optional) final gradient descent refinement.
 *
 * @param {Array<Array<number>>} points - The input point set.
 * @param {string} circleType - "MIC", "MCC", or "MZC".
 * @param {Object} leastSquaresCircle - Initial estimate.
 * @param {number} initialTemperature - Starting temperature (default: 1000).
 * @param {Object} coolingType - Cooling schedule object (default: { type: 'logarithmic', rate: 0.1 }).
 * @param {number} maxIterations - Number of temperature drops (default: 1000).
 * @param {number} equilibriumSteps - Number of neighbors tried per temperature (default: 20).
 * @param {number} maxNeighborIterations - Max tries to generate valid neighbor (default: 5000).
 * @param {number} stepSize - Perturbation scale.
 * @param {boolean} refineWithGradientDescent - Use gradient descent to refine the solution once a low temperature is reached (default: false).
 * @returns {Object} - Best solution found.
 */
function simulatedAnnealing(points, circleType, leastSquaresCircle, initialTemperature = 1000, coolingType = { type: 'logarithmic', rate: 0.1 }, maxIterations = 1000, equilibriumSteps = 20, maxNeighborIterations = 5000, stepSize, refineWithGradientDescent = false) {
  const convexHull = computeConvexHull(points);
  const initialSolution = getInitialSolution(points, circleType, leastSquaresCircle, convexHull);

  let currentSolution = initialSolution;
  let bestSolution = initialSolution;
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

    // Gradient Descent Refinement at low temp
    /*
    if (refineWithGradientDescent) {
      if (i > 0.9 * maxIterations && !didRefine) {
        console.log("Switching to gradient descent refinement...");

        const refinedParams = gradientDescent(points, circleType, {
          center: bestSolution.center || bestSolution.outerCircle?.center,
          radius: bestSolution.radius || (bestSolution.outerCircle?.radius ?? 0)
        });

        if (circleType === "MIC" || circleType === "MCC") {
          bestSolution = {
            center: [refinedParams[0], refinedParams[1]],
            radius: refinedParams[2]
          };
        } else if (circleType === "MZC") {
          const center = [refinedParams[0], refinedParams[1]];
          const dists = points.map(p => Math.sqrt(distanceSquared(p, center)));
          const rMin = Math.min(...dists);
          const rMax = Math.max(...dists);
          
          bestSolution = {
            innerCircle: {
              center,
              radius: rMin
            },
            outerCircle: {
              center,
              radius: rMax
            }
          };
        }

        didRefine = true;
      }
    }
    */
    if (i + 1 === maxIterations) {
      console.log("Reached max iterations: " + i);
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
 * Generates a neighbor solution for the current solution using random perturbations.
 *
 * @param {Array<Array<number>>} points - The set of 2D points.
 * @param {Array<Array<number>>} convexHull - The convex hull of the points.
 * @param {string} circleType - The type of circle optimization ("MIC", "MCC", or "MZC").
 * @param {Object} currentSolution - The current solution state.
 * @param {number} maxNeighborIterations - The maximum attempts to generate a valid neighbor.
 * @param {number} stepSize - The maximum step size for generating neighbor solutions.
 * @returns {Object|null} - A valid neighbor solution or `null` if none could be found.
 */
function generateNeighborOld(points, convexHull, circleType, currentSolution, maxNeighborIterations, stepSize) {
  let foundValidNeighbor = false;
  for (let i = 0; i < maxNeighborIterations; i++) {
    let neighborCandidate;

    if (circleType !== "MZC") {
      const radiusOrCenterOrBoth = getRandomBetween(0, 100);
      const xOrY = getRandomBetween(0, 100);
      if (radiusOrCenterOrBoth < 33) {
        // 33% of the time we try generating a new neighbor by changing the center location.
        const centerXRandomizerExp = getRandomBetween(-8, 0);
        const centerYRandomizerExp = getRandomBetween(-8, 0);
        neighborCandidate = {
          center: [currentSolution.center[0] + getRandomBetween(-stepSize, stepSize) * (10 ** centerXRandomizerExp) * (xOrY < 50 ? 1 : 0), currentSolution.center[1] + getRandomBetween(-stepSize * stepSize) * (10 ** centerYRandomizerExp) * (xOrY > 50 ? 1 : 0)],
          radius: currentSolution.radius
        };
      } else if (radiusOrCenterOrBoth < 66) {
        // 33% of the time we change the radius (but never going smaller for the MIC and never going bigger for the MCC).
        const radiusRandomizerExp = getRandomBetween(-8, 0);
        let newRadius;
        if (circleType === "MIC") {
          newRadius = currentSolution.radius + getRandomBetween(0, stepSize) * (10 ** radiusRandomizerExp);
        } else {
          newRadius = currentSolution.radius - getRandomBetween(0, stepSize) * (10 ** radiusRandomizerExp);
        }
        neighborCandidate = {
          center: [currentSolution.center[0], currentSolution.center[1]],
          radius: newRadius
        };
      } else {
        // 33% of the time change both.
        const centerXRandomizerExp = getRandomBetween(-8, 0);
        const centerYRandomizerExp = getRandomBetween(-8, 0);
        const radiusRandomizerExp = getRandomBetween(-8, 0);
        let newRadius;
        if (circleType === "MIC") {
          newRadius = currentSolution.radius + getRandomBetween(0, stepSize) * (10 ** radiusRandomizerExp);
        } else {
          newRadius = currentSolution.radius - getRandomBetween(0, stepSize) * (10 ** radiusRandomizerExp);
        }
        neighborCandidate = {
          center: [currentSolution.center[0] + getRandomBetween(-stepSize, stepSize) * (10 ** centerXRandomizerExp) * (xOrY < 50 ? 1 : 0), currentSolution.center[1] + getRandomBetween(-stepSize * stepSize) * (10 ** centerYRandomizerExp) * (xOrY > 50  ? 1 : 0)],
          radius: newRadius
        };
      }
    } else {
      // When finding the neighbor for MZC we want to share the same center because concentricity must be maintained.
      if (getRandomBetween(0, 100) > 50) {
        // Half of the time we try generating a new neighbor by changing the center location.
        let center = [currentSolution.center[0] + getRandomBetween(-stepSize, stepSize), currentSolution.center[1] + getRandomBetween(-stepSize, stepSize)];
        neighborCandidate = {
          outerCircle: {
            center: center,
            radius: currentSolution.radius
          }, innerCircle: {
            center: center,
            radius: currentSolution.radius
          }
        };
      } else {
        // TODO: Should we only change one radius at a time?
        // The other half of the time we try changing the radius (TODO: but never going smaller for the inner circle and never going larger for the outer circle?).
        neighborCandidate = {
          outerCircle: {
            center: [currentSolution.center[0], currentSolution.center[1]],
            radius: currentSolution.radius + getRandomBetween(-stepSize, stepSize)
          }, innerCircle: {
            center: [currentSolution.center[0], currentSolution.center[1]],
            radius: currentSolution.radius + getRandomBetween(-stepSize, stepSize)
          }
        };
      }
    }

    foundValidNeighbor = true;

    switch (circleType) {
      case "MIC":
        // Make sure no point is strictly inside (can be on) the circle.
        // FIXME: This is not working correctly because sometimes the center moves outside the circle trace completely and thus no points are inside but that's not a valid candidate...
        for (const point of points) {
          const distance = distanceSquared(point, neighborCandidate.center);
          if (distance < neighborCandidate.radius * neighborCandidate.radius) {
            foundValidNeighbor = false;
            break;
          }
        }
        // Need to also make sure that we didn't move center outside the circle trace. We do this by making sure the center is still inside the convex hull of the point cloud.
        if (foundValidNeighbor) {
          if (!isPointInPolygon(neighborCandidate.center, convexHull)) {
            foundValidNeighbor = false;
            break;
          }
        }
        break;
      case "MCC":
        // Make sure no point is strictly outside (can be on) the circle.
        for (const point of points) {
          const distance = distanceSquared(point, neighborCandidate.center);
          if (distance > neighborCandidate.radius * neighborCandidate.radius) {
            foundValidNeighbor = false;
            break;
          }
        }
        // Need to also make sure that we didn't move center outside the circle trace. We do this by making sure the center is still inside the convex hull of the point cloud.
        if (foundValidNeighbor) {
          if (!isPointInPolygon(neighborCandidate.center, convexHull)) {
            foundValidNeighbor = false;
            break;
          }
        }
        break;
      case "MZC":
        // Make sure no point is strictly inside (can be on) the inner circle.
        for (const point of points) {
          const distance = distanceSquared(point, neighborCandidate.outerCircle.center);
          if (distance < neighborCandidate.radius * neighborCandidate.radius) {
            foundValidNeighbor = false;
            break;
          }
        }
        if (!foundValidNeighbor) {
          return null;
        }
        // Make sure no point is strictly outside (can be on) the outer circle.
        for (const point of points) {
          const distance = distanceSquared(point, neighborCandidate.innerCircle.center);
          if (distance > neighborCandidate.radius * neighborCandidate.radius) {
            foundValidNeighbor = false;
            break;
          }
        }
        // Need to also make sure that we didn't move center outside the circle trace. We do this by making sure the center is still inside the convex hull of the point cloud.
        if (foundValidNeighbor) {
          if (!isPointInPolygon(neighborCandidate.center, convexHull)) {
            foundValidNeighbor = false;
            break;
          }
        }
        // Make sure both circles are concentric.
        if (foundValidNeighbor) {
          foundValidNeighbor = areConcentric(neighborCandidate.outerCircle, neighborCandidate.innerCircle);
        }
        break;
    }

    if (foundValidNeighbor) {
      return neighborCandidate;
    }
  }
  if (!foundValidNeighbor) {
    // We hit maxNeighborIterations and couldn't find a good next neighbor cantidate so we must terminate.
    return null;
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

/**
 * Calculates the energy difference between two solutions.
 *
 * @param {Object} newSolution - The candidate solution.
 * @param {Object} oldSolution - The current solution.
 * @param {string} circleType - The type of circle optimization ("MIC", "MCC", or "MZC").
 * @param {Array<Array<number>>} points - The set of 2D points.
 * @returns {number} - The energy difference between the two solutions.
 */
function calculateEnergyDifference(newSolution, oldSolution, circleType, points) {
  switch (circleType) {
    case "MIC": {
      // If the radius is larger, the energy should be negative (meaning it's a better solution). Also if the center is farther away from the nearest point, the energy should be negative. This is because
      // being farther away from the nearest point means we can expand the radius (on the next step).
      const nearestPointDistance = findNearestPoint(points, oldSolution.center).distance - findNearestPoint(points, newSolution.center).distance;
      return (-newSolution.radius - oldSolution.radius) + nearestPointDistance;
    }
    case "MCC": {
      // If the radius is smaller, energy is negative. If the center is farther away from the farthest point, the energy should be negative. This is because being farther away from the farthest point means
      // we can shrink the radius (on the next step).
      const farthestPointDistance = findFarthestPoint(points, oldSolution.center).distance - findFarthestPoint(points, newSolution.center).distance;
      return (newSolution.radius - oldSolution.radius) + farthestPointDistance;
    }
    case "MZC":
      // In this case we use the center distance energy calculation of the MIC for the inner circle and the MCC for the outer circle and then the difference of the radii for the rest of the energy calculation.
      // If the radius got smaller, the energy is negative.
      const nearestPointDistance = findNearestPoint(points, oldSolution.innerCircle.center).distance - findNearestPoint(points, newSolution.innerCircle.center).distance;
      const farthestPointDistance = findFarthestPoint(points, oldSolution.outerCircle.center).distance - findFarthestPoint(points, newSolution.outerCircle.center).distance;
      const newRadiusDiff = newSolution.outerCircle.radius - newSolution.innerCircle.radius;
      const oldRadiusDiff = oldSolution.outerCircle.radius - oldSolution.innerCircle.radius;
      return (newRadiusDiff - oldRadiusDiff) + nearestPointDistance + farthestPointDistance;
  }
}

/**
 * Evaluates the quality of a solution based on the optimization criteria.
 *
 * @param {Object} solution - The solution to evaluate.
 * @param {string} circleType - The type of circle optimization ("MIC", "MCC", or "MZC").
 * @returns {number} - The evaluation score of the solution.
 */
function evaluateSolution(solution, circleType) {
  switch (circleType) {
    case "MIC":
      return solution.radius;
    case "MCC":
      return -solution.radius;
    case "MZC":
      return -(solution.outerCircle.radius - solution.innerCircle.radius);
  }
}

export { simulatedAnnealing };
