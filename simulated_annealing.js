import { areConcentric, distanceSquared, findFarthestPoint, findFarthestPoint, getInitialSolution, getRandomBetween } from './utils.js';

/**
 * Implements simulated annealing to optimize the parameters of one or more circles
 * to fit a given point set.
 *
 * @param {Array<Array<number>>} points - The set of 2D points to optimize the circle(s) for.
 * @param {string} circleType - The type of circle optimization ("MIC", "MCC", or "MZC").
 * @param {Object} leastSquaresCircle - The least-squares circle, with properties `a` and `b` for its center.
 * @param {number} initialTemperature - The starting temperature for the annealing process.
 * @param {Object} coolingType - The cooling schedule with properties `type` ("exponential", "linear", "logarithmic") and `rate`.
 * @param {number} maxIterations - The maximum number of iterations for annealing.
 * @param {number} maxNeighborIterations - The maximum attempts to generate a valid neighbor solution.
 * @param {number} stepSize - The maximum step size for generating neighbor solutions.
 * @returns {Object} - The best solution found during the process.
 */
function simulatedAnnealing(points, circleType, leastSquaresCircle, initialTemperature, coolingType, maxIterations, maxNeighborIterations, stepSize) {
  const convexHull = computeConvexHull(points);
  const initialSolution = getInitialSolution(points, circleType, leastSquaresCircle, convexHull);
  let currentSolution = initialSolution;
  let bestSolution = initialSolution;

  for (let i = 0; i < maxIterations; i++) {
    const temperature = temperatureSchedule(initialTemperature, i, coolingType);

    // We could repeat the following until equilibrium is approached sufficiently closely (for this temperature value - would need to define that carefully.)
    const neighbor = generateNeighbor(points, convexHull, circleType, currentSolution, maxNeighborIterations, stepSize);
    if (neighbor === null) {
      console.log("Could not find next neighbor at iteration: " + i);
      break;
    }
    const deltaEnergy = calculateEnergyDifference(neighbor, currentSolution, circleType, points);

    // Metropolis criterion
    if (deltaEnergy < 0 || Math.exp(-deltaEnergy / temperature) > Math.random()) {
      currentSolution = neighbor;
    }

    if (evaluateSolution(currentSolution, circleType) > evaluateSolution(bestSolution, circleType)) {
      bestSolution = currentSolution;
    }

    if (i + 1 === maxIterations) {
      console.log("Got to max iterations: " + i);
    }
  }

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
function generateNeighbor(points, convexHull, circleType, currentSolution, maxNeighborIterations, stepSize) {
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

/**
 * Checks if a point is inside a polygon using the ray casting algorithm.
 *
 * @param {Array<number>} point - The point to check [x, y].
 * @param {Array<Array<number>>} polygon - The vertices of the polygon in order.
 * @returns {boolean} - `true` if the point is inside the polygon, `false` otherwise.
 */
function isPointInPolygon(point, polygon) {
  let inside = false;
  for (let i = 0, j = polygon.length - 1; i < polygon.length; j = i++) {
    const xi = polygon[i][0], yi = polygon[i][1];
    const xj = polygon[j][0], yj = polygon[j][1];
    const intersect = ((yi > point[1]) != (yj > point[1]))
      && (point[0] < (xj - xi) * (point[1] - yi) / (yj - yi) + xi);
    if (intersect) {
      inside = !inside;
    }
  }
  return inside;
}

export { simulatedAnnealing };
