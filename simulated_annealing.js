function simulatedAnnealing(points, circleType, initialSolution, initialTemperature, coolingType, maxIterations, maxNeighborIterations, stepSize) {
  let currentSolution = initialSolution;
  let bestSolution = initialSolution;

  for (let i = 0; i < maxIterations; i++) {
    const temperature = temperatureSchedule(initialTemperature, i, coolingType);

    // We could repeat the following until equilibrium is approached sufficiently closely (for this temperature value - would need to define that carefully.)
    const neighbor = generateNeighbor(points, circleType, currentSolution, maxNeighborIterations, stepSize);
    if (neighbor === null) {
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
  }

  return bestSolution;
}

function findNearestPoint(points, point) {
  let minDistance = Infinity;
  let nearestPoint = null;

  for (const p of points) {
    const distance = distanceSquared(p, point);
    if (distance < minDistance) {
      minDistance = distance;
      nearestPoint = p;
    }
  }
  return { point: nearestPoint, distance: Math.sqrt(minDistance) };
}

function findFarthestPoint(points, point) {
  let maxDistance = Infinity;
  let farthestPoint = null;

  for (const p of points) {
    const distance = distanceSquared(p, point);
    if (distance > maxDistance) {
      maxDistance = distance;
      farthestPoint = p;
    }
  }

  return { point: farthestPoint, distance: Math.sqrt(maxDistance) };
}

function getInitialSolution(points, circleType, leastSquaresCircle) {
  // TODO: We could experiment with using the LSC center as the initial solution.
  if (circleType === "MIC") {
    const nearestPointDistance = findNearestPoint(points, [leastSquaresCircle.a, leastSquaresCircle.b]).distance;
    const maxRadius = Math.sqrt(Math.max(...points.map(p => distanceSquared(p, [leastSquaresCircle.a, leastSquaresCircle.b]))));
    const minRadius = Math.sqrt(Math.min(...points.map(p => distanceSquared(p, [leastSquaresCircle.a, leastSquaresCircle.b]))));
    const offset = ((maxRadius - minRadius) / 100);
    // Find the point closest to the origin and use that as the radius of a circle centered at (0, 0) for the initial solution.
    return { radius: nearestPointDistance - offset, center: [leastSquaresCircle.a, leastSquaresCircle.b] };
  } else if (circleType === "MCC") {
    // Find the point farthest from the origin and use that as the radius of a circle centered at (0, 0) for the initial solution.
    return { radius: findFarthestPoint(points, [0, 0]).distance, center: [0, 0] };
  } else if (circleType === "MZC") {
    // Use the MIC initial guess as the inner circle and the MCC guess as the outer circle.
    return { outerCircle: { radius: findFarthestPoint(points, [0, 0]).distance, center: [0, 0] }, innerCircle: { radius: findNearestPoint(points, [0, 0]).distance, center: [0, 0] } };
  }
}

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

function generateNeighbor(points, circleType, currentSolution, maxNeighborIterations, stepSize) {
  let foundValidNeighbor = false;
  for (let i = 0; i < maxNeighborIterations; i++) {
    let neighborCandidate;
    if (circleType !== "MZC") {
      if (getRandomBetween(0, 100) > 50) {
        neighborCandidate = {
          center: [currentSolution.center[0] + getRandomBetween(-stepSize, stepSize), currentSolution.center[1] + getRandomBetween(-stepSize, stepSize)],
          radius: currentSolution.radius
        };
      } else {
        neighborCandidate = {
          center: [currentSolution.center[0], currentSolution.center[1]],
          radius: currentSolution.radius + getRandomBetween(-stepSize, stepSize)
        };
      }

    } else {
      // When finding the neighbor for MZC we want to share the same center because concentricity must be maintained.
      if (getRandomBetween(0, 100) > 50) {
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
        for (const point of points) {
          const distance = distanceSquared(point, neighborCandidate.center);
          if (distance < neighborCandidate.radius * neighborCandidate.radius) {
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

function getRandomBetween(min, max) {
  return Math.random() * (max - min) + min;
}

function distanceSquared(a, b) {
  return (a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]);
}

function areConcentric(circle1, circle2) {
  const dx = circle1.center[0] - circle2.center[0];
  const dy = circle1.center[1] - circle2.center[1];

  // Check if the centers are the same.
  return dx === 0 && dy === 0;
}

export { simulatedAnnealing, getInitialSolution };
