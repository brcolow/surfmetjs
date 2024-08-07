function simulatedAnnealing(points, circleType, initialTemperature, coolingType, maxIterations, maxNeighborIterations, stepSize) {
  let currentSolution = getInitialSolution(points, circleType);
  let bestSolution = currentSolution;

  for (let i = 0; i < maxIterations; i++) {
    const temperature = temperatureSchedule(initialTemperature, i, coolingType);

    // We could repeat the following until equilibrium is approached sufficiently closely (for this temperature value - would need to define that carefully.)
    const neighbor = generateNeighbor(points, circleType, currentSolution, maxNeighborIterations, stepSize);
    if (neighbor === null) {
      console.log("Couldn't find next neighbor, at iteration: " + i);
      break;
    }
    const deltaEnergy = calculateEnergyDifference(neighbor, currentSolution, circleType, points);

    // Metropolis criterion
    if (deltaEnergy < 0 || Math.exp(-deltaEnergy / temperature) > Math.random()) {
      currentSolution = neighbor;
    }

    if (evaluateSolution(bestSolution, circleType) > evaluateSolution(currentSolution, circleType)) {
      bestSolution = currentSolution;
    }

    if (i === maxIterations - 1) {
      console.log("Made it to last iteration: " + i);
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

function getInitialSolution(points, circleType) {
  switch (circleType) {
    case "MIC":
      // Find the point closest to the origin and use that as the radius of a circle centered at (0, 0) for the initial solution.
      return { radius: findNearestPoint(points, [0, 0]).distance, center: [0, 0] };
    case "MCC":
      // Find the point farthest from the origin and use that as the radius of a circle centered at (0, 0) for the initial solution.
      return { radius: findFarthestPoint(points, [0, 0]).distance, center: [0, 0] };
    case "MZC":
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
      // FIXME: It doesn't make sense to only change the center because then this neighbor will never have less energy, unless we can include in the energy function
      // a measure of how close the center is to the centroid of all the points? Something like that (what would be more likely for the MIC to be...) bascially
      // if we move the center and the new center is farther away from every point (by some measure) then it's a good one to keep as the next solution because
      // then the chance that the radius can be increased next time is higher. For the MCC we would rate centers that are closer to every point higher.
      if (getRandomBetween(0, 100) > 50) {
        neighborCandidate = {
          center: [currentSolution.center[0] += getRandomBetween(-stepSize, stepSize), currentSolution.center[1] += getRandomBetween(-stepSize, stepSize)],
          radius: currentSolution.radius
        };
      } else {
        neighborCandidate = {
          center: currentSolution.center,
          radius: currentSolution.radius += getRandomBetween(-stepSize, stepSize)
        };
      }

    } else {
      // When finding the neighbor for MZC we want to share the same center because concentricity must be maintained.
      if (getRandomBetween(0, 100) > 50) {
        let center = [currentSolution.center[0] += getRandomBetween(-stepSize, stepSize), currentSolution.center[1] += getRandomBetween(-stepSize, stepSize)];
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
        neighborCandidate = {
          outerCircle: {
            center: currentSolution.center,
            radius: currentSolution.radius += getRandomBetween(-stepSize, stepSize)
          }, innerCircle: {
            center: currentSolution.center,
            radius: currentSolution.radius += getRandomBetween(-stepSize, stepSize)
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

export { simulatedAnnealing };
