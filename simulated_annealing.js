function simulatedAnnealing(points, circleType, leastSquaresCircle, initialTemperature, coolingType, maxIterations, maxNeighborIterations, stepSize) {
  const convexHull = computeConvexHull(points);
  const initialSolution = getInitialSolution(points, "MIC", leastSquaresCircle, convexHull);
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

function getInitialSolution(points, circleType, leastSquaresCircle, convexHull) {
  if (circleType === "MIC") {
    const nearestPointDistance = findNearestPoint(points, [leastSquaresCircle.a, leastSquaresCircle.b]).distance;
    const maxRadius = Math.sqrt(Math.max(...points.map(p => distanceSquared(p, [leastSquaresCircle.a, leastSquaresCircle.b]))));
    const minRadius = Math.sqrt(Math.min(...points.map(p => distanceSquared(p, [leastSquaresCircle.a, leastSquaresCircle.b]))));

    // We want to use the LSC center plus or minus a small amount to randomize the starting center. Go until we find one that is inside the convex hull of the point cloud.
    let center = [leastSquaresCircle.a + getRandomBetween(-(maxRadius - minRadius) / 10000, (maxRadius - minRadius) / 10000), leastSquaresCircle.b + getRandomBetween(-(maxRadius - minRadius) / 10000, (maxRadius - minRadius) / 10000)];
    while (!isPointInPolygon(center, convexHull)) {
      center = [leastSquaresCircle.a + getRandomBetween(-(maxRadius - minRadius) / 10000, (maxRadius - minRadius) / 10000), leastSquaresCircle.b + getRandomBetween(-(maxRadius - minRadius) / 10000, (maxRadius - minRadius) / 10000)];
    }

    // We want to use the radius corresponding to the distance from the LSC center to the nearest point in the point cloud, plus or minus a small amount. Go until we find one such that no points are (strictly) inside the initial solution circle.
    let radius = nearestPointDistance + getRandomBetween(-(maxRadius - minRadius) / 10000, 0);
    let tries = 0;
    while (!noPointsStrictlyInside(points, center, radius)) {
      radius = nearestPointDistance + getRandomBetween(-(maxRadius - minRadius) / 10000, 0);
      tries++;

      if (tries > 1000) {
        console.log("We could not find a suitable radius for the center even after 1000 tries!");
        return { radius: nearestPointDistance - 0.00001, center: [leastSquaresCircle.a, leastSquaresCircle.b] };
      }
    }
    
    return { radius: radius, center: [center[0], center[1]] };
  } else if (circleType === "MCC") {
    // Find the point farthest from the origin and use that as the radius of a circle centered at (0, 0) for the initial solution.
    return { radius: findFarthestPoint(points, [0, 0]).distance, center: [0, 0] };
  } else if (circleType === "MZC") {
    // Use the MIC initial guess as the inner circle and the MCC guess as the outer circle.
    return { outerCircle: { radius: findFarthestPoint(points, [0, 0]).distance, center: [0, 0] }, innerCircle: { radius: findNearestPoint(points, [0, 0]).distance, center: [0, 0] } };
  }
}

function noPointsStrictlyInside(points, center, radius) {
  for (const point of points) {
    const distance = distanceSquared(point, center);
    if (distance < radius * radius) {
      return false;
    }
  }

  return true;
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

/**
 * Computes and returns the convex hull of the given points using the Graham scan.
 */
function computeConvexHull(points) {
  // Find the lowest point
  let minY = Infinity;
  let minIndex = 0;
  for (let i = 0; i < points.length; i++) {
    if (points[i][1] < minY || (points[i][1] === minY && points[i][0] < points[minIndex][0])) {
      minY = points[i][1];
      minIndex = i;
    }
  }

  // Swap the lowest point to the first position
  [points[0], points[minIndex]] = [points[minIndex], points[0]];

  // Sort points by polar angle
  // Note: We may be able to skip this because the circlular trace is always ordered this way...
  points.sort((a, b) => {
    const angleA = Math.atan2(a[1] - points[0][1], a[0] - points[0][0]);
    const angleB = Math.atan2(b[1] - points[0][1], b[0] - points[0][0]);
    return angleA - angleB;
  });

  const stack = [];
  stack.push(points[0]);
  stack.push(points[1]);

  // Build the convex hull.
  for (let i = 2; i < points.length; i++) {
    let top = stack.pop();
    while (orientation(stack[stack.length - 1], top, points[i]) <= 0) {
      top = stack.pop();
    }
    stack.push(top);
    stack.push(points[i]);
  }

  return stack;
}

/**
 * Takes an ordered triplet and returns 0, 1, or 2 depending on if the points are orientated colinearly, clockwise, or counter-clockwise respectively.
 */
function orientation(p, q, r) {
  // To find orientation of ordered triplet (p, q, r).
  // The function returns following values
  // 0 : Colinear points
  // 1 : Clockwise points
  // 2 : Counterclockwise points
  const val = (q[1] - p[1]) * (r[0] - q[0]) - (q[0] - p[0]) * (r[1] - q[1]);
  if (val === 0) {
    return 0;  // colinear
  } else {
    return (val > 0) ? 1 : 2; // clockwise or counterclock wise
  }
}

/**
 * Returns true if the point is inside the given polygon; false otherwise. Uses the ray casting algorithm (i.e. the crossing number algorithm).
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

export { simulatedAnnealing, getInitialSolution, computeConvexHull };
