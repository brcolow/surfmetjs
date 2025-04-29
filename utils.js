/**
 * Calculates the squared distance between two points.
 *
 * @param {Array<number>} a - The first point [x, y].
 * @param {Array<number>} b - The second point [x, y].
 * @returns {number} - The squared distance between the points.
 */
export function distanceSquared(a, b) {
  return (a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]);
}

/**
 * Generates a random number between the given range.
 *
 * @param {number} min - The minimum value (inclusive).
 * @param {number} max - The maximum value (exclusive).
 * @returns {number} - A random number in the specified range.
 */
export function getRandomBetween(min, max) {
  return Math.random() * (max - min) + min;
}

/**
 * Checks whether two circles are concentric.
 *
 * @param {Object} circle1 - The first circle, with `center` and `radius`.
 * @param {Object} circle2 - The second circle, with `center` and `radius`.
 * @returns {boolean} - `true` if the circles are concentric, `false` otherwise.
 */
export function areConcentric(circle1, circle2) {
  const dx = circle1.center[0] - circle2.center[0];
  const dy = circle1.center[1] - circle2.center[1];

  // Check if the centers are the same.
  return dx === 0 && dy === 0;
}

/**
 * Computes the mean of a set of points.
 *
 * @param {Array<Array<number>>} points - The set of 2D points.
 * @returns {Array<number>} - The mean of the points
*/
export function mean(points) {
  let mean = points.reduce(([x, y], p) => [x + p[0], y + p[1]], [0, 0]);
  mean[0] /= points.length;
  mean[1] /= points.length;
  return mean;
}

/**
 * Apply the transformation p |--> (p - center) / ||p - center||^2
 *
 * @param {Array<Array<number>>} points - The set of 2D points to apply the transformation to
 * @param {Array<number>} center - The center of the transformation
 * @returns {Array<Array<number>>} - The transformed points
 */
export function invertAtUnitCircle(points, center) {
  return points.map(p => {
    let pointsOffset = [p[0] - center[0], p[1] - center[1]];
    let distSquared = distanceSquared(p, center);
    return [
      center[0] + pointsOffset[0] / distSquared,
      center[1] + pointsOffset[1] / distSquared,
    ]
  });
}

/**
 * Computes the convex hull of a set of 2D points using the Graham scan algorithm.
 *
 * @param {Array<Array<number>>} points - The set of 2D points.
 * @returns {Array<Array<number>>} - The vertices of the convex hull in counterclockwise order.
 */
export function computeConvexHull(points) {
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
 * Generates an initial solution for the circle optimization problem.
 *
 * @param {Array<Array<number>>} points - The set of 2D points.
 * @param {string} circleType - The type of circle optimization ("MIC", "MCC", or "MZC").
 * @param {Object} leastSquaresCircle - The least-squares circle, with properties `a` and `b` for its center.
 * @returns {Object} - The initial solution, with format depending on `circleType`.
 */
export function getInitialSolution(points, circleType, leastSquaresCircle) {
  if (circleType === "MIC") {
    const nearestPointDistance = findNearestPoint(points, [leastSquaresCircle.a, leastSquaresCircle.b]).distance;
    const maxRadius = Math.sqrt(Math.max(...points.map(p => distanceSquared(p, [leastSquaresCircle.a, leastSquaresCircle.b]))));
    const minRadius = Math.sqrt(Math.min(...points.map(p => distanceSquared(p, [leastSquaresCircle.a, leastSquaresCircle.b]))));

    // We want to use the LSC center plus or minus a small amount to randomize the starting center. Go until we find one that is inside the convex hull of the point cloud.
    let center = [leastSquaresCircle.a + getRandomBetween(-(maxRadius - minRadius) / 10000, (maxRadius - minRadius) / 10000), leastSquaresCircle.b + getRandomBetween(-(maxRadius - minRadius) / 10000, (maxRadius - minRadius) / 10000)];
    const convexHull = computeConvexHull(points);
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
    // For MCC (Minimum Circumscribed Circle), we want an initial circle that covers all points.
    // A simple strategy is:
    // 1. Find the axis-aligned bounding box (AABB) of the points.
    // 2. Use the center of the bounding box as the initial circle center.
    // 3. Set the radius to the maximum distance from this center to any point (plus a small epsilon).
    //
    // This guarantees:
    // - All points are initially enclosed.
    // - The radius is finite and reasonably tight (no infinity or absurdly large values).
    // - A good neutral starting point for gradient descent optimization.
    const xs = points.map(p => p[0]);
    const ys = points.map(p => p[1]);
    const minX = Math.min(...xs);
    const maxX = Math.max(...xs);
    const minY = Math.min(...ys);
    const maxY = Math.max(...ys);

    const centerX = (minX + maxX) / 2;
    const centerY = (minY + maxY) / 2;

    const maxDist = Math.max(...points.map(([px, py]) => Math.hypot(px - centerX, py - centerY)));

    return { radius: maxDist - 1e-6, center: [centerX, centerY] };
  } else if (circleType === "MZC") {
    // Balance outer and inner points influence
    // Compute outer convex hull
    const outerConvexHull = computeConvexHull(points);

    // Compute weighted mean: outer hull points (weight 2), inner points (weight 1)
    const weights = new Array(points.length).fill(1);
    const outerPointSet = new Set(outerConvexHull.map(pt => pt.toString()));
    points.forEach((pt, i) => {
      if (outerPointSet.has(pt.toString())) {
        weights[i] = 2;
      }
    });

    const weightedMean = meanWeighted(points, weights);
    // Invert all points around weighted mean with normalized inversion
    const pointsInverted = invertAtUnitCircleNormalized(points, weightedMean);
    // Compute convex hull of inverted points
    const invertedInnerHull = computeConvexHull(pointsInverted);
    // Invert back to original coordinate system
    const innerConvexHull = invertAtUnitCircleNormalized(invertedInnerHull, weightedMean);
    const meanInnerPoints = mean(innerConvexHull);

    return {
      outerCircle: {
        center: weightedMean,
        radius: findFarthestPoint(points, weightedMean).distance
      },
      innerCircle: {
        center: meanInnerPoints,
        radius: findNearestPoint(points, meanInnerPoints).distance
      }
    };
  }
}

/**
 * Computes the weighted mean (centroid) of a set of 2D points.
 * 
 * @param {Array<Array<number>>} points - Array of points [[x1, y1], [x2, y2], ...].
 * @param {Array<number>} weights - Array of weights corresponding to each point.
 * @returns {Array<number>} The weighted mean [x, y].
 */
function meanWeighted(points, weights) {
  let sumX = 0, sumY = 0, sumW = 0;
  for (let i = 0; i < points.length; i++) {
    sumX += points[i][0] * weights[i];
    sumY += points[i][1] * weights[i];
    sumW += weights[i];
  }
  return [sumX / sumW, sumY / sumW];
}

/**
 * Inverts a set of 2D points relative to a given center using normalized inversion.
 * 
 * Each point is mapped such that points closer to the center are pushed farther away,
 * and points farther from the center are pulled closer, using a 1/distance² scaling.
 * 
 * @param {Array<Array<number>>} points - Array of points [[x1, y1], [x2, y2], ...].
 * @param {Array<number>} center - The center [x, y] to invert points around.
 * @returns {Array<Array<number>>} Array of inverted points.
 */
function invertAtUnitCircleNormalized(points, center) {
  return points.map(p => {
    const dx = p[0] - center[0];
    const dy = p[1] - center[1];
    const dist = Math.hypot(dx, dy);
    if (dist === 0) return [center[0], center[1]]; // Avoid division by zero
    const factor = 1 / (dist * dist); // Normalized inversion
    return [center[0] + dx * factor, center[1] + dy * factor];
  });
}

/**
 * Checks if no points lie strictly inside the given circle.
 *
 * @param {Array<Array<number>>} points - The set of 2D points.
 * @param {Array<number>} center - The center of the circle [x, y].
 * @param {number} radius - The radius of the circle.
 * @returns {boolean} - `true` if no points are strictly inside, `false` otherwise.
 */
function noPointsStrictlyInside(points, center, radius) {
  for (const point of points) {
    const distance = distanceSquared(point, center);
    if (distance < radius * radius) {
      return false;
    }
  }

  return true;
}

/**
 * Determines the orientation of three points.
 *
 * @param {Array<number>} p - The first point [x, y].
 * @param {Array<number>} q - The second point [x, y].
 * @param {Array<number>} r - The third point [x, y].
 * @returns {number} - 0 for collinear, 1 for clockwise, 2 for counterclockwise.
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
 * Finds the nearest point to the given reference point from a set of points.
 *
 * @param {Array<Array<number>>} points - The set of points to search through.
 * @param {Array<number>} point - The reference point [x, y].
 * @returns {Object} - The nearest point and its distance: `{ point, distance }`.
 */
export function findNearestPoint(points, point) {
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

/**
 * Finds the farthest point from the given reference point in a set of points.
 *
 * @param {Array<Array<number>>} points - The set of points to search through.
 * @param {Array<number>} point - The reference point [x, y].
 * @returns {Object} - The farthest point and its distance: `{ point, distance }`.
 */
export function findFarthestPoint(points, point) {
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

/**
 * Checks if a point is inside a polygon using the ray casting algorithm.
 *
 * @param {Array<number>} point - The point to check [x, y].
 * @param {Array<Array<number>>} polygon - The vertices of the polygon in order.
 * @returns {boolean} - `true` if the point is inside the polygon, `false` otherwise.
 */
export function isPointInPolygon(point, polygon) {
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
