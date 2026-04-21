const GEOMETRY_EPSILON = 1e-12;

/**
 * Calculates the squared distance between two points.
 *
 * @param {Array<number>} a - The first point [x, y].
 * @param {Array<number>} b - The second point [x, y].
 * @returns {number} The squared distance between the points.
 */
export function distanceSquared(a, b) {
  const dx = a[0] - b[0];
  const dy = a[1] - b[1];
  return dx * dx + dy * dy;
}

/**
 * Generates a random number between the given range.
 *
 * @param {number} min - The minimum value (inclusive).
 * @param {number} max - The maximum value (exclusive).
 * @returns {number} A random number in the specified range.
 */
export function getRandomBetween(min, max) {
  return Math.random() * (max - min) + min;
}

/**
 * Checks whether two circles are concentric.
 *
 * @param {Object} circle1 - The first circle, with `center`.
 * @param {Object} circle2 - The second circle, with `center`.
 * @param {number} [epsilon=1e-12] - Center comparison tolerance.
 * @returns {boolean} True if the circles are concentric.
 */
export function areConcentric(circle1, circle2, epsilon = GEOMETRY_EPSILON) {
  return distanceSquared(circle1.center, circle2.center) <= epsilon * epsilon;
}

/**
 * Computes the mean of a set of points.
 *
 * @param {Array<Array<number>>} points - The set of 2D points.
 * @returns {Array<number>} The mean of the points.
 */
export function mean(points) {
  assertNonEmptyPoints(points);

  const [sumX, sumY] = points.reduce(
    ([x, y], point) => [x + point[0], y + point[1]],
    [0, 0]
  );

  return [sumX / points.length, sumY / points.length];
}

/**
 * Applies the inversion p -> center + (p - center) / ||p - center||^2.
 *
 * @param {Array<Array<number>>} points - Points to invert.
 * @param {Array<number>} center - Inversion center.
 * @returns {Array<Array<number>>} Inverted points.
 */
export function invertAtUnitCircle(points, center) {
  return points.map(point => invertPointAtUnitCircle(point, center));
}

function invertPointAtUnitCircle(point, center) {
  const dx = point[0] - center[0];
  const dy = point[1] - center[1];
  const radiusSquared = dx * dx + dy * dy;

  if (radiusSquared <= GEOMETRY_EPSILON) {
    return [center[0], center[1]];
  }

  return [
    center[0] + dx / radiusSquared,
    center[1] + dy / radiusSquared
  ];
}

/**
 * Computes the convex hull of a set of 2D points with Andrew's monotone chain
 * algorithm. The input array is not mutated.
 *
 * Collinear interior points are omitted. If all points are collinear, the two
 * endpoints are returned.
 *
 * @param {Array<Array<number>>} points - The set of 2D points.
 * @returns {Array<Array<number>>} Hull vertices in counterclockwise order.
 */
export function computeConvexHull(points) {
  if (!Array.isArray(points) || points.length === 0) {
    return [];
  }

  const sorted = uniquePoints(points)
    .sort((a, b) => a[0] - b[0] || a[1] - b[1]);

  if (sorted.length <= 2) {
    return sorted.map(point => [...point]);
  }

  const lower = [];
  for (const point of sorted) {
    while (lower.length >= 2 && cross(lower[lower.length - 2], lower[lower.length - 1], point) <= 0) {
      lower.pop();
    }
    lower.push(point);
  }

  const upper = [];
  for (let i = sorted.length - 1; i >= 0; i--) {
    const point = sorted[i];
    while (upper.length >= 2 && cross(upper[upper.length - 2], upper[upper.length - 1], point) <= 0) {
      upper.pop();
    }
    upper.push(point);
  }

  lower.pop();
  upper.pop();
  return lower.concat(upper).map(point => [...point]);
}

function uniquePoints(points) {
  const seen = new Set();
  const result = [];

  for (const point of points) {
    const key = `${point[0]},${point[1]}`;
    if (!seen.has(key)) {
      seen.add(key);
      result.push([point[0], point[1]]);
    }
  }

  return result;
}

function cross(origin, a, b) {
  return (a[0] - origin[0]) * (b[1] - origin[1]) -
    (a[1] - origin[1]) * (b[0] - origin[0]);
}

/**
 * Generates an initial solution for the circle optimization problem.
 *
 * This helper is intentionally deterministic and finite. It is a starting
 * point generator, not an exact MIC/MCC/MZC solver.
 *
 * @param {Array<Array<number>>} points - The set of 2D points.
 * @param {string} circleType - "MIC", "MCC", or "MZC".
 * @param {Object} leastSquaresCircle - Optional estimate with `a`, `b`, and `r`.
 * @returns {Object} Initial solution with format depending on `circleType`.
 */
export function getInitialSolution(points, circleType, leastSquaresCircle = {}) {
  assertNonEmptyPoints(points);

  const hull = computeConvexHull(points);
  const centerDomain = circleType === "MIC" ? points : hull;
  const candidates = initialCenterCandidates(points, leastSquaresCircle)
    .filter(center => centerDomain.length < 3 || isPointInPolygon(center, centerDomain));
  const centers = candidates.length > 0 ? candidates : [mean(points)];

  switch (circleType) {
    case "MIC":
      return centers
        .map(center => ({
          center,
          radius: boundaryClearance(points, center),
          domain: "POLYGON_INTERIOR",
          profileReconstruction: "CHORDAL_POLYGON_FROM_ORDERED_SAMPLES"
        }))
        .reduce((best, candidate) => candidate.radius > best.radius ? candidate : best);

    case "MCC": {
      const center = boundingBoxCenter(points);
      return { center, radius: maxDistance(points, center) };
    }

    case "MZC": {
      const center = centers
        .map(candidate => {
          const innerRadius = minDistance(points, candidate);
          const outerRadius = maxDistance(points, candidate);
          return { center: candidate, innerRadius, outerRadius };
        })
        .reduce((best, candidate) => {
          return candidate.outerRadius - candidate.innerRadius < best.outerRadius - best.innerRadius
            ? candidate
            : best;
        });

      return {
        innerCircle: { center: center.center, radius: center.innerRadius },
        outerCircle: { center: center.center, radius: center.outerRadius }
      };
    }

    default:
      throw new Error("Invalid circleType");
  }
}

function initialCenterCandidates(points, estimate) {
  const candidates = [
    mean(points),
    boundingBoxCenter(points)
  ];

  if (Number.isFinite(estimate.a) && Number.isFinite(estimate.b)) {
    candidates.unshift([estimate.a, estimate.b]);
  }

  if (estimate.center && Number.isFinite(estimate.center[0]) && Number.isFinite(estimate.center[1])) {
    candidates.unshift([estimate.center[0], estimate.center[1]]);
  }

  return uniquePoints(candidates);
}

function boundingBoxCenter(points) {
  const xs = points.map(([x]) => x);
  const ys = points.map(([, y]) => y);

  return [
    (Math.min(...xs) + Math.max(...xs)) / 2,
    (Math.min(...ys) + Math.max(...ys)) / 2
  ];
}

function minDistance(points, center) {
  return findNearestPoint(points, center).distance;
}

function maxDistance(points, center) {
  return findFarthestPoint(points, center).distance;
}

function boundaryClearance(points, center) {
  let clearanceSquared = Infinity;
  for (let i = 0, j = points.length - 1; i < points.length; j = i++) {
    clearanceSquared = Math.min(
      clearanceSquared,
      distanceToSegmentSquared(center, points[j], points[i])
    );
  }
  return Math.sqrt(clearanceSquared);
}

/**
 * Finds the nearest point to a reference point.
 *
 * @param {Array<Array<number>>} points - Points to search.
 * @param {Array<number>} point - Reference point [x, y].
 * @returns {{point: Array<number>, distance: number}} Nearest point and distance.
 */
export function findNearestPoint(points, point) {
  assertNonEmptyPoints(points);

  let minDistanceSquared = Infinity;
  let nearestPoint = null;

  for (const candidate of points) {
    const candidateDistance = distanceSquared(candidate, point);
    if (candidateDistance < minDistanceSquared) {
      minDistanceSquared = candidateDistance;
      nearestPoint = candidate;
    }
  }

  return { point: nearestPoint, distance: Math.sqrt(minDistanceSquared) };
}

/**
 * Finds the farthest point from a reference point.
 *
 * @param {Array<Array<number>>} points - Points to search.
 * @param {Array<number>} point - Reference point [x, y].
 * @returns {{point: Array<number>, distance: number}} Farthest point and distance.
 */
export function findFarthestPoint(points, point) {
  assertNonEmptyPoints(points);

  let maxDistanceSquared = -Infinity;
  let farthestPoint = null;

  for (const candidate of points) {
    const candidateDistance = distanceSquared(candidate, point);
    if (candidateDistance > maxDistanceSquared) {
      maxDistanceSquared = candidateDistance;
      farthestPoint = candidate;
    }
  }

  return { point: farthestPoint, distance: Math.sqrt(maxDistanceSquared) };
}

/**
 * Checks if a point is inside a polygon. Points on the boundary are considered
 * inside by default.
 *
 * @param {Array<number>} point - The point to check [x, y].
 * @param {Array<Array<number>>} polygon - Polygon vertices.
 * @param {boolean} [includeBoundary=true] - Whether boundary points count.
 * @returns {boolean} True if the point is inside the polygon.
 */
export function isPointInPolygon(point, polygon, includeBoundary = true) {
  if (!Array.isArray(polygon) || polygon.length === 0) {
    return false;
  }

  if (includeBoundary && isPointOnPolygonBoundary(point, polygon)) {
    return true;
  }

  if (polygon.length < 3) {
    return false;
  }

  let inside = false;
  for (let i = 0, j = polygon.length - 1; i < polygon.length; j = i++) {
    const xi = polygon[i][0];
    const yi = polygon[i][1];
    const xj = polygon[j][0];
    const yj = polygon[j][1];
    const intersects = (yi > point[1]) !== (yj > point[1]) &&
      point[0] < ((xj - xi) * (point[1] - yi)) / (yj - yi) + xi;

    if (intersects) {
      inside = !inside;
    }
  }

  return inside;
}

export function isPointOnPolygonBoundary(point, polygon, epsilon = GEOMETRY_EPSILON) {
  for (let i = 0, j = polygon.length - 1; i < polygon.length; j = i++) {
    if (distanceToSegmentSquared(point, polygon[j], polygon[i]) <= epsilon * epsilon) {
      return true;
    }
  }

  return false;
}

export function distanceToSegmentSquared(point, start, end) {
  const dx = end[0] - start[0];
  const dy = end[1] - start[1];
  const lengthSquared = dx * dx + dy * dy;

  if (lengthSquared <= GEOMETRY_EPSILON) {
    return distanceSquared(point, start);
  }

  const t = Math.max(
    0,
    Math.min(1, ((point[0] - start[0]) * dx + (point[1] - start[1]) * dy) / lengthSquared)
  );
  const projection = [start[0] + t * dx, start[1] + t * dy];

  return distanceSquared(point, projection);
}

function assertNonEmptyPoints(points) {
  if (!Array.isArray(points) || points.length === 0) {
    throw new Error("At least one point is required.");
  }
}
