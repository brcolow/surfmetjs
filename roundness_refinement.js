import {
  computeConvexHull,
  distanceSquared,
  distanceToSegmentSquared,
  isPointInPolygon
} from './utils.js';

const DEFAULT_REFINEMENT_OPTIONS = {
  directionCount: 64,
  initialStepDivisor: 50,
  minStep: 1e-8,
  maxIterations: 1000
};

export function refineRoundnessCenter(points, circleType, initialCenter, options = {}) {
  const {
    directionCount,
    initialStepDivisor,
    minStep,
    maxIterations
  } = {
    ...DEFAULT_REFINEMENT_OPTIONS,
    ...options
  };

  if (!initialCenter || !Number.isFinite(initialCenter[0]) || !Number.isFinite(initialCenter[1])) {
    return initialCenter;
  }

  const hull = computeConvexHull(points);
  const directions = compassDirections(directionCount);
  let center = [...initialCenter];
  let objective = exactRoundnessObjective(points, circleType, center, hull);
  let stepSize = boundingBoxDiagonal(points) / initialStepDivisor;

  if (!Number.isFinite(objective) || stepSize <= 0) {
    return center;
  }

  for (let iteration = 0; iteration < maxIterations && stepSize > minStep; iteration++) {
    let bestCandidate = center;
    let bestObjective = objective;

    for (const [dx, dy] of directions) {
      const candidate = [
        center[0] + dx * stepSize,
        center[1] + dy * stepSize
      ];
      const candidateObjective = exactRoundnessObjective(points, circleType, candidate, hull);

      if (candidateObjective < bestObjective) {
        bestCandidate = candidate;
        bestObjective = candidateObjective;
      }
    }

    if (bestCandidate !== center) {
      center = bestCandidate;
      objective = bestObjective;
    } else {
      stepSize *= 0.5;
    }
  }

  return center;
}

export function roundnessResultFromCenter(points, circleType, center) {
  const radii = distancesFromCenter(points, center);
  const innerRadius = Math.min(...radii);
  const outerRadius = Math.max(...radii);

  switch (circleType) {
    case "MIC": {
      const radius = boundaryClearance(points, center);
      return {
        center,
        radius,
        roundness: outerRadius - radius,
        domain: "POLYGON_INTERIOR",
        profileReconstruction: "CHORDAL_POLYGON_FROM_ORDERED_SAMPLES",
        circleType: "MIC"
      };
    }

    case "MCC":
      return {
        center,
        radius: outerRadius,
        roundness: outerRadius - innerRadius,
        domain: "POINT_SET_EQUIVALENT_TO_POLYGON",
        circleType: "MCC"
      };

    case "MZC": {
      const zoneWidth = outerRadius - innerRadius;
      return {
        center,
        outerRadius,
        innerRadius,
        zoneWidth,
        roundness: zoneWidth / 2,
        innerCircle: { center, radius: innerRadius },
        outerCircle: { center, radius: outerRadius },
        domain: "POINT_SET",
        objective: "RADIAL_WIDTH",
        circleType: "MZC"
      };
    }

    default:
      throw new Error("Invalid circleType");
  }
}

function exactRoundnessObjective(points, circleType, center, hull) {
  if (!isAllowedCenter(points, circleType, center, hull)) {
    return Infinity;
  }

  const result = roundnessResultFromCenter(points, circleType, center);

  switch (circleType) {
    case "MIC":
      return -result.radius;
    case "MCC":
      return result.radius;
    case "MZC":
      return result.zoneWidth;
    default:
      throw new Error("Invalid circleType");
  }
}

function isAllowedCenter(points, circleType, center, hull) {
  const domain = circleType === "MIC" ? points : hull;
  return isPointInPolygon(center, domain);
}

function compassDirections(directionCount) {
  return Array.from({ length: directionCount }, (_, index) => {
    const angle = 2 * Math.PI * index / directionCount;
    return [Math.cos(angle), Math.sin(angle)];
  });
}

function distancesFromCenter(points, center) {
  return points.map(point => Math.sqrt(distanceSquared(point, center)));
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

function boundingBoxDiagonal(points) {
  const xs = points.map(([x]) => x);
  const ys = points.map(([, y]) => y);
  return Math.hypot(
    Math.max(...xs) - Math.min(...xs),
    Math.max(...ys) - Math.min(...ys)
  );
}
