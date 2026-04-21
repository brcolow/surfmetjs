import { describe, expect, it } from 'vitest';
import fixtureData from './fixtures/roundness/roundness_profiles.json';
import { adaptiveGradientDescent } from '../gradient_descent.js';
import { simulatedAnnealing } from '../simulated_annealing.js';

const profiles = fixtureData.profiles;
const ORACLE_ERROR_TOLERANCE = 1e-8;
const GRADIENT_DESCENT_TOLERANCE = ORACLE_ERROR_TOLERANCE;
const SIMULATED_ANNEALING_STEP = 0.2;
const SIMULATED_ANNEALING_TOLERANCE = ORACLE_ERROR_TOLERANCE;
const SIMULATED_ANNEALING_OPTIONS = {
  initialTemperature: 0.02,
  coolingType: { type: 'exponential', rate: 0.995 },
  maxIterations: 2000,
  equilibriumSteps: 20,
  maxNeighborIterations: 2000
};

function radiiFromCenter(points, center) {
  return points.map(([x, y]) => Math.hypot(x - center[0], y - center[1]));
}

function zoneWidth(points, center) {
  const radii = radiiFromCenter(points, center);
  return Math.max(...radii) - Math.min(...radii);
}

function maxRadius(points, center) {
  return Math.max(...radiiFromCenter(points, center));
}

function polylineClearance(points, center) {
  let clearanceSquared = Infinity;

  for (let i = 0; i < points.length; i++) {
    const start = points[i];
    const end = points[(i + 1) % points.length];
    clearanceSquared = Math.min(
      clearanceSquared,
      distanceToSegmentSquared(center, start, end)
    );
  }

  return Math.sqrt(clearanceSquared);
}

function distanceToSegmentSquared(point, start, end) {
  const dx = end[0] - start[0];
  const dy = end[1] - start[1];
  const lengthSquared = dx * dx + dy * dy;

  if (lengthSquared === 0) {
    return distanceSquared(point, start);
  }

  const t = Math.max(
    0,
    Math.min(1, ((point[0] - start[0]) * dx + (point[1] - start[1]) * dy) / lengthSquared)
  );
  return distanceSquared(point, [start[0] + t * dx, start[1] + t * dy]);
}

function distanceSquared(a, b) {
  return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2;
}

function meanPoint(points) {
  const [sumX, sumY] = points.reduce(
    ([sx, sy], [x, y]) => [sx + x, sy + y],
    [0, 0]
  );
  return [sumX / points.length, sumY / points.length];
}

function seededRandom(seed) {
  return () => {
    seed += 0x6D2B79F5;
    let value = seed;
    value = Math.imul(value ^ value >>> 15, value | 1);
    value ^= value + Math.imul(value ^ value >>> 7, value | 61);
    return ((value ^ value >>> 14) >>> 0) / 4294967296;
  };
}

function seedFromString(text) {
  let hash = 2166136261;

  for (let i = 0; i < text.length; i++) {
    hash ^= text.charCodeAt(i);
    hash = Math.imul(hash, 16777619);
  }

  return hash >>> 0;
}

function transformedPoint(point, transform) {
  const cos = Math.cos(transform.rotation);
  const sin = Math.sin(transform.rotation);
  return [
    transform.translation[0] + transform.scale * (point[0] * cos - point[1] * sin),
    transform.translation[1] + transform.scale * (point[0] * sin + point[1] * cos)
  ];
}

describe('CGAL roundness fixtures', () => {
  it('records metrology-domain goldens with passing CGAL certificates', () => {
    expect(fixtureData.schemaVersion).toBe(1);
    expect(fixtureData.profileReconstruction).toBe('CHORDAL_POLYGON_FROM_ORDERED_SAMPLES');
    expect(profiles.length).toBeGreaterThan(0);

    for (const profile of profiles) {
      expect(profile.cgal.mic.valid).toBe(true);
      expect(profile.cgal.mic.domain).toBe('POLYGON_INTERIOR');
      expect(profile.cgal.mcc.valid).toBe(true);
      expect(profile.cgal.mcc.domain).toBe('POINT_SET_EQUIVALENT_TO_POLYGON');
      expect(profile.cgal.mzc.valid).toBe(true);
      expect(profile.cgal.mzc.domain).toBe('POINT_SET');
      expect(profile.cgal.mzc.objective).toBe('RADIAL_WIDTH');
      expect(profile.cgal.mzc.solverExactFor).toBe('FINITE_POINT_SET_RADIAL_WIDTH');
      expect(profile.cgal.certs.mic).toBe('OK');
      expect(profile.cgal.certs.mcc).toBe('OK');
      expect(profile.cgal.certs.mzc).toBe('OK');
    }
  });

  it.each(profiles.map(profile => [profile.name, profile]))(
    'matches the CGAL MIC radius against the reconstructed profile boundary: %s',
    (_name, profile) => {
      const radius = polylineClearance(profile.points, profile.cgal.mic.center);
      expect(radius).toBeCloseTo(profile.cgal.mic.radius, 10);
    }
  );

  it.each(profiles.map(profile => [profile.name, profile]))(
    'matches the CGAL MCC radius when evaluated on extracted profile points: %s',
    (_name, profile) => {
      const radius = maxRadius(profile.points, profile.cgal.mcc.center);
      expect(radius).toBeCloseTo(profile.cgal.mcc.radius, 10);
    }
  );

  it.each(profiles.map(profile => [profile.name, profile]))(
    'matches the CGAL MZC width when evaluated in JS at the oracle center: %s',
    (_name, profile) => {
      const width = zoneWidth(profile.points, profile.cgal.mzc.center);
      expect(width).toBeCloseTo(profile.cgal.mzc.width, 10);
    }
  );

  it.each(profiles.map(profile => [profile.name, profile]))(
    'keeps adaptive gradient descent close to the CGAL MIC oracle: %s',
    (_name, profile) => {
      const result = adaptiveGradientDescent(
        profile.points,
        'MIC',
        meanPoint(profile.points),
        {
          gamma: 800,
          initialStepSize: 0.2,
          learningRateDecay: 0.7,
          slopeDamper: 0.1,
          maxIterations: 2500,
          maxIterationsLineSearch: 100,
          epsilon: 1e-14
        }
      );

      expect(result.domain).toBe('POLYGON_INTERIOR');
      expect(result.profileReconstruction).toBe('CHORDAL_POLYGON_FROM_ORDERED_SAMPLES');
      expect(result.radius).toBeLessThanOrEqual(profile.cgal.mic.radius + ORACLE_ERROR_TOLERANCE);
      expect(profile.cgal.mic.radius - result.radius).toBeLessThan(GRADIENT_DESCENT_TOLERANCE);
    }
  );

  it.each(profiles.map(profile => [profile.name, profile]))(
    'keeps adaptive gradient descent close to the CGAL MCC oracle: %s',
    (_name, profile) => {
      const result = adaptiveGradientDescent(
        profile.points,
        'MCC',
        meanPoint(profile.points),
        {
          gamma: 800,
          initialStepSize: 0.2,
          learningRateDecay: 0.7,
          slopeDamper: 0.1,
          maxIterations: 2500,
          maxIterationsLineSearch: 100,
          epsilon: 1e-14
        }
      );

      expect(result.domain).toBe('POINT_SET_EQUIVALENT_TO_POLYGON');
      expect(result.radius).toBeGreaterThanOrEqual(profile.cgal.mcc.radius - ORACLE_ERROR_TOLERANCE);
      expect(result.radius - profile.cgal.mcc.radius).toBeLessThan(GRADIENT_DESCENT_TOLERANCE);
    }
  );

  it.each(profiles.map(profile => [profile.name, profile]))(
    'keeps adaptive gradient descent close to the CGAL MZC oracle: %s',
    (_name, profile) => {
      const result = adaptiveGradientDescent(
        profile.points,
        'MZC',
        meanPoint(profile.points),
        {
          gamma: 800,
          initialStepSize: 0.2,
          learningRateDecay: 0.7,
          slopeDamper: 0.1,
          maxIterations: 2500,
          maxIterationsLineSearch: 100,
          epsilon: 1e-14
        }
      );

      expect(result.zoneWidth).toBeGreaterThanOrEqual(profile.cgal.mzc.width - ORACLE_ERROR_TOLERANCE);
      expect(result.zoneWidth - profile.cgal.mzc.width).toBeLessThan(GRADIENT_DESCENT_TOLERANCE);
    }
  );

  it.each(profiles.map(profile => [profile.name, profile]))(
    'keeps simulated annealing close to the CGAL MIC oracle: %s',
    (_name, profile) => {
      const result = simulatedAnnealing(
        profile.points,
        'MIC',
        meanPoint(profile.points),
        SIMULATED_ANNEALING_STEP,
        {
          ...SIMULATED_ANNEALING_OPTIONS,
          random: seededRandom(seedFromString(`${profile.name}:MIC`))
        }
      );

      expect(result.domain).toBe('POLYGON_INTERIOR');
      expect(result.profileReconstruction).toBe('CHORDAL_POLYGON_FROM_ORDERED_SAMPLES');
      expect(result.radius).toBeLessThanOrEqual(profile.cgal.mic.radius + ORACLE_ERROR_TOLERANCE);
      expect(profile.cgal.mic.radius - result.radius).toBeLessThan(SIMULATED_ANNEALING_TOLERANCE);
    }
  );

  it.each(profiles.map(profile => [profile.name, profile]))(
    'keeps simulated annealing close to the CGAL MCC oracle: %s',
    (_name, profile) => {
      const result = simulatedAnnealing(
        profile.points,
        'MCC',
        meanPoint(profile.points),
        SIMULATED_ANNEALING_STEP,
        {
          ...SIMULATED_ANNEALING_OPTIONS,
          random: seededRandom(seedFromString(`${profile.name}:MCC`))
        }
      );

      expect(result.domain).toBe('POINT_SET_EQUIVALENT_TO_POLYGON');
      expect(result.radius).toBeGreaterThanOrEqual(profile.cgal.mcc.radius - ORACLE_ERROR_TOLERANCE);
      expect(result.radius - profile.cgal.mcc.radius).toBeLessThan(SIMULATED_ANNEALING_TOLERANCE);
    }
  );

  it.each(profiles.map(profile => [profile.name, profile]))(
    'keeps simulated annealing close to the CGAL MZC oracle: %s',
    (_name, profile) => {
      const result = simulatedAnnealing(
        profile.points,
        'MZC',
        meanPoint(profile.points),
        SIMULATED_ANNEALING_STEP,
        {
          ...SIMULATED_ANNEALING_OPTIONS,
          random: seededRandom(seedFromString(`${profile.name}:MZC`))
        }
      );
      const width = result.outerCircle.radius - result.innerCircle.radius;

      expect(result.domain).toBe('POINT_SET');
      expect(result.objective).toBe('RADIAL_WIDTH');
      expect(width).toBeGreaterThanOrEqual(profile.cgal.mzc.width - ORACLE_ERROR_TOLERANCE);
      expect(width - profile.cgal.mzc.width).toBeLessThan(SIMULATED_ANNEALING_TOLERANCE);
    }
  );

  it('preserves MIC/MCC/MZC under the generated scale/rotate/translate transform', () => {
    const base = profiles.find(profile => profile.name === 'mixed_harmonics_64');
    const transformed = profiles.find(profile => profile.name === 'mixed_harmonics_translated_scaled');
    const expectedMicCenter = transformedPoint(base.cgal.mic.center, transformed.source);
    const expectedMccCenter = transformedPoint(base.cgal.mcc.center, transformed.source);
    const expectedMzcCenter = transformedPoint(base.cgal.mzc.center, transformed.source);

    expect(transformed.cgal.mic.radius).toBeCloseTo(
      base.cgal.mic.radius * transformed.source.scale,
      10
    );
    expect(transformed.cgal.mic.center[0]).toBeCloseTo(expectedMicCenter[0], 10);
    expect(transformed.cgal.mic.center[1]).toBeCloseTo(expectedMicCenter[1], 10);
    expect(transformed.cgal.mcc.radius).toBeCloseTo(
      base.cgal.mcc.radius * transformed.source.scale,
      10
    );
    expect(transformed.cgal.mcc.center[0]).toBeCloseTo(expectedMccCenter[0], 10);
    expect(transformed.cgal.mcc.center[1]).toBeCloseTo(expectedMccCenter[1], 10);
    expect(transformed.cgal.mzc.width).toBeCloseTo(
      base.cgal.mzc.width * transformed.source.scale,
      10
    );
    expect(transformed.cgal.mzc.center[0]).toBeCloseTo(expectedMzcCenter[0], 10);
    expect(transformed.cgal.mzc.center[1]).toBeCloseTo(expectedMzcCenter[1], 10);
  });
});
