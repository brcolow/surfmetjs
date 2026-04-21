import { describe, expect, it } from 'vitest';
import { getCurvePoints, getCurvePoints2 } from '../cardinal_spline.js';
import {
  computeConvexHull,
  findFarthestPoint,
  getInitialSolution,
  isPointInPolygon
} from '../utils.js';

function asPointKeys(points) {
  return points.map(point => point.join(',')).sort();
}

describe('cardinal spline helpers', () => {
  it('returns identical output from both exported spline helpers', () => {
    const points = [0, 0, 1, 2, 3, 1, 4, 0];

    expect(Array.from(getCurvePoints2(points, 0.5, 8, false)))
      .toEqual(Array.from(getCurvePoints(points, 0.5, 8, false)));
    expect(Array.from(getCurvePoints2(points, 0.5, 8, true)))
      .toEqual(Array.from(getCurvePoints(points, 0.5, 8, true)));
  });

  it('preserves expected endpoints and output length', () => {
    const points = [0, 0, 1, 2, 3, 1];
    const open = getCurvePoints(points, 0.5, 4, false);
    const closed = getCurvePoints(points, 0.5, 4, true);

    expect(open.length).toBe(((points.length / 2 - 1) * 4 + 1) * 2);
    expect(open[0]).toBeCloseTo(0);
    expect(open[1]).toBeCloseTo(0);
    expect(open[open.length - 2]).toBeCloseTo(3);
    expect(open[open.length - 1]).toBeCloseTo(1);

    expect(closed.length).toBe(((points.length / 2) * 4 + 1) * 2);
    expect(closed[closed.length - 2]).toBeCloseTo(0);
    expect(closed[closed.length - 1]).toBeCloseTo(0);
  });

  it('rejects malformed point arrays', () => {
    expect(() => getCurvePoints([0, 0, 1], 0.5, 4)).toThrow(/x,y pairs/);
    expect(() => getCurvePoints([0, 0, 1, 1], 0.5, 0)).toThrow(/segments/);
  });
});

describe('geometry utilities', () => {
  it('computes a convex hull without mutating the input points', () => {
    const points = [[1, 1], [0, 0], [1, 0], [0, 1], [0.5, 0.5], [1, 1]];
    const original = points.map(point => [...point]);
    const hull = computeConvexHull(points);

    expect(points).toEqual(original);
    expect(asPointKeys(hull)).toEqual(['0,0', '0,1', '1,0', '1,1']);
  });

  it('returns endpoints for a collinear hull', () => {
    const hull = computeConvexHull([[1, 0], [0, 0], [2, 0], [0.5, 0]]);

    expect(hull).toEqual([[0, 0], [2, 0]]);
  });

  it('finds the actual farthest point', () => {
    const result = findFarthestPoint([[1, 0], [2, 0], [-1, 0]], [0, 0]);

    expect(result.point).toEqual([2, 0]);
    expect(result.distance).toBeCloseTo(2);
  });

  it('treats polygon boundary points as inside by default', () => {
    const square = [[0, 0], [1, 0], [1, 1], [0, 1]];

    expect(isPointInPolygon([0, 0.5], square)).toBe(true);
    expect(isPointInPolygon([0.5, 0.5], square)).toBe(true);
    expect(isPointInPolygon([2, 0.5], square)).toBe(false);
  });

  it('produces finite initial circle estimates', () => {
    const points = Array.from({ length: 16 }, (_, i) => [
      3 + 5 * Math.cos(2 * Math.PI * i / 16),
      -2 + 5 * Math.sin(2 * Math.PI * i / 16)
    ]);

    const mcc = getInitialSolution(points, 'MCC', { a: 3, b: -2, r: 5 });
    const mzc = getInitialSolution(points, 'MZC', { a: 3, b: -2, r: 5 });

    expect(mcc.radius).toBeCloseTo(5);
    expect(mzc.innerCircle.radius).toBeCloseTo(5);
    expect(mzc.outerCircle.radius).toBeCloseTo(5);
    expect(Number.isFinite(mzc.outerCircle.radius)).toBe(true);
  });
});
