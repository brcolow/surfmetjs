/*! Curve calc function for canvas 3.0.1
 *  (c) Epistemex 2013-2018, 2024
 *  License: MIT
 */

'use strict';

/**
 * Calculates a flattened point array representing a cardinal spline through
 * the given flattened point array: [x1, y1, x2, y2, ..., xn, yn].
 *
 * The returned points include each segment start and the final endpoint. For
 * closed splines, the final endpoint is the first input point.
 *
 * @param {Array<number>|TypedArray} points - Flattened point array.
 * @param {number} [tension=0.5] - Tangent scale.
 * @param {number} [numOfSeg=25] - Number of samples per segment.
 * @param {boolean} [closed=false] - Whether to close the spline.
 * @returns {Float32Array} Flattened spline points.
 */
function getCurvePoints(points, tension = 0.5, numOfSeg = 25, closed = false) {
  if (!points || points.length < 2) {
    return new Float32Array(0);
  }

  validateSplineInput(points, numOfSeg);

  const pointCount = points.length / 2;
  if (pointCount < 2) {
    return new Float32Array(0);
  }

  const segmentCount = closed ? pointCount : pointCount - 1;
  const result = new Float32Array((segmentCount * numOfSeg + 1) * 2);
  let resultIndex = 0;

  for (let segment = 0; segment < segmentCount; segment++) {
    const p0 = getControlPoint(points, segment - 1, pointCount, closed);
    const p1 = getControlPoint(points, segment, pointCount, closed);
    const p2 = getControlPoint(points, segment + 1, pointCount, closed);
    const p3 = getControlPoint(points, segment + 2, pointCount, closed);
    const tangent1 = [
      (p2[0] - p0[0]) * tension,
      (p2[1] - p0[1]) * tension
    ];
    const tangent2 = [
      (p3[0] - p1[0]) * tension,
      (p3[1] - p1[1]) * tension
    ];

    for (let step = 0; step < numOfSeg; step++) {
      const t = step / numOfSeg;
      const t2 = t * t;
      const t3 = t2 * t;
      const h00 = 2 * t3 - 3 * t2 + 1;
      const h10 = t3 - 2 * t2 + t;
      const h01 = -2 * t3 + 3 * t2;
      const h11 = t3 - t2;

      result[resultIndex++] = h00 * p1[0] + h10 * tangent1[0] + h01 * p2[0] + h11 * tangent2[0];
      result[resultIndex++] = h00 * p1[1] + h10 * tangent1[1] + h01 * p2[1] + h11 * tangent2[1];
    }
  }

  const finalPoint = closed
    ? getControlPoint(points, 0, pointCount, closed)
    : getControlPoint(points, pointCount - 1, pointCount, closed);
  result[resultIndex++] = finalPoint[0];
  result[resultIndex] = finalPoint[1];

  return result;
}

function validateSplineInput(points, numOfSeg) {
  if (points.length % 2 !== 0) {
    throw new Error('Point array must contain x,y pairs.');
  }

  if (!Number.isInteger(numOfSeg) || numOfSeg < 1) {
    throw new Error('Number of segments cannot be less than one.');
  }
}

function getControlPoint(points, index, pointCount, closed) {
  const pointIndex = closed
    ? positiveModulo(index, pointCount)
    : Math.max(0, Math.min(pointCount - 1, index));
  const offset = pointIndex * 2;

  return [points[offset], points[offset + 1]];
}

function positiveModulo(value, divisor) {
  return ((value % divisor) + divisor) % divisor;
}

function getCurvePoints2(points, tension = 0.5, numOfSeg = 25, closed = false) {
  return getCurvePoints(points, tension, numOfSeg, closed);
}

export { getCurvePoints, getCurvePoints2 };
