import { add, identity, matrix, multiply, pinv, transpose } from 'mathjs';

const MIN_CIRCLE_POINTS = 3;
const MIN_RADIUS = 1e-12;
const DISTANCE_EPSILON = 1e-12;

/*
 * FitCircle Java class for fitting circles to 2D coordinate data
 *
 * Copyright 2009 2010 Michael Doube
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 */

/**
 * Kasa circle fit.
 *
 * Adapted from Java https://github.com/mdoube/BoneJ/blob/1c5a90873f5d2dda0294dabc9698f002670aa816/src/org/doube/geometry/FitCircle.java#L53
 *
 * @param {Array<Array<number>>} points 2D array containing x,y points.
 * @returns {{a: number, b: number, r: number}} Circle centered at (a, b) with radius r.
 */
function kasaFit(points) {
  assertValidPoints(points);

  const z = [];
  const xy1 = [];

  for (const [x, y] of points) {
    z.push([x * x + y * y]);
    xy1.push([x, y, 1]);
  }

  const coefficients = multiply(pinv(matrix(xy1)), matrix(z));
  const p0 = coefficients.get([0, 0]);
  const p1 = coefficients.get([1, 0]);
  const p2 = coefficients.get([2, 0]);
  const radiusSquared = (p0 * p0 + p1 * p1) / 4 + p2;

  return {
    a: p0 / 2,
    b: p1 / 2,
    r: Math.sqrt(Math.max(0, radiusSquared))
  };
}

/**
 * Fits a circle to a set of 2D points using Levenberg-Marquardt geometric
 * least squares. The objective is the sum of squared radial residuals:
 * distance(point, center) - radius.
 *
 * @param {Array<Array<number>>} points - An array of [x, y] coordinate pairs.
 * @param {number} [lambdaIni=1] - Initial damping parameter.
 * @returns {{a: number, b: number, r: number}} Best-fit circle parameters.
 */
function circleFitOrthogonal(points, lambdaIni = 1) {
  assertValidPoints(points);

  const guess = kasaFit(points);
  let params = [
    guess.a,
    guess.b,
    Number.isFinite(guess.r) && guess.r > 0 ? guess.r : meanRadius(points, [guess.a, guess.b])
  ];

  let lambda = Math.max(lambdaIni, 1e-12);
  let state = evaluateCircle(points, params);
  const maxIterations = 100;
  const maxDampingSteps = 25;
  const tolerance = 1e-12;

  for (let iteration = 0; iteration < maxIterations; iteration++) {
    let accepted = false;

    for (let dampingStep = 0; dampingStep < maxDampingSteps; dampingStep++) {
      const step = solveLevenbergMarquardtStep(state.jacobian, state.residuals, lambda);
      const candidate = [
        params[0] - step[0],
        params[1] - step[1],
        params[2] - step[2]
      ];

      if (candidate[2] <= MIN_RADIUS || !candidate.every(Number.isFinite)) {
        lambda *= 10;
        continue;
      }

      const candidateState = evaluateCircle(points, candidate);
      if (candidateState.cost < state.cost) {
        const stepNorm = Math.hypot(step[0], step[1], step[2]);
        const paramsNorm = Math.hypot(params[0], params[1], params[2]);

        params = candidate;
        state = candidateState;
        lambda = Math.max(lambda / 10, 1e-12);
        accepted = true;

        if (stepNorm <= tolerance * (paramsNorm + tolerance)) {
          return circleObject(params);
        }

        break;
      }

      lambda *= 10;
    }

    if (!accepted) {
      break;
    }
  }

  return circleObject(params);
}

function solveLevenbergMarquardtStep(jacobian, residuals, lambda) {
  const jt = transpose(jacobian);
  const normalMatrix = add(multiply(jt, jacobian), multiply(lambda, identity(3)));
  const rhs = multiply(jt, residuals);
  const step = multiply(pinv(normalMatrix), rhs);

  return [
    step.get([0, 0]),
    step.get([1, 0]),
    step.get([2, 0])
  ];
}

function evaluateCircle(points, [a, b, r]) {
  const jacobian = [];
  const residuals = [];
  let cost = 0;

  for (const [px, py] of points) {
    const dx = a - px;
    const dy = b - py;
    const distance = Math.hypot(dx, dy);
    const residual = distance - r;
    const safeDistance = Math.max(distance, DISTANCE_EPSILON);

    jacobian.push([
      distance <= DISTANCE_EPSILON ? 0 : dx / safeDistance,
      distance <= DISTANCE_EPSILON ? 0 : dy / safeDistance,
      -1
    ]);
    residuals.push([residual]);
    cost += residual * residual;
  }

  return {
    cost,
    jacobian: matrix(jacobian),
    residuals: matrix(residuals)
  };
}

function meanRadius(points, center) {
  const sum = points.reduce((total, point) => {
    return total + Math.hypot(point[0] - center[0], point[1] - center[1]);
  }, 0);

  return Math.max(sum / points.length, MIN_RADIUS);
}

function circleObject([a, b, r]) {
  return { a, b, r };
}

function assertValidPoints(points) {
  if (!Array.isArray(points) || points.length < MIN_CIRCLE_POINTS) {
    throw new Error("Too few points: must have at least 3 but was: " + (points?.length ?? 0));
  }

  for (const point of points) {
    if (!Array.isArray(point) || point.length < 2 || !point.every(Number.isFinite)) {
      throw new Error("Circle fitting points must be finite [x, y] pairs.");
    }
  }
}

export { circleFitOrthogonal, kasaFit };
