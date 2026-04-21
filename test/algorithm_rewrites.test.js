import { describe, expect, it } from 'vitest';
import { FFT } from '../fft.js';
import { transform } from '../fft2.js';
import { circleFitOrthogonal, kasaFit } from '../circlefit.js';
import { simulatedAnnealing } from '../simulated_annealing.js';
import { gradientDescent } from '../gradient_descent.js';

function naiveDft(real, imag) {
  const n = real.length;
  const outReal = new Array(n).fill(0);
  const outImag = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = -2 * Math.PI * t * k / n;
      outReal[k] += real[t] * Math.cos(angle) - imag[t] * Math.sin(angle);
      outImag[k] += real[t] * Math.sin(angle) + imag[t] * Math.cos(angle);
    }
  }

  return { real: outReal, imag: outImag };
}

function exactCirclePoints(center, radius, count = 24) {
  return Array.from({ length: count }, (_, i) => [
    center[0] + radius * Math.cos(2 * Math.PI * i / count),
    center[1] + radius * Math.sin(2 * Math.PI * i / count)
  ]);
}

describe('FFT implementations', () => {
  it('fft.js round-trips real input without extra scaling', () => {
    const input = [1, 0.5, -2, 3, 0, 4, -1, 2];
    const fft = new FFT(input.length);
    const spectrum = fft.createComplexArray();
    const output = fft.createComplexArray();

    fft.realTransform(spectrum, input);
    fft.completeSpectrum(spectrum);
    fft.inverseTransform(output, spectrum);

    for (let i = 0; i < input.length; i++) {
      expect(output[2 * i]).toBeCloseTo(input[i], 12);
      expect(output[2 * i + 1]).toBeCloseTo(0, 12);
    }
  });

  it('fft2.js matches a naive DFT for non-power-of-two lengths', () => {
    const real = [1, -2, 0.5, 3, 4];
    const imag = [0.25, -1, 2, 0, 1.5];
    const expected = naiveDft(real, imag);

    transform(real, imag);

    for (let i = 0; i < real.length; i++) {
      expect(real[i]).toBeCloseTo(expected.real[i], 12);
      expect(imag[i]).toBeCloseTo(expected.imag[i], 12);
    }
  });
});

describe('circle fitting', () => {
  it('fits an exact circle with both algebraic and orthogonal fits', () => {
    const points = exactCirclePoints([3, -2], 5);

    for (const fit of [kasaFit(points), circleFitOrthogonal(points)]) {
      expect(fit.a).toBeCloseTo(3, 12);
      expect(fit.b).toBeCloseTo(-2, 12);
      expect(fit.r).toBeCloseTo(5, 12);
    }
  });
});

describe('circle optimizers', () => {
  it('derives finite MZC radii from the center instead of preserving invalid initial radii', () => {
    const points = exactCirclePoints([3, -2], 5);
    const result = simulatedAnnealing(points, 'MZC', { a: 3, b: -2, r: 5 }, 0.01, {
      maxIterations: 0
    });

    expect(result.innerCircle.radius).toBeCloseTo(5, 12);
    expect(result.outerCircle.radius).toBeCloseTo(5, 12);
    expect(Number.isFinite(result.outerCircle.radius)).toBe(true);
  });

  it('returns geometrically valid MIC and MCC radii from gradient descent', () => {
    const count = 24;
    const radius = 5;
    const points = exactCirclePoints([3, -2], radius, count);
    const chordalMicRadius = radius * Math.cos(Math.PI / count);

    const mic = gradientDescent(points, 'MIC', { a: 3, b: -2, r: radius });
    const mcc = gradientDescent(points, 'MCC', { a: 3, b: -2, r: radius });

    expect(mic.center[0]).toBeCloseTo(3, 12);
    expect(mic.center[1]).toBeCloseTo(-2, 12);
    expect(mic.radius).toBeCloseTo(chordalMicRadius, 12);
    expect(mic.domain).toBe('POLYGON_INTERIOR');
    expect(mic.profileReconstruction).toBe('CHORDAL_POLYGON_FROM_ORDERED_SAMPLES');
    expect(mcc.center[0]).toBeCloseTo(3, 12);
    expect(mcc.center[1]).toBeCloseTo(-2, 12);
    expect(mcc.radius).toBeCloseTo(radius, 12);
    expect(mcc.domain).toBe('POINT_SET_EQUIVALENT_TO_POLYGON');
  });
});
