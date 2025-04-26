import { describe, it, expect, beforeEach } from 'vitest';
import { adaptiveGradientDescent, objectiveFunctionMZCSmooth, gradient } from '../gradient_descent';

describe('gradient descent', () => {

  const P = [ 
    [1, 0], [0, 1], [-1, 0], [0, -1], 
    [2, 0], [0, 2], [-2, 0], [0, -2], 
  ];
  const gamma = 200;

  describe('functions for MZC', () => {
    it('should compute objective function for MZC correctly', () => {
      expect(objectiveFunctionMZCSmooth([0, 0, 0], P, gamma)).toEqual(1.0);
      expect(objectiveFunctionMZCSmooth([0.5, 0, 0], P, gamma)).toEqual(2.0);
      // Result value taken from python
      expect(objectiveFunctionMZCSmooth([.5, .5, 0], P, gamma)).toBeCloseTo(1.8424029756098448, 20);
      // Result value taken from python 
      expect(objectiveFunctionMZCSmooth([.5, .75, 0], P, gamma)).toBeCloseTo(2.23606797749979, 20);
    });

    it('should compute gradient function for MZC correctly', () => {
      let dx, dy, _;
      // We don't work with dr, hence just _
      [dx, dy, _] = gradient([0, 0, 0], P, 'MZC');
      expect(dx).toBeCloseTo(0, 20);
      expect(dy).toBeCloseTo(0, 20);
      
      // Again, results taken from python
      [dx, dy, _] = gradient([0.5, 0, 0], P, 'MZC');
      expect(dx).toBeCloseTo(2., 20);
      expect(dy).toBeCloseTo(0, 20);
      
      // Again, results taken from python - this time be less accurate
      [dx, dy, _] = gradient([0.5, 0.75, 0], P, 'MZC');
      expect(dx).toBeCloseTo(-0.71554175, 7);
      expect(dy).toBeCloseTo(1.43108351, 7);
    });
  });

  describe('adaptive gradient descent', () => {
    it('should descent towards [0, 0]', () => {
      let adGD;

      adGD = adaptiveGradientDescent(P, 'MZC', [0.5, 0.25, 0], gamma);
      expect(adGD.center[0]).toBeCloseTo(0, 6);
      expect(adGD.center[1]).toBeCloseTo(0, 6);
      
      adGD = adaptiveGradientDescent(P, 'MZC', [10, 10, 0], gamma);
      expect(adGD.center[0]).toBeCloseTo(0, 6);
      expect(adGD.center[1]).toBeCloseTo(0, 6);
    });
  });

});

describe('adaptive gradient descent w/ more random like numbers', () => {

  // Values generated with python
  const P = [
    [0.63696169, 0.26978671],
    [0.04097352, 0.01652764],
    [0.81327024, 0.91275558],
    [0.60663578, 0.72949656],
    [0.54362499, 0.93507242],
  ];
  const gamma = 200;
  const minimizer = [0.08995785, 0.68466782];

  describe('adaptive gradient descent', () => {
    it('should descent towards [0, 0]', () => {
      let adGD;

      adGD = adaptiveGradientDescent(
        P, 
        'MZC', 
        [0.5, 0.25, 0], 
        gamma,
        0.5,
        0.25
      );
      expect(adGD.center[0]).toBeCloseTo(minimizer[0], 4);
      expect(adGD.center[1]).toBeCloseTo(minimizer[1], 4);
    });
  });

});