import { FFT } from './fft.js';
import { transform } from './fft2.js';
import { levenMarqFull } from './circlefit.js';
import { simulatedAnnealing } from './simulated_annealing.js';
import { gradientDescent } from './gradient_descent.js';
import * as Plotly from 'plotly.js-dist-min';
import Chart from 'chart.js/auto';
import nistData from './nist/cir2d22.ds';
import { getCurvePoints, getCurvePoints2 } from './cardinal_spline.js';

// const roundnessData = generateRandomRoundnessProfile(numPoints, minRadius, maxRadius);

function createPolarChart(points) {
  // We want to also be able to draw the LSC, MIC, etc but plotly only accepts data as [r, theta]. We could either create a bunch of fake points on those circles and draw in lines only mode (no markers)
  // or extend plotly ourselves or make our own polar chart with the canvas (most likely).
  const rawPoints = {
    r: cartesianToPolar(points).map(point => point[0]),
    theta: cartesianToPolar(points).map(point => point[1] * (180 / Math.PI)),
    mode: 'lines',
    name: 'Traced Profile',
    line: {
      color: 'green',
      size: 2
    },
    marker: {
      color: 'green',
      size: 2
    },
    mode: 'lines+markers',
    type: 'scatterpolar'
  };
  const layout = {
    font: {
      family: 'Arial, sans-serif;',
      size: 12,
    },
    showlegend: true,
    orientation: -90
  };

  Plotly.newPlot(document.getElementById("polarChart"), [rawPoints], layout, { displaylogo: false, scrollZoom: true });
}

function cartesianToPolar(cartesianCoords) {
  return cartesianCoords.map(([x, y]) => {
    const r = Math.sqrt(x * x + y * y);
    const theta = Math.atan2(y, x);
    return [r, theta];
  });
}

function polarToCartesian(polarCoords) {
  return polarCoords.map(([radius, angle]) => {
    const x = radius * Math.cos(angle);
    const y = radius * Math.sin(angle);
    return [x, y];
  });
}

function analyzeHarmonics(roundnessData) {
  const f = new FFT(roundnessData.length);
  const out = f.createComplexArray();

  f.realTransform(out, roundnessData);

  for (let i = 0; i < out.length; i += 2) {
    const real = out[i]; // a_n (cos)
    const imag = out[i + 1]; // b_n (sin)
    // amplitude = c_n = sqrt(a_n^2 + b_n^2)
    // phase = γ_n = arctan(b_n / a_n)
    // a_n = c_n cos(γ_n) 
    // b_n = c_n sin(γ_n)
    const amplitude = Math.sqrt(real * real + imag * imag) * ((i === 0 ? 1 : 2) / roundnessData.length); // This multiplies by the scaling factor so that the amplitudes stay in-tact.
    const phase = Math.atan2(imag, real);

    if (i === 0) {
      console.log("Average radius (DC component): " + amplitude);
    } else {
      console.log(`Harmonic ${i / 2}: Amplitude = ${amplitude}, Phase = ${phase}`);
    }
  }
}

/**
 * Calculates the Discrete Fourier Transform sample frequencies.
 *
 * @param {number} n - The length of the input signal (number of samples).
 * @param {number} d - The sample spacing (inverse of the sampling rate). Defaults to 1.
 * @returns {Array<number>} An array containing the frequency bin centers in cycles per unit of the sample spacing.
 */
function fftfreq(n, d = 1) {
  const N = n;
  const vals = new Array(N);

  if (N % 2 === 0) {
    // Even number of samples
    for (let i = 0; i < N / 2; i++) {
      vals[i] = i / (d * N);
    }
    for (let i = N / 2; i < N; i++) {
      vals[i] = (i - N) / (d * N);
    }
  } else {
    // Odd number of samples
    for (let i = 0; i < N; i++) {
      vals[i] = (i - (N - 1) / 2) / (d * N);
    }
  }

  return vals;
}

function analyzeHarmonics2(roundnessData) {
  const amplitudes = [];
  const out = roundnessData.slice();
  transform(out, new Array(roundnessData.length).fill(0));

  for (let i = 0; i < out.length; i += 2) {
    const real = out[i];
    const imag = out[i + 1];
    const amplitude = Math.sqrt(real * real + imag * imag) * (1 / roundnessData.length); // This multiplies by the scaling factor so that the amplitudes stay in-tact.
    const phase = Math.atan2(imag, real);


    if (i === 0) {
      console.log("Unscaled amplitude: " + Math.sqrt(real * real + imag * imag));
      console.log("Average radius (DC component): " + amplitude);
      amplitudes.push({ x: String(i / 2), y: amplitude });
    } else {
      if (i < 100) {
        console.log(`Harmonic ${i / 2}: Amplitude = ${amplitude}, Phase = ${phase}`);
        amplitudes.push({ x: String(i / 2), y: amplitude });

      }
    }
  }

  return amplitudes;
}

function verifyMic(points, center, radius) {
  for (let point of points) {
    const distance = distanceSquared(point, center);
    if (Math.abs(distance - radius * radius) < 0.0000000001) {
      console.log("POINT " + point + " IS ON CIRCLE!");
    }
    if (distance < radius * radius) {
      throw new Error("THIS IS NOT A MIC!");
    }
  }
}

function generateRandomRoundnessProfile(numPoints, minRadius, maxRadius) {
  const data = [];
  for (let i = 0; i < numPoints; i++) {
    // Generate random value between min and max radius
    const randomRadius = Math.random() * (maxRadius - minRadius) + minRadius;
    const theta = 360 * (i / numPoints);
    data.push([randomRadius, theta]);
  }
  return data;
}

async function readPointsFromURL(url, pointsPerLine = 2) {
  try {
    const response = await fetch(url);
    if (!response.ok) {
      throw new Error(`HTTP error! status: ${response.status}`);
    }
    const data = await response.text();
    const points = [];
    for (const line of data.split('\n')) {
      if (pointsPerLine === 2) {
        const [x, y] = line.trim().split(/\s+/);
        if (x && y) { // Check for empty lines
          points.push([parseFloat(x), parseFloat(y)]);
        }
      } else if (pointsPerLine === 3) {
        const [x, y, z] = line.trim().split(/\s+/);
        if (x && y && z) { // Check for empty lines
          points.push([parseFloat(x), parseFloat(y), parseFloat(z)]);
        }
      }
    }
    return points;
  } catch (error) {
    console.error("Error fetching data:", error);
    return null;
  }
}

/**
 * Calculates the (squared) distance between two points a and b.
 *
 * @param {*} a First point.
 * @param {*} b Second point.
 * @returns The squared distance between points a and b.
 */
function distanceSquared(a, b) {
  return (a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]);
}

function calculateCircle(p1, p2, p3) {
  // Reference: http://www.faqs.org/faqs/graphics/algorithms-faq/
  // Subject 1.04: How do I generate a circle through three points?
  const A = p2[0] - p1[0];
  const B = p2[1] - p1[1];
  const C = p3[0] - p1[0];
  const D = p3[1] - p1[1];

  const E = A * (p1[0] + p2[0]) + B * (p1[1] + p2[1]);
  const F = C * (p1[0] + p3[0]) + D * (p1[1] + p3[1]);

  const G = 2.0 * (A * (p3[1] - p2[1]) - B * (p3[0]- p2[0]));
  if (G === 0) {
    throw new Error("The three points given to calculateCircumcenter must not be co-linear.");
  }

  const cx = (D * E - B * F) / G
  const cy = (A * F - C * E) / G

  const dx = cx - p1[0];
  const dy = cy - p1[1];
  const radiusSquared = dx * dx + dy * dy;

  return { center: [cx, cy], radiusSquared: radiusSquared};
}

function bruteForceMic(points) {
  let biggestMic = null;

  for (let i = 0; i < points.length - 2; i++) {
    for (let j = i + 1; j < points.length - 1; j++) {
      for (let k = j + 1; k < points.length; k++) {
        let potentialMic = calculateCircle(points[i], points[j], points[k]);
        let isMicCandidate = true;
        for (let point of points) {
          const distance = distanceSquared(point, potentialMic.center);
          if (distance < potentialMic.radiusSquared) {
            isMicCandidate = false;
            break;
          }
        }

        if (isMicCandidate) {
          if (biggestMic == null || potentialMic.radiusSquared > biggestMic.radiusSquared) {
            biggestMic = potentialMic;
          }
        }
      }
    }
  }

  return { center: [biggestMic.center[0], biggestMic.center[1]], radius: Math.sqrt(biggestMic.radiusSquared) };
}

/**
 * Uses the Box–Muller transform to return a random number between (min, max) where the probability of picking 
 * a number from that range follows the normal distribution. From https://stackoverflow.com/a/36481059.
 * 
 * See: https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
 */
function randomNormalDistribution(min, max) {
  let u, v;

  // Generate random numbers until both are non-zero
  do {
    u = Math.random();
    v = Math.random();
  } while (u === 0 || v === 0);

  let num = Math.sqrt(-2.0 * Math.log(u)) * Math.cos(2.0 * Math.PI * v);
  
  // Scale and shift to the desired range
  num = (num + 1) / 2 * (max - min) + min;

  // Resample if outside the range
  if (num < min || num > max) {
    return randomNormalDistribution(min, max);
  }

  return num;
}

// const url = "http://localhost:5000/nist/cir2d22.ds";
readPointsFromURL(nistData)
  .then(points => {
    if (points) {
      console.log(randomNormalDistribution(0, 50));
      const curvePoints = getCurvePoints(points.flatMap(pair => pair), 0.5, 25, true);
      const curvePoints2 = getCurvePoints2(points.flatMap(pair => pair), 0.5, 25, true);

      console.log(curvePoints);
      console.log(curvePoints2);
      const leastSquaresCircle = levenMarqFull(points);
      let timeStart = performance.now();
      const bruceForceMic = bruteForceMic(points);
      console.log("Brute force time: " + (performance.now() - timeStart) + " milliseconds");
      console.log("Bruce force MIC:");
      console.log(bruceForceMic);
      console.log(points);

      console.log(gradientDescent(points, "MIC", leastSquaresCircle));

      // Subtract the LSC center from the points to center the point-cloud around the origin (0,0).
      // const centeredPoints = points.slice().map(point => [point[0] - leastSquaresCircle.a, point[1] - leastSquaresCircle.b]);

      // We should use the order of max radius - min radius to determine the stepSize for neighbors.
      const maxRadius = Math.sqrt(Math.max(...points.map(p => distanceSquared(p, [leastSquaresCircle.a, leastSquaresCircle.b]))));
      const minRadius = Math.sqrt(Math.min(...points.map(p => distanceSquared(p, [leastSquaresCircle.a, leastSquaresCircle.b]))));
      console.log("maxRadius - minRadius: " + (maxRadius - minRadius));

      console.log("NIST MIC: center = [-600.5093622580, -428.7134351928], radius = 169.4623601410");
      verifyMic(points, [-600.5093622669094, -428.713435109361], 169.46236020895796);

      timeStart = performance.now();
      let biggestMic = null;
      for (let i = 0; i < 100; i++) {
        const mic = simulatedAnnealing(points.slice(), "MIC", leastSquaresCircle, 1000, { type: 'logarithmic', rate: 0.1 }, 1000, 5000, 0.001);
        if (biggestMic == null || mic.radius > biggestMic.radius) {
          biggestMic = mic;
        }
      }
  
      console.log("SA time: " + (performance.now() - timeStart) + " milliseconds");
      console.log("Simulated Annealing MIC:");
      console.log({ center: [biggestMic.center[0], biggestMic.center[1]], radius: biggestMic.radius});
      if (true) {
        return;
      }

      let data3 = polarToCartesian(generateRandomRoundnessProfile(3, 5, 6));
      console.log("Random data of 3 points around circle: " + data3);
      let MIC = maximumIscribedCircle(data3);
      console.log(MIC);

      let data4 = polarToCartesian(generateRandomRoundnessProfile(4, 5, 6));
      console.log("Random data of 4 points around circle: " + data4);
      MIC = maximumIscribedCircle(data4);
      console.log(MIC);
      createPolarChart(data4);

      if (true) {
        return;
      }

      // Strip out the constant z-value of the 3d points from the NIST file.
      points = points.map(point => [point[0], point[1]]);
      leastSquaresCircle = levenMarqFull(points);
      console.log(leastSquaresCircle);
      // Subtract the LSC center from the points to center the point-cloud around the origin (0,0).
      centeredPoints = points.map(point => [point[0] - leastSquaresCircle.a, point[1] - leastSquaresCircle.b]);
      // Note: The NIST profiles are not full circles so the last harmonic is huge compared to all the others (because of the discontinuity).
      // We could explore trying to use a window function for partial profiles such as Hanning, Hamming, and Blackman-Harris.
      const amplitudes = analyzeHarmonics2(cartesianToPolar(centeredPoints).map(point => point[0]));
      console.log(maximumIscribedCircle(points));
      createPolarChart(centeredPoints);

      const ctx = document.getElementById('harmonicChart');
      new Chart(ctx, {
        type: 'bar',
        data: {
          datasets: [{
            data: amplitudes,
          }]
        },
        options: {
          scales: {
            x: {
              display: true,
            },
            y: {
              display: true,
              type: 'logarithmic',
              beginAtZero: true,
              ticks: {
                callback: function (val, index) {
                  // If tick value is a whole number, don't show decimal places.
                  return Number.isInteger(val) ? Math.trunc(val) : val;
                },
              }
            }
          },
          plugins: {
            title: {
              display: true,
              text: 'Harmonic Amplitudes',
            },
            legend: {
              display: false
            }
          }
        }
      });
    }
  })
  .catch(error => {
    console.error("Error reading data:", error);
  });