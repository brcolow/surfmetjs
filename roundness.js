import { FFT } from './fft.js';
import { transform } from './fft2.js';
import { levenMarqFull } from './circlefit.js';

const radius = 1.5; // Target radius of the circle (assuming perfect circle) in micrometers
const tolerance = 0.1; // Adjust tolerance for desired randomness
const minRadius = radius - tolerance;
const maxRadius = radius + tolerance;
const numPoints = 512;

// const roundnessData = generateRandomRoundnessProfile(numPoints, minRadius, maxRadius);

const url = "http://localhost:5000/nist/cir2d30.ds";
readPointsFromURL(url)
  .then(points => {
    if (points) {
      // Strip out the constant z-value of the 3d points from the NIST file.
      points = points.map(point => [point[0], point[1]]);
      const leastSquaresCircle = levenMarqFull(points);
      console.log(leastSquaresCircle);

      // Note: The NIST profiles are not full circles so the last harmonic is huge compared to all the others (because of the discontinuity).
      // We could explore trying to use a window function for partial profiles such as Hanning, Hamming, and Blackman-Harris.
      analyzeHarmonics2(cartesianToPolar(points).map(point => point[0]));
    }
  })
  .catch(error => {
    console.error("Error reading data:", error);
  });

function cartesianToPolar(cartesianCoords) {
  return cartesianCoords.map(([x, y]) => {
    const r = Math.sqrt(x * x + y * y);
    const theta = Math.atan2(y, x);
    return [r, theta];
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
    const amplitude = Math.sqrt(real * real + imag * imag) * ((i === 0 ? 1 : 2) / numPoints); // This multiplies by the scaling factor so that the amplitudes in-tact.
    const phase = Math.atan2(imag, real);

    if (i === 0) {
      console.log("Average radius: " + amplitude);
    } else {
      console.log(`Harmonic ${i / 2}: Amplitude = ${amplitude}, Phase = ${phase}`);
    }
  }
}

function analyzeHarmonics2(roundnessData) {
  const out = roundnessData.slice();
  transform(out, new Array(roundnessData.length).fill(0));

  for (let i = 0; i < out.length; i += 2) {
    const real = out[i];
    const imag = out[i + 1];
    const amplitude = Math.sqrt(real * real + imag * imag) * ((i === 0 ? 1 : 2) / numPoints); // This multiplies by the scaling factor so that the amplitudes in-tact.
    const phase = Math.atan2(imag, real);

    if (i === 0) {
      console.log("Average radius: " + amplitude);
    } else {
      if (amplitude > 1) {
        console.log(`Harmonic ${i / 2}: Amplitude = ${amplitude}, Phase = ${phase}`);
      }
    }
  }
}

function generateRandomRoundnessProfile(numPoints, minRadius, maxRadius) {
  const data = [];
  for (let i = 0; i < numPoints; i++) {
    // Generate random value between min and max radius
    const randomValue = Math.random() * (maxRadius - minRadius) + minRadius;
    data.push(randomValue);
  }
  return data;
}

// Need to convert our equally-spaced N data-points around the circle of radius measurements into (θ, r) and then (x, y) coordinates i.e. polar => cartesian (or can you use polar coordinates with the Levenberg–Marquardt algorithm?) which would be done using:
// x = r cos θ , y = r sin θ
const dataPoints = [
  { x: 1, y: 2 },
  { x: 3, y: 4 },
  { x: 5, y: 3 },
];

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
