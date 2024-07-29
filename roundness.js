import { FFT } from './fft.js';
import { transform } from './fft2.js';

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
      points = points.map(point => [ point[0], point[1] ]);
      const leastSquaresCircle = levenMarqFull(points);
      console.log(leastSquaresCircle);
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

function twoDArray(rows, cols) {
  let x = new Array(rows);

  for (let i = 0; i < rows; i++) {
    x[i] = new Array(cols).fill(0);
  }
  return x;
}

/**
 * Kåsa circle fit.
 * 
 * Adapted from https://github.com/mdoube/BoneJ/blob/1c5a90873f5d2dda0294dabc9698f002670aa816/src/org/doube/geometry/FitCircle.java#L53
 * 
 * @param {*} points 2d array containing x,y points (double[n][2])
 * @returns circle object centered at (a, b) with radius r
 */
function kasaFit(points) {
  const nPoints = points.length;
  if (nPoints < 3) {
    throw new Error("Too few points: must have at least 3");
  }
  const z = twoDArray(nPoints, 1);
  const xy1 = twoDArray(nPoints, 3);
  for (let n = 0; n < nPoints; n++) {
    z[n][0] = points[n][0] * points[n][0] + points[n][1] * points[n][1];
    xy1[n][0] = points[n][0];
    xy1[n][1] = points[n][1];
    xy1[n][2] = 1;
  }

  const XY1 = math.matrix(xy1);
  const Z = math.matrix(z);
  const P = XY1.size()[0] === XY1.size()[1] ? math.multiply(math.inv(XY1), Z) : math.multiply(math.pinv(XY1), Z);

  const p0 = P.get([0, 0]);
  const p1 = P.get([1, 0]);
  const p2 = P.get([2, 0]);
  return {"a": p0 / 2, "b": p1 / 2, "r": Math.sqrt((p0 * p0 + p1 * p1) / 4 + p2)};
}

/**
 * Levenberg-Marquardt fit in the "full" (a, b, R) space.
 * 
 * Adapted from https://github.com/mdoube/BoneJ/blob/1c5a90873f5d2dda0294dabc9698f002670aa816/src/org/doube/geometry/FitCircle.java#L500
 *
 * @param points 2d array containing x,y points (double[n][2])
 * @param lambdaIni the initial value of lambda which controls the dampening - larger values means it's more like gradient descent. Defaults to 1.
 * @return circle object centered at (a, b) with radius r
 */
function levenMarqFull(points, lambdaIni = 1) {
  const nPoints = points.length;
  if (nPoints < 3) {
    throw new Error("Too few points: must have at least 3");
  }
  const guess = kasaFit(points);
  let x = guess.a
  let y = guess.b;
  let r = guess.r;
  const par = [[ x, y, r ]];
  let Par = math.matrix(par);
  let ParTemp = math.matrix(par);
  const epsilon = 1e-6;
  let progress = epsilon;
  const iterMax = 50;
  let lambda_sqrt = Math.sqrt(lambdaIni);

  let f = 0;
  let j = twoDArray(nPoints + 3, 3);
  let g = twoDArray(nPoints + 3, 1);
  for (let i = 0; i < nPoints; i++) {
    const dX = points[i][0] - x;
    const dY = points[i][1] - y;
    const d = Math.sqrt(dX * dX + dY * dY);
    j[i][0] = -dX / d;
    j[i][1] = -dY / d;
    j[i][2] = -1;
    g[i][0] = d - r;
    f += (d - r) * (d - r);
  }

  const J = math.matrix(j);
  const G = math.matrix(g);

  let fTemp = 0;
  const jTemp = twoDArray(nPoints + 3, 3);
  const gTemp = twoDArray(nPoints + 3, 1);

  for (let iter = 0; iter < iterMax; iter++) {
    let safety = 0;
    while (safety < 100) {
      safety++;
      J.set([nPoints, 0], lambda_sqrt);
      J.set([nPoints + 1, 1], lambda_sqrt);
      J.set([nPoints + 2, 2], lambda_sqrt);
      G.set([nPoints, 0], 0);
      G.set([nPoints + 1, 0], 0);
      G.set([nPoints + 2, 0], 0);

      const DelPar = J.size()[0] === J.size()[1] ? math.multiply(math.inv(J), G) : math.multiply(math.pinv(J), G);
      progress = math.norm(DelPar, 'fro') / (math.norm(Par, 'fro') + epsilon);
      if (progress < epsilon) {
        break;
      }
      ParTemp = math.subtract(Par, math.transpose(DelPar));
      x = ParTemp.get([0, 0]);
      y = ParTemp.get([0, 1]);
      r = ParTemp.get([0, 2]);
      for (let i = 0; i < nPoints; i++) {
        const dX = points[i][0] - x;
        const dY = points[i][1] - y;
        const d = Math.sqrt(dX * dX + dY * dY);
        jTemp[i][0] = -dX / d;
        jTemp[i][1] = -dY / d;
        jTemp[i][2] = -1;
        gTemp[i][0] = d - r;
        fTemp += (d - r) * (d - r);
      }
      if (fTemp < f && ParTemp.get([0, 2]) > 0) {
        lambda_sqrt /= 2;
        break;
      } 
      lambda_sqrt *= 2;
      continue;
    }
    if (progress < epsilon) {
      break;
    }
    Par = ParTemp;
    j = jTemp;
    g = gTemp;
    f = fTemp;
  }
  return {"a": Par.get([0, 0]), "b": Par.get([0, 1]), "r": Par.get([0, 2])};
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
