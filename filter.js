import { FFT } from './fft.js';

/**
 * Applies a Gaussian filter to roundness profile data, compliant with ISO 16610-21.
 * The Gaussian filter has 50% transmission at the specified cutoff UPR.
 *
 * @param {number[]} profile - Array of radial values (1D roundness profile).
 * @param {'lowpass'|'highpass'|'bandpass'} passType - Type of filter to apply.
 * @param {number} cutoffLow - Lower cutoff UPR (used in highpass and bandpass).
 * @param {number} cutoffHigh - Upper cutoff UPR (used in lowpass and bandpass).
 * @returns {number[]} Filtered roundness profile.
 */
function gaussianFilter(profile, passType, cutoffLow = 0, cutoffHigh = Infinity) {
  const N = profile.length;

  // Create filter kernels based on the UPR cutoffs
  const { kernel: gHigh, radius: rH } = makeKernel(sigmaFromUPR(cutoffHigh));
  const { kernel: gLow, radius: rL } = makeKernel(sigmaFromUPR(cutoffLow));

  // Apply the appropriate filter operation
  let result;
  if (passType === 'lowpass') {
    result = convolve(profile, gHigh, rH);
  } else if (passType === 'highpass') {
    const lowpassed = convolve(profile, gLow, rL);
    result = profile.map((v, i) => v - lowpassed[i]);
  } else if (passType === 'bandpass') {
    const low = convolve(profile, gLow, rL);
    const high = convolve(profile, gHigh, rH);
    result = low.map((v, i) => v - high[i]);
  } else {
    throw new Error("Invalid passType: " + passType);
  }

  return result;
}

/**
 * Applies a robust Gaussian filter to roundness profile data.
 * Supports low-pass, high-pass, and band-pass modes.
 * Based on ISO 16610-31 with Tukey’s biweight M-estimator.
 *
 * @param {number[]} profile - Roundness profile (radial values).
 * @param {'lowpass'|'highpass'|'bandpass'} passType - Filter mode.
 * @param {number} cutoffLow - Lower cutoff UPR (used in high-pass/band-pass).
 * @param {number} cutoffHigh - Upper cutoff UPR (used in low-pass/band-pass).
 * @param {number} [maxIter=10] - Max iterations for convergence.
 * @param {number} [tune=4.685] - Tuning constant for outlier suppression.
 * @returns {number[]} - Filtered roundness profile.
 */
function robustFilter(profile, passType, cutoffLow = 0, cutoffHigh = Infinity, maxIter = 10, tune = 4.685) {
  const N = profile.length;

  function applyRobustSmoother(input, cutoffUPR) {
    const sigma = sigmaFromUPR(cutoffUPR, N);
    const { kernel, radius } = makeGaussianKernel(sigma);
    let estimate = [...input];

    for (let iter = 0; iter < maxIter; iter++) {
      const residuals = input.map((v, i) => v - estimate[i]);

      const sorted = residuals.slice().sort((a, b) => a - b);
      const median = sorted[Math.floor(N / 2)];
      const absDev = residuals.map(r => Math.abs(r - median));
      const mad = absDev.sort((a, b) => a - b)[Math.floor(N / 2)] || 1e-6;

      const weights = residuals.map(r => {
        const u = r / (tune * mad);
        return Math.abs(u) >= 1 ? 0 : (1 - u ** 2) ** 2;
      });

      const weighted = input.map((v, i) => v * weights[i]);
      const smoothWeighted = convolveCircular(weighted, kernel, radius);
      const weightSum = convolveCircular(weights, kernel, radius);

      estimate = smoothWeighted.map((v, i) =>
        weightSum[i] !== 0 ? v / weightSum[i] : input[i]
      );
    }

    return estimate;
  }

  if (passType === 'lowpass') {
    return applyRobustSmoother(profile, cutoffHigh);
  }

  if (passType === 'highpass') {
    const smooth = applyRobustSmoother(profile, cutoffLow);
    // Subtract smoothed signal from original to extract high-frequency components.
    return profile.map((v, i) => v - smooth[i]);
  }

  if (passType === 'bandpass') {
    const smoothLow = applyRobustSmoother(profile, cutoffLow);
    const smoothHigh = applyRobustSmoother(profile, cutoffHigh);
    // Subtract the high-cutoff smoothed profile (removes roughness) from the low-cutoff smoothed profile (removes form)
    // Result is the band-pass filtered profile that retains only mid-scale features (e.g., waviness).
    return smoothLow.map((v, i) => v - smoothHigh[i]);
  }

  throw new Error(`Invalid passType: ${passType}`);
}

/**
 * Applies a 2CR (Two-Cutoff-Radii) Gaussian filter to roundness profile data.
 * This filter isolates a band of spatial frequencies (UPR) by subtracting a high-pass
 * Gaussian filter from a low-pass Gaussian filter, consistent with ISO 12181-2.
 *
 * Effectively:
 *   filtered = LowPass(cutoffLow) − LowPass(cutoffHigh)
 *
 * This removes both form and fine roughness, preserving mid-scale waviness features.
 *
 * @param {number[]} profile - The 1D roundness profile data (radial values), uniformly sampled.
 * @param {number} [cutoffLow=0] - Lower cutoff frequency in Undulations Per Revolution (UPR).
 *                                 Removes low-frequency content (e.g., form); must be > 0 to apply.
 * @param {number} [cutoffHigh=Infinity] - Upper cutoff frequency in UPR.
 *                                         Removes high-frequency content (e.g., roughness).
 * @returns {number[]} - The band-pass filtered profile isolating the waviness zone between the two cutoffs.
 *
 * @example
 * // Isolate 15–150 UPR waviness band
 * const filtered = twoCRFilter(roundnessData, 15, 150);
 */
function twoCRFilter(profile, cutoffLow = 0, cutoffHigh = Infinity) {
  const N = profile.length;

  const sigmaLow = sigmaFromUPR(cutoffLow, N);   // long wavelength (removes form)
  const sigmaHigh = sigmaFromUPR(cutoffHigh, N); // short wavelength (removes roughness)

  const { kernel: kLow, radius: rLow } = makeGaussianKernel(sigmaLow);
  const { kernel: kHigh, radius: rHigh } = makeGaussianKernel(sigmaHigh);

  const smoothLow = convolveCircular(profile, kLow, rLow);
  const smoothHigh = convolveCircular(profile, kHigh, rHigh);

  return smoothLow.map((v, i) => v - smoothHigh[i]); // band-pass result
}

/**
 * Performs a Fourier filter (1D) on the input roundness profile.
 * 
 * @param {number[]} profile - Real-valued input data.
 * @param {'lowpass'|'highpass'|'bandpass'} passType - Type of filter.
 * @param {number} cutoffLow - Lower frequency in units of UPR (cycles/rev).
 * @param {number} cutoffHigh - Upper frequency in UPR.
 * @returns {number[]} - Filtered real profile.
 */
function fourierFilter(profile, passType, cutoffLow = 0, cutoffHigh = Infinity) {
  const N = profile.length;
  const fft = new FFT(N);

  const input = new Float64Array(N);
  profile.forEach((v, i) => input[i] = v);

  const spectrum = fft.createComplexArray();
  fft.realTransform(spectrum, input);
  fft.completeSpectrum(spectrum);

  const nyquist = N / 2;

  // Zero out unwanted frequencies
  for (let i = 0; i < nyquist; i++) {
    const upr = i;

    let keep = false;
    if (passType === 'lowpass') {
      keep = upr <= cutoffHigh;
    } else if (passType === 'highpass') {
      keep = upr >= cutoffLow;
    } else if (passType === 'bandpass') {
      keep = (upr >= cutoffLow && upr <= cutoffHigh);
    }

    if (!keep) {
      const reIdx = 2 * i;
      const imIdx = 2 * i + 1;
      spectrum[reIdx] = 0;
      spectrum[imIdx] = 0;

      // Zero out the mirrored frequency too
      const mirrorIdx = 2 * (N - i);
      spectrum[mirrorIdx] = 0;
      spectrum[mirrorIdx + 1] = 0;
    }
  }

  // Inverse transform
  const output = new Float64Array(N);
  fft.inverseTransform(output, spectrum);

  // Normalize result (fft.js does not normalize)
  for (let i = 0; i < N; i++) {
    output[i] /= N;
  }

  return Array.from(output);
}

/**
 * Converts Undulations Per Revolution (UPR) to Gaussian sigma in samples,
 * based on ISO 16610-21: σ = λ / (2π), where λ = N / UPR.
 */
function sigmaFromUPR(upr, N) {
  const wavelength = N / upr;
  return wavelength / (2 * Math.PI);
}

/**
 * Creates a normalized Gaussian kernel for convolution.
 * Truncated at ±3σ for efficiency and ISO-compliant shape.
 */
function makeGaussianKernel(sigma) {
  const radius = Math.ceil(3 * sigma);
  const kernel = new Float64Array(2 * radius + 1);
  let sum = 0;

  for (let i = -radius; i <= radius; i++) {
    const value = Math.exp(-(i * i) / (2 * sigma * sigma));
    kernel[i + radius] = value;
    sum += value;
  }

  for (let i = 0; i < kernel.length; i++) kernel[i] /= sum;

  return { kernel, radius };
}

/**
 * Performs circular convolution on 1D periodic data.
 */
function convolveCircular(data, kernel, radius) {
  const N = data.length;
  const out = new Array(N);

  for (let i = 0; i < N; i++) {
    let sum = 0;
    for (let k = -radius; k <= radius; k++) {
      const idx = (i + k + N) % N; // wrap around
      sum += data[idx] * kernel[k + radius];
    }
    out[i] = sum;
  }

  return out;
}

/**
 * Unified filter dispatcher for roundness profile filtering.
 * Selects filter type and pass-band mode based on input parameters.
 *
 * @param {number[]} profile - 1D roundness profile (uniform radial values).
 * @param {{
*   filterType: 'gaussian' | 'fourier' | 'robust' | '2cr',
*   passType: 'lowpass' | 'highpass' | 'bandpass',
*   cutoffLow?: number,
*   cutoffHigh?: number
* }} options
* @returns {number[]} - Filtered profile.
*
* @example
* const filtered = filterRoundnessProfile(profile, {
*   filterType: 'robust',
*   passType: 'bandpass',
*   cutoffLow: 15,
*   cutoffHigh: 150
* });
*/
export function filterRoundnessProfile(profile, {
 filterType,
 passType,
 cutoffLow = 0,
 cutoffHigh = Infinity,
}) {
 switch (filterType) {
   case 'gaussian':
     return gaussianFilter(profile, passType, cutoffLow, cutoffHigh);

   case 'fourier':
     return fourierFilter(profile, passType, cutoffLow, cutoffHigh);

   case 'robust':
     return robustFilter(profile, passType, cutoffLow, cutoffHigh);

   case '2cr':
     if (passType !== 'bandpass') {
       throw new Error(`2CR filter only supports bandpass mode`);
     }
     return twoCRFilter(profile, cutoffLow, cutoffHigh);

   default:
     throw new Error(`Unknown filterType: ${filterType}`);
 }
}