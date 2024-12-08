/*!	Curve calc function for canvas 3.0.1
 *	(c) Epistemex 2013-2018, 2024
 *	License: MIT
 */

 'use strict';

 /**
  * Calculates an array containing points representing a cardinal spline through given point array.
  * Points must be arranged as: [x1, y1, x2, y2, ..., xn, yn].
  *
  * There must be a minimum of two points in the input array but the function
  * is only useful where there are three points or more.
  *
  * The points for the cardinal spline are returned as a new array.
  *
  * @param {Array} points - point array
  * @param {Number} [tension=0.5] - tension. Typically between [0.0, 1.0] but can be exceeded
  * @param {Number} [numOfSeg=25] - number of segments between two points (line resolution)
  * @param {Boolean} [closed=false] - Close the ends making the line continuous
  * @returns {Float32Array} New array with the calculated points that was added to the path
  */
 function getCurvePoints(points, tension = 0.5, numOfSeg = 25, closed = false) {
   if (typeof points === "undefined" || points.length < 2) {
     return new Float32Array(0);
   }
 
   if (numOfSeg < 1) {
     throw new Error('Number of segments cannot be less than one.');
   }
 
   let ptsClone;                                                       // for cloning point array
   let i = 1;
   let l = points.length;
   let resPos = 0;
   let cachePtr = 4;
 
   const resLength = (l - 2) * numOfSeg + 2 + (closed ? 2 * numOfSeg : 0);
   const result = new Float32Array(resLength);
   const cache = new Float32Array((numOfSeg + 2) << 2);
 
   ptsClone = points.slice(0);
 
   if (closed) {
     ptsClone.unshift(points[l - 1]);                                // insert end point as first point
     ptsClone.unshift(points[l - 2]);
     ptsClone.push(points[0], points[1]);                            // first point as last point
   }
   else {
     ptsClone.unshift(points[1]);                                    // copy 1. point and insert at beginning
     ptsClone.unshift(points[0]);
     ptsClone.push(points[l - 2], points[l - 1]);                    // duplicate end-points
   }
 
   // cache inner-loop calculations as they are based on t alone
   cache[0] = 1;                                                       // 1,0,0,0
 
   for(; i < numOfSeg; i++) {
     const st = i / numOfSeg;
     const st2 = st * st;
     const st3 = st2 * st;
     const st23 = st3 * 2;
     const st32 = st2 * 3;
 
     cache[cachePtr++] = st23 - st32 + 1;                            // c1
     cache[cachePtr++] = st32 - st23;                                // c2
     cache[cachePtr++] = st3 - 2 * st2 + st;                         // c3
     cache[cachePtr++] = st3 - st2;                                  // c4
   }
 
   cache[++cachePtr] = 1;                                              // 0,1,0,0
 
   // calc. points
 
   fill(ptsClone, cache, l, tension);
 
   if (closed) {
     //l = points.length;
     ptsClone = [];
     ptsClone.push(points[l - 4], points[l - 3],
       points[l - 2], points[l - 1],                               // second last and last
       points[0], points[1],
       points[2], points[3]);                                      // first and second
     fill(ptsClone, cache, 4, tension);
   }
 
   function fill(pts, cache, l, tension) {
     for(let i = 2; i < l; i += 2) {
       const pt1 = pts[i];                                         // x1
       const pt2 = pts[i + 1];                                     // y1
       const pt3 = pts[i + 2];                                     // x2
       const pt4 = pts[i + 3];                                     // y2
 
       const t1x = (pt3 - pts[i - 2]) * tension;                   // x2-x0
       const t1y = (pt4 - pts[i - 1]) * tension;                   // y2-y0
       const t2x = (pts[i + 4] - pt1) * tension;                   // x3-x1
       const t2y = (pts[i + 5] - pt2) * tension;                   // y3-y1
       let c = 0, c1, c2, c3, c4;
 
       for(let t = 0; t < numOfSeg; t++) {
         c1 = cache[c++];
         c2 = cache[c++];
         c3 = cache[c++];
         c4 = cache[c++];
 
         result[resPos++] = c1 * pt1 + c2 * pt3 + c3 * t1x + c4 * t2x;
         result[resPos++] = c1 * pt2 + c2 * pt4 + c3 * t1y + c4 * t2y;
       }
     }
   }
 
   // append last point
   l = closed ? 0 : points.length - 2;
   result[resPos++] = points[l++];
   result[resPos] = points[l];
 
   return result
 }

 function getCurvePoints2(points, tension = 0.5, numOfSeg = 25, closed = false) {
  if (points.length < 2) {
    return new Float32Array(0);
  }

  // Pre-calculate constants
  const numSegments = (points.length - 2) * numOfSeg + 2 + (closed ? 2 * numOfSeg : 0);
  const cacheSize = (numOfSeg + 2) << 2;
  const cache = new Float32Array(cacheSize);

  // Initialize cache
  cache[0] = 1;
  for (let i = 1; i < numOfSeg; i++) {
    const st = i / numOfSeg;
    const st2 = st * st;
    const st3 = st2 * st;
    const st23 = st3 * 2;
    const st32 = st2 * 3;
    cache[i * 4] = st23 - st32 + 1;
    cache[i * 4 + 1] = st32 - st23;
    cache[i * 4 + 2] = st3 - 2 * st2 + st;
    cache[i * 4 + 3] = st3 - st2;
  }
  cache[cacheSize - 4] = 1;

  // Calculate control points and fill result array
  const result = new Float32Array(numSegments);
  let resPos = 0;
  const adjustedPoints = closed
  ? [...points, ...points.slice(0, 2)] // Duplicate the first two points
  : [...points.slice(0, 2), ...points, ...points.slice(-2)];
  
  for (let i = 2; i < adjustedPoints.length; i += 2) {
    const pt1 = adjustedPoints[i - 2];
    const pt2 = adjustedPoints[i - 1];
    const pt3 = adjustedPoints[i];
    const pt4 = adjustedPoints[i + 1];
    const t1x = (pt3 - pt1) * tension;
    const t1y = (pt4 - pt2) * tension;
    const t2x = (adjustedPoints[i + 2] - pt1) * tension;
    const t2y = (adjustedPoints[i + 3] - pt2) * tension;

    for (let j = 0; j < numOfSeg; j++) {
      const c1 = cache[j * 4];
      const c2 = cache[j * 4 + 1];
      const c3 = cache[j * 4 + 2];
      const c4 = cache[j * 4 + 3];
      result[resPos++] = c1 * pt1 + c2 * pt3 + c3 * t1x + c4 * t2x;
      result[resPos++] = c1 * pt2 + c2 * pt4 + c3 * t1y + c4 * t2y;
    }
  }

  return result;
}
 
 export { getCurvePoints, getCurvePoints2 };
 