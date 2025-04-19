importScripts('../utils/utils.js');

self.onmessage = function ({
  data: {
    points,
    polygonVertices,
    polygonSides,
    maxRadius,
    maxCenter,
    subDivide,
    level,
  },
}) {
  ({ maxRadius, maxCenter } = computeAllDistancesToPolygon(
    points,
    polygonVertices,
    polygonSides,
    maxRadius,
    maxCenter
  ));
  self.postMessage({ points, maxRadius, maxCenter });
  const newPoints = subDivide
    ? computeNextLevelGrid(points, maxRadius, level)
    : [];
  self.postMessage({ newPoints });
};

function computeAllDistancesToPolygon(
  interiorPixels,
  polygonVertices,
  polygonSides,
  maxRadius,
  maxCenter
) {
  for (const pixel of interiorPixels) {
    const [x0, y0, r0] = pixel;
    let r;
    if (r0) {
      r = r0;
    } else {
      r = 1 / 0;
      for (const side of polygonSides) {
        const ri = _distanceToSide({ x: x0, y: y0 }, side);
        r = Math.min(r, ri);
      }
      pixel[2] = r;
    }
    if (r > maxRadius) {
      maxCenter = [x0, y0];
      maxRadius = r;
    }
  }
  return { maxRadius, maxCenter };
}

function computeNextLevelGrid(interiorPixels, maxRadius, level) {
  const res = 1 / 3 ** level;
  let relativeRadii = interiorPixels.map(([_1, _2, r], i) => [
    r / maxRadius,
    i,
  ]);
  const intendedCount = (600 * 600) / 9;
  let relativeRadiusThreshold = 0.7;
  while (relativeRadii.length > intendedCount) {
    relativeRadii = relativeRadii.filter(
      ([rFrac]) => rFrac > relativeRadiusThreshold
    );
    relativeRadiusThreshold =
      relativeRadiusThreshold + (1 - relativeRadiusThreshold) / 5;
  }

  let topPixels = relativeRadii.map(([_, i]) => interiorPixels[i]);
  const n = topPixels.length;
  for (let i = 0; i < n; i++) {
    const [x, y] = topPixels[i];
    const delta = res / 3;
    topPixels.push(
      [x - delta, y, null],
      [x + delta, y, null],
      [x, y - delta, null],
      [x, y + delta, null],
      [x - delta, y - delta, null],
      [x + delta, y - delta, null],
      [x - delta, y + delta, null],
      [x + delta, y + delta, null]
    );
  }

  return topPixels;
}
