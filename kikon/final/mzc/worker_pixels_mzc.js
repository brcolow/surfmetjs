importScripts('../utils/utils.js');

self.onmessage = function ({
  data: {
    points,
    polygonVertices,
    polygonSides,
    minZone,
    minZoneCenter,
    minZoneRExt,
    minZoneRInt,
    subDivide,
    level,
  },
}) {
  ({ minZone, minZoneCenter, minZoneRInt, minZoneRExt } =
    computeAllDistancesToPolygon(
      points,
      polygonVertices,
      polygonSides,
      minZone,
      minZoneCenter,
      minZoneRExt,
      minZoneRInt
    ));
  self.postMessage({
    points,
    minZone,
    minZoneCenter,
    minZoneRInt,
    minZoneRExt,
  });
  const newPoints = subDivide
    ? computeNextLevelGrid(points, minZone, level)
    : [];
  self.postMessage({ newPoints });
};

function computeAllDistancesToPolygon(
  interiorPixels,
  polygonVertices,
  polygonSides,
  minZone,
  minZoneCenter,
  minZoneRExt,
  minZoneRInt
) {
  for (const pixel of interiorPixels) {
    const [x0, y0] = pixel;

    let r_int = 1 / 0;
    for (const side of polygonSides) {
      const ri = _distanceToSide({ x: x0, y: y0 }, side);
      r_int = Math.min(r_int, ri);
    }
    let r_ext = 0;
    for (const vertex of polygonVertices) {
      const ri = _distance2({ x: x0, y: y0 }, vertex);
      r_ext = Math.max(r_ext, ri);
    }
    const dZone = r_ext - r_int;
    pixel[2] = dZone;

    if (dZone < minZone) {
      minZoneCenter = [x0, y0];
      minZone = dZone;
      minZoneRExt = r_ext;
      minZoneRInt = r_int;
    }
  }
  return { minZone, minZoneCenter, minZoneRInt, minZoneRExt };
}

function computeNextLevelGrid(interiorPixels, minZone, level) {
  const res = 1 / 3 ** level;
  let relativeZones = interiorPixels.map(([_1, _2, zone], i) => [
    minZone / zone,
    i,
  ]);
  const intendedCount = (600 * 600) / 9;
  let relativeThreshold = 0.7;
  while (relativeZones.length > intendedCount) {
    relativeZones = relativeZones.filter(
      ([zoneFrac]) => zoneFrac > relativeThreshold
    );
    relativeThreshold = relativeThreshold + (1 - relativeThreshold) / 5;
  }

  let topPixels = relativeZones.map(([_, i]) => interiorPixels[i]);
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
