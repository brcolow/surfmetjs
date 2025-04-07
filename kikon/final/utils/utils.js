function numberToFixed(x, n) {
  return x.toFixed(n).replace(/\.?0+$/, '');
}

///////////////////
function _normalize({ x, y }) {
  const mag = Math.sqrt(x ** 2 + y ** 2);
  return { x: x / mag, y: y / mag };
}

function _distance2({ x: x1, y: y1 }, { x: x2, y: y2 }) {
  return Math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2);
}

function _distanceToSide({ x: x0, y: y0 }, side, onlyNormal = false) {
  // distance from point to segment
  const { x1, y1, x12, y12 } = side;
  const x01 = x0 - x1,
    y01 = y0 - y1,
    d12 = x12 ** 2 + y12 ** 2;
  const t = -(x01 * x12 + y01 * y12) / d12;
  if (t >= 0 && t <= 1) {
    return Math.abs(x01 * y12 - y01 * x12) / Math.sqrt(d12);
  } else if (onlyNormal) {
    return 1 / 0;
  }
  const { x2, y2 } = side;
  const x02 = x0 - x2,
    y02 = y0 - y2;
  return Math.min(
    Math.sqrt(x01 ** 2 + y01 ** 2),
    Math.sqrt(x02 ** 2 + y02 ** 2)
  );
}

function _normalized({ x, y }) {
  const mag = Math.sqrt(x ** 2 + y ** 2);
  return { x: x / mag, y: y / mag };
}

function _dotProduct({ x: x1, y: y1 }, { x: x2, y: y2 }) {
  return x1 * x2 + y1 * y2;
}

function _unitDir({ from: { x: x0, y: y0 }, to: { x, y } }) {
  return _normalized({ x: x - x0, y: y - y0 });
}

function _addToXY(o, { x, y }) {
  o.x += x;
  o.y += y;
}

function _setNormTo(v, mag) {
  const old_mag = Math.sqrt(v.x ** 2 + v.y ** 2);
  v.x = (v.x / old_mag) * mag;
  v.y = (v.y / old_mag) * mag;
}

function _rotateTo(v, theta) {
  const mag = Math.sqrt(v.x ** 2 + v.y ** 2);
  v.x = mag * Math.cos(theta);
  v.y = mag * Math.sin(theta);
}
