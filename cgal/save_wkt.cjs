// save_wkt.js
const fs = require("fs");

function polygonToWKT(points) {
  if (points.length < 3) throw new Error("Need at least 3 points.");
  const first = points[0];
  const last  = points[points.length - 1];
  const closed = (first[0] === last[0] && first[1] === last[1])
    ? points
    : [...points, first];
  const ring = closed.map(([x, y]) => `${x} ${y}`).join(", ");
  return `POLYGON((${ring}))\n`;
}

// Example usage:
// node save_wkt.js '[[-1,0],[2,0],[2,1],[0,2],[-1,1]]' poly.wkt
const pts = JSON.parse(process.argv[2]);
const out = process.argv[3] || "polygon.wkt";
fs.writeFileSync(out, polygonToWKT(pts), "utf8");
console.log(`Wrote ${out}`);
