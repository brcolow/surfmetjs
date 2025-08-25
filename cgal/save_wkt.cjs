// save_wkt.js
const fsprom = require('fs').promises;
const fs = require('fs')

/**
 * Reads point data from a text file. Each line in the file 
 * is expected to contain coordinates separated by whitespace.
 *
 * @param {string} file - T text file containing point data.
 * @param {number} [pointsPerLine=2] - The number of coordinates per line. Supported values:
 *   - 2: For 2D points (x, y).
 *   - 3: For 3D points (x, y, z).
 *
 * @returns {Promise<Array>} A promise that resolves to an array of points, where 
 * each point is represented as an array of numbers ([x, y] for 2D or [x, y, z] for 3D).
 * Returns `null` if an error occurs.
 *
 * @throws {Error} If an unsupported number of points per line is specified.
 */
async function readPointsFromFile(file, pointsPerLine = 2) {
  try {
    const response = await fsprom.readFile(file, { encoding: 'utf8' });
    const points = [];
    for (const line of response.split('\n')) {
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
      } else {
        throw new Error('Unsupported number of points per line: ' + pointsPerLine);
      }
    }
    return points;
  } catch (error) {
    console.error("Error fetching data:", error);
    return null;
  }
}

function polygonToWKT(points) {
  if (points.length < 3) throw new Error("Need at least 3 points.");
  const first = points[0];
  const last = points[points.length - 1];
  const closed = (first[0] === last[0] && first[1] === last[1])
    ? points
    : [...points, first];
  const ring = closed.map(([x, y]) => `${x} ${y}`).join(", ");
  return `POLYGON((${ring}))\n`;
}

// Example usage:
// node save_wkt.js points.txt poly.wkt
readPointsFromFile(process.argv[2])
  .then(points => {
    if (points) {
      const out = process.argv[3] || "polygon.wkt";
      fs.writeFileSync(out, polygonToWKT(points), "utf8");
      console.log(`Wrote ${out}`);
    }
  });
