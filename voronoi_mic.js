import { computeConvexHull, isPointInPolygon } from './utils.js';

/**
 * Computes the Maximum Inscribed Circle (MIC) for a set of 2D points.
 * 
 * Method:
 * - Compute Delaunay triangulation
 * - Build Voronoi diagram
 * - Filter Voronoi vertices inside the convex hull
 * - Find the vertex farthest from the point set
 * 
 * @param {Array<Array<number>>} points - Array of [x, y] points
 * @returns {{ center: Vertex, radius: number }} Center and radius of MIC
 */
function computeMIC(points) {
  const vertices = points.map(([x, y]) => new Vertex(x, y, 0));
  const delaunayTriangles = bowyerWatson(vertices);
  const voronoi = voronoi(delaunayTriangles);

  const hull = computeConvexHull(points);

  let bestCenter = null;
  let bestRadius = -Infinity;

  for (const centerVertex of voronoi.vertices) {
    const center = [centerVertex.x, centerVertex.y];

    if (isPointInPolygon(center, hull)) {
      let minDist = Infinity;
      for (const p of points) {
        const dx = p[0] - center[0];
        const dy = p[1] - center[1];
        const dist = Math.hypot(dx, dy);
        if (dist < minDist) {
          minDist = dist;
        }
      }

      if (minDist > bestRadius) {
        bestRadius = minDist;
        bestCenter = center;
      }
    }
  }

  return { center: bestCenter, radius: bestRadius };
}