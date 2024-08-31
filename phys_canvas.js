addEventListener("DOMContentLoaded", (event) => {
  // const curvePoints = getCurvePoints(points.flatMap(pair => pair), 0.5, 25, true);
  const canvas = document.getElementById('myCanvas');
  const ctx = canvas.getContext('2d');

  // Set canvas dimensions
  canvas.width = 800;
  canvas.height = 600;
  // Polygon vertices
  const polygonVertices = ensureClosed(regularPolygonVertices(12, 100, {x: 100, y: 100}));
  const centroid = calculateCentroid(polygonVertices);

  const vAngle = Math.random() * 2 * Math.PI;
  const vMag = 30;
  const maxVelComp = 30;
  let previousTime;
  let msSinceLastInflation = 0;

  const ball = {
    x: centroid.x,
    y: centroid.y,
    radius: 5,
    vx: vMag * Math.cos(vAngle),
    vy: vMag * Math.sin(vAngle),
    lastLineCollidedWith: -1
  };
  

  function draw() {
    ctx.clearRect(0, 0, canvas.width, canvas.height);

    // Draw the polygon (will become the circle trace)
    ctx.beginPath();
    ctx.moveTo(polygonVertices[0].x, polygonVertices[0].y);
    for (let i = 1; i < polygonVertices.length; i++) {
      ctx.lineTo(polygonVertices[i].x, polygonVertices[i].y);
    }
    ctx.closePath();
    ctx.strokeStyle = 'black';
    ctx.lineWidth = 2;
    ctx.stroke();

    // Draw the ball
    ctx.fillStyle = 'red';
    ctx.beginPath();
    ctx.arc(ball.x, ball.y, ball.radius, 0, 2 * Math.PI);
    ctx.closePath();
    ctx.fill();
  }

  function handleCollisions() {
    // We need to add the first vertex to the end to create a full loop for intersection.
    const intersection = circlePolylineIntersection(ball, polygonVertices);
    // console.log(performance.measure('collision', 'start-collision', 'end-collision').duration);
    if (intersection != null && intersection.lineIndex != ball.lastLineCollidedWith) {

      const normal = { x: intersection.line.end.y - intersection.line.start.y, y: intersection.line.start.x - intersection.line.end.x };
      const normalMagnitude = Math.sqrt(normal.x * normal.x + normal.y * normal.y);
      const normalizedNormal = {
        x: normal.x / normalMagnitude,
        y: normal.y / normalMagnitude
      };

      const relativeVelocity = ball.vx * normalizedNormal.x + ball.vy * normalizedNormal.y;

      ball.lastLineCollidedWith = intersection.lineIndex;
      ball.vx = ball.vx - (2 * relativeVelocity * normalizedNormal.x);
      ball.vy = ball.vy - (2 * relativeVelocity * normalizedNormal.y);
      ball.vx = Math.max(-maxVelComp, Math.min(maxVelComp, ball.vx));
      ball.vy = Math.max(-maxVelComp, Math.min(maxVelComp, ball.vy));
    }
  }

  function animate(currentTime) {
    if (previousTime === undefined) {
      previousTime = currentTime;
    } else {
      const timePerFrame = currentTime - previousTime;
      msSinceLastInflation += timePerFrame;
      if (msSinceLastInflation > 1000) {
        // See if the ball will collide with a new larger radius.
        const biggerBallCollides = circlePolylineIntersection({ ...ball, radius: ball.radius + 1 }, polygonVertices) != null;
        if (biggerBallCollides != null) {
          ball.radius += 1;
          msSinceLastInflation = 0;
          console.log(biggerBallCollides);
        }
      }

      for (let i = 0; i < 1; i++) {
        handleCollisions();
        ball.x += ball.vx * (timePerFrame / 1000);
        ball.y += ball.vy * (timePerFrame / 1000);
      }

    }
    previousTime = currentTime;

    draw();
    requestAnimationFrame(animate);
  }

  animate();
});

function regularPolygonVertices(sides, radius, center = { x: 0, y: 0 }) {
  const vertices = [];
  const angle = 360 / sides;

  for (let i = 0; i < sides; i++) {
    const x = radius * Math.cos(Math.PI * i * angle / 180) + center.x;
    const y = radius * Math.sin(Math.PI * i * angle / 180) + center.y;
    vertices.push({x: x, y: y});
  }

  return vertices;
}


function circlePolylineIntersection(circle, polyline) {
  performance.mark("start-collision");
  let closestIntersection = null;
  let closestDistance = Infinity;

  for (let i = 0; i < polyline.length - 1; i++) {
    if (circle.lastLineCollidedWith == i) {
      continue;
    }
    const p1 = polyline[i];
    const p2 = polyline[i + 1];

    // Calculate distance from circle center to segment
    const distance = pointLineDistance(circle.x, circle.y, p1.x, p1.y, p2.x, p2.y);

    // Check if the distance is less than or equal to the circle's radius
    if (distance <= circle.radius && distance < closestDistance) {
      closestDistance = distance;
      closestIntersection = {
        lineIndex: i,
        line: { start: p1, end: p2 },
        point: pointOnLineClosestTo(circle.x, circle.y, p1.x, p1.y, p2.x, p2.y)
      };
    }
  }

  performance.mark("end-collision");
  return closestIntersection ? closestIntersection : null;
}

function pointLineDistance(px, py, x1, y1, x2, y2) {
  // Check if the point lies on the line segment
  if (isPointOnLineSegment(px, py, x1, y1, x2, y2)) {
    return 0;
  }

  // Calculate the perpendicular distance
  const dx = x2 - x1;
  const dy = y2 - y1;
  const u = ((px - x1) * dx + (py - y1) * dy) / (dx * dx + dy * dy);

  if (u < 0 || u > 1) {
    // The closest point is one of the endpoints
    const d1 = Math.sqrt((px - x1) * (px - x1) + (py - y1) * (py - y1));
    const d2 = Math.sqrt((px - x2) * (px - x2) + (py - y2) * (py - y2));
    return Math.min(d1, d2);
  }

  const ix = x1 + u * dx;
  const iy = y1 + u * dy;
  return Math.sqrt((px - ix) * (px - ix) + (py - iy) * (py - iy));
}

function pointOnLineClosestTo(px, py, x1, y1, x2, y2) {
  const dx = x2 - x1;
  const dy = y2 - y1;
  const u = ((px - x1) * dx + (py - y1) * dy) / (dx * dx + dy * dy);

  if (u < 0) {
    return { x: x1, y: y1 };
  } else if (u > 1) {
    return { x: x2, y: y2 };
  } else {
    return { x: x1 + u * dx, y: y1 + u * dy };
  }
}

function isPointOnLineSegment(px, py, x1, y1, x2, y2) {
  const v1x = px - x1;
  const v1y = py - y1;
  const v2x = x2 - x1;
  const v2y = y2 - y1;

  const cross = v1x * v2y - v1y * v2x;
  if (cross !== 0) {
    return false; // Point is not on the line
  }

  const dot = v1x * v2x + v1y * v2y;
  return dot >= 0 && dot <= v2x * v2x + v2y * v2y;
}

function calculateCentroid(polygon) {
  let signedArea = 0;
  let cx = 0;
  let cy = 0;

  for (let i = 0; i < polygon.length - 1; i++) {
    const xi = polygon[i].x;
    const yi = polygon[i].y;
    const xi1 = polygon[i + 1].x;
    const yi1 = polygon[i + 1].y;

    const factor = xi * yi1 - xi1 * yi;
    signedArea += factor;
    cx += (xi + xi1) * factor;
    cy += (yi + yi1) * factor;
  }

  signedArea *= 0.5;
  cx /= (6 * signedArea);
  cy /= (6 * signedArea);

  return { x: cx, y: cy };
}

function ensureClosed(polygon) {
  const clonedPolygon = [...polygon]; // Create a copy of the polygon
  // Ensure polygon is a closed polygon
  if (polygon[0] !== polygon[polygon.length - 1]) {
    clonedPolygon.push(polygon[0]);
  }

  return clonedPolygon;
}