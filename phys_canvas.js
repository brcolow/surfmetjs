addEventListener("DOMContentLoaded", (event) => {
  // const curvePoints = getCurvePoints(points.flatMap(pair => pair), 0.5, 25, true);
  const canvas = document.getElementById('myCanvas');
  const ctx = canvas.getContext('2d');

  // Set canvas dimensions
  canvas.width = 800;
  canvas.height = 600;
  // Polygon vertices
  const polygonVertices = ensureClosed(regularPolygonVertices(12, 100, { x: 100, y: 100 }));

  const insertAt = 5;
  const xNew = polygonVertices[insertAt].x / 3 + polygonVertices[insertAt + 1].x / 3 + 100 / 3,
    yNew = polygonVertices[insertAt].y / 3 + polygonVertices[insertAt + 1].y / 3 + 100 / 3;
  polygonVertices.splice(insertAt, 0, { x: xNew, y: yNew });

  const centroid = calculateCentroid(polygonVertices);

  const vAngle = Math.random() * 2 * Math.PI;
  const vMag = 30;
  const maxVelComp = 30;
  let previousTime;
  let msSinceLastInflation = 0;
  const ball = {
    // Start ball at the centroid of the polyline.
    x: centroid.x,
    y: centroid.y,
    radius: 5,
    vx: vMag * Math.cos(vAngle),
    vy: vMag * Math.sin(vAngle),
    lastLineCollidedWith: -1,
    vRadius: 1,
  };
  let impactedLineIndex = -1;
  let speedFactor = 1;


  function draw() {
    ctx.clearRect(0, 0, canvas.width, canvas.height);

    // Draw the polygon (will become the circle trace)
    for (let i = 1; i < polygonVertices.length; i++) {
      ctx.beginPath();
      ctx.moveTo(polygonVertices[i - 1].x, polygonVertices[i - 1].y);
      ctx.lineTo(polygonVertices[i].x, polygonVertices[i].y);
      ctx.closePath();
      if (impactedLineIndex === i) {
        ctx.strokeStyle = 'green';
      } else {
        ctx.strokeStyle = 'black';
      }
      ctx.lineWidth = 2;
      ctx.stroke();
    }


    // Draw the ball
    ctx.fillStyle = 'red';
    ctx.beginPath();
    ctx.arc(ball.x, ball.y, ball.radius, 0, 2 * Math.PI);
    ctx.closePath();
    ctx.fill();

    if (xComputedCollision !== null && yComputedCollision !== null) {
      ctx.strokeStyle = 'green';
      ctx.lineWidth = 2;
      ctx.beginPath();
      ctx.arc(xComputedCollision, yComputedCollision, 2, 0, 2 * Math.PI);
      ctx.closePath();
      ctx.stroke();
    }
  }

  function handleCollisions() {
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
      return true;
    }
    return false;
  }

  let xComputedCollision = null;
  let yComputedCollision = null;

  function computeNextImpact() {
    // At the beginning when the ball is at the centroid of the polyline or after it has collided with a line segment of the polyline we want to know
    // which line it will hit next, where it will hit it and when it will it.
    // We can analyze this problem and come up with a solution by utilizing parametric equations:

    // Let the ball radius be: r = r0 + rr * t
    //     the ball's current center be (xc, yc)
    //     the ball's velocity be (vx, vy)
    //     Then the ball's center after time t is: (xc + vx * t, yc + vy * t)

    // Parametric equation of the circumference of inflating ball:
    // x = xc + vx * t + (r0 + rr * t) * cos(phi),   (1)
    // y = yc + vy * t + (r0 + rr * t) * sin(phi)    (2)
    // t > 0; -pi < phi < pi

    // Parametric equation of the current line segment of the polyline:
    // x = x1 + (x2 - x1) * s,  (3)
    // y = y1 + (y2 - y1) * s,  (4)
    // 0 <= s <= 1

    // We want to solve for the moment when the ball hits the line segment and is thus tangent with it - in this case the x and y correspond
    // to the same points of the circle and the line segment so we get a system of two equations in s, t, and phi. We can solve this for s an t
    // keeping phi a free parameter. Also because the polyline is not necessarily convex we also need to solve the equations s(phi) = 0 and s(phi) = 1
    // and find a formula for the values of phi for the extrema of t(phi). The minimum of t gives us the moment when the contact happens. (If we were
    // to ignore the contact and continue to draw the ball going in the same direction, the contact will continue while the ball traverses the segment;
    // those contacts will also obey our equations. However, we are only interested in the first moment of the contact which is defined by min t(phi)).
    // Once we have phi for the contact, we can compute s(phi) and use it to determine the exact point of the contact.

    // All of these can be solved in closed form - numerical approximation is not necessary. We can find the solutons using sympy:

    /*
    import sympy as sym

    (xc, yc, vx, vy, r0, rr, x1, y1, x2, y2, t, s) = sym.symbols("xc yc vx vy r0 rr x1 y1 x2 y2 t s")
    phi = sym.Symbol("phi", real=True)

    eqx1 = xc + vx * t + (r0 + rr * t) * sym.cos(phi)
    eqy1 = yc + vy * t + (r0 + rr * t) * sym.sin(phi)

    eqx2 = x1 + (x2 - x1) * s
    eqy2 = y1 + (y2 - y1) * s

    sol = sym.solve([eqx1 - eqx2, eqy1 - eqy2], [t, s])
    solS = sol[s]
    print('s(phi) =', sym.together(solS))
    solT = sol[t]
    print('t(phi) =', sym.together(solT))

    print(80 * '-')
    print('Solutions for s(phi) = 0 ...')
    sol_s0 = sym.solve(solS, phi)
    print(sol_s0)

    print('Solutions for s(phi) = 1 ...')
    sol_s1 = sym.solve(solS - 1, phi)
    print(sol_s1)

    print(80 * '-')

    print('The values of phi for the extremes of t(phi) ...')
    dt = sym.diff(solT, phi)
    sol_t = sym.solve(dt, phi)
    print(sol_t)
    */

    // Running this code prints the following (after a few minutes!):

    // s(phi) = (r0*vx*sin(phi) - r0*vy*cos(phi) + rr*x1*sin(phi) - rr*xc*sin(phi) - rr*y1*cos(phi) + rr*yc*cos(phi) - vx*y1 + vx*yc + vy*x1 - vy*xc)/(rr*x1*sin(phi) - rr*x2*sin(phi) - rr*y1*cos(phi) + rr*y2*cos(phi) - vx*y1 + vx*y2 + vy*x1 - vy*x2)
    // t(phi) = (-r0*x1*sin(phi) + r0*x2*sin(phi) + r0*y1*cos(phi) - r0*y2*cos(phi) + x1*y2 - x1*yc - x2*y1 + x2*yc + xc*y1 - xc*y2)/(rr*x1*sin(phi) - rr*x2*sin(phi) - rr*y1*cos(phi) + rr*y2*cos(phi) - vx*y1 + vx*y2 + vy*x1 - vy*x2)
    // --------------------------------------------------------------------------------
    // Solutions for s(phi) = 0 ...
    // [-2*atan((r0*vx + rr*x1 - rr*xc - sqrt(r0**2*vx**2 + r0**2*vy**2 + 2*r0*rr*vx*x1 - 2*r0*rr*vx*xc + 2*r0*rr*vy*y1 - 2*r0*rr*vy*yc + rr**2*x1**2 - 2*rr**2*x1*xc + rr**2*xc**2 + rr**2*y1**2 - 2*rr**2*y1*yc + rr**2*yc**2 - vx**2*y1**2 + 2*vx**2*y1*yc - vx**2*yc**2 + 2*vx*vy*x1*y1 - 2*vx*vy*x1*yc - 2*vx*vy*xc*y1 + 2*vx*vy*xc*yc - vy**2*x1**2 + 2*vy**2*x1*xc - vy**2*xc**2))/(r0*vy + rr*y1 - rr*yc - vx*y1 + vx*yc + vy*x1 - vy*xc)), -2*atan((r0*vx + rr*x1 - rr*xc + sqrt(r0**2*vx**2 + r0**2*vy**2 + 2*r0*rr*vx*x1 - 2*r0*rr*vx*xc + 2*r0*rr*vy*y1 - 2*r0*rr*vy*yc + rr**2*x1**2 - 2*rr**2*x1*xc + rr**2*xc**2 + rr**2*y1**2 - 2*rr**2*y1*yc + rr**2*yc**2 - vx**2*y1**2 + 2*vx**2*y1*yc - vx**2*yc**2 + 2*vx*vy*x1*y1 - 2*vx*vy*x1*yc - 2*vx*vy*xc*y1 + 2*vx*vy*xc*yc - vy**2*x1**2 + 2*vy**2*x1*xc - vy**2*xc**2))/(r0*vy + rr*y1 - rr*yc - vx*y1 + vx*yc + vy*x1 - vy*xc))]
    // Solutions for s(phi) = 1 ...
    // [-2*atan((r0*vx + rr*x2 - rr*xc - sqrt(r0**2*vx**2 + r0**2*vy**2 + 2*r0*rr*vx*x2 - 2*r0*rr*vx*xc + 2*r0*rr*vy*y2 - 2*r0*rr*vy*yc + rr**2*x2**2 - 2*rr**2*x2*xc + rr**2*xc**2 + rr**2*y2**2 - 2*rr**2*y2*yc + rr**2*yc**2 - vx**2*y2**2 + 2*vx**2*y2*yc - vx**2*yc**2 + 2*vx*vy*x2*y2 - 2*vx*vy*x2*yc - 2*vx*vy*xc*y2 + 2*vx*vy*xc*yc - vy**2*x2**2 + 2*vy**2*x2*xc - vy**2*xc**2))/(r0*vy + rr*y2 - rr*yc - vx*y2 + vx*yc + vy*x2 - vy*xc)), -2*atan((r0*vx + rr*x2 - rr*xc + sqrt(r0**2*vx**2 + r0**2*vy**2 + 2*r0*rr*vx*x2 - 2*r0*rr*vx*xc + 2*r0*rr*vy*y2 - 2*r0*rr*vy*yc + rr**2*x2**2 - 2*rr**2*x2*xc + rr**2*xc**2 + rr**2*y2**2 - 2*rr**2*y2*yc + rr**2*yc**2 - vx**2*y2**2 + 2*vx**2*y2*yc - vx**2*yc**2 + 2*vx*vy*x2*y2 - 2*vx*vy*x2*yc - 2*vx*vy*xc*y2 + 2*vx*vy*xc*yc - vy**2*x2**2 + 2*vy**2*x2*xc - vy**2*xc**2))/(r0*vy + rr*y2 - rr*yc - vx*y2 + vx*yc + vy*x2 - vy*xc))]
    // --------------------------------------------------------------------------------
    // The values of phi for the extremes of t(phi) ...
    // [-2*atan((-y1 + y2 + sqrt(x1**2 - 2*x1*x2 + x2**2 + y1**2 - 2*y1*y2 + y2**2))/(x1 - x2)), 2*atan((y1 - y2 + sqrt(x1**2 - 2*x1*x2 + x2**2 + y1**2 - 2*y1*y2 + y2**2))/(x1 - x2))]

    console.log(`ball = ${JSON.stringify(ball)}; impactedLineIndex = ${impactedLineIndex};`)
    // this data (and the Pause button) allows us to restart the animation from the previous
    // impact, if an erroneous computation is detected visually

    xComputedCollision = null; yComputedCollision = null;
    let shortestTimeTillImpact = 1 / 0;
    let nextImpactedLineIndex = -1;
    const r0 = ball.radius;
    const rr = ball.vRadius;
    const vx = ball.vx;
    const vy = ball.vy;
    for (let i = 1; i < polygonVertices.length; i++) {
      if (i === impactedLineIndex) { // exclude current line
        continue;
      }
      const p1 = polygonVertices[i - 1];
      const p2 = polygonVertices[i];
      const x1 = p1.x;
      const y1 = p1.y;
      const x2 = p2.x;
      const y2 = p2.y;
      const xc = ball.x;
      const yc = ball.y;

      const x1c = x1 - xc, y1c = y1 - yc,
        x2c = x2 - xc, y2c = y2 - yc,
        x12 = x1 - x2, y12 = y1 - y2;

      const sA = vy * x1c - vx * y1c,
        sB = - r0 * vy - rr * y1c,
        sC = r0 * vx + rr * x1c,
        sD = vy * x12 - vx * y12,
        sE = - rr * y12,
        sF = rr * x12;

      const getS = phi => {
        const cos_phi = Math.cos(phi),
          sin_phi = Math.sin(phi);
        return (sA + sB * cos_phi + sC * sin_phi) / (sD + sE * cos_phi + sF * sin_phi);
      }

      const tA = x1 * y2c - x2 * y1c + xc * y12,
        tB = r0 * y12,
        tC = - r0 * x12,
        tD = vy * x12 - vx * y12,
        tE = -rr * y12,
        tF = rr * x12;

      const getT = phi => {
        const cos_phi = Math.cos(phi),
          sin_phi = Math.sin(phi);
        return (tA + tB * cos_phi + tC * sin_phi) / (tD + tE * cos_phi + tF * sin_phi);
      }

      const phi_t1 = -2 * Math.atan((-y12 + Math.sqrt(x12 ** 2 + y12 ** 2)) / x12),
        phi_t2 = 2 * Math.atan((y12 + Math.sqrt(x12 ** 2 + y12 ** 2)) / x12),
        t1 = getT(phi_t1),
        t2 = getT(phi_t2);
      let tMin12 = Math.min(...[t1, t2]);//.filter(t => t > 0));
      if (tMin12 < shortestTimeTillImpact) {
        const phi4Min = (tMin12 === t1) ? phi_t1 : phi_t2;
        let s4Min = getS(phi4Min);
        if (s4Min < 0) {
          const expr1 = r0 ** 2 * (vx ** 2 + vy ** 2) + 2 * r0 * rr * (vx * x1c + vy * y1c) +
            rr ** 2 * (x1c ** 2 + y1c ** 2) - (vy * x1c - vx * y1c) ** 2;
          let phi_s0_1 = -2 * Math.atan((r0 * vx + rr * x1c - Math.sqrt(expr1)) / ((rr - vx) * y1c + vy * (r0 + x1c)));
          let phi_s0_2 = -2 * Math.atan((r0 * vx + rr * x1c + Math.sqrt(expr1)) / ((rr - vx) * y1c + vy * (r0 + x1c)));
          let t0_1 = 1 / 0, t0_2 = 1 / 0;
          if (Number.isFinite(phi_s0_1)) {
            const s0_1 = getS(phi_s0_1);
            if (Math.abs(s0_1) < 1e-6) { // s0_1 ~= 0
              // slightly change phi to help compare the two lines that share a point
              phi_s0_1 += 1e-3;
              const s0_1a = getS(phi_s0_1);
              if (s0_1a < s0_1) {
                phi_s0_1 -= 2e-3;
              }
              t0_1 = getT(phi_s0_1);
              if (t0_1 < 0) {
                t0_1 = 1 / 0;
              }
            }
          }

          if (Number.isFinite(phi_s0_2)) {
            const s0_2 = getS(phi_s0_2);
            if (Math.abs(s0_2) < 1e-6) { // s0_2 ~= 0
              phi_s0_2 += 1e-3;
              const s0_2a = getS(phi_s0_2);
              if (s0_2a < s0_2) {
                phi_s0_2 -= 2e-3;
              }
              t0_2 = getT(phi_s0_2);
              if (t0_2 < 0) {
                t0_2 = 1 / 0;
              }
            }
          }
          const tCorrected = Math.min(t0_1, t0_2);
          if (tCorrected < shortestTimeTillImpact) {
            shortestTimeTillImpact = tCorrected;
            xComputedCollision = x1;
            yComputedCollision = y1;
            nextImpactedLineIndex = i;
          }
        } else if (s4Min > 1) {
          const expr2 = r0 ** 2 * (vx ** 2 + vy ** 2) + 2 * r0 * rr * (vx * x2c + vy * y2c) + rr ** 2 * (x2c ** 2 + y2c ** 2) - (vy * x2c - vx * y2c) ** 2;
          let phi_s1_1 = -2 * Math.atan((r0 * vx + rr * x2c - Math.sqrt(expr2)) / ((rr - vx) * y2c + vy * (r0 + x2c)));
          let phi_s1_2 = -2 * Math.atan((r0 * vx + rr * x2c + Math.sqrt(expr2)) / ((rr - vx) * y2c + vy * (r0 + x2c)));
          let t1_1 = 1 / 0, t1_2 = 1 / 0;
          if (Number.isFinite(phi_s1_1)) {
            const s1_1 = getS(phi_s1_1);
            if (Math.abs(s1_1 - 1) < 1e-6) { // s1_1 ~= 1
              phi_s1_1 += 1e-3;
              const s1_1a = getS(phi_s1_1);
              if (s1_1a > s1_1) {
                phi_s1_1 -= 2e-3;
              }
              t1_1 = getT(phi_s1_1);
              if (t1_1 < 0) {
                t1_1 = 1 / 0;
              }
            }
          }

          if (Number.isFinite(phi_s1_2)) {
            const s1_2 = getS(phi_s1_2);
            if (Math.abs(s1_2 - 1) < 1e-6) {  // s1_2 ~= 1
              phi_s1_2 += 1e-3;
              const s1_2a = getS(phi_s1_2);
              if (s1_2a > s1_2) {
                phi_s1_2 -= 2e-3;
              }
              t1_2 = getT(phi_s1_2);
              if (t1_2 < 0) {
                t1_2 = 1 / 0;
              }
            }
          }
          const tCorrected = Math.min(t1_1, t1_2);
          if (tCorrected < shortestTimeTillImpact) {
            shortestTimeTillImpact = tCorrected;
            xComputedCollision = x2;
            yComputedCollision = y2;
            nextImpactedLineIndex = i;
          }
        } else if (tMin12 > 0) {
          shortestTimeTillImpact = tMin12;
          xComputedCollision = x1 + (x2 - x1) * s4Min;
          yComputedCollision = y1 + (y2 - y1) * s4Min;
          nextImpactedLineIndex = i;
        }
      }
    }

    impactedLineIndex = nextImpactedLineIndex;
    console.log('next impactedLineIndex =', impactedLineIndex);
  }

  draw();
  requestAnimationFrame(animate);
  let collision = false;
  
  function animate(currentTime) {
    if (previousTime === undefined) {
      setTimeout(computeNextImpact);
    } else {
      const timePerFrame = currentTime - previousTime;
      msSinceLastInflation += timePerFrame;

      ball.radius += ball.vRadius * (timePerFrame / 1000 / speedFactor);

      collision = handleCollisions();
      ball.x += ball.vx * (timePerFrame / 1000 / speedFactor);
      ball.y += ball.vy * (timePerFrame / 1000 / speedFactor);
    }
    previousTime = currentTime;
    draw();
    if (collision) {
      // next impact should only be computed after a collision
      // setTimeout(computeNextImpact) // asynchronous
      computeNextImpact();
      if (ball.radius < 100 * Math.cos(Math.PI / (polygonVertices.length - 1)) - 1) {
        requestAnimationFrame(animate);
      }
    } else {
      requestAnimationFrame(animate);
    }
  }
});

function regularPolygonVertices(sides, radius, center = { x: 0, y: 0 }) {
  const vertices = [];
  const angle = 360 / sides;

  for (let i = 0; i < sides; i++) {
    const x = radius * Math.cos(Math.PI * i * angle / 180) + center.x;
    const y = radius * Math.sin(Math.PI * i * angle / 180) + center.y;
    vertices.push({ x: x, y: y });
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