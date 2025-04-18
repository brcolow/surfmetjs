importScripts('../utils/CircularZones.js', '../utils/utils.js');

let polygonVertices,
  polygonSides,
  previousCollisionBall,
  vMag,
  collisionCounter,
  size0,
  RExt,
  ballMoving,
  ballFinalized;

self.onmessage = function ({ data }) {
  if (data.type === 'init') {
    ({ polygonVertices, previousCollisionBall, vMag, size0, RExt } = data);
    polygonSides = buildPolygonData(polygonVertices);
    collisionCounter = 0;
    ballMoving = true;
    ballFinalized = false;
    self.postMessage({ type: 'postInit', polygonVertices, polygonSides });
  } else if (data.type === 'computeNextCollision') {
    if (!ballMoving) {
      return; // don't answer
    }
    self.postMessage({
      type: 'nextCollisionBall',
      nextCollisionBall: computeNextImpact(),
    });
  } else if (data.type === 'finalizeBall') {
    if (ballMoving || ballFinalized) {
      return;
    }
    ballFinalized = true;
    self.postMessage({ type: 'ballFinally', ballFinally: finalizeBall() });
  }
};

function buildPolygonData(polygonVertices) {
  // precomputes all static (that don't depend on the ball) properties of the polygon vertices and sides
  const nVertices = polygonVertices.length;
  for (let i = 1; i <= nVertices; i++) {
    const vertex = polygonVertices[i - 1],
      nextVertex = polygonVertices[i === nVertices ? 0 : i];
    Object.assign(vertex, {
      next: nextVertex,
      index: i - 1,
    });
    nextVertex.previous = vertex;
  }
  const sides = [];
  for (let i = 1; i <= nVertices; i++) {
    const p1 = polygonVertices[i - 1],
      p2 = p1.next,
      x1 = p1.x,
      y1 = p1.y,
      x2 = p2.x,
      y2 = p2.y,
      x12 = x1 - x2,
      y12 = y1 - y2,
      normal = _normalize({ x: y12, y: -x12 });
    sides.push({
      index: i - 1,
      p1,
      p2,
      x1,
      y1,
      x2,
      y2,
      x12,
      y12,
      normal,
    });
  }
  for (let i = 1; i <= nVertices; i++) {
    const side = sides[i - 1],
      sideNext = sides[i === nVertices ? 0 : i];
    side.next = sideNext;
    sideNext.previous = side;
  }
  return sides;
}

function computeAllImpacted(
  { x: xcNext, y: ycNext, radius: rNext },
  threshold
) {
  // called after computeNextImpact has found the nearest impact, to find all
  // other impact points (vertices/sides) within threshold (radius-wise)
  let nTotal = 0;
  const impactedVertexIndices = [];
  for (const vertex of polygonVertices) {
    const ri = Math.sqrt((vertex.x - xcNext) ** 2 + (vertex.y - ycNext) ** 2);
    if (Math.abs(ri - rNext) < threshold) {
      impactedVertexIndices.push(vertex.index);
      nTotal++;
    }
  }
  const impactedSideIndices = [];
  const center = { x: xcNext, y: ycNext };
  for (const side of polygonSides) {
    const ri = _distanceToSide(center, side, true);
    if (/*Number.isFinite(ri) &&*/ Math.abs(ri - rNext) < threshold) {
      impactedSideIndices.push(side.index);
      nTotal++;
    }
  }

  let allowedDirections = null;
  if (nTotal > 1) {
    allowedDirections = CircularZones.fullCircle();
    for (const vertexIndex of impactedVertexIndices) {
      const vertex = polygonVertices[vertexIndex];
      const dirCenter = Math.atan2(-vertex.y + ycNext, -vertex.x + xcNext),
        dir1 = Math.atan2(vertex.x - xcNext, -vertex.y + ycNext),
        dir2 = Math.atan2(-vertex.x + xcNext, vertex.y - ycNext);
      allowedDirections.intersectWith(
        CircularZones.betweenDirectionsContainingDirection(
          dir1,
          dir2,
          dirCenter
        )
      );
    }
    for (const sideIndex of impactedSideIndices) {
      const side = polygonSides[sideIndex];
      const p1 = side.p1,
        p2 = side.p2,
        p0 = { x: (p1.x + p2.x) / 2, y: (p1.y + p2.y) / 2 };
      const dir1 = Math.atan2(p1.y - p0.y, p1.x - p0.x),
        dir2 = Math.atan2(p2.y - p0.y, p2.x - p0.x),
        dirCenter = Math.atan2(center.y - p0.y, center.x - p0.x);
      allowedDirections.intersectWith(
        CircularZones.betweenDirectionsContainingDirection(
          dir1,
          dir2,
          dirCenter
        )
      );
    }
  }

  return {
    impactedVertexIndices,
    impactedSideIndices,
    nTotal,
    allowedDirections,
  };
}

function computeNextImpact() {
  let xComputedCollision = null,
    yComputedCollision = null;
  let timeTillNextImpact = 1 / 0;
  let phi4Impact = -1;
  let impactedSideIndex = -1;
  let impactedSideNormal = null;
  let impactedVertexIndex = -1;
  const r0 = previousCollisionBall.radius;
  const rr = previousCollisionBall.vRadius;
  const vx = previousCollisionBall.vx;
  const vy = previousCollisionBall.vy;
  for (const side of polygonSides) {
    if (
      previousCollisionBall.impactedSideIndex === side.index ||
      [side.p1.index, side.p2.index].includes(
        previousCollisionBall.impactedVertexIndex
      )
    ) {
      // exclude current line
      continue;
    }
    const { x1, y1, x2, y2, x12, y12 } = side;
    const xc = previousCollisionBall.x;
    const yc = previousCollisionBall.y;

    const x1c = x1 - xc,
      y1c = y1 - yc,
      x2c = x2 - xc,
      y2c = y2 - yc;

    const sA = vy * x1c - vx * y1c,
      sB = -r0 * vy - rr * y1c,
      sC = r0 * vx + rr * x1c,
      sD = vy * x12 - vx * y12,
      sE = -rr * y12,
      sF = rr * x12;

    const getS = (phi) => {
      const cos_phi = Math.cos(phi),
        sin_phi = Math.sin(phi);
      return (
        (sA + sB * cos_phi + sC * sin_phi) / (sD + sE * cos_phi + sF * sin_phi)
      );
    };

    const tA = x1 * y2c - x2 * y1c + xc * y12,
      tB = r0 * y12,
      tC = -r0 * x12,
      tD = vy * x12 - vx * y12,
      tE = -rr * y12,
      tF = rr * x12;

    const getT = (phi) => {
      const cos_phi = Math.cos(phi),
        sin_phi = Math.sin(phi);
      return (
        (tA + tB * cos_phi + tC * sin_phi) / (tD + tE * cos_phi + tF * sin_phi)
      );
    };

    const phi_t1 =
        -2 * Math.atan((-y12 + Math.sqrt(x12 ** 2 + y12 ** 2)) / x12),
      phi_t2 = 2 * Math.atan((y12 + Math.sqrt(x12 ** 2 + y12 ** 2)) / x12),
      t1 = getT(phi_t1),
      t2 = getT(phi_t2);

    //let [tMin12, phi12] = t1 < t2 ? [t1, phi_t1] : [t2, phi_t2];//Math.min(...[t1, t2]);//.filter(t => t > 0));
    for (const [t12, phi12] of [
      [t1, phi_t1],
      [t2, phi_t2],
    ]) {
      if (t12 < timeTillNextImpact) {
        let s12 = getS(phi12);
        if (s12 < 0) {
          const expr1 =
            r0 ** 2 * (vx ** 2 + vy ** 2) +
            2 * r0 * rr * (vx * x1c + vy * y1c) +
            rr ** 2 * (x1c ** 2 + y1c ** 2) -
            (vy * x1c - vx * y1c) ** 2;
          let t0_1 = 1 / 0,
            t0_2 = 1 / 0;
          if (expr1 >= 0) {
            const phiA1 = r0 * vx + rr * x1c,
              phiB1 = (rr - vx) * y1c + vy * (r0 + x1c);
            const phi_s0_1 = -2 * Math.atan((phiA1 - Math.sqrt(expr1)) / phiB1);
            const phi_s0_2 = -2 * Math.atan((phiA1 + Math.sqrt(expr1)) / phiB1);
            t0_1 = getT(phi_s0_1);
            t0_2 = getT(phi_s0_2);
            if (t0_1 < 0) {
              t0_1 = 1 / 0;
            }
            if (t0_2 < 0) {
              t0_2 = 1 / 0;
            }
            const [tCorner, phi_corner] =
              t0_1 < t0_2 ? [t0_1, phi_s0_1] : [t0_2, phi_s0_2];
            if (tCorner < timeTillNextImpact) {
              timeTillNextImpact = tCorner;
              phi4Impact = phi_corner;
              xComputedCollision = x1;
              yComputedCollision = y1;
              impactedSideIndex = -1;
              impactedSideNormal = null;
              impactedVertexIndex = side.index;
            }
          }
        } else if (s12 > 1) {
          let t1_1 = 1 / 0,
            t1_2 = 1 / 0;
          const expr2 =
            r0 ** 2 * (vx ** 2 + vy ** 2) +
            2 * r0 * rr * (vx * x2c + vy * y2c) +
            rr ** 2 * (x2c ** 2 + y2c ** 2) -
            (vy * x2c - vx * y2c) ** 2;
          if (expr2 >= 0) {
            const phiA2 = r0 * vx + rr * x2c,
              phiB2 = (rr - vx) * y2c + vy * (r0 + x2c);
            const phi_s1_1 = -2 * Math.atan((phiA2 - Math.sqrt(expr2)) / phiB2);
            const phi_s1_2 = -2 * Math.atan((phiA2 + Math.sqrt(expr2)) / phiB2);

            t1_1 = getT(phi_s1_1);
            t1_2 = getT(phi_s1_2);
            if (t1_1 < 0) {
              t1_1 = 1 / 0;
            }
            if (t1_2 < 0) {
              t1_2 = 1 / 0;
            }
            const [tCorner, phi_corner] =
              t1_1 < t1_2 ? [t1_1, phi_s1_1] : [t1_2, phi_s1_2];
            if (tCorner < timeTillNextImpact) {
              timeTillNextImpact = tCorner;
              phi4Impact = phi_corner;
              xComputedCollision = x2;
              yComputedCollision = y2;
              impactedSideIndex = -1;
              impactedSideNormal = null;
              impactedVertexIndex = side.next.index;
            }
          }
        } else {
          if (t12 > 0) {
            timeTillNextImpact = t12;
            phi4Impact = phi12;
            xComputedCollision = x1 + (x2 - x1) * s12;
            yComputedCollision = y1 + (y2 - y1) * s12;
            impactedSideIndex = side.index;
            impactedSideNormal = side.normal;
            impactedVertexIndex = -1;
          }
        }
      }
    }
  }

  const rAtImpact = r0 + rr * timeTillNextImpact;
  const computedBallAtNextImpact = {
    radius: rAtImpact,
    vRadius: previousCollisionBall.vRadius,
    x: xComputedCollision - rAtImpact * Math.cos(phi4Impact),
    y: yComputedCollision - rAtImpact * Math.sin(phi4Impact),
    impactedSideIndex,
    impactedVertexIndex,
  };

  let vx_next = 0,
    vy_next = 0;
  const ballEscaping = isBallEscaping(computedBallAtNextImpact);
  // if(ballEscaping){
  //  	console.warn('ballEscaping', collisionCounter, ballEscaping);
  // }

  ballMoving = !ballEscaping && timeTillNextImpact > 1e-8;
  let allImpacted;
  if (ballMoving) {
    allImpacted = computeAllImpacted(computedBallAtNextImpact, 1e-6);
    if (allImpacted.allowedDirections?.isVoid()) {
      ballMoving = false;
    }
  }
  if (ballMoving) {
    if (impactedVertexIndex === -1) {
      const normalizedNormal = impactedSideNormal;
      let v_normal = vx * normalizedNormal.x + vy * normalizedNormal.y;
      let reverseSign = v_normal > 0 ? -1 : 1;
      v_normal = reverseSign * v_normal;

      const delta_v_normal = 2 * v_normal - computedBallAtNextImpact.vRadius;
      vx_next = vx - /*reverseSign **/ delta_v_normal * normalizedNormal.x;
      vy_next = vy - /*reverseSign **/ delta_v_normal * normalizedNormal.y;
    } else {
      // collision with an angle
      //({vx_next, vy_next} = angleCollision_model1({vx, vy}, impactedVertexIndex));
      ({ vx_next, vy_next } = angleCollision_model2(
        computedBallAtNextImpact,
        { vx, vy },
        impactedVertexIndex
      ));
    }
    let magSet = false;
    if (allImpacted.allowedDirections) {
      const dir_v_next = Math.atan2(vy_next, vx_next);
      const insideOrClosest =
        allImpacted.allowedDirections.isDirectionInsideOrClosest(dir_v_next);
      let dir_v_next1 = dir_v_next;
      if (!insideOrClosest.inside) {
        //console.warn("not inside")
        const [dMin, dMax] = insideOrClosest.interval;
        const deltaDir = Math.min((5 * Math.PI) / 180, (dMax - dMin) / 10);
        dir_v_next1 = insideOrClosest.toMinInterval
          ? dMin + deltaDir
          : dMax - deltaDir;
        vx_next = vMag * Math.cos(dir_v_next1);
        vy_next = vMag * Math.sin(dir_v_next1);
        magSet = true;
      }
    }

    if (!magSet) {
      const vNextMag = Math.sqrt(vx_next ** 2 + vy_next ** 2);
      vx_next *= vMag / vNextMag;
      vy_next *= vMag / vNextMag;
    }
    computedBallAtNextImpact.moving = true;
  } else {
    computedBallAtNextImpact.escaping =
      ballEscaping || !verifyBall(computedBallAtNextImpact);
    computedBallAtNextImpact.moving = false;
  }
  computedBallAtNextImpact.vx = vx_next;
  computedBallAtNextImpact.vy = vy_next;
  computedBallAtNextImpact.t = previousCollisionBall.t + timeTillNextImpact;
  computedBallAtNextImpact.dt = timeTillNextImpact;

  if (!computedBallAtNextImpact.escaping) {
    collisionCounter++;
    previousCollisionBall = computedBallAtNextImpact;
  }

  return computedBallAtNextImpact;
}

function finalizeBall() {
  const ballFinally = { ...previousCollisionBall, dt: 0, final: true };
  return ballFinally;
  // let r = 1/0;
  // for(const side of polygonSides){
  // 	r = Math.min(r, _distanceToSide(ballFinally, side));
  // }
  // ballFinally.impactedVertexIndex = null;
  // ballFinally.impactedSideIndices = [];
  // ballFinally.radius = r;
  // return ballFinally;
}

// function angleCollision_model1({vx, vy}, impactedVertexIndex){
// 	// the bisecting line of the angle takes the place of the normal in the flat case
// 	// this is equivalent to taking the average of the normals of the two polygon lines meeting at that corner
//
// 	// computed as bisector; note that this can be precomputed for each polygon point
// 	const a = polygonVertices[impactedVertexIndex],
// 		b = polygonVertices[impactedVertexIndex === polygonVertices.length-2 ? 0 : impactedVertexIndex+1],
// 		c = polygonVertices[impactedVertexIndex === 0 ? polygonVertices.length- 2 : impactedVertexIndex-1],
// 		ba_normalized = _normalize({x: b.x - a.x, y: b.y - a.y}),
// 		ca_normalized = _normalize({x: c.x - a.x, y: c.y - a.y}),
// 		b1 = {x: a.x+ ba_normalized.x, y: a.y + ba_normalized.y},
// 		c1 = {x: a.x + ca_normalized.x, y: a.y + ca_normalized.y},
// 		bis_normalized = _normalize({x: b1.x + c1.x - 2 * a.x, y: b1.y + c1.y - 2 * a.y});
//
// 	// computed as the normalized average of the normalized normals
// 	// const normal1_normalized = _normalize({x: b.y - a.y, y: a.x - b.x}),
// 	// 	normal2_normalized = _normalize({x: a.y - c.y, y: c.x - a.x}),
// 	// 	normal_avg_normalized = _normalize({
// 	// 		x: normal1_normalized.x + normal2_normalized.x,
// 	// 		y: normal1_normalized.y + normal2_normalized.y
// 	// 	}); // dividing the sum by 2 (for average) isn't necessary because of the normalization
// 	// normal_avg_normalized is approx equal to bis_normalized or - bis_normalized
//
// 	const v_bis = vx * bis_normalized.x + vy * bis_normalized.y;
// 	const vx_next = vx - (2 * v_bis * bis_normalized.x),
// 		vy_next = vy - (2 * v_bis * bis_normalized.y);
//
// 	return {vx_next, vy_next};
// }

function angleCollision_model2(ball, { vx, vy }, impactedVertexIndex) {
  // the point of impact is taken as a very small sphere
  // then the "surface of contact" is the tangent to the ball at the point of impact
  // the normal from the flat case is now the radius of the ball at the point of impact

  const a = polygonVertices[impactedVertexIndex],
    normal_normalized = _normalize({ x: a.x - ball.x, y: a.y - ball.y });

  let v_normal = vx * normal_normalized.x + vy * normal_normalized.y;
  if (v_normal < 0) {
    v_normal = -v_normal;
  }
  const delta_v_normal = 2 * v_normal - ball.vRadius;
  const vx_next = vx - delta_v_normal * normal_normalized.x,
    vy_next = vy - delta_v_normal * normal_normalized.y;

  return { vx_next, vy_next };
}

function isBallEscaping(ball, delta = 1) {
  for (const { x, y, index } of polygonVertices) {
    const d = Math.sqrt((ball.x - x) ** 2 + (ball.y - y) ** 2);
    if (d < ball.radius - delta) {
      return 'V' + index;
    }
  }
  for (const side of polygonSides) {
    if (_distanceToSide(ball, side) < ball.radius - delta) {
      return 'S' + side.index;
    }
  }
  return false;
}

function verifyBall(ball) {
  return !(
    Math.sqrt((ball.x - size0 / 2) ** 2 + (ball.y - size0 / 2) ** 2) > RExt ||
    _distanceToSide(ball, {
      x1: 0,
      y1: 0,
      x2: size0,
      y2: 0,
      x12: -size0,
      y12: 0,
    }) < ball.radius ||
    _distanceToSide(ball, {
      x1: 0,
      y1: 0,
      x2: 0,
      y2: size0,
      x12: 0,
      y12: -size0,
    }) < ball.radius ||
    _distanceToSide(ball, {
      x1: size0,
      y1: size0,
      x2: size0,
      y2: 0,
      x12: 0,
      y12: size0,
    }) < ball.radius ||
    _distanceToSide(ball, {
      x1: size0,
      y1: size0,
      x2: 0,
      y2: size0,
      x12: size0,
      y12: 0,
    }) < ball.radius
  );
}
