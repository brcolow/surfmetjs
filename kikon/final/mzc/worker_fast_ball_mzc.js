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
let rExt0, rInt0;

self.onmessage = function ({ data }) {
  if (data.type === 'init') {
    ({ polygonVertices, previousCollisionBall, vMag, size0, RExt } = data);
    rExt0 = previousCollisionBall.rExt;
    rInt0 = previousCollisionBall.rInt;
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
      type: 'nextCollisionBalls',
      nextCollisionBalls: computeNextImpact(),
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

let deltaReverseExt = -1,
  deltaReverseInt = -1;
function computeNextImpact() {
  const ballInt = computeNextImpactInterior();
  const ballExt = computeNextImpactExterior();
  const nextBall =
    (ballInt?.dt ?? 1 / 0) < (ballExt?.dt ?? 1 / 0) ? ballInt : ballExt;
  let suppBall = null;
  nextBall.t = previousCollisionBall.t + nextBall.dt;
  if (!ballInt.moving) {
    if (deltaReverseInt < 0) {
      deltaReverseInt = (ballInt.rInt - rInt0) / 2; // initialization
    } else {
      deltaReverseInt = deltaReverseInt / 2;
    }
    if (deltaReverseInt < 1) {
      nextBall.moving = false;
    } else {
      nextBall.moving = true;
      suppBall = { ...nextBall };
      suppBall.rInt -= deltaReverseInt;
      suppBall.dt = deltaReverseInt / suppBall.vrInt / 5;
      suppBall.t += suppBall.dt;
      suppBall.vx = -previousCollisionBall.vx;
      suppBall.vy = -previousCollisionBall.vy;
      suppBall.moving = true;
    }
  }
  if (!ballExt.moving) {
    if (deltaReverseExt < 0) {
      deltaReverseExt = (rExt0 - ballExt.rExt) / 2; // initialization
    } else {
      deltaReverseExt = deltaReverseExt / 2;
    }
    if (deltaReverseExt < 1) {
      nextBall.moving = false;
    } else {
      nextBall.moving = true;
      suppBall = { ...nextBall };
      suppBall.rExt += deltaReverseExt;
      suppBall.dt = deltaReverseExt / suppBall.vrExt / 5;
      suppBall.t += suppBall.dt;
      suppBall.vx = -previousCollisionBall.vx;
      suppBall.vy = -previousCollisionBall.vy;
      suppBall.moving = true;
    }
  }
  previousCollisionBall = suppBall ?? nextBall;
  ballMoving = nextBall.moving;
  return suppBall ? [nextBall, suppBall] : [nextBall];
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Interior collision (from MIC)

function computeAllImpactedInterior(
  { x: xcNext, y: ycNext, rInt: rNext },
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

function computeNextImpactInterior() {
  let xComputedCollision = null,
    yComputedCollision = null;
  let timeTillNextImpact = 1 / 0;
  let phi4Impact = -1;
  let impactedSideIndex = -1;
  let impactedSideNormal = null;
  let impactedVertexIndex = -1;
  const r0 = previousCollisionBall.rInt;
  const rr = previousCollisionBall.vrInt;
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

  const rIntAtImpact = r0 + rr * timeTillNextImpact;
  const computedBallAtNextImpact = {
    rInt: rIntAtImpact,
    rExt:
      previousCollisionBall.rExt -
      previousCollisionBall.vrExt * timeTillNextImpact,
    vrInt: previousCollisionBall.vrInt,
    vrExt: previousCollisionBall.vrExt,
    x: xComputedCollision - rIntAtImpact * Math.cos(phi4Impact),
    y: yComputedCollision - rIntAtImpact * Math.sin(phi4Impact),
    impactedSideIndex,
    impactedVertexIndex,
  };

  let vx_next = 0,
    vy_next = 0;
  const ballEscaping = isBallEscaping(computedBallAtNextImpact);
  // if(ballEscaping){
  //  	console.warn('ballEscaping', collisionCounter, ballEscaping);
  // }

  let ballIntMoving = !ballEscaping && timeTillNextImpact > 1e-8;
  let allImpacted;
  if (ballIntMoving) {
    allImpacted = computeAllImpactedInterior(computedBallAtNextImpact, 1e-6);
    if (allImpacted.allowedDirections?.isVoid()) {
      ballIntMoving = false;
    }
  }
  if (ballIntMoving) {
    if (impactedVertexIndex === -1) {
      const normalizedNormal = impactedSideNormal;
      let v_normal = vx * normalizedNormal.x + vy * normalizedNormal.y;
      let reverseSign = v_normal > 0 ? -1 : 1;
      v_normal = reverseSign * v_normal;

      const delta_v_normal = 2 * v_normal - computedBallAtNextImpact.vrInt;
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
    computedBallAtNextImpact.invalid =
      ballEscaping || !verifyBallInt(computedBallAtNextImpact);
    computedBallAtNextImpact.moving = false;
  }
  computedBallAtNextImpact.vx = vx_next;
  computedBallAtNextImpact.vy = vy_next;
  computedBallAtNextImpact.t = previousCollisionBall.t + timeTillNextImpact;
  computedBallAtNextImpact.dt = timeTillNextImpact;

  //if(!computedBallAtNextImpact.invalid){
  //collisionCounter++;
  //previousCollisionBall = computedBallAtNextImpact;
  //}

  return computedBallAtNextImpact;
}

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
  const delta_v_normal = 2 * v_normal - ball.vrInt;
  const vx_next = vx - delta_v_normal * normal_normalized.x,
    vy_next = vy - delta_v_normal * normal_normalized.y;

  return { vx_next, vy_next };
}

function isBallEscaping(ball, delta = 1) {
  for (const { x, y, index } of polygonVertices) {
    const d = Math.sqrt((ball.x - x) ** 2 + (ball.y - y) ** 2);
    if (d < ball.rInt - delta) {
      return 'V' + index;
    }
  }
  for (const side of polygonSides) {
    if (_distanceToSide(ball, side) < ball.rInt - delta) {
      return 'S' + side.index;
    }
  }
  return false;
}

function verifyBallInt(ball) {
  return !(
    Math.sqrt((ball.x - size0 / 2) ** 2 + (ball.y - size0 / 2) ** 2) > RExt ||
    _distanceToSide(ball, {
      x1: 0,
      y1: 0,
      x2: size0,
      y2: 0,
      x12: -size0,
      y12: 0,
    }) < ball.rInt ||
    _distanceToSide(ball, {
      x1: 0,
      y1: 0,
      x2: 0,
      y2: size0,
      x12: 0,
      y12: -size0,
    }) < ball.rInt ||
    _distanceToSide(ball, {
      x1: size0,
      y1: size0,
      x2: size0,
      y2: 0,
      x12: 0,
      y12: size0,
    }) < ball.rInt ||
    _distanceToSide(ball, {
      x1: size0,
      y1: size0,
      x2: 0,
      y2: size0,
      x12: size0,
      y12: 0,
    }) < ball.rInt
  );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Exterior collision (from MCC)

function computeAllImpactedExterior(xcNext, ycNext, rNext, threshold) {
  const impactedVertexIndices = [];
  for (const { x, y, index } of polygonVertices) {
    const ri = Math.sqrt((x - xcNext) ** 2 + (y - ycNext) ** 2);
    if (Math.abs(ri - rNext) < threshold) {
      impactedVertexIndices.push(index);
    }
  }
  return impactedVertexIndices;
}

function computeNextImpactExterior() {
  const r0 = previousCollisionBall.rExt;
  const vr = previousCollisionBall.vrExt;
  const vx = previousCollisionBall.vx;
  const vy = previousCollisionBall.vy;
  const xc = previousCollisionBall.x;
  const yc = previousCollisionBall.y;

  const theta_v0 = Math.acos(vr / vMag);

  let tFirstImpact = 1 / 0;
  for (const { x, y, index } of polygonVertices) {
    //if(!previousCollisionBall.impactedVertexIndices.includes(index)){
    const a = vx ** 2 + vy ** 2 - vr ** 2,
      b = vr * r0 - (x - xc) * vx - (y - yc) * vy,
      c = (x - xc) ** 2 + (y - yc) ** 2 - r0 ** 2,
      d = b ** 2 - a * c;
    if (d >= 0) {
      const dSqrt = Math.sqrt(d),
        t1 = (-b + dSqrt) / a,
        t2 = (-b - dSqrt) / a;
      for (const tContact of [t1, t2]) {
        if (tContact >= 0 && tContact < tFirstImpact && tContact >= 1e-6) {
          tFirstImpact = tContact;
        }
      }
    }
    //}
  }

  const xcNext = xc + vx * tFirstImpact,
    ycNext = yc + vy * tFirstImpact,
    rNext = r0 - vr * tFirstImpact;

  // get the exact number of vertices that are touching the circle
  const impactedVertexIndices = computeAllImpactedExterior(
    xcNext,
    ycNext,
    rNext,
    1e-5
  );

  let computedBallAtNextImpact = {
    rExt: rNext,
    rInt:
      previousCollisionBall.rInt + previousCollisionBall.vrInt * tFirstImpact,
    vrExt: vr,
    vrInt: previousCollisionBall.vrInt,
    x: xcNext,
    y: ycNext,
    impactedVertexIndices,
  };

  const vNext = { x: 0, y: 0 };

  const ballEnteringPolygon = isBallEnteringPolygon(computedBallAtNextImpact);
  if (ballEnteringPolygon) {
    computedBallAtNextImpact = { ...previousCollisionBall, dt: 0 };
    if (isBallEnteringPolygon(computedBallAtNextImpact)) {
      computedBallAtNextImpact.invalid = true; // typically when the initial solution is inside the polygon
    }
  }
  let stopBall = ballEnteringPolygon || tFirstImpact <= 1e-5;
  if (!ballEnteringPolygon) {
    let impactNormal = { x: 0, y: 0 };
    for (const ii of impactedVertexIndices) {
      const { x: xi, y: yi } = polygonVertices[ii];
      _addToXY(
        impactNormal,
        /*_normalized*/ { x: xcNext - xi, y: ycNext - yi }
      );
    }

    const normalizedNormal = _normalized(impactNormal);
    const v_normal = vx * normalizedNormal.x + vy * normalizedNormal.y;
    const delta_v_normal = 2 * v_normal;
    vNext.x =
      vx - delta_v_normal * normalizedNormal.x - 0.1 * vr * normalizedNormal.x;
    vNext.y =
      vy - delta_v_normal * normalizedNormal.y - 0.1 * vr * normalizedNormal.y;
    const vx1 = vNext.x,
      vy1 = vNext.y;
    _setNormTo(vNext, vMag);
    let theta_v = Math.atan2(vNext.y, vNext.x);

    let allowedDirections = CircularZones.fullCircle();
    for (const vertexIndex of impactedVertexIndices) {
      const vertex = polygonVertices[vertexIndex];
      const dirOutCenter = Math.atan2(vertex.y - ycNext, vertex.x - xcNext);
      allowedDirections.intersectWith(
        CircularZones.directionPlusMinus(dirOutCenter, theta_v0)
      );
    }

    stopBall = allowedDirections.isVoid();
    if (!stopBall) {
      const ret = allowedDirections.isDirectionInsideOrClosest(theta_v);
      if (!ret.inside) {
        const theta_v = (ret.interval[0] + ret.interval[1]) / 2; // sefer solution
        _rotateTo(vNext, theta_v);
      }
    }
  }

  computedBallAtNextImpact.vPrev = { vx: vx, vy: vy }; // for illustration purposes in draw
  computedBallAtNextImpact.vx = stopBall ? 0 : vNext.x;
  computedBallAtNextImpact.vy = stopBall ? 0 : vNext.y;
  computedBallAtNextImpact.moving = !stopBall;

  computedBallAtNextImpact.dt = tFirstImpact;

  collisionCounter++;

  return computedBallAtNextImpact;
}

function isBallEnteringPolygon(ball, delta = 0.1) {
  for (const { x, y, index } of polygonVertices) {
    const ri = Math.sqrt((x - ball.x) ** 2 + (y - ball.y) ** 2);
    if (ball.rExt < ri - delta) {
      return 'V' + index;
    }
  }
  return false;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

function finalizeBall() {
  const ballFinally = { ...previousCollisionBall, dt: 0, final: true };
  return ballFinally;
}
