importScripts("../utils/CircularZones.js", "../utils/utils.js");

let polygonVertices, polygonSides, previousCollisionBall, vMag, collisionCounter, ballMoving;

self.onmessage = function({data}){
	if(data.type === 'init'){
		({polygonVertices, ball: previousCollisionBall, vMag} = data);
		polygonSides = buildPolygonData(polygonVertices);
		collisionCounter = 0;
		ballMoving = true;
		self.postMessage({type: 'postInit', polygonVertices, polygonSides, ball: previousCollisionBall});
	}
	else if(data.type === 'computeNextCollision'){
		if(!ballMoving){
			return; // don't answer
		}
		const nextCollisionBall = computeNextCollision();
		self.postMessage({type: 'nextCollisionBall', nextCollisionBall});
	}
	else if(data.type === 'finalizeBall'){
		if(ballMoving){
			return;
		}
		const ballFinally = finalizeBall();
		self.postMessage({type: 'ballFinally', ballFinally});
	}
};

function buildPolygonData(polygonVertices){
	// precomputes all static (that don't depend on the ball) properties of the polygon vertices and sides
	const nVertices = polygonVertices.length;
	for(let i = 1; i <= nVertices; i++){
		const current = polygonVertices[i - 1],
			next = polygonVertices[i === nVertices ? 0 : i],
			prev = polygonVertices[i === 1 ? nVertices - 1 : i - 2];
		Object.assign(current, {
			next,
			prev,
			index: i - 1,
			unitDirections: {
				prev: _unitDir({from: current, to: prev}),
				next: _unitDir({from: current, to: next}),
			}
		});
	}
	
	const sides = [];
	for(let i = 1; i <= nVertices; i++){
		const p1 = polygonVertices[i-1],
			p2 = p1.next
		sides.push({
			index: i-1,
			p1, p2
		});
	}
	for(let i = 1; i <= nVertices; i++){
		const side = sides[i - 1],
			sideNext = sides[i === nVertices ? 0 : i];
		side.next = sideNext;
		sideNext.previous = side;
	}
	return sides;
}

function computeImpactedVertices(xcNext, ycNext, rNext, threshold){
	const impactedVertexIndices = [];
	for(const {x, y, index} of polygonVertices){
		const ri = Math.sqrt((x - xcNext) ** 2 + (y - ycNext) ** 2);
		if(Math.abs(ri - rNext) < threshold){
			impactedVertexIndices.push(index);
		}
	}
	return impactedVertexIndices;
}

function computeNextCollision(){
	const r0 = previousCollisionBall.radius;
	const vr = previousCollisionBall.vRadius;
	const vx = previousCollisionBall.vx;
	const vy = previousCollisionBall.vy;
	const xc = previousCollisionBall.x;
	const yc = previousCollisionBall.y;
	
	const theta_v0 = Math.acos(vr / vMag);
	
	let tFirstImpact = 1/0;
	for(const {x, y, index} of polygonVertices){
		//if(!previousCollisionBall.impactedVertexIndices.includes(index)){
			const a = vx ** 2 + vy ** 2 - vr ** 2,
				b = (vr * r0 - (x - xc) * vx - (y - yc) * vy),
				c = (x - xc) ** 2 + (y - yc) ** 2 - r0 ** 2,
				d = b ** 2 - a * c;
			if(d >= 0){
				const dSqrt = Math.sqrt(d),
					t1 = (-b + dSqrt) / a,
					t2 = (-b - dSqrt) / a;
				for(const tContact of [t1, t2]){
					if(tContact >= 0 && tContact < tFirstImpact && tContact >= 1e-6){
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
	const impactedVertexIndices = computeImpactedVertices(xcNext, ycNext, rNext, 1e-5);
	
	let computedBallAtNextImpact = {
		radius: rNext,
		vRadius: vr,
		x: xcNext,
		y: ycNext,
		impactedVertexIndices
	}
	
	const vNext = {x: 0, y: 0};
	
	const ballEnteringPolygon = isBallEnteringPolygon(computedBallAtNextImpact);
	if(ballEnteringPolygon){
		computedBallAtNextImpact = {...previousCollisionBall, dt: 0};
		if(isBallEnteringPolygon(computedBallAtNextImpact)){
			computedBallAtNextImpact.invalid = true; // typically when the initial solution is inside the polygon
		}
	}
	let stopBall = ballEnteringPolygon || tFirstImpact <= 1e-5;
	if(!ballEnteringPolygon){
		let impactNormal = {x: 0, y: 0};
		for(const ii of impactedVertexIndices){
			const {x: xi, y: yi} = polygonVertices[ii];
			_addToXY(impactNormal, /*_normalized*/({x: xcNext - xi, y: ycNext - yi}));
		}
		
		const normalizedNormal = _normalized(impactNormal);
		const v_normal = vx * normalizedNormal.x + vy * normalizedNormal.y;
		const delta_v_normal = 2 * v_normal;
		vNext.x = vx - (delta_v_normal * normalizedNormal.x) - 0.1 * vr * normalizedNormal.x;
		vNext.y = vy - (delta_v_normal * normalizedNormal.y) - 0.1 * vr * normalizedNormal.y;
		const vx1 = vNext.x, vy1 = vNext.y;
		_setNormTo(vNext, vMag);
		let theta_v = Math.atan2(vNext.y, vNext.x);
		
		let allowedDirections = CircularZones.fullCircle();
		for(const vertexIndex of impactedVertexIndices){
			const vertex = polygonVertices[vertexIndex];
			const dirOutCenter = Math.atan2(vertex.y - ycNext, vertex.x - xcNext);
			allowedDirections.intersectWith(CircularZones.directionPlusMinus(dirOutCenter, theta_v0));
		}
		
		stopBall = allowedDirections.isVoid();
		if(!stopBall){
			const ret = allowedDirections.isDirectionInsideOrClosest(theta_v);
			if(!ret.inside){
				const theta_v = (ret.interval[0] + ret.interval[1]) / 2; // sefer solution
				_rotateTo(vNext, theta_v);
			}
		}
	}
	
	computedBallAtNextImpact.vPrev = {vx: vx, vy: vy}; // for illustration purposes in draw
	computedBallAtNextImpact.vx = stopBall ? 0 : vNext.x;
	computedBallAtNextImpact.vy = stopBall ? 0 : vNext.y;
	computedBallAtNextImpact.running = !stopBall;
	ballMoving = !stopBall;
	
	computedBallAtNextImpact.t = previousCollisionBall.t + tFirstImpact;
	computedBallAtNextImpact.dt = tFirstImpact;
	
	collisionCounter++;
	
	previousCollisionBall = computedBallAtNextImpact;
	return computedBallAtNextImpact;
}

function isBallEnteringPolygon(ball, delta = 0.1){
	for(const {x, y, index} of polygonVertices){
		const ri = Math.sqrt((x - ball.x)**2 + (y - ball.y)**2);
		if(ball.radius < ri - delta){
			return 'V'+index;
		}
	}
	return false;
}

function finalizeBall(){
	const ballFinally = {...previousCollisionBall, dt: 0, final: true};
	return ballFinally;
	
	// const r0 = previousCollisionBall.radius;
	// const xc0 = previousCollisionBall.x;
	// const yc0 = previousCollisionBall.y;
	//
	// const closestContactPoints = Array(3).fill(-1),
	// 	minDistancesToBall = Array(3).fill(1/0);
	// for(const {x, y, index} of polygonVertices){
	// 	const dist = Math.abs(Math.sqrt((xc0 - x) ** 2 + (yc0 - y) ** 2) - r0);
	// 	for(let i = 0; i < minDistancesToBall.length; i++){
	// 		if(dist < minDistancesToBall[i]){
	// 			for(let j = minDistancesToBall.length - 1; j >= i + 1; j--){
	// 				minDistancesToBall[j] = minDistancesToBall[j - 1];
	// 				closestContactPoints[j] = closestContactPoints[j - 1];
	// 			}
	// 			minDistancesToBall[i] = dist;
	// 			closestContactPoints[i] = index;
	// 			break;
	// 		}
	// 	}
	// }
	//
	// const {x: x1, y: y1} = polygonVertices[closestContactPoints[0]],
	// 	{x: x2, y: y2} = polygonVertices[closestContactPoints[1]],
	// 	{x: x3, y: y3} = polygonVertices[closestContactPoints[2]];
	// //using https://stackoverflow.com/a/78328096/16466946
	// const d = 2 * (x1 * (y2-y3) + x2 * (y3-y1) + x3 * (y1-y2));
	// if (Math.abs(d) < 1e-6){
	// 	return ballFinally;
	// }
	// const [t1, t2, t3] = [x1*x1+y1*y1, x2*x2+y2*y2, x3*x3+y3*y3],
	// 	xc1 = (t1 * (y2 - y3) + t2 * (y3 - y1) + t3 * (y1 - y2)) / d,
	// 	yc1 = (t1 * (x3 - x2) + t2 * (x1 - x3) + t3 * (x2 - x1)) / d,
	// 	r1 = Math.sqrt((x1-xc1)**2 + (y1-yc1)**2);
	//
	// ballFinally.x = xc1;
	// ballFinally.y = yc1;
	// ballFinally.radius = r1;
	// ballFinally.dt = Math.sqrt((xc0 - xc1)**2 + (yc0 - yc1)**2)/vMag;
	// ballFinally.t += ballFinally.dt;
	// ballFinally.impactedVertexIndices = computeImpactedVertices(xc1, yc1, r1, 0.1);
	//
	// return ballFinally;
}

