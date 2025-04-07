importScripts("../utils/utils.js");

self.onmessage = function({data: {points, polygonVertices, minRadius, minCenter, subDivide, level}}){
	({minRadius, minCenter} = computeAllDistancesToPolygon(points, polygonVertices, minRadius, minCenter));
	self.postMessage({points, minRadius, minCenter});
	const newPoints = subDivide ? computeNextLevelGrid(points, minRadius, level) : [];
	self.postMessage({newPoints});
};

function computeAllDistancesToPolygon(points, polygonVertices, minRadius, minCenter){
	for(const point of points){
		const [x0, y0, r0] = point;
		let r;
		if(r0){
			r = r0;
		}
		else{
			r = 0;
			for(const vertex of polygonVertices){
				const ri = _distance2({x: x0, y: y0}, vertex);
				r = Math.max(r, ri);
			}
			point[2] = r;
		}
		if(r < minRadius){
			minCenter = [x0, y0];
			minRadius = r;
		}
	}
	return {minRadius, minCenter};
}

function computeNextLevelGrid(interiorPixels, maxRadius, level){
	const res = 1 / 3**level;
	let relativeRadii = interiorPixels.map(([_1, _2, r], i) => [maxRadius/r, i]);
	const intendedCount = 610 * 610 / 9;
	let relativeRadiusThreshold = 0.7;
	while(relativeRadii.length > intendedCount){
		relativeRadii = relativeRadii.filter(([rFrac])=>rFrac > relativeRadiusThreshold);
		relativeRadiusThreshold = relativeRadiusThreshold + (1-relativeRadiusThreshold)/5;
	}
	
	let topPixels = relativeRadii.map(([_, i]) => interiorPixels[i]);
	const n = topPixels.length;
	for(let i = 0; i < n; i++){
		const [x, y] = topPixels[i];
		const delta = res / 3;
		topPixels.push(
			[x - delta, y, null], [x + delta, y, null], [x, y - delta, null], [x, y + delta, null],
			[x - delta, y - delta, null], [x + delta, y - delta, null], [x - delta, y + delta, null],
			[x + delta, y + delta, null]
		);
	}
	
	return topPixels;
}

