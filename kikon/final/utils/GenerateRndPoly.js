if (!globalThis.document) {
  // if used as worker
  self.onmessage = function ({ data }) {
    if (data.type === 'generate') {
      const { nInitialVertices, RExt, center, insertAtIndices } = data;
      const polygonVertices = generateRndPolygonVertices(
        nInitialVertices,
        RExt,
        center,
        insertAtIndices
      );
      self.postMessage({ type: 'generated', polygonVertices });
    }
  };
}

function generateRndPolygonVertices(
  nInitialVertices,
  RExt,
  center,
  { atIndices: { insertAt, moveAt } } = {
    atIndices: { insertAt: [], moveAt: [] },
  },
  { radial: { insertR, moveR } } = {
    radial: {
      insertR: { from: 0.5, to: 1.25 },
      moveR: { from: 0.5, to: 1.25 },
    },
  },
  { angular: { insertAng, moveAng } } = {
    angular: {
      insertAng: { from: -0.4, to: 0.4 },
      moveAng: { from: -0.4, to: 0.4 },
    },
  }
) {
  if (nInitialVertices.from && nInitialVertices.to) {
    nInitialVertices = Math.floor(
      Math.random() * (nInitialVertices.to - nInitialVertices.from + 1) +
        nInitialVertices.from
    );
  }
  if (RExt.from && RExt.to) {
    RExt = Math.random() * (RExt.to - RExt.from) + RExt.from;
  }
  if (center.x.from && center.x.to) {
    center.x = Math.random() * (center.x.to - center.x.from) + center.x.from;
  }
  if (center.y.from && center.y.to) {
    center.y = Math.random() * (center.y.to - center.y.from) + center.y.from;
  }
  const polygonVertices = regularPolygonVertices(
    nInitialVertices,
    RExt,
    center
  );

  if (moveAt && moveAt.length > 0) {
    if (Number.isFinite(moveR)) {
      moveR = { from: moveR, to: moveR };
    }
    if (Number.isFinite(moveAng)) {
      moveAng = { from: moveAng, to: moveAng };
    }
    const dAngle = (2 * Math.PI) / nInitialVertices;
    for (const move of moveAt) {
      let index =
        move.from && move.to
          ? Math.floor(Math.random() * (move.to - move.from + 1) + move.from)
          : move;
      if (Number.isInteger(index)) {
        if (index < 0) {
          index = nInitialVertices - index;
        }
        const r = (Math.random() * (moveR.to - moveR.from) + moveR.from) * RExt,
          angleFrom = dAngle * index + moveAng.from * dAngle,
          angleTo = dAngle * index + moveAng.to * dAngle,
          angle = Math.random() * (angleTo - angleFrom) + angleFrom;
        polygonVertices[index].x = center.x + r * Math.cos(angle);
        polygonVertices[index].y = center.y + r * Math.sin(angle);
      }
    }
  }

  if (insertAt && insertAt.length > 0) {
    const insertedVertices = [];
    if (Number.isFinite(insertR)) {
      insertR = { from: insertR, to: insertR };
    }
    if (Number.isFinite(insertAng)) {
      insertAng = { from: insertAng, to: insertAng };
    }
    const dAngle = (2 * Math.PI) / nInitialVertices;
    for (const insert of insertAt) {
      let index =
        insert.from && insert.to
          ? Math.floor(
              Math.random() * (insert.to - insert.from + 1) + insert.from
            )
          : insert;
      if (Number.isInteger(index)) {
        if (index < 0) {
          index = nInitialVertices + index;
        }
        const r =
            (Math.random() * (insertR.to - insertR.from) + insertR.from) * RExt,
          angleFrom =
            dAngle * (index + 0 * insertedVertices.length) +
            insertAng.from * dAngle,
          angleTo =
            dAngle * (index + 0 * insertedVertices.length) +
            insertAng.to * dAngle,
          angle = Math.random() * (angleTo - angleFrom) + angleFrom,
          xi = center.x + r * Math.cos(angle),
          yi = center.y + r * Math.sin(angle);
        insertedVertices.unshift({ index, x: xi, y: yi });
      }
    }

    insertedVertices.forEach(({ index, x, y }) => {
      polygonVertices.splice(index, 0, { x, y });
    });
  }

  return polygonVertices;
}

function regularPolygonVertices(sides, radius, center = { x: 0, y: 0 }) {
  const vertices = [];
  const angle = 360 / sides;

  for (let i = 0; i < sides; i++) {
    const x = radius * Math.cos((Math.PI * i * angle) / 180) + center.x;
    const y = radius * Math.sin((Math.PI * i * angle) / 180) + center.y;
    vertices.push({ x, y });
  }

  return vertices;
}
