#!/usr/bin/env node

const fs = require('fs');
const path = require('path');
const { execFileSync } = require('child_process');

const rootDir = path.resolve(__dirname, '..');
const cgalDir = path.join(rootDir, 'cgal');
const cgalSource = path.join(cgalDir, 'cgal.cpp');
const cgalExe = path.join(cgalDir, process.platform === 'win32' ? 'cgal' : 'cgal');
const fixtureDir = path.join(rootDir, 'test', 'fixtures', 'roundness');
const wktDir = path.join(fixtureDir, 'wkt');
const fixturePath = path.join(fixtureDir, 'roundness_profiles.json');

main();

function main() {
  fs.mkdirSync(wktDir, { recursive: true });
  ensureCgalBuilt();

  const profiles = buildProfiles();
  const fixtures = profiles.map(profile => {
    const wktPath = path.join(wktDir, `${profile.name}.wkt`);
    fs.writeFileSync(wktPath, toWkt(profile.points));

    const cgalOutput = runCgal(wktPath);
    const cgal = parseCgalOutput(cgalOutput);

    assertCgalCerts(profile.name, cgal);

    return {
      name: profile.name,
      description: profile.description,
      source: profile.source,
      points: profile.points.map(([x, y]) => [roundFixtureNumber(x), roundFixtureNumber(y)]),
      cgal
    };
  });

  const payload = {
    schemaVersion: 1,
    generatedBy: 'scripts/generate_roundness_fixtures.cjs',
    cgalSource: 'cgal/cgal.cpp',
    profileReconstruction: 'CHORDAL_POLYGON_FROM_ORDERED_SAMPLES',
    note: 'CGAL MZC/MCC goldens are exact for extracted finite points; MIC is exact for the chordal polygon reconstructed from the ordered samples.',
    profiles: fixtures
  };

  fs.writeFileSync(fixturePath, `${JSON.stringify(payload, null, 2)}\n`);
  console.log(`Wrote ${fixtures.length} roundness fixtures to ${path.relative(rootDir, fixturePath)}`);
}

function ensureCgalBuilt() {
  if (!fs.existsSync(cgalSource)) {
    throw new Error(`Missing ${cgalSource}`);
  }

  const shouldBuild =
    !fs.existsSync(cgalExe) ||
    fs.statSync(cgalSource).mtimeMs > fs.statSync(cgalExe).mtimeMs;

  if (!shouldBuild) {
    return;
  }

  console.log('Building cgal/cgal...');
  runInCgalDir(
    'g++ -std=c++17 -O3 -g -fno-omit-frame-pointer cgal.cpp ' +
    '-I./build/include -L./build/lib -lmpfr -lgmp -o cgal'
  );
}

function runCgal(wktPath) {
  return runInCgalDir(`./cgal ${shellQuote(toExecutionPath(wktPath))}`);
}

function runInCgalDir(command) {
  if (process.platform === 'win32') {
    return execFileSync(
      'bash',
      ['-lc', `cd ${shellQuote(toExecutionPath(cgalDir))} && ${command}`],
      { encoding: 'utf8', stdio: ['ignore', 'pipe', 'pipe'] }
    );
  }

  return execFileSync(
    'sh',
    ['-lc', `cd ${shellQuote(cgalDir)} && ${command}`],
    { encoding: 'utf8', stdio: ['ignore', 'pipe', 'pipe'] }
  );
}

function toExecutionPath(filePath) {
  if (process.platform !== 'win32') {
    return filePath;
  }

  const resolved = path.resolve(filePath);
  const drive = resolved[0].toLowerCase();
  const rest = resolved.slice(2).replace(/\\/g, '/');
  return `/mnt/${drive}${rest}`;
}

function shellQuote(value) {
  return `'${String(value).replace(/'/g, `'\\''`)}'`;
}

function buildProfiles() {
  const mixed = radialProfile({
    name: 'mixed_harmonics_64',
    description: 'Centered polar profile with several low-order harmonics.',
    count: 64,
    baseRadius: 10,
    harmonics: [
      { order: 3, amplitude: 0.055, phase: 0.35 },
      { order: 5, amplitude: -0.028, phase: 1.1 },
      { order: 8, amplitude: 0.012, phase: -0.4 }
    ]
  });

  return [
    circleProfile({
      name: 'perfect_circle_32',
      description: 'All sampled points lie on one circle.',
      count: 32,
      radius: 10,
      center: [0, 0]
    }),
    radialProfile({
      name: 'four_lobe_48',
      description: 'Centered four-lobe profile with a small second harmonic.',
      count: 48,
      baseRadius: 10,
      harmonics: [
        { order: 4, amplitude: 0.06, phase: 0.2 },
        { order: 8, amplitude: 0.015, phase: -0.6 }
      ]
    }),
    mixed,
    radialProfile({
      name: 'localized_dent_64',
      description: 'Centered profile with one broad bump and one local dent.',
      count: 64,
      baseRadius: 10,
      harmonics: [
        { order: 2, amplitude: 0.025, phase: 0.3 },
        { order: 6, amplitude: 0.018, phase: -0.8 }
      ],
      lobes: [
        { center: 1.15, amplitude: -0.07, sigma: 0.16 },
        { center: 4.2, amplitude: 0.045, sigma: 0.24 }
      ]
    }),
    transformProfile(mixed, {
      name: 'mixed_harmonics_translated_scaled',
      description: 'Scaled, rotated, and translated copy used for metamorphic checks.',
      scale: 1.7,
      rotation: 0.73,
      translation: [25, -13]
    })
  ];
}

function circleProfile({ name, description, count, radius, center }) {
  return {
    name,
    description,
    source: { type: 'circle', count, radius, center },
    points: Array.from({ length: count }, (_, i) => {
      const theta = 2 * Math.PI * i / count;
      return [
        center[0] + radius * Math.cos(theta),
        center[1] + radius * Math.sin(theta)
      ];
    })
  };
}

function radialProfile({ name, description, count, baseRadius, harmonics = [], lobes = [] }) {
  return {
    name,
    description,
    source: { type: 'polar', count, baseRadius, harmonics, lobes },
    points: Array.from({ length: count }, (_, i) => {
      const theta = 2 * Math.PI * i / count;
      const radius = radiusAt(theta, baseRadius, harmonics, lobes);
      return [
        radius * Math.cos(theta),
        radius * Math.sin(theta)
      ];
    })
  };
}

function transformProfile(profile, { name, description, scale, rotation, translation }) {
  const cos = Math.cos(rotation);
  const sin = Math.sin(rotation);

  return {
    name,
    description,
    source: {
      type: 'transform',
      baseProfile: profile.name,
      scale,
      rotation,
      translation
    },
    points: profile.points.map(([x, y]) => [
      translation[0] + scale * (x * cos - y * sin),
      translation[1] + scale * (x * sin + y * cos)
    ])
  };
}

function radiusAt(theta, baseRadius, harmonics, lobes) {
  let radius = baseRadius;

  for (const { order, amplitude, phase } of harmonics) {
    radius += amplitude * Math.cos(order * theta + phase);
  }

  for (const { center, amplitude, sigma } of lobes) {
    const delta = angularDistance(theta, center);
    radius += amplitude * Math.exp(-(delta * delta) / (2 * sigma * sigma));
  }

  if (radius <= 0) {
    throw new Error(`Generated non-positive radius ${radius}`);
  }

  return radius;
}

function angularDistance(a, b) {
  let delta = a - b;
  while (delta <= -Math.PI) delta += 2 * Math.PI;
  while (delta > Math.PI) delta -= 2 * Math.PI;
  return delta;
}

function toWkt(points) {
  const closed = [...points, points[0]];
  return `POLYGON((${closed.map(([x, y]) => `${formatWktNumber(x)} ${formatWktNumber(y)}`).join(', ')}))\n`;
}

function formatWktNumber(value) {
  return value.toFixed(12).replace(/\.?0+$/, '');
}

function roundFixtureNumber(value) {
  return Number(value.toFixed(12));
}

function parseCgalOutput(output) {
  const values = {};
  for (const line of output.split(/\r?\n/)) {
    const match = /^([A-Z0-9_]+)=(.*)$/.exec(line.trim());
    if (match) {
      values[match[1]] = match[2];
    }
  }

  return {
    input: {
      vertexCount: integer(values.INPUT_VERTEX_COUNT),
      uniquePointCount: integer(values.INPUT_UNIQUE_POINT_COUNT),
      simple: values.INPUT_SIMPLE === '1',
      pointSetSource: values.INPUT_POINT_SET_SOURCE
    },
    mic: circle(values, 'MIC'),
    mcc: circle(values, 'MCC'),
    mzc: {
      valid: values.MZC_VALID === '1',
      domain: values.MZC_DOMAIN,
      objective: values.MZC_OBJECTIVE,
      solverExactFor: values.MZC_SOLVER_EXACT_FOR,
      center: point(values.MZC_CENTER_DEC),
      innerRadius: number(values.MZC_INNER_RADIUS_DEC),
      outerRadius: number(values.MZC_OUTER_RADIUS_DEC),
      width: number(values.MZC_WIDTH_DEC),
      centerExact: [values.MZC_CENTER_EXACT_X, values.MZC_CENTER_EXACT_Y],
      innerRadiusSquaredExact: values.MZC_R_IN2_EXACT,
      outerRadiusSquaredExact: values.MZC_R_OUT2_EXACT,
      widthExact: values.MZC_WIDTH_EXACT,
      candidatesEvaluated: integer(values.MZC_CANDIDATES_EVALUATED)
    },
    polylineDiagnostic: values.MZC_POLYLINE_DIAGNOSTIC ? {
      diagnostic: values.MZC_POLYLINE_DIAGNOSTIC,
      innerRadius: number(values.MZC_POLYLINE_INNER_RADIUS_DEC),
      outerRadius: number(values.MZC_POLYLINE_OUTER_RADIUS_DEC),
      width: number(values.MZC_POLYLINE_WIDTH_DEC),
      innerRadiusSquaredExact: values.MZC_POLYLINE_R_IN2_EXACT,
      outerRadiusSquaredExact: values.MZC_POLYLINE_R_OUT2_EXACT,
      widthExact: values.MZC_POLYLINE_WIDTH_EXACT
    } : null,
    certs: {
      mic: values.MIC_CERT,
      mcc: values.MCC_CERT,
      mzc: values.MZC_CERT
    }
  };
}

function circle(values, label) {
  return {
    valid: values[`${label}_VALID`] === '1',
    domain: values[`${label}_DOMAIN`],
    center: point(values[`${label}_CENTER_DEC`]),
    radius: number(values[`${label}_RADIUS_DEC`]),
    centerExact: [
      values[`${label}_CENTER_EXACT_X`],
      values[`${label}_CENTER_EXACT_Y`]
    ],
    radiusSquaredExact: values[`${label}_RADIUS2_EXACT`],
    radiusExact: values[`${label}_RADIUS_EXACT`]
  };
}

function assertCgalCerts(name, cgal) {
  for (const label of ['mic', 'mcc', 'mzc']) {
    if (!cgal[label].valid) {
      throw new Error(`${name}: CGAL ${label.toUpperCase()} result was invalid`);
    }
  }

  for (const [label, status] of Object.entries(cgal.certs)) {
    if (status !== 'OK') {
      throw new Error(`${name}: CGAL ${label.toUpperCase()} cert failed (${status})`);
    }
  }
}

function number(value) {
  const parsed = Number(value);
  if (!Number.isFinite(parsed)) {
    throw new Error(`Expected finite number, got ${value}`);
  }
  return parsed;
}

function integer(value) {
  return Math.trunc(number(value));
}

function point(value) {
  const match = /^\((.*),\s*(.*)\)$/.exec(value || '');
  if (!match) {
    throw new Error(`Expected point, got ${value}`);
  }
  return [number(match[1]), number(match[2])];
}
