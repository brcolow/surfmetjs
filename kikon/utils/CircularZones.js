class CircularZones {
  constructor() {
    this._intervals = [];
  }

  // For the two static factory functions below: dir1 and dir2 are two directions in the plane, defined
  // by the angles measured anti-clockwise from the x-axis (thus dir1 and dir2 are two angles in radians, in the
  // interval [-pi, pi]). We may however imagine the directions as versors starting from the origin. No relation
  // can be assumed between dir1 and dir2 (e.g., either one may be the largest of the two angles).
  //
  // The two functions construct a CircularZones object corresponding to the angle between the two directions.
  // The CircularZones object contains a "canonical" array of "correct" intervals. An interval is "correct", if the
  // first value is smaller than the second, and they are both in the interval [-pi, pi]. The array is "canonical" if no
  // two intervals overlap and they are sorted ascending, i.e., the max value of an interval is strictly smaller than the
  // min value of the next.
  //
  // Why a set of intervals and not just one interval that just sorts in ascending order
  // the set [dir1, dir2]? This in fact is the solution in many cases, that is [[min(dir1, dir2), max(dir1, dir2)]],
  // the array od intervals consists of just one interval. However, if the angle contains the -pi=pi line (the direction
  // along the negative x-axis, or towards West), it should be split in two "correct" intervals: for instance, the angle from
  // 3*pi/2 (NW) to -3*pi/2 (SW) will become [[-pi, -3*pi/2], [3*pi/2, pi]].
  //
  // And, of course, there's an ambiguity about the actual interval: there are always two intervals between the two directions,
  // the one from dir1 to dir2 and the one from dir2 to dir1.
  // The first function, betweenDirectionsExteriorAngle always computes the exterior angle, that is the one that is larger than
  // 180 degrees. The second function, betweenDirectionsContainingDirection may be called with two directions dir1 and dir2 that
  // are completely opposed to each other, so there are two angles of 180 degrees; it chooses that one that contains a third
  // direction, dir0

  static fullCircle() {
    const newCZs = new CircularZones();
    newCZs._intervals = [[-Math.PI, Math.PI]];
    return newCZs;
  }

  static betweenDirectionsExteriorAngle(dir1, dir2) {
    // returns the angle intervals corresponding to the exterior angle between the directions (angles) dir1 and dir2
    // the exterior angle is the angle that is > 180 deg
    const newCZs = new CircularZones();
    const dirMin = Math.min(dir1, dir2),
      dirMax = Math.max(dir1, dir2);
    if (dirMax - dirMin > Math.PI) {
      newCZs._intervals = [[dirMin, dirMax]];
    } else {
      newCZs._intervals = [
        [-Math.PI, dirMin],
        [dirMax, Math.PI],
      ];
    }
    return newCZs;
  }

  static betweenDirectionsContainingDirection(dir1, dir2, dir0) {
    // returns the angle intervals corresponding to the angle between the directions (angles) dir1 and dir2 that contains dir0
    const newCZs = new CircularZones();
    const dirMin = Math.min(dir1, dir2),
      dirMax = Math.max(dir1, dir2);
    if (dirMin < dir0 && dir0 < dirMax) {
      newCZs._intervals = [[dirMin, dirMax]];
    } else {
      newCZs._intervals = [
        [-Math.PI, dirMin],
        [dirMax, Math.PI],
      ];
    }
    return newCZs;
  }

  static directionPlusMinus(dir, deltaDir) {
    // returns the intervals corresponding to the angle from dir - deltaDir to dir + deltaDir
    // assumes deltaDir < Math.PI
    const newCZs = new CircularZones();
    const dirMin = dir - deltaDir,
      dirMax = dir + deltaDir;
    if (dirMin < -Math.PI) {
      newCZs._intervals = [
        [-Math.PI, dirMax],
        [dirMin + 2 * Math.PI, Math.PI],
      ];
    } else if (dirMax > Math.PI) {
      newCZs._intervals = [
        [-Math.PI, dirMax - 2 * Math.PI],
        [dirMin, Math.PI],
      ];
    } else {
      newCZs._intervals = [[dirMin, dirMax]];
    }
    return newCZs;
  }

  unionWith(otherIntervals) {
    for (const newInterval of otherIntervals._intervals) {
      this._addInterval(newInterval);
    }
  }

  //   |-------|
  //...........       dMin
  //    ...........   dMax
  intersectWith(otherIntervals) {
    const newCZs = new CircularZones();
    //console.warn('intersect', JSON.stringify(this._intervals), JSON.stringify(otherIntervals._intervals))
    for (const [dMin, dMax] of otherIntervals._intervals) {
      const intersect_intervals = [];
      for (const [d0Min, d0Max] of this._intervals) {
        if (dMin <= d0Max && dMax >= d0Min) {
          const minInterior = dMin >= d0Min,
            maxInterior = dMax <= d0Max;
          if (minInterior && maxInterior) {
            intersect_intervals.push([dMin, dMax]);
            break;
          } else if (minInterior) {
            intersect_intervals.push([dMin, d0Max]);
            // no break here
          } else if (maxInterior) {
            intersect_intervals.push([d0Min, dMax]);
            break;
          } else {
            intersect_intervals.push([d0Min, d0Max]);
            // no break here
          }
        } else if (intersect_intervals.length > 0) {
          break;
        }
      }
      for (const interval of intersect_intervals) {
        newCZs._addInterval(interval);
      }
    }
    //console.warn('==>', JSON.stringify(newCZs._intervals))
    this._intervals = newCZs._intervals;
  }

  _addInterval([dMin, dMax]) {
    for (let i = 0; i < this._intervals.length; i++) {
      const [d0Min, d0Max] = this._intervals[i];
      if (dMax < d0Min) {
        this._intervals.splice(i, 0, [dMin, dMax]);
        return;
      }
      const minInterior = d0Min <= dMin && dMin <= d0Max,
        maxInterior = d0Min <= dMax && dMax <= d0Max;
      if (maxInterior) {
        if (minInterior) {
          return;
        }
        this._intervals[i][0] = dMin;
        return;
      }
      if (minInterior) {
        let newMax = dMax,
          nToDelete = 0;
        for (let j = i + 1; j < this._intervals.length; j++) {
          const [d1Min, d1Max] = this._intervals[j];
          if (dMax < d1Min) {
            break;
          }
          if (dMax > d1Max) {
            nToDelete++;
          } else {
            nToDelete++;
            newMax = d1Max;
            break;
          }
        }
        this._intervals.splice(i, nToDelete);
        this._intervals[i][1] = newMax;
        return;
      }
    }
    this._intervals.push([dMin, dMax]);
  }

  static _angularAbsDiff(dir1, dir2) {
    return Math.min(
      Math.abs(dir1 - dir2),
      Math.abs(dir1 - Math.PI) + Math.abs(dir2 - Math.PI)
    );
  }

  isDirectionInsideOrClosest(dir) {
    const closest = {
      delta: 1 / 0,
      interval: null,
      toMinInterval: null,
      inside: true,
    };
    for (const [dMin, dMax] of this._intervals) {
      if (dir >= dMin && dir <= dMax) {
        return closest;
      }
      const deltaMin = CircularZones._angularAbsDiff(dMin, dir);
      if (deltaMin < closest.delta) {
        closest.delta = deltaMin;
        closest.interval = [dMin, dMax];
        closest.toMinInterval = true;
      }
      const deltaMax = CircularZones._angularAbsDiff(dMax, dir);
      if (deltaMax < closest.delta) {
        closest.delta = deltaMax;
        closest.interval = [dMin, dMax];
        closest.toMinInterval = false;
      }
    }
    closest.inside = false;
    return closest;
  }

  isVoid() {
    return this._intervals.length === 0;
  }
}
