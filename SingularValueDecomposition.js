// Copyright 2013-2022, University of Colorado Boulder

/**
 * SVD decomposition, based on Jama (http://math.nist.gov/javanumerics/jama/)
 *
 * @author Jonathan Olson <jonathan.olson@colorado.edu>
 */

import { Matrix } from './matrix.js';

const ArrayType = window.Float64Array || Array;

class SingularValueDecomposition {
  U;
  V;
  s;
  m;
  n;

  /**
   * @param {Matrix} matrix
   */
  constructor(matrix) {
    this.matrix = matrix;
    console.log(this.matrix.toString());
    // Derived from LINPACK code.
    // Initialize.
    const A = structuredClone(matrix.data);
    this.m = matrix.rows;
    this.n = matrix.cols;
    const m = this.m;
    const n = this.n;

    /* Apparently the failing cases are only a proper subset of (m<n),
     so let's not throw error.  Correct fix to come later?
     if (m<n) {
     throw new IllegalArgumentException("Jama SVD only works for m >= n"); }
     */
    const nu = Math.min( m, n );
    this.s = new ArrayType( Math.min( m + 1, n ) ).fill(0);
    const s = this.s;
    this.U = this.#twoDArray(m, nu);
    const U = this.U;
    this.V = this.#twoDArray(n, n);
    const V = this.V;
    const e = new ArrayType( n ).fill(0);
    const work = new ArrayType( m ).fill(0);
    const wantu = true;
    const wantv = true;

    let cs;
    let sn;

    // Reduce A to bidiagonal form, storing the diagonal elements
    // in s and the super-diagonal elements in e.

    const nct = Math.min( m - 1, n );
    const nrt = Math.max( 0, Math.min( n - 2, m ) );
    for (let k = 0; k < Math.max( nct, nrt ); k++ ) {
      if ( k < nct ) {

        // Compute the transformation for the k-th column and
        // place the k-th diagonal in s[k].
        // Compute 2-norm of k-th column without under/overflow.
        s[ k ] = 0;
        for (let i = k; i < m; i++ ) {
          s[ k ] = SingularValueDecomposition.hypot( s[ k ], A[ i * n + k ] );
        }
        if ( s[ k ] !== 0.0 ) {
          if ( A[ k * n + k ] < 0.0 ) {
            s[ k ] = -s[ k ];
          }
          for (let i = k; i < m; i++ ) {
            A[i][k] /= s[k];
          }
          A[k][k] += 1.0;
        }
        s[ k ] = -s[ k ];
      }
      for (let j = k + 1; j < n; j++ ) {
        if ( ( k < nct ) && ( s[ k ] !== 0.0 ) ) {

          // Apply the transformation.

          let t = 0;
          for (let i = k; i < m; i++ ) {
            t += A[i][k] * A[i][j];
        }
        t = -t / A[k][k];
        for (let i = k; i < m; i++ ) {
            A[i][j] += t * A[i][k];
          }
        }

        // Place the k-th row of A into e for the
        // subsequent calculation of the row transformation.

        e[ j ] = A[k][j];
      }
      if ( wantu && ( k < nct ) ) {

        // Place the transformation in U for subsequent back
        // multiplication.

        for (let i = k; i < m; i++ ) {
          U[i][k] = A[i][k];
        }
      }
      if ( k < nrt ) {

        // Compute the k-th row transformation and place the
        // k-th super-diagonal in e[k].
        // Compute 2-norm without under/overflow.
        e[ k ] = 0;
        for (let i = k + 1; i < n; i++ ) {
          e[ k ] = SingularValueDecomposition.hypot( e[ k ], e[ i ] );
        }
        if ( e[ k ] !== 0.0 ) {
          if ( e[ k + 1 ] < 0.0 ) {
            e[ k ] = -e[ k ];
          }
          for (let i = k + 1; i < n; i++ ) {
            e[ i ] /= e[ k ];
          }
          e[ k + 1 ] += 1.0;
        }
        e[ k ] = -e[ k ];
        if ( ( k + 1 < m ) && ( e[ k ] !== 0.0 ) ) {

          // Apply the transformation.

          for (let i = k + 1; i < m; i++ ) {
            work[ i ] = 0.0;
          }
          for (let j = k + 1; j < n; j++ ) {
            for (let i = k + 1; i < m; i++ ) {
                work[i] += e[j] * A[i][j];
            }
          }
          for (let j = k + 1; j < n; j++ ) {
            let t = -e[ j ] / e[ k + 1 ];
            for (let i = k + 1; i < m; i++ ) {
                A[i][j] += t * work[i];
            }
          }
        }
        if ( wantv ) {

          // Place the transformation in V for subsequent
          // back multiplication.

          for (let i = k + 1; i < n; i++ ) {
            V[i][k] = e[i];
          }
        }
      }
    }

    // Set up the final bidiagonal matrix or order p.

    let p = Math.min( n, m + 1 );
    if ( nct < n ) {
        s[nct] = A[nct][nct];
    }
    if ( m < p ) {
      s[ p - 1 ] = 0.0;
    }
    if ( nrt + 1 < p ) {
        e[nrt] = A[nrt][p - 1];
    }
    e[ p - 1 ] = 0.0;

    // If required, generate U.

    if ( wantu ) {
      for (let j = nct; j < nu; j++ ) {
        for (let i = 0; i < m; i++ ) {
          U[i][j] = 0.0;
        }
        U[j][j] = 1.0;
      }
      for (let k = nct - 1; k >= 0; k-- ) {
        if ( s[ k ] !== 0.0 ) {
          for (let j = k + 1; j < nu; j++ ) {
            let t = 0;
            for (let i = k; i < m; i++ ) {
              t += U[i][k] * U[i][j];
            }
            t = -t / U[k][k];
            for (let i = k; i < m; i++ ) {
                U[i][j] += t * U[i][k];
            }
          }
          for (let i = k; i < m; i++ ) {
            U[i][k] = -U[i][k];
          }
          U[k][k] = 1.0 + U[k][k];
          for (let i = 0; i < k - 1; i++ ) {
            U[i][k] = 0.0;
          }
        }
        else {
          for (let i = 0; i < m; i++ ) {
            U[i][k] = 0.0;
          }
          U[k][k] = 1.0;
        }
      }
    }

    // If required, generate V.

    if ( wantv ) {
      for (let k = n - 1; k >= 0; k-- ) {
        if ( ( k < nrt ) && ( e[ k ] !== 0.0 ) ) {
          for (let j = k + 1; j < nu; j++ ) {
            let t = 0;
            for (let i = k + 1; i < n; i++ ) {
                t += V[i][k] * V[i][j];
            }
            t = -t / V[k + 1][k];
            for (let i = k + 1; i < n; i++ ) {
                V[i][j] += t * V[i][k];
            }
          }
        }
        for (let i = 0; i < n; i++ ) {
          V[i][k] = 0.0;
        }
        V[k][k] = 1.0;
      }
    }

    // Main iteration loop for the singular values.

    const pp = p - 1;
    let iter = 0;
    const eps = Math.pow( 2.0, -52.0 );
    const tiny = Math.pow( 2.0, -966.0 );
    while ( p > 0 ) {
      let k, kase;

      // Here is where a test for too many iterations would go.
      if ( iter > 10 ) {
        break;
      }

      // This section of the program inspects for
      // negligible elements in the s and e arrays.  On
      // completion the variables kase and k are set as follows.

      // kase = 1   if s(p) and e[k-1] are negligible and k<p
      // kase = 2   if s(k) is negligible and k<p
      // kase = 3   if e[k-1] is negligible, k<p, and
      //        s(k), ..., s(p) are not negligible (qr step).
      // kase = 4   if e(p-1) is negligible (convergence).

      for (k = p - 2; k >= -1; k-- ) {
        console.log(k);
        if ( k === -1 ) {
          break;
        }
        if ( Math.abs( e[ k ] ) <=
             tiny + eps * ( Math.abs( s[ k ] ) + Math.abs( s[ k + 1 ] ) ) ) {
          e[ k ] = 0.0;
          break;
        }
      }
      if ( k === p - 2 ) {
        kase = 4;
      }
      else {
        let ks;
        for (let ks = p - 1; ks >= k; ks-- ) {
          if ( ks === k ) {
            break;
          }
          const t = ( ks !== p ? Math.abs( e[ ks ] ) : 0 ) +
              ( ks !== k + 1 ? Math.abs( e[ ks - 1 ] ) : 0 );
          if ( Math.abs( s[ ks ] ) <= tiny + eps * t ) {
            s[ ks ] = 0.0;
            break;
          }
        }
        if ( ks === k ) {
          kase = 3;
        } else if ( ks === p - 1 ) {
          kase = 1;
        } else {
          kase = 2;
          k = ks;
        }
      }
      k++;

      // Perform the task indicated by kase.

      switch( kase ) {

        // Deflate negligible s(p).

        case 1: {
          f = e[ p - 2 ];
          e[ p - 2 ] = 0.0;
          for (let j = p - 2; j >= k; j-- ) {
            t = SingularValueDecomposition.hypot( s[ j ], f );
            cs = s[ j ] / t;
            sn = f / t;
            s[ j ] = t;
            if ( j !== k ) {
              f = -sn * e[ j - 1 ];
              e[ j - 1 ] = cs * e[ j - 1 ];
            }
            if ( wantv ) {
              for (let i = 0; i < n; i++ ) {
                t = cs * V[i][j] + sn * V[i][p - 1];
                V[i][p - 1] = -sn * V[i][j] + cs * V[i][p - 1];
                V[i][j] = t;
              }
            }
          }
        }
          break;

        // Split at negligible s(k).

        case 2: {
          let f = e[ k - 1 ];
          e[ k - 1 ] = 0.0;
          for (let j = k; j < p; j++ ) {
            t = SingularValueDecomposition.hypot( s[ j ], f );
            cs = s[ j ] / t;
            sn = f / t;
            s[ j ] = t;
            f = -sn * e[ j ];
            e[ j ] = cs * e[ j ];
            if ( wantu ) {
              for (let i = 0; i < m; i++ ) {
                t = cs * U[i][j] + sn * U[i][k - 1];
                U[i][k - 1] = -sn * U[i][j] + cs * U[i][k - 1];
                U[i][j] = t;
              }
            }
          }
        }
          break;

        // Perform one qr step.

        case 3: {

          // Calculate the shift.

          const scale = Math.max(Math.max(Math.max(Math.max(Math.abs( s[ p - 1 ] ), Math.abs( s[ p - 2 ] ) ), Math.abs( e[ p - 2 ] ) ), Math.abs( s[ k ] ) ), Math.abs( e[ k ] ) );
          const sp = s[ p - 1 ] / scale;
          const spm1 = s[ p - 2 ] / scale;
          const epm1 = e[ p - 2 ] / scale;
          const sk = s[ k ] / scale;
          const ek = e[ k ] / scale;
          const b = ( ( spm1 + sp ) * ( spm1 - sp ) + epm1 * epm1 ) / 2.0;
          const c = ( sp * epm1 ) * ( sp * epm1 );
          let shift = 0.0;
          if ( ( b !== 0.0 ) || ( c !== 0.0 ) ) {
            shift = Math.sqrt( b * b + c );
            if ( b < 0.0 ) {
              shift = -shift;
            }
            shift = c / ( b + shift );
          }
          f = ( sk + sp ) * ( sk - sp ) + shift;
          let g = sk * ek;

          // Chase zeros.

          for (let j = k; j < p - 1; j++ ) {
            t = SingularValueDecomposition.hypot( f, g );
            cs = f / t;
            sn = g / t;
            if ( j !== k ) {
              e[ j - 1 ] = t;
            }
            f = cs * s[ j ] + sn * e[ j ];
            e[ j ] = cs * e[ j ] - sn * s[ j ];
            g = sn * s[ j + 1 ];
            s[ j + 1 ] = cs * s[ j + 1 ];
            if ( wantv ) {
              for (let i = 0; i < n; i++ ) {
                t = cs * V[i][j] + sn * V[i][j + 1];
                V[i][j + 1] = -sn * V[i][j] + cs * V[i][j + 1];
                V[i][j] = t;
              }
            }
            t = SingularValueDecomposition.hypot( f, g );
            cs = f / t;
            sn = g / t;
            s[ j ] = t;
            f = cs * e[ j ] + sn * s[ j + 1 ];
            s[ j + 1 ] = -sn * e[ j ] + cs * s[ j + 1 ];
            g = sn * e[ j + 1 ];
            e[ j + 1 ] = cs * e[ j + 1 ];
            if ( wantu && ( j < m - 1 ) ) {
              for (let i = 0; i < m; i++ ) {
                t = cs * U[i][j] + sn * U[i][j + 1];
                U[i][j + 1] = -sn * U[i][j] + cs * U[i][j + 1];
                U[i][j] = t;
              }
            }
          }
          e[ p - 2 ] = f;
          iter = iter + 1;
        }
          break;

        // Convergence.

        case 4: {

          // Make the singular values positive.

          if ( s[ k ] <= 0.0 ) {
            s[ k ] = ( s[ k ] < 0.0 ? -s[ k ] : 0.0 );
            if ( wantv ) {
              for (let i = 0; i <= pp; i++ ) {
                V[i][k] = -V[i][k];
              }
            }
          }

          // Order the singular values.

          while ( k < pp ) {
            if ( s[ k ] >= s[ k + 1 ] ) {
              break;
            }
            t = s[ k ];
            s[ k ] = s[ k + 1 ];
            s[ k + 1 ] = t;
            if ( wantv && ( k < n - 1 ) ) {
              for (let i = 0; i < n; i++ ) {
                t = V[i][k + 1];
                V[i][k + 1] = V[i][k];
                V[i][k] = t;
              }
            }
            if ( wantu && ( k < m - 1 ) ) {
              for (let i = 0; i < m; i++ ) {
                t = U[i][k + 1];
                U[i][k + 1] = U[i][k];
                U[i][k] = t;
              }
            }
            k++;
          }
          iter = 0;
          p--;
        }
          break;

        default:
          throw new Error( `invalid kase: ${kase}` );
      }

    }
  }

 #twoDArray(rows, cols) {
    let x = new Array(rows);
  
    for (let i = 0; i < rows; i++) {
      x[i] = new ArrayType(cols).fill(0);
    }
    return x;
  }

  /**
   * @public
   *
   * @returns {Matrix}
   */
  getU() {
    return new Matrix( this.m, Math.min( this.m + 1, this.n ), this.U);
  }

  /**
   * @public
   *
   * @returns {Matrix}
   */
  getV() {
    return new Matrix( this.n, this.n, this.V);
  }

  /**
   * @public
   *
   * @returns {Array.<number>}
   */
  getSingularValues() {
    return this.s;
  }

  /**
   * @public
   *
   * @returns {Matrix}
   */
  getS() {
    const result = new Matrix( this.n, this.n );
    for ( let i = 0; i < this.n; i++ ) {
      for ( let j = 0; j < this.n; j++ ) {
        result.data[i][j] = 0.0;
      }
      result.data[i, i] = this.s[ i ];
    }
    return result;
  }

  /**
   * @public
   *
   * @returns {number}
   */
  norm2() {
    return this.s[ 0 ];
  }

  /**
   * @public
   *
   * @returns {number}
   */
  cond() {
    return this.s[ 0 ] / this.s[ Math.min( this.m, this.n ) - 1 ];
  }

  /**
   * @public
   *
   * @returns {number}
   */
  rank() {
    // changed to 23 from 52 (bits of mantissa), since we are using floats here!
    const eps = Math.pow( 2.0, -23.0 );
    const tol = Math.max( this.m, this.n ) * this.s[ 0 ] * eps;
    let r = 0;
    for ( let i = 0; i < this.s.length; i++ ) {
      if ( this.s[ i ] > tol ) {
        r++;
      }
    }
    return r;
  }

  /**
   * sqrt(a^2 + b^2) without under/overflow.
   * @public
   *
   * @param {number} a
   * @param {number} b
   * @returns {number}
   */
  static hypot( a, b ) {
    let r;
    if ( Math.abs( a ) > Math.abs( b ) ) {
      r = b / a;
      r = Math.abs( a ) * Math.sqrt( 1 + r * r );
    }
    else if ( b !== 0 ) {
      r = a / b;
      r = Math.abs( b ) * Math.sqrt( 1 + r * r );
    }
    else {
      r = 0.0;
    }
    return r;
  }

  /**
   * Constructs the Moore-Penrose pseudoinverse of the specified matrix, using the SVD construction.
   * @public
   *
   * See https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_pseudoinverse for details. Helpful for
   * linear least-squares regression.
   *
   * @param {Matrix} matrix, m x n
   * @returns {Matrix} - n x m
   */
  static pseudoinverse( matrix ) {
    console.log("computing pseudoinverse of matrix: " + matrix.toString());
    const svd = new SingularValueDecomposition( matrix );
    console.log("singular values: " + svd.getSingularValues());
    const sigmaPseudoinverse = Matrix.diagonalMatrix( svd.getSingularValues().map( value => {
      if ( Math.abs( value ) < 1e-300 ) {
        console.log(0);
        return 0;
      }
      else {
        console.log((1 / value));
        return 1 / value;
      }
    } ) );
    console.log("sigmaPseudoinverse rows: " + sigmaPseudoinverse);
    return svd.getV().multiply( sigmaPseudoinverse ).multiply( svd.getU().transpose() );
  }
}

export default SingularValueDecomposition;