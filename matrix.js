import SingularValueDecomposition from './SingularValueDecomposition.js';

class Matrix {
  rows;
  cols;
  data;

  constructor(rows, cols, data = []) {
    this.rows = rows;
    this.cols = cols;
    if (data.length === 0) {
      this.data = new Array(rows);
      for (let i = 0; i < rows; i++) {
        this.data[i] = new Array(cols).fill(0);
      }
    } else {
      this.data = new Array(rows);
      for (let i = 0; i < rows; i++) {
        this.data[i] = new Array(cols);
        for (let j = 0; j < cols; j++) {
          this.data[i][j] = data[i][j];
        }
      }
    }
  }

  static create(data) {
    // Assuming a 2D array
    const rows = data.length;
    const cols = data[0].length;
    return new Matrix(rows, cols, data);
  }

  toString(){
    return this.rows + "x" + this.cols + " matrix";
    // :\n" + this.data;
  }

  /**
   * Returns a square diagonal matrix, whose entries along the diagonal are specified by the passed-in array, and the
   * other entries are 0.
   * @public
   *
   * @param {Array.<number>} diagonalValues
   * @returns {Matrix}
   */
  static diagonalMatrix(diagonalValues) {
    const n = diagonalValues.length;
    const result = new Matrix(n, n);
    for (let i = 0; i < n; i++) {
      result.data[i][i] = diagonalValues[i];
    }
    return result;
  }

  get(row, col) {
    return this.data[row][col];
  }

  set(row, col, value) {
    this.data[row][col] = value;
  }

  add(otherMatrix) {
    if (this.rows !== otherMatrix.rows || this.cols !== otherMatrix.cols) {
      throw new Error("Matrix dimensions must match for addition");
    }

    const result = new Matrix(this.rows, this.cols);
    for (let i = 0; i < this.rows; i++) {
      for (let j = 0; j < this.cols; j++) {
        result.data[i][j] = this.data[i][j] + otherMatrix.data[i][j];
      }
    }
    return result;
  }

  subtract(otherMatrix) {
    if (this.rows !== otherMatrix.rows || this.cols !== otherMatrix.cols) {
      throw new Error("Matrix dimensions must match for subtraction");
    }

    const result = new Matrix(this.rows, this.cols);
    for (let i = 0; i < this.rows; i++) {
      for (let j = 0; j < this.cols; j++) {
        result.data[i][j] = this.data[i][j] - otherMatrix.data[i][j];
      }
    }
    return result;
  }

  scale(scalar) {
    for (let i = 0; i < this.rows; i++) {
      for (let j = 0; j < this.cols; j++) {
        this.data[i][j] *= scalar;
      }
    }
    return this;
  }

  multiply(otherMatrix) {
    if (this.cols !== otherMatrix.rows) {
      throw new Error("Number of columns in first matrix must match number of rows in second matrix");
    }

    const result = new Matrix(this.rows, otherMatrix.cols);
    for (let i = 0; i < this.rows; i++) {
      for (let j = 0; j < otherMatrix.cols; j++) {
        let sum = 0;
        for (let k = 0; k < this.cols; k++) {
          sum += this.data[i][k] * otherMatrix.data[k][j];
        }
        result.data[i][j] = sum;
      }
    }
    return result;
  }

  transpose() {
    const transposed = new Matrix(this.cols, this.rows);

    for (let i = 0; i < this.cols; i++) {
      for (let j = 0; j < this.rows; j++) {
        transposed.data[i][j] = this.data[j][i];
      }
    }

    return transposed;
  }

  fNorm() {
    let sum = 0;
    for (let i = 0; i < this.rows; i++) {
      for (let j = 0; j < this.cols; j++) {
        sum += this.data[i][j] * this.data[i][j];
      }
    }
    return Math.sqrt(sum);
  }

	static identity(rows, cols) {
		const A = new Matrix(rows, cols);
		for (let i = 0; i < rows; i++) {
			for (let j = 0; j < cols; j++) {
				A[i][j] = (i == j ? 1.0 : 0.0);
			}
		}
		return A;
	}

  svd() {
    const { u, q, v } = SVD(this.data);
    return { u, q, v };
  }

  inverse() {
    if (this.rows === this.cols) {
      const n = this.rows;
      const augmented = new Matrix(n, 2 * n);

      // Create augmented matrix
      for (let i = 0; i < n; i++) {
        for (let j = 0; j < n; j++) {
          augmented.data[i][j] = this.data[i][j];
        }
        augmented.data[i][i + n] = 1;
      }

      // Gaussian elimination with partial pivoting
      for (let i = 0; i < n - 1; i++) {
        // Find pivot row
        let maxRow = i;
        for (let j = i + 1; j < n; j++) {
          if (Math.abs(augmented.data[j][i]) > Math.abs(augmented.data[maxRow][i])) {
            maxRow = j;
          }
        }

        // Swap rows
        if (maxRow !== i) {
          augmented.#swapRows(i, maxRow);
        }

        // Elimination
        for (let j = i + 1; j < n; j++) {
          const factor = augmented.data[j][i] / augmented.data[i][i];
          for (let k = i; k < 2 * n; k++) {
            augmented.data[j][k] -= factor * augmented.data[i][k];
          }
        }
      }

      // Back substitution
      for (let i = n - 1; i >= 0; i--) {
        for (let j = i + 1; j < n; j++) {
          augmented.data[i][n + i] -= augmented.data[i][j] * augmented.data[j][n + j];
        }
        augmented.data[i][n + i] /= augmented.data[i][i];
      }

      // Extract inverse
      const inverse = new Matrix(n, n);
      for (let i = 0; i < n; i++) {
        for (let j = 0; j < n; j++) {
          inverse.data[i][j] = augmented.data[i][j + n];
        }
      }
      return inverse;
    } else {
      return SingularValueDecomposition.pseudoinverse(this);
    }
  }

  // Helper function to swap rows
  #swapRows(row1, row2) {
    const temp = this.data[row1];
    this.data[row1] = this.data[row2];
    this.data[row2] = temp;
  }
}

export { Matrix }
