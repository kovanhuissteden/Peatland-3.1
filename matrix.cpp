/***************************************************************************
                          matrix.cpp  -  description
                             -------------------
    begin                : |02-01-2003|
    copyright            : (C) |2003| by |J. van Huissteden|
    email                : |ko.van.huissteden@geo.falw.vu.nl|
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <cstring>
using namespace std;
#include "matrix.h"
#include "general.h"

#define ERROR1    "Matrix memory allocation error"
#define ERROR2    "Matrix index out of range "
#define ERROR3    "This insert procedure only works with vectors"
#define ERROR4    "Addition: matrices must be of equal size"
#define ERROR5    "Cannot copy data, source matrix exceeds size of destination matrix"
#define ERROR6    "Matrix multiplication: size of matrices incompatible."
#define ERROR7    " Offending index: "
#define ERROR8    "Zero or negative row/col size."

Matrix::Matrix()                 // initializes a 1 x 1 matrix
{
  data = new double;
  if (data == NULL)
  {
    cout << ERROR1 << endl;
    exit(EXIT_FAILURE);
  }
  *data = 0.0;
  cols = 1;
  rows = 1;
  is_scalar = TRUE;
}

Matrix::Matrix(int r)            // initializes a vector with r elements
{
  int i;
  double *p;

  data = new double[r];
  if (data == NULL)
  {
    cout << ERROR1 << endl;
    exit(EXIT_FAILURE);
  }
  cols = r;
  rows = 1;
  p = data;
    for (i = 0; i < r; i++) {*p = 0.0; p++;}
  if (cols == 1) is_scalar = TRUE; else is_scalar = FALSE;
}

Matrix::Matrix(int r, int c)     // initializes matrix with r rows and c columns
{
  int i;
  double *p;

  data = new double[r*c];
  if (data == NULL)
  {
    cout << ERROR1 << endl;
    exit(EXIT_FAILURE);
  }
  rows = r;
  cols = c;
  p = data;
  for (i = 0; i < (r * c); i++) {*p = 0.0; p++;}
  if ((cols == 1) && (rows == 1)) is_scalar = TRUE; else is_scalar = FALSE;
}


Matrix::Matrix(Matrix &m)   // copy constructor
{
  rows = m.rows;
  cols = m.cols;
  is_scalar = m.is_scalar;
  data = new double[rows * cols];
  if (data == NULL)
  {
    cout << ERROR1 << endl;
    exit(EXIT_FAILURE);
  }
  memcpy(data, m.data, rows * cols * sizeof(double));
}
 


Matrix::Matrix(int r, double *d)  // initializes avector with r elements and data specified in *d
{

  data = new double[r];
  if (data == NULL)
  {
    cout << ERROR1 << endl;
    exit(EXIT_FAILURE);
  }
  cols = r;
  rows = 1;
  memcpy(data, d, cols * sizeof(double));
  if (cols == 1) is_scalar = TRUE; else is_scalar = FALSE;
}

Matrix::Matrix(int r, int c, double *d)     // initializes matrix with r rows and c columns and data specified in *d
{

  data = new double[r*c];
  if (data == NULL)
  {
    cout << ERROR1 << endl;
    exit(EXIT_FAILURE);
  }
  rows = r;
  cols = c;
  memcpy(data, d, rows * cols * sizeof(double));
  if ((cols == 1) && (rows == 1)) is_scalar = TRUE; else is_scalar = FALSE;
}

double &Matrix::operator()(int i)           // Matlab-style index operator
{
  if ((i > cols) || (i < 1))
  {
     cout << ERROR2 << cols << ERROR7 << i << endl;
  }
  return *(data + (i - 1) * rows);
}


double &Matrix::operator()(int i, int j)    // Matlab-style index operator
{
  if ((i > rows) || (j > cols) || (i < 1) || (j < 1))
  {
    cout << ERROR2 << rows << " " << cols << ERROR7 << i << " " << j << endl;
  }
  return *(data + (i - 1) * cols + j - 1);
}



Matrix &Matrix::operator=(const Matrix &m)   // assignment operator
{
  delete data;
  rows = m.rows;
  cols = m.cols;
  is_scalar = m.is_scalar;
  data = new double[rows * cols];
  if (data == NULL)
  {
    cout << ERROR1 << endl;
    exit(EXIT_FAILURE);
  }
  memcpy(data, m.data, rows * cols * sizeof(double));
  return *this;
}

Matrix &Matrix::operator+=(const Matrix &m)    // addition operator
{
  int i, j;
  double *p, *d;

  d = data;                                  // assign pointers to data and do the element-by element addition
  p = m.data;
  if (m.is_scalar)                          // Matrix with one element is treated as addition with scalar
  {
      for (i = 0; i < rows; i++) for (j = 0; j < cols; j++) {*d = *d + *p; d++;}
    return *this;
  }
  if ((cols != m.cols) || (rows != m.rows))  // check size if matrix is larger
  {
    cout << ERROR4 << endl;
    exit(EXIT_FAILURE);
  }
    for (i = 0; i < rows; i++) for (j = 0; j < cols; j++) {*d = *d + *p; d++; p++;}
  return *this;
}


Matrix &Matrix::operator-=(const Matrix &m)    // substraction operator
{
  int i, j;
  double *p, *d;

  d = data;                                  // assign pointers to data and do the element-by element substraction
  p = m.data;
  if (m.is_scalar)                          // Matrix with one element is treated as substraction with scalar
  {
      for (i = 0; i < rows; i++) for (j = 0; j < cols; j++) {*d = *d - *p; d++;}
    return *this;
  }
  if ((cols != m.cols) || (rows != m.rows))  // check size if matrix is larger
  {
    cout << ERROR4 << endl;
    exit(EXIT_FAILURE);
  }
    for (i = 0; i < rows; i++) for (j = 0; j < cols; j++) {*d = *d - *p; d++; p++;}
  return *this;
}

Matrix &Matrix::operator+=(const double m)    // addition of scalar
{
  int i, j;
  double *d;

  d = data;
    for (i = 0; i < rows; i++) for (j = 0; j < cols; j++) {*d = *d + m; d++;}
  return *this;
}

Matrix &Matrix::operator-=(const double m)    // substraction of scalar
{
  int i, j;
  double *d;

  d = data;
    for (i = 0; i < rows; i++) for (j = 0; j < cols; j++) {*d = *d - m; d++;}
  return *this;
}


Matrix &Matrix::operator*=(const Matrix &m)   // element by element multiplication operator - if m is a1 x 1 matrix it is treated as a scalar
{
  int i, j;
  double *p, *d;

  d = data;                                  // assign pointers to data and do the element-by element multiplication
  p = m.data;
  if (m.is_scalar)                           // Matrix with one element is treated as multiplication with scalar
  {
      for (i = 0; i < rows; i++) for (j = 0; j < cols; j++) {*d = *d * *p; d++;}
    return *this;
  }
  if ((cols != m.cols) || (rows != m.rows))  // check size if matrix is larger
  {
    cout << ERROR4 << endl;
    exit(EXIT_FAILURE);
  }
    for (i = 0; i < rows; i++) for (j = 0; j < cols; j++) {*d = *d * *p; d++; p++;}
  return *this;
}


Matrix &Matrix::operator*=(const double m)    // multiplication by scalar
{
  int i, j;
  double *d;

  d = data;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            *d = *d * m;
            d++;
        }
    }
  return *this;
}

Matrix &Matrix::operator/=(const Matrix &m)   // element by element division operator - if m is a1 x 1 matrix it is treated as a scalar
{
  int i, j;
  double *p, *d;

  d = data;                                  // assign pointers to data and do the element-by element division
  p = m.data;
  if (m.is_scalar)                           // Matrix with one element is treated as division with scalar
  {
      for (i = 0; i < rows; i++) for (j = 0; j < cols; j++) {*d = *d / *p; d++;}
    return *this;
  }
  if ((cols != m.cols) || (rows != m.rows))  // check size if matrix is larger
  {
    cout << ERROR4 << endl;
    exit(EXIT_FAILURE);
  }
    for (i = 0; i < rows; i++) for (j = 0; j < cols; j++) {*d = *d / *p; d++; p++;}
  return *this;
}


Matrix &Matrix::operator/=(const double m)    // division by scalar
{
  int i, j;
  double *d;

  d = data;
    for (i = 0; i < rows; i++) for (j = 0; j < cols; j++) {*d = *d / m; d++;}
  return *this;
}


Matrix &Matrix::operator-()       // unary minus
{
  int i;
  double *d;

  d = data;
  for (i = 0; i < (rows * cols); i++)
  {
    *d = -*d;
    d++;
  }
  return *this;
}


Matrix &operator+(Matrix &a, Matrix &b) // addition operator 2 matrices
// this has to be a non-member function of class Matrix
{
  int i, j, r, c;
  double *pa, *pb, *pc;
  Matrix *m;

  if (!CheckSize(&a, &b, &r, &c)) exit(EXIT_FAILURE);  // check sizes
  pa = a.data;                                  // assign pointers to data and do the element-by element addition
  pb = b.data;
  m = new Matrix;
  m->Resize(r, c);
  pc = m->Data();
  if (a.is_scalar)        // Matrix with one element is treated as addition with scalar
  {
      for (i = 0; i < r; i++) for (j = 0; j < c; j++) {*pc = *pb + *pa; pc++; pb++;}
  } else
  {
    if (b.is_scalar)
    {
        for (i = 0; i < r; i++) for (j = 0; j < c; j++) {*pc = *pa + *pb; pc++; pa++;}
    } else for (i = 0; i < r; i++) for (j = 0; j < c; j++) {*pc = *pa + *pb; pc++; pa++; pb++;}
  }
  return *m;
}

Matrix &operator-(Matrix &a, Matrix &b) // substraction operator 2 matrices
// this has to be a non-member function of class Matrix
{
  int i, j, r, c;
  double *pa, *pb, *pc;
  Matrix *m;

  if (!CheckSize(&a, &b, &r, &c)) exit(EXIT_FAILURE);  // check sizes
  pa = a.data;                                  // assign pointers to data and do the element-by element addition
  pb = b.data;
  m = new Matrix;
  m->Resize(r, c);
  pc = m->Data();
  if (a.is_scalar)                       // Matrix with one element is treated as substraction with scalar
  {
      for (i = 0; i < r; i++) for (j = 0; j < c; j++) {*pc =  *pa - *pb; pc++; pb++;}
  } else
  {
    if (b.is_scalar)
    {
        for (i = 0; i < r; i++) for (j = 0; j < c; j++) {*pc = *pa - *pb; pc++; pa++;}
    } else for (i = 0; i < r; i++) for (j = 0; j < c; j++) {*pc = *pa - *pb; pc++; pa++; pb++;}
  }
  return *m;
}

Matrix &operator*(Matrix &a, Matrix &b) // element-by-element multiplication operator 2 matrices
// this has to be a non-member function of class Matrix
{
  int i, j, r, c;
  double *pa, *pb, *pc;
  Matrix *m;

  if (!CheckSize(&a, &b, &r, &c)) exit(EXIT_FAILURE);  // check sizes
  pa = a.data;                                  // assign pointers to data of matirces
  pb = b.data;
  m = new Matrix;
  m->Resize(r, c);
  pc = m->Data();
  if (a.is_scalar)        // Matrix with one element is treated as addition with scalar
  {
      for (i = 0; i < r; i++) for (j = 0; j < c; j++) {*pc = *pb * *pa; pc++; pb++;}
  } else
  {
    if (b.is_scalar)
    {
        for (i = 0; i < r; i++) for (j = 0; j < c; j++) {*pc = *pa * *pb; pc++; pa++;}
    } else for (i = 0; i < r; i++) for (j = 0; j < c; j++) {*pc = *pa * *pb; pc++; pa++; pb++;}
  }
  return *m;
}

Matrix &operator/(Matrix &a, Matrix &b) // element-by-element division operator 2 matrices
// this has to be a non-member function of class Matrix
{
  int i, j, r, c;
  double *pa, *pb, *pc;
  Matrix *m;

  if (!CheckSize(&a, &b, &r, &c)) exit(EXIT_FAILURE);  // check sizes
  pa = a.data;                                  // assign pointers to data
  pb = b.data;
  m = new Matrix;
  m->Resize(r, c);
  pc = m->Data();
  if (a.is_scalar)        // Matrix with one element is treated as addition with scalar
  {
      for (i = 0; i < r; i++) for (j = 0; j < c; j++) {*pc =  *pa / *pb; pc++; pb++;}
  } else
  {
    if (b.is_scalar)
    {
        for (i = 0; i < r; i++) for (j = 0; j < c; j++) {*pc = *pa / *pb; pc++; pa++;}
    } else for (i = 0; i < r; i++) for (j = 0; j < c; j++) {*pc = *pa / *pb; pc++; pa++; pb++;}
  }
  return *m;
}

Matrix &operator+(Matrix &a, const double b) // addition operator matrix and scalar
// this has to be a non-member function of class Matrix
{
  int i, j;
  double *pa, *pm;
  Matrix *m;

  pa = a.data;                                  // assign pointers to data and do the element-by element addition
  m = new Matrix(a.rows, a.cols, a.data);
  pm = m->Data();
    for (i = 0; i < a.rows; i++) for (j = 0; j < a.cols; j++){ *pm = *pa + b; pm++; pa++;}
  return *m;
}

Matrix &operator-(Matrix &a, const double b) // substraction operator matrix and scalar
// this has to be a non-member function of class Matrix
{
  int i, j;
  double *pa, *pm;
  Matrix *m;

  pa = a.data;                                  // assign pointers to data and do the element-by element addition
  m = new Matrix(a.rows, a.cols, a.data);
  pm = m->Data();
  for (i = 0; i < a.rows; i++) for (j = 0; j < a.cols; j++) { *pm = *pa - b; pm++; pa++;}
  return *m;
}

Matrix &operator*(Matrix &a, const double b) // multiplication operator matrix and scalar
// this has to be a non-member function of class Matrix
{
  int i, j;
  double *pa, *pm;
  Matrix *m;

  pa = a.data;                                  // assign pointers to data and do the element-by element addition
  m = new Matrix(a.rows, a.cols, a.data);
  pm = m->Data();
  for (i = 0; i < a.rows; i++) for (j = 0; j < a.cols; j++) { *pm = *pa * b; pm++; pa++;}
  return *m;
}

Matrix &operator/(Matrix &a, const double b) // division operator matrix and scalar
// this has to be a non-member function of class Matrix
{
  int i, j;
  double *pa, *pm;
  Matrix *m;

  pa = a.data;                                  // assign pointers to data and do the element-by element addition
  m = new Matrix(a.rows, a.cols, a.data);
  pm = m->Data();
  for (i = 0; i < a.rows; i++) for (j = 0; j < a.cols; j++) { *pm = *pa / b; pm++; pa++;}
  return *m;
}

Matrix &operator-(const double a, Matrix &b)         // substraction operator scalar and matrix
{
  int i, j;
  double *pb, *pm;
  Matrix *m;

  pb = b.data;                                  // assign pointers to data and do the element-by element addition
  m = new Matrix(b.rows, b.cols, b.data);
  pm = m->Data();
    for (i = 0; i < b.rows; i++) for (j = 0; j < b.cols; j++) {*pm = a - *pb; pm++; pb++;}
  return *m;
}


Matrix &operator/(const double a, Matrix &b)         // element-by-element division operator scalar and matrix
{
  int i, j;
  double *pb, *pm;
  Matrix *m;

  pb = b.data;                                  // assign pointers to data and do the element-by element addition
  m = new Matrix(b.rows, b.cols, b.data);
  pm = m->Data();
    for (i = 0; i < b.rows; i++) for (j = 0; j < b.cols; j++) {*pm = a / *pb; pm++; pb++;}
  return *m;
}


Matrix &Matrix::Row(const int r)             // returns row of a matrix as a new matrix
{
  double *d;
  Matrix *m;

  if ((r > rows) || (r < 1))                // check out of range
  {
    cout << ERROR2 << rows << ERROR7 << r << endl;
    exit(EXIT_FAILURE);
  }
  d = data + (r - 1) * cols;                 // put data pointer on right location
  m = new Matrix(1, cols, d);                // make new matrix and insert the row rightaway
  return *m;
}


Matrix &Matrix::Col(const int c)               // returns column of a matrix as a new matrix
{
  double *d,*n;
  Matrix *m;
  int i;

  if ((c > cols) || (c < 1))                // check out of range
  {
    cout << ERROR2 << cols << ERROR7 << c <<endl;
    exit(EXIT_FAILURE);
  }
  d = data + c - 1;                         // put data pointer on right location
  m = new Matrix(rows, 1);
  n = m->Data();
  for (i = 0; i < rows; i++)
  {
     *n++ = *d;
     d+= cols;
  }
  return *m;
}


Matrix &Matrix::Range(const int r1, const int r2, const int c1, const int c2)        // returns a subrange as a new matrix
{
  Matrix *m;
  double *d, *n;
  int i;

  if ((c1 > cols) || (c1 < 1) || (c2 > cols) || (c2 < 1) || (r1 > rows) || (r1 < 1) || (r2 > rows) || (r2 < 1))
// check out of range
  {
    cout << ERROR2 << endl;
    exit(EXIT_FAILURE);
  }
  m = new Matrix(r2 - r1 + 1, c2 - c1 + 1);
  n = m->Data();
  d = data + (r1 - 1) * cols + c1 - 1;
  for (i = 0; i < m->rows; i++)
  {
    memcpy(n, d, (c2 - c1 + 1) * sizeof(double));   // copy rows
    d = d + cols;
    n = n + (c2 - c1 + 1);
  }
  return *m;
}

Matrix &Matrix::Range(const Matrix &m, const int r1, const int r2, const int c1, const int c2)        // Copies a subrange of m into current matrix
{

  double *d, *n;
  int i;

  if ((c1 > m.cols) || (c1 < 1) || (c2 > m.cols) || (c2 < 1) || (r1 > m.rows) || (r1 < 1) || (r2 > m.rows) || (r2 < 1))
// check out of range
  {
    cout << ERROR2 << endl;
    exit(EXIT_FAILURE);
  }
  Resize(r2 - r1 + 1, c2 - c1 + 1);
  n = m.data + (r1 - 1) * m.cols + c1 - 1;
  d = data;
  for (i = 0; i < rows; i++)
  {
    memcpy(d, n, cols * sizeof(double));   // copy rows
    d = d + cols;
    n = n + m.cols;
  }
  return *this;
}

Matrix &Matrix::Fliplr(const Matrix &m)                      // Flip matrix from left to right and return result in new matrix
{
  double *d,*n;
  int i, j;

  Resize(m.rows, m.cols);
  n = data;
  d = m.data;
  for (i = 0; i < rows; i++)
  {
      for (j = (cols - 1); j >= 0; j--) {*n = *(d + j); n++;}
    d+= cols;
  }
  return *this;
}


Matrix &Matrix::Flipud(const Matrix &m)                        // Flip matrix up-down  and return result in new matrix
{
  double *d,*n;
  int i;

  Resize(m.rows, m.cols);
  n = data;
  d = m.data + (rows - 1) * cols;
  for (i = (rows - 1); i >= 0; i--)
  {
    memcpy(n, d, cols * sizeof(double));
    d -= cols;
    n += cols;
  }
  return *this;
}


Matrix &Matrix::Transp(const Matrix &m)                               // Transpose matrix and return the result as a new matrix
{
  double *d,*n;
  int i, j;

  Resize(m.cols, m.rows);
  n = data;
  for (i = 1; i <= m.cols; i++)
  {
    d = m.data + (i - 1);
    for (j = 1; j <= m.rows; j++)
    {
      *n = *d;
        n++;
      d += m.cols;
    }
  }
  return *this;
}


Matrix &Matrix::Fliplr()                      // Flip matrix from left to right
{
  double *d, *n, *p;
  int i, j;

  n = new double[rows * cols];
  p = n;
  d = data;
  for (i = 0; i < rows; i++)
  {
      for (j = (cols - 1); j >= 0; j--) {*n = *(d + j); n++;}
    d+= cols;
  }
  memcpy(data, p, rows * cols * sizeof(double));
  delete p;
  return *this;
}


Matrix &Matrix::Flipud()                        // Flip matrix up-down
{
  double *d, *n, *p;
  int i;


  n = new double[rows * cols];
  p = n;
  d = data + (rows - 1) * cols;
  for (i = (rows - 1); i >= 0; i--)
  {
    memcpy(n, d, cols * sizeof(double));
    d -= cols;
    n += cols;
  }
  memcpy(data, p, rows * cols * sizeof(double));
  delete p;
  return *this;
}



Matrix &Matrix::Transp()                               // Transpose matrix
{
  double *d,*n, *p;
  int i, j;

  n = new double[rows * cols];
  p = n;                  // copy columns into rows
  for (i = 1; i <= cols; i++)
  {
    d = data + (i - 1);
    for (j = 1; j <= rows; j++)
    {
        *n = *d;
        n++;
        d += cols;
    }
  }
  memcpy(data, p, rows * cols * sizeof(double));
  i = rows;               // swap rows and columns
  rows = cols;
  cols = i;
  delete p;
  return *this;
}

Matrix &Matrix::Mult(const Matrix &m, const Matrix &n)
// Matrix multiplicaton of m and n, result is stored in current matrix
{
  int i, j, k;
  double *d;

  if (m.cols != n.rows)
  {
    cout << ERROR6 << endl;
    exit(EXIT_FAILURE);
  }
  Resize(m.rows, n.cols);
  d = data;
  for (i = 1; i <= m.rows; i++)
  {
    for (j = 1; j <= n.cols; j++)
    {
      for (k = 1; k <= m.cols; k++) *d += *(m.data + (i - 1) * m.cols + k - 1) * *(n.data + (k - 1) * n.cols + j - 1);
      d++;
    }
  }
  return *this;
}

void Matrix::Resize(int c)                   // resizes matrix to row vector
{
  int i;
  double *p;

  if (c < 1) cout << ERROR8 << endl;
  if (data != NULL) delete data;
  data = new double[c];
  if (data == NULL)
  {
    cout << ERROR1 << endl;
    exit(EXIT_FAILURE);
  }
  cols = c;
  rows = 1;
  p = data;
    for (i = 0; i < c; i++) {*p = 0.0; p++;}
  if (cols == 1) is_scalar = TRUE; else is_scalar = FALSE;
}

void Matrix::Resize(int r, int c)             // resizes to matrix with r rows and c columns
{
  int i;
  double *p;

  if ((c < 1) || (r < 1)) cout << ERROR8 << endl;
  if (data != NULL) delete data;
  data = new double[r*c];
  if (data == NULL)
  {
    cout << ERROR1 << endl;
    exit(EXIT_FAILURE);
  }
  rows = r;
  cols = c;
  p = data;
    for (i = 0; i < (r * c); i++) {*p = 0.0; p++;}
  if ((cols == 1) && (rows == 1)) is_scalar = TRUE; else is_scalar = FALSE;
}


void Matrix::Resize(int r, double *d)          // resizes to vector with r elements and read new data
{
  if (r < 1) cout << ERROR8 << endl;
  if (data != NULL) delete data;
  data = new double[r];
  if (data == NULL)
  {
    cout << ERROR8 << endl;
    exit(EXIT_FAILURE);
  }
  cols = r;
  rows = 1;
  memcpy(data, d, cols * sizeof(double));
  if ((cols == 1) && (rows == 1)) is_scalar = TRUE; else is_scalar = FALSE;
}

void Matrix::Resize(int r, int c, double *d)    // resizes to matrix with r rows and c columns and read new data
{
  if (data != NULL) delete data;

  if ((c < 1) || (r < 1)) cout << ERROR8 << endl;
  data = new double[r*c];
  if (data == NULL)
  {
    cout << ERROR1 << endl;
    exit(EXIT_FAILURE);
  }
  rows = r;
  cols = c;
  memcpy(data, d, rows * cols * sizeof(double));
  if ((cols == 1) && (rows == 1)) is_scalar = TRUE; else is_scalar = FALSE;
}


void Matrix::Insert(const int index, const double d)
// inserts new value d after position index in a vector
{
  double *tmp, *p, *t;

  if (rows > 1)
  {
    cout << ERROR3 << endl;
    exit(EXIT_FAILURE);
  }
  tmp = new double[cols];             // temporary array
  if (tmp == NULL)
  {
    cout << ERROR1 << endl;
    exit(EXIT_FAILURE);
  }
  t = tmp;
  p = data;
  memcpy(t, p, cols*sizeof(double));   // copy data;
  delete data;
  data = new double[cols + 1];         // resize the data
  if (data == NULL)
  {
    cout << ERROR1 << endl;
    exit(EXIT_FAILURE);
  }
  t = tmp;
  p = data;
  memcpy(p, t, index * sizeof(double));
  *(p + index) = d;
  memcpy(p + index +1, t + index, (cols - index) * sizeof(double));
  cols++;
  delete tmp;
  is_scalar = FALSE;
}


void Matrix::Disp()                           // displays matrix on screen
{
  int i, j;
  double *d;

  d = data;
  for (i = 0; i < rows; i++)
  {
    for (j = 0; j < cols; j++)
    {
        cout << *d << "  ";
        d++;
    }
    cout << endl;
  }
  cout << endl;
}


void Matrix::Write(ofstream *outfile)       //writes content of the matrix to output file stream
{
  int i, j;
  double *d;

  outfile->setf(ios::scientific, ios::floatfield);  // set scientific format
  d = data;
  for (i = 0; i < rows; i++)
  {
    for (j = 0; j < cols; j++)

    {
        *outfile << *d << " ";
        d++;
    }
    *outfile << endl;
  }
}

int CheckSize(Matrix *a, Matrix *b, int *r, int *c)    // checks the size of the matrix for validity of arithmetic operations
{
  if (a->is_scalar)
  {
    *r = b->rows;
    *c = b->cols;
  } else
  {
    if (b->is_scalar)
    {
      *r = a->rows;
      *c = a->cols;
    } else
    {
      if ((a->cols != b->cols) || (a->rows != b->rows))
      {
        cout << ERROR4 << endl;
        return FALSE;
      } else
      {
        *r = a->rows;
        *c = a->cols;
      }
    }
  }
  return TRUE;
}

double Matrix::Min()                                               // returns minimum value
{
  double m = DBL_MAX;
  double *d;
  int i, j;

  d = data;
  for (i = 0; i < rows; i++)
  {
    for (j = 0; j < cols; j++)
    {
      if (*d < m) m = *d;
      d++;
    }
  }
  return m;
}

double Matrix::Max()                                               // returns maximum value
{
  double m = -DBL_MAX;
  double *d;
  int i, j;

  d = data;
  for (i = 0; i < rows; i++)
  {
    for (j = 0; j < cols; j++)
    {
      if (*d > m) m = *d;
      d++;
    }
  }
  return m;
}


double Matrix::Sum()                                               // returns grand total of the matrix
{
  double s = 0;
  double *d;
  int i, j;

  d = data;
  for (i = 0; i < rows; i++)

  {
    for (j = 0; j < cols; j++)
    {
      s += *d;
      d++;
    }
  }
  return s;
}

double Matrix::Mean()                                               // returns grand total of the matrix
{
  double s = 0;
  double *d;
  int i, j;

  d = data;
  for (i = 0; i < rows; i++)

  {
    for (j = 0; j < cols; j++)
    {
      s += *d;
      d++;
    }
  }
  return (s / (rows * cols));
}


double Matrix::SumRow(int r)                                       // returns sum of selected row
{
  double s = 0;
  double *d;
  int i;

  if ((r > rows) || (r < 1))
  {
    cout << ERROR2 << endl;
    exit(EXIT_FAILURE);
  }
  d = data + (r - 1) * cols;
    for (i = 0; i < cols; i++) {s += *d; d++;}
  return s;
}

double Matrix::SumCol(int c)                                       // returns sum of selceted column
{
  double s = 0;
  double *d;
  int i;

  if ((c > cols) || (c < 1))
  {
    cout << ERROR2 << endl;
    exit(EXIT_FAILURE);
  }
  d = data + c - 1;
  for (i = 0; i < rows; i++)
  {
      s += *d;
      d+= cols;
  }
  return s;
}



void Matrix::PutData(int i, int nr, double *d)    // index operator to put more nr data from array d at index i
{

  if (((i + nr - 1) > (rows * cols)) || (i < 1))
  {
    cout << ERROR2 << endl;
    exit(EXIT_FAILURE);
  }
  memcpy(data + i - 1, d, nr * sizeof(double));
}

void Matrix::PutData(const int r, const int c, Matrix &m)              // copies all elements of Matrix m into destination matrix starting at position row r and column c
{
  int l;

  l = m.Length();
  if ((r > rows) || (c > cols) || (r < 1) || (c < 1) || ((r - 1) * cols + c - 1 + l > rows * cols))
  {
    cout << ERROR5 << endl;
    exit(EXIT_FAILURE);
  }
  memcpy(data + cols * (r - 1) + c - 1, m.data, l * sizeof(double));
}


void Matrix::PutData(const int r, const int c, Matrix &m, const int sr, const int sc, const int n)
// copies n elements of Matrix m, starting at row sr and column sc into destination matrix starting at position row r and column c
{
  int mcols, mrows;

  mcols = m.Cols();
  mrows = m.Rows();
  if ((r > rows) || (c > cols) || (r < 1) || (c < 1) || (sr > mrows) || (sc > mcols) || (sr < 1) || (sc < 1) || ((r - 1) * cols + c - 1 + n > rows * cols) || ((sr - 1) * mcols + sc - 1 + n > mrows * mcols))
  {
    cout << ERROR5 << endl;
    exit(EXIT_FAILURE);
  }
  memcpy(data + cols * (r - 1) + c - 1, m.data + mcols * (sr - 1) + sc - 1, n * sizeof(double));
}


Matrix::~Matrix()  // destructor
{
  if (data != NULL) delete data;
}

void Matrix::Fill(double v)                    // fills the matrix with value v;
{
  double *d;
  int i, j;

  d = data;
    for (i = 0; i < rows; i++) for (j = 0; j < cols; j++) {*d = v; d++;}
}


void Matrix::Fill(double first, double step, double last)
// fills the matrix with values starting with first, stepwise with step
// filling is row-wise, 
// a warning is generated when the number of fill values do not match the size of the matrix
{
  double *d;
  double v;
  int i, j, nsteps, s, stepcount;

  nsteps = ceil(abs(last - first) / step) + 1;  // number of steps error checking
  s = rows * cols;
  if (nsteps != s) cout << "Matrix fill: sizes do not match, matrix entries: " << s << ", nr of fill values: " << nsteps << endl;
  d = data;								// initialization
  v = first;
  stepcount = 0;
  for (i = 0; i < rows; i++)
  {
	for (j = 0; j < cols; j++)
	{
		stepcount += 1;					// fill; stop when the size of the matrix is exceeded
		if (stepcount > s) break;
		if (stepcount == nsteps) v = last;  // the last value is always equal to the stop value, whatever the step size is
		*d = v;
        d++;
		v += step;
		if (stepcount == nsteps) break;   // also stop when the step counter is exceeded
	}
	if (stepcount > s) break;
	if (stepcount == nsteps) break;
  }
}


int Matrix::Contains(double v)                    // checks whether v occurs in the matrix
{
  int yes = FALSE;
  double *d;
  int i, j;

  d = data;
  for (i = 0; i < rows; i++)
  {
    for (j = 0; j < cols; j++)
    {
      if (*d == v) yes = TRUE;
      d++;
    }
  }
  return yes;
}
