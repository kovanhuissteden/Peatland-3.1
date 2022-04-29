/***************************************************************************
       matrix.h  -  matrix class to facilitate translation from Matlab
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


/**
  *@author J. van Huissteden
  */

#include <iostream>
#include <fstream>
using namespace std;


class Matrix
{
public: 
	Matrix();
  Matrix(int r);                                              // initializes a vector with r elements
  Matrix(int r, int c);                                       // initializes matrix with r rows and c columns
  Matrix(int r, double *d);                                   // initializes avector with r elements and data specified in *d
  Matrix(int r, int c, double *d);                            // initializes matrix with r rows and c columns and data specified in *d
  Matrix(Matrix &m);                                          // copy constructor
  ~Matrix();                                                  // destructor
  double &operator()(int i);                                  // index operator
  double &operator()(int i, int j);                           // index operator row  column
  Matrix &operator=(const Matrix &m);                         // assignment operator
  Matrix &operator+=(const Matrix &m);                        // addition operator - if m is a1 x 1 matrix it is treated as a scalar
  Matrix &operator-=(const Matrix &m);                        // substraction operator - if m is a1 x 1 matrix it is treated as a scalar
  Matrix &operator+=(const double m);                         // addition of scalar
  Matrix &operator-=(const double m);                         // substraction of scalar
  Matrix &operator*=(const Matrix &m);                        // element by element multiplication operator - if m is a1 x 1 matrix it is treated as a scalar
  Matrix &operator*=(const double m);                         // multiplication by scalar
  Matrix &operator/=(const Matrix &m);                        // element by element division operator - if m is a1 x 1 matrix it is treated as a scalar
  Matrix &operator/=(const double m);                         // division by scalar
  Matrix &operator-();                                        // unary minus
  friend Matrix &operator+(Matrix &a, Matrix &b);             // addition operator 2 matrices
  friend Matrix &operator-(Matrix &a, Matrix &b);             // substraction operator 2 matrices
  friend Matrix &operator*(Matrix &a, Matrix &b);             // element-by-element multiplication operator 2 matrices
  friend Matrix &operator/(Matrix &a, Matrix &b);             // element-by-element division operator 2 matrices
  friend Matrix &operator+(Matrix &a, const double b);        // addition operator matrix and scalar
  friend Matrix &operator-(Matrix &a, const double b);        // substraction operator matrix and scalar
  friend Matrix &operator*(Matrix &a, const double b);        // multiplication operator matrix and scalar
  friend Matrix &operator/(Matrix &a, const double b);        // division operator matrix and scalar
  friend Matrix &operator-(const double a, Matrix &b);        // substraction operator scalar and matrix
  friend Matrix &operator/(const double a, Matrix &b);        // element-by-element division operator scalar and matrix
  Matrix &Row(const int r);                                   // returns row of a matrix as a new matrix
  Matrix &Col(const int c);                                   // returns column of a matrix as a new matrix
  Matrix &Range(const int r1, const int r2, const int c1, const int c3);  // returns a subrange as a new matrix
  Matrix &Range(const Matrix &m, const int r1, const int r2, const int c1, const int c2);        // Copies a subrange of m into current matrix
  Matrix &Fliplr(const Matrix &m);                            // Flip matrix m from left to right and copy in current matrix
  Matrix &Fliplr();                                           // Flip matrix from left to right
  Matrix &Flipud(const Matrix &m);                            // Flip matrix up down and  copy in current matrix
  Matrix &Flipud();                                           // Flip matrix up down
  Matrix &Transp(const Matrix &m);                            // Transpose matrix m and  and copy in current matrix
  Matrix &Transp();                                           // Transpose matrix
  Matrix &Mult(const Matrix &m, const Matrix &n);             // Matrix multiplicaton of m and n, result is stored in current matrix
  inline double *Data() {return (data);}                      // returns pointer to data
  inline int Rows() {return (rows);}                          // returns number of rows
  inline int Cols() {return (cols);}                          // returns number of columns
  inline int Length() {return (rows * cols);}                 // returns number of matrix entries
  inline int IsScalar() {return (is_scalar);}                 // returns TRUE if matrix is 1 * 1
  void Resize(int c);                                         // resizes matrix to vector
  void Resize(int r, int c);                                  // resizes matrix
  void Resize(int r, double *d);                              // resizes to vector and read new data
  void Resize(int r, int c, double *d);                       // resizes and read new data
  void Insert(const int index, const double d);               // inserts new value d after position index in a vector
  void Disp();                                                // displays content of the matrix on screen
  void Write(ofstream *outfile);                              // writes content of the matrix to output file stream
  double Min();                                               // returns minimum value
  double Max();                                               // returns maximum value
  double Mean();
  double Sum();                                               // returns grand total of the matrix
  double SumRow(int r);                                       // returns sum of selected row
  double SumCol(int c);                                       // returns sum of selceted column
  void PutData(int i, int nr, double *d);                     // index operator to put more nr data from array d at index i
  void PutData(const int r, const int c, Matrix &m);          // copies all elements of Matrix m into destination matrix starting at position row r and column c
  void PutData(const int r, const int c, Matrix &m, const int sr, const int sc, const int n);
  // copies n elements of Matrix m, starting at row sr and column sc into destination matrix starting at position row r and column c
  void Fill(double v);                                        // fills the matrix with value v
  void Fill(double first, double step, double last);          // fills the matrix with values starting with first, stepwise with step
															  // filling is row-wise, a warning is generated when the number of fill values do not match the size of the matrix
  int Contains(double v);                                     // checks whether v occurs in the matrix

private:
  double *data;                                               // pointer to data
  int rows;                                                   // number of rows
  int cols;                                                   // number of columns
  int is_scalar;                                              // TRUE for a 1 * 1 matrix
  friend int CheckSize(Matrix *a, Matrix *b, int *r, int *c); // checks the size of the matrix for validity of arithmetic operations
};
