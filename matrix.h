//=================================
// include guard
#ifndef __MATRIX_H_INCLUDED__
#define __MATRIX_H_INCLUDED__

//=================================
// included dependencies
#include <iostream>
#include <string>
#include <stdlib.h>
#include <string.h>
using namespace std;


class CMatrix {
	int nR , nC ;
	double ** values ;
	string name;

public :
	// constructors and destructors
	CMatrix () ;
	CMatrix (int nR , int nC) ;
	void reset ();
	~CMatrix () ;
	CMatrix(string s);
	CMatrix(double d);
	CMatrix(const CMatrix& m);
	//functions
	void copy(const CMatrix& m);
	void copy(double d);
	void copy(string s);
	void addColumn(CMatrix& m);
	void addRow(CMatrix& m);
	void setSubMatrix(int iR, int iC, CMatrix& m);
	void add(const CMatrix& m);
	void sub(const CMatrix& m);
	void mul(const CMatrix& m);
	void set_name(string n);
	string get_name();
	void PrintMatrix();
	void transpose(const CMatrix &m);
	void dotSlash(int num);
	CMatrix getCofactor(int r, int c);
	double getDeterminant();
	CMatrix getInverse();
	CMatrix getAdjoint(const CMatrix &m);

	//operators
	CMatrix operator=(const CMatrix& m);
	CMatrix operator=(double d);
	void operator+=(const CMatrix& m);
	CMatrix operator+(const CMatrix& m);
	void operator-=(const CMatrix& m);
	CMatrix operator-(const CMatrix& m);
	void operator*=(const CMatrix& m);
	CMatrix operator*(const CMatrix& m);
	CMatrix operator/(const CMatrix& m);
	void operator/=(const CMatrix& m);
	
	//Trig Functions
	CMatrix matsin(CMatrix &M);
	CMatrix matcos(CMatrix &M);
	CMatrix mattan(CMatrix &M);

};

#endif
