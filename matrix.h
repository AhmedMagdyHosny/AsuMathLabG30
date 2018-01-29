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
#include <ctime>
#include <vector>
#include <math.h>
using namespace std;


class CMatrix {
	int nR , nC ;
	double ** values ;
	string name;

	static bool isNumber(string s);

public :
	enum MI { MI_ZEROS, MI_ONES, MI_EYE, MI_RAND, MI_VALUE };
	
	//=================================
	//Constructors and destructors
	CMatrix () ;
	CMatrix (int nR , int nC) ;
	CMatrix(int nR, int nC, int initialization, double initializationValue=0);
	CMatrix(const CMatrix& m);
	CMatrix(string s);
	CMatrix(double d);
	void reset ();
	~CMatrix () ;
	
	//=================================
	//Functions
	void copy(const CMatrix& m);
	void copy(double d);
	void copy(string s);
	void addColumn(CMatrix& m);
	void addRow(CMatrix& m);
	void setSubMatrix(int iR, int iC, CMatrix& m);
	void add(const CMatrix& m);
	void sub(const CMatrix& m);
	void mul(const CMatrix& m);
	void dotSlash(int num);
	CMatrix power(double p);
	static string MatrixCat(string s);
	static CMatrix trigStrtoDouble(string trigExpression ,vector<CMatrix> &matrix_List);
	static CMatrix matrixInit(string expression);
	void set_name(string n);
	string get_name();
	string get_string();
	string getString();
	void PrintMatrix();
	int get_nR();
	int get_nC();
	void transpose(const CMatrix &m);
	CMatrix getCofactor(int r, int c);
	double getDeterminant();
	CMatrix getInverse();
	CMatrix getAdjoint(const CMatrix &m);
	//=================================
	// Trigonometric functions (sin ,cos ,tan ,sinh ,cosh ,tanh ,sin^-1 ,cos^-1 ,tan^-1 ,sinh^-1 ,cosh^-1 ,tanh^-1);
    static CMatrix matsin(CMatrix &M);  //sin
    static CMatrix matcos(CMatrix &M);  //cos
    static CMatrix mattan(CMatrix &M);  //tan
    static CMatrix matsinInv(CMatrix &M);  //asin2
    static CMatrix matcosInv(CMatrix &M);  //acos2
    static CMatrix mattanInv(CMatrix &M);  //atan2
    static CMatrix matcosh(CMatrix &M);    //cosh
    static CMatrix matsinh(CMatrix &M);    //sinh
    static CMatrix mattanh(CMatrix &M);    //tanh
    /*static CMatrix matcoshInv(CMatrix &M);   //acosh
    static CMatrix matsinhInv(CMatrix &M);   //asinh
    static CMatrix mattanhInv(CMatrix &M);   //atanh*/
	
	//=================================
	//operators
	CMatrix operator=(const CMatrix& m);
	CMatrix operator=(double d);
	void operator+=(const CMatrix& m);
	CMatrix operator+(const CMatrix& m);
	CMatrix operator+(double d);
	void operator+=(double d);
	CMatrix operator-(const CMatrix& m);
	void operator-=(const CMatrix& m);
	CMatrix operator-(double d);
	void operator-=(double d);
	CMatrix operator*(const CMatrix& m);
	void operator*=(const CMatrix& m);
	CMatrix operator/(const CMatrix& m);
	void operator/=(const CMatrix& m);
};

#endif
