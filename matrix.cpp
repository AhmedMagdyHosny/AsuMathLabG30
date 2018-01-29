#include "matrix.h"

string CMatrix::getString() {
	string s;
	for (int iR = 0; iR<nR; iR++) {
		for (int iC = 0; iC<nC; iC++) {
			char buffer[50];
			sprintf_s(buffer, 50, "%g ", values[iR][iC]);
			s += buffer;
		}
		if (iR != nR - 1)s += ";";
	}
	return s;
}

string CMatrix::MatrixCat(string s) {
	string buffer = s;
	int count = 0;
	buffer = buffer.substr(buffer.find("[") + 1, buffer.length() - buffer.find("[") - 2);
	string temp;
	temp = buffer.substr(buffer.find("["), buffer.rfind("]") - buffer.find("[") + 1);

	char* input = new char[temp.length() + 1];

	strcpy(input, temp.c_str());
	char*token;

	token = strtok(input, "],");
	CMatrix A;

	while (token) {

		CMatrix B(token);

		A.addColumn(B);
		if (A.nR != B.nR) throw("invalid concatenation");
		token = strtok(NULL, "],");

	}
	
	string output = A.getString();
	//cout << output << endl;
	s.replace(s.find(temp), temp.length(), output);

	return string(s);
}

CMatrix  CMatrix::getInverse(){
	int i,j,k,n;
	double a[100][200],t;
	n=nR ;
	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
			a[i][j]=values[i][j];
	for(i=0;i<n;i++)
	{
		for(j=n;j<2*n;j++)
		{
			if(i==j-n)
				a[i][j]=1;
			else
				a[i][j]=0;
		}
	}
	for(i=0;i<n;i++)
	{
		t=a[i][i];
		for(j=i;j<2*n;j++){
			if (t == 0){
				throw((string)"Singular matrix, can't find its inverse");
			}
			a[i][j]=a[i][j]/t;
		}
		for(j=0;j<n;j++)
		{
			if(i!=j)
			{
				t=a[j][i];
				for(k=0;k<2*n;k++)
					a[j][k]=a[j][k]-t*a[i][k];
			}
		}
	}
	double b [100][100];
	for(i=0;i<n;i++)
	{
		for(j=n;j<2*n;j++)
			b[i][j-n]= a[i][j];

	}
	nR=nC=n ;
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++){values[i][j]=b[i][j];}
	}
	return *this ;
}

CMatrix CMatrix::getAdjoint(const CMatrix &m){

	CMatrix temp, adj(m.nC, m.nR),cof;
	if(m.nC * m.nR == 1){
		temp = m;
		return temp;
	}
	int sign = 1;
	temp = m;
	for(int i=0; i<m.nC; i++){
		for(int j=0; j<m.nC; j++){
			cof = temp.getCofactor(i, j);
			sign = ((i+j)%2==0)? 1: -1;
			adj.values[j][i] = (sign)*(cof.getDeterminant());
		}
	}
	return adj;
}

CMatrix CMatrix::getCofactor(int r, int c)
{
	if(nR<=1 && nC<=1)
		throw("Invalid matrix dimension");
	CMatrix m(nR-1, nC-1);
	for(int iR=0;iR<m.nR;iR++)
		for(int iC=0;iC<m.nC;iC++)
		{
			int sR = (iR<r)?iR:iR+1;
			int sC = (iC<c)?iC:iC+1;
			m.values[iR][iC] = values[sR][sC];
		}
		return m;
}
double CMatrix::getDeterminant()
{
	if(nR!=nC)
		throw("Invalid matrix dimension");
	if(nR==1 && nC==1)return values[0][0];
	double value = 0, m = 1;
	for(int iR=0;iR<nR;iR++)
	{
		value+= m * values[0][iR] * getCofactor(0, iR).getDeterminant();
		m *= -1;
	}
	return value;
}

CMatrix CMatrix::operator/(const CMatrix& m){
	CMatrix r = *this;
	try{
		r /= m;
	}catch(string error){
		throw(error);
	}
	return r;
}

void CMatrix::operator/=(const CMatrix& m){
	CMatrix r = *this;
	CMatrix temp = m;
	try{
		r *= temp.getInverse();
	}catch(string error){
		throw(error);
		return;
	}
	*this = r;
}

void CMatrix::sub(const CMatrix& m){
	if(nR!=m.nR||nC!=m.nC)
		throw("Invalid matrix dimension");
	for(int iR=0;iR<nR;iR++)
		for(int iC=0;iC<nC;iC++)
			values[iR][iC] -= m.values[iR][iC];
}

void CMatrix::operator-=(const CMatrix& m){
	sub(m);
}
CMatrix CMatrix::operator-(const CMatrix& m){
	CMatrix r = *this;
	r-=m;
	return r;
}

void CMatrix::mul(const CMatrix& m){
	if(nC != m.nR)
		throw("Invalid matrix dimension");
	CMatrix r(nR, m.nC);
	for(int iR=0;iR<r.nR;iR++)
		for(int iC=0;iC<r.nC;iC++){
			r.values[iR][iC] = 0;
			for(int k=0;k<nC;k++)
				r.values[iR][iC] += values[iR][k]*m.values[k][iC];
		}
		copy(r);
}

void CMatrix::operator*=(const CMatrix& m){
	mul(m);
}

CMatrix CMatrix::operator*(const CMatrix& m){
	CMatrix r = *this;
	r*=m;
	return r;
}

CMatrix::CMatrix(const CMatrix& m){
	nR = nC = 0;
	values = NULL;
	copy(m);
}

void CMatrix::set_name(string n){
	this->name = n;
}

string CMatrix::get_name(){
	return this->name;
}

void CMatrix::transpose(const CMatrix &m){

	reset();
	this->nC = m.nR;
	this->nR = m.nC;
	if((nR*nC)==0){values=NULL; return;}
	values = new double*[nR];
	for(int iR=0;iR<nR;iR++)
	{
		values[iR] = new double[nC];
		for(int iC=0;iC<nC;iC++)
		{
			values[iR][iC] = m.values[iC][iR];
		}
	}

}

void CMatrix::dotSlash(int num){

	for (int i=0; i<this->nC; i++){
		for(int j=0; j<this->nR; j++){
			this->values[i][j] = num / this->values[i][j];
		}
	}

}

CMatrix CMatrix::operator+(const CMatrix &m)
{
	CMatrix r = *this;
	r+=m;
	return r;
}

void CMatrix::operator+=(const CMatrix &m)
{
	add(m);
}

void CMatrix::add(const CMatrix &m)
{
	if(nR!=m.nR||nC!=m.nC)
		throw("Invalid matrix dimension");
	for(int iR=0;iR<nR;iR++)
		for(int iC=0;iC<nC;iC++)
CMatrix::CMatrix(int nR, int nC, int initialization , double initializationValue) { 
	this->nR = nR;
	this->nC = nC; 
	if ((nR*nC) == 0) 
	{ values = NULL; return; }
	values = new double*[nR];
	for (int iR = 0; iR<nR; iR++) 
	{ 
		values[iR] = new double[nC]; 

		switch (initialization)
		{

		case MI_ZEROS:
			{for (int iC = 0; iC < nC; iC++)
				values[iR][iC] = 0; break; }
		case MI_ONES: 
			{for (int iC = 0; iC < nC; iC++)
				values[iR][iC] = 1; break; }
		case MI_EYE: 
			{ if(nC!=nR) throw("Invalid matrix dimension");

			for (int iC = 0; iC < nC; iC++)
				values[iR][iC] = (iR == iC) ? 1 : 0; break; }
		case MI_RAND: 
			{
				srand(time(0));
				for (int iC = 0; iC < nC; iC++)
					values[iR][iC] = (rand() % 1000000) / 1000.0; break; }
		case MI_VALUE:
			{for (int iC = 0; iC < nC; iC++)
				values[iR][iC] = initializationValue; break; }

		} 
	} 
}			values[iR][iC] += m.values[iR][iC];
}

CMatrix CMatrix::operator=(const CMatrix &m)
{
	copy(m);
	return *this;
}

CMatrix CMatrix::operator= (double d)
{
	copy(d);
	return *this;
}

CMatrix::CMatrix(string s){
	nR = nC = 0;
	values = NULL;
	copy(s);
}
CMatrix::CMatrix(double d){
	nR = nC = 0;
	values = NULL;
	copy(d);
}

CMatrix::CMatrix() {
	nR=nC=0 ;
	values = NULL ;
}

CMatrix::CMatrix(int nR , int nC) {
	this-> nR=nR ;
	this-> nC=nC;
	if((nR*nC)==0){values=NULL; return;}
	values = new double * [nR] ;
	//va_list va ;
	//va_start (va , nC);
	for(int i = 0 ; i<nR ; i++){
		values[i]= new double [nC];
		/*for(int j=0 ; j<nC ; j++){
		values[i][j] = va_arg(va,double) ;
		}*/
	}
	//va_end(va) ;
}

void CMatrix::reset(){
	if(values){
		for(int i=0;i<nR;i++)
			delete[] values[i];
		delete[] values;
	}
	nR = nC = 0;
	values = NULL;}

CMatrix::~CMatrix(){
	this->reset();
}


void CMatrix::copy(const CMatrix& m)
{
	reset();
	this->name = m.name;
	this->nR = m.nR;
	this->nC = m.nC;
	if((nR*nC)==0){values=NULL; return;}
	values = new double*[nR];
	for(int iR=0;iR<nR;iR++)
	{
		values[iR] = new double[nC];
		for(int iC=0;iC<nC;iC++)
		{
			values[iR][iC] = m.values[iR][iC];
		}
	}
}


void CMatrix::copy(string s){
	reset();
	char* buffer = new char[s.length()+1];
	strcpy(buffer, s.c_str());
	char* lineContext;
	const char* lineSeparators = ";\r\n";
	char* line = strtok_r(buffer, lineSeparators, &lineContext);
	while(line){
		CMatrix row;
		char* context;
		const char* separators = " []";
		char* token = strtok_r(line, separators, &context);
		while(token){
			CMatrix item = atof(token);
			row.addColumn(item);
			token = strtok_r(NULL, separators, &context);
		}
		if(row.nC>0 && (row.nC==nC || nR==0))addRow(row);
		line = strtok_r(NULL, lineSeparators, &lineContext);
	}
	delete[] buffer;
}



void CMatrix::copy(double d){
	reset();
	this->nR = 1;
	this->nC = 1;
	values = new double*[1];
	values[0] = new double[1];
	values[0][0] = d;}

void CMatrix::addColumn(CMatrix& m){
	CMatrix n(max(nR, m.nR), nC+m.nC);
	n.setSubMatrix(0, 0, *this);
	n.setSubMatrix(0, nC, m);
	*this = n;
}
void CMatrix::addRow(CMatrix& m){
	CMatrix n(nR+m.nR, max(nC, m.nC));
	n.setSubMatrix(0, 0, *this);
	n.setSubMatrix(nR, 0, m);
	*this = n;
}
void CMatrix::setSubMatrix(int r, int c, CMatrix& m){
	if((r+m.nR)>nR || (c+m.nC)>nC)
		throw("Invalid matrix dimension");
	for(int iR=0;iR<m.nR;iR++)
		for(int iC=0;iC<m.nC;iC++)
			values[r+iR][c+iC] = m.values[iR][iC];
}

void CMatrix::PrintMatrix()
{
	cout << this -> get_name() << " = [";
	for (int i = 0; i < nR; i++)
	{
		for (int j = 0; j < nC; j++){
			if (j < nC-1)
				cout << values[i][j] << " ";
			else
				cout << values[i][j];
		}
		if (i < nR-1)
			cout << "; ";
	}
	cout << "]" << endl;
}


// Trig functions....
CMatrix CMatrix::matsin(CMatrix &M){
	CMatrix Result(M.nR, M.nC);
	for (int i=0; i<Result.nR; i++) {
		for (int j=0; j<Result.nC; j++)
			Result.values[i][j] = sin(M.values[i][j]);
	}
	return Result;
}

CMatrix CMatrix::matcos(CMatrix &M){
	CMatrix Result(M.nR, M.nC);
	for (int i=0; i<Result.nR; i++) {
		for (int j=0; j<Result.nC; j++)
			Result.values[i][j] = cos(M.values[i][j]);
	}
	return Result;
}

CMatrix CMatrix::mattan(CMatrix &M){
	CMatrix Result(M.nR, M.nC);
	for (int i=0; i<Result.nR; i++) {
		for (int j=0; j<Result.nC; j++)
			Result.values[i][j] = tan(M.values[i][j]);
	}
	return Result;
}


CMatrix CMatrix::trigStrtoDouble(string trigExpression ,vector<CMatrix> &matrix_List){

	string trigFunction = trigExpression.substr(0,trigExpression.find('('));
	string angle = trigExpression.substr(trigExpression.find('(')+1 , trigExpression.find(')')-trigExpression.find('(')-1);

	if(isNumber(angle)){
		CMatrix Result(1,1);

		if (trigFunction == "sin") {
			Result.values[0][0] = sin(atof(angle.c_str()));
		}else if (trigFunction == "cos"){
			Result.values[0][0] = cos(atof(angle.c_str()));
		}else if (trigFunction == "tan"){
			Result.values[0][0] = tan(atof(angle.c_str()));
		}else if (trigFunction == "sinh"){
			Result.values[0][0] = sinh(atof(angle.c_str()));
		}else if (trigFunction == "cosh"){
			Result.values[0][0] = cosh(atof(angle.c_str()));
		}else if (trigFunction == "tanh"){
			Result.values[0][0] = tanh(atof(angle.c_str()));
		}else if (trigFunction == "sinInv") {
			Result.values[0][0] = asin(atof(angle.c_str()));
		}else if (trigFunction == "cosInv"){
			Result.values[0][0] = acos(atof(angle.c_str()));
		}else if (trigFunction == "tanInv"){
			Result.values[0][0] = atan(atof(angle.c_str()));
		}/*else if (trigFunction == "sinhInv"){
		 Result.values[0][0] = asinh(atof(angle.c_str()));
		 }else if (trigFunction == "coshInv"){
		 Result.values[0][0] = acosh(atof(angle.c_str()));
		 }else if (trigFunction == "tanhInv"){
		 Result.values[0][0] = atanh(atof(angle.c_str()));
		 }*/else if (trigFunction == "sqrt"){
			 Result.values[0][0] = pow(atof(angle.c_str()), 0.5);
		}
		return Result;

	}else{
		for (int i=0; i<matrix_List.size(); i++){
			if (matrix_List[i].get_name() == angle){
				CMatrix temp = matrix_List[i];
				CMatrix ResultMat;
				if (trigFunction == "sin") {
					ResultMat = matsin(temp);
				}else if (trigFunction == "cos"){
					ResultMat = matcos(temp);
				}else if (trigFunction == "tan"){
					ResultMat = mattan(temp);
				}else if (trigFunction == "sinh"){
					ResultMat = matsinh(temp);
				}else if (trigFunction == "cosh"){
					ResultMat = matcosh(temp);
				}else if (trigFunction == "tanh"){
					ResultMat = mattanh(temp);
				}else if (trigFunction == "sinInv") {
					ResultMat = matsinInv(temp);
				}else if (trigFunction == "cosInv"){
					ResultMat = matcosInv(temp);
				}else if (trigFunction == "tanInv"){
					ResultMat = mattanInv(temp);
				}/*else if (trigFunction == "sinhInv"){
				 ResultMat = matsinhInv(temp);
				 }else if (trigFunction == "coshInv"){
				 ResultMat = matcoshInv(temp);
				 }else if (trigFunction == "tanhInv"){
				 ResultMat = mattanhInv(temp);
				 }*/else if (trigFunction == "sqrt"){
					 ResultMat = temp;
					 ResultMat.power(0.5);
				}
				return ResultMat;
			}
		}
	}
	throw string("invalid Expression");
}


CMatrix CMatrix::power(double p) {

	CMatrix A(nR,nC);
	for (int ir = 0; ir < nR; ir++)
		for (int ic = 0; ic < nC; ic++)
			A.values[ir][ic] = pow(values[ir][ic], p);

	return A;
}


string CMatrix::get_string(){

	string matrixData = "";
	if(!(nR == 1 && nC == 1))
		matrixData += "[";
	for (int i = 0; i < nR; i++)
	{
		for (int j = 0; j < nC; j++){
			if (j < nC-1){
				matrixData += to_string(values[i][j]); 
				matrixData += " ";
			}
			else
				matrixData += to_string(values[i][j]);
		}
		if (i < nR-1)
			matrixData += "; ";
	}
	if(!(nR == 1 && nC == 1))
		matrixData += "]";
	return matrixData;

}


CMatrix CMatrix::matrixInit(string expression){

	string Function = expression.substr(0,expression.find('('));
	string param = expression.substr(expression.find('(')+1 , expression.find(')')-expression.find('(')-1);
	CMatrix ResultMat;

	if (Function == "rand"){
		int a,b,comma;
		string a_str, b_str;
		comma = param.find(",");
		a_str = param.substr(0,comma);
		if(param[comma+1] == ' ') comma++;
		b_str = param.substr(comma+1);
		a = atof(a_str.c_str());
		b = atof(b_str.c_str());
		ResultMat = CMatrix(a, b, CMatrix::MI_RAND);
	}else if (Function == "eye"){
		int a,b,comma;
		string a_str, b_str;
		comma = param.find(",");
		a_str = param.substr(0,comma);
		if(param[comma+1] == ' ') comma++;
		b_str = param.substr(comma+1);
		a = atof(a_str.c_str());
		b = atof(b_str.c_str());
		ResultMat = CMatrix(a, b, CMatrix::MI_EYE);
	}else if (Function == "zeros"){
		int a,b,comma;
		string a_str, b_str;
		comma = param.find(",");
		a_str = param.substr(0,comma);
		if(param[comma+1] == ' ') comma++;
		b_str = param.substr(comma+1);
		a = atof(a_str.c_str());
		b = atof(b_str.c_str());
		ResultMat = CMatrix(a, b, CMatrix::MI_ZEROS);
	}else if (Function == "ones"){
		int a,b,comma;
		string a_str, b_str;
		comma = param.find(",");
		a_str = param.substr(0,comma);
		if(param[comma+1] == ' ') comma++;
		b_str = param.substr(comma+1);
		a = atof(a_str.c_str());
		b = atof(b_str.c_str());
		ResultMat = CMatrix(a, b, CMatrix::MI_ONES);
	}
	return ResultMat;
}
