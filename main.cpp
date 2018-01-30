#include"matrix.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

vector<CMatrix> matrix_list;

//=================================
//Functions defenition
string check_command(string command);
void to_upperCase(string &s);
bool isNumber(string s);
bool isName(string s);
void formatStringSpaces(string &s);
string stringToken(string s, string separators, int &currentIndex);
void checkOperation(vector<string> &command_list);
string do_operation(string A, string op, string B);
double evaluate(double A, string op, double B);
CMatrix evaluate(CMatrix &A, string op, CMatrix &B);
CMatrix evaluate(CMatrix &A, string op, double B, bool doubleFirst);
string addMatrix(CMatrix &m);
void assignAndPrint(vector<string> &command_list);
string do_transpose(string A);



int main(int argc, char* argv[]){
	string inputMode = "console";
	string command;
	ifstream inputFile;
	if (argc == 2){
		inputFile.open(argv[1]);
		if (inputFile.bad() || inputFile.fail()){
			cout << "Failed to open input file" << endl;
			return 1;
		}
		inputMode = "file";
	}
	while(true){
		if (inputMode == "console"){
			getline(cin, command);
			if (command.find("\r") != -1){
				command.erase(command.find("\r"),2);
			}
			if (command == "exit"){
				break;
			}
			if(command == "")	continue;
			//formatStringSpaces(command);//Make spaces consistent
			try{
				check_command(command);
			}catch(string error){
				cout << error << endl;
			}
		}else if (inputMode == "file"){
			getline(inputFile,command);
			if (command.find("\r") != -1){
				command.erase(command.find("\r"),2);
			}
			if (command.find('[') != -1 && command.find(']') == -1){
				string temp;
				while (true){
					getline(inputFile, temp);
					if (temp.find("\r") != -1){
						temp.erase(temp.find("\r"),2);
					}
					command += temp;
					if (temp.find(']') != -1) break;
				}
			}
			if(command == "")	continue;
			//formatStringSpaces(command);//Make spaces consistent
			try{
				check_command(command);
			}catch(string error){
				cout << error << endl;
			}
			if (inputFile.eof()){
				inputFile.close();
				break;
			}
		}
	}
	matrix_list.clear();
	return 0;
}

string check_command(string command){
	string name, matrixData;
	int bracketIndex = command.find('[');
	if (bracketIndex != -1){
		if(command.find('[', bracketIndex+1) == -1){
			name = command.substr(0,command.find(" "));
			to_upperCase(name);
			matrixData = command.substr(bracketIndex);
			CMatrix temp(matrixData);
			temp.set_name(name);
			matrix_list.push_back(temp);
		}else{
			//=================================
			//Concatination
			vector<string> command_list;
			string separators = " [];,";
			int currentIndex = 0;
			string token = stringToken(command, separators, currentIndex);
			while(token != ""){
				command_list.push_back(token);
				token = stringToken(command, separators, currentIndex);
			}
			for(int i=2; i<command_list.size(); i++){
				if(isName(command_list[i])){
					name = command_list[i];
					CMatrix temp;
					for(int j=0; j<matrix_list.size(); j++){
						if(matrix_list[j].get_name() == command_list[i]){
							temp = matrix_list[j];
							break;
						}
					}
					int index = command.find(name);
					command.erase(index,name.length());
					command.insert(index, temp.get_string());
				}else if(!isNumber(command_list[i])){
					string result = check_command(command_list[i]);
					int index = command.find(command_list[i]);
					command.erase(index,command_list[i].length());
					command.insert(index, result);
				}
			}
			try{
				string matrixData = CMatrix::MatrixCat(command);
				CMatrix temp(matrixData.substr(matrixData.find("[")));
				string name = addMatrix(temp);
				command_list.erase(command_list.begin()+2,command_list.end());
				command_list.push_back(name);
				assignAndPrint(command_list);
			}catch (string error){
				throw(error);
			}
		}
	}else{
		formatStringSpaces(command);//Make spaces consistent
		vector<string> command_list;
		string separators = " ";
		int currentIndex = 0;
		string token = stringToken(command, separators, currentIndex);
		while(token != ""){
			command_list.push_back(token);
			token = stringToken(command, separators, currentIndex);
		}
		//=================================
		//Check if trasnpose is the only operation
		if(command_list.size() == 3 && command_list[2].find("'") != -1){
			string name = do_transpose(command_list[2]);
			command_list[2] = name;
		}
		try{
			checkOperation(command_list);
			if(command_list.size() == 1)	return command_list[0];
			assignAndPrint(command_list);
			command_list.clear();
		}catch(string error){
			command_list.clear();
			throw(error);
		}
	}
	return "";
}

void checkOperation(vector<string> &command_list){

	for (int i=0; i<command_list.size(); i++){
		if( command_list[i][0] == '('){
			command_list[i].erase(0,1);
			command_list[i].erase(command_list[i].end()-1);
			string result = check_command(command_list[i]);
			std::vector<string>::iterator it;
			it = command_list.insert(command_list.begin()+i, result);
			command_list.erase(it+1);
		}
	}

	for (int i=0; i<command_list.size(); i++){
		int bracket = command_list[i].find('(');
		if (bracket == -1)	continue;
		string trigfuncs = "sin cos tan sqrt";
		string initfuncs = "rand eye zeros ones";
		string funcName = command_list[i].substr(0,bracket);
		if (trigfuncs.find(funcName) != -1){
			string angle = command_list[i].substr(bracket+1, command_list[i].rfind(")")-bracket-1);
			if(!isNumber(angle) && !isName(angle)){
				string ang = check_command(angle);
				command_list[i].erase(bracket+1, angle.length());
				command_list[i].insert(bracket+1, ang);
			}
			CMatrix result =  CMatrix::trigStrtoDouble(command_list[i], matrix_list);
			string s = addMatrix(result);
			std::vector<string>::iterator it;
			it = command_list.insert(command_list.begin()+i, s);
			command_list.erase(it+1);
		}else if(initfuncs.find(funcName) != -1){
			CMatrix result =  CMatrix::matrixInit(command_list[i]);
			string s = addMatrix(result);
			std::vector<string>::iterator it;
			it = command_list.insert(command_list.begin()+i, s);
			command_list.erase(it+1);
		}

	}

	for (int i=0; i<command_list.size(); i++){
		if(command_list[i] == "^" || command_list[i] == ".^"){
			std::vector<string>::iterator it;
			string result;
			try{
				result = do_operation(command_list[i-1], command_list[i], command_list[i+1]);
			}catch(string error){
				throw(error);
			}
			it = command_list.insert(command_list.begin()+i-1, result);
			for(int j=0; j<3; j++) command_list.erase(it+1);
			i--;
		}
	}

	for (int i=0; i<command_list.size(); i++){
		if(command_list[i] == "*" || command_list[i] == "/" || command_list[i] == "./"){
			std::vector<string>::iterator it;
			string result;
			try{
				result = do_operation(command_list[i-1], command_list[i], command_list[i+1]);
			}catch(string error){
				throw(error);
			}
			it = command_list.insert(command_list.begin()+i-1, result);
			for(int j=0; j<3; j++) command_list.erase(it+1);
			i--;
		}
	}

	for (int i=0; i<command_list.size(); i++){
		if(command_list[i] == "+" || command_list[i] == "-" || command_list[i] == ".+" || command_list[i] == ".-"){
			std::vector<string>::iterator it;
			string result;
			try{
				result = do_operation(command_list[i-1], command_list[i], command_list[i+1]);
			}catch(string error){
				throw(error);
			}
			it = command_list.insert(command_list.begin()+i-1, result);
			for(int j=0; j<3; j++) command_list.erase(it+1);
			i--;
		}
	}

}

//=================================
//Convert a whole string into upper case.
void to_upperCase(string &s){

	for (int i=0, n=s.size(); i<n; i++){
		s[i] = toupper(s[i]);
	}
}

//=================================
//Check if a string is a number.
bool isNumber(string s){

	for(int i=0, n=s.length(); i<n; i++){
		if (!isdigit(s[i]) && s[i] != '.' && s[i] != '-')	return false;
	}
	return true;

}

//=================================
//Check if a string is a name.
bool isName(string s){

	for(int i=0, n=s.length(); i<n; i++){
		if (!isalpha(s[i]))	return false;
	}
	return true;

}


//=================================
//Make the spaces consistent in the input string
void formatStringSpaces(string &s){

	if (s == "")	return;
	string operators = "+-*/^=";
	if (s[s.size()-1] == ';'){
		s.erase(s.size()-1);
	}
	for(int i=0, n=s.length(); i<n; i++){
		int index = operators.find(s[i]);
		if (index != -1){
			if (s[i-1] == '.'){
				if (s[i-2] != ' '){
					s.insert(i-1, " ");
					n++;
					i++;
				}
				if (s[i+1] != ' '){
					s.insert(i+1, " ");
					n++;
					i++;
				}
			}else{
				if (s[i-1] != ' '){
					s.insert(i, " ");
					n++;
					i++;
				}
				if (s[i+1] != ' '){
					s.insert(i+1, " ");
					n++;
					i++;
				}
			}
		}
		if (s[i] == ','){
			if (s[i-1] == ' '){
				s.erase(i-1,1);
				n--;
				i--;
			}
			if (s[i+1] == ' '){
				s.erase(i+1,1);
				n--;
				i--;
			}
		}
	}

}

//=================================
//Splits the string into blocks, it returns an empty string when the end of the string is reached
//Works like strtok_s
string stringToken(string s, string separators, int &currentIndex){

	string toReturn;
	bool parnths = true, trig = true;;
	for(int i=currentIndex, n=s.length(); i<n; i++){
		if(s[i] == '(' && i != 0 && !isalpha(s[i-1]) && parnths){
			int parnthsclose = s.find(')',i);
			int parnthsopen = s.find('(', i+1);
			while (parnthsopen < parnthsclose && parnthsopen != -1){
				parnthsclose = s.find(')',parnthsclose+1);
				parnthsopen = s.find('(', parnthsopen+1);
			}
			toReturn = "(";
			toReturn += s.substr(i+1, parnthsclose-i-1);
			toReturn += ")";
			currentIndex = parnthsclose+1;
			return toReturn;
		}else if(s[i] == '(' && i == 0 && parnths){
			int parnthsclose = s.find(')',i);
			int parnthsopen = s.find('(', i+1);
			while (parnthsopen < parnthsclose && parnthsopen != -1){
				parnthsclose = s.find(')',parnthsclose+1);
				parnthsopen = s.find('(', parnthsopen+1);
			}
			toReturn = "(";
			toReturn += s.substr(i+1, parnthsclose-i-1);
			toReturn += ")";
			currentIndex = parnthsclose+1;
			return toReturn;
		}else if(s[i] == '(' && i != 0 && isalpha(s[i-1]) && trig){
			int parnthsclose = s.find(')',i);
			toReturn = "";
			toReturn += s.substr(currentIndex, parnthsclose-currentIndex+1);
			currentIndex = parnthsclose+1;
			return toReturn;
		}else{
			parnths = false;
			if(!isalpha(s[i]))	trig = false;
			if(separators.find(s[i]) != -1){
				if(i-currentIndex == 0){
					currentIndex = i+1;
					continue;
				}
				toReturn = s.substr(currentIndex, i-currentIndex);
				currentIndex = i+1;
				return toReturn;
			}
		}
		if(i == n-1){
			toReturn = s.substr(currentIndex);
			currentIndex = i+1;
			return toReturn;
		}
	}
	return (string)"";
}

string do_operation(string A, string op, string B){
	bool AIsNum, BIsNum;
	bool A_trans = false, B_trans = false;
	AIsNum = isNumber(A);
	BIsNum = isNumber(B);
	double A_double, B_double;
	CMatrix A_mat, B_mat;
	CMatrix temp;
	//=================================
	//Check for transpose
	if (A.find("'") != -1){
		A.erase(A.find("'"),1);
		A_trans = true;
	}
	if (B.find("'") != -1){
		B.erase(B.find("'"),1);
		B_trans = true;
	}

	//if(op.find('.') != -1){//one operand is a number, the other is a matrix
	if(AIsNum && !BIsNum){
		A_double = atof(A.c_str());
		for (int i=0; i<matrix_list.size(); i++){
			if (matrix_list[i].get_name() == B){
				if(B_trans)
					B_mat.transpose(matrix_list[i]);
				else
					B_mat = matrix_list[i];
				break;
			}
		}
		try{
			temp = evaluate(B_mat, op, A_double, AIsNum);
		}catch(string error){
			throw(error);
		}
		return addMatrix(temp);
	}else if (!AIsNum && BIsNum){
		B_double = atof(B.c_str());
		for (int i=0; i<matrix_list.size(); i++){
			if (matrix_list[i].get_name() == A){
				if(A_trans)
					A_mat.transpose(matrix_list[i]);
				else
					A_mat = matrix_list[i];
				break;
			}
		}
		try{
			temp = evaluate(A_mat, op, B_double, AIsNum);
		}catch(string error){
			throw(error);
		}
		return addMatrix(temp);
	}else if(AIsNum && BIsNum){
		A_double = atof(A.c_str());
		B_double = atof(B.c_str());
		return to_string(evaluate(A_double, op, B_double));
	}else{
		for (int i=0; i<matrix_list.size(); i++){
			if (matrix_list[i].get_name() == A){
				if(A_trans)
					A_mat.transpose(matrix_list[i]);
				else
					A_mat = matrix_list[i];
			}
			if (matrix_list[i].get_name() == B){
				if(B_trans)
					B_mat.transpose(matrix_list[i]);
				else
					B_mat = matrix_list[i];
			}
		}
		try{
			temp = evaluate(A_mat, op, B_mat);
		}catch(string error){
			throw(error);
		}
		return addMatrix(temp);
	}

}

string do_transpose(string A){

	int trans = A.find("'");
	CMatrix temp;
	if(trans != -1){
		A.erase(trans, 1);
		for (int i=0; i<matrix_list.size(); i++){
			if (matrix_list[i].get_name() == A){
				temp.transpose(matrix_list[i]);
			}
		}
		string name = addMatrix(temp);
		return name;
	}

}

double evaluate(double A, string op, double B){

	if (op == "+"){
		return A+B;
	}else if (op == "-"){
		return A-B;
	}else if (op == "*"){
		return A*B;
	}else if (op == "/"){
		try{
			return A/B;
		}catch(string err){
			throw(err);
		}
	}else if(op == "^"){
		return pow(A,B);
	}else{
		throw((string)"Wrong operation!");
	}

}

CMatrix evaluate(CMatrix &A, string op, CMatrix &B){

	if (op == "+"){
		return A+B;
	}else if (op == "-"){
		return A-B;
	}else if (op == "*"){
		return A*B;
	}else if (op == "/"){
		try{
			return A/B;
		}catch(string err){
			throw(err);
		}
	}else{
		throw((string)"Wrong operation!");
	}

}

CMatrix evaluate(CMatrix &A, string op, double B, bool doubleFirst){

	if (op == ".+"){
		return A+B;
	}else if (op == ".-"){
		return A-B;
	}else if (op == ".^" || op == "^"){
		return A.power(B);
	}
	else if (op == "./"){
		CMatrix temp = A;
		temp.dotSlash(B, doubleFirst);
		return temp;
	}else{
		throw((string)"Wrong operation!");
	}

}

string addMatrix(CMatrix &m){

	static int index = 0;
	if(m.get_nR() == 1 && m.get_nC() == 1){
		return m.get_string();
	}
	char* buffer = new char[10];
	sprintf(buffer,"mat%d", index);
	CMatrix temp = m;
	temp.set_name((string)buffer);
	matrix_list.push_back(temp);
	string toreturn = (string)buffer;
	delete[] buffer;
	index++;
	return toreturn;

}

void assignAndPrint(vector<string> &command_list){

	if (command_list.size() != 3){
		checkOperation(command_list);
	}
	if (command_list.size() != 3){
		throw((string)"Error in assigning");
	}else{
		CMatrix temp;
		if(isNumber(command_list[2])){
			temp =  CMatrix(command_list[2]);
		}else{
			for (int i=0; i<matrix_list.size(); i++){
				if(matrix_list[i].get_name() == command_list[2]){
					temp = matrix_list[i];
					break;
				}
			}
		}
		temp.set_name(command_list[0]);
		matrix_list.push_back(temp);
		temp.PrintMatrix();
	}

}
