//
//  symdiff.c
//  physher
//
//  Created by Mathieu Fourment on 2/01/2018.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#include "symdiff.h"

#include <string.h>

#include "mstring.h"
#include "matrix.h"

ExpressionItem* new_ExpressionItem(char* expression){
	ExpressionItem* item = malloc(sizeof(ExpressionItem));
	item->expr = expression;
	item->dev = NULL;
	item->op = 0;
	return item;
}

ExpressionStack* new_ExpressionStack(size_t capacity){
	ExpressionStack* stack = malloc(sizeof(ExpressionStack));
	stack->expressions = calloc(capacity, sizeof(ExpressionItem*));
	stack->count = 0;
	stack->capacity = capacity;
	return stack;
}

void free_Stack(ExpressionStack* stack){
	for (int i = 0; i < stack->count; i++) {
		free(stack->expressions[i]->expr);
		free(stack->expressions[i]->dev);
		free(stack->expressions[i]);
	}
	free(stack->expressions);
	free(stack);
}

void ExpressionStack_insert(ExpressionStack* stack, size_t index, ExpressionItem* item){
	if (stack->capacity == stack->count) {
		stack->capacity *= 2;
		stack->expressions = realloc(stack->expressions, sizeof(ExpressionItem*)*stack->capacity);
	}
	
	if(stack->count != index){
		memmove(stack->expressions+index+1, stack->expressions+index, sizeof(ExpressionItem*)*(stack->count-index));
		stack->expressions[index] = item;
	}
	else{
		stack->expressions[stack->count] = item;
	}
	stack->count++;
}

void ExpressionStack_add(ExpressionStack* stack, ExpressionItem* item){
	if (stack->capacity == stack->count) {
		stack->capacity *= 2;
		stack->expressions = realloc(stack->expressions, sizeof(ExpressionItem*)*stack->capacity);
	}
	stack->expressions[stack->count] = item;
	stack->count++;
}

#define DIGITS			"0123456789."
#define IsDigit(c)		(c && strchr(DIGITS, c) != NULL)

bool IsNumeric(const char* lpcs){
	char* p = (char*)lpcs;
	if(*p == '-' || *p == '+')
		p++;
	if(*p == 'p' && *(p+1) == 'i' && *(p+2) == 0)
		return true;
	while (*p){
		if(IsDigit(*p) == false)
			return false;
		p++;
	}
	return true;
}

bool IsNumericRange(const char* lpcs, size_t start, size_t end){
	const char* p = lpcs+start;
	if(*p == '-' || *p == '+'){
		p++;
		start++;
	}
	if(start < end+1 && *p == 'p' && *(p+1) == 'i' && *(p+2) == 0)
		return true;
	while (start != end){
		if(IsDigit(*p) == false)
			return false;
		p++;
		start++;
	}
	return true;
}

bool isVariable(const char* lpcs){
	const char* p = lpcs;
	if(*p == '-' || *p == '+')
		p++;
	// must start with letter or underscore
	if (!(*p >= 65 && *p <=90) && !(*p >= 97 && *p <=122) && *p != '_') {
		return false;
	}
	p++;
	while (*p) {
		if (!((*p >= 65 && *p <=90) || (*p >= 97 && *p <=122)  || (*p >= 48 && *p <=57) || *p == '_')) {
			return false;
		}
		p++;
	}
	return true;
}

bool isVariableRange(const char* lpcs, size_t start, size_t end){
	const char* p = lpcs+start;
	if(*p == '-' || *p == '+')
		p++;
	// must start with letter or underscore
	if (!(*p >= 65 && *p <=90) && !(*p >= 97 && *p <=122) && *p != '_') {
		return false;
	}
	p++;
	while (start != end) {
		if (!((*p >= 65 && *p <=90) || (*p >= 97 && *p <=122)  || (*p >= 48 && *p <=57) || *p == '_')) {
			return false;
		}
		start++;
		p++;
	}
	return true;
}

bool IsRightSign(char c, const char* lpcsOperators[], int nIndex, size_t len){
	for(; nIndex < len; nIndex++)
		if(strchr(lpcsOperators[nIndex], c) != NULL)
			return false;
	return true;
}

int GetOperator(const char* lpcs, const char* lpcsOperators[], size_t n){
	for(int nIndex = 0; nIndex < n; nIndex++){
		int open = 0;
		// scan the expression from its end
		const char* p = lpcs+strlen(lpcs)-1;
		// loop tell reach expression start
		while(p != lpcs){
			// check for close
			if(*p == ')'){
				open--;
			}
			// check for open
			else	if(*p == '('){
				open++;
			}
			// check for operator
			else if(open == 0 && strchr(lpcsOperators[nIndex], *p) != NULL)
				// check if the operator in not a sign mark
				if((*p != '-' && *p != '+') || (p != lpcs && IsRightSign(*(p-1), lpcsOperators, nIndex+1, n)))
					// return operator index
					return (int)(p-lpcs);
			p--;
		}
	}
	// operator not found
	return -1;
}

ExpressionStack* FillStack(const char* lpcsInput){
	// operators array from high to low priority
	const char* lpcsOperators[3] = { "+-", "*/", "^%" };
	ExpressionStack* stack = new_ExpressionStack(100);
	
	// insert first input into the stack
	ExpressionStack_add(stack, new_ExpressionItem(String_clone(lpcsInput)));
//printf("bbb\n");
	// loop in Expression stack to check if any Expression can be divided to two queries
	for(int i = 0; i < stack->count; i++){
		// check if Expression item is operator
		if(stack->expressions[i]->op == 0){
			char* str = stack->expressions[i]->expr;
			// parse expression to find operators
			int nOpIndex = GetOperator(str, lpcsOperators, 3);
//			printf("aaa\n");
			if(nOpIndex != -1){	// split the Expression into two queries at the operator index
				stack->expressions[i]->op = str[nOpIndex];

				char* substr = malloc(sizeof(char)*(nOpIndex+1));
				memcpy(substr, str, sizeof(char)*nOpIndex);
				substr[nOpIndex] = '\0';
				ExpressionStack_insert(stack, i+1, new_ExpressionItem(substr));

				size_t len = strlen(str);
				char* substr2 = malloc(sizeof(char)*(len-nOpIndex));
				memcpy(substr2, str+nOpIndex+1, sizeof(char)*(len-nOpIndex));// include \0
				ExpressionStack_insert(stack, i+2, new_ExpressionItem(substr2));
			}
			else if(str[0] == '('){
				int open = 1;
				int j = 1;
				while (open != 0) {
					if (str[j] == '(') {
						open++;
					}
					else if (str[j] == ')'){
						open--;
					}
					j++;
				}
				j--;
				memmove(str+j, str+j+1, sizeof(char)*(strlen(str)-j));//remove extra )
				memmove(str, str+1, sizeof(char)*(strlen(str))); // remove extra (
				i--;
			}
			// string starts with + (positive): remove it
			else if(str[0] == '+'){
				memmove(str, str+1, sizeof(char)*(strlen(str)));
			}
			// expression starts with -(: flip the sign of expression inside the parentethis
			else if(strlen(str) > 1 && str[0] == '-' && str[1] == '('){
//				printf("{%s}\n", str);
				memmove(str, str+1, sizeof(char)*(strlen(str))); // remove -
				
				int open = 0;
				for (int j = 0; j < strlen(str); j++) {
					if(str[j] == ')'){
						open--;
					}
					else if(str[j] == '('){
						open++;
					}
					else if(open == 1){
						if (str[j] == '-') {
							if(str[j-1] == '('){
								memmove(str+j, str+j+1, sizeof(char)*(strlen(str)-j));
							}
							else{
								str[j] = '+';
							}
						}
						else if (str[j] == '+') {
							str[j] = '-';
						}
						else if(str[j-1] == '('){
							str = stack->expressions[i]->expr = realloc(stack->expressions[i]->expr, sizeof(char)*(strlen(str)+2));
							memmove(str+j+1, str+j, sizeof(char)*(strlen(str)-j+1));
							str[j] = '-';
						}
					}
				}
				stack->expressions[i]->expr = str;
				i--;
//				printf("[%s]\n", str);
			}
		}
	}
	return stack;
}

void OptimizeSign(char* str){
	// replace "--" with "" or "+"
	for (int i = 0; i < strlen(str)-1; i++) {
		if (str[i] == '-' && str[i+1] == '-') {
			if(i == 0 || strchr("(+-/*^", str[i-1]) != NULL){
				memmove(str+i, str+i+2, sizeof(char)*(strlen(str)-i-1));
			}
			else{
				memmove(str+i, str+i+1, sizeof(char)*(strlen(str)-i));
				str[i] = '+';
			}
			i--;
		}
	}
	
	// replace "+-" with "-"
	for (int i = 0; i < strlen(str)-1; i++) {
		if (str[i] == '+' && str[i+1] == '-') {
			memmove(str+i, str+i+1, sizeof(char)*(strlen(str)-i));
		}
	}
}

void Optimize(char* str){
	size_t nLength = strlen(str);
	
	OptimizeSign(str);
	
	// replace "((....))"  with "(....)"
	for (int i = 0; i < strlen(str)-1; i++) {
		if (str[i] == '(' && str[i+1] == '(') {
			int open = 2;
			int j = i+2;
			bool flag = false;
			while (open != 0) {
				if (str[j] == '(') {
					open++;
				}
				else if (str[j] == ')'){
					open--;
					if(strchr("+-*/", str[j+1]) != NULL){
						flag = true;
					}
				}
				j++;
			}
			j--;
			if(str[j-1] == ')' && flag == false){
				memmove(str+j-1, str+j, sizeof(char)*(strlen(str)-j+1));//remove extra )
				memmove(str+i, str+i+1, sizeof(char)*(strlen(str)-i)); // remove extra (
			}
		}
	}

	// remove any 1*
	for (int i = 0; i < strlen(str)-1; i++) {
		if (str[i] == '1' && str[i+1] == '*' && (i == 0 || strchr("+-*(", str[i-1]) != NULL)) {
			memmove(str+i, str+i+2, sizeof(char)*(strlen(str)-i-1));
		}
	}
	
	// remove any *1
	for (int i = 0; i < strlen(str)-1; i++) {
		if (str[i] == '*' && str[i+1] == '1' && (i+2 == strlen(str) || strchr("+-*(", str[i+2]) == NULL)) {
			memmove(str+i, str+i+2, sizeof(char)*(strlen(str)-i-1));
		}
	}
	
	// remove any exponent equal 1
	for (int i = 0; i < strlen(str)-1; i++) {
		if (str[i] == '^' && str[i+1] == '1' && (i+2 == strlen(str) || strchr("+-*(", str[i+2]) != NULL)) {
			memmove(str+i, str+i+2, sizeof(char)*(strlen(str)-i-1));
		}
	}
	
	// remove unneeded parentheses
	for (int i = 0; i < strlen(str)-1; i++) {
		if (str[i] == '('){
			// find the parenthesis close
			int open = 1;
			int j = i+1;
			while (open != 0) {
				if (str[j] == '(') {
					open++;
				}
				else if (str[j] == ')'){
					open--;
				}
				j++;
			}
			j--;
			
			// check if the parentheses in the start and the end of the string
			if((i == 0 && j == strlen(str)-1) ||
			   j == i+2 ||
			   // check if the string doesn't include any operator
			   IsNumericRange(str, i+1, j-1) || isVariableRange(str, i+1, j-1)){
				// delete the far index of ')'
				memmove(str+j, str+j+1, sizeof(char)*(strlen(str)-j));
				// delete the near index of '('
				memmove(str+i, str+i+1, sizeof(char)*(strlen(str)-i));
			}
		}
	}
	
	if(nLength != strlen(str)){
		Optimize(str);
	}
}

//CString TrimFloat(double f){
//	CString str;
//	str.Format("%f", f);
//	if(str.Find('.') != -1)
//		str.TrimRight('0');
//	str.TrimRight('.');
//	return str;
//}

char* GetInput(ExpressionItem* item){
	if(IsNumeric(item->expr) == false && isVariable(item->expr) == false &&
	   (String_contains(item->expr, '+')|| String_contains(item->expr, '-') || String_contains(item->expr, '*') ||
		String_contains(item->expr, '/'))){
		
		StringBuffer* buffer = new_StringBuffer(strlen(item->expr)+3);
		StringBuffer_append_format(buffer, "(%s)", item->expr);
		char* temp = StringBuffer_tochar(buffer);
		free_StringBuffer(buffer);
		return temp;
	}
	return String_clone(item->expr);
}

char* GetDifferentiation(ExpressionItem* item, const char* dx){
	StringBuffer* buffer = new_StringBuffer(10);
	int nIndex;
//	if(item->expr[0] == 'd' && item->expr[1] == '('){
//		nIndex = m_strInput.Find('(');
//		CString str = m_strInput.Mid(nIndex+1);
//		// get the string between function parentheses
//		str = str.Left(str.ReverseFind(')'));
//		switch(m_nFunction)
//		{
//			case 0:		m_strOutput = d(str, m_strStack);			break;
//			case 1:		m_strOutput = d_sin(str, m_strStack);		break;
//			case 2:		m_strOutput = d_cos(str, m_strStack);		break;
//			case 3:		m_strOutput = d_tan(str, m_strStack);		break;
//			case 4:		m_strOutput = d_sec(str, m_strStack);		break;
//			case 5:		m_strOutput = d_cosec(str, m_strStack);		break;
//			case 6:		m_strOutput = d_cot(str, m_strStack);		break;
//			case 7:		m_strOutput = d_sinh(str, m_strStack);		break;
//			case 8:		m_strOutput = d_cosh(str, m_strStack);		break;
//			case 9:		m_strOutput = d_tanh(str, m_strStack);		break;
//			case 10:	m_strOutput = d_sech(str, m_strStack);		break;
//			case 11:	m_strOutput = d_cosech(str, m_strStack);	break;
//			case 12:	m_strOutput = d_coth(str, m_strStack);		break;
//			case 13:	m_strOutput = d_asin(str, m_strStack);		break;
//			case 14:	m_strOutput = d_acos(str, m_strStack);		break;
//			case 15:	m_strOutput = d_atan(str, m_strStack);		break;
//			case 16:	m_strOutput = d_asec(str, m_strStack);		break;
//			case 17:	m_strOutput = d_acosec(str, m_strStack);	break;
//			case 18:	m_strOutput = d_acot(str, m_strStack);		break;
//			case 19:	m_strOutput = d_asinh(str, m_strStack);		break;
//			case 20:	m_strOutput = d_acosh(str, m_strStack);		break;
//			case 21:	m_strOutput = d_atanh(str, m_strStack);		break;
//			case 22:	m_strOutput = d_asech(str, m_strStack);		break;
//			case 23:	m_strOutput = d_acosech(str, m_strStack);	break;
//			case 24:	m_strOutput = d_acoth(str, m_strStack);		break;
//			case 25:	m_strOutput = d_sqrt(str, m_strStack);		break;
//			case 26:	m_strOutput = d_log10(str, m_strStack);		break;
//			case 27:	m_strOutput = d_log(str, m_strStack);		break;
//			case 28:	m_strOutput = d_ln(str, m_strStack);		break;
//			case 29:	m_strOutput = d_sign(str, m_strStack);		break;
//			case 30:	m_strOutput = d_abs(str, m_strStack);		break;
//		}
//		m_strOutput = (m_nSign == -1?"-":"")+m_strOutput;
//		Differentiate(u,stack)
//	}
//	else{
		// dx/dx = 1
	char* m_strInput = item->expr;
		if(strcmp(m_strInput, dx) == 0 || (strlen(m_strInput) > 1 && m_strInput[0] == '+' && strcmp(m_strInput+1, dx) == 0)){
			StringBuffer_append_string(buffer, "1");
		}
		else if(strlen(m_strInput) > 1 && m_strInput[0] == '-' && strcmp(m_strInput+1, dx) == 0){
			StringBuffer_append_string(buffer, "-1");
		}
		else if(IsNumeric(m_strInput)){
			// dc/dx = 0, where c is constant
			StringBuffer_append_string(buffer, "0");
		}
		else if(isVariable(m_strInput)){
			// dc/dx = 0, where c is a variable but not x
			StringBuffer_append_string(buffer, "0");
		}
		else{
			printf("\n%s\n", m_strInput);
			// du/dx, where u is a function of x
			StringBuffer_append_strings(buffer, 4, "d", m_strInput, "/d", dx);exit(3);
		}
//	}
	char* temp = StringBuffer_tochar(buffer);
	free_StringBuffer(buffer);
	return temp;
}

char* DifferentiateStack(ExpressionStack* stack, int* nExpression, const char* dx){
	ExpressionItem *expression = stack->expressions[*nExpression];
	(*nExpression)++;
	
	if(expression->op) {
		StringBuffer* buffer = new_StringBuffer(10);
		// get left operand
		char* u = GetInput(stack->expressions[*nExpression]);
		// get left operand differentiation
		char* du = DifferentiateStack(stack, nExpression, dx);
		//		printf("u: %s du: %s\n", u, du);
		// get right operand
		char* v = GetInput(stack->expressions[*nExpression]);
		// get right operand differentiation
		char* dv = DifferentiateStack(stack, nExpression, dx);
		//		printf("v: %s dv: %s\n", v, dv);
		
		switch(expression->op){
			// d(u-v) = du-dv
			case '-':{
				if(strcmp(du, dv) == 0){
					StringBuffer_append_string(buffer, "0");
				}
				else if(strcmp(du, "0") == 0){
					StringBuffer_append_format(buffer, "(-%s)", dv);
				}
				else if(strcmp(dv, "0") == 0){
					StringBuffer_append_format(buffer, "%s", du);
				}
				else{
					StringBuffer_append_format(buffer, "(%s-%s)", du, dv);
				}
				break;
			}
			// d(u+v) = du+dv
			case '+':{
				if((du[0] == '-' && dv[0] != '-' && strcmp(du+1, dv) == 0) ||
				   (du[0] != '-' && dv[0] == '-' && strcmp(du, dv+1) == 0) ){
					StringBuffer_append_string(buffer, "0");
				}
				else if(strcmp(du, "0") == 0){
					StringBuffer_append_format(buffer, "%s", dv);
				}
				else if(strcmp(dv, "0") == 0){
					StringBuffer_append_format(buffer, "%s", du);
				}
				else{
					StringBuffer_append_format(buffer, "(%s%c%s)", du, expression->op, dv);
				}
				break;
			}
			// d(u*v) = u*dv+du*v
			case '*':{
				bool zero_term1 = false;
				bool zero_term2 = false;
				if(strcmp(u, "0") == 0 || strcmp(dv, "0") == 0){
					zero_term1 = true;
				}
				if(strcmp(du, "0") == 0 || strcmp(v, "0") == 0){
					zero_term2 = true;
				}
				if(zero_term1 && zero_term2){
					StringBuffer_append_string(buffer, "0");
				}
				else if(zero_term1){
					StringBuffer_append_format(buffer, "(%s*%s)", du, v);
				}
				else if(zero_term2){
					StringBuffer_append_format(buffer, "(%s*%s)", u, dv);
				}
				else{
					StringBuffer_append_format(buffer, "(%s*%s+%s*%s)", u, dv, du, v);
				}
				break;
			}
			// d(u/v) = (du*v-u*dv)/v^2
			case '/':{
				if( (strcmp(du, "0") == 0 && strcmp(u, "0") == 0)){
					StringBuffer_append_string(buffer, "0");
				}
				else if( (strcmp(du, "0") == 0)){
					StringBuffer_append_format(buffer, "(-%s*%s)/(%s)^2", u, dv, v);
				}
				else if( (strcmp(u, "0") == 0)){
					StringBuffer_append_format(buffer, "(%s*%s)/(%s)^2", du, v, v);
				}
				else{
					StringBuffer_append_format(buffer, "(%s*%s-%s*%s)/(%s)^2", du, v, u, dv, v);
				}
				break;
			}
			// d(u^v) = v*u^(v-1)*du+u^v*ln(u)*dv
			case '^':{
				// v constant
				// d(u^v) = v*u^(v-1)*du
				if(strcmp(dv, "1") == 0){
					StringBuffer_append_format(buffer, "%s", du);
				}
				else if(strcmp(dv, "0") == 0){
					double exponent = atof(v)-1;
					if(exponent == 1){
						StringBuffer_append_format(buffer, "%s*%s*%s", v, u, du);
					}
					else{
						StringBuffer_append_format(buffer, "%s*%s^%f*%s", v, u, atof(v)-1, du);
					}
				}
				else{
					StringBuffer_append_format(buffer, "(%s*%s^(%s-1)*%s+%s^%s*ln(%s)*%s)", v, u, v, du, u, v, u, dv);
				}
				break;
			}
			default:{
				fprintf(stderr, "du!=0 d Operator %c\n", expression->op);
				exit(1);
			}
		}
		expression->dev = StringBuffer_tochar(buffer);
		
		free(v);
		free(u);
		free_StringBuffer(buffer);
	}
	else{
		// get Expression differentiation
		expression->dev = GetDifferentiation(expression, dx);
	}
	
	return expression->dev;
}



double evalueStack(ExpressionStack* stack, int* nExpression, Parameters* list){
	ExpressionItem *expression = stack->expressions[*nExpression];
	(*nExpression)++;
	if(expression->op) {

		double du = evalueStack(stack, nExpression, list);
		double dv = evalueStack(stack, nExpression, list);
		//		printf("v: %s dv: %s\n", v, dv);
		
		switch(expression->op){
				// d(u-v) = du-dv
			case '-':{
				expression->value = du - dv;
				break;
			}
				// d(u+v) = du+dv
			case '+':{
				expression->value = du + dv;
				break;
			}
				// d(u*v) = u*dv+du*v
			case '*':{
				expression->value = du * dv;
				break;
			}
				// d(u/v) = (du*v-u*dv)/v^2
			case '/':{
				expression->value = du/dv;
				break;
			}
				// d(u^v) = v*u^(v-1)*du+u^v*ln(u)*dv
			case '^':{
				expression->value = pow(du, dv);
				break;
			}
			default:{
				fprintf(stderr, "du!=0 d Operator %c\n", expression->op);
				exit(1);
			}
		}
	}
	else if(IsNumeric(expression->expr)){
		expression->value = atof(expression->expr);
	}
	else{
		if(isVariable(expression->expr)){
			int index = 0;
			if (expression->expr[0] == '-') {
				index = 1;
			}
			
			for (int i = 0; i < Parameters_count(list); i++) {
				if (strcmp(Parameters_name(list, i), expression->expr+index) == 0) {
					expression->value = Parameters_value(list, i);
					break;
				}
			}
			if(index!=0)expression->value =  -expression->value;
		}
		else{
			fprintf(stderr, "sdfsn");
			exit(11);
		}
	}
	
	return expression->value;
}

double differentiate2(ExpressionStack* stack, Parameters* list){
	int nExpression = 0;
	return evalueStack(stack, &nExpression, list);
}

char* differentiate(char* lpcsInput, const char* dx){
	char* strInput = lpcsInput;
	// remove spaces
//	strInput.Remove(' ');
	// make all characters lower case
//	strInput.MakeLower();
	// Optimize "--"
	OptimizeSign(strInput);
	
	ExpressionStack* stack = FillStack(strInput);
	
//	for (int i = 0; i < stack->count; i++) {
//		if(stack->expressions[i]->op != 0){
//			printf("%c\n", stack->expressions[i]->op);
//		}
//		else{
//			printf("%s\n", stack->expressions[i]->expr);
//		}
//	}
	
	int nExpression = 0;
	// apply operators to operands
	char* strOutput = DifferentiateStack(stack, &nExpression, dx);
	char* out = String_clone(strOutput);
	Optimize(strOutput);
	free_Stack(stack);
	return out;
}

void test_symdiff(){
	char** vs = malloc(16*sizeof(char*));
	/*
	 pa = phi1
	 pc = (1-phi1)*phi2
	 pg = (1 - (pa+pc))*phi3
	 pt = 1 - (pa+pc+pg)*/
	char* pia = "phi1";
	char* pic = "(1-phi1)*phi2";
	char* pig = "(1-(phi1+((1-phi1)*phi2)))*phi3";
	char* pit = "(1-(phi1+(1-phi1)*phi2+(1-(phi1+((1-phi1)*phi2)))*phi3))";
	
	//	vs[0] = String_clone("-(a*phi2+b*phi3+c)");
	//	vs[1] = String_clone("a*phi2");
	//	vs[2] = String_clone("b*phi3");
	//	vs[3] = String_clone("c");
	//
	//	vs[4] = String_clone("a*phi1");
	//	vs[5] = String_clone("-(a*phi1+d*phi3+e)");
	//	vs[6] = String_clone("d*phi3");
	//	vs[7] = String_clone("e");
	//
	//	vs[8] = String_clone("b*phi1");
	//	vs[9] = String_clone("d*phi2");
	//	vs[10] = String_clone("-(b*phi1+d*phi2+1)");
	//	vs[11] = String_clone("1");
	//
	//	vs[12] = String_clone("c*phi1");
	//	vs[13] = String_clone("e*phi2");
	//	vs[14] = String_clone("phi3");
	//	vs[15] = String_clone("-(c*phi1+e*phi2+phi3)");
	
	vs[0] = String_clone("-(a*(1-phi1)*phi2+b*(1-(phi1+((1-phi1)*phi2)))*phi3+c*(1-(phi1+(1-phi1)*phi2+(1-(phi1+((1-phi1)*phi2)))*phi3)))");
	//	vs[0] = String_clone("-(a+x^(-2))");
	vs[1] = String_clone("a*(1-phi1)*phi2");
	vs[2] = String_clone("b*(1-(phi1+((1-phi1)*phi2)))*phi3");
	vs[3] = String_clone("c*(1-(phi1+(1-phi1)*phi2+(1-(phi1+((1-phi1)*phi2)))*phi3))");
	
	vs[4] = String_clone("a*phi1");
	vs[5] = String_clone("-(a*phi1+d*(1-(phi1+((1-phi1)*phi2)))*phi3+e*(1-(phi1+(1-phi1)*phi2+(1-(phi1+((1-phi1)*phi2)))*phi3)))");
	vs[6] = String_clone("d*(1-(phi1+((1-phi1)*phi2)))*phi3");
	vs[7] = String_clone("e*(1-(phi1+(1-phi1)*phi2+(1-(phi1+((1-phi1)*phi2)))*phi3))");
	
	vs[8] = String_clone("b*phi1");
	vs[9] = String_clone("d*(1-phi1)*phi2");
	vs[10] = String_clone("-(b*phi1+d*(1-phi1)*phi2+(1-(phi1+(1-phi1)*phi2+(1-(phi1+((1-phi1)*phi2)))*phi3)))");
	vs[11] = String_clone("(1-(phi1+(1-phi1)*phi2+(1-(phi1+((1-phi1)*phi2)))*phi3))");
	
	vs[12] = String_clone("c*phi1");
	vs[13] = String_clone("e*(1-phi1)*phi2");
	vs[14] = String_clone("(1-(phi1+((1-phi1)*phi2)))*phi3");
	vs[15] = String_clone("-(c*phi1+e*(1-phi1)*phi2+(1-(phi1+((1-phi1)*phi2)))*phi3)");
	
	//	char*freqs[4] = {"phi1", "phi2", "phi3", "1"};
	char*freqs[4] = {pia, pic, pig, pit};
	
	StringBuffer*buffer = new_StringBuffer(100);
	
	StringBuffer_append_string(buffer, vs[0]);
	StringBuffer_append_string(buffer, "/(");
	int j =0;
	for (int i = 0; i < 4; i++) {
		j = 4*i+ i;
		if(freqs[i][0] == '1'){
			StringBuffer_append_format(buffer, "%s", vs[j]+1);
		}
		else{
			StringBuffer_append_format(buffer, "%s*%s", vs[j]+1, freqs[i]);
		}
		if (j < 16-1) {
			StringBuffer_append_char(buffer, '+');
		}
	}
	//StringBuffer_append_string(buffer, ")/(1+phi1+phi2+phi3)^2");
	StringBuffer_append_string(buffer, ")");
	printf("%s\n", buffer->c);
	
	char* output = differentiate(buffer->c, "a");
	//	char* output = differentiate("a+(2+b)*(a+c)", "a");
	
	printf("output: ==%s==\n", output);
	
	Parameters* list = new_Parameters(8);
	Parameters_add(list, new_Parameter("a", 1, NULL));
	Parameters_add(list, new_Parameter("b", 1, NULL));
	Parameters_add(list, new_Parameter("c", 1, NULL));
	Parameters_add(list, new_Parameter("d", 1, NULL));
	Parameters_add(list, new_Parameter("e", 1, NULL));
	
	Parameters_add(list, new_Parameter("phi1", 0.5, NULL));
	Parameters_add(list, new_Parameter("phi2", 0.5, NULL));
	Parameters_add(list, new_Parameter("phi3", 0.5, NULL));
	
	ExpressionStack* stack = FillStack(output);
	int nExpression = 0;
	//	DifferentiateStack(stack, &nExpression, "a");
	
	double df = differentiate2(stack, list);
	printf("%f\n", df);
	
	free_cmatrix(vs, 16);
	free_StringBuffer(buffer);
	//	free(output);
}

//cVector* infixToRPN(char* str){
//
//	cVector* Q = new_cVector(100); // queue
//	cVector* S = new_cVector(100); // stack
//
//	for (size_t i = 0 ; i < strlen(str); i++) {
//		if (str[i] == '('){
//			cVector_push(S, "(");
//		}
//		else if (str[i] == ')'){
//			while (strcmp(cVector_peek(S), "(") != 0) {
//				char* temp = cVector_pop(S);
//				cVector_push(Q, temp);
//				free(temp);
//			}
//			cVector_pop(S);
//		}
//		// an operator
//		else if (strchr("+-/*^", str[i]) != NULL) {
//			while (!S.empty() && prededence.get(token) <= prededence.get(S.peek())) {
//				cVector_push(Q, cVector_pop(S));
//			}
//			cVector_push(S, token);
//			continue;
//		}
//		// variable or number
//		else {
//			// loop here until the end of variable or number
////			cVector_push(Q, token);
//		}
//	}
//	// at the end, pop all the elements in S to Q
//	while (cVector_length(S) != 0) {
//		cVector_push(Q, cVector_pop(S));
//	}
//	free_cVector(S);
//
//	return Q;
//}
