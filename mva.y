%{
/* User Code Block*/ 
#include <stdio.h>
#include <string.h>
#define YYSTYPE char*
//buffer to store classification expression
char fbuf[1000]		;
//store the result
char expVal[1000]	;
int  fbufHasRead	;
//yyinput function will be used by flex
int readBuff(char*buff, int maxlen){
	int ret;
	if(fbufHasRead==1)
		return 0;
	else
		fbufHasRead=1;
	if(strlen(fbuf)>maxlen-1){
		ret=maxlen-1;
		strncpy (buff,fbuf,ret);
	}
	else{
		ret=strlen(fbuf);
		strcpy(buff,fbuf);
	}
	return ret;	
}
//procedure when error happen during parsing the classification expression.
void yyerror(const char *str){
	fprintf(stderr,"error: %s\n",str);
}
//Default yacc yywrap() function in order to support processing of multiple input files as one. We don't need it here. 
int yywrap()
{
	return 1;
}
//Parse a classifiction expression.
//input variable	: expr   --> classification expression
//output variable:	: outpur --> Reverse Polish notation which will be used to calculate the value of classification expression
//return 		: 0      --> classification expression is valid
//			: !=0	 --> classification expression is not valid
int parseExpression(const char* expr, char* output)
{
	int ret;
	fbufHasRead=0;
	sprintf(fbuf,expr);
	ret=yyparse();
	sprintf(output,"%s",expVal);
	return ret;
}
%}

%{ 
	/*RULE BLOCK*/
%}
%token VARIABLE
%%
exp: comExp '+' comExp 
{
	int l1=strlen($1)		; 
	int l2=strlen($3)		; 
	$$=(char*)malloc(l1+l2+4)	;
	sprintf($$,"%s|%s|+",$3,$1)	;
	free($1)			;
	free($3)			;
	sprintf(expVal,"%s",$$)		;
}	/*classification expression*/
|
'(' exp ')' 
{
	$$=strdup($2); free($2);
}
;
comExp:
var  {$$=strdup($1); free($1);}
|
'(' exp ')' {$$=strdup($2);  free($2);}
;
var:VARIABLE 
{
	$$=strdup($1);
	free($1);
}	/*variable name*/
%%
