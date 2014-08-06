%{
/* User Code Block
   will be cpoied directly to the .c source file by flex
*/
#define YYSTYPE char*
#include <stdio.h>
#include "y.tab.h"
extern YYSTYPE yylval;
extern int readBuff(char *buff, int buffsize);
#define YY_INPUT(buff, res, buffsize) (res = readBuff(buff, buffsize))
%}

%%
[A-Za-z][A-Z_a-z0-9]+	{
				yylval=strdup(yytext);
				return VARIABLE;
			  /* Recognize valid varibale name. A valid variable name should start with character "A-Za-z" and can be followed with character "A-Za-z0-9" */
			}
[+()]			{
				return *yytext ;
			  /* Recognize valid operator. In this package we only use three operators, namely "(",")" and "+" */
			}
\n			  /* Ignore end of line */;
[ \t]+			  /* Ignore whitespace */;
%%
