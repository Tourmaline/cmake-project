#ifndef lint
static const char yysccsid[] = "@(#)yaccpar	1.9 (Berkeley) 02/21/93";
#endif

#define YYBYACC 1
#define YYMAJOR 1
#define YYMINOR 9
#define YYPATCH 20140101

#define YYEMPTY        (-1)
#define yyclearin      (yychar = YYEMPTY)
#define yyerrok        (yyerrflag = 0)
#define YYRECOVERING() (yyerrflag != 0)

#define YYPREFIX "yy"

#define YYPURE 0

#line 2 "ergo_input_processor.y"

  /** @file ergo_input_processor.c Parses the input.
      Uses bison code generator to generate the parses.
  */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ergo_scripted.h"
#define YYERROR_VERBOSE
int yylex(void);
int yyerror(const char *s);
static const char *last_token = NULL;

#line 18 "ergo_input_processor.y"
#ifdef YYSTYPE
#undef  YYSTYPE_IS_DECLARED
#define YYSTYPE_IS_DECLARED 1
#endif
#ifndef YYSTYPE_IS_DECLARED
#define YYSTYPE_IS_DECLARED 1
typedef union {
  double num;     /* for returning numbers */
  char str[256];  /* for returning strings */
  struct variable *var; /* for returning lvalues */
} YYSTYPE;
#endif /* !YYSTYPE_IS_DECLARED */
#line 47 "ergo_input_processor.c"

/* compatibility with bison */
#ifdef YYPARSE_PARAM
/* compatibility with FreeBSD */
# ifdef YYPARSE_PARAM_TYPE
#  define YYPARSE_DECL() yyparse(YYPARSE_PARAM_TYPE YYPARSE_PARAM)
# else
#  define YYPARSE_DECL() yyparse(void *YYPARSE_PARAM)
# endif
#else
# define YYPARSE_DECL() yyparse(void)
#endif

/* Parameters sent to lex. */
#ifdef YYLEX_PARAM
# define YYLEX_DECL() yylex(void *YYLEX_PARAM)
# define YYLEX yylex(YYLEX_PARAM)
#else
# define YYLEX_DECL() yylex(void)
# define YYLEX yylex()
#endif

/* Parameters sent to yyerror. */
#ifndef YYERROR_DECL
#define YYERROR_DECL() yyerror(const char *s)
#endif
#ifndef YYERROR_CALL
#define YYERROR_CALL(msg) yyerror(msg)
#endif

extern int YYPARSE_DECL();

#define NUMBER 257
#define DOT 258
#define SYMBOL 259
#define EQUAL 260
#define STRING 261
#define EOFTAG 262
#define GETEXC 263
#define GETPOL 264
#define K_ALL 265
#define HELP 266
#define MOLTAG 267
#define GHOSTTAG 268
#define MOLDAL 269
#define QUIT 270
#define RUNTAG 271
#define SYSTEM 272
#define GHOST 273
#define ANGSTROM 274
#define PRECISION 275
#define RANGE 276
#define WARRANTY 277
#define SET_NTHREADS 278
#define PLUS 279
#define MINUS 280
#define TIMES 281
#define DIVIDE 282
#define POWER 283
#define LEFT_PARENTHESIS 284
#define RIGHT_PARENTHESIS 285
#define EOL 286
#define NEG 287
#define YYERRCODE 256
static const short yylhs[] = {                           -1,
    0,    0,    0,    3,    3,    3,    3,    4,    4,    2,
    2,    5,    5,    5,    5,    5,    5,    5,    5,    5,
    5,    5,    5,    5,    5,    5,    5,    5,    5,    5,
    5,    6,    6,    7,    1,    1,    1,    1,    1,    1,
    1,    1,
};
static const short yylen[] = {                            2,
    0,    1,    3,    0,    1,    1,    1,    3,    3,    1,
    3,    3,    4,    3,    4,    1,    2,    1,    2,    3,
    2,    4,    4,    2,    1,    3,    4,    4,    1,    6,
    7,    1,    2,    5,    1,    3,    3,    3,    3,    2,
    3,    3,
};
static const short yydefred[] = {                         0,
    7,   10,    0,    0,    0,    0,    0,    0,   29,    0,
    0,   18,    0,   25,    0,    0,    0,    0,    5,    6,
    0,    0,    0,    0,    0,    0,    0,   19,    0,   21,
   24,    0,    0,    0,    0,    0,    0,   26,    0,    0,
    0,    0,   32,   12,    0,    0,   14,   20,    0,    0,
   35,    0,    0,    0,    0,   11,    9,    0,    3,   27,
   28,   13,    0,   33,   15,    0,    0,   23,    0,    0,
    0,    0,    0,    0,    0,   22,    0,    0,    0,   42,
    0,    0,    0,    0,    0,    0,   30,    0,   34,   31,
};
static const short yydgoto[] = {                         16,
   55,   17,   18,   19,   20,   44,   45,
};
static const short yysindex[] = {                      -214,
    0,    0, -243, -229, -225, -264, -258, -237,    0, -217,
 -213,    0, -252,    0, -238,    0, -247, -227,    0,    0,
 -197, -218, -182, -203, -169, -195, -169,    0, -167,    0,
    0, -168, -162, -257, -163, -249, -214,    0, -160, -159,
 -169, -158,    0,    0, -169, -169,    0,    0, -157, -156,
    0, -184, -254, -254, -208,    0,    0, -194,    0,    0,
    0,    0, -155,    0,    0, -154, -152,    0, -177, -201,
 -254, -254, -254, -254, -254,    0, -150, -153, -148,    0,
 -215, -215, -177, -177, -177, -176,    0, -149,    0,    0,
};
static const short yyrindex[] = {                         6,
    0,    0,    0,    0,    7,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,  111,    0,    0,
    0,    0,   14,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    6,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,   15,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    1,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
   29,   33,    9,   17,   25,    0,    0,    0,    0,    0,
};
static const short yygindex[] = {                        76,
  -34,  109,    0,    0,    0,   24,    0,
};
#define YYTABLESIZE 319
static const short yytable[] = {                         51,
   40,   58,   51,   52,   32,    1,   16,   51,   38,   24,
   35,   57,   36,   17,    8,   26,   39,   21,   69,   70,
   33,   25,   53,   28,   41,   53,   54,   27,   36,   54,
   53,   22,   37,    2,   54,   29,   81,   82,   83,   84,
   85,    1,   39,   30,    2,   34,   40,   31,    3,    4,
   47,    5,    6,    7,    8,    9,   10,   11,   37,   38,
   12,   13,   14,   15,   62,   73,   74,   75,   64,   65,
   71,   72,   73,   74,   75,   35,   76,   71,   72,   73,
   74,   75,   41,   80,   71,   72,   73,   74,   75,   42,
   46,   49,   43,   48,   50,   56,   60,   61,   63,   66,
   68,   77,   78,   67,   79,   75,   86,   87,   88,   89,
    2,   90,   59,   23,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,   40,
   40,   40,   40,    0,    0,   40,   40,   38,   38,   38,
   38,    4,   16,   38,   38,   39,   39,   39,   39,   17,
    8,   39,   39,   41,   41,   41,   41,   36,   36,   41,
   41,   37,   37,   36,   36,    0,    0,   37,   37,
};
static const short yycheck[] = {                        257,
    0,   36,  257,  261,  257,    0,    0,  257,    0,  274,
  258,  261,  260,    0,    0,  274,    0,  261,   53,   54,
  273,  286,  280,  261,    0,  280,  284,  286,    0,  284,
  280,  261,    0,  259,  284,  273,   71,   72,   73,   74,
   75,  256,  261,  261,  259,  284,  265,  261,  263,  264,
   27,  266,  267,  268,  269,  270,  271,  272,  286,  257,
  275,  276,  277,  278,   41,  281,  282,  283,   45,   46,
  279,  280,  281,  282,  283,  258,  285,  279,  280,  281,
  282,  283,  286,  285,  279,  280,  281,  282,  283,  259,
  286,  260,  262,  261,  257,  259,  257,  257,  257,  257,
  285,  257,  257,  260,  257,  283,  257,  261,  257,  286,
    0,  261,   37,    5,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,  279,
  280,  281,  282,   -1,   -1,  285,  286,  279,  280,  281,
  282,  286,  286,  285,  286,  279,  280,  281,  282,  286,
  286,  285,  286,  279,  280,  281,  282,  279,  280,  285,
  286,  279,  280,  285,  286,   -1,   -1,  285,  286,
};
#define YYFINAL 16
#ifndef YYDEBUG
#define YYDEBUG 0
#endif
#define YYMAXTOKEN 287
#define YYTRANSLATE(a) ((a) > YYMAXTOKEN ? (YYMAXTOKEN + 1) : (a))
#if YYDEBUG
static const char *yyname[] = {

"end-of-file",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"NUMBER","DOT","SYMBOL","EQUAL",
"STRING","EOFTAG","GETEXC","GETPOL","K_ALL","HELP","MOLTAG","GHOSTTAG","MOLDAL",
"QUIT","RUNTAG","SYSTEM","GHOST","ANGSTROM","PRECISION","RANGE","WARRANTY",
"SET_NTHREADS","PLUS","MINUS","TIMES","DIVIDE","POWER","LEFT_PARENTHESIS",
"RIGHT_PARENTHESIS","EOL","NEG","illegal-symbol",
};
static const char *yyrule[] = {
"$accept : Input",
"Input :",
"Input : Line",
"Input : Line EOL Input",
"Line :",
"Line : Assignment",
"Line : Command",
"Line : error",
"Assignment : Lvalue EQUAL Expression",
"Assignment : Lvalue EQUAL STRING",
"Lvalue : SYMBOL",
"Lvalue : Lvalue DOT SYMBOL",
"Command : MOLTAG EOL Molinput",
"Command : MOLTAG ANGSTROM EOL Molinput",
"Command : GHOSTTAG EOL Molinput",
"Command : GHOSTTAG ANGSTROM EOL Molinput",
"Command : HELP",
"Command : HELP Lvalue",
"Command : PRECISION",
"Command : MOLDAL STRING",
"Command : MOLDAL GHOST STRING",
"Command : RUNTAG STRING",
"Command : SET_NTHREADS LEFT_PARENTHESIS Expression RIGHT_PARENTHESIS",
"Command : SET_NTHREADS LEFT_PARENTHESIS STRING RIGHT_PARENTHESIS",
"Command : SYSTEM STRING",
"Command : WARRANTY",
"Command : GETEXC STRING NUMBER",
"Command : GETPOL STRING STRING NUMBER",
"Command : GETPOL STRING K_ALL NUMBER",
"Command : QUIT",
"Command : RANGE NUMBER EQUAL NUMBER NUMBER STRING",
"Command : RANGE GHOST NUMBER EQUAL NUMBER NUMBER STRING",
"Molinput : EOFTAG",
"Molinput : Molline Molinput",
"Molline : SYMBOL NUMBER NUMBER NUMBER EOL",
"Expression : NUMBER",
"Expression : Expression PLUS Expression",
"Expression : Expression MINUS Expression",
"Expression : Expression TIMES Expression",
"Expression : Expression DIVIDE Expression",
"Expression : MINUS Expression",
"Expression : Expression POWER Expression",
"Expression : LEFT_PARENTHESIS Expression RIGHT_PARENTHESIS",

};
#endif

int      yydebug;
int      yynerrs;

int      yyerrflag;
int      yychar;
YYSTYPE  yyval;
YYSTYPE  yylval;

/* define the initial stack-sizes */
#ifdef YYSTACKSIZE
#undef YYMAXDEPTH
#define YYMAXDEPTH  YYSTACKSIZE
#else
#ifdef YYMAXDEPTH
#define YYSTACKSIZE YYMAXDEPTH
#else
#define YYSTACKSIZE 10000
#define YYMAXDEPTH  10000
#endif
#endif

#define YYINITSTACKSIZE 200

typedef struct {
    unsigned stacksize;
    short    *s_base;
    short    *s_mark;
    short    *s_last;
    YYSTYPE  *l_base;
    YYSTYPE  *l_mark;
} YYSTACKDATA;
/* variables for the parser stack */
static YYSTACKDATA yystack;
#line 128 "ergo_input_processor.y"

YYSTYPE yylval;
int ergo_scanner_lineno = 1;
int ergo_scanner_reading_stdin = 0;

int yyerror(const char *s) {
  if (last_token) {
    printf("line %d: %s at '%s'\n",ergo_scanner_lineno, s, last_token);
    last_token = NULL;
  } else {
    printf("line %d: %s\n",ergo_scanner_lineno, s);
  }
  return !ergo_scanner_reading_stdin;
}

#ifdef SCANNER_TEST
int main(void) {
  yyparse();
}
#endif
#line 356 "ergo_input_processor.c"

#if YYDEBUG
#include <stdio.h>		/* needed for printf */
#endif

#include <stdlib.h>	/* needed for malloc, etc */
#include <string.h>	/* needed for memset */

/* allocate initial stack or double stack size, up to YYMAXDEPTH */
static int yygrowstack(YYSTACKDATA *data)
{
    int i;
    unsigned newsize;
    short *newss;
    YYSTYPE *newvs;

    if ((newsize = data->stacksize) == 0)
        newsize = YYINITSTACKSIZE;
    else if (newsize >= YYMAXDEPTH)
        return -1;
    else if ((newsize *= 2) > YYMAXDEPTH)
        newsize = YYMAXDEPTH;

    i = (int) (data->s_mark - data->s_base);
    newss = (short *)realloc(data->s_base, newsize * sizeof(*newss));
    if (newss == 0)
        return -1;

    data->s_base = newss;
    data->s_mark = newss + i;

    newvs = (YYSTYPE *)realloc(data->l_base, newsize * sizeof(*newvs));
    if (newvs == 0)
        return -1;

    data->l_base = newvs;
    data->l_mark = newvs + i;

    data->stacksize = newsize;
    data->s_last = data->s_base + newsize - 1;
    return 0;
}

#if YYPURE || defined(YY_NO_LEAKS)
static void yyfreestack(YYSTACKDATA *data)
{
    free(data->s_base);
    free(data->l_base);
    memset(data, 0, sizeof(*data));
}
#else
#define yyfreestack(data) /* nothing */
#endif

#define YYABORT  goto yyabort
#define YYREJECT goto yyabort
#define YYACCEPT goto yyaccept
#define YYERROR  goto yyerrlab

int
YYPARSE_DECL()
{
    int yym, yyn, yystate;
#if YYDEBUG
    const char *yys;

    if ((yys = getenv("YYDEBUG")) != 0)
    {
        yyn = *yys;
        if (yyn >= '0' && yyn <= '9')
            yydebug = yyn - '0';
    }
#endif

    yynerrs = 0;
    yyerrflag = 0;
    yychar = YYEMPTY;
    yystate = 0;

#if YYPURE
    memset(&yystack, 0, sizeof(yystack));
#endif

    if (yystack.s_base == NULL && yygrowstack(&yystack)) goto yyoverflow;
    yystack.s_mark = yystack.s_base;
    yystack.l_mark = yystack.l_base;
    yystate = 0;
    *yystack.s_mark = 0;

yyloop:
    if ((yyn = yydefred[yystate]) != 0) goto yyreduce;
    if (yychar < 0)
    {
        if ((yychar = YYLEX) < 0) yychar = 0;
#if YYDEBUG
        if (yydebug)
        {
            yys = yyname[YYTRANSLATE(yychar)];
            printf("%sdebug: state %d, reading %d (%s)\n",
                    YYPREFIX, yystate, yychar, yys);
        }
#endif
    }
    if ((yyn = yysindex[yystate]) && (yyn += yychar) >= 0 &&
            yyn <= YYTABLESIZE && yycheck[yyn] == yychar)
    {
#if YYDEBUG
        if (yydebug)
            printf("%sdebug: state %d, shifting to state %d\n",
                    YYPREFIX, yystate, yytable[yyn]);
#endif
        if (yystack.s_mark >= yystack.s_last && yygrowstack(&yystack))
        {
            goto yyoverflow;
        }
        yystate = yytable[yyn];
        *++yystack.s_mark = yytable[yyn];
        *++yystack.l_mark = yylval;
        yychar = YYEMPTY;
        if (yyerrflag > 0)  --yyerrflag;
        goto yyloop;
    }
    if ((yyn = yyrindex[yystate]) && (yyn += yychar) >= 0 &&
            yyn <= YYTABLESIZE && yycheck[yyn] == yychar)
    {
        yyn = yytable[yyn];
        goto yyreduce;
    }
    if (yyerrflag) goto yyinrecovery;

    yyerror("syntax error");

    goto yyerrlab;

yyerrlab:
    ++yynerrs;

yyinrecovery:
    if (yyerrflag < 3)
    {
        yyerrflag = 3;
        for (;;)
        {
            if ((yyn = yysindex[*yystack.s_mark]) && (yyn += YYERRCODE) >= 0 &&
                    yyn <= YYTABLESIZE && yycheck[yyn] == YYERRCODE)
            {
#if YYDEBUG
                if (yydebug)
                    printf("%sdebug: state %d, error recovery shifting\
 to state %d\n", YYPREFIX, *yystack.s_mark, yytable[yyn]);
#endif
                if (yystack.s_mark >= yystack.s_last && yygrowstack(&yystack))
                {
                    goto yyoverflow;
                }
                yystate = yytable[yyn];
                *++yystack.s_mark = yytable[yyn];
                *++yystack.l_mark = yylval;
                goto yyloop;
            }
            else
            {
#if YYDEBUG
                if (yydebug)
                    printf("%sdebug: error recovery discarding state %d\n",
                            YYPREFIX, *yystack.s_mark);
#endif
                if (yystack.s_mark <= yystack.s_base) goto yyabort;
                --yystack.s_mark;
                --yystack.l_mark;
            }
        }
    }
    else
    {
        if (yychar == 0) goto yyabort;
#if YYDEBUG
        if (yydebug)
        {
            yys = yyname[YYTRANSLATE(yychar)];
            printf("%sdebug: state %d, error recovery discards token %d (%s)\n",
                    YYPREFIX, yystate, yychar, yys);
        }
#endif
        yychar = YYEMPTY;
        goto yyloop;
    }

yyreduce:
#if YYDEBUG
    if (yydebug)
        printf("%sdebug: state %d, reducing by rule %d (%s)\n",
                YYPREFIX, yystate, yyn, yyrule[yyn]);
#endif
    yym = yylen[yyn];
    if (yym)
        yyval = yystack.l_mark[1-yym];
    else
        memset(&yyval, 0, sizeof yyval);
    switch (yyn)
    {
case 7:
#line 52 "ergo_input_processor.y"
	{ if(!ergo_scanner_reading_stdin) { yyerror("Aborted."); YYABORT; } }
break;
case 8:
#line 56 "ergo_input_processor.y"
	{ es_assign_num(yystack.l_mark[-2].var, yystack.l_mark[0].num);}
break;
case 9:
#line 57 "ergo_input_processor.y"
	{ es_assign_str(yystack.l_mark[-2].var, yystack.l_mark[0].str);}
break;
case 10:
#line 61 "ergo_input_processor.y"
	{ yyval.var=es_find_var(NULL, yystack.l_mark[0].str); 
                              if(!yyval.var) { last_token = yystack.l_mark[0].str;
                                        yyerror("Unknown variable");YYERROR; }}
break;
case 11:
#line 64 "ergo_input_processor.y"
	{ yyval.var=es_find_var(yystack.l_mark[-2].var,   yystack.l_mark[0].str);
                              if(!yyval.var) { last_token = yystack.l_mark[0].str;
                                        yyerror("Unknown variable");YYERROR;}}
break;
case 12:
#line 70 "ergo_input_processor.y"
	{ es_mol_commit(); }
break;
case 13:
#line 71 "ergo_input_processor.y"
	{ es_mol_commit(); }
break;
case 14:
#line 72 "ergo_input_processor.y"
	{ es_mol_commit(); }
break;
case 15:
#line 73 "ergo_input_processor.y"
	{ es_mol_commit(); }
break;
case 16:
#line 74 "ergo_input_processor.y"
	{ es_print_help(); }
break;
case 17:
#line 75 "ergo_input_processor.y"
	{ es_print_help_var(yystack.l_mark[0].var); }
break;
case 18:
#line 76 "ergo_input_processor.y"
	{ es_print_precision(); }
break;
case 19:
#line 77 "ergo_input_processor.y"
	{ if(es_mol_read_molecule(yystack.l_mark[0].str,MOL_MAIN)) { 
                              yyerror("Reading MOLECULE failed"); YYERROR;} }
break;
case 20:
#line 79 "ergo_input_processor.y"
	{ if(es_mol_read_molecule(yystack.l_mark[0].str, MOL_GHOST)) { 
                              yyerror("Reading GHOST MOLECULE failed"); YYERROR;} }
break;
case 21:
#line 81 "ergo_input_processor.y"
	{ if(es_run(yystack.l_mark[0].str, 0)) {
                              yyerror("RUN failed"); YYERROR;} }
break;
case 22:
#line 83 "ergo_input_processor.y"
	{
          if(es_set_nthreads(yystack.l_mark[-1].num)) { yyerror("setNThreads failed"); YYERROR;} }
break;
case 23:
#line 85 "ergo_input_processor.y"
	{
          if(es_set_nthreads_string(yystack.l_mark[-1].str)) { yyerror("setNThreads failed"); YYERROR;} }
break;
case 24:
#line 87 "ergo_input_processor.y"
	{ puts(yystack.l_mark[0].str); 
                   if(system(yystack.l_mark[0].str) != 0) {yyerror("system() failed"); YYERROR;} }
break;
case 25:
#line 89 "ergo_input_processor.y"
	{ es_warranty(); }
break;
case 26:
#line 90 "ergo_input_processor.y"
	{ if(es_getexc(yystack.l_mark[-1].str, yystack.l_mark[0].num)) {
                              yyerror("get_excited_state failed"); YYERROR;} }
break;
case 27:
#line 92 "ergo_input_processor.y"
	{ if(es_get_polarisability(yystack.l_mark[-2].str, yystack.l_mark[-1].str, yystack.l_mark[0].num)) {
                              yyerror("get_polarisability failed"); YYERROR;} }
break;
case 28:
#line 94 "ergo_input_processor.y"
	{ if(es_get_polarisability(yystack.l_mark[-2].str, NULL,yystack.l_mark[0].num)) {
                              yyerror("get_polarisability failed"); YYERROR;} }
break;
case 29:
#line 96 "ergo_input_processor.y"
	{ YYACCEPT; }
break;
case 30:
#line 97 "ergo_input_processor.y"
	{
                         if(!es_assign_range(MOL_MAIN,yystack.l_mark[-4].num, yystack.l_mark[-2].num, yystack.l_mark[-1].num, yystack.l_mark[0].str)) {
                            yyerror("Invalid main basis set range");YYERROR;}
                         }
break;
case 31:
#line 101 "ergo_input_processor.y"
	{
                         if(!es_assign_range(MOL_GHOST,yystack.l_mark[-4].num, yystack.l_mark[-2].num, yystack.l_mark[-1].num, yystack.l_mark[0].str)) {
                            yyerror("Invalid ghost basis set range");YYERROR;}
                         }
break;
case 34:
#line 113 "ergo_input_processor.y"
	{ es_add_atom(yystack.l_mark[-4].str, yystack.l_mark[-3].num, yystack.l_mark[-2].num, yystack.l_mark[-1].num); }
break;
case 35:
#line 117 "ergo_input_processor.y"
	{ yyval.num=yystack.l_mark[0].num;  }
break;
case 36:
#line 118 "ergo_input_processor.y"
	{ yyval.num=yystack.l_mark[-2].num+yystack.l_mark[0].num; }
break;
case 37:
#line 119 "ergo_input_processor.y"
	{ yyval.num=yystack.l_mark[-2].num-yystack.l_mark[0].num; }
break;
case 38:
#line 120 "ergo_input_processor.y"
	{ yyval.num=yystack.l_mark[-2].num*yystack.l_mark[0].num; }
break;
case 39:
#line 121 "ergo_input_processor.y"
	{ yyval.num=yystack.l_mark[-2].num/yystack.l_mark[0].num; }
break;
case 40:
#line 122 "ergo_input_processor.y"
	{ yyval.num=-yystack.l_mark[0].num; }
break;
case 41:
#line 123 "ergo_input_processor.y"
	{ yyval.num=pow(yystack.l_mark[-2].num,yystack.l_mark[0].num); }
break;
case 42:
#line 124 "ergo_input_processor.y"
	{ yyval.num=yystack.l_mark[-1].num; }
break;
#line 713 "ergo_input_processor.c"
    }
    yystack.s_mark -= yym;
    yystate = *yystack.s_mark;
    yystack.l_mark -= yym;
    yym = yylhs[yyn];
    if (yystate == 0 && yym == 0)
    {
#if YYDEBUG
        if (yydebug)
            printf("%sdebug: after reduction, shifting from state 0 to\
 state %d\n", YYPREFIX, YYFINAL);
#endif
        yystate = YYFINAL;
        *++yystack.s_mark = YYFINAL;
        *++yystack.l_mark = yyval;
        if (yychar < 0)
        {
            if ((yychar = YYLEX) < 0) yychar = 0;
#if YYDEBUG
            if (yydebug)
            {
                yys = yyname[YYTRANSLATE(yychar)];
                printf("%sdebug: state %d, reading %d (%s)\n",
                        YYPREFIX, YYFINAL, yychar, yys);
            }
#endif
        }
        if (yychar == 0) goto yyaccept;
        goto yyloop;
    }
    if ((yyn = yygindex[yym]) && (yyn += yystate) >= 0 &&
            yyn <= YYTABLESIZE && yycheck[yyn] == yystate)
        yystate = yytable[yyn];
    else
        yystate = yydgoto[yym];
#if YYDEBUG
    if (yydebug)
        printf("%sdebug: after reduction, shifting from state %d \
to state %d\n", YYPREFIX, *yystack.s_mark, yystate);
#endif
    if (yystack.s_mark >= yystack.s_last && yygrowstack(&yystack))
    {
        goto yyoverflow;
    }
    *++yystack.s_mark = (short) yystate;
    *++yystack.l_mark = yyval;
    goto yyloop;

yyoverflow:
    yyerror("yacc stack overflow");

yyabort:
    yyfreestack(&yystack);
    return (1);

yyaccept:
    yyfreestack(&yystack);
    return (0);
}
