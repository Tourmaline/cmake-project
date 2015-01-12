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
extern YYSTYPE yylval;
