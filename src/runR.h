#define VARTOOLS_R_MESSAGE_ENDPROCESS 0
#define VARTOOLS_R_MESSAGE_READDATA 1


typedef struct {
  void *dataptr;
  int datatype;
} _RArrayData;

typedef struct {
  int Nfunc;
  int Nvars;
  _RArrayData *data;
} _RObjectContainer;
