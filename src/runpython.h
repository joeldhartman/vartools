#define VARTOOLS_PYTHON_MESSAGE_ENDPROCESS 0
#define VARTOOLS_PYTHON_MESSAGE_READDATA 1


typedef struct {
  void *dataptr;
  int datatype;
} _PythonArrayData;

typedef struct {
  int Nfunc;
  PyObject *CompiledUserCode;
  PyObject *UserModule;
  PyObject **UserFunctionToRun;
  PyObject **Variables;
  PyObject *VariableList;
  int Nvars;
  _PythonArrayData *data;
} _PythonObjectContainer;

