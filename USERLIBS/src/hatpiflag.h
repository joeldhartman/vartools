typedef struct {
  char ***fiphotflagvals;
  double **rejbadframemaskvals;
  double **tfaoutliermaskvals;
  double **pointingoutlierflagvals;
  short **outflagvals;

  _Variable *fiphotflagvar;
  _Variable *rejbadframemaskvar;
  _Variable *tfaoutliermaskvar;
  _Variable *pointingoutlierflagvar;
  _Variable *outflagvar;

} _Hatpiflag;
