#ifndef ARBEVAL8_H
#define ARBEVAL8_H

#include <red_black_tree.h>
#include <doubly_linked_list.h>
#include <math.h>

#if defined(_WIN32) || defined(_WIN64)
#define ARBEVAL_EXPORTS 1
#ifdef ARBEVAL_EXPORTS
#define ARBEVAL_API __declspec(dllexport)
#else
#define ARBEVAL_API __declspec(dllimport)
#endif
#else
#define ARBEVAL_API
#endif



typedef struct arbenv arbenv;

extern void *arbeval_globals;

ARBEVAL_API char *GetLocalPath();
ARBEVAL_API char *GetLocalOutPath();
ARBEVAL_API void SetPlatformCode(long platcode);
ARBEVAL_API void SetLocalPath(char *lpath);
ARBEVAL_API void SetLocalOutPath(char *lpath);
ARBEVAL_API void SetPlatformName(char *lpath);
ARBEVAL_API void SetInstructions(char *init,long isFile);
ARBEVAL_API void SetGraphingInstructions(char *init,long isFile);
ARBEVAL_API void SetFourierInstructions(char *init,long isFile);
ARBEVAL_API void SetProgrammingInstructions(char *init,long isFile);
ARBEVAL_API void SetInitScript(char *init,long isFile);
ARBEVAL_API char *ConstructReadObjectAssignmentFromFile(char *fname_s);
ARBEVAL_API void AddArbevalEnvironment(arbenv *state, arbenv *ae, int overwrite);
ARBEVAL_API char *ArbevalHistoryToString(arbenv *ae);
ARBEVAL_API void ClearArbevalHistoryEntry(arbenv *ae, long eind);
ARBEVAL_API void ClearArbevalHistory(arbenv *ae);
ARBEVAL_API void WriteArbevalSettingsHeader(FILE *ff, char *hname, double version);
ARBEVAL_API void WriteArbevalSettingsHeaderToStream(void *ff, char *hname, double version,
  stream_type st);
ARBEVAL_API rbtree *LoadFontFromStream(MEM_FILE *file);
ARBEVAL_API void AppendHistoryItem(arbenv *ae, char *in, char *out);
ARBEVAL_API void SetLastHistoryItemOutput(arbenv *ae, char *out);
ARBEVAL_API unsigned char *DummyBitmap(int cols, int rows);
ARBEVAL_API uint32_t **LoadBitmapToRGBMatrix(char *fname);
ARBEVAL_API  uint32_t **DownSampleRGBMatrix(uint32_t **rgbmat, long dsfact);
ARBEVAL_API void RemoveHistoryItem(arbenv *ae, long indx);
ARBEVAL_API char SetUpdateHistoryFlag(arbenv *ae, char val);
ARBEVAL_API size_t ArbenvSize();
ARBEVAL_API long AtoL(char *nums);
ARBEVAL_API double AtoF(char *nums);
ARBEVAL_API long TestNumericalString(char *, double *val);
ARBEVAL_API void FreeVariable(void *val);
ARBEVAL_API char **GetPredefinedFunctionInfo(char *flabel);
ARBEVAL_API char **GetEnvironmentFunctionInfo(char *flabel);
ARBEVAL_API char **GetProgrammingCommandInfo(char *flabel);
ARBEVAL_API arbenv *InitializeTopLevelArbevalEnvironment(long defaults);
ARBEVAL_API void TerminateArbevalEnvironment(arbenv *ae);
ARBEVAL_API void *copyArbenv(void *v);
ARBEVAL_API void freeArbenv(void *ae);
ARBEVAL_API void writeArbenv(void *v, void *f, stream_type st);
ARBEVAL_API char *FormatRBST(char **sin, long expand, long sd, char * prevnum, char * prevdenom, arbenv *ae);
ARBEVAL_API void *readTopArbenv(void *f, stream_type st);
ARBEVAL_API char *ArbevalScript(char *script1, arbenv *ae, long delete_literal_objects);
ARBEVAL_API char *RestoreDefaultAssignments(arbenv *ae, long delete_lits);
ARBEVAL_API char *PrintEnvironment(arbenv *ae, long indent, long print_gen);
ARBEVAL_API void TerminateArbevalGlobals(void);
ARBEVAL_API void InitializeArbevalGlobals(char *init, long init_is_file);
ARBEVAL_API void ReinitializeHistory(arbenv *ae);
ARBEVAL_API void HistoryIteratorSeekEnd(arbenv *ae);
ARBEVAL_API dl_list *GetHistory(arbenv *ae);
ARBEVAL_API double ArbevalSettingsFileVersionCode(char *fname, size_t offset);
ARBEVAL_API double ArbevalSettingsFileVersionCodeFromStream(void *f1, stream_type st);
ARBEVAL_API char *Instructions(void);
ARBEVAL_API char *InstructionsGraphing(void);
ARBEVAL_API char *InstructionsProgramming(void);
ARBEVAL_API char *InstructionsFourier(void);
ARBEVAL_API char *About(void);
ARBEVAL_API long GetSignificantDigits(arbenv *ae);
ARBEVAL_API long GetCancelFlag(arbenv *ae);
ARBEVAL_API void SetSignificantDigits(arbenv *ae, long nsig );
ARBEVAL_API void SetCancelFlag(arbenv *ae, char val);
ARBEVAL_API arbenv *CopyArbevalEnvironment(arbenv *source, arbenv *outer,
  long cpyhist, long ischild);
ARBEVAL_API void SetSeed(arbenv *ae, long val);
ARBEVAL_API char *ReadFileToString(char *path);
ARBEVAL_API double GetArbevalVersion();
ARBEVAL_API char **GetPredefinedFunctionLabels();
ARBEVAL_API char **GetEnvironmentFunctionLabels();
ARBEVAL_API char **GetProgrammingCommandLabels();
ARBEVAL_API long GetProgrammingCommandListSize();
ARBEVAL_API long GetPredefinedFunctionListSize();
ARBEVAL_API long GetEnvironmentFunctionListSize();
ARBEVAL_API char *GetPlatformName();


ARBEVAL_API unsigned char **GeneratePlotBitmaps(char *expr, double xmin, double xmax, uint32_t nx,
  cmplx *ymin, cmplx *ymax,
  uint32_t dimx, uint32_t dimy, rbtree *font,
  arbenv *ae);
ARBEVAL_API unsigned char *GenerateParametricPlotBitmap(char *expr11, double tmin,
  double tmax, long nt, cmplx ymin,
  cmplx ymax, uint32_t dimx,
  uint32_t dimy, rbtree *font,
  arbenv *ae);
ARBEVAL_API unsigned char **GenerateFourierPlotBitmaps(char *expr11, double xmin, double xmax, long pow2,
  cmplx yminsig, cmplx ymaxsig, cmplx yminspec, cmplx ymaxspec,
  uint32_t dimx, uint32_t dimy,
  double ***ft, char **dnames,
  rbtree *font, arbenv *ae);
ARBEVAL_API unsigned char **GenerateSpectrallyFilteredPlotBitmaps(char *expr11,
  double **tf1d1,
  double xmin, double xmax, uint32_t dimx,
  cmplx ymin, cmplx ymax,
  uint32_t dimy,
  rbtree *font, arbenv *ae);
ARBEVAL_API unsigned char **GeneratePlotBitmaps2D(char *expr11, double xmin, double xmax,
  double ymin, double ymax, uint32_t dimx,
  uint32_t dimy, uint32_t nx, uint32_t ny, rbtree *font,
  arbenv *ae);
ARBEVAL_API unsigned char **GeneratePlotBitmapsComplex(char *expr11, double xmin,
  double xmax, double ymin, double ymax, uint32_t dimx,
  uint32_t dimy, uint32_t nx, uint32_t ny,
  rbtree *font, arbenv *ae);
ARBEVAL_API rbtree *LoadFontFromFile(char *file);
ARBEVAL_API void FreeFont(rbtree *font);
ARBEVAL_API void writeRBNodeValueBitmap(void *val, void *f, stream_type st);
ARBEVAL_API void *copyRBNodeValueBitmap(void *v);
ARBEVAL_API void *readRBNodeValueBitmap(void *f, stream_type st); 
ARBEVAL_API void freeRBNodeValueBitmap(void *v);
ARBEVAL_API unsigned char *getBitmap(char *key, rbtree *rbt);


#endif // ARBEVAL8_H
