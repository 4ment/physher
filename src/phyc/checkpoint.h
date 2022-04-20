
#ifndef checkpoint_h
#define checkpoint_h

#include "parameters.h"

void checkpoint_apply(const char* file_path, Parameters* parameters);

void checkpoint_save(const char* file_path, Parameters* parameters);

#endif /* checkpoint_h */