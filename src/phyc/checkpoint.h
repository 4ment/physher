// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef checkpoint_h
#define checkpoint_h

#include "parameters.h"

// Restore parameter values from a JSON checkpoint written by checkpoint_save.
// Parameters are matched by name; any parameter absent from the file, of the
// wrong shape, or with a mismatched dimension is left untouched (a message is
// printed). See checkpoint.c for the on-disk format.
void checkpoint_apply(const char* file_path, Parameters* parameters);

// Atomically write all parameter values to a JSON checkpoint file (written to a
// temporary sibling then renamed, so an existing checkpoint is never left
// partially overwritten).
void checkpoint_save(const char* file_path, Parameters* parameters);

#endif /* checkpoint_h */