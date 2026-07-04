// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef checkpoint_h
#define checkpoint_h

#include "parameters.h"

void checkpoint_apply(const char* file_path, Parameters* parameters);

void checkpoint_save(const char* file_path, Parameters* parameters);

#endif /* checkpoint_h */