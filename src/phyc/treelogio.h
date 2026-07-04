// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef PhyC_treelogio_h
#define PhyC_treelogio_h

int TreeLog_define_range( const char *filename, int *start, int *end, double alpha );

void TreeLog_define_range_nclasses( const char *filename, int nclasses, int *start, int *end );

#endif
