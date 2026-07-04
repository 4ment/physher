// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef PhyC_combinatorics_h
#define PhyC_combinatorics_h

double bell_number( unsigned int n );

int choose ( unsigned int n, unsigned int k );

void combination_at_index( unsigned int n,  unsigned int k, long m, unsigned *ans );

#endif
