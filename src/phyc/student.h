// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef _STUDENT_DISTRIBUTION_H_
#define _STUDENT_DISTRIBUTION_H_

/**
 * Returns the quantile for the two-tailed student's t distribution.
 * This is an implementation of the algorithm in
 * G. W. Hill. "Algorithm 396: Student's t-Quantiles." Communications
 * of the ACM 13(10):619--620.  ACM Press, October, 1970.
 */

double qt(double p, long n);

#endif
