/*
 *  lh.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 10/12/2015.
 *  Copyright (C) 2016 Mathieu Fourment. All rights reserved.
 *
 *  This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along with this program; if not,
 *  write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

#ifndef PhyC_lg_h
#define PhyC_lg_h

#include "substmodel.h"

static const double AMINO_ACID_MODEL_LG[20][20] = {
    {             0,     2.4890837,    0.39514435,     1.0385443,    0.25370059,      2.066039,    0.35885784,    0.14982985,    0.53651783,    0.39533733,     1.1240347,    0.27681832,     1.1776509,    0.96989392,    0.42509245,     4.7271804,     2.1395006,     2.5478697,    0.18071662,    0.21895908},
    {     2.4890837,             0,   0.062555906,   0.003499274,     1.1052511,    0.56926458,    0.64054259,    0.32062657,   0.013265853,    0.59400724,    0.89367979,    0.52876813,   0.075381891,   0.084808359,    0.53455088,     2.7844772,     1.1434797,       1.95929,    0.67012807,     1.1655316},
    {    0.39514435,   0.062555906,             0,      5.243868,    0.01741582,     0.8449254,     0.9271141,    0.01069023,     0.2829585,   0.015075551,   0.025548356,     5.0761471,    0.39445605,    0.52338552,    0.12395413,     1.2402744,    0.42585978,   0.037967113,   0.029890093,    0.13510725},
    {     1.0385443,   0.003499274,      5.243868,             0,   0.018811039,    0.34884696,    0.42388087,   0.044265021,     1.8071764,   0.069672559,    0.17373539,    0.54171192,    0.41940891,     4.1285897,    0.36396992,    0.61197284,    0.60454457,    0.24503411,   0.077852072,    0.12003742},
    {    0.25370059,     1.1052511,    0.01741582,   0.018811039,             0,   0.089585483,    0.68213848,     1.1127262,   0.023918071,     2.5926917,     1.7988525,   0.089525087,   0.094463956,   0.035855056,   0.052721535,    0.36181849,    0.16500105,    0.65468292,     2.4571204,     7.8038994},
    {      2.066039,    0.56926458,     0.8449254,    0.34884696,   0.089585483,             0,    0.31148391,  0.0087054516,    0.29663617,    0.04426053,    0.13953789,     1.4376449,    0.19696132,    0.26795876,    0.39019229,     1.7399898,    0.12983619,    0.07670085,    0.26849133,   0.054678753},
    {    0.35885784,    0.64054259,     0.9271141,    0.42388087,    0.68213848,    0.31148391,             0,    0.10888195,    0.69726378,     0.3663167,    0.44247196,     4.5092364,    0.50885055,     4.8135036,     2.4266008,     0.9900121,    0.58426218,    0.11901337,    0.59705339,     5.3068322},
    {    0.14982985,    0.32062657,    0.01069023,   0.044265021,     1.1127262,  0.0087054516,    0.10888195,             0,    0.15906847,     4.1450659,     4.2736062,    0.19150282,    0.07828126,   0.072854177,    0.12699054,   0.064105457,     1.0337382,     10.649104,    0.11166008,    0.23252269},
    {    0.53651783,   0.013265853,     0.2829585,     1.8071764,   0.023918071,    0.29663617,    0.69726378,    0.15906847,             0,    0.13749964,    0.65660428,     2.1450773,    0.39032178,     3.2342926,     6.3260657,    0.74868243,     1.1368625,    0.18520165,   0.049905553,    0.13193188},
    {    0.39533733,    0.59400724,   0.015075551,   0.069672559,     2.5926917,    0.04426053,     0.3663167,     4.1450659,    0.13749964,             0,      6.312356,    0.06842664,    0.24906001,    0.58245637,    0.30184773,     0.1822869,    0.30293571,      1.702745,    0.61963194,    0.29964768},
    {     1.1240347,    0.89367979,   0.025548356,    0.17373539,     1.7988525,    0.13953789,    0.44247196,     4.2736062,    0.65660428,      6.312356,             0,    0.37100375,   0.099848702,     1.6725683,    0.48413321,    0.34695977,     2.0203655,     1.8987176,    0.69617497,    0.48130595},
    {    0.27681832,    0.52876813,     5.0761471,    0.54171192,   0.089525087,     1.4376449,     4.5092364,    0.19150282,     2.1450773,    0.06842664,    0.37100375,             0,    0.16178707,     1.6957511,    0.75187754,     4.0083573,     2.0006785,    0.08368748,   0.045375965,    0.61202438},
    {     1.1776509,   0.075381891,    0.39445605,    0.41940891,   0.094463956,    0.19696132,    0.50885055,    0.07828126,    0.39032178,    0.24906001,   0.099848702,    0.16178707,             0,    0.62429374,    0.33253337,     1.3381318,    0.57146773,    0.29650076,   0.095130632,   0.089613415},
    {    0.96989392,   0.084808359,    0.52338552,     4.1285897,   0.035855056,    0.26795876,     4.8135036,   0.072854177,     3.2342926,    0.58245637,     1.6725683,     1.6957511,    0.62429374,             0,      2.807907,     1.2238279,      1.080136,    0.21033229,    0.23619896,    0.25733597},
    {    0.42509245,    0.53455088,    0.12395413,    0.36396992,   0.052721535,    0.39019229,     2.4266008,    0.12699054,     6.3260657,    0.30184773,    0.48413321,    0.75187754,    0.33253337,      2.807907,             0,    0.85815097,    0.57898665,    0.17088738,     0.5936072,    0.31443952},
    {     4.7271804,     2.7844772,     1.2402744,    0.61197284,    0.36181849,     1.7399898,     0.9900121,   0.064105457,    0.74868243,     0.1822869,    0.34695977,     4.0083573,     1.3381318,     1.2238279,    0.85815097,             0,     6.4722771,   0.098368473,    0.24886188,    0.40054664},
    {     2.1395006,     1.1434797,    0.42585978,    0.60454457,    0.16500105,    0.12983619,    0.58426218,     1.0337382,     1.1368625,    0.30293571,     2.0203655,     2.0006785,    0.57146773,      1.080136,    0.57898665,     6.4722771,             0,     2.1881575,    0.14082482,    0.24584051},
    {     2.5478697,       1.95929,   0.037967113,    0.24503411,    0.65468292,    0.07670085,    0.11901337,     10.649104,    0.18520165,      1.702745,     1.8987176,    0.08368748,    0.29650076,    0.21033229,    0.17088738,   0.098368473,     2.1881575,             0,    0.18950974,    0.24931339},
    {    0.18071662,    0.67012807,   0.029890093,   0.077852072,     2.4571204,    0.26849133,    0.59705339,    0.11166008,   0.049905553,    0.61963194,    0.69617497,   0.045375965,   0.095130632,    0.23619896,     0.5936072,    0.24886188,    0.14082482,    0.18950974,             0,     3.1518138},
    {    0.21895908,     1.1655316,    0.13510725,    0.12003742,     7.8038994,   0.054678753,     5.3068322,    0.23252269,    0.13193188,    0.29964768,    0.48130595,    0.61202438,   0.089613415,    0.25733597,    0.31443952,    0.40054664,    0.24584051,    0.24931339,     3.1518138,             0}
    
};

static const double AMINO_ACID_MODEL_LG_FREQUENCIES[20] = {0.079066, 0.0129369, 0.0530516, 0.0715863, 0.0423017, 0.0573372, 0.0223546, 0.0621565, 0.0646003, 0.099081, 0.0229506, 0.0419774, 0.0440395, 0.0407668, 0.0559413, 0.0611971, 0.0532871, 0.0691469, 0.0120656, 0.0341554};


SubstitutionModel * new_LG();

SubstitutionModel * new_LG_with_values( const double *freqs );

#endif
