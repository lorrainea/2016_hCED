/**
    hCED: Heuristic Cyclic Edit Distance
    Copyright (C) 2016 Solon P. Pissis, Lorraine A. K. Ayad 

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/

/**
This matrix was created by Todd Lowe 12/10/92
Uses ambiguous nucleotide codes, probabilities rounded to
nearest integer
Lowest score = -4, Highest score = 5
  A T G C S W R Y K M B V H D N
A 5 -4 -4 -4 -4 1 1 -4 -4 1 -4 -1 -1 -1 -2
T -4 5 -4 -4 -4 1 -4 1 1 -4 -1 -4 -1 -1 -2
G -4 -4 5 -4 1 -4 1 -4 1 -4 -1 -1 -4 -1 -2
C -4 -4 -4 5 1 -4 -4 1 -4 1 -1 -1 -1 -4 -2
S -4 -4 1 1 -1 -4 -2 -2 -2 -2 -1 -1 -3 -3 -1
W 1 1 -4 -4 -4 -1 -2 -2 -2 -2 -3 -3 -1 -1 -1
R 1 -4 1 -4 -2 -2 -1 -4 -2 -2 -3 -1 -3 -1 -1
Y -4 1 -4 1 -2 -2 -4 -1 -2 -2 -1 -3 -1 -3 -1
K -4 1 1 -4 -2 -2 -2 -2 -1 -4 -1 -3 -3 -1 -1
M 1 -4 -4 1 -2 -2 -2 -2 -4 -1 -3 -1 -1 -3 -1
B -4 -1 -1 -1 -1 -3 -3 -1 -1 -3 -1 -2 -2 -2 -1
V -1 -4 -1 -1 -1 -3 -1 -3 -3 -1 -2 -1 -2 -2 -1
H -1 -1 -4 -1 -3 -1 -3 -1 -3 -1 -2 -2 -1 -2 -1
D -1 -1 -1 -4 -3 -1 -1 -3 -1 -3 -2 -2 -2 -1 -1
N -2 -2 -2 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
Availability: <ftp://ftp.ncbi.nih.gov/blast/matrices/NUC.4.4>
**/

double EDNAFULL_matrix[15][15] = {
{ 5, -4, -4, -4, -4, 1, 1, -4, -4, 1, -4, -1, -1, -1, -2},
{-4, 5, -4, -4, -4, 1, -4, 1, 1, -4, -1, -4, -1, -1, -2},
{-4, -4, 5, -4, 1, -4, 1, -4, 1, -4, -1, -1, -4, -1, -2},
{-4, -4, -4, 5, 1, -4, -4, 1, -4, 1, -1, -1, -1, -4, -2},
{-4, -4, 1, 1, -1, -4, -2, -2, -2, -2, -1, -1, -3, -3, -1},
{1, 1, -4, -4, -4, -1, -2, -2, -2, -2, -3, -3, -1, -1, -1},
{1, -4, 1, -4, -2, -2, -1, -4, -2, -2, -3, -1, -3, -1, -1},
{-4, 1, -4, 1, -2, -2, -4, -1, -2, -2, -1, -3, -1, -3, -1},
{-4, 1, 1, -4, -2, -2, -2, -2, -1, -4, -1, -3, -3, -1, -1},
{1, -4, -4, 1, -2, -2, -2, -2, -4, -1, -3, -1, -1, -3, -1},
{-4, -1, -1, -1, -1, -3, -3, -1, -1, -3, -1, -2, -2, -2, -1},
{-1, -4, -1, -1, -1, -3, -1, -3, -3, -1, -2, -1, -2, -2, -1},
{-1, -1, -4, -1, -3, -1, -3, -1, -3, -1, -2, -2, -1, -2, -1},
{-1, -1, -1, -4, -3, -1, -1, -3, -1, -3, -2, -2, -2, -1, -1},
{-2, -2, -2, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1} };
