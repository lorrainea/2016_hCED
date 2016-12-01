/**
    CSC: Circular Sequence Comparison
    Copyright (C) 2015 Solon P. Pissis, Ahmad Retha, Fatima Vayani 

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

#ifndef __CSC__
#define __CSC__

#define ALLOC_SIZE              1048576
#define DEL                     '$'
#define DEL_STR                 "$"
#define no_edit_distance	"N"
#define edit_distance		"Y"

#define MAX2(a,b) ((a) > (b)) ? (a) : (b)
#define MIN2(a,b) ((a) < (b)) ? (a) : (b)

struct TSwitch
{
    char *               input_filename;         // the input file name
    char *               output_filename;        // the output file name
    unsigned int         l;                      // block length (min. required)
    unsigned int         q;                      // q-gram size (min. required)
    double               P;                      // Percent Sequence to align at ends
    char * 		 e;			 // edit distance method
    unsigned int         S, I, D;                // EDIT DISTANCE costs for substitution, insertion, deletion
    int                  m, r, f, g, O, E;       // SIMILARITY costs for edit distance match, substitution, insertion, deletion, gap open, gap extend
    int 		 R;			 // computes edit distance for x,y and y,x
};

double gettime( void );
int decode_switches ( int argc, char * argv [], struct TSwitch * sw , char * eD );
void usage ( void );
void create_rotation ( unsigned char * x, unsigned int offset, unsigned char * rotation );
void create_backward_rotation ( unsigned char * x, unsigned int offset, unsigned char * rotation );

#endif
