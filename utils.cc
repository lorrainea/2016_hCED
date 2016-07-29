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

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include <limits.h>
#include <sys/time.h>
#include "hced.h"


static struct option long_options[] =
 {
   { "input-file",              required_argument, NULL, 'i' },
   { "output-file",             required_argument, NULL, 'o' },
   { "q-length",            	required_argument, NULL, 'q' },
   { "block-length",        	required_argument, NULL, 'l' },
   { "refine-blocks",          	optional_argument, NULL, 'P' },
   { "edit-distance",           optional_argument, NULL, 'e' },
   { "cost-substitution",	optional_argument, NULL, 'S' },
   { "cost-insertion",		optional_argument, NULL, 'I' },
   { "cost-deletion",		optional_argument, NULL, 'D' },
   { "score-match",     	optional_argument, NULL, 'm' },
   { "score-mismatch",     	optional_argument, NULL, 'r' },
   { "score-insertion",     	optional_argument, NULL, 'f' },
   { "score-deletion",     	optional_argument, NULL, 'g' },
   { "gap-open",     		optional_argument, NULL, 'O' },
   { "gap-extend",     		optional_argument, NULL, 'E' },
   { "help",                    no_argument,       NULL, 'h' },
   {  NULL,                     0,                 NULL,  0  }
 };


/* 
Decode the input switches 
*/
int decode_switches ( int argc, char * argv [], struct TSwitch * sw, char * eD )
 {
   int          oi;
   int          opt;
   double       val;
   char       * ep;
   int          args;

   /* initialisation */
   sw -> input_filename                 = NULL;
   sw -> output_filename                = NULL;
   sw -> q                              = 5;
   sw -> l                              = 50;
   sw -> P                              = 1.0;
   sw -> e                              = eD;
   sw -> e[0]      			= 'Y';
   sw -> e[1]      			= '\0';
   sw -> S				= 1;
   sw -> I 				= 1;
   sw -> D				= 1;
   sw -> m				= 1;
   sw -> r				= -1;
   sw -> f				= -1;
   sw -> g				= -1;
   sw -> O				= 0;
   sw -> E				= 0;
   args = 0;

   while ( ( opt = getopt_long ( argc, argv, "i:o:q:l:P:e:S:I:D:m:r:f:g:O:E:h", long_options, &oi ) ) != - 1 )
    {
      switch ( opt )
       {
         case 'i':
           sw -> input_filename = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> input_filename, optarg );
           args ++;
           break;

         case 'o':
           sw -> output_filename = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> output_filename, optarg );
           args ++;
           break;

         case 'q':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> q = val;
           args ++;
           break;

         case 'l':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> l = val;
           args ++;
           break;

         case 'P':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> P = val;
           args ++;
           break;

	 case 'e':
           sw -> e = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> e, optarg );
           args ++;
           break;

	case 'S':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> S = val;
           break;

	case 'I':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> I = val;
           break;

	case 'D':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> D = val;
           break;

	case 'm':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> m = val;
           break;

	case 'r':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> r = val;
           break;

	case 'f':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> f = val;
           break;

	case 'g':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> g = val;
           break;

	case 'O':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> O = val;
           break;

	case 'E':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> E = val;
           break;

         case 'h':
           return ( 0 );
       }
    }

   if ( args < 4 )
     {
       usage ();
       exit ( 1 );
     }
   else
     return ( optind );
 }


/* 
Usage of the tool 
*/
void usage ( void )
 {
   fprintf ( stdout, " hCED <options>\n" );
   fprintf ( stdout, " Required arguments for saCSC:\n" );
   fprintf ( stdout, "  -i, --input-file            <str>     (Multi)FASTA input filename.\n" );
   fprintf ( stdout, "  -o, --output-file           <str>     Output filename for the rotated sequences.\n" );
   fprintf ( stdout, "  -q, --q-length              <int>     The q-gram length.\n");
   fprintf ( stdout, "  -l, --block-length          <int>     The length of each block.\n\n");
   fprintf ( stdout, " Optional for the Refinement stage:\n" );
   fprintf ( stdout, "  -P, --refine-blocks         <dbl>     Refine the alignment of saCSC by\n"
                     "                                        checking P blocks of the ends. Default: 1.\n" );
   fprintf ( stdout, "  -m, --score-match           <int>     Score of match for refinement. Default: 1.\n" );
   fprintf ( stdout, "  -r, --score-mismatch        <int>     Score of mismatch for refinement. Default: -1.\n" );
   fprintf ( stdout, "  -f, --score-insertion       <int>     Score of insertion for refinement. Default: -1.\n" );
   fprintf ( stdout, "  -g, --score-deletion        <int>     Score of deletion for refinement. Default: -1.\n" );
   fprintf ( stdout, "  -O, --gap-open              <int>     Score of gap opening for refinement. Default: NOT USED.\n" );
   fprintf ( stdout, "  -E, --gap-extend            <int>     Score of gap extension for refinement. Default: NOT USED.\n\n" );
   fprintf ( stdout, " Optional for edit distance model:\n" );
   fprintf ( stdout, "  -e, --edit-distance         <str>     Choose 'Y' to calculate edit distance and 'N' to output\n"
                     "                                        rotation only. Default: Y.\n" );
   fprintf ( stdout, "  -S, --cost-substitution     <int>     Cost of substitution. Default: 1.\n" );
   fprintf ( stdout, "  -I, --cost-insertion        <int>     Cost of insertion. Default: 1.\n" );
   fprintf ( stdout, "  -D, --cost-deletion         <int>     Cost of deletion. Default: 1.\n" );
 }

double gettime( void )
 {
    struct timeval ttime;
    gettimeofday( &ttime , 0 );
    return ttime.tv_sec + ttime.tv_usec * 0.000001;
 }

void create_rotation ( unsigned char * x, unsigned int offset, unsigned char * rotation )
 {
    unsigned int m = strlen ( ( char * ) x );
    memmove ( &rotation[0], &x[offset], m - offset );
    memmove ( &rotation[m - offset], &x[0], offset );
    rotation[m] = '\0';
 }
