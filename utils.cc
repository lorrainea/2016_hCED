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
#include "csc.h"


static struct option long_options[] =
 {
   { "alphabet",                required_argument, NULL, 'a' },
   { "input-file",              required_argument, NULL, 'i' },
   { "output-file",             required_argument, NULL, 'o' },
   { "q-length-min",            required_argument, NULL, 'q' },
   { "q-length-max",            optional_argument, NULL, 'Q' },
   { "block-length-min",        required_argument, NULL, 'l' },
   { "block-length-max",        optional_argument, NULL, 'L' },
   { "percent-refine",          optional_argument, NULL, 'P' },
   { "edit-distance",           optional_argument, NULL, 'e' },
   { "cost-match",     		optional_argument, NULL, 'M' },
   { "cost-substitution",	optional_argument, NULL, 'S' },
   { "cost-insertion",		optional_argument, NULL, 'I' },
   { "cost-deletion",		optional_argument, NULL, 'D' },
   { "help",                    no_argument,       NULL, 'h' },
   {  NULL,                     0,                 NULL,  0  }
 };


/* 
Decode the input switches 
*/
int decode_switches ( int argc, char * argv [], struct TSwitch * sw )
 {
   int          oi;
   int          opt;
   double       val;
   char       * ep;
   int          args;

   /* initialisation */
   sw -> alphabet                       = NULL;
   sw -> input_filename                 = NULL;
   sw -> output_filename                = NULL;
   sw -> q                              = 5;
   sw -> Q                              = 0;
   sw -> l                              = 10;
   sw -> L                              = 0;
   sw -> P                              = 0.0;
   sw -> e      			= NULL;
   sw -> M				= 0;
   sw -> S				= 1;
   sw -> I 				= 1;
   sw -> D				= 1;
   args = 0;

   while ( ( opt = getopt_long ( argc, argv, "a:i:o:q:Q:l:L:P:e:M:S:I:D:h", long_options, &oi ) ) != - 1 )
    {
      switch ( opt )
       {
         case 'a':
           sw -> alphabet = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> alphabet, optarg );
           args ++;
           break;

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

	 case 'e':
           sw -> e = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> e, optarg );
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

         case 'Q':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> Q = val;
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

         case 'L':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> L = val;
           break;

         case 'P':
           val = (double) atof ( optarg );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> P = val;
           break;
	
	case 'M':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> M = val;
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

         case 'h':
           return ( 0 );
       }
    }

   if ( sw -> Q == 0 )
     {
       sw -> Q = sw -> q;
     }
   if ( sw -> L == 0 )
     {
       sw -> L = sw -> l;
     }

   if ( args < 6 )
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
   fprintf ( stdout, " Required:\n" );
   fprintf ( stdout, "  -a, --alphabet              <str>     `DNA' for nucleotide  sequences or `PROT'\n"
                     "                                        for protein  sequences. \n" );
   fprintf ( stdout, "  -i, --input-file            <str>     (Multi)FASTA input filename.\n" );
   fprintf ( stdout, "  -o, --output-file           <str>     Output filename for the rotated sequences.\n" );
   fprintf ( stdout, "  -q, --q-length-min          <int>     The q-gram length.\n");
   fprintf ( stdout, "  -l, --block-length-min      <int>     The length of each block.\n\n");
   fprintf ( stdout, "  -Q, --q-length-max          <int>     The maximum q-gram length. The program\n"
                     "                                        will try all in range min .. max.\n" );
   fprintf ( stdout, "  -L, --block-length-max      <int>     The maximum length of each block. The\n"
                     "                                        program will try all in range min .. max.\n" );
   fprintf ( stdout, "  -P, --percent-refine        <float>   Refine the alignment of hCSC/saCSC by\n"
                     "                                        checking a percentage of the ends (e.g. 2.5)\n" );
   fprintf ( stdout, "  -e, --edit distance method  <str>     Choose 'Y' for Myers or 'V' for standard edit distance\n\n" );
   fprintf ( stdout, " Optional:\n\n" );
   fprintf ( stdout, "  -M, --Cost of match         <int>     Cost of match when 'V' is chosen for edit\n"
                     "                                        distance method. Default: M = 0.\n" );
   fprintf ( stdout, "  -S, --Cost of substitution  <int>     Cost of substitution when 'V is chosen for edit\n"
                     "                                        distance method. Default: S = 1.\n" );
   fprintf ( stdout, "  -I, --Cost of insertion     <int>     Cost of insertion when 'V' is chosen for edit\n"
                     "                                        distance method. Default: I = 1.\n" );
   fprintf ( stdout, "  -D, --Cost of deletion      <int>     Cost of deletion when 'V' is chosen for edit\n"
                     "                                        distance method. Default: D = 1.\n" );
	
	
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
