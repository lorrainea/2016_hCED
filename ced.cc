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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <climits>
#include <algorithm>
#include "sacsc.h"
#include "hced.h"
#include "ced.h"
#include "edlib.h"

using namespace std;

#define MIN3(a, b, c) ((a) < (b) ? ((a) < (c) ? (a) : (c)) : ((b) < (c) ? (b) : (c)))
#define MAX3(a, b, c) ((a) > (b) ? ((a) > (c) ? (a) : (c)) : ((b) > (c) ? (b) : (c)))


int delta ( char a, char b, struct TSwitch sw )
 {
    if ( a == DEL || b == DEL ) 
    {
	return 0;
    }
    if ( a == b )
    {
	return sw . m;
    }
    else
    {
	return sw . r;
    }
 }


unsigned int nw_ag_allocation( unsigned int m, unsigned int n, int ** &I, int **& D, int ** &T )
{
	int i;

	if ( ( T = ( int ** ) calloc ( ( m + 1 ) , sizeof( int * ) ) ) == NULL )
        {
                fprintf( stderr, " Error: T could not be allocated!\n");
                return ( 0 );
        }
        for ( i = 0; i < m + 1; i ++ )
        {
                if ( ( T[i] = ( int * ) calloc ( ( n + 1 ) , sizeof( int ) ) ) == NULL )
                {
                        fprintf( stderr, " Error: T could not be allocated!\n");
                        return ( 0 );
                }
        }
	if ( ( I = ( int ** ) calloc ( ( m + 1 ) , sizeof( int * ) ) ) == NULL )
        {
                fprintf( stderr, " Error: I could not be allocated!\n");
                return ( 0 );
        }
        for ( i = 0; i < m + 1; i ++ )
        {
                if ( ( I[i] = ( int * ) calloc ( ( n + 1 ) , sizeof( int ) ) ) == NULL )
                {
                        fprintf( stderr, " Error: I could not be allocated!\n");
                        return ( 0 );
                }
        }
	if ( ( D = ( int ** ) calloc ( ( m + 1 ) , sizeof( int * ) ) ) == NULL )
        {
                fprintf( stderr, " Error: D could not be allocated!\n");
                return ( 0 );
        }
        for ( i = 0; i < m + 1; i ++ )
        {
                if ( ( D[i] = ( int * ) calloc ( ( n + 1 ) , sizeof( int ) ) ) == NULL )
                {
                        fprintf( stderr, " Error: D could not be allocated!\n");
                        return ( 0 );
                }
        }

	return EXIT_SUCCESS;
}


unsigned int nw_ag ( unsigned char * p, unsigned int m, unsigned char * t, unsigned int n, struct TSwitch  sw, int * score, int **&I, int **&D, int **&T)
{
	int i, j;
	int g = sw . O;
	int h = sw . E;
	int u, v, w;
        
        for ( i = 0; i < m + 1; i++ )
	{
		D[i][0] = m * sw . r;
		I[i][0] = m * sw . r;
	}
	for ( j = 0; j < n + 1; j++ )
	{
		D[0][j] = n * sw . r;
		I[0][j] = n * sw . r;
	}

	T[0][0] = 0;
	if ( m > 0 )
		T[1][0] = g;
	for ( i = 2; i < m + 1; i++ )
	T[i][0] = T[i - 1][0] + h;
	if ( n > 0 )
	T[0][1] = g;
	for ( j = 2; j < n + 1; j++ )
	T[0][j] = T[0][j - 1] + h;

	for( i = 1; i < m + 1; i++ )
	{
        	for( j = 1; j < n + 1; j++ )
        	{
			D[i][j] = MAX2 ( D[i - 1][j] + h, T[i - 1][j] + g );
			u = D[i][j];

			I[i][j] = MAX2 ( I[i][j - 1] + h, T[i][j - 1] + g );
			v = I[i][j];

			w = T[i - 1][j - 1] + delta ( t[j - 1], p[i - 1], sw );

			T[i][j] = MAX3 ( w, u, v );
        	}
    	}

	( * score ) = T[m][n];
	
	return EXIT_SUCCESS;
}


unsigned int nw_allocation( unsigned int m, unsigned int n, int ** &T )
{

	if ( ( T = ( int ** ) calloc ( ( m + 1 ) , sizeof( int * ) ) ) == NULL )
        {
                fprintf( stderr, " Error: T could not be allocated!\n");
                return ( 0 );
        }
        for ( int i = 0; i < m + 1; i ++ )
        {
                if ( ( T[i] = ( int * ) calloc ( ( n + 1 ) , sizeof( int ) ) ) == NULL )
                {
                        fprintf( stderr, " Error: T could not be allocated!\n");
                        return ( 0 );
                }
        }
	
	return EXIT_SUCCESS;
}


unsigned int nw ( unsigned char * p, unsigned int m, unsigned char * t, unsigned int n, struct TSwitch  sw, int * score, int **& T )
{
	int ins = sw . f;
	int del = sw . g;
	int i, j;
	
        for ( i = 0; i < m + 1; i++ )
	{
		T[i][0] = del * i;
	}
	for ( j = 1; j < n + 1; j++ )
	{
		T[0][j] = ins * j;
	}

        for ( i = 1; i < m + 1; i++ )
	{
		for ( j = 1; j < n + 1; j++ )
			T[i][j] = MAX3(T[i][j-1] + ins, T[i-1][j] + del, T[i-1][j-1] + delta ( t[j - 1], p[i - 1], sw ));
	}

    	( * score ) = T[m][n];

	return EXIT_SUCCESS;
}

unsigned int sacsc_refinement ( unsigned char * x, unsigned char * y, struct TSwitch  sw, unsigned int * rotation, unsigned int * distance ) 
{
	unsigned int rot;
	unsigned int dist;
	int ** I;
	int ** D;
	int ** T;
	circular_sequence_comparison ( x, y, sw, &rot, &dist );

	( * distance ) = dist;

	unsigned int m = strlen ( ( char * ) x );
	unsigned int n = strlen ( ( char * ) y );
	unsigned char * xr;	
	
	xr = ( unsigned char * ) calloc( ( m + 1 ) , sizeof( unsigned char ) );
	create_rotation ( x, rot, xr );

	unsigned char * X;
	unsigned char * Y;

	unsigned int sl = sw . P * ( sw . l ); //section length
	sl = MIN3 ( sl, m/2, n/2 );

	X = ( unsigned char * ) calloc( ( 3 * sl + 1 ) , sizeof( unsigned char ) );
	Y = ( unsigned char * ) calloc( ( 3 * sl + 1 ) , sizeof( unsigned char ) );

	memcpy ( &X[0], &xr[0], sl );
	for ( int i = 0; i < sl; i++ )
		X[sl + i] = DEL;
	memcpy ( &X[sl + sl], &xr[m - sl], sl );
	X[3 * sl] = '\0';
	
	memcpy ( &Y[0], &y[0], sl );
	for ( int i = 0; i < sl; i++ )
		Y[sl + i] = DEL;
	memcpy ( &Y[sl + sl], &y[n - sl], sl );
	Y[3 * sl] = '\0';

	unsigned int mm = sl + sl + sl;
	unsigned int nn = sl + sl + sl;

	int score = -INT_MAX;
	int max_score = score;
	unsigned int rrot = 0;
	unsigned char * Xr = ( unsigned char * ) calloc( ( 3 * sl + 1 ) , sizeof( unsigned char ) );

	if( sw . O < 0 )
		nw_ag_allocation ( mm, nn, I, D, T );
	else 
		nw_allocation( mm, nn, T );

	for ( int i = 0; i < mm; i++ )
	{
		if ( i >= sl && i < 2 * sl )
			continue;

		memcpy ( &Xr[0], &X[i], ( 3 * sl ) - i );
    		memcpy ( &Xr[( 3 * sl ) - i], &X[0], i );
    		Xr[3 * sl] = '\0';

		if ( sw . O < 0 )
			nw_ag ( Xr, mm , Y, nn, sw, &score, I, D, T );
		else
			nw ( Xr, mm, Y, nn, sw, &score, T );

		if ( score > max_score )
		{
			max_score = score;
			rrot = i;
		}	 
	}

	for ( int j = 0; j < mm + 1; j ++ )
	{
		free ( T[j] );
		if ( sw . O < 0 )
		{
			free ( D[j] );
			free ( I[j] );
		}
	}
	free ( Xr );

	if( sw . O < 0 )
	{
		free ( I );
		free ( D );
	}
	free ( T );

	int final_rot;
        if ( rrot < sl )
        {
                final_rot = rot + rrot;
        }
        else
        {
                final_rot = rot - ( 3 * sl - rrot );
        }

        if ( final_rot > ( int ) m )
        {
                ( * rotation ) = final_rot % m;
        }
        else if ( final_rot < 0 )
        {
                ( * rotation ) = m + final_rot;
        }
        else
                ( * rotation ) = final_rot;
	
	unsigned char * x_final_rotation = ( unsigned char * ) calloc( ( m + 1 ) , sizeof( unsigned char ) );
	
	create_rotation( x, ( *rotation ), x_final_rotation );

	int sub = sw . S;
	int ins = sw . I; 
	int del = sw . D;

	if ( strcmp ( sw . e, edit_distance ) == 0 && ins == 1 && del == 1 && sub == 1 )
	{
		editDistanceMyers( x_final_rotation, y, m, n, distance );
	}
	else if ( strcmp ( sw . e, edit_distance ) == 0 )
	{
		if( m < n )
			editDistance( x_final_rotation, y, m, n, distance, sub, ins, del ); 
		else editDistance( y, x_final_rotation, n, m, distance,  sub, ins, del );	
	}

	free ( xr );
	free ( X );
	free ( Y );
	free( x_final_rotation );

	return EXIT_SUCCESS;
}

/*
Myers Bit-Vector algorithm implemented using edlib Library
*/
int editDistanceMyers( unsigned char * xInput, unsigned char * yInput, int mInput, int nInput, unsigned int * distance )
{
	int score = edlibAlign( (const char*) xInput, strlen( (char*) xInput ), (const char*) yInput, strlen( (char*) yInput ), edlibDefaultAlignConfig()).editDistance;

	( * distance ) = score;

	return EXIT_SUCCESS;
}




unsigned int editDistance(unsigned char * x, unsigned char * y, int xSize, int ySize, unsigned int * distance, int sub, int ins, int del )
{
	
	int pds = 0;
	int mn = min( xSize, ySize );
	int mx = max( xSize, ySize );
	
	unsigned int * ed =  ( unsigned int * ) calloc ( mn + 1 , sizeof(unsigned int));
    	ed[0] = 0;

	for (int i = 1; i < mn + 1; ++i) 
	{
		ed[i] = ed[i - 1] + del;
	}

	int prev_diag = 0;
	for (int j = 1; j < mx + 1; ++j) 
	{
		prev_diag = ed[0], 
		pds = 0;
        	ed[0] = ed[0] + ins;
	

		for (int i = 1; i < mn + 1; ++i) 
		{
		    	pds = ed[i];
			if (x[j - 1] == y[i - 1])	
			{
		        	ed[i] = prev_diag;
			} 
			else 
			{
				ed[i] = MIN3( ed[i - 1] + del, ed[i] + ins ,  prev_diag + sub );
			}
			prev_diag = pds;
		}
	}
  

	( * distance ) = ed[mn];

	free( ed );
	return EXIT_SUCCESS;
}

