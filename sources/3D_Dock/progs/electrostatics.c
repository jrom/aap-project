/*
This file is part of ftdock, a program for rigid-body protein-protein docking 
Copyright (C) 1997-2000 Gidon Moont

Biomolecular Modelling Laboratory
Imperial Cancer Research Fund
44 Lincoln's Inn Fields
London WC2A 3PX

+44 (0)20 7269 3348
http://www.bmm.icnet.uk/

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

*/

#include "structures.h"

void assign_charges( struct Structure This_Structure ) {

/************/

  /* Counters */

  int	residue , atom ;

/************/

  for( residue = 1 ; residue <= This_Structure.length ; residue ++ ) {
    for( atom = 1 ; atom <= This_Structure.Residue[residue].size ; atom ++ ) {

      This_Structure.Residue[residue].Atom[atom].charge = 0.0 ;

      /* peptide backbone */

      if( strcmp( This_Structure.Residue[residue].Atom[atom].atom_name , " N  " ) == 0 ) {
        if( strcmp( This_Structure.Residue[residue].res_name , "PRO" ) == 0 ) {
          This_Structure.Residue[residue].Atom[atom].charge = -0.10 ;
        } else {
          This_Structure.Residue[residue].Atom[atom].charge =  0.55 ;
          if( residue == 1 ) This_Structure.Residue[residue].Atom[atom].charge = 1.00 ;
        }
      }

      if( strcmp( This_Structure.Residue[residue].Atom[atom].atom_name , " O  " ) == 0 ) {
        This_Structure.Residue[residue].Atom[atom].charge = -0.55 ;
        if( residue == This_Structure.length  ) This_Structure.Residue[residue].Atom[atom].charge = -1.00 ;
      }

      /* charged residues */

      if( ( strcmp( This_Structure.Residue[residue].res_name , "ARG" ) == 0 ) && ( strncmp( This_Structure.Residue[residue].Atom[atom].atom_name , " NH" , 3 ) == 0 ) ) This_Structure.Residue[residue].Atom[atom].charge =  0.50 ;
      if( ( strcmp( This_Structure.Residue[residue].res_name , "ASP" ) == 0 ) && ( strncmp( This_Structure.Residue[residue].Atom[atom].atom_name , " OD" , 3 ) == 0 ) ) This_Structure.Residue[residue].Atom[atom].charge = -0.50 ;
      if( ( strcmp( This_Structure.Residue[residue].res_name , "GLU" ) == 0 ) && ( strncmp( This_Structure.Residue[residue].Atom[atom].atom_name , " OE" , 3 ) == 0 ) ) This_Structure.Residue[residue].Atom[atom].charge = -0.50 ;
      if( ( strcmp( This_Structure.Residue[residue].res_name , "LYS" ) == 0 ) && ( strcmp( This_Structure.Residue[residue].Atom[atom].atom_name , " NZ " ) == 0 ) ) This_Structure.Residue[residue].Atom[atom].charge =  1.00 ;

    }
  }

/************/

}



/************************/



void electric_field( struct Structure This_Structure , float grid_span , int grid_size , fftw_real *grid )
{

/************/

  /* Counters */

  int	residue , atom ;

  /* Co-ordinates */

  int	x , y , z ;
  float		x_centre , y_centre , z_centre ;

  /* Variables */

  float		distance ;
  float		phi , epsilon ;

/************/

  for( x = 0 ; x < grid_size ; x ++ ) {
    for( y = 0 ; y < grid_size ; y ++ ) {
      for( z = 0 ; z < grid_size ; z ++ ) {

        grid[gaddress(x,y,z,grid_size)] = (fftw_real)0 ;

      }
    }
  }

/************/



  // Pre-process needed floats

  unsigned int cnt = 0;
  for( residue = 1 ; residue <= This_Structure.length ; residue ++ )
    cnt += residue * This_Structure.Residue[residue].size;

  unsigned int atoms = 0;
  float *atom_coords, *aux_coord;
  float *atom_charges, *aux_charge;
  float *atom_distances, *aux_distance;
  posix_memalign((void**) &atom_coords, 16, cnt * 4 * sizeof(float));
  posix_memalign((void**) &atom_charges, 16, cnt * sizeof(float));

  for( aux_coord = atom_coords, aux_charge = atom_charges, residue = 1 ; residue <= This_Structure.length ; residue ++ )
  {
    for( atom = 1 ; atom <= This_Structure.Residue[residue].size ; atom ++ )
    {
      if( This_Structure.Residue[residue].Atom[atom].charge != 0 )
      {
        *aux_coord = This_Structure.Residue[residue].Atom[atom].coord[1];
        *(aux_coord + 1) = This_Structure.Residue[residue].Atom[atom].coord[2];
        *(aux_coord + 2) = This_Structure.Residue[residue].Atom[atom].coord[3];
        aux_coord += 4;
        *aux_charge = This_Structure.Residue[residue].Atom[atom].charge;
        aux_charge++;
        atoms++;
      }
    }
  }

  atom_distances = malloc(atoms * sizeof(float));
  // End pre-process

  setvbuf( stdout , (char *)NULL , _IONBF , 0 ) ;

  printf( "  electric field calculations ( one dot / grid sheet ) " ) ;
  for( x = 0 ; x < grid_size ; x ++ )
  {
    printf( "." ) ;
    x_centre  = gcentre( x , grid_span , grid_size ) ;
    for( y = 0 ; y < grid_size ; y ++ )
    {
      y_centre  = gcentre( y , grid_span , grid_size ) ;
      for( z = 0 ; z < grid_size ; z ++ )
      {
        z_centre  = gcentre( z , grid_span , grid_size ) ;
        phi = 0 ;
        aux_coord = atom_coords;
        aux_distance = atom_distances;
        __m128 _centers = _mm_setr_ps(x_centre, y_centre, z_centre, 0.0);
        for ( atom = 1; atom < atoms-3; atom += 4)
        {
          float pyth1, pyth2, pyth3, pyth4;
          __m128 pyths;
          __m128 *coords = (__m128*) aux_coord;
          __m128 aux = _mm_sub_ps(*coords, _centers);
          aux = _mm_mul_ps(aux, aux);
          *((float*) &pyths) = *((float*) &aux) + *((float*) (&aux)+1) + *((float*) (&aux)+2);
          aux_coord += 4;

          coords = (__m128*) aux_coord;
          aux = _mm_sub_ps(*coords, _centers);
          aux = _mm_mul_ps(aux, aux);
          *((float*) &pyths+1) = *((float*) &aux) + *((float*) (&aux)+1) + *((float*) (&aux)+2);
          aux_coord += 4;

          coords = (__m128*) aux_coord;
          aux = _mm_sub_ps(*coords, _centers);
          aux = _mm_mul_ps(aux, aux);
          *((float*) &pyths+2) = *((float*) &aux) + *((float*) (&aux)+1) + *((float*) (&aux)+2);
          aux_coord += 4;

          coords = (__m128*) aux_coord;
          aux = _mm_sub_ps(*coords, _centers);
          aux = _mm_mul_ps(aux, aux);
          *((float*) &pyths+3) = *((float*) &aux) + *((float*) (&aux)+1) + *((float*) (&aux)+2);
          aux_coord += 4;

          *((__m128*) aux_distance) = _mm_rsqrt_ps(pyths);

          aux_distance +=4 ;
        }

        for ( ; atom < atoms; atom++)
        {
          __m128 *coords = (__m128*) aux_coord;
          __m128 aux = _mm_sub_ps(*coords, _centers);
          aux = _mm_mul_ps(aux, aux);
          *aux_distance = sqrt( *((float*) &aux) + *((float*) (&aux)+1) + *((float*) (&aux)+2) );
          aux_coord += 4;
          aux_distance++;
        }

        aux_charge = atom_charges;
        aux_distance = atom_distances;
        for ( atom = 1; atom < atoms; atom++)
        {
          if( *aux_distance < 2.0 ) *aux_distance = 2.0 ;
          if( *aux_distance >= 8.0 )
          {
            epsilon = 80 ;
          }
          else
          { 
            if( *aux_distance <= 6.0 )
            { 
              epsilon = 4 ;
            }
            else
            {
              epsilon = ( 38 * *aux_distance ) - 224 ;
            }
          }
          phi += ( *aux_charge / ( epsilon * *aux_distance ) ) ;
          aux_charge++;
          aux_distance++;
        }

        grid[gaddress(x,y,z,grid_size)] = (fftw_real)phi;
      }
    }
  }

  printf( "\n" ) ;

/************/

  return ;

}



/************************/



void electric_point_charge( struct Structure This_Structure , float grid_span , int grid_size , fftw_real *grid ) {

/************/

  /* Counters */

  int	residue , atom ;

  /* Co-ordinates */

  int	x , y , z ;
  int	x_low , x_high , y_low , y_high , z_low , z_high ;

  float		a , b , c ;
  float		x_corner , y_corner , z_corner ;
  float		w ;

  /* Variables */

  float		one_span ;

/************/

  for( x = 0 ; x < grid_size ; x ++ ) {
    for( y = 0 ; y < grid_size ; y ++ ) {
      for( z = 0 ; z < grid_size ; z ++ ) {

        grid[gaddress(x,y,z,grid_size)] = (fftw_real)0 ;

      }
    }
  }

/************/

  one_span = grid_span / (float)grid_size ;

  for( residue = 1 ; residue <= This_Structure.length ; residue ++ ) {
    for( atom = 1 ; atom <= This_Structure.Residue[residue].size ; atom ++ ) {

      if( This_Structure.Residue[residue].Atom[atom].charge != 0 ) {

        x_low = gord( This_Structure.Residue[residue].Atom[atom].coord[1] - ( one_span / 2 ) , grid_span , grid_size ) ;
        y_low = gord( This_Structure.Residue[residue].Atom[atom].coord[2] - ( one_span / 2 ) , grid_span , grid_size ) ;
        z_low = gord( This_Structure.Residue[residue].Atom[atom].coord[3] - ( one_span / 2 ) , grid_span , grid_size ) ;

        x_high = x_low + 1 ;
        y_high = y_low + 1 ;
        z_high = z_low + 1 ;

        a = This_Structure.Residue[residue].Atom[atom].coord[1] - gcentre( x_low , grid_span , grid_size ) - ( one_span / 2 ) ;
        b = This_Structure.Residue[residue].Atom[atom].coord[2] - gcentre( y_low , grid_span , grid_size ) - ( one_span / 2 ) ;
        c = This_Structure.Residue[residue].Atom[atom].coord[3] - gcentre( z_low , grid_span , grid_size ) - ( one_span / 2 ) ;

        for( x = x_low ; x <= x_high  ; x ++ ) {
 
          x_corner = one_span * ( (float)( x - x_high ) + .5 ) ;

          for( y = y_low ; y <= y_high  ; y ++ ) {

            y_corner = one_span * ( (float)( y - y_high ) + .5 ) ;

            for( z = z_low ; z <= z_high  ; z ++ ) {

              z_corner = one_span * ( (float)( z - z_high ) + .5 ) ;

              w = ( ( x_corner + a ) * ( y_corner + b ) * ( z_corner + c ) ) / ( 8.0 * x_corner * y_corner * z_corner ) ;

              grid[gaddress(x,y,z,grid_size)] += (fftw_real)( w * This_Structure.Residue[residue].Atom[atom].charge ) ;

            }
          }
        }

      }

    }
  }

/************/

  return ;

}



/************************/



void electric_field_zero_core( int grid_size , fftw_real *elec_grid , fftw_real *surface_grid , float internal_value ) {

/************/

  /* Co-ordinates */

  int	x , y , z ;

/************/

  for( x = 0 ; x < grid_size ; x ++ ) {
    for( y = 0 ; y < grid_size ; y ++ ) {
      for( z = 0 ; z < grid_size ; z ++ ) {

        if( surface_grid[gaddress(x,y,z,grid_size)] == (fftw_real)internal_value ) elec_grid[gaddress(x,y,z,grid_size)] = (fftw_real)0 ;

      }
    }
  }

/************/

  return ;

}
