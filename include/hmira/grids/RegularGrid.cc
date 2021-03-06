/*===========================================================================*\
 *                                                                           *
 *                                IsoEx                                      *
 *        Copyright (C) 2002 by Computer Graphics Group, RWTH Aachen         *
 *                         www.rwth-graphics.de                              *
 *                                                                           *
 *---------------------------------------------------------------------------*
 *                                                                           *
 *                                License                                    *
 *                                                                           *
 *  This library is free software; you can redistribute it and/or modify it  *
 *  under the terms of the GNU Library General Public License as published   *
 *  by the Free Software Foundation, version 2.                              *
 *                                                                           *
 *  This library is distributed in the hope that it will be useful, but      *
 *  WITHOUT ANY WARRANTY; without even the implied warranty of               *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU        *
 *  Library General Public License for more details.                         *
 *                                                                           *
 *  You should have received a copy of the GNU Library General Public        *
 *  License along with this library; if not, write to the Free Software      *
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                *
 *                                                                           *
\*===========================================================================*/

//=============================================================================
//
//  CLASS RegularGrid - IMPLEMENTATION
//
//=============================================================================

//== INCLUDES =================================================================

#include "RegularGrid.hh"
#include <algorithm>

//== NAMESPACES ===============================================================

namespace IsoEx
{

//== IMPLEMENTATION ==========================================================


void
RegularGrid::
initialize( const Vec3&  _origin,
            const Vec3&  _x_axis,
            const Vec3&  _y_axis,
            const Vec3&  _z_axis,
            unsigned int            _x_res,
            unsigned int            _y_res,
            unsigned int            _z_res )
{
    origin_ = _origin;
    x_axis_ = _x_axis;
    y_axis_ = _y_axis;
    z_axis_ = _z_axis;
    x_res_ = _x_res;
    y_res_ = _y_res;
    z_res_ = _z_res;

    // CubeIdx and PointIdx have 32 bits
    // -> its 3 components should require <= 10 bits
    assert( x_res_ < 1024 );
    assert( y_res_ < 1024 );
    assert( z_res_ < 1024 );

    n_cubes_  = ( _x_res - 1 ) * ( _y_res - 1 ) * ( _z_res - 1 );
    n_points_ = _x_res * _y_res * _z_res;

    dx_ = x_axis_ / ( real )( x_res_ -1 );
    dy_ = y_axis_ / ( real )( y_res_ -1 );
    dz_ = z_axis_ / ( real )( z_res_ -1 );

    Vec3 x_axisn = x_axis_; x_axisn.normalize();
    Vec3 y_axisn = y_axis_; y_axisn.normalize();
    Vec3 z_axisn = z_axis_; z_axisn.normalize();


    offsets_[0] = 0;
    offsets_[1] = 1;
    offsets_[2] = 1 + x_res_;
    offsets_[3] =     x_res_;
    offsets_[4] =              x_res_ * y_res_;
    offsets_[5] = 1          + x_res_ * y_res_;
    offsets_[6] = 1 + x_res_ + x_res_ * y_res_;
    offsets_[7] =     x_res_ + x_res_ * y_res_;
}


//-----------------------------------------------------------------------------

RegularGrid::PointIdx
RegularGrid::
nearest_point( const Vec3 &_p )
{
    Vec3 pl = ( _p );
    
    // grid spacing
    real dxf = dx().norm();
    real dyf = dy().norm();
    real dzf = dz().norm();

    // get gird cell
    int x = ( int )( pl[0] / dxf );
    int y = ( int )( pl[1] / dyf );
    int z = ( int )( pl[2] / dzf );

    // transform to point coordinates
    return x + y * x_res_ + z * x_res_ * y_res_;
}

RegularGrid::CubeIdx
RegularGrid::
get_proper_cube( const Vec3 &_p)
{
	unsigned int X( x_res_ -1 ), Y( y_res_ -1 );

	Vec3 pl = const_cast<Vec3&>( _p );
	
	real dxf = dx().norm();
	real dyf = dy().norm();
	real dzf = dz().norm();
	
	int x = ( int )( pl[0] / dxf );
	int y = ( int )( pl[1] / dyf );
	int z = ( int )( pl[2] / dzf );

	return x + y * X + z * X * Y;
}

//-----------------------------------------------------------------------------


RegularGrid::PointIdx
RegularGrid::
point_idx( CubeIdx _idx, unsigned char _corner ) const
{
    assert( _corner < 8 );

    // get cube coordinates
    unsigned int X( x_res_ -1 ), Y( y_res_ -1 );
    unsigned int x = _idx % X;  _idx /= X;
    unsigned int y = _idx % Y;  _idx /= Y;
    unsigned int z = _idx;

    // transform to point coordinates
    _idx = x + y * x_res_ + z * x_res_ * y_res_;

    // add offset
    return _idx + offsets_[_corner];
}


//-----------------------------------------------------------------------------


RegularGrid::Vec3
RegularGrid::
point( PointIdx _idx ) const
{
    unsigned int x = _idx % x_res_;  _idx /= x_res_;
    unsigned int y = _idx % y_res_;  _idx /= y_res_;
    unsigned int z = _idx;

    return origin_ + dx_*x + dy_*y + dz_*z;
}

//-----------------------------------------------------------------------------

RegularGrid::Vec3
RegularGrid::
point( int x, int y, int z ) const
{
    return origin_ + dx_*x + dy_*y + dz_*z;
}

//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------

bool
RegularGrid::
ray_intersect_local( const Vec3 &_o, const Vec3 &_d,
                     Vec3  &_entry,  Vec3 &_exit )
{
    // ray cube intersections
    std::vector< Intersection > intersections_;

    // upper,right,front cube corner
    Vec3 cc( x_axis().norm(), y_axis().norm(), z_axis().norm() );

    // find ray-cube intersections
    for ( unsigned int i = 0; i < 3; ++i )
    {
        if ( fabs( _d[i] ) > 1E-9 )
        {
            // enter
            intersections_.resize( intersections_.size() + 1 );
            intersections_.back().lambda = -_o[i]/_d[i];
            if ( _d[i] > 0.0 )
                intersections_.back().enter  = true;
            else
                intersections_.back().enter  = false;

            // leaving
            intersections_.resize( intersections_.size() + 1 );
            intersections_.back().lambda = ( _o[i]-cc[i] )/( -_d[i] );
            if ( _d[i] > 0.0 )
                intersections_.back().enter  = false;
            else
                intersections_.back().enter  = true;
        }
    }

    // find last entry and first exit point
    real l_en = -1000000.f;
    real l_ex = 1000000.f;
    for ( unsigned int i = 0; i < intersections_.size(); ++i )
    {
        if ( intersections_[i].enter )
        {
            if ( l_en < intersections_[i].lambda )
                l_en = intersections_[i].lambda;
        }
        else
        {
            if ( l_ex > intersections_[i].lambda )
                l_ex = intersections_[i].lambda;
        }
    }

    // does the ray intersect the cube
    if ( l_en == -1000000.f || l_ex == 1000000.f || l_ex < l_en )
        return false;

    _entry = _o + _d*l_en;
    _exit  = _o + _d*l_ex;

    return true;
}


//=============================================================================
} // namespace IsoEx
//=============================================================================
