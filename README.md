Program: fflo

Copyright (C) by Sherryl Manalo

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

Contact the author at: sm@antitachyon.com

The program fflo calculates the FFLO critical field for both s- and d-wave
superconductors as a function of temperature and for various film
thicknesses.


----------------------- compile ----------------------
To compile the program: 
make fflo

----------------------- run --------------------------
To run the program:
fflo fflo.cfg temperaturefile outputname

fflo.cfg: 		config file
temperaturefile:	one-column file with temperature points,
			preferably from high temperatures to lower
			i.e. dtemp.conf
outputname:		main name of files in which the data are stored 
			outputname.lst and outputname.orig are created
			(if calcgap=1: outputname.gap)

---------------------- fflo.cfg ----------------------
sm:		parameter sm

mucon:		constant multiple of mu, normal value: 1.0

calcmu:		1 ... calculate mu-equivalent to sm 
		0 ... do not calclate mu-equivalent to sm 

mudamp:		damping of mu (only active if calcmu = 1)

tolmu:		tolerance of mu (only active if calcmu = 1)

hkapt_0:	first value of H_fflo

ntemp:		number of temperatures to calculate

maxiter:	maximum number of temperatures to calculate

TOL:		tolerance of iteration of temperatures

TOLq:		tolerance of iteration of q

sdamp:		damping of iteration of temperatures

mesh:		mesh of integration

calcgap:	0 ... normal iteration of temperatures
		1 ... calculation of H as function of q at 
		      temperature specified in temperature-file

stepsize:	stepsize of q (only active if calcgap = 1)

rangehmax:	maximum number of q-points to be calculated 
		(only active if calcgap = 1)

qstart:		starting value of q

qdamp:		damping of q

MODUS:		0 ... s-wave
		1 ... d-wave with \phi_q = 0
		2 ... d-wave with variable \phi_q

angle_q:	\phi_q in degrees

angle_plane:	\phi_plane in degrees

qstep:		adds this constant to q at the beginning of the iteration
