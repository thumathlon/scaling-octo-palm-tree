/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#pragma once

#include "Outputter.h"

using namespace std;

//! Material base class
class CMaterial
{
public:
	unsigned int nset;	//!< Material property set number
	double E;			//!< Young's modulus

public:

//!	Constructor
	CMaterial() : nset(0), E(0.0) {}

//! Virtual deconstructor
	virtual ~CMaterial() {}

//!	Read material data from stream Input
	virtual bool Read(std::ifstream& Input, unsigned int nset) = 0;

//!	Write material data to Stream
	virtual void Write(COutputter& output) = 0;

};

//! Material class for bar element
class CBarMaterial : public CMaterial
{
public:

	double Area;	//!< Sectional area of a bar element

public:
	
//!	Constructor
	CBarMaterial() : CMaterial(), Area(0.0) {}

//!	Read material data from stream Input
	virtual bool Read(std::ifstream& Input, unsigned int mset) override;

//!	Write material data to Stream
	virtual void Write(COutputter& output) override;
};
