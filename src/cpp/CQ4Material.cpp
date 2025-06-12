/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "CQ4Material.h"
#include "Outputter.h" // For COutputter, assuming it's in the same include path context
#include <iomanip>    // For setw
#include <iostream>   // For std::cerr
#include <string>     // For string and getline

//! Read material data from stream Input
//! For Q4 Material, expects: NSet E nu
bool CQ4Material::Read(std::ifstream& Input, unsigned int mset)
{
    nset = mset;
    Input >> E >> nu;
    
    // After reading the required values, we clear the failbit
    // in case the line has extra data and to ensure future reads work.
    if (Input.fail()) {
        Input.clear();
    }
    
    // Consume the rest of the line to avoid issues with subsequent file parsing
    std::string dummy;
    std::getline(Input, dummy);
    
    return true;
}

//! Write material data to Stream
void CQ4Material::Write(COutputter& output)
{
    output << std::setw(5) << nset 
           << std::setw(12) << E 
           << std::setw(12) << nu 
           << std::setw(12) << thickness 
           << std::endl;
} 