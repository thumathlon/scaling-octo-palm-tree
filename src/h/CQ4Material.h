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

#include "Material.h" // 包含基类 CMaterial
#include <fstream>   // For ifstream

// Forward declaration
class COutputter;

enum class AnalysisType2D { // 使用 enum class 增强类型安全
    PlaneStress,
    PlaneStrain
};

//! Material class for Q4 (4-node quadrilateral) plane stress/strain element
class CQ4Material : public CMaterial
{
public:
    double nu;        //!< Poisson's ratio
    double thickness; //!< Element thickness
    AnalysisType2D analysisType; 

public:
    //! Constructor
    CQ4Material() : CMaterial(), nu(0.0), thickness(1.0), analysisType(AnalysisType2D::PlaneStress) {}

    //! Read material data from stream Input
    virtual bool Read(std::ifstream& Input, unsigned int mset) override;

    //! Write material data to Stream
    virtual void Write(COutputter& output) override;
}; 