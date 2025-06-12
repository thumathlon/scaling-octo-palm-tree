/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "ElementGroup.h"
#include "Domain.h"
#include "../h/CQ4Element.h"
#include "../h/CQ4Material.h"
#include "../h/Bar.h"
#include "../h/Material.h"
#include <iostream>

CNode* CElementGroup::NodeList_ = nullptr;

//! Constructor
CElementGroup::CElementGroup()
{
    if (!NodeList_)
    {
        CDomain* FEMData = CDomain::GetInstance();
        NodeList_ = FEMData->GetNodeList();
    }
    ElementType_ = ElementTypes::UNDEFINED;
    NUME_ = 0;
    ElementList_ = nullptr;
    NUMMAT_ = 0;
    MaterialList_ = nullptr;
}

//! Deconstructor
CElementGroup::~CElementGroup()
{
    if (ElementList_) {
        for (unsigned int i = 0; i < NUME_; ++i) {
            delete ElementList_[i];
        }
        delete [] ElementList_;
    }
    if (MaterialList_) {
        for (unsigned int i = 0; i < NUMMAT_; ++i) {
            delete MaterialList_[i];
        }
        delete [] MaterialList_;
    }
}

//! operator []
CElement& CElementGroup::operator[](unsigned int i)
{
    return *ElementList_[i];
}

//! Return index-th material in this element group
CMaterial& CElementGroup::GetMaterial(unsigned int i)
{
    return *MaterialList_[i];
}

//! Read element group data from stream Input
bool CElementGroup::Read(ifstream& Input)
{
    Input >> (int&)ElementType_ >> NUME_ >> NUMMAT_;
    
    if (Input.fail()) {
        std::cerr << "*** Error *** Failed to read element group header." << std::endl;
        return false;
    }

    // Allocate pointer arrays
    if (NUMMAT_ > 0) {
        MaterialList_ = new CMaterial*[NUMMAT_];
    }
    if (NUME_ > 0) {
        ElementList_ = new CElement*[NUME_];
    }

    // Read material properties
    for (unsigned int m = 0; m < NUMMAT_; ++m) {
        // Read material set number first
        unsigned int mset;
        Input >> mset;
        if (Input.fail() || mset != m + 1) {
             std::cerr << "*** Error *** Reading material set " << m + 1 << std::endl;
             return false;
        }

        switch(ElementType_) {
            case ElementTypes::Bar:
                MaterialList_[m] = new CBarMaterial();
                break;
            case ElementTypes::Q4:
                MaterialList_[m] = new CQ4Material();
                break;
            default:
                std::cerr << "Material allocation failed: Unsupported element type " << ElementType_ << std::endl;
                return false;
        }
        if (!MaterialList_[m]->Read(Input, mset)) {
            std::cerr << "*** Error *** Reading material data for set " << mset << std::endl;
            return false;
        }
    }

    // Read element data
    for (unsigned int e = 0; e < NUME_; ++e) {
        unsigned int N;
        Input >> N;
        if (Input.fail() || N != e + 1) {
            std::cerr << "*** Error *** Reading element number " << e + 1 << std::endl;
            return false;
        }
        
        switch(ElementType_) {
            case ElementTypes::Bar:
                ElementList_[e] = new CBar();
                break;
            case ElementTypes::Q4:
                ElementList_[e] = new CQ4Element();
                break;
            default:
                std::cerr << "Element allocation failed: Unsupported element type " << ElementType_ << std::endl;
                return false;
        }
        if (!ElementList_[e]->Read(Input, MaterialList_, NodeList_)) {
            std::cerr << "*** Error *** Reading element data for element " << N << std::endl;
            return false;
        }
    }

    return true;
}
