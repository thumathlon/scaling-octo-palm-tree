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

#include "Element.h"
#include "CQ4Material.h" // Include the Q4 material definition
#include "Node.h"      // For CNode
#include <fstream>   // For ifstream

// Forward declaration
class COutputter;

//! Q4 Element class (4-node isoparametric quadrilateral)
class CQ4Element : public CElement
{
public:
    //! Constructor
    CQ4Element();

    //! Destructor
    virtual ~CQ4Element();

    //! Read element data from stream Input
    virtual bool Read(std::ifstream& Input, CMaterial** MaterialSets, CNode* NodeList) override;

    //! Write element data to stream
    virtual void Write(COutputter& output) override;

    //! Generate location matrix for Q4 element (2 DOFs per node: UX, UY)
    virtual void GenerateLocationMatrix() override;

    //! Calculate element stiffness matrix (Upper triangular matrix, stored as an array column by column)
    virtual void ElementStiffness(double* stiffness) override;

    //! Calculate element stress/strain
    //! stress_strain array: [sigma_x, sigma_y, tau_xy, epsilon_x, epsilon_y, gamma_xy] for each Gauss point (ordered by Gauss point index)
    //! For a 2x2 Gauss integration, this would be 6 values * 4 points = 24 doubles.
    //! The output format needs to be coordinated with COutputter.
    virtual void ElementStress(double* stress_strain, double* Displacement) override;

    //! Return the number of double values this element will write into stress_array
    //! For Q4: 4 Gauss points * (sig_x, sig_y, tau_xy, eps_x, eps_y, gam_xy) = 4 * 6 = 24 values
    //! + 1 for average von Mises stress
    virtual unsigned int GetNumberOfStressStrainOutputs() const override 
    { 
        return 25; // 4*6 + 1
    }

private:
    // Helper methods for stiffness matrix calculation
    //! Calculate shape functions N at point (xi, eta)
    void ShapeFunctions(double xi, double eta, double N[4]);

    //! Calculate derivatives of shape functions dN/dxi and dN/deta at point (xi, eta)
    void ShapeFunctionDerivatives(double xi, double eta, double dNdxi[4], double dNdeta[4]);

    //! Calculate Jacobian matrix, its determinant, and its inverse at point (xi, eta)
    void Jacobian(double xi, double eta, const CNode* nodes[4], double J[2][2], double& detJ, double invJ[2][2]);

    //! Calculate B matrix (strain-displacement matrix) at point (xi, eta)
    void BMatrix(double xi, double eta, const CNode* nodes[4], double B[3][8]); // 3 strains, 8 DOFs

    //! Calculate D matrix (material constitutive matrix for plane stress)
    void DMatrixPlaneStress(double D[3][3], double E, double nu);

    //! Calculate D matrix (material constitutive matrix for plane strain)
    void DMatrixPlaneStrain(double D[3][3], double E, double nu);
}; 