/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "CQ4Element.h"
#include "Outputter.h" // For COutputter
#include <iostream>    // For cerr
#include <iomanip>     // For setw
#include <cmath>       // For sqrt, fabs

// Function to clear an array (already in CElement.h, but for self-containment if needed)
// template <class type> void clear(type* a, unsigned int N) {
//     for (unsigned int i = 0; i < N; i++)
//         a[i] = 0;
// }

// Constructor
CQ4Element::CQ4Element()
{
    NEN_ = 4; // 4 nodes per element
    nodes_ = new CNode*[NEN_];
    for (unsigned int i = 0; i < NEN_; ++i) nodes_[i] = nullptr;

    // Assuming 2 DOFs per node (UX, UY) for plane problems
    // CNode::NDF should reflect this (e.g. might need to be configured or checked)
    // If CNode::NDF is 3 (for 3D problems like Bar), this needs adjustment or
    // the Q4 element needs to handle 3D coordinates but use only 2D DOFs for plane problems.
    // For now, assume CNode::NDF is appropriate for 2D plane problem (e.g. 2).
    // Let's assume node has X, Y, Z, but for Q4 we only use X, Y.
    // And each node has 2 DOFs for this element: horizontal and vertical displacement in XY plane.
    ND_ = NEN_ * 2; // 4 nodes * 2 DOFs/node = 8 DOFs
    LocationMatrix_ = new unsigned int[ND_];
    for (unsigned int i = 0; i < ND_; ++i) LocationMatrix_[i] = 0;

    ElementMaterial_ = nullptr;
}

// Destructor
CQ4Element::~CQ4Element()
{
    // Base class destructor handles deletion of nodes_, ElementMaterial_, LocationMatrix_
    // if they were allocated by base. Here nodes_ and LocationMatrix_ are allocated in derived.
    // CElement::~CElement() handles nodes_ (if not null), ElementMaterial_ (if not null), LocationMatrix_ (if not null)
    // This might lead to double deletion if CElement also allocates. Let's check CElement constructor/destructor.
    // CElement constructor: NEN_(0), nodes_(nullptr), ElementMaterial_(nullptr)
    // CElement destructor: delete [] nodes_; delete [] ElementMaterial_; delete [] LocationMatrix_;
    // So, if derived class allocates these, base class should not, or derived should nullify them before base destructor is called.
    // The current CElement destructor will delete them. This is fine as long as they are allocated here.
}

// Read element data from stream Input
bool CQ4Element::Read(std::ifstream& Input, CMaterial** MaterialSets, CNode* NodeList)
{
    unsigned int MSet;
    unsigned int N[4];
    double element_thickness;
    int type_code;

    Input >> N[0] >> N[1] >> N[2] >> N[3] >> MSet >> element_thickness >> type_code;
    if (Input.fail()) {
        return false;
    }

    for (unsigned int i = 0; i < NEN_; ++i)
    {
        nodes_[i] = &NodeList[N[i] - 1];
    }

    ElementMaterial_ = MaterialSets[MSet - 1];
    if (!ElementMaterial_) {
        return false;
    }

    CQ4Material* mat = dynamic_cast<CQ4Material*>(ElementMaterial_);
    if (mat) {
        mat->thickness = element_thickness;
        if (type_code == 1) {
            mat->analysisType = AnalysisType2D::PlaneStrain;
        } else {
            mat->analysisType = AnalysisType2D::PlaneStress;
        }
    } else {
        return false;
    }

    return true;
}

// Write element data to stream
void CQ4Element::Write(COutputter& output)
{
    output << std::setw(5) << "Q4" // Element type identifier
           << std::setw(9) << nodes_[0]->NodeNumber
           << std::setw(9) << nodes_[1]->NodeNumber
           << std::setw(9) << nodes_[2]->NodeNumber
           << std::setw(9) << nodes_[3]->NodeNumber
           << std::setw(12) << (dynamic_cast<CQ4Material*>(ElementMaterial_))->nset << std::endl;
}

// Shape functions N(xi, eta)
void CQ4Element::ShapeFunctions(double xi, double eta, double N[4])
{
    N[0] = 0.25 * (1.0 - xi) * (1.0 - eta);
    N[1] = 0.25 * (1.0 + xi) * (1.0 - eta);
    N[2] = 0.25 * (1.0 + xi) * (1.0 + eta);
    N[3] = 0.25 * (1.0 - xi) * (1.0 + eta);
}

// Derivatives of shape functions dN/dxi, dN/deta
void CQ4Element::ShapeFunctionDerivatives(double xi, double eta, double dNdxi[4], double dNdeta[4])
{
    dNdxi[0] = -0.25 * (1.0 - eta);
    dNdxi[1] =  0.25 * (1.0 - eta);
    dNdxi[2] =  0.25 * (1.0 + eta);
    dNdxi[3] = -0.25 * (1.0 + eta);

    dNdeta[0] = -0.25 * (1.0 - xi);
    dNdeta[1] = -0.25 * (1.0 + xi);
    dNdeta[2] =  0.25 * (1.0 + xi);
    dNdeta[3] =  0.25 * (1.0 - xi);
}

// Jacobian matrix J, its determinant detJ, and inverse invJ
void CQ4Element::Jacobian(double xi, double eta, const CNode* elem_nodes[4], double J[2][2], double& detJ, double invJ[2][2])
{
    double dNdxi[4], dNdeta[4];
    ShapeFunctionDerivatives(xi, eta, dNdxi, dNdeta);

    J[0][0] = J[0][1] = J[1][0] = J[1][1] = 0.0;

    for (unsigned int i = 0; i < NEN_; ++i)
    {
        J[0][0] += dNdxi[i] * elem_nodes[i]->XYZ[0]; // dX/dxi
        J[0][1] += dNdxi[i] * elem_nodes[i]->XYZ[1]; // dY/dxi
        J[1][0] += dNdeta[i] * elem_nodes[i]->XYZ[0]; // dX/deta
        J[1][1] += dNdeta[i] * elem_nodes[i]->XYZ[1]; // dY/deta
    }

    detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];

    if (std::fabs(detJ) < 1.0e-12) { // Check for zero or very small determinant
        std::cerr << "*** Error: Jacobian determinant is zero or too small for Q4 element. detJ = " << detJ << std::endl;
        // Handle error, possibly by setting detJ to a small non-zero value or throwing an exception
        // For now, let it proceed, which might lead to NaN/inf if detJ is actually zero.
        // A robust implementation would stop or use a very small number to avoid division by zero.
    }

    double inv_detJ = 1.0 / detJ;
    invJ[0][0] =  J[1][1] * inv_detJ;
    invJ[0][1] = -J[0][1] * inv_detJ;
    invJ[1][0] = -J[1][0] * inv_detJ;
    invJ[1][1] =  J[0][0] * inv_detJ;
}

// B matrix (Strain-Displacement matrix)
void CQ4Element::BMatrix(double xi, double eta, const CNode* elem_nodes[4], double B[3][8])
{
    // Clear B matrix before filling it
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 8; ++j) {
            B[i][j] = 0.0;
        }
    }

    double dNdxi[4], dNdeta[4];
    ShapeFunctionDerivatives(xi, eta, dNdxi, dNdeta);

    double J[2][2], invJ[2][2], detJ;
    Jacobian(xi, eta, elem_nodes, J, detJ, invJ);

    double dNdx[4], dNdy[4];
    for (unsigned int i = 0; i < NEN_; ++i)
    {
        dNdx[i] = invJ[0][0] * dNdxi[i] + invJ[0][1] * dNdeta[i];
        dNdy[i] = invJ[1][0] * dNdxi[i] + invJ[1][1] * dNdeta[i];
    }

    for (unsigned int i = 0; i < NEN_; ++i)
    {
        B[0][2 * i] = dNdx[i]; B[0][2 * i + 1] = 0.0;
        B[1][2 * i] = 0.0;     B[1][2 * i + 1] = dNdy[i];
        B[2][2 * i] = dNdy[i]; B[2][2 * i + 1] = dNdx[i];
    }
}

// D matrix (Material Constitutive matrix for Plane Stress)
void CQ4Element::DMatrixPlaneStress(double D[3][3], double E, double nu)
{
    double factor = E / (1.0 - nu * nu);

    D[0][0] = factor * 1.0;
    D[0][1] = factor * nu;
    D[0][2] = 0.0;

    D[1][0] = factor * nu;
    D[1][1] = factor * 1.0;
    D[1][2] = 0.0;

    D[2][0] = 0.0;
    D[2][1] = 0.0;
    D[2][2] = factor * (1.0 - nu) / 2.0;
}

// D matrix for Plane Strain
void CQ4Element::DMatrixPlaneStrain(double D[3][3], double E, double nu)
{
    double factor = E / ((1.0 + nu) * (1.0 - 2.0 * nu));

    D[0][0] = factor * (1.0 - nu);
    D[0][1] = factor * nu;
    D[0][2] = 0.0;

    D[1][0] = factor * nu;
    D[1][1] = factor * (1.0 - nu);
    D[1][2] = 0.0;

    D[2][0] = 0.0;
    D[2][1] = 0.0;
    D[2][2] = factor * (1.0 - 2.0 * nu) / 2.0;
}

// Calculate element stiffness matrix
void CQ4Element::ElementStiffness(double* stiffness)
{
    clear(stiffness, SizeOfStiffnessMatrix());
    double D[3][3];
    CQ4Material* mat = static_cast<CQ4Material*>(ElementMaterial_);

    if (mat->analysisType == AnalysisType2D::PlaneStrain) {
        DMatrixPlaneStrain(D, mat->E, mat->nu);
    } else {
        DMatrixPlaneStress(D, mat->E, mat->nu);
    }

    std::cout << "\n*** DEBUG: Stiffness Calculation for Element (nodes: " 
              << nodes_[0]->NodeNumber << "," << nodes_[1]->NodeNumber << "," << nodes_[2]->NodeNumber << "," << nodes_[3]->NodeNumber << ") ***" << std::endl;

    // Gauss quadrature points (2x2)
    const int NGPs = 4;
    double GPx[NGPs] = {-0.577350269189626, 0.577350269189626, 0.577350269189626, -0.577350269189626};
    double GPy[NGPs] = {-0.577350269189626, -0.577350269189626, 0.577350269189626, 0.577350269189626};
    double GPw[NGPs] = {1.0, 1.0, 1.0, 1.0};

    double Ke[8][8]; // Local full stiffness matrix
    for(int i=0; i<8; ++i) for(int j=0; j<8; ++j) Ke[i][j] = 0.0;

    for (int gp = 0; gp < NGPs; ++gp)
    {
        double xi = GPx[gp];
        double eta = GPy[gp];
        double weight = GPw[gp];

        double B[3][8];
        BMatrix(xi, eta, (const CNode**)nodes_, B);

        double J_mat[2][2], invJ[2][2], detJ;
        Jacobian(xi, eta, (const CNode**)nodes_, J_mat, detJ, invJ);

        if (gp == 0) { // Only print for the first GP to reduce clutter
             std::cout << "  Element (nodes: " << nodes_[0]->NodeNumber 
                       << "," << nodes_[1]->NodeNumber << "," << nodes_[2]->NodeNumber << "," << nodes_[3]->NodeNumber << ")"
                       << ", detJ at GP1: " << std::scientific << std::setprecision(5) << detJ << std::endl;
        }

        if (detJ <= 0) {
             std::cerr << "*** WARNING: Non-positive Jacobian determinant in element " 
                       << nodes_[0]->NodeNumber << " at GP " << gp << ". detJ = " << detJ << std::endl;
        }

        // Perform B^T * D * B
        double BD[8][3];
        for (int i = 0; i < 8; ++i) {
            for (int j = 0; j < 3; ++j) {
                BD[i][j] = 0.0;
                for (int k = 0; k < 3; ++k) {
                    BD[i][j] += B[k][i] * D[k][j]; // B is transposed
                }
            }
        }
        
        for (int i = 0; i < 8; ++i) {
            for (int j = 0; j < 8; ++j) {
                double temp = 0.0;
                for (int k = 0; k < 3; ++k) {
                    temp += BD[i][k] * B[k][j];
                }
                Ke[i][j] += temp * detJ * weight * mat->thickness;
            }
        }
    }
    
    std::cout << "  Final Ke matrix for Element (nodes: " << nodes_[0]->NodeNumber << "," << nodes_[1]->NodeNumber << "," << nodes_[2]->NodeNumber << "," << nodes_[3]->NodeNumber << "):" << std::endl;
    for (int i = 0; i < 8; ++i) {
        std::cout << "    ";
        for (int j = 0; j < 8; ++j) {
            std::cout << std::scientific << std::setprecision(3) << std::setw(12) << Ke[i][j];
        }
        std::cout << std::endl;
    }

    // Store upper triangular part of Ke into stiffness array (column-wise)
    unsigned int l = 0;
    for (unsigned int j = 0; j < 8; ++j) { // Columns
        for (unsigned int i = 0; i <= j; ++i) { // Rows
            stiffness[l++] = Ke[i][j];
        }
    }
}

// Calculate element stress/strain
// stress_strain_out array will store: 
// GP1_sig_x, GP1_sig_y, GP1_tau_xy, GP1_eps_x, GP1_eps_y, GP1_gam_xy, (repeated for GP2, GP3, GP4)
// then finally, Average_Von_Mises_Stress
void CQ4Element::ElementStress(double* stress_strain_out, double* Displacement_global)
{
    CQ4Material* mat = static_cast<CQ4Material*>(ElementMaterial_); // Safe static_cast based on design
    if (!mat) {
        // This should ideally not happen if Read function and material setup is correct
        for(unsigned int i=0; i < GetNumberOfStressStrainOutputs(); ++i) stress_strain_out[i] = 0.0;
        return;
    }

    double D[3][3];
    // Choose D matrix based on analysis type
    if (mat->analysisType == AnalysisType2D::PlaneStrain) {
        DMatrixPlaneStrain(D, mat->E, mat->nu); // Corrected call
    } else { // Default to PlaneStress
        DMatrixPlaneStress(D, mat->E, mat->nu); // Corrected call
    }

    // Element nodal displacements (8 values: u1x, u1y, u2x, u2y, ...)
    double Ue[8];
    for (unsigned int i = 0; i < ND_; ++i) {
        if (LocationMatrix_[i] > 0) {
            Ue[i] = Displacement_global[LocationMatrix_[i] - 1];
        } else {
            Ue[i] = 0.0; // Fixed DOF
        }
    }
    
    const int NGPs = 4; // Number of Gauss points
    double GP_coords[NGPs][2] = {
        {-0.577350269189626, -0.577350269189626},
        { 0.577350269189626, -0.577350269189626},
        { 0.577350269189626,  0.577350269189626},
        {-0.577350269189626,  0.577350269189626}
    };

    const CNode* current_nodes[4] = {nodes_[0], nodes_[1], nodes_[2], nodes_[3]};
    double sum_von_mises_stress = 0.0;

    for (int gp = 0; gp < NGPs; ++gp) {
        double xi = GP_coords[gp][0];
        double eta = GP_coords[gp][1];

        double B[3][8];
        BMatrix(xi, eta, current_nodes, B); // BMatrix is already implemented

        // Strain: epsilon = B * Ue  (epsilon is [eps_x, eps_y, gamma_xy]^T)
        double epsilon[3] = {0.0, 0.0, 0.0};
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 8; ++j) {
                epsilon[i] += B[i][j] * Ue[j];
            }
        }

        // Stress: sigma = D * epsilon (sigma is [sigma_x, sigma_y, tau_xy]^T)
        double sigma[3] = {0.0, 0.0, 0.0};
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                sigma[i] += D[i][j] * epsilon[j];
            }
        }

        // Store results for this Gauss point
        int base_idx = gp * 6;
        stress_strain_out[base_idx + 0] = sigma[0];   // sigma_x
        stress_strain_out[base_idx + 1] = sigma[1];   // sigma_y
        stress_strain_out[base_idx + 2] = sigma[2];   // tau_xy
        stress_strain_out[base_idx + 3] = epsilon[0]; // epsilon_x
        stress_strain_out[base_idx + 4] = epsilon[1]; // epsilon_y
        stress_strain_out[base_idx + 5] = epsilon[2]; // gamma_xy

        // Accumulate for average Von Mises stress
        double s_x = sigma[0];
        double s_y = sigma[1];
        double t_xy = sigma[2];
        sum_von_mises_stress += std::sqrt(s_x*s_x - s_x*s_y + s_y*s_y + 3*t_xy*t_xy);
    }

    // Store average Von Mises stress at the end
    stress_strain_out[NGPs * 6] = sum_von_mises_stress / NGPs;
}

// Generate location matrix for Q4 element (2 DOFs per node: UX, UY)
void CQ4Element::GenerateLocationMatrix()
{
    unsigned int i = 0;
    for (unsigned int N = 0; N < NEN_; N++)
    {
        // For Q4, we only consider 2 DOFs (X and Y) per node
        for (unsigned int D = 0; D < 2; D++) 
        {
            LocationMatrix_[i++] = nodes_[N]->bcode[D];
        }
    }
}
