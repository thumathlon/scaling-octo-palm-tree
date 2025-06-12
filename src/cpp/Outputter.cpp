/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include <ctime>
#include <vector>

#include "Domain.h"
#include "Outputter.h"
#include "SkylineMatrix.h"
#include "CQ4Element.h"

using namespace std;

//	Output current time and date
void COutputter::PrintTime(const struct tm* ptm, COutputter &output)
{
	const char* weekday[] = {"Sunday", "Monday", "Tuesday", "Wednesday",
							 "Thursday", "Friday", "Saturday"};
	const char* month[] = {"January", "February", "March", "April", "May", "June",
						   "July", "August", "September", "October", "November", "December"};

	output << "        (";
	output << ptm->tm_hour << ":" << ptm->tm_min << ":" << ptm->tm_sec << " on ";
	output << month[ptm->tm_mon] << " " << ptm->tm_mday << ", " << ptm->tm_year + 1900 << ", "
		   << weekday[ptm->tm_wday] << ")" << endl
		   << endl;
}

COutputter* COutputter::_instance = nullptr;

//	Constructor
COutputter::COutputter(string FileName)
{
	OutputFile.open(FileName);

	if (!OutputFile)
	{
		cerr << "*** Error *** File " << FileName << " does not exist !" << endl;
		exit(3);
	}
}

//	Return the single instance of the class
COutputter* COutputter::GetInstance(string FileName)
{
	if (!_instance)
		_instance = new COutputter(FileName);
    
	return _instance;
}

//	Print program logo
void COutputter::OutputHeading()
{
	CDomain* FEMData = CDomain::GetInstance();

	*this << "TITLE : " << FEMData->GetTitle() << endl;

	time_t rawtime;
	struct tm timeinfo_data;
	struct tm* timeinfo = &timeinfo_data;

	time(&rawtime);
	localtime_s(timeinfo, &rawtime);

	PrintTime(timeinfo, *this);
}

//	Print nodal data
void COutputter::OutputNodeInfo()
{
	CDomain* FEMData = CDomain::GetInstance();

	CNode* NodeList = FEMData->GetNodeList();

	*this << "C O N T R O L   I N F O R M A T I O N" << endl
		  << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	unsigned int NUMNP = FEMData->GetNUMNP();
	unsigned int NUMEG = FEMData->GetNUMEG();
	unsigned int NLCASE = FEMData->GetNLCASE();
	unsigned int MODEX = FEMData->GetMODEX();

	*this << "      NUMBER OF NODAL POINTS . . . . . . . . . . (NUMNP)  =" << setw(6) << NUMNP << endl;
	*this << "      NUMBER OF ELEMENT GROUPS . . . . . . . . . (NUMEG)  =" << setw(6) << NUMEG << endl;
	*this << "      NUMBER OF LOAD CASES . . . . . . . . . . . (NLCASE) =" << setw(6) << NLCASE << endl;
	*this << "      SOLUTION MODE  . . . . . . . . . . . . . . (MODEX)  =" << setw(6) << MODEX << endl;
	*this << "         EQ.0, DATA CHECK" << endl
		  << "         EQ.1, EXECUTION" << endl
		  << endl;

	*this << " N O D A L   P O I N T   D A T A" << endl << endl;
	*this << "    NODE       BOUNDARY                         NODAL POINT" << endl
		  << "   NUMBER  CONDITION  CODES                     COORDINATES" << endl;

	for (unsigned int np = 0; np < NUMNP; np++)
		NodeList[np].Write(*this);

	*this << endl;
}

//	Output equation numbers
void COutputter::OutputEquationNumber()
{
	CDomain* FEMData = CDomain::GetInstance();
	unsigned int NUMNP = FEMData->GetNUMNP();

	CNode* NodeList = FEMData->GetNodeList();

	*this << " EQUATION NUMBERS" << endl
		  << endl;
	*this << "   NODE NUMBER   DEGREES OF FREEDOM" << endl;
	*this << "        N           X    Y    Z" << endl;

	for (unsigned int np = 0; np < NUMNP; np++) // Loop over for all node
		NodeList[np].WriteEquationNo(*this);

	*this << endl;
}

//	Output element data
void COutputter::OutputElementInfo()
{
	//	Print element group control line

	CDomain* FEMData = CDomain::GetInstance();

	unsigned int NUMEG = FEMData->GetNUMEG();

	*this << " E L E M E N T   G R O U P   D A T A" << endl
		  << endl
		  << endl;

	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
		*this << " E L E M E N T   D E F I N I T I O N" << endl
			  << endl;

		ElementTypes ElementType = FEMData->GetEleGrpList()[EleGrp].GetElementType();
		unsigned int NUME = FEMData->GetEleGrpList()[EleGrp].GetNUME();

		*this << " ELEMENT TYPE  . . . . . . . . . . . . .( NPAR(1) ) . . =" << setw(5)
			  << ElementType << endl;
		*this << "     EQ.1, TRUSS ELEMENTS" << endl
			  << "     EQ.2, ELEMENTS CURRENTLY" << endl
			  << "     EQ.3, NOT AVAILABLE" << endl
			  << endl;

		*this << " NUMBER OF ELEMENTS. . . . . . . . . . .( NPAR(2) ) . . =" << setw(5) << NUME
			  << endl
			  << endl;

		switch (ElementType)
		{
			case ElementTypes::Bar: // Bar element
				OutputBarElements(EleGrp);
				break;
            case ElementTypes::Q4: // Q4 element
                // Implemented, but no specific output routine here yet.
                break;
		    default:
		        *this << ElementType << " has not been implemented yet." << endl;
		        break;
		}
	}
}
//	Output bar element data
void COutputter::OutputBarElements(unsigned int EleGrp)
{
	CDomain* FEMData = CDomain::GetInstance();

	CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
	unsigned int NUMMAT = ElementGroup.GetNUMMAT();

	*this << " M A T E R I A L   D E F I N I T I O N" << endl
		  << endl;
	*this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << endl;
	*this << " AND CROSS-SECTIONAL  CONSTANTS  . . . .( NPAR(3) ) . . =" << setw(5) << NUMMAT
		  << endl
		  << endl;

	*this << "  SET       YOUNG'S     CROSS-SECTIONAL" << endl
		  << " NUMBER     MODULUS          AREA" << endl
		  << "               E              A" << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	//	Loop over for all property sets
	for (unsigned int mset = 0; mset < NUMMAT; mset++)
    {
        *this << setw(5) << mset+1;
		dynamic_cast<CBarMaterial&>(ElementGroup.GetMaterial(mset)).Write(*this);
    }

	*this << endl << endl
		  << " E L E M E N T   I N F O R M A T I O N" << endl;
    
	*this << " ELEMENT     NODE     NODE       MATERIAL" << endl
		  << " NUMBER-N      I        J       SET NUMBER" << endl;

	unsigned int NUME = ElementGroup.GetNUME();

	//	Loop over for all elements in group EleGrp
	for (unsigned int Ele = 0; Ele < NUME; Ele++)
    {
        *this << setw(5) << Ele+1;
		ElementGroup[Ele].Write(*this);
    }

	*this << endl;
}

//	Print load data
void COutputter::OutputLoadInfo()
{
	CDomain* FEMData = CDomain::GetInstance();

	for (unsigned int lcase = 1; lcase <= FEMData->GetNLCASE(); lcase++)
	{
		CLoadCaseData* LoadData = &FEMData->GetLoadCases()[lcase - 1];

		*this << setiosflags(ios::scientific);
		*this << " L O A D   C A S E   D A T A" << endl
			  << endl;

		*this << "     LOAD CASE NUMBER . . . . . . . =" << setw(6) << lcase << endl;
		*this << "     NUMBER OF CONCENTRATED LOADS . =" << setw(6) << LoadData->nloads << endl
			  << endl;
		*this << "    NODE       DIRECTION      LOAD" << endl
			  << "   NUMBER                   MAGNITUDE" << endl;

		LoadData->Write(*this);

		*this << endl;
	}
}

//	Print nodal displacement
void COutputter::OutputNodalDisplacement()
{
	CDomain* FEMData = CDomain::GetInstance();
	CNode* NodeList = FEMData->GetNodeList();
	double* Displacement = FEMData->GetDisplacement();

	*this << setiosflags(ios::scientific);

	*this << " D I S P L A C E M E N T S" << endl
		  << endl;
	*this << "  NODE           X-DISPLACEMENT    Y-DISPLACEMENT    Z-DISPLACEMENT" << endl;

	for (unsigned int np = 0; np < FEMData->GetNUMNP(); np++)
		NodeList[np].WriteNodalDisplacement(*this, Displacement);

	*this << endl;
}

//	Calculate stresses
void COutputter::OutputElementStress()
{
	CDomain* FEMData = CDomain::GetInstance();
	double* Displacement = FEMData->GetDisplacement();
	unsigned int NUMEG = FEMData->GetNUMEG();

	*this << setiosflags(ios::scientific) << setprecision(5);

	for (unsigned int EleGrpIndex = 0; EleGrpIndex < NUMEG; EleGrpIndex++)
	{
		*this << " S T R E S S  C A L C U L A T I O N S  F O R  E L E M E N T  G R O U P" << setw(5)
			  << EleGrpIndex + 1 << endl
			  << endl;

		CElementGroup& EleGrp = FEMData->GetEleGrpList()[EleGrpIndex];
		unsigned int NUME = EleGrp.GetNUME();
		ElementTypes ElementType = EleGrp.GetElementType(); // Keep for now, might be useful for type-specific headers
		CElement& firstElement = EleGrp[0]; // Get a representative element to query output size
		unsigned int numStressValues = firstElement.GetNumberOfStressStrainOutputs();
		std::vector<double> stress_strain_values(numStressValues);

		// Output headers based on element type
		// This part remains somewhat specific, or needs more generic header generation
		if (dynamic_cast<CBar*>(&firstElement)) // Check if it's a Bar element
		{
			*this << "  ELEMENT             FORCE            STRESS" << endl
				  << "  NUMBER" << endl;
		}
		else if (dynamic_cast<CQ4Element*>(&firstElement)) // Check if it's a Q4 element
		{
			*this << "  ELEMENT    GP    Sigma_X    Sigma_Y    Tau_XY     Eps_X      Eps_Y      Gam_XY" << endl;
			*this << "  NUMBER" << endl;
			*this << "  ----------------------------------------------------------------------------------" << endl;
			// Specific header for the last value (Average Von Mises)
			// Will be printed after GP data for each element
		}
		else
		{
			*this << "  ELEMENT    STRESS/STRAIN VALUES (Count: " << numStressValues << ")" << endl
				  << "  NUMBER" << endl;
		}

		for (unsigned int Ele = 0; Ele < NUME; Ele++)
		{
			CElement& Element = EleGrp[Ele];
			Element.ElementStress(stress_strain_values.data(), Displacement);

			*this << setw(5) << Ele + 1;

			if (dynamic_cast<CBar*>(&Element)) // Bar specific output
			{
				CBarMaterial* material = dynamic_cast<CBarMaterial*>(Element.GetElementMaterial());
				if (material && numStressValues == 1) {
					double stress = stress_strain_values[0];
					*this << setw(22) << stress * material->Area << setw(18)
						  << stress << endl;
				}
			}
			else if (dynamic_cast<CQ4Element*>(&Element)) // Q4 specific output
			{
				if (numStressValues == 25) { // 4 GPs * 6 values/GP + 1 avg VonMises
					// Output for each Gauss Point (4 GPs, 6 values each)
					for (int gp = 0; gp < 4; ++gp) {
						if (gp > 0) *this << setw(5) << " "; // Indent subsequent GP lines
						*this << setw(7) << gp + 1; // Gauss Point number
						for (int val_idx = 0; val_idx < 6; ++val_idx) {
							*this << setw(12) << stress_strain_values[gp * 6 + val_idx];
						}
						*this << endl;
					}
					// Output Average Von Mises Stress
					*this << setw(5) << " " << " AvgVM: " << setw(12) << stress_strain_values[24] << endl;
				} else {
					 // Fallback for unexpected number of values from Q4
					for (unsigned int i = 0; i < numStressValues; ++i) {
						*this << setw(15) << stress_strain_values[i];
					}
					*this << endl;
				}
			}
			else // Generic output for other element types
			{
				for (unsigned int i = 0; i < numStressValues; ++i) {
					*this << setw(15) << stress_strain_values[i];
				}
				*this << endl;
			}
		}
		*this << endl;
	}
}

//	Print total system data
void COutputter::OutputTotalSystemData()
{
	CDomain* FEMData = CDomain::GetInstance();

	*this << "	TOTAL SYSTEM DATA" << endl
		  << endl;

	*this << "     NUMBER OF EQUATIONS . . . . . . . . . . . . . .(NEQ) = " << FEMData->GetNEQ()
		  << endl
		  << "     NUMBER OF MATRIX ELEMENTS . . . . . . . . . . .(NWK) = " << FEMData->GetStiffnessMatrix()->size()
		  << endl
		  << "     MAXIMUM HALF BANDWIDTH  . . . . . . . . . . . .(MK ) = " << FEMData->GetStiffnessMatrix()->GetMaximumHalfBandwidth()
		  << endl
		  << "     MEAN HALF BANDWIDTH . . . . . . . . . . . . . .(MM ) = " << FEMData->GetStiffnessMatrix()->size() / FEMData->GetNEQ() << endl
		  << endl
		  << endl;
}

#ifdef _DEBUG_

//	Print column heights for debuging
void COutputter::PrintColumnHeights()
{
	*this << "*** _Debug_ *** Column Heights" << endl;

	CDomain* FEMData = CDomain::GetInstance();

	unsigned int NEQ = FEMData->GetNEQ();
	CSkylineMatrix<double> *StiffnessMatrix = FEMData->GetStiffnessMatrix();
	unsigned int* ColumnHeights = StiffnessMatrix->GetColumnHeights();

	for (unsigned int col = 0; col < NEQ; col++)
	{
		if (col + 1 % 10 == 0)
		{
			*this << endl;
		}

		*this << setw(8) << ColumnHeights[col];
	}

	*this << endl
		  << endl;
}

//	Print address of diagonal elements for debuging
void COutputter::PrintDiagonalAddress()
{
	*this << "*** _Debug_ *** Address of Diagonal Element" << endl;

	CDomain* FEMData = CDomain::GetInstance();

	unsigned int NEQ = FEMData->GetNEQ();
	CSkylineMatrix<double> *StiffnessMatrix = FEMData->GetStiffnessMatrix();
	unsigned int* DiagonalAddress = StiffnessMatrix->GetDiagonalAddress();

	for (unsigned int col = 0; col <= NEQ; col++)
	{
		if (col + 1 % 10 == 0)
		{
			*this << endl;
		}

		*this << setw(8) << DiagonalAddress[col];
	}

	*this << endl
		  << endl;
}

//	Print banded and full stiffness matrix for debuging
void COutputter::PrintStiffnessMatrix()
{
	*this << "*** _Debug_ *** Banded stiffness matrix" << endl;

	CDomain* FEMData = CDomain::GetInstance();

	unsigned int NEQ = FEMData->GetNEQ();
	CSkylineMatrix<double> *StiffnessMatrix = FEMData->GetStiffnessMatrix();
	unsigned int* DiagonalAddress = StiffnessMatrix->GetDiagonalAddress();

	*this << setiosflags(ios::scientific) << setprecision(5);

	for (unsigned int i = 0; i < DiagonalAddress[NEQ] - DiagonalAddress[0]; i++)
	{
		*this << setw(14) << (*StiffnessMatrix)(i);

		if ((i + 1) % 6 == 0)
		{
			*this << endl;
		}
	}

	*this << endl
		  << endl;

	*this << "*** _Debug_ *** Full stiffness matrix" << endl;

	for (int I = 1; I <= NEQ; I++)
	{
		for (int J = 1; J <= NEQ; J++)
		{
			int J_new = (J > I) ? J : I;
			int I_new = (J > I) ? I : J;
			int H = DiagonalAddress[J_new] - DiagonalAddress[J_new - 1];
			if (J_new - I_new - H >= 0)
			{
				*this << setw(14) << 0.0;
			}
			else
			{
				*this << setw(14) << (*StiffnessMatrix)(I_new, J_new);
			}
		}

		*this << endl;
	}

	*this << endl;
}

//	Print displacement vector for debuging
void COutputter::PrintDisplacement()
{
	*this << "*** _Debug_ *** Displacement vector" << endl;

	CDomain* FEMData = CDomain::GetInstance();

	unsigned int NEQ = FEMData->GetNEQ();
	double* Force = FEMData->GetForce();

	*this << setiosflags(ios::scientific) << setprecision(5);

	for (unsigned int i = 0; i < NEQ; i++)
	{
		if ((i + 1) % 6 == 0)
		{
			*this << endl;
		}

		*this << setw(14) << Force[i];
	}

	*this << endl
		  << endl;
}

#endif

//! Write results to a legacy VTK file for post-processing
void COutputter::WriteLegacyVTKFile(const std::string& filename)
{
    CDomain* domain = CDomain::GetInstance();
    if (!domain) return;

    std::ofstream vtk_file(filename);
    if (!vtk_file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }

    unsigned int num_nodes = domain->GetNUMNP();
    CNode* nodes = domain->GetNodeList();
    double* displacements = domain->GetDisplacement();

    // --- Calculate Nodal Stresses (temporary, not stored in CNode) ---
    struct NodalStressData {
        double sxx = 0.0, syy = 0.0, sxy = 0.0;
        int count = 0;
    };
    std::vector<NodalStressData> nodal_stresses(num_nodes);
    std::vector<double> von_mises_stresses(num_nodes, 0.0);

    for (unsigned int i = 0; i < domain->GetNUMEG(); ++i) {
        CElementGroup& group = domain->GetEleGrpList()[i];
        if (dynamic_cast<CQ4Element*>(&group[0])) {
            for (unsigned int j = 0; j < group.GetNUME(); ++j) {
                CQ4Element* elem = static_cast<CQ4Element*>(&group[j]);
                double stress_strain_out[25];
                elem->ElementStress(stress_strain_out, displacements);

                double avg_sxx = 0, avg_syy = 0, avg_sxy = 0;
                for (int gp = 0; gp < 4; ++gp) {
                    avg_sxx += stress_strain_out[gp * 6 + 0];
                    avg_syy += stress_strain_out[gp * 6 + 1];
                    avg_sxy += stress_strain_out[gp * 6 + 2];
                }
                avg_sxx /= 4.0;
                avg_syy /= 4.0;
                avg_sxy /= 4.0;

                for (unsigned int k = 0; k < elem->GetNEN(); ++k) {
                    unsigned int node_idx = elem->GetNodes()[k]->NodeNumber - 1;
                    nodal_stresses[node_idx].sxx += avg_sxx;
                    nodal_stresses[node_idx].syy += avg_syy;
                    nodal_stresses[node_idx].sxy += avg_sxy;
                    nodal_stresses[node_idx].count++;
                }
            }
        }
    }

    for (unsigned int i = 0; i < num_nodes; ++i) {
        if (nodal_stresses[i].count > 0) {
            double sxx = nodal_stresses[i].sxx / nodal_stresses[i].count;
            double syy = nodal_stresses[i].syy / nodal_stresses[i].count;
            double sxy = nodal_stresses[i].sxy / nodal_stresses[i].count;
            von_mises_stresses[i] = std::sqrt(sxx*sxx - sxx*syy + syy*syy + 3*sxy*sxy);
        }
    }

    // --- Write VTK File ---
    vtk_file << "# vtk DataFile Version 3.0" << std::endl;
    vtk_file << "STAP++ Post-processing file" << std::endl;
    vtk_file << "ASCII" << std::endl;
    vtk_file << "DATASET UNSTRUCTURED_GRID" << std::endl;

    // POINTS (Deformed)
    vtk_file << "POINTS " << num_nodes << " double" << std::endl;
    for (unsigned int i = 0; i < num_nodes; ++i) {
        double dx = nodes[i].bcode[0] > 0 ? displacements[nodes[i].bcode[0] - 1] : 0.0;
        double dy = nodes[i].bcode[1] > 0 ? displacements[nodes[i].bcode[1] - 1] : 0.0;
        vtk_file << nodes[i].XYZ[0] + dx << " " << nodes[i].XYZ[1] + dy << " " << nodes[i].XYZ[2] << std::endl;
    }

    // CELLS
    unsigned int num_q4_elements = 0;
    for (unsigned int i = 0; i < domain->GetNUMEG(); ++i) {
        CElementGroup& group = domain->GetEleGrpList()[i];
        if (dynamic_cast<CQ4Element*>(&group[0])) {
            num_q4_elements += group.GetNUME();
        }
    }
    vtk_file << "CELLS " << num_q4_elements << " " << num_q4_elements * 5 << std::endl;
    for (unsigned int i = 0; i < domain->GetNUMEG(); ++i) {
        CElementGroup& group = domain->GetEleGrpList()[i];
        if (dynamic_cast<CQ4Element*>(&group[0])) {
            for (unsigned int j = 0; j < group.GetNUME(); j++) {
                CElement& elem = group[j];
                vtk_file << "4";
                for (unsigned int k = 0; k < elem.GetNEN(); k++) {
                    vtk_file << " " << elem.GetNodes()[k]->NodeNumber - 1;
                }
                vtk_file << std::endl;
            }
        }
    }

    // CELL_TYPES
    vtk_file << "CELL_TYPES " << num_q4_elements << std::endl;
    for (unsigned int i = 0; i < num_q4_elements; ++i) {
        vtk_file << "9" << std::endl; // VTK_QUAD
    }

    // POINT_DATA
    vtk_file << "POINT_DATA " << num_nodes << std::endl;
    // Displacements
    vtk_file << "VECTORS Displacements double" << std::endl;
    for (unsigned int i = 0; i < num_nodes; ++i) {
        double dx = nodes[i].bcode[0] > 0 ? displacements[nodes[i].bcode[0] - 1] : 0.0;
        double dy = nodes[i].bcode[1] > 0 ? displacements[nodes[i].bcode[1] - 1] : 0.0;
        vtk_file << dx << " " << dy << " 0.0" << std::endl;
    }
    // Von Mises Stress
    vtk_file << "SCALARS VonMisesStress double 1" << std::endl;
    vtk_file << "LOOKUP_TABLE default" << std::endl;
    for (unsigned int i = 0; i < num_nodes; ++i) {
        vtk_file << von_mises_stresses[i] << std::endl;
    }

    vtk_file.close();
}
