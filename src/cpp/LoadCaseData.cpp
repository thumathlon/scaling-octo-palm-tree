/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "LoadCaseData.h"

#include <iomanip>
#include <iostream>

using namespace std;

CLoadCaseData :: ~CLoadCaseData()
{
	delete [] node;
	delete [] dof;
	delete [] load;
}

void CLoadCaseData :: Allocate(unsigned int num)
{
	nloads = num;
	node = new unsigned int[nloads];
	dof = new unsigned int[nloads];
	load = new double[nloads];
}; 

//	Read load case data from stream Input
bool CLoadCaseData :: Read(ifstream& Input)
{
//	Load case number (LL) is read by CDomain::ReadLoadCases() prior to this call.
//	This function reads number of concentrated loads in this load case (NL)
	
	unsigned int NL; // Number of loads in *this* specific load case

	Input >> NL;
    if (Input.fail()) {
        cerr << "*** Error (CLoadCaseData::Read): Failed to read number of concentrated loads (NL)." << endl;
        return false;
    }

    if (NL > 0) { // Only allocate if there are loads
	    Allocate(NL); // Allocate calls new, could fail
        if (!node || !dof || !load) { // Basic check for allocation failure
            cerr << "*** Error (CLoadCaseData::Read): Memory allocation failed for loads (NL = " << NL << ")." << endl;
            // Cleanup already allocated memory if any (tricky if Allocate does partial alloc)
            // For simplicity, assuming Allocate handles its own failures or nloads remains 0 if it fails. Or it fully allocates or nothing.
            // A more robust Allocate would throw or return status.
            nloads = 0; // Prevent access to bad pointers
            delete [] node; node = nullptr;
            delete [] dof; dof = nullptr;
            delete [] load; load = nullptr;
            return false;
        }

	    for (unsigned int i = 0; i < NL; i++)
        {
		    Input >> node[i] >> dof[i] >> load[i];
            if (Input.fail()) {
                cerr << "*** Error (CLoadCaseData::Read): Failed to read data for load entry " << i + 1 << " (out of " << NL << ")." << endl;
                return false;
            }
            // Add validation for node[i], dof[i] if necessary
            // e.g., if (node[i] == 0 || node[i] > MAX_NODE_NUMBER) { cerr << "Invalid node " << node[i] << endl; return false; }
            // e.g., if (dof[i] == 0 || dof[i] > MAX_DOF_PER_NODE) { cerr << "Invalid dof " << dof[i] << endl; return false; }
        }
        cout << "    Successfully read " << NL << " load entries for this load case." << endl;
    }
    else {
        nloads = 0; // Ensure nloads is 0 if NL read from file is 0
        cout << "    No concentrated loads (NL = 0) for this load case." << endl;
    }

	return true;
}

//	Write load case data to stream
void CLoadCaseData::Write(COutputter& output)
{
	for (unsigned int i = 0; i < nloads; i++)
		output << setw(7) << node[i] << setw(13) << dof[i]  << setw(19) << load[i] << endl;
}
