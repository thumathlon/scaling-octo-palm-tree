/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include <iostream>
#include <iomanip>

#include "Node.h"

CNode::CNode(double X, double Y, double Z)
{
    XYZ[0] = X;		// Coordinates of the node
    XYZ[1] = Y;
    XYZ[2] = Z;
    
    bcode[0] = 0;	// Boundary codes
    bcode[1] = 0;
    bcode[2] = 0;
};

//	Read element data from stream Input
bool CNode::Read(ifstream& Input)
{
    cout << "DEBUG: Starting CNode::Read" << endl;
    cout << "DEBUG: Current stream state before reading NodeNumber: failbit=" << Input.fail() 
         << ", badbit=" << Input.bad() << ", eofbit=" << Input.eof() << endl;
    
    // 打印下一个字符
    char nextChar = Input.peek();
    cout << "DEBUG: Next character to read: '" << nextChar << "' (ASCII: " << (int)nextChar << ")" << endl;

    Input >> NodeNumber;
    cout << "DEBUG: Attempted to read NodeNumber: " << NodeNumber << endl;
    cout << "DEBUG: Stream state after reading NodeNumber: failbit=" << Input.fail() 
         << ", badbit=" << Input.bad() << ", eofbit=" << Input.eof() << endl;

    if (Input.fail() || NodeNumber == 0) {
        cerr << "*** Error (CNode::Read): Failed to read NodeNumber or NodeNumber is 0." << endl;
        cerr << "NodeNumber: " << NodeNumber << endl;
		return false;
	}
	
    cout << "DEBUG: Attempting to read boundary codes and coordinates..." << endl;
	Input >> bcode[0] >> bcode[1] >> bcode[2]
		  >> XYZ[0] >> XYZ[1] >> XYZ[2];
    
    cout << "DEBUG: Read values: bcode=[" << bcode[0] << "," << bcode[1] << "," << bcode[2] 
         << "], XYZ=[" << XYZ[0] << "," << XYZ[1] << "," << XYZ[2] << "]" << endl;
    
    cout << "DEBUG: Stream state after reading all values: failbit=" << Input.fail() 
         << ", badbit=" << Input.bad() << ", eofbit=" << Input.eof() << endl;

    if (Input.fail()) {
        cerr << "*** Error (CNode::Read): Failed to read boundary codes or coordinates for Node " << NodeNumber << endl;
        return false;
    }

    cout << "DEBUG: Successfully read node " << NodeNumber << endl;
	return true;
}

//	Output nodal point data to stream
void CNode::Write(COutputter& output)
{
	output << setw(9) << NodeNumber << setw(5) << bcode[0] << setw(5) << bcode[1] << setw(5) << bcode[2]
		   << setw(18) << XYZ[0] << setw(15) << XYZ[1] << setw(15) << XYZ[2] << endl;
}

//	Output equation numbers of nodal point to stream
void CNode::WriteEquationNo(COutputter& output)
{
	output << setw(9) << NodeNumber << "       ";

	for (unsigned int dof = 0; dof < CNode::NDF; dof++)	// Loop over for DOFs of node np
	{
		output << setw(5) << bcode[dof];
	}

	output << endl;
}

//	Write nodal displacement
void CNode::WriteNodalDisplacement(COutputter& output, double* Displacement)
{
	output << setw(5) << NodeNumber << "        ";

	for (unsigned int j = 0; j < NDF; j++)
	{
		if (bcode[j] == 0)
		{
			output << setw(18) << 0.0;
		}
		else
		{
			output << setw(18) << Displacement[bcode[j] - 1];
		}
	}

	output << endl;
}
