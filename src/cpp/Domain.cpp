/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Domain.h"
#include "Material.h"
#include <iostream>
#include <string>

using namespace std;

//	Clear an array
//  THIS IS THE BUG. The implementation of a template function must be in the header file.
//  I will move this to Element.h
/*
template <class type> void clear( type* a, unsigned int N )
{
	for (unsigned int i = 0; i < N; i++)
		a[i] = (type) 0;
}
*/

CDomain* CDomain::_instance = nullptr;

//	Constructor
CDomain::CDomain()
{
	Title[0] = '0';
	MODEX = 0;

	NUMNP = 0;
	NodeList = nullptr;
	
	NUMEG = 0;
	EleGrpList = nullptr;
	
	NLCASE = 0;
	NLOAD = nullptr;
	LoadCases = nullptr;
	
	NEQ = 0;

	Force = nullptr;
	StiffnessMatrix = nullptr;
}

//	Desconstructor
CDomain::~CDomain()
{
	delete [] NodeList;

	delete [] EleGrpList;

	delete [] NLOAD;
	delete [] LoadCases;

	delete [] Force;
	delete StiffnessMatrix;
}

//	Return pointer to the instance of the Domain class
CDomain* CDomain::GetInstance()
{
	if (!_instance) 
		_instance = new CDomain();
	
	return _instance;
}

//	Read domain data from the input data file
bool CDomain::ReadData(string FileName, string OutFile)
{
	Input.open(FileName);

	if (!Input.is_open()) // Explicitly check if file is open
	{
		cerr << "*** Error *** Failed to open input file: " << FileName << endl;
		return false; // Return false instead of exit, allows caller to handle
	}

	COutputter* Output = COutputter::GetInstance(OutFile);
    if (!Output) { // Ensure Outputter instance is successfully obtained
        cerr << "*** Error *** Failed to initialize outputter. Output file might be invalid: " << OutFile << endl;
        Input.close(); // Close already opened input file
        return false;
    }

//	Read the heading line
	Input.getline(Title, 256);
    // DEBUG LINES START
    cout << "DEBUG: Read Title line: [\\\"" << Title << "\\\"]" << endl; 
    cout << "DEBUG: Input stream state after getline (before Output->OutputHeading): failbit=" << Input.fail() 
         << ", badbit=" << Input.bad() << ", eofbit=" << Input.eof() << endl;
    // DEBUG LINES END

    if (Input.fail()) { // MODIFIED IF BLOCK START
        cerr << "*** Error *** Failed to read title line from input file. Stream state will be cleared to attempt control line read." << endl;
        Input.clear(); // Clear stream state, very important
    } // MODIFIED IF BLOCK END
    
    // If Output->OutputHeading() itself operates on the input stream, that could complicate things,
    // but typically it only handles output.
	Output->OutputHeading();
    // DEBUG LINES START
    cout << "DEBUG: Input stream state after Output->OutputHeading (before control line read): failbit=" << Input.fail() 
         << ", badbit=" << Input.bad() << ", eofbit=" << Input.eof() << endl;
    // DEBUG LINES END


//	Read the control line
    // DEBUG LINES START
	cout << "DEBUG: Attempting to read control line (NUMNP, NUMEG, NLCASE, MODEX)..." << endl; 
    // DEBUG LINES END
	Input >> NUMNP >> NUMEG >> NLCASE >> MODEX; // Attempt to read
    // DEBUG LINES START
    cout << "DEBUG: Values after attempting to read control line: " 
         << "NUMNP=" << NUMNP << ", NUMEG=" << NUMEG 
         << ", NLCASE=" << NLCASE << ", MODEX=" << MODEX << endl; // Print values after read attempt
    cout << "DEBUG: Input stream state after control line read attempt: failbit=" << Input.fail() 
         << ", badbit=" << Input.bad() << ", eofbit=" << Input.eof() << endl; // Print stream state
    // DEBUG LINES END

    if (Input.fail()) { // This is the existing error check for the control line
        cerr << "*** Error *** Failed to read control line (NUMNP, NUMEG, NLCASE, MODEX) from input file." << endl;
        Input.close();
        return false;
    }
    if (NUMNP == 0) { // This is an existing warning
        cerr << "*** Warning *** Number of nodal points (NUMNP) is 0. Check input file." << endl;
        // Decide if this is a fatal error, warning only for now
    }


//	Read nodal point data
	if (ReadNodalPoints())
        Output->OutputNodeInfo();
    else
    {
        cerr << "*** Error *** Failed during reading nodal points." << endl; // Pinpoint error location
        Input.close();
        return false;
    }

//	Update equation number
	CalculateEquationNumber();
	Output->OutputEquationNumber();

//	Read load data
	if (ReadLoadCases())
        Output->OutputLoadInfo();
    else
    {
        cerr << "*** Error *** Failed during reading load cases." << endl; // Pinpoint error location
        Input.close();
        return false;
    }

//	Read element data
	if (ReadElements())
        Output->OutputElementInfo();
    else
    {
        cerr << "*** Error *** Failed during reading element data." << endl; // Pinpoint error location
        Input.close();
        return false;
    }

    Input.close(); // Ensure file is closed after successful read
	return true;
}

//	Read nodal point data
bool CDomain::ReadNodalPoints()
{
    cout << "DEBUG: Starting ReadNodalPoints" << endl;
    if (NUMNP == 0) {
        cout << "DEBUG: NUMNP is 0, skipping nodal point reading" << endl;
        return true;
    }

    string dummyLine;
    cout << "DEBUG: Attempting to read first line before node data" << endl;
    
    // Skip the "NODAL_POINTS" line
    if (!getline(Input, dummyLine)) { 
        cerr << "*** Error *** Failed to read/skip 'NODAL_POINTS' header line before node data." << endl;
        return false; 
    }
    cout << "DEBUG: Skipped line 1 (expected NODAL_POINTS header): [\\\"" << dummyLine << "\\\"]" << endl;
    cout << "DEBUG: Stream state after reading first line: failbit=" << Input.fail() 
         << ", badbit=" << Input.bad() << ", eofbit=" << Input.eof() << endl;

    // Skip the comment line "// NodeID..."
    if (!getline(Input, dummyLine)) {
        cerr << "*** Error *** Failed to read/skip comment line before node data." << endl;
        return false;
    }
    cout << "DEBUG: Skipped line 2 (expected comment line): [\\\"" << dummyLine << "\\\"]" << endl;
    cout << "DEBUG: Stream state after reading second line: failbit=" << Input.fail() 
         << ", badbit=" << Input.bad() << ", eofbit=" << Input.eof() << endl;

    // 清除可能的错误状态
    Input.clear();
    
    // 跳过可能的空白字符
    while (Input.peek() == ' ' || Input.peek() == '\t' || Input.peek() == '\n') {
        char c = Input.get();
        cout << "DEBUG: Skipped whitespace character: '" << c << "' (ASCII: " << (int)c << ")" << endl;
    }

    cout << "DEBUG: Allocating NodeList for " << NUMNP << " nodes" << endl;
	NodeList = new CNode[NUMNP];
    if (!NodeList) {
        cerr << "*** Error *** Memory allocation failed for NodeList (NUMNP = " << NUMNP << ")" << endl;
        return false;
    }

    cout << "DEBUG: Starting to read " << NUMNP << " nodes" << endl;
	for (unsigned int np = 0; np < NUMNP; np++)
    {
        cout << "DEBUG: Attempting to read node " << np + 1 << endl;
		if (!NodeList[np].Read(Input))
        {
            cerr << "*** Error *** Failed to read data for nodal point " << np + 1 << endl;
			return false;
        }
    
        if (NodeList[np].NodeNumber != np + 1)
        {
            cerr << "*** Error *** Nodes must be inputted in order ! (At node " << np + 1 << ")" << endl
            << "   Expected node number : " << np + 1 << endl
            << "   Provided node number : " << NodeList[np].NodeNumber << endl;
            return false;
        }
    }
    cout << "DEBUG: Successfully read all " << NUMNP << " nodes" << endl;
	return true;
}

//	Calculate global equation numbers corresponding to every degree of freedom of each node
void CDomain::CalculateEquationNumber()
{
	NEQ = 0;
	for (unsigned int np = 0; np < NUMNP; np++)	// Loop over for all node
	{
		for (unsigned int dof = 0; dof < CNode::NDF; dof++)	// Loop over for DOFs of node np
		{
			if (NodeList[np].bcode[dof]) 
				NodeList[np].bcode[dof] = 0;
			else
			{
				NEQ++;
				NodeList[np].bcode[dof] = NEQ;
			}
		}
	}
}

//	Read load case data
bool CDomain::ReadLoadCases()
{
    if (NLCASE == 0) {
        // cerr << "Info: NLCASE is 0, skipping load case reading." << endl;
        return true;
    }
//	Read load data lines
	LoadCases = new CLoadCaseData[NLCASE];	// List all load cases
    if (!LoadCases) {
        cerr << "*** Error *** Memory allocation failed for LoadCases (NLCASE = " << NLCASE << ")" << endl;
        return false;
    }

//	Loop over for all load cases
	for (unsigned int lcase = 0; lcase < NLCASE; lcase++)
    {
        unsigned int LL;
        // --- 新增调试和预处理代码 ---
        cout << "DEBUG: CDomain::ReadLoadCases - Reading Load Case ID for iteration " << lcase + 1 << endl;
        cout << "DEBUG: CDomain::ReadLoadCases - Stream state before reading LL: failbit=" << Input.fail()
             << ", badbit=" << Input.bad() << ", eofbit=" << Input.eof() << endl;
        Input >> ws; // 跳过前导的空白字符（非常重要！）
        char peek_char_LL = '#'; // Default
        if(Input.good()) peek_char_LL = Input.peek();
        cout << "DEBUG: CDomain::ReadLoadCases - Peek next char AFTER ws before reading LL: '" << peek_char_LL << "' (ASCII: " << (int)peek_char_LL << ")" << endl;
        // --- 结束新增代码 ---
        Input >> LL;
        // --- 新增调试代码 ---
        cout << "DEBUG: CDomain::ReadLoadCases - Attempted to read LL: " << LL << endl;
        cout << "DEBUG: CDomain::ReadLoadCases - Stream state after reading LL: failbit=" << Input.fail()
             << ", badbit=" << Input.bad() << ", eofbit=" << Input.eof() << endl;
        // --- 结束新增代码 ---

        if (Input.fail()){
            cerr << "*** Error *** Failed to read load case number for load case " << lcase + 1 << endl;
            // 新增更详细的错误输出
            Input.clear(); 
            char error_snippet_buffer_LL[100] = {0}; 
            streampos error_pos_LL = Input.tellg();
            Input.read(error_snippet_buffer_LL, 99);
            Input.clear(); 
            Input.seekg(error_pos_LL); 
            cerr << "       CDomain::ReadLoadCases - Failed when trying to interpret data. Input near error: [\"" << error_snippet_buffer_LL << "\"] as LL." << endl;
            // printf("LL: %d\n", LL); // LL的值此时是不可靠的
            return false;
        }
        
        if (LL != lcase + 1)
        {
            cerr << "*** Error *** Load case must be inputted in order ! (At load case " << lcase + 1 << ")" << endl
            << "   Expected load case : " << lcase + 1 << endl
            << "   Provided load case : " << LL << endl;
            
            return false;
        }
        // --- 新增调试代码 ---
        cout << "DEBUG: CDomain::ReadLoadCases - Successfully read LL = " << LL << ". Calling LoadCases[" << lcase << "].Read(Input)." << endl;
        // --- 结束新增代码 ---

        if (!LoadCases[lcase].Read(Input)) // Check return value of Read
        {
            cerr << "*** Error *** Failed to read data for load case " << LL << endl;
            return false;
        }
    }
    cout << "Successfully read " << NLCASE << " load cases." << endl; // Success message
	return true;
}

// Read element data
bool CDomain::ReadElements()
{
    if (NUMEG == 0 && NUMNP > 0) { // If nodes exist but no element groups, may need a warning
        cerr << "*** Warning *** NUMEG (Number of Element Groups) is 0, but nodes exist. No elements will be defined." << endl;
        return true; // May be a valid case
    }
    if (NUMEG == 0) {
        return true; // No element groups to read
    }

    EleGrpList = new CElementGroup[NUMEG];
    if (!EleGrpList) {
        cerr << "*** Error *** Memory allocation failed for EleGrpList (NUMEG = " << NUMEG << ")" << endl;
        return false;
    }

//	Loop over for all element group
	for (unsigned int eg = 0; eg < NUMEG; eg++) // Use different loop variable name 'eg'
    {
        cout << "Reading Element Group " << eg + 1 << "..." << endl;
        if (!EleGrpList[eg].Read(Input))
        {
            cerr << "*** Error *** Failed to read data for element group " << eg + 1 << endl;
            return false;
        }
        if (EleGrpList[eg].GetNUME() == 0 && EleGrpList[eg].GetElementType() != ElementTypes::UNDEFINED) { // If type is defined but no elements
            cerr << "*** Warning *** Element group " << eg + 1 << " (Type: " << EleGrpList[eg].GetElementType() 
                 << ") has 0 elements defined. Check input file if this is unintended." << endl;
        }
    }
    cout << "Successfully read " << NUMEG << " element groups." << endl; // Success message
    return true;
}

//	Calculate column heights
void CDomain::CalculateColumnHeights()
{
#ifdef _DEBUG_
    COutputter* Output = COutputter::GetInstance();
    *Output << setw(9) << "Ele = " << setw(22) << "Location Matrix" << endl;
#endif

	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)		//	Loop over for all element groups
    {
        CElementGroup& ElementGrp = EleGrpList[EleGrp];
        unsigned int NUME = ElementGrp.GetNUME();
        
		for (unsigned int Ele = 0; Ele < NUME; Ele++)	//	Loop over for all elements in group EleGrp
        {
            CElement& Element = ElementGrp[Ele];

            // Generate location matrix
            Element.GenerateLocationMatrix();
            
#ifdef _DEBUG_
            unsigned int* LocationMatrix = Element.GetLocationMatrix();
            
            *Output << setw(9) << Ele+1;
            for (int i=0; i<Element.GetND(); i++)
                *Output << setw(5) << LocationMatrix[i];
            *Output << endl;
#endif

            StiffnessMatrix->CalculateColumnHeight(Element.GetLocationMatrix(), Element.GetND());
        }
    }
    
    StiffnessMatrix->CalculateMaximumHalfBandwidth();
    
#ifdef _DEBUG_
    *Output << endl;
	Output->PrintColumnHeights();
#endif

}

/*!	Allocate Force, ColumnHeights, DiagonalAddress and StiffnessMatrix and 
    calculate the column heights and address of diagonal elements */
void CDomain::AllocateMatrices()
{
	StiffnessMatrix = new CSkylineMatrix<double>(NEQ);
	Force = new double[NEQ];

	CalculateColumnHeights();

	StiffnessMatrix->CalculateDiagnoalAddress();
	StiffnessMatrix->Allocate();
}

//	Assemble the banded gloabl stiffness matrix
void CDomain::AssembleStiffnessMatrix()
{
	COutputter* Output = COutputter::GetInstance();
    *Output << "\n*** DEBUG: Entering AssembleStiffnessMatrix... ***\n";

//	Loop over for all element groups
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
		*Output << "*** DEBUG: Processing Element Group " << EleGrp + 1 << "/" << NUMEG << " ***\n";
		CElementGroup& ElementGroup = EleGrpList[EleGrp];
		unsigned int NUME = ElementGroup.GetNUME();
		*Output << "*** DEBUG: Number of elements in this group: " << NUME << " ***\n";

//		Loop over for all elements in group EleGrp
		for (unsigned int Ele = 0; Ele < NUME; Ele++)
		{
			*Output << "*** DEBUG: Processing Element " << Ele + 1 << "/" << NUME << " in group " << EleGrp + 1 << " ***\n";
			// Get the element
			CElement& Element = ElementGroup[Ele];

			// Get the size of the element stiffness matrix
			unsigned int size = Element.SizeOfStiffnessMatrix();
			*Output << "*** DEBUG: Element stiffness matrix size: " << size << " ***\n";

			// Allocate memory for the element stiffness matrix
			double* kes = new double[size];
			*Output << "*** DEBUG: Allocated temporary kes array. ***\n";

			// Calculate element stiffness matrix
			*Output << "*** DEBUG: Calling ElementStiffness for element " << Ele + 1 << "... ***\n";
			Element.ElementStiffness(kes);
			*Output << "*** DEBUG: ElementStiffness for element " << Ele + 1 << " returned. ***\n";

			// Assemble element stiffness matrix to global stiffness matrix
			*Output << "*** DEBUG: Calling Assembly for element " << Ele + 1 << "... ***\n";
			StiffnessMatrix->Assembly(kes, Element.GetLocationMatrix(), Element.GetND());
			*Output << "*** DEBUG: Assembly for element " << Ele + 1 << " returned. ***\n";

			// Deallocate memory
			delete[] kes;
			*Output << "*** DEBUG: Deallocated temporary kes array for element " << Ele + 1 << ". ***\n";
		}
	}

#ifdef _DEBUG_
	*Output << "\n\nS T I F F N E S S   M A T R I X\n\n";
	StiffnessMatrix->Print(*Output);
#endif
	*Output << "\n*** DEBUG: Exiting AssembleStiffnessMatrix. ***\n";
}

//	Assemble the global nodal force vector for load case LoadCase
bool CDomain::AssembleForce(unsigned int LoadCase)
{
	if (LoadCase > NLCASE) 
		return false;

	CLoadCaseData* LoadData = &LoadCases[LoadCase - 1];

    clear(Force, NEQ);

//	Loop over for all concentrated loads in load case LoadCase
	for (unsigned int lnum = 0; lnum < LoadData->nloads; lnum++)
	{
		unsigned int dof = NodeList[LoadData->node[lnum] - 1].bcode[LoadData->dof[lnum] - 1];
        
        if(dof) // The DOF is activated
            Force[dof - 1] += LoadData->load[lnum];
	}

	return true;
}

