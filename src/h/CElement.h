//!	Calculate element stress 
	virtual void ElementStress(double* stress_array, double* Displacement) = 0; 

//! Return the number of double values this element will write into stress_array
    virtual unsigned int GetNumberOfStressStrainOutputs() const 
    { 
        return 1; // 默认值为1，以兼容CBar单元 (输出一个轴力/应力)
    }

//! Return number of nodes per element 