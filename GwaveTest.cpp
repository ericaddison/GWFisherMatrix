// GwaveTest.cpp -- a test function for the Gwave namespace


#include "fisher.h"
#include <iostream>
#include "GravWave.h"
#include "gslWrappers.h"

int main()
{
	using std::cout;
	using GravWave::GW_Test;
	using GravWave::Ftest;
	using gslWrappers::MatrixTest;
	using gslWrappers::IntTest;
	using gslWrappers::rootTest;

	cout << "Running GW_Test()...\n\n";
	GW_Test();

	cout << "Running Ftest()...\n\n";
	Ftest();

	//cout << "Running MatrixTest()...\n\n";
	//MatrixTest();	

//	cout << "Running IntTest()...\n\n";
//	IntTest();	
	
	//cout << "Running rootTest()...\n\n";
	//rootTest();
	return 0;
}
