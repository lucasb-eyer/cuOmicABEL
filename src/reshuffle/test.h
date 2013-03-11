/*
 * test.h
 *
 *  Created on: 07.03.2013
 *      Author: lima
 */

#ifndef TEST_H_
#define TEST_H_
#include "parameters.h"
#include "iout_file.h"
#include "reshuffle.h"
//#include <dir.h>
using namespace std;

class test{
public:
	string name;
	string result;
	string check;
	Parameters Params;

	Parameters *Params_test=&Params;
	test(string,string,string,string);
	void run();
};
#endif /* TEST_H_ */
