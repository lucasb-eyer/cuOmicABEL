/*
 * test.cpp
 *
 *  Created on: 07.03.2013
 *      Author: lima
 */
#include "test.h"
#include <string.h>
#include <stdlib.h>

using namespace std;

ofstream test_txt("test.txt");
test::test(string Name,string cmd,string Result,string Check){

	name=Name;
	result=Result;
	check=Check;
	Params_test->iout_fname="data_4test/4test.iout";
	Params_test->out_fname="data_4test/4test.out";
	if(name=="datadims")
		Params_test->datadims = Parameter(cmd,"datadims");
	if(name=="traitnamesdef"||name=="traitnames")
		Params_test->traitnames = Parameter(cmd,"traitnames");
	if(name=="snpnamesdef"||name=="snpnames")
		Params_test->snpnames = Parameter(cmd,"snpnames");
	if(name=="traitsnumbers"||name=="traitsdef"||name=="traitsbynames")
		Params_test->traits = Parameter(cmd,"traits");
	if(name=="snpsnumbers"||name=="snpsdef"||name=="snpsbynames")
		Params_test->snps = Parameter(cmd,"snps");
	if(name=="traitssnpscombo"){
		Params_test->traits=Parameter(cmd,"traits");
		Params_test->snps=Parameter(cmd,"snps");
	}
	if(name=="chiall"||name=="chiover")
		Params_test->chi = Parameter(cmd,"chi");
	if(name=="traitssnpschicombo"){
		Params_test->traits=Parameter(cmd,"traits");
		Params_test->snps=Parameter(cmd,"snps");
		Params_test->chi = Parameter(cmd,"chi");
	}
	if(name=="estdef"||name=="estnumbers"||name=="estnames")
		Params_test->heritabilities = Parameter(cmd,"heritabilities");
	if(name=="dataslim"){
		Params_test->dataslim = Parameter(cmd,"dataslim");
		Params_test->chi = Parameter(cmd,"chi");
	}

}

void test::run(){

	test_txt<<"START TEST "<<name<<"\t";
	iout_file iout_F(*Params_test);
	Reshuffle reshh(iout_F,*Params_test);
	reshh.run();
	//result="datadims/datadims.txt";
	//check="data_4test_check/datadims/datadims.txt";
	ifstream result_f(result.c_str());
	ifstream check_f(check.c_str());
	string str_res="";
	string str_che="";
	int checker=0;
	while (getline(result_f,str_res)){
		getline(check_f,str_che);
		if(strcmp(str_che.c_str(),str_res.c_str()))
			checker++;
	}
	if(checker!=0){
		test_txt<<"Test "<<name<<" FAILED!!!"<<endl;
	}else
		test_txt<<"Test "<<name<<" OK"<<endl;

	//TODO: remove doesn't work from test.cpp
	/*if(remove(result.c_str()) != 0 )
	  cout<<"Error deleting file";*/

}
