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
	Params_test= new Parameters;
	name=Name;
	result=Result;
	check=Check;
	Params_test->iout_fname="data_4test/4test.iout";
	Params_test->out_fname="data_4test/4test.out";
	if(name=="datadims")
		Params_test->datadims = Parameter(cmd,"datadims","datadims.txt");
	if(name=="traitnamesdef"||name=="traitnames")
		Params_test->traitnames = Parameter(cmd,"traitnames","traitnames.txt");
	if(name=="snpnamesdef"||name=="snpnames","snpnames.txt")
		Params_test->snpnames = Parameter(cmd,"snpnames","snpnames.txt");
	if(name=="traitsnumbers"||name=="traitsdef"||name=="traitsbynames")
		Params_test->traits = Parameter(cmd,"traits","data.txt");
	if(name=="snpsnumbers"||name=="snpsdef"||name=="snpsbynames")
		Params_test->snps = Parameter(cmd,"snps","data.txt");
	if(name=="traitssnpscombo"){
		Params_test->traits=Parameter(cmd,"traits","data.txt");
		Params_test->snps=Parameter(cmd,"snps","data.txt");
	}
	if(name=="chiall"||name=="chiover")
		Params_test->chi = Parameter(cmd,"chi","chi_data.txt");
	if(name=="traitssnpschicombo"){
		Params_test->traits=Parameter(cmd,"traits","data.txt");
		Params_test->snps=Parameter(cmd,"snps","data.txt");
		Params_test->chi = Parameter(cmd,"chi","chi_data.txt");
	}
	if(name=="estdef"||name=="estnumbers"||name=="estnames")
		Params_test->heritabilities = Parameter(cmd,"heritabilities","estimates.txt");
	if(name=="dataslim"){
		Params_test->dataslim = Parameter(cmd,"dataslim","slim_data.txt");
		Params_test->chi = Parameter(cmd,"chi","chi_data.txt");
	}
}

void test::run(){
	test_txt<<"START TEST "<<name<<"\t";
	iout_file iout_F(*Params_test);
	if(Params_test->traits.use)
		Params_test->traits.setbynames(*(iout_F.labels.trait_names));
	if(Params_test->snps.use)
		Params_test->snps.setbynames(*(iout_F.labels.snp_names));
	if(Params_test->heritabilities.use)
		Params_test->heritabilities.setbynames(*(iout_F.labels.trait_names));
	Reshuffle reshh(iout_F,*Params_test);
	reshh.run();
	ifstream result_f(result.c_str());
	ifstream check_f(check.c_str());
	string str_res="";
	string str_che="";
	int checker=0;
	while (getline(check_f,str_che)){
		getline(result_f,str_res);
		if(strcmp(str_che.c_str(),str_res.c_str())!=0)
			checker++;
	}
	if(getline(result_f,str_res))
		checker++;
	if(checker!=0){
		test_txt<<"Test "<<name<<" FAILED!!!"<<endl;
	}else
		test_txt<<"Test "<<name<<" OK"<<endl;

	//TODO: remove doesn't work from test.cpp
	/*if(remove(result.c_strf()) != 0 )
	  cout<<"Error deleting file";*/

}
