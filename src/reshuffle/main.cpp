/*
 * main.cpp
 *
 *  Created on: 31.01.2013
 *      Author: Sodbo
 */
#include <iostream>
#include <ctime>
#include "iout_file.h"
#include "Parameters.h"
#include "reshuffle.h"
#include <iterator>
#include "test.h"

using namespace std;

int main(int argc, char* argv[]) {
	//cout<<"Start reshuffeling"<<endl;
	Parameters Params(argc, argv);
	if(Params.h.use||Params.help.use||(Params.iout_fname==".iout"&&Params.out_fname==".out")){
		cout<<"Available commands"<<endl<<endl;
		cout<<" --datadims | to get data's dimension"<<endl<<endl;
		cout<<" --snpnames=<indexes> | to get names of SNP by indexes"<<endl<<endl;
		cout<<" --traitnames=<indexes> | to get names of trait by indexes"<<endl<<endl;
		cout<<" --traits=<indexes OR/AND names OR/AND regexp> | to get data"<<endl;
		cout<<"   by trait's indexes OR/AND names OR/AND regexp"<<endl<<endl;
		cout<<" --snp=<indexes OR/AND names OR/AND focus> | to get data"<<endl;
		cout<<"   by SNP's indexes OR/AND names"<<endl<<endl;
		cout<<" --heritabilities=<indexes OR/AND names OR/AND regexp>"<<endl;
		cout<<"   to get estimates of trait's heritability, sigma, res_sigma and betas"<<endl<<endl;
		cout<<" --chi=<number> | to get data by snp's, which chi2>number"<<endl<<endl;
		cout<<" --dataslim | to get slim data [You should set --chi=<number>]"<<endl<<endl;
		exit(1);
	}

	if(Params.test.use){
		test t_datadims("datadims","--datadims--","datadims.txt","data_4test_check/datadims.txt");
		t_datadims.run();
		remove(t_datadims.result.c_str());
		test t_snpnamesdef("snpnamesdef","--snpnames--","snpnames.txt","data_4test_check/snpnamesdef.txt");
		t_snpnamesdef.run();
		remove(t_snpnamesdef.result.c_str());

		test t_traitnamesdef("traitnamesdef","--traitnames--","traitnames.txt","data_4test_check/traitnamesdef.txt");
		t_traitnamesdef.run();
		remove(t_traitnamesdef.result.c_str());

		test t_snpnames("snpnames","--snpnames=10-27,1,20-35,55--","snpnames.txt","data_4test_check/snpnames.txt");
		t_snpnames.run();
		remove(t_snpnames.result.c_str());

		test t_traitnames("traitnames","--traitnames=1-3,5--","traitnames.txt","data_4test_check/traitnames.txt");
		t_traitnames.run();
		remove(t_traitnames.result.c_str());

		test t_snpsdef("snpsdef","--snps--","data.txt","data_4test_check/data_def.txt");
		t_snpsdef.run();
		remove(t_snpsdef.result.c_str());

		test t_traitsdef("traitsdef","--traits--","data.txt","data_4test_check/data_def.txt");
		t_traitsdef.run();
		remove(t_traitsdef.result.c_str());

		test t_traitsnumbers("traitsnumbers","--traits=1-3,5--","data.txt","data_4test_check/data_traitnumbers.txt");
		t_traitsnumbers.run();
		remove(t_traitsnumbers.result.c_str());

		test t_traitsbynames("traitsbynames","--traits=1-2,tca--","data.txt","data_4test_check/data_traitsbynames.txt");
		t_traitsbynames.run();
		remove(t_traitsbynames.result.c_str());

		test t_snpsnumbers("snpsnumbers","--snps=1-20,66,15-30,13--","data.txt","data_4test_check/data_snpsnumbers.txt");
		t_snpsnumbers.run();
		remove(t_snpsnumbers.result.c_str());

		test t_snpsbynames("snpsbynames","--snps=3-8,rs3121561,10,rs6687776--","data.txt","data_4test_check/data_snpsbynames.txt");
		t_snpsbynames.run();
		remove(t_snpsbynames.result.c_str());

		test t_traitssnpscombo("traitssnpscombo","--snps=1-20,66,15-30,13--traits=1-3,5--","data.txt","data_4test_check/data_ts_combo.txt");
		t_traitssnpscombo.run();
		remove(t_traitssnpscombo.result.c_str());

		test t_chiall("chiall","--chi--","chi_data.txt","data_4test_check/data_chiall.txt");
		t_chiall.run();
		remove(t_chiall.result.c_str());

		test t_chiover("chiover","--chi=2--","chi_data.txt","data_4test_check/data_chiover.txt");
		t_chiover.run();
		remove(t_chiover.result.c_str());

		test t_traitssnpschicombo("traitssnpschicombo","--chi=2--traits--1-3,5--snps=1-20,66,15-30,13--","chi_data.txt","data_4test_check/data_chi_trsnp.txt");
		t_traitssnpschicombo.run();
		remove(t_traitssnpschicombo.result.c_str());

		test t_estdef("estdef","--heritabilities--","estimates.txt","data_4test_check/estdef.txt");
		t_estdef.run();
		remove(t_estdef.result.c_str());

		test t_estnum("estnumbers","--heritabilities=1-3,5--","estimates.txt","data_4test_check/estnum.txt");
		t_estnum.run();
		remove(t_estnum.result.c_str());

		test t_estnames("estnames","--heritabilities=1-2,tca,hdla--","estimates.txt","data_4test_check/estnames.txt");
		t_estnames.run();
		remove(t_estnames.result.c_str());

		test t_dataslim("dataslim","--chi=5--dataslim--","slim_data.txt","data_4test_check/data_slim.txt");
		t_dataslim.run();
		remove(t_dataslim.result.c_str());
	}
	if(Params.info.use)
		cout << Params;
	iout_file iout_F(Params);
	cout << "Finish iout_file read\t" << double(clock()) / CLOCKS_PER_SEC <<" sec" << endl;
	if(Params.info.use){
		cout<<iout_F.header;
		cout<<iout_F.labels;
	}
	if(Params.traits.use)
		Params.traits.setbynames(*(iout_F.labels.trait_names));
	if(Params.snps.use)
		Params.snps.setbynames(*(iout_F.labels.snp_names));
	if(Params.heritabilities.use)
		Params.heritabilities.setbynames(*(iout_F.labels.trait_names));
	Reshuffle reshh(iout_F,Params);
	reshh.run();
	cout << "Finish reshuffling " << double(clock()) / CLOCKS_PER_SEC <<" sec" << endl;
	return (0);
}
