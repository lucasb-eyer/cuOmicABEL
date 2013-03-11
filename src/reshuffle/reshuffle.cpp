/*
 * reshuffle.cpp
 *
 *  Created on: 01.03.2013
 *      Author: lima
 */

#include "reshuffle.h"
//#include <dir.h>
#include <ctime>
#include <math.h>
#include <iterator>
#include <list>

using namespace std;

#define PRECISION_DOUBLE 15//Precision of double

Reshuffle::Reshuffle(iout_file &iout,Parameters &Params){
	p_iout_file = &iout;
	p_Parameters = &Params;
	per_trait_per_snp = (*p_iout_file).header.p + (*p_iout_file).header.p * ((*p_iout_file).header.p + 1) / 2;
	herest_startpos = (p_iout_file->header.m * p_iout_file->header.t*
			(p_iout_file->header.p + p_iout_file->header.p * (p_iout_file->header.p + 1) / 2))
					* sizeof(double);
}

string Reshuffle::create_filename(string prename, string name) {
	string path = prename + "_"+name+".txt";
	return path;
}
string Reshuffle::create_filename(string name) {
	string path = name+".txt";
	return path;
}

void Reshuffle::write_datadims(ofstream& txt_datadims){

	txt_datadims << "Number of traits\t" << (*p_iout_file).header.t << endl;
	txt_datadims << "Number of SNP\t" << (*p_iout_file).header.m << endl;
	txt_datadims << "Number of covariates\t" << ((*p_iout_file).header.p - 2);
}

void Reshuffle::write_snpnames(ofstream& txt_snpnames){

	if ((*p_Parameters).snpnames.value == "all"){
		for (int i=0;i<(*(*p_iout_file).labels.snp_names).size();i++)
			(*p_Parameters).snpnames.valueset.insert(i);
		cout<<"SNPNAMES VALUE SET CHANGED TO ALL"<<endl;
	}
	for(set<int>::iterator it= (*p_Parameters).snpnames.valueset.begin();it!=(*p_Parameters).snpnames.valueset.end();++it)
		txt_snpnames << "SNP #"<<(*it+1)<<"\t"<<(*(*p_iout_file).labels.snp_names)[*it]<<endl;
	cout<<"END WRITE SNPNAMES"<<endl;
}

void Reshuffle::write_traitnames(ofstream& txt_traitnames){

	if ((*p_Parameters).traitnames.value == "all"){
		for (int i=0;i<(*(*p_iout_file).labels.trait_names).size();i++)
			(*p_Parameters).traitnames.valueset.insert(i);
		cout<<"TRAITNAMES VALUE SET CHANGED TO ALL"<<endl;
	}
	for(std::set<int>::iterator it= (*p_Parameters).traitnames.valueset.begin();it!=p_Parameters->traitnames.valueset.end();++it)
		txt_traitnames<<"TRAIT #"<<(*it+1)<<"\t"<<(*(*p_iout_file).labels.trait_names)[*it]<<endl;
	cout<<"END WRITE TRAITNAMES"<<endl;
}

void Reshuffle::write_data(ifstream& out_file){
	out_file.seekg(0, ios_base::beg);
	cout << "startwritetxt=" << double(clock()) / CLOCKS_PER_SEC << endl;
	ofstream txt_trait(create_filename("data").c_str());
	for (set<int>::iterator trait= (*p_Parameters).traits.valueset.begin();trait!=(*p_Parameters).traits.valueset.end();trait++) {

		//Set precision of double
		cout<<(*(*p_iout_file).labels.trait_names)[*trait]<<endl;
		txt_trait.precision(PRECISION_DOUBLE);
		txt_trait<<	(*(*p_iout_file).labels.trait_names)[*trait]<<endl;
		for (int beta = 0;	beta < (*(*p_iout_file).labels.beta).size(); beta++)
			txt_trait << (*(*p_iout_file).labels.beta)[beta] << "\t";
		for (int se = 0;se < (*(*p_iout_file).labels.se).size(); se++)
			txt_trait << (*(*p_iout_file).labels.se)[se] << "\t";
		for (int cov = 0;cov < (*(*p_iout_file).labels.cov).size(); cov++)
			txt_trait << (*(*p_iout_file).labels.cov)[cov] << "\t";
		txt_trait << endl;
		cout << "endwritetrait_colnames" << *trait << clock() << endl;
		double* buf = new double[per_trait_per_snp];
		int oldPos = 0;
		char s[30];
		for (set<int>::iterator snp= (*p_Parameters).snps.valueset.begin();snp!=(*p_Parameters).snps.valueset.end();snp++) {
			txt_trait << (*(*p_iout_file).labels.snp_names)[*snp] << "\t";
			int pos = (*p_iout_file).tilecoordinates(*trait, *snp);
			//cout << oldPos << "-" << pos << endl;
			if(pos != oldPos)
			{
				out_file.seekg(pos,ios_base::beg);
			}
			oldPos=pos+sizeof(double)*per_trait_per_snp;
			out_file.read((char *)buf, sizeof(double)*per_trait_per_snp);
			for (int i = 0; i < per_trait_per_snp; i++) {
				sprintf(s, "%.15g", buf[i]);
				txt_trait << s << "\t";
			}
			txt_trait << endl;
		}
		delete buf;
		//txt_trait.close();
		cout << "endwritetrait " << (*(*p_iout_file).labels.trait_names)[*trait] << " "<< double(clock()) / CLOCKS_PER_SEC << endl;
	}
	cout << "finishwritetxt=" << double(clock()) / CLOCKS_PER_SEC << endl;
}

void Reshuffle::write_data_chi(ifstream& out_file){
	out_file.seekg(0, ios_base::beg);
	double chi = 0;
	double CheckChi = (*p_Parameters).chi.value == "all" ? -1.0 : atof((*p_Parameters).chi.value.c_str());
	cout << "startwritetxt=" << double(clock()) / CLOCKS_PER_SEC << endl;
	ofstream txt_chi(create_filename("chi_data").c_str());
	for (set<int>::iterator trait= (*p_Parameters).traits.valueset.begin();trait!=(*p_Parameters).traits.valueset.end();trait++) {
		//ofstream txt_chi(create_filename("chi_data//chi", (*(*p_iout_file).labels.trait_names)[*trait]).c_str());
		//Set precision of double
		txt_chi.precision(PRECISION_DOUBLE);
		txt_chi<<	(*(*p_iout_file).labels.trait_names)[*trait]<<endl;
		for (int beta = 0;	beta < (*(*p_iout_file).labels.beta).size(); beta++)
			txt_chi << (*(*p_iout_file).labels.beta)[beta] << "\t";
		for (int se = 0;se < (*(*p_iout_file).labels.se).size(); se++)
			txt_chi << (*(*p_iout_file).labels.se)[se] << "\t";
		for (int cov = 0;cov < (*(*p_iout_file).labels.cov).size(); cov++)
			txt_chi << (*(*p_iout_file).labels.cov)[cov] << "\t";
		txt_chi << "Chi2" << endl;
		double* buf = new double[per_trait_per_snp];
		int oldPos = 0;
		char s[30];
		for (set<int>::iterator snp= (*p_Parameters).snps.valueset.begin();snp!=(*p_Parameters).snps.valueset.end();snp++) {
			int pos = (*p_iout_file).tilecoordinates(*trait, *snp);
			//cout << oldPos << "-" << pos << endl;
			if(pos != oldPos)
			{
				out_file.seekg(pos,ios_base::beg);
			}
			oldPos=pos+sizeof(double)*per_trait_per_snp;
			out_file.read((char *)buf, sizeof(double)*per_trait_per_snp);
			chi=pow((buf[(*(*p_iout_file).labels.beta).size()-1]/buf[(*(*p_iout_file).labels.beta).size()+(*(*p_iout_file).labels.se).size()-1]),2);
			if(chi>CheckChi){
				txt_chi << (*(*p_iout_file).labels.snp_names)[*snp] << "\t";
				for (int i = 0; i < per_trait_per_snp; i++) {
					sprintf(s, "%.15g", buf[i]);
					txt_chi << s << "\t";
				}
				txt_chi << chi << endl;
			}
		}
		delete buf;
		//txt_chi.close();
		cout << "endwritechitrait " << (*(*p_iout_file).labels.trait_names)[*trait] << " "<< double(clock()) / CLOCKS_PER_SEC << endl;
	}
	cout << "finishwritechitxt=" << double(clock()) / CLOCKS_PER_SEC << endl;
}

void Reshuffle::write_slim_data(ifstream& out_file){
	out_file.seekg(0, ios_base::beg);
	set<int> goodtraits;
	set<int> goodsnps;
	double chi = 0;
	if((*p_Parameters).chi.value=="all"||(*p_Parameters).chi.value=="None"){
		cout << "ERROR" << "Chi value doesn't set"<<endl;
		cout << "Please, set Chi value to get slim data" << endl;
		exit(1);
	}
	double CheckChi = atof((*p_Parameters).chi.value.c_str());
	for (set<int>::iterator trait= (*p_Parameters).traits.valueset.begin();trait!=(*p_Parameters).traits.valueset.end();trait++) {
		double* buf = new double[per_trait_per_snp];
		int oldPos = 0;
		char s[30];
		for (set<int>::iterator snp= (*p_Parameters).snps.valueset.begin();snp!=(*p_Parameters).snps.valueset.end();snp++) {
			int pos = (*p_iout_file).tilecoordinates(*trait, *snp);
			//cout << oldPos << "-" << pos << endl;
			if(pos != oldPos)
			{
				out_file.seekg(pos,ios_base::beg);
			}
			oldPos=pos+sizeof(double)*per_trait_per_snp;
			out_file.read((char *)buf, sizeof(double)*per_trait_per_snp);
			chi=pow((buf[(*(*p_iout_file).labels.beta).size()-1]/buf[(*(*p_iout_file).labels.beta).size()+(*(*p_iout_file).labels.se).size()-1]),2);
			if(chi>CheckChi){
				goodtraits.insert(*trait);
				goodsnps.insert(*snp);
			}
		}
		delete buf;
	}
	out_file.seekg(0, ios_base::beg);
	ofstream txt_slim(create_filename("slim_data").c_str());

	for (set<int>::iterator trait= goodtraits.begin();trait!=goodtraits.end();trait++) {
		//Set precision of double
		txt_slim.precision(PRECISION_DOUBLE);
		txt_slim<<	(*(*p_iout_file).labels.trait_names)[*trait]<<endl;
		for (int beta = 0;	beta < (*(*p_iout_file).labels.beta).size(); beta++)
			txt_slim << (*(*p_iout_file).labels.beta)[beta] << "\t";
		for (int se = 0;se < (*(*p_iout_file).labels.se).size(); se++)
			txt_slim << (*(*p_iout_file).labels.se)[se] << "\t";
		for (int cov = 0;cov < (*(*p_iout_file).labels.cov).size(); cov++)
			txt_slim << (*(*p_iout_file).labels.cov)[cov] << "\t";
		txt_slim << "Chi2" << endl;
		double* buf = new double[per_trait_per_snp];
		int oldPos = 0;
		char s[30];
		for (set<int>::iterator snp= goodsnps.begin();snp!=goodsnps.end();snp++) {
			txt_slim << (*(*p_iout_file).labels.snp_names)[*snp] << "\t";
			int pos = (*p_iout_file).tilecoordinates(*trait, *snp);
			//cout << oldPos << "-" << pos << endl;
			if(pos != oldPos)
			{
				out_file.seekg(pos,ios_base::beg);
			}
			oldPos=pos+sizeof(double)*per_trait_per_snp;
			out_file.read((char *)buf, sizeof(double)*per_trait_per_snp);
			chi=pow((buf[(*(*p_iout_file).labels.beta).size()-1]/buf[(*(*p_iout_file).labels.beta).size()+(*(*p_iout_file).labels.se).size()-1]),2);
			for (int i = 0; i < per_trait_per_snp; i++) {
				sprintf(s, "%.15g", buf[i]);
				txt_slim << s << "\t";
			}
			txt_slim << chi << endl;
		}
		delete buf;
		//txt_slim.close();
	}
}

int Reshuffle::est_shift(int counter){
	int shift = ((p_iout_file->header.m * p_iout_file->header.t* (p_iout_file->header.p +
			p_iout_file->header.p * (p_iout_file->header.p + 1) / 2))
			+counter*p_iout_file->header.t)* sizeof(double);
	return shift;
}

int Reshuffle::est_beta_shift(int counter){
	int shift = est_shift(3)+counter*(p_iout_file->header.p-1)*sizeof(double);
	return shift;
}

void Reshuffle::write_herest(ifstream& out_file){
	ofstream txt_est("estimates.txt");
	out_file.seekg(herest_startpos, ios_base::beg);
	if (p_Parameters->heritabilities.value == "all")
		for(int i=0;i<(*(p_iout_file->labels.trait_names)).size();i++)
			p_Parameters->heritabilities.valueset.insert(i);
	txt_est.precision(PRECISION_DOUBLE);
	txt_est<<"\t";
	for (set<int>::iterator trait= p_Parameters->heritabilities.valueset.begin();trait!=p_Parameters->heritabilities.valueset.end();trait++)
		txt_est << (*(p_iout_file->labels.trait_names))[*trait] << "\t";
	txt_est << endl;
	list<string> est_names;
	est_names.insert(est_names.end(), "heritabilities");
	est_names.insert(est_names.end(), "sigma");
	est_names.insert(est_names.end(), "res_sigma");
	double tmp_number = 0;
	int counter=0;
	for (list<string>::iterator name = est_names.begin();name != est_names.end(); ++name) {
		txt_est << *name << "\t";
		for (std::set<int>::iterator trait= p_Parameters->heritabilities.valueset.begin();trait!=p_Parameters->heritabilities.valueset.end();++trait) {
			out_file.seekg(*trait*sizeof(double),ios_base::cur);
			out_file.read((char *) &tmp_number, sizeof(double));
			txt_est << tmp_number << "\t";
			out_file.seekg(est_shift(counter), ios_base::beg);
		}
		counter++;
		txt_est << "\n";
		out_file.seekg(est_shift(counter), ios_base::beg);
	}
	out_file.seekg(est_shift(3), ios_base::beg);
	counter=0;
	for (int beta=0;beta<(*(p_iout_file->labels.beta)).size();beta++) {
		beta++;
		if (beta != (*(p_iout_file->labels.beta)).size()) {
			beta--;
			txt_est << (*(p_iout_file->labels).beta)[beta] << "\t";
			for (std::set<int>::iterator trait= p_Parameters->heritabilities.valueset.begin();trait!=p_Parameters->heritabilities.valueset.end();++trait) {
				out_file.seekg(*trait*sizeof(double),ios_base::cur);
				out_file.read((char *) &tmp_number, sizeof(double));
				txt_est << tmp_number << "\t";
				out_file.seekg(est_beta_shift(counter), ios_base::beg);
			}
			counter++;
			out_file.seekg(est_beta_shift(counter), ios_base::beg);
			txt_est << "\n";
			beta++;
		}
	}
}

void Reshuffle::run(){
	if((*p_Parameters).datadims.use){
		//mkdir("datadims");
		ofstream datadims("datadims.txt");
		write_datadims(datadims);
	}
	if((*p_Parameters).snpnames.use){
		//mkdir("snpnames");
		ofstream snpnames("snpnames.txt");
		write_snpnames(snpnames);
	}

	if((*p_Parameters).traitnames.use){
		//mkdir("traitnames");
		ofstream* traitnames = new ofstream("traitnames.txt");
		write_traitnames(*traitnames);
		delete traitnames;
	}

	//Open *.out to read data
	ifstream out_file((*p_Parameters).out_fname.c_str(), ios::binary | ios::in);
	if (!out_file) {
		cout << "Error open " << (*p_Parameters).out_fname << " file";
		cout << "Maybe, file doesn't exist" << endl;
		exit(1);
	}

	//If any of parameters traits||snps||chi use, this block fill traits.valueset and snps.valueset
	//(if their values are default all)
	if((*p_Parameters).traits.use||(*p_Parameters).snps.use||(*p_Parameters).chi.use){

			if((*p_Parameters).traits.value=="all"||(*p_Parameters).traits.value=="None"){
				for(int i=0;i<(*(*p_iout_file).labels.trait_names).size();i++)
				(*p_Parameters).traits.valueset.insert(i);
			cout<<"TRAITS VALUE SET CHANGED TO ALL"<<endl;
		}

		if((*p_Parameters).snps.value=="all"||(*p_Parameters).snps.value=="None"){
			for(int i=0;i<(*(*p_iout_file).labels.snp_names).size();i++)
				(*p_Parameters).snps.valueset.insert(i);
			cout<<"SNPS VALUE SET CHANGED TO ALL"<<endl;
		}
	}

	if(((*p_Parameters).traits.use||(*p_Parameters).snps.use)&&!(*p_Parameters).chi.use){
		//mkdir("data");
		write_data(out_file);
	}

	if((*p_Parameters).chi.use&&!(*p_Parameters).dataslim.use){
		//mkdir("chi_data");
		write_data_chi(out_file);

	}

	if((*p_Parameters).dataslim.use){
		//mkdir("slimdata");
		write_slim_data(out_file);

	}

	if((*p_Parameters).heritabilities.use){
		//mkdir("estimates");
		write_herest(out_file);
	}
}
