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
	cout<<"End_write_data_dimension\t"<< double(clock()) / CLOCKS_PER_SEC <<" sec" << endl;
}

void Reshuffle::write_snpnames(ofstream& txt_snpnames){

	if ((*p_Parameters).snpnames.value == "all"){
		for (int i=0;i<(*(*p_iout_file).labels.snp_names).size();i++)
			(*p_Parameters).snpnames.numbersset.insert(i);
		cout<<"SNPNAMES VALUE SET CHANGED TO ALL"<<endl;
	}
	for(set<int>::iterator it= (*p_Parameters).snpnames.numbersset.begin();it!=(*p_Parameters).snpnames.numbersset.end();++it)
		txt_snpnames << "SNP #"<<(*it+1)<<"\t"<<(*(*p_iout_file).labels.snp_names)[*it]<<endl;
	cout<<"END WRITE SNPNAMES\t" << double(clock()) / CLOCKS_PER_SEC <<" sec" << endl;
}

void Reshuffle::write_traitnames(ofstream& txt_traitnames){

	if ((*p_Parameters).traitnames.value == "all"){
		for (int i=0;i<(*(*p_iout_file).labels.trait_names).size();i++)
			(*p_Parameters).traitnames.numbersset.insert(i);
		cout<<"TRAITNAMES VALUE SET CHANGED TO ALL"<<endl;
	}
	for(std::set<int>::iterator it= (*p_Parameters).traitnames.numbersset.begin();it!=p_Parameters->traitnames.numbersset.end();++it)
		txt_traitnames<<"TRAIT #"<<(*it+1)<<"\t"<<(*(*p_iout_file).labels.trait_names)[*it]<<endl;
	cout<<"END WRITE TRAITNAMES\t" << double(clock()) / CLOCKS_PER_SEC <<" sec" << endl;
}

void Reshuffle::write_data(ifstream& out_file,ofstream& data){
	out_file.seekg(0, ios_base::beg);
	cout << "Start_write_data\t" << double(clock()) / CLOCKS_PER_SEC <<" sec" << endl;
	data.precision(PRECISION_DOUBLE);
	data<<"SNP\t";
	data<<	"Trait\t";
	for (int beta = 0;	beta < (*(*p_iout_file).labels.beta).size(); beta++)
		data << (*(*p_iout_file).labels.beta)[beta] << "\t";
	for (int se = 0;se < (*(*p_iout_file).labels.se).size(); se++)
		data << (*(*p_iout_file).labels.se)[se] << "\t";
	for (int cov = 0;cov < (*(*p_iout_file).labels.cov).size(); cov++)
		data << (*(*p_iout_file).labels.cov)[cov] << "\t";
	data << endl;
	for (set<int>::iterator trait= (*p_Parameters).traits.numbersset.begin();trait!=(*p_Parameters).traits.numbersset.end();trait++) {
		double* buf = new double[per_trait_per_snp];
		int oldPos = 0;
		char s[30];
		for (set<int>::iterator snp= (*p_Parameters).snps.numbersset.begin();snp!=(*p_Parameters).snps.numbersset.end();snp++) {
			data << (*(*p_iout_file).labels.snp_names)[*snp] << "\t";
			data << (*(*p_iout_file).labels.trait_names)[*trait]<<"\t";
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
				data << s << "\t";
			}
			data << endl;
		}
		delete buf;
		//txt_trait.close();
		cout << "End_write_trait\t" << (*(*p_iout_file).labels.trait_names)[*trait] << " "<< double(clock()) / CLOCKS_PER_SEC <<" sec" << endl;
	}
	cout << "Finish_write_data\t" << double(clock()) / CLOCKS_PER_SEC <<" sec" << endl;
}

void Reshuffle::write_data_chi(ifstream& out_file,ofstream& txt_chi){
	out_file.seekg(0, ios_base::beg);
	double chi = 0;
	double CheckChi = (*p_Parameters).chi.value == "all" ? -1.0 : atof((*p_Parameters).chi.value.c_str());
	cout << "Start_write_chi_data=" << double(clock()) / CLOCKS_PER_SEC <<" sec" << endl;
	txt_chi.precision(PRECISION_DOUBLE);
	txt_chi << "SNP\t";
	txt_chi << "Trait\t";
	for (int beta = 0;	beta < (*(*p_iout_file).labels.beta).size(); beta++)
		txt_chi << (*(*p_iout_file).labels.beta)[beta] << "\t";
	for (int se = 0;se < (*(*p_iout_file).labels.se).size(); se++)
		txt_chi << (*(*p_iout_file).labels.se)[se] << "\t";
	for (int cov = 0;cov < (*(*p_iout_file).labels.cov).size(); cov++)
		txt_chi << (*(*p_iout_file).labels.cov)[cov] << "\t";
	txt_chi << "Chi2" << endl;
	for (set<int>::iterator trait= (*p_Parameters).traits.numbersset.begin();trait!=(*p_Parameters).traits.numbersset.end();trait++) {
		//ofstream txt_chi(create_filename("chi_data//chi", (*(*p_iout_file).labels.trait_names)[*trait]).c_str());
		double* buf = new double[per_trait_per_snp];
		int oldPos = 0;
		char s[30];
		for (set<int>::iterator snp= (*p_Parameters).snps.numbersset.begin();snp!=(*p_Parameters).snps.numbersset.end();snp++) {
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
				txt_chi << (*(*p_iout_file).labels.trait_names)[*trait]<<"\t";
				for (int i = 0; i < per_trait_per_snp; i++) {
					sprintf(s, "%.15g", buf[i]);
					txt_chi << s << "\t";
				}
				txt_chi << chi << endl;
			}
		}
		delete buf;
		//txt_chi.close();
		cout << "End_write_chi_trait\t" << (*(*p_iout_file).labels.trait_names)[*trait] << " "<< double(clock()) / CLOCKS_PER_SEC <<" sec" << endl;
	}
	cout << "Finish_write_chi_data\t" << double(clock()) / CLOCKS_PER_SEC <<" sec" << endl;
}

void Reshuffle::write_slim_data(ifstream& out_file, ofstream& txt_slim){
	out_file.seekg(0, ios_base::beg);
	set<int> goodtraits;
	set<int> goodsnps;
	double chi = 0;
	if((*p_Parameters).chi.value=="all"||(*p_Parameters).chi.value=="None"){
		cout << "ERROR: " << "Chi value doesn't set"<<endl;
		cout << "Please, set Chi value to get slim data" << endl;
		exit(1);
	}
	double CheckChi = atof((*p_Parameters).chi.value.c_str());
	for (set<int>::iterator trait= (*p_Parameters).traits.numbersset.begin();trait!=(*p_Parameters).traits.numbersset.end();trait++) {
		double* buf = new double[per_trait_per_snp];
		int oldPos = 0;
		char s[30];
		for (set<int>::iterator snp= (*p_Parameters).snps.numbersset.begin();snp!=(*p_Parameters).snps.numbersset.end();snp++) {
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
	txt_slim << "SNP\t";
	txt_slim << "Trait\t";
	for (int beta = 0;	beta < (*(*p_iout_file).labels.beta).size(); beta++)
		txt_slim << (*(*p_iout_file).labels.beta)[beta] << "\t";
	for (int se = 0;se < (*(*p_iout_file).labels.se).size(); se++)
		txt_slim << (*(*p_iout_file).labels.se)[se] << "\t";
	for (int cov = 0;cov < (*(*p_iout_file).labels.cov).size(); cov++)
		txt_slim << (*(*p_iout_file).labels.cov)[cov] << "\t";
	txt_slim << "Chi2" << endl;
	for (set<int>::iterator trait= goodtraits.begin();trait!=goodtraits.end();trait++) {
		double* buf = new double[per_trait_per_snp];
		int oldPos = 0;
		char s[30];
		for (set<int>::iterator snp= goodsnps.begin();snp!=goodsnps.end();snp++) {
			txt_slim << (*(*p_iout_file).labels.snp_names)[*snp] << "\t";
			txt_slim << (*(*p_iout_file).labels.trait_names)[*trait]<<"\t";
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
	}
	cout <<"End_write_slim_data\t" << double(clock()) / CLOCKS_PER_SEC <<" sec" << endl;
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

void Reshuffle::write_herest(ifstream& out_file, ofstream& herest){
	ofstream txt_est("estimates.txt");
	out_file.seekg(herest_startpos, ios_base::beg);
	if (p_Parameters->heritabilities.value == "all")
		for(int i=0;i<(*(p_iout_file->labels.trait_names)).size();i++)
			p_Parameters->heritabilities.numbersset.insert(i);
	txt_est.precision(PRECISION_DOUBLE);
	txt_est<<"\t";
	for (set<int>::iterator trait= p_Parameters->heritabilities.numbersset.begin();trait!=p_Parameters->heritabilities.numbersset.end();trait++)
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
		for (std::set<int>::iterator trait= p_Parameters->heritabilities.numbersset.begin();trait!=p_Parameters->heritabilities.numbersset.end();++trait) {
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
			for (std::set<int>::iterator trait= p_Parameters->heritabilities.numbersset.begin();trait!=p_Parameters->heritabilities.numbersset.end();++trait) {
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
	cout << "End_write_estimates\t" << double(clock()) / CLOCKS_PER_SEC <<" sec" << endl;
}

void Reshuffle::run(){
	if((*p_Parameters).datadims.use){
		ofstream datadims((*p_Parameters).datadims.outfile.c_str());
		write_datadims(datadims);
	}
	if((*p_Parameters).snpnames.use){
		ofstream snpnames((*p_Parameters).snpnames.outfile.c_str());
		write_snpnames(snpnames);
	}

	if((*p_Parameters).traitnames.use){
		ofstream traitnames((*p_Parameters).traitnames.outfile.c_str());
		write_traitnames(traitnames);
		//delete traitnames;
	}

	//Open *.out to read data
	ifstream out_file((*p_Parameters).out_fname.c_str(), ios::binary | ios::in);
	if (!out_file) {
		cout << "Error open " << (*p_Parameters).out_fname << " file";
		cout << "Maybe, file doesn't exist" << endl;
		exit(1);
	}

	//If any of parameters traits||snps||chi use, this block fill traits.numbersset and snps.numbersset
	//(if their values are default all)
	if((*p_Parameters).traits.use||(*p_Parameters).snps.use||(*p_Parameters).chi.use||!(*p_Parameters).defaultstate){

			if((*p_Parameters).traits.value=="all"||(*p_Parameters).traits.value=="None"){
				for(int i=0;i<(*(*p_iout_file).labels.trait_names).size();i++)
				(*p_Parameters).traits.numbersset.insert(i);
			cout<<"TRAITS VALUE SET CHANGED TO ALL"<<endl;
		}

		if((*p_Parameters).snps.value=="all"||(*p_Parameters).snps.value=="None"){
			for(int i=0;i<(*(*p_iout_file).labels.snp_names).size();i++)
				(*p_Parameters).snps.numbersset.insert(i);
			cout<<"SNPS VALUE SET CHANGED TO ALL"<<endl;
		}
	}

	if((((*p_Parameters).traits.use||(*p_Parameters).snps.use)&&!(*p_Parameters).chi.use)||!(*p_Parameters).defaultstate){
		ofstream data;
		if((*p_Parameters).traits.outfile!="data.txt"){
			data.open((*p_Parameters).traits.outfile.c_str());
		}else
			if((*p_Parameters).snps.outfile!="data.txt"){
				data.open((*p_Parameters).snps.outfile.c_str());
			}else
				data.open("data.txt");
		write_data(out_file,data);
	}

	if((*p_Parameters).chi.use&&!(*p_Parameters).dataslim.use){
		ofstream chi_data((*p_Parameters).chi.outfile.c_str());
		write_data_chi(out_file,chi_data);

	}

	if((*p_Parameters).dataslim.use){
		ofstream dataslim((*p_Parameters).dataslim.outfile.c_str());
		write_slim_data(out_file,dataslim);

	}

	if((*p_Parameters).heritabilities.use){
		ofstream herest((*p_Parameters).heritabilities.outfile.c_str());
		write_herest(out_file,herest);
	}
}
