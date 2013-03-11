/*
 * iout_file.cpp
 *
 *  Created on: 28.02.2013
 *      Author: Sodbo
 */
#include "iout_file.h"

iout_file::iout_file(Parameters& Params){
	//Fill iout_header
	ifstream iout_f(Params.iout_fname.c_str(), ios::binary | ios::in);
	if (!iout_f) {
		cout << "Error opening " << Params.iout_fname << " file" << endl;
		cout << "Maybe, file doesn't exist" << endl;
		exit(1);
	}
	iout_f.read((char *)&header.type, sizeof(int));
	iout_f.read((char *)&header.nbytes, sizeof(int));
	iout_f.read((char *)&header.p, sizeof(int));
	iout_f.read((char *)&header.m, sizeof(int));
	iout_f.read((char *)&header.t, sizeof(int));
	iout_f.read((char *)&header.tile_m, sizeof(int));
	iout_f.read((char *)&header.tile_t, sizeof(int));
	iout_f.read((char *)&header.namelength, sizeof(int));

	//Fill labels_data
	iout_f.seekg(sizeof(header), ios_base::beg); //Change position of read after iout_header
	labels.number =(header.p * 2 + header.m + header.t + header.p * (header.p - 1) / 2);
	char *labels_char = new char[labels.number * header.namelength];
	iout_f.read(labels_char, labels.number * header.namelength);
	int beta_end = header.p;
	int se_end = beta_end+header.p;
	int cov_end = se_end+header.p * (header.p - 1) / 2;
	int snp_end = cov_end+header.m;
	int traits_end = snp_end+header.t;
	int i=0;

	labels.beta=new vector<string>;
	labels.se=new vector<string>;
	labels.cov=new vector<string>;
	labels.snp_names=new vector<string>;
	labels.trait_names=new vector<string>;

	for (; i < beta_end; i++)
		(*labels.beta).push_back(string(labels_char + i * header.namelength));
	for (; i < se_end; i++)
		(*labels.se).push_back(string(labels_char + i * header.namelength));
	for (; i < cov_end;i++)
		(*labels.cov).push_back(string(labels_char + i * header.namelength));
	for (; i < snp_end; i++)
		(*labels.snp_names).push_back(string(labels_char + i * header.namelength));
	for (; i< traits_end; i++)
		(*labels.trait_names).push_back(string(labels_char + i * header.namelength));

	delete labels_char;
	iout_f.close();
}

//Overloading operator cout for iout_header
ostream &operator <<(ostream &os,iout_header header) {
	os << "IOUT_HEADER [ type: " << header.type << "]" << endl;
	os << "IOUT_HEADER [ nbytes: " << header.nbytes << "]" << endl;
	os << "IOUT_HEADER [ p: " << header.p << "]" << endl;
	os << "IOUT_HEADER [ m: " << header.m << "]" << endl;
	os << "IOUT_HEADER [ t: " << header.t << "]" << endl;
	os << "IOUT_HEADER [ tile_m: " << header.tile_m << "]" << endl;
	os << "IOUT_HEADER [ tile_t: " << header.tile_t << "]" << endl;
	os << "IOUT_HEADER [ namelength: " << header.namelength << "]" << endl;

	return os;
}

ostream &operator <<(ostream &os,labels_data labels) {
	os <<"NUMBER OF BETAS\t"<<(*labels.beta).size()<<endl;
	os <<"NUMBER OF SE\t"<<(*labels.se).size()<<endl;
	os <<"NUMBER OF COVS\t"<<(*labels.cov).size()<<endl;
	os <<"NUMBER OF SNPNAMES\t"<<(*labels.snp_names).size()<<endl;
	os <<"NUMBER OF TRAITNAMES\t"<<(*labels.trait_names).size()<<endl;
}

int iout_file::tilecoordinates(int traitNo, int snpNo) {
	int tileCoor = 0;
	int t_tile = traitNo / header.tile_t;
	int t_off = traitNo % header.tile_t;
	int m_tile = snpNo / header.tile_m;
	int m_off = snpNo % header.tile_m;
	int per_trait_per_snp = header.p + header.p * (header.p + 1) / 2;
	tileCoor = (t_tile * (header.m * header.tile_t)
			+ m_tile* (header.tile_m* min(header.tile_t,header.t - header.tile_t * t_tile))
			+ t_off * (min(header.tile_m, header.m - header.tile_m * m_tile))
			+ m_off) * per_trait_per_snp * sizeof(double);
	return tileCoor;
}
