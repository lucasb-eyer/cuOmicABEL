/*
 * iout_file.h
 *
 *  Created on: 28.02.2013
 *      Author: Sodbo
 */

#ifndef IOUT_FILE_H_
#define IOUT_FILE_H_

#include <fstream>
#include <iostream>
#include <vector>
#include "Parameters.h"
#include <stdlib.h>
#include <stdio.h>


using namespace std;

//Contain meta data
class iout_header {
public:
	int type; // Double -> 6 as in Databel
	int nbytes; // Double -> 8
	int p; // # of covariates + 2 (intercept and SNP)
	int m; // # of SNPs
	int t; // # of traits
	int tile_m; // More on these two later
	int tile_t;
	int namelength; // 32 as in Databel

};

//Overloading "cout" operator
ostream &operator <<(ostream &,iout_header);

//Contain labels
class labels_data {
public:
	int number; //number of labels
	vector<string>* beta;
	vector<string>* se;
	vector<string>* cov;
	vector<string>* snp_names;
	vector<string>* trait_names;
};

ostream &operator <<(ostream &,labels_data);

class iout_file{
public:
	iout_header header;
	labels_data labels;
	iout_file(Parameters &);
	int tilecoordinates(int, int);
};
#endif /* IOUT_FILE_H_ */
