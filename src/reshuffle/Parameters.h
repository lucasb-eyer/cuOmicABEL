/*
 * Parameters.h
 *
 *  Created on: 28.02.2013
 *      Author: Sodbo
 */

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <iostream>
#include <string>
#include <set>
#include <vector>
using namespace std;

class Parameter {

public:
	string name; 	//name of parametr
	bool use; 	// Parameter use or not (Is paameter in command line?)
	string value; 	// value of parametr,chars after "=" symbol
	set<int> numbersset;
	set<string> namesset;
	string delfromcmdline(string);
	void setbynames(vector<string>);
	string outfile;
	Parameter(string,string,string); 	//constructor
	Parameter();		//default constructor
};

ostream &operator <<(ostream &, Parameter);

class Parameters {
public:
	Parameter info; // Write info about progeamm's run
	string iout_fname; //iout_file_name
	string out_fname;
	Parameter datadims;
	Parameter snpnames;
	Parameter traitnames;
	Parameter traits;
	Parameter snps;
	Parameter heritabilities;
	Parameter chi;
	Parameter dataslim;
	Parameter test;
	Parameters();
	Parameters(int, char*[]);		//	Constructor from cmdline
	static string get_cmd_line(int,char*[]);
	bool defaultstate;//Check: are there parameters instead of infile's path

};

ostream &operator <<(ostream &, Parameters);

#endif /* PARAMETERS_H_ */
