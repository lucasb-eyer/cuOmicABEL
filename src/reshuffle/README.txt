Reshuffle ver 0.000

Parameters:

For the first it should be file name without iout and out extensions:
	(--B_eigen_1000) 
	if your files are B_eigen.iout and B_eigen.out

Data dimensions 
	(--datadims)  Gives back t, m, p
output: "datadims//datadims.txt"

SNP names
	default: (--snpnames) : all
	by index (--snpnames=27) : name of snp#27
	by index range, combination (--snpnames=27,2-12,4-20) : name of snp #2-20,27
output: "snpnames//snpnames.txt"
	
Trait names 
	default (--traitnames): all
	by index (--traitnames=27) : name of trait #27
	by index range, combination (--traitnames=27,2-12,15)
output: "traitnames//traitnames.txt"

Heritabilities, sigma, res_sigma, estimates
	Default: (--heritabilities) : all traits
	Reange,indexes (--heritabilities=1-10,4,5-12) : for traits #1-12
output: "estimates//estimates.txt"	
Results - association
	by SNP
		Dafault (--snp) : for all snp for all traits
		indexes, ranges (--snp=12,100-1000,500-10000,22000) :  for all traits
	output: "data//[trait].txt"
	by trait
		Dafault (--trait) : for all snp for all traits
		indexes,range (--trait=1-10,12) : for snp #12 for all traits
	output: "data//[trait].txt"
	Chi2 more than some threshold
		(--chi=20) for snps,which chi>20 for all traits
	output: "chi_data//[trait].txt"
	All combinations of traits,snps, Chi2
		(--traits=1-2--snps=1-1000--chi=15) for traits #1-2 for snps#1-1000, which chi>15
		output: "chi_data//[trait].txt"
		
Just write data with chi
	(--chi) : write data(with chi column) for all snps for all traits
	all combinations are supported
	output: "chi_data//[trait].txt"
Dataslim
	Creating a sub-matrix by Chi2>X (--chi=X--dataslim)
	output: "slim_data//slim_[trait].txt"
	
Test
	(--test) run all tests
	output: test.txt

Examples:
Reshuffle.exe B_1112_NA_clear_RNA_nocovar --traitnames=1-4,1--snpnames=1-10,20-30,46--traits=1--snps=1-100

outputs from B_1112_NA_clear_RNA_nocovar.iout and B_1112_NA_clear_RNA_nocovar.out:

	traitnames//traitnames.txt : with trait's names hgta,tga,tca,ldla
	snpnames//snpnames.txt Names of snps #1-10,20-30,46
	data//trait_hgta.txt : result of association for trait "hgta" for snps #1-100