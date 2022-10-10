/*
 * smoothT.cpp
 *
 *  Created on: Oct 7, 2019
 *      Author,Copyright: Rene Staritzbichler
 */

#include <iostream>
#include <vector>
#include <unistd.h>
#include <memory>
#include <map>
#include <math.h>
#include <fstream>
#include <algorithm>
#include <cctype>
#include <iterator>
#include <sstream>
#include <iomanip>
#include <string>
#include <list>

#include "string_functions.h"
#include "reader.h"
#include "math_functions.h"
#include "path.h"  // node and edges


// g++ -O3 smoothT.cpp -o smoothT  -std=c++11


struct SortPair
{
	bool operator()( const std::pair< float, float> &A , const std::pair< float, float> &B) const
	{
		if( A.first == B.first)
		{
			return A.second < B.second;
		}
		return A.first < B.first;
	}
};



void
Options()
{
	std::cout << "\t\t\t=====================================" << std::endl;
	std::cout << "\t\t\t=====================================" << std::endl;
	std::cout << "\t\t\t===========    SmoothT    ===========" << std::endl;
	std::cout << "\t\t\t=====================================" << std::endl;
	std::cout << "\t\t\t=====================================\n\n" << std::endl;
	std::cout << "smoothT allows for the determination of low energy pathways from an ensemble of PDB models/poses.";
	std::cout << "The ensemble has to be sampled densely enough to have sufficient neighbors within the desired RMSD cutoff (specified by '-d' flag).\n\n\n";


	std::cout << "######     FLAGS OVERVIEW     ######\n";
	std::cout << "###  REQUIRED FLAGS:  ###\nINPUT:\n";
	std::cout <<
			"\t-b FILENAME of first node in pathway\n"
			"\t-e FILENAME of last node in pathway\n"
			"\t-l FILENAME containing list of pdb files to create pathway\n"
			"\t-d RMSD maximum distance for structures considered to be neighbors(=connected)\n"
			"\t-i ENERGY_IDENTIFYER string in PDB file followed by energy value. IF not set, energies have to be passed in FILENAME.\n"
//			"OR:\n"
//			"\t-m read energies of all nodes from '-l LIST' instead from PDBs\n"
			"OUTPUT:\n"
			"\t-o DIR output directory for pdb, energy\n" << std::endl;
	std::cout << "###  OPTIONAL FLAGS:  ###\n"
			"\t-a CA, C, N, ... atomtypes as in PDB used for RMSD calculation\n"
			"\t-c CHAIN1 ... chains\n"
//			"\t-z ZERO_ENERGY interaction energy when molecules are separated (default 0.0)\n"
//			"\t-w WIDTH of transition for estimation of energy barriers (default: 5.0)\n"
			"\t-n NR of top scoring pathways that are written to output directory (default: 200)\n"
			"\t-k NR of barriers :: first sorting criterium (default: 5)\n"
			"\t-z NR of pathways per barrier :: second sorting criterium (default: 20000)\n"
			"\t-f FILENAME of alignment\n\n" << std::endl;
}


void
Details()
{

	std::cout << "######     DETAILS     ######\n";
	std::cout << "The user has to provide a list of all PDB files to consider for pathway construction: specified by '-l' flag\n";
	std::cout << "Additionally, from that list, start and end point of the pathway has to be specified: '-b' and '-e' flag\n\n";
	std::cout << "PDBs need to be structured identically!!\n";
	std::cout << "Each PDB needs the energy value being added: 'ENERGY_IDENTIFIER -123.455'" << std::endl;
	std::cout << "The ENERGY_IDENTIFIER string is specified with the '-i' flag. \n\n";
	std::cout << "smoothT constructs pathways by determining all neighbors for each PDB that are within the RMSD cutoff, specified by '-d' flag. \n";
	std::cout << "For smaller ligands binding to larger receptors, it is useful to restrict the RMSD calculation to the chain of the ligand by applying the '-c' flag.\n";
	std::cout << "Mostly, it is sufficient to use C-alpha atoms only for RMSD calculation.\n";
	std::cout << "If no atom types is specified using the '-a' flag, all atoms are used." << std::endl;
	std::cout << "\n";
	std::cout << "##   OUTPUT:  ##\n";
	std::cout << "All output files are written into the directory specified by '-o' flag.\n";
	std::cout << "Three kind of files are provided: \n";
	std::cout << "\tpath_N_I.pdb: PDB file containing all models forming a pathway\n";
	std::cout << "\t\tcan be used directly in e.g. VMD or MDsrv to visualize as dynamic transition\n";
	std::cout << "\tenergy_N_I.txt: text file containing the energies of all PDBs contained in the pathway.\n";
	std::cout << "\t\tcan be plotted using e.g. gnuplot ('plot \"energy_0_0.txt\" us 3:4 w l') x is accumulated RMSD, y the energy. \n";
	std::cout << "\t\taccumulated RMSD allows to visualize them subsequently, but actual meaning has only the difference in RMSD of neighboring models.\n";
	std::cout << "\ttransition_N_I.txt: text file containing additionally a simple estimation of the energy barrier between subsequent models.\n";
	std::cout << "\t\tthe transition is based on the assumption that each model represents a local energy minimum\n";
	std::cout << "\t\tmeaning any step away from this minimum should increase the energy,\n";
	std::cout << "\t\treference is the unbound state, specified by 'z' flag\n";
	std::cout << "\t\ta sinoid transition is assumed, its maximum is reached after a RMSD value that is specified via the '-w' flag.\n";
	std::cout << "\n";
	std::cout << "All output indices follow this schema:\n";
	std::cout << "\tN is the ranking in energy barrier (the first sorting criteria)\n";
	std::cout << "\tI is the ranking in sum of energies in pathway (the second sorting criteria)\n";
	std::cout << "\tExample: energy_0_0.txt contains the energies of all nodes in the best scoring pathway.\n";
	std::cout << "\n";
	std::cout << "##  ALIGNMENTS:  ##\n";
	std::cout << "In case there are two ensembles of PDB models that are similar but not identical, they can be used for pathway construction by passing an alignment using the '-f' flag. \n";
	std::cout << "An example: there are PDB structures solved for distinct binding states of two molecules and you want to model the transition between them. It is rather common that the molecules are modified for crystallization, thus may differ in sequence although they belong to the same protein.\n";
	std::cout << "Here independent docking can be performed for both for sampling the conformational space.\n";
	std::cout << "Using the alignment, the RMSD is calculate based on aligned positions.\n";
	std::cout << "\n";
	std::cout << "If an alignment is passed, the sequences have to match the coordinate section of the PDB file! Sometimes the sequences in the coordinate section (ATOM lines) differ from the SEQRES section or the ones in FASTA files.\n\n";
	std::cout << "\n";
	std::cout << "Do not hesitate to ask further questions or to send us feedback!\n";
	std::cout << "\n";
	std::cout << "SUCCESSFUL MODELING !!!\n";
	std::cout << "\n";
}


void Help()
{
	Options();
	Details();
}


void Author()
{
	std::cout << "Author: \nRené Staritzbichler, rene@staritzbichler.com, http://proteinformatics.org , 10.10.2019\n" << std::endl;
	std::cout << "Please cite: \nRené Staritzbichler, Nikola Ristic, Peter W. Hildebrand, 2022, \n\"SmoothT – a server constructing basic transition pathways from conformational ensembles for interactive visualization and enhanced sampling\"\n\n\n" << std::endl;
}



int main(  int ARGC, char ** ARGV)
{
	const clock_t start_time = clock();


	for( int i = 0; i < ARGC; ++i)
		std::cout << ARGV[i] << "  ";
	std::cout << "\n\n" << std::endl;
	Author();
	if( ARGC == 1)
	{
		Options();
		return 0;
	}

	//  ===== DECLARATIONS =====
	std::vector< std::string>
    	atom_types;
	std::vector< char>
		chains;
	std::string
		energy_identifyer,
		file,
		outdir,
		first_file,
		last_file,
		list_file,
		alignment_file;
	std::vector< std::shared_ptr< Node > >
		all;

	std::vector< std::map<char,std::string> >
		alignments;
	float
		offset = 0.0,
		width = 5.0,
		zero_energy = std::numeric_limits<float>::max(),
		max_dist;
	int
		nr_output = 200,
		max_nr_paths = 20000,
		max_nr_barriers = 5,
		c, id = 0;
//	bool
//		read_energies_from_list = false;

	//enum eSortType {maxE,sumE,maxRMSD};

	while( ( c = getopt( ARGC, ARGV, "hb:e:l:d:p:o:t:a:c:f:i:n:w:z:k:") ) != -1 )
    {
		switch(c)
		{
		case 'h':
			Help();
			return 0;
		case 'b':
			if( optarg) first_file = optarg;
			break;
		case 'e':
			if( optarg) last_file = optarg;
			break;
		case 'l':
			if( optarg) list_file = optarg;
			break;
		case 'd':
			if( optarg) max_dist = std::atof( optarg);
			break;
		case 'o':
			if( optarg) outdir = optarg;
			break;
		case 'a':
			if( optarg) atom_types.push_back( optarg);
			break;
		case 'c':
			if( optarg) chains.push_back( optarg[0]);
			break;
		case 'n':
			if(optarg) nr_output = std::atoi(optarg);
			break;
		case 'k':
			if(optarg) max_nr_paths = std::atoi(optarg);
			break;
		case 'z':
			if(optarg) max_nr_barriers = std::atoi(optarg);
			break;
//		case 'm':
//			read_energies_from_list = true;
//			break;
		case 'f':
			if( optarg) alignment_file = optarg;
			alignments = ReadAlignments( alignment_file);
			id = 1;
			break;
		case 'i':
			if( optarg) energy_identifyer = optarg;
			break;
		case 'w':
			if( optarg ) width = atof( optarg );
			break;
//		case 'z':
//			if( optarg)	offset = atof(optarg);
//			break;
		}
    }

	if( alignments.size() != 2){
		alignments.push_back( std::map<char,std::string>() );}

	Edge
		edge;
	std::shared_ptr<Node>
		node,
		first_node( new Node( first_file, energy_identifyer, atom_types, chains, alignments[ 0])),
		last_node( new Node( last_file, energy_identifyer, atom_types, chains, alignments[id]));   // last node needed - or just name?

	std::cout << *first_node << std::endl;
	std::cout << *last_node << std::endl;

	std::cout << "rmsd of start and endpoint: " << RMSD( first_node->GetPos(), last_node->GetPos()) << std::endl;

//	std::map<std::string,int>
//		grouping;

	std::ifstream
		in( list_file);
	if( !in)
	{ std::cerr << "ERROR: not opened: " << list_file << std::endl; exit(1);}
	std::cout << "read list of nodes" << std::endl;
	int loc = 0;

	//while( in >> file >> id)
	while( in >> file )
    {
		if( ++loc % 10000 == 0) std::cout << loc << std::endl;
		//grouping[file] = id;
		std::shared_ptr<Node> tmp( new Node( file, energy_identifyer, atom_types, chains, alignments[id]));
		if( tmp->GetPos().size() > 0){
			all.push_back( tmp ); }
		zero_energy = std::min( zero_energy, tmp->GetEnergy());
    }

	all.push_back( last_node); // !!! ???

	in.close();
	in.clear();

	std::cout << all.size() + 1 << " nodes read" << std::endl;
	std::cout << "zero energy: " << zero_energy << std::endl;

	std::vector< std::shared_ptr< Node> >
		current;  // jeweiliger Startpunkt bei jeder iteration
	current.push_back( first_node);

//// TODO CHECK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//	//////  switch first (and last??) node 'off'  /////////////////////////////
//	float
//		initial_first_energy = first_node->GetEnergy();
//	current[0]->SetEnergy( 0.0);
//	std::cout << "switch off first node: " << first_node->GetEnergy()  << std::endl;
//	//////////////////////////////////////////////////////////////

	clock_t now = clock();
	std::cout << "TIMER: reading: " << float( now - start_time) / CLOCKS_PER_SEC << std::endl;

	std::cout << "STATUS: constructing graph ... " << std::endl;

	std::vector< std::vector< std::shared_ptr< Node > > >
		generations;


	std::cout << "unshifted energies, first node: " << first_node->GetEnergy() << " last node: " << last_node->GetEnergy() << std::endl;
//	std::cout << "unshifted barriers, first node: " << first_node->GetBarrier() << " last node: " << last_node->GetBarrier() << std::endl;
//	std::cout << "unshifted integrals, first node: " << first_node->GetSum() << " last node: " << last_node->GetSum() << std::endl;

	for( auto a : all)
	{
		Shift( a, zero_energy);  //  TODO: shift sum, barrier in nodes??
	}
	Shift( first_node, zero_energy);

	std::cout << "shifted energies, first node: " << first_node->GetEnergy() << " last node: " << last_node->GetEnergy() << std::endl;
//	std::cout << "shifted barriers, first node: " << first_node->GetBarrier() << " last node: " << last_node->GetBarrier() << std::endl;
//	std::cout << "shifted integrals, first node: " << first_node->GetSum() << " last node: " << last_node->GetSum() << std::endl;
	////////////////////////////////////////////////////////
	////////////////      BUILD GRAPH    ///////////////////
	// network/graph building loop
	while( Iterate( current, all, max_dist, last_node, generations) ){}
	////////////////////////////////////////////////////////
	clock_t now2 = clock();
	std::cout << "TIMER: graph construction: " << float( now2 - now) / CLOCKS_PER_SEC << std::endl;
	std::cout << "STATUS: graph constructed" << std::endl;

	generations.push_back( std::vector< std::shared_ptr< Node > >( 1, last_node));

//	int cc = 0;
//	for( std::vector< std::vector< std::shared_ptr< Node > > >::const_iterator itr = generations.begin(); itr != generations.end(); ++itr, ++cc)
//		for( std::vector< std::shared_ptr< Node > >::const_iterator jtr = itr->begin(); jtr != itr->end(); ++jtr)
//		{
//			if( *jtr == last_node)
//			{
//				std::cout << "CHECK: last node found in generation " << cc <<  " / " << generations.size() << std::endl;
//			}
//			for( std::vector< Edge>::const_iterator etr = (*jtr)->GetParentEdges().begin(); etr != (*jtr)->GetParentEdges().end(); ++etr)
//				if( etr->GetNode() == first_node)
//				{
//					std::cout << "CHECK: first node found as parent in generation " << cc <<  " / " << generations.size() << " " << (*jtr)->GetParentEdges().size() <<  std::endl;
//				}
//		}

//	std::cout << "\nBEFORE GENERATIONWALK:\nfirst node: " << *first_node << " \n\nlast node: " << *last_node << std::endl;

	GenerationWalk( generations);

//	std::cout << "\nAFTER GENERATIONWALK: \nfirst node: " << *first_node << " \n\nlast node: " << *last_node << std::endl;

	now = clock();
	std::cout << "TIMER: generation walk: " << float( now - now2) / CLOCKS_PER_SEC << std::endl;

	std::vector< std::shared_ptr< Node> >
		final_path;

//	Shift( last_node, -zero_energy);  // UNCOMMENT THIS AFTER TEST ROUND !!!!

	Backtrace( last_node, final_path);

	now2 = clock();
	std::cout << "TIMER: backtrace: " << float( now2 - now) / CLOCKS_PER_SEC << std::endl;

	std::cout << "path length: " << final_path.size() << std::endl;

	std::reverse( final_path.begin(), final_path.end());

//	for( auto a : final_path)
//	{
//		Shift( a, -zero_energy);  //  TODO: shift sum, barrier in nodes??
//	}

	if( final_path[0] == first_node)
	{
		std::cout << "first node found in final path" << std::endl;
	}
	if( final_path.back() == last_node)
	{
		std::cout << "last node found in final path" << std::endl;
	}

	std::cout << "Best path: barrier: " << last_node->GetBarrier() << " integral: " << last_node->GetSum() << std::endl;

//	Write( final_path, 0.0, last_node->GetBarrier(), outdir);
	Write( final_path, zero_energy, last_node->GetBarrier() + zero_energy, outdir); // UNCOMMENT THIS AFTER TEST ROUND !!!!
	now = clock();
	std::cout << "TIMER: writing: " << float( now - now2) / CLOCKS_PER_SEC << std::endl;

	std::cout << "STATUS: finished" << std::endl;
	return 0;
};


