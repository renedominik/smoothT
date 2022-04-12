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


// g++ -O3 smoothT.cpp -o smoothT


void
Help()
{
	std::cout << "\t\t=====================================" << std::endl;
	std::cout << "\t\t=====================================" << std::endl;
	std::cout << "\t\t===========    smoothT    ===========" << std::endl;
	std::cout << "\t\t=====================================" << std::endl;
	std::cout << "\t\t=====================================\n\n" << std::endl;
	std::cout << "smoothT allows for the determination of low energy pathways from an ensemble of PDB models/poses.";
	std::cout << "The ensemble can be created either by means of Monte-Carlo algorithms or Molecular dynamics simulations.\n";
	std::cout << "The ensemble has to be sampled densely enough to have sufficient neighbors within the desired RMSD cutoff (specified by '-d' flag).\n\n\n";


	std::cout << "######     FLAGS OVERVIEW     ######\n";
	std::cout << "###  REQUIRED FLAGS:  ###\nINPUT:\n";
	std::cout << "\t-b FILENAME of first node in pathway\n"
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
			"\t-z ZERO_ENERGY interaction energy when molecules are separated (default 0.0)\n"
			"\t-w WIDTH of transition for estimation of energy barriers (default: 5.0)\n"
			"\t-n NR of top scoring pathways that are written to output directory (default: 200)\n"
			"\t-f FILENAME of alignment\n\n" << std::endl;


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

void Author()
{
	std::cout << "Author: Rene Staritzbichler, rene@staritzbichler.com, http://proteinformatics.org , 10.10.2019\n" << std::endl;
	std::cout << "please cite: \nStaritzbichler R, Hildebrand P, 2020\n\"smoothT\"\n\n\n" << std::endl;
}



int main(  int ARGC, char ** ARGV)
{
	for( int i = 0; i < ARGC; ++i)
		std::cout << ARGV[i] << "  ";
	std::cout << "\n\n" << std::endl;
	Author();
	if( ARGC == 1)
	{
		Help();
		return 0;
	}

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
		max_dist;
	int
		nr_output = 200,
		c, id = 0;
//	bool
//		read_energies_from_list = false;

	//enum eSortType {maxE,sumE,maxRMSD};

	while( ( c = getopt( ARGC, ARGV, "hb:e:l:d:p:o:t:a:c:f:i:n:w:z:") ) != -1 )
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
		case 'z':
			if( optarg)	offset = atof(optarg);
			break;
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
		if( ++loc % 1000 == 0) std::cout << loc << std::endl;
		//grouping[file] = id;
		std::shared_ptr<Node> tmp( new Node( file, energy_identifyer, atom_types, chains, alignments[id]));
		if( tmp->GetPos().size() > 0){
			all.push_back( tmp ); }
    }

	all.push_back( last_node); // !!! ???

	in.close();
	in.clear();

	std::cout << all.size() + 1 << " nodes read" << std::endl;

	std::vector< std::shared_ptr< Node> >
		current;
	current.push_back( first_node);


	//////  switch first node 'off'  /////////////////////////////
	float
		initial_first_energy = first_node->GetEnergy();
	current[0]->SetEnergy( 0.0);
	std::cout << "check: " << first_node->GetEnergy() << " " << current[0]->GetEnergy() << " (ideally not the same)" << std::endl;
	//////////////////////////////////////////////////////////////


	// network/graph building loop
	while( Iterate( current, all, max_dist, last_node) ){}

	// score sorted paths
	std::multimap< float, std::vector<int> >
		score_sorted_paths;
	std::vector<int>
		path;
	node = last_node;

	// go through all possible pathways through the graph, search best solution
	std::pair< float, float>
		min_max = std::make_pair( 1e10, -1e10);
	Backtrace( score_sorted_paths, path, min_max, node, first_node);
	min_max.first = std::min( min_max.first, initial_first_energy);
	min_max.second = std::max( min_max.second, initial_first_energy);
	std::cout << "energy min: " << min_max.first << " max: " << min_max.second << std::endl;

	std::cout << score_sorted_paths.size() << " complete paths" << std::endl;
	std::cout << "minimum energy path: " << score_sorted_paths.begin()->first << std::endl;

	std::cout << "all paths, sorted by energy threshold: " << std::endl;
	for( std::multimap< float, std::vector<int> >::const_iterator itr = score_sorted_paths.begin(); itr != score_sorted_paths.end(); ++itr)
	{
		std::cout << itr->first + initial_first_energy << ": ";
		for( int i : itr->second)
		{ std::cout << i << " ";}
		std::cout << std::endl;
    }

	std::cout << "nr parents of first node (should be 0): " << first_node->GetParentEdges().size() << std::endl;


	//////  switch first node back 'on'  /////////////////////////////
	first_node->SetEnergy( initial_first_energy);
	std::cout << "reset first node: " << first_node->GetEnergy() << std::endl;
	//////////////////////////////////////////////////////////////


	//////////////////////////////
	////   second sorting     ////
	//////////////////////////////
	std::cout << "sort paths belonging to same energy barrier by sum of energies " << std::endl;
	std::multimap< float, std::vector<int> >
		sub_sorted;
	std::vector<float>
		energies,
		rmsds;
	std::vector<std::string>
		names;
	float
		//integral,
		previous = score_sorted_paths.begin()->first;
	int
		count = 0,
		cc = 1,
		i1 = 0;
	std::cout << "print only first " << nr_output << std::endl;
	for( std::multimap< float, std::vector<int> >::const_iterator itr = score_sorted_paths.begin(); itr != score_sorted_paths.end() && count < nr_output; ++itr, ++cc)
	{

		if( itr->first != previous || cc == (signed) score_sorted_paths.size() || count == nr_output - 1)
		{
			previous = itr->first;
			int i2 = 0;

			// for plotting energy profiles along pathways with gnuplot
			std::ofstream
				gnu( outdir + "/gnu_" + std::to_string( i1) + ".txt");
			gnu << "plot ";

			// write sub_sorted
			for( std::multimap< float, std::vector<int> >::const_iterator jtr = sub_sorted.begin(); jtr != sub_sorted.end() && count < nr_output; ++jtr, ++i2, ++count )
			{
				std::string
					tail = "_" + std::to_string( i1) + "_" + std::to_string( i2) + ".txt",
					energy_file = "energy" + tail,
					transition_file = "transition" + tail,
					path_file = "path" + tail;
				gnu << "\"" + energy_file + "\" us 3:4 w l ti \"\"";   //  for plotting of profiles
//				gnu << "\"" + energy_file + "\" us 3:4 w l ti \"" + std::to_string( i2) + "\"";
				if( count + 1 < nr_output && i2 + 1 < (signed) sub_sorted.size())
				{ gnu << ", ";}

				std::cout << count << " " << itr->first << "  " << jtr->first << ":  ";
				for( int i3 : jtr->second)
				{ std::cout << i3 << " ";}
				std::cout << std::endl;

				// calc sum of energies of all nodes in path
				names.clear();
				energies.clear();
				rmsds.clear();
				node = last_node;
				for( std::vector<int>::const_iterator ktr = jtr->second.begin(); ktr != jtr->second.end(); ++ktr)
				{
					energies.push_back( node->GetEnergy() );
					names.push_back( node->GetName());
					edge = node->GetParentEdges()[*ktr];
					rmsds.push_back( edge.GetDistance() );
					node = edge.GetNode();
				}
				energies.push_back( node->GetEnergy() );
				names.push_back( node->GetName());

				float
					min_energy = *std::min_element( energies.begin(), energies.end() ),
					max_energy = *std::max_element( energies.begin(), energies.end() );

				std::cout << "sum: " << jtr->first << " barrier: " << itr->first << " energy range: " << min_energy << " " << max_energy << std::endl;

				std::ofstream
					out( outdir + "/" + energy_file),
					trance( outdir + "/" + transition_file);

				// write energies to file
				float
					rmsd = 0, sum = 0;
				int
					count = 1;

				if( !out){ std::cerr << "ERROR: not opened energy output file: " << std::endl;}
				for( unsigned int i = 1; i <= energies.size(); ++i)
				{
					sum += rmsd;
					out    << i << "\t" << rmsd << "\t" << sum << "\t" << energies[energies.size()-i] << "\t" << names[energies.size()-i] << std::endl;
					trance << count << "\t" << rmsd << "\t" << sum << "\t" << energies[energies.size()-i] << "\t" << names[energies.size()-i] << std::endl;
					rmsd = rmsds[rmsds.size() - i];
					if( i < energies.size())
					{ EnergyTransition( trance , rmsd , sum , offset , energies[energies.size()-i] , energies[energies.size()-i-1] , 20 , count, 0.2 , width );}
				}
				out << "# rmsd: " << rmsd << " barrier: " << itr->first << " sum: " << jtr->first << std::endl;
				out.close(); out.clear();
				trance.close(); trance.clear();

//					// TODO: energy barrier estimation // extra script???
//					if( out_transition_file != "")
//					{
//						out.open( out_transition_file.c_str());
//						if(!out){ std::cerr << "ERROR: not opened transition energy file: <" << out_transition_file << ">" << std::endl; }
//
//						for( unsigned int i = 0; i < energies.size(); ++i)
//						{ out << i << "\t" << energies[i] << "\t" << names[i] << std::endl;}
//
//						out.close(); out.clear();
//					}


				// collect PDBs for creating nice movie
				out.open( ( outdir + "/path_" + std::to_string( i1) + "_" + std::to_string( i2) + ".pdb").c_str() );
				if(!out){ std::cerr << "ERROR: not opened pdb output file: " << std::endl; }

				// pdb with all models for visualization
				for( unsigned int i = 0; i < names.size(); ++i)
				{
					int id = names.size() - i - 1;
					out << "MODEL " << i << std::endl;
					out << "HEADER" << std::endl;
					out << names[id] << "\t" << energies[id] << std::endl;
					WritePDB( out, names[id], min_energy, energies[id], max_energy );
					out << "ENDMDL" << std::endl;
				}
				out.close();
				out.clear();
			}
			++i1;
			sub_sorted.clear();
			gnu.close();
			gnu.clear();
		}  // primary score improved

		// calc sum of energies of all nodes in path
		energies.clear();
		node = last_node;
		for( std::vector<int>::const_iterator jtr = itr->second.begin(); jtr != itr->second.end(); ++jtr)
		{
			energies.push_back( node->GetEnergy() );
			edge = node->GetParentEdges()[*jtr];
			node = edge.GetNode();
		}
		energies.push_back( node->GetEnergy() );

		// secondary score calculated here:
//		float score = 0.0;
//		for( auto& e : energies)
//		{ score += e;}
		//
		float score = energies.size() * min_max.first;
		for( auto& e : energies)
		{ score += e;}


		sub_sorted.insert( std::make_pair( score, itr->second) );
	}  // iterating score_sorted_paths

};
