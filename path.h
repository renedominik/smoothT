/*
 * path.h
 *
 *  Created on: Oct 7, 2019
 *      Author: Rene Staritzbichler
 */

#ifndef PATH_H_
#define PATH_H_


// path = graph consists of nodes and edges


// forward declaration
class Node;


class Edge
{
private:
	std::shared_ptr<Node> m_Node;   // pointer to parent
	float                 m_Dist;   // RMSD of current and parent node
public:
	// default constructur
	Edge()
	: m_Node(), m_Dist(-1.0){}

	// construct from member
	Edge( const std::shared_ptr<Node> &NODE, const float & DIST)
	: m_Node( NODE), m_Dist( DIST){}

	// copy constructor
	Edge( const Edge &A)
	: m_Node( A.m_Node), m_Dist( A.m_Dist)
	{}

	// destructor
	~Edge(){}

	const float &
	GetDistance() const
	{return m_Dist;}

	const std::shared_ptr<Node> &
	GetNode() const
	{ return m_Node;}

};


// node contains protein conformations
class Node
{
private:
  std::string m_Name;                             // file name
  float m_Energy;                                 // energy/score describing stored conformation
  std::vector< std::vector< float > > m_Pos;      // list of atom positions  // only information stored from PDB
  std::vector< Edge> m_Parents;                   // list of edges to parent nodes
  float m_Sum;                                    // sum over energies along (sub)path  (initialized undefined in constructor as max float)
  float m_Barrier;                                // max over previous (sub)path
  Edge  m_Best;                  // best parent ever

public:
  // constructor
	Node( const std::string &NAME, const std::string &ENERGY_IDENTIFIER, const std::vector< std::string> &ATOM_TYPES, const std::vector<char> &CHAINS,  const std::map<char,std::string> &ALIGNMENT)
    : m_Name( NAME), m_Energy(), m_Pos(), m_Parents(), m_Sum( std::numeric_limits<float>::max()), m_Barrier(), m_Best()
	{
		//std::cout << __FUNCTION__ << " construct from: <" << NAME << "> <" << ENERGY_IDENTIFIER << ">" << std::endl;
		auto all = ReadPDBPositionsAndEnergy( NAME, ENERGY_IDENTIFIER, ATOM_TYPES, CHAINS, ALIGNMENT);
		m_Energy = all.first;
		m_Pos = all.second;
		m_Barrier = m_Energy;
	}

	~Node(){} // std::cout << __FUNCTION__ << std::endl;}

	void
	AddParent( const std::shared_ptr<Node> &NODE, const float & DIST)
	{
		m_Parents.push_back( Edge( NODE, DIST) );
	}


	const std::vector< Edge> &
	GetParentEdges() const
	{ return m_Parents;}


	const std::vector< std::vector< float > > &
	GetPos() const
	{ return m_Pos;}


	const std::string &
	GetName() const
	{ return m_Name;}

	const float &
	GetEnergy() const
	{ return m_Energy;}

	void
	SetEnergy( const float &VALUE)
	{ m_Energy = VALUE;}

	const Edge &
	GetBest() const
	{ return m_Best;}

	void
	SetBest( const Edge &PTR)
	{  m_Best = PTR;}

	const float &
	GetSum() const
	{ return m_Sum;}

	void
	SetSum( const float &VALUE)
	{  m_Sum = VALUE;}

	const float &
	GetBarrier() const
	{ return m_Barrier;}

	void
	SetBarrier( const float &VALUE)
	{ m_Barrier = VALUE;}



	std::ostream & Write( std::ostream &OUT) const
	{
		OUT << m_Name << std::endl;
		OUT << "energy: "<< m_Energy << std::endl;
		OUT << "nr positions: " << m_Pos.size() << std::endl;
		OUT << "nr edges: " << m_Parents.size() << std::endl;
		//		OUT << m_Pos[0][0] << "\t" << m_Pos[0][1] << "\t" << m_Pos[0][2] << std::endl;
		//		OUT << m_Pos.back()[0] << "\t" << m_Pos.back()[1] << "\t" << m_Pos.back()[2] << std::endl;
		OUT << "barrier: " << m_Barrier << std::endl;
		OUT << "sum: " << m_Sum << std::endl;
		OUT << "best edge distance: " << m_Best.GetNode().GetDistance() << std::endl;
		return OUT;
	}
};


bool operator == ( const Node &N1, const Node &N2 )
{ return N1.GetName() == N2.GetName();}

std::ostream & operator << ( std::ostream &OUT, const Node &N1)
{
	return N1.Write( OUT);
}



// function for graph construction
bool Iterate
(
		std::vector< std::shared_ptr< Node > > & LATEST,
		std::vector< std::shared_ptr< Node > > & REMAINING,
		const float & MAX_DIST,
		const std::shared_ptr< Node > & LAST,
		std::vector< std::vector< std::shared_ptr< Node > > > & GENERATIONS
 )
{
	static int
		count = 0;
    ++count;

//    static int
//		prev_nr_connected = -1;

//    std::cout << __FUNCTION__ << " " << LATEST.size() << " " << REMAINING.size() << std::endl;
    std::list< std::shared_ptr< Node > >
    	connected;
    float
		distance;

//    std::ofstream
//		out( "pathway/rmsds" + std::to_string(count) + ".txt");
    int i = 0, j = 0;
    // find remaining that are connected to latest
    for( std::vector< std::shared_ptr< Node > >::const_iterator rtr = REMAINING.begin(); rtr != REMAINING.end(); ++rtr, ++i)
      {
	j = 0;
    	for( std::vector< std::shared_ptr< Node > >::iterator ltr = LATEST.begin(); ltr != LATEST.end(); ++ltr, ++j)
        {
    		if( *rtr == *ltr && *rtr != LAST )
    		{
    			std::cout << __FUNCTION__ << ": WARNING: nodes are equal!" << std::endl;
    			break;
    		}

    		distance = RMSD( (*rtr)->GetPos(), (*ltr)->GetPos() );

//    		out << i << "\t" << j << "\t" << distance << "\t" << (*ltr)->GetName() << "\t" << (*rtr)->GetName() << std::endl;

    		if( distance < MAX_DIST)
    		{
    			(*rtr)->AddParent( *ltr , distance);  // add pointer to parent, for backtrace construction of pathways
    			connected.push_back( *rtr );
    		}
    		//if( distance < 100 ){ std::cout << "d: " << distance << std::endl;}
        }
      }
//    out.clear(); out.close();

//    std::cout << "nr connected " << connected.size() << std::endl;

    connected.unique();


    // stop if no new connections found
    if( connected.size() == 0 )
    {
    	std::cout << "nothing was found to be within RMSD cutoff, terminate" << std::endl;
    	return false;
    }

    // remove connected from remaining
    for( std::vector< std::shared_ptr< Node > >::iterator itr = REMAINING.begin(); itr != REMAINING.end(); ++itr)
    {
    	if( find( connected.begin(), connected.end(), *itr) != connected.end() && *itr != LAST)
    	{
    		REMAINING.erase( itr);
    		--itr;
    	}
    }

    std::cout << __FUNCTION__ << " " << count << ": unique nr connected " << connected.size() << " remaining: " << REMAINING.size() << std::endl;

    // stop if no nodes are in REMAINING
    if( REMAINING.size() == 0 )
    {
    	std::cout << "no remaining nodes, iteration is terminated" << std::endl;
    	return false;
    }
    else if( REMAINING.size() == 1 && REMAINING[0] == LAST)
    {
    	std::cout << "last node remaining, iteration is terminated" << std::endl;
        return false;
    }

//    if( connected.size() == prev_nr_connected)
//    {
//    	std::cout  << __FUNCTION__<< " no new connections found in this iteration, terminate" << std::endl;
//    	return false;
//    }
//    prev_nr_connected = connected.size();


    // pass newly connected to latest for next iteration
    //std::cout << "node destructors:" << std::endl;
    LATEST.clear(); // probably egal
    //std::cout << "end destruction" << std::endl;
    LATEST = std::vector<std::shared_ptr< Node >  >( connected.begin(), connected.end());

    GENERATIONS.push_back( LATEST);

    // stop if LAST node is in LATEST
    if( std::find( LATEST.begin(), LATEST.end(), LAST) != LATEST.end())
    {
        std::cout << "hooked up last node, having " << LAST->GetParentEdges().size() << " parent edges" << std::endl;
//        std::cout << REMAINING.size() << " remaining nodes (these nodes are not linked to graph)." << std::endl;
	// return false;  //////    VERSION 1.0  === TODO: block and find better strategy! ////////////
    }

    // continue
    return true;
}



void Insert
(
		const std::vector<int> & PATH,
		std::map< float, std::multimap< float, std::vector<int> > > & POOL,
		const std::shared_ptr< Node> &LAST_NODE,
		const float & ZERO_ENERGY,
		const int & MAX_NR_BARRIERS,
		const int & MAX_NR_PATHS
)
{
	std::vector<float>
		rmsds,
		energies;
	std::shared_ptr< Node>
		node = LAST_NODE;
	Edge
		edge;
	for( std::vector<int>::const_iterator itr = PATH.begin(); itr != PATH.end(); ++itr)
	{
		energies.push_back( node->GetEnergy() - ZERO_ENERGY);
		edge = node->GetParentEdges()[*itr];
		rmsds.push_back( edge.GetDistance() );
		node = edge.GetNode();
	}
	energies.push_back( node->GetEnergy() - ZERO_ENERGY);
	float
		barrier = *std::max_element( energies.begin(), energies.end()),
		integral = 0.0;
	for( int i = 0; i < rmsds.size(); ++i)
	{
		integral += 0.5* rmsds[i] * (energies[i] + energies[i+1]);
	}
	auto loc = POOL.find( barrier);
	if( loc  != POOL.end())
	{
		// barrier exists in POOL
		loc->second.insert( std::make_pair( integral, PATH ) );
		if( loc->second.size() > MAX_NR_PATHS + 500 )
		{
			auto x = loc->second.begin();
			for( int i = 0; i < MAX_NR_PATHS; ++i){ ++x;}
			loc->second = std::multimap< float, std::vector<int> >( loc->second.begin(), x);
			std::cout << "trimmed: " << loc->second.size() << std::endl;
		}
	}
	else
	{
		// 	new barrier
		std::multimap< float, std::vector<int> >  m;
		m.insert( std::make_pair( integral, PATH ) );
		POOL[barrier] = m;
		if( POOL.size() > MAX_NR_BARRIERS)
		{
			POOL.erase( std::prev( POOL.end() ) );
		}
	}
	// FILTER POOL

}



void Write
(
		std::map< float, std::multimap< float, std::vector<int> > > & POOL,
		const std::shared_ptr< Node> &LAST_NODE,
		const std::string & OUTDIR
)
{
	int
		c1 = 0;
	for( std::map< float, std::multimap< float, std::vector<int> > >::const_iterator itr = POOL.begin(); itr != POOL.end(); ++itr, ++c1)
	{
		int c2 = 0;
		for( std::multimap< float, std::vector<int> >::const_iterator jtr = itr->second.begin(); jtr != itr->second.end(); ++jtr, ++c2)
		{
			std::ofstream
				out1( OUTDIR + "/energies_" + std::to_string( c1) + "_" + std::to_string(c2) + ".txt"),
				out2( OUTDIR + "/path_" + std::to_string( c1) + "_" + std::to_string(c2) + ".pdb");

			std::vector<float>
				energies,
				rmsds;
			std::vector<std::string>
				names;
			std::shared_ptr< Node>
				node = LAST_NODE;
			Edge
				edge;

			for( std::vector<int>::const_iterator ktr = jtr->second.begin(); ktr != jtr->second.end(); ++ktr)
			{
				energies.push_back( node->GetEnergy());
				names.push_back( node->GetName());
				edge = node->GetParentEdges()[*ktr];
				rmsds.push_back( edge.GetDistance() );
				node = edge.GetNode();
			}
			energies.push_back( node->GetEnergy());
			names.push_back( node->GetName());
			rmsds.push_back( 0.0);

			std::reverse( energies.begin(), energies.end());
			std::reverse( rmsds.begin(), rmsds.end());
			std::reverse( names.begin(), names.end());


			float
				max_energy = *std::max_element( energies.begin(), energies.end()),
				min_energy = *std::min_element( energies.begin(), energies.end()),
				rmsd = 0;

			// pdb with all models for visualization
			for( unsigned int i = 0; i < names.size(); ++i)
			{
				rmsd += rmsds[i];
				out1 << rmsd << "\t" << energies[i] << std::endl;
				out2 << "MODEL " << i + 1 << std::endl;
				out2 << "HEADER" << std::endl;
				out2 << names[i] << "\t" << energies[i] << std::endl;
				WritePDB( out2, names[i], min_energy, energies[i], max_energy );
				out2 << "ENDMDL" << std::endl;
			}
			out1.close(); out1.clear();
			out2.close(); out2.clear();
		}

	}

}



void Write
(
		const std::vector< std::shared_ptr<Node> >  &PATH,
		const float &MIN,
		const float &MAX,
		const std::string & OUTDIR
)
{

	std::ofstream
		out1( OUTDIR + "/energies.txt"),
		out2( OUTDIR + "/path.pdb");


	float
		rmsd,
		min_energy,
		max_energy;

	// pdb with all models for visualization
	for( unsigned int i = 0; i < PATH.size(); ++i)
	{
		if( i == 0){
			rmsd = 0.0;
		}else{
			rmsd += PATH[i]->GetBest().GetDistance();
		}
		out1 << rmsd << "\t" << PATH[i]->GetEnergy() << std::endl;
		out2 << "MODEL " << i + 1 << std::endl;
		out2 << "HEADER" << std::endl;
		out2 << PATH[i]->GetName() << "\t" << PATH[i]->GetEnergy() << std::endl;
		WritePDB( out2, PATH[i]->GetName(), MIN, PATH[i]->GetEnergy(), MAX );
		out2 << "ENDMDL" << std::endl;
	}

	out1.close(); out1.clear();
	out2.close(); out2.clear();
}





void
Backtrace( const std::shared_ptr<Node> &NODE, std::vector< std::shared_ptr<Node> >  &PATH)
{
	PATH.push_back( NODE);
	std::cout << __FUNCTION__ << " " << PATH.size() << "  " << NODE << std::endl;
	exit(1);
	if( NODE->GetBest().GetDistance() >= 0)
	{
		Backtrace( NODE->GetBest().GetNode(), PATH);
	}
}


// PATH here is more a navigation instruction, which edge to follow in a specific node
void Backtrace
(
		std::map< float, std::multimap< float, std::vector<int> > > & POOL,
		std::vector<int> PATH,
		const float & ZERO_ENERGY,
		const std::shared_ptr< Node> & NODE,
		const std::shared_ptr< Node> &FIRST_NODE,
		const std::shared_ptr< Node> &LAST_NODE,
		const int & MAX_NR_BARRIERS,
		const int & MAX_NR_PATHS
)
{
    // EITHER YOU HAVE A COMPLETE PATH
    if( *NODE == *FIRST_NODE){

    	Insert( PATH, POOL, LAST_NODE, ZERO_ENERGY, MAX_NR_BARRIERS, MAX_NR_PATHS);
	
        return; // ENDS FUNCTION
    }

    // OR YOU STILL ARE ON THE WAY TO THE FIRST NODE

    auto edges = NODE->GetParentEdges();

//    std::cout << edges.size() << " edges" << std::endl;
    if( edges.size() == 0)
    {
        std::cout <<  "WARNING: node has no parents!!" << std::endl;
    }

    for( unsigned int i = 0; i < edges.size(); ++i)
    {
        auto path = PATH;
        path.push_back( i );  // extending the path with current edge ID
        // ID is representing the local structure of the graph
        Backtrace( POOL, path, ZERO_ENERGY, edges[i].GetNode(), FIRST_NODE, LAST_NODE, MAX_NR_BARRIERS, MAX_NR_PATHS);  // recursive call of itself
    }
}


void Shift( std::shared_ptr< Node> & NODE, const float &SHIFT)
{
	NODE->SetEnergy( NODE->GetEnergy() - SHIFT);
//	for( auto itr = NODE->GetParentEdges().begin(); itr != NODE->GetParentEdges().end(); ++itr)
//	{
//		auto node = itr->GetNode();
//		Shift( node,  SHIFT);
//	}
}


void GenerationWalk( const std::vector< std::vector< std::shared_ptr< Node > > > &GENERATIONS)
{
	for( auto gen = GENERATIONS.begin(); gen != GENERATIONS.end(); ++gen)
		for( auto node = gen->begin(); node != gen->end(); ++node)
		{
			float
				prev_sum,
				prev_barrier,
				best_tmp_area = std::numeric_limits<float>::max(),
				best_area = std::numeric_limits<float>::max(),
				sum = 0.0,
				prev_area,
				current_area,
				barrier = std::numeric_limits<float>::max();
			std::shared_ptr< Node>
				prev_node,
				prev_tmp_node;
			bool
				first_time = true;

			for( auto edge = (*node)->GetParentEdges().begin(); edge != (*node)->GetParentEdges().end(); ++edge)
			{
				prev_node = edge->GetNode();
			       
				prev_barrier = prev_node->GetBarrier();
				prev_sum = prev_node->GetSum();
				current_area = 0.5 * ( prev_node->GetEnergy() + (*node)->GetEnergy() ) * edge->GetDistance();

				if( barrier > (*node)->GetEnergy() )
				{
					if( prev_barrier < barrier)
					{
						barrier = prev_barrier;
						best_area = current_area;
//						best_prev_sum = prev_sum;
						(*node)->SetBest( *edge);
					}
					else if( prev_barrier == barrier && current_area < best_area)
					{
						best_area = current_area;
//						best_prev_sum = prev_sum;
						(*node)->SetBest( *edge);
					}
				}
				else if( prev_barrier <= (*node)->GetEnergy() )
				{
					if( first_time)
					{
						best_area =  std::numeric_limits<float>::max();
						first_time = false;
					}
					if( current_area < best_area)
					{
						best_area = current_area;
						barrier = prev_barrier;
//						best_prev_sum = prev_sum;
						(*node)->SetBest( *edge);
					}
				}


			}
			// barrier
			if( (*node)->GetBest().GetNode())
			{
				(*node)->SetBarrier( std::max( (*node)->GetEnergy(), barrier ));
				(*node)->SetSum( (*node)->GetBest().GetNode()->GetSum() + best_area);
			}
			else
			{
				std::cerr << "ERROR: Generation Walk did not find a path" << std::endl;
				return;
			}
			// sum
		}
}






#endif /* PATH_H_ */
