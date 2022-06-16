/*
 * path.h
 *
 *  Created on: Oct 7, 2019
 *      Author: Rene Staritzbichler
 */

#ifndef PATH_H_
#define PATH_H_




// forward declaration
class Node;


class Edge
{
private:
	std::shared_ptr<Node> m_Node;
	float                 m_Dist;
public:
	Edge()
	: m_Node(), m_Dist(){}

	Edge( const std::shared_ptr<Node> &NODE, const float & DIST)
	: m_Node( NODE), m_Dist( DIST){}

	const float &
	GetDistance() const {return m_Dist;}

	const std::shared_ptr<Node> &
	GetNode() const
	{ return m_Node;}

};


class Node
{
private:
  std::string m_Name;
  float m_Energy;
  std::vector< std::vector< float > > m_Pos; // list of atom positions
  // list of atoms containing a list of x,y,z coordinates
  /* python 
  m_Pos = []
  with open( filename) as r:
     for l in r:
         if "ATOM" = l[:4]:
	     x = float( l[46:54]
	     y = ...
	     z = ....
	     m_Pos.append( [x,y,z] )
  */  
  std::vector< Edge> m_Parents;

public:
  // constructor
	Node( const std::string &NAME, const std::string &ENERGY_IDENTIFIER, const std::vector< std::string> &ATOM_TYPES, const std::vector<char> &CHAINS,  const std::map<char,std::string> &ALIGNMENT)
    : m_Name( NAME), m_Energy(), m_Parents()
	{
		//std::cout << __FUNCTION__ << " construct from: <" << NAME << "> <" << ENERGY_IDENTIFIER << ">" << std::endl;
		auto all = ReadPDBPositionsAndEnergy( NAME, ENERGY_IDENTIFIER, ATOM_TYPES, CHAINS, ALIGNMENT);
		m_Energy = all.first;
		m_Pos = all.second;
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


	std::ostream & Write( std::ostream &OUT) const
	{
		OUT << m_Name << std::endl;
		OUT << "energy: "<< m_Energy << std::endl;
		OUT << "nr positions: " << m_Pos.size() << std::endl;
		OUT << "nr edges: " << m_Parents.size() << std::endl;
		OUT << m_Pos[0][0] << "\t" << m_Pos[0][1] << "\t" << m_Pos[0][2] << std::endl;
		OUT << m_Pos.back()[0] << "\t" << m_Pos.back()[1] << "\t" << m_Pos.back()[2] << std::endl;
		return OUT;
	}
};


bool operator == ( const Node &N1, const Node &N2 )
{ return N1.GetName() == N2.GetName();}

std::ostream & operator << ( std::ostream &OUT, const Node &N1)
{
	return N1.Write( OUT);
}




bool Iterate
(
		std::vector< std::shared_ptr< Node > > & LATEST,
		std::vector< std::shared_ptr< Node > > & REMAINING,
		const float & MAX_DIST,
		const std::shared_ptr< Node > & LAST
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

    std::ofstream
		out( "pathway/rmsds" + std::to_string(count) + ".txt");
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

    		out << i << "\t" << j << "\t" << distance << "\t" << (*ltr)->GetName() << "\t" << (*rtr)->GetName() << std::endl;

    		if( distance < MAX_DIST)
    		{
    			(*rtr)->AddParent( *ltr , distance);
    			connected.push_back( *rtr );
    		}
    		//if( distance < 100 ){ std::cout << "d: " << distance << std::endl;}
        }
      }
    out.clear(); out.close();

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
    if( REMAINING.size() == 0 || (REMAINING.size() == 1 && REMAINING[0] == LAST))
    {
    	std::cout << "no remaining nodes, iteration is terminated" << std::endl;
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

    // stop if LAST node is in LATEST
    if( std::find( LATEST.begin(), LATEST.end(), LAST) != LATEST.end())
    {
        std::cout << "hooked up last node (first node), having " << LAST->GetParentEdges().size() << " parent edges" << std::endl;
//        std::cout << REMAINING.size() << " remaining nodes (these nodes are not linked to graph)." << std::endl;
	// return false;  //////    VERSION 1.0  === TODO: block and find better strategy! ////////////
    }

    // continue
    return true;
}


// PATH here is more a navigation instruction, which edge to follow in a specific node
void Backtrace(  std::multimap< float, std::vector<int> > & POOL, std::vector<int> PATH, std::pair<float,float> & MIN_MAX_ENERGY, float PATH_MAX , const std::shared_ptr< Node> & NODE, std::shared_ptr< Node> &FIRST_NODE)
{
//    std::cout << __FUNCTION__ << " " << POOL.size() << " " << PATH.size() << " " << MIN_MAX_ENERGY.first << " " << MIN_MAX_ENERGY.second << " " << NODE->GetEnergy() << std::endl;
    MIN_MAX_ENERGY.first  = std::min( MIN_MAX_ENERGY.first , NODE->GetEnergy());
    MIN_MAX_ENERGY.second = std::max( MIN_MAX_ENERGY.second, NODE->GetEnergy());
    PATH_MAX = std::max( PATH_MAX, NODE->GetEnergy());
    // EITHER YOU HAVE A COMPLETE PATH
    if( *NODE == *FIRST_NODE){      // ptr vs obj
//        std::cout << __FUNCTION__ << " connected to first node: " << std::endl;
//        std::cout << NODE->GetName() << " ";
//        std::cout << FIRST_NODE->GetName() << std::endl;
//        std::copy( PATH.begin(), PATH.end(), std::ostream_iterator<int>( std::cout , " ") );
//        std::cout << std::endl;
        POOL.insert( std::make_pair( PATH_MAX , PATH));

		if(POOL.size() > 25000)
		{
			auto b = POOL.begin();
			for( int i = 0; i < 20000; ++i)
			{ ++b;}
			POOL =  std::multimap< float, std::vector<int> >( POOL.begin(), POOL.begin() + 20000);
		}
	
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
        Backtrace( POOL, path, MIN_MAX_ENERGY, PATH_MAX, edges[i].GetNode(), FIRST_NODE);  // recursive call of itself
    }
}








#endif /* PATH_H_ */
