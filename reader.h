/*
 * reader.h
 *
 *   read and write
 *   pdb and alignments
 *
 *  Created on: Oct 7, 2019
 *      Author,Copyright: Rene Staritzbichler
 */

#ifndef READER_H_
#define READER_H_


///////////      ALIGNMENT    ///////////////////////////////

void
InvertGapLogic( std::string &S1, std::string &S2)
{
  std::string
    s1, s2;

  for( std::string::iterator i1 = S1.begin(), i2 = S2.begin(); i1 != S1.end(); ++i1, ++i2)
    {
      if( *i1 == '-'){ s2.push_back( '-');}
      else if( *i2 == '-'){ s1.push_back( '-');}
      else{ s1.push_back(*i1); s2.push_back(*i2);}
    }

  S1 = s1;
  S2 = s2;
}


std::vector< std::map<char,std::string> >
ReadAlignments( const std::string &FILE){
    std::cout << __FUNCTION__ << std::endl;
    std::ifstream in( FILE);
    if( !in){ std::cerr << "ERROR: not opened: " << FILE << std::endl; exit(1);}
    std::vector< std::map<char,std::string> >
    	alignments(2);
    char chain; std::string seq1, seq2;
    int nr;
    in >> nr;
    for( int i = 0; i < nr; ++i){
        in >> chain >> seq1;
        in >> chain >> seq2;
        std::cout << chain << std::endl;
        std::cout << seq1 << std::endl;
        std::cout << seq2 << std::endl;
        InvertGapLogic( seq1, seq2);
        alignments[0][chain] = seq1;
        alignments[1][chain] = seq2;
        std::cout << chain << std::endl;
        std::cout << seq1 << std::endl;
        std::cout << seq2 << std::endl;
    }
    return alignments;
}


/////////////////////////      PDB      ///////////////////////////////////////





float X( const std::string &LINE)
{ return std::atof( LINE.substr(30,8).c_str());}
float Y( const std::string &LINE)
{ return std::atof( LINE.substr(38,8).c_str());}
float Z( const std::string &LINE)
{ return std::atof( LINE.substr(46,8).c_str());}
std::string AtomName( const std::string &LINE)
{
  std::string str = LINE.substr( 12,4);
  return Strip( str);
}
std::string Resname( const std::string &LINE)
{
  std::string str = LINE.substr( 17,4);
  return Strip( str);
}
char Chain( const std::string &LINE)
{ return LINE[21];}
int Resid( const std::string &LINE)
{ return std::atoi( LINE.substr(22,4).c_str());}

std::string &
WriteBFactor( std::string &LINE, float X)
{
  std::stringstream stream;
  stream << std::fixed << std::setprecision(2) << std::setw(6) << std::right << X;
  LINE = LINE.substr(0,60) + stream.str() + LINE.substr( 66, LINE.size() - 66 );
  return LINE;
}


std::string &
WritePos( std::string &LINE, std::vector<float> POS)
{
  std::stringstream stream;
  stream << std::fixed << std::setprecision(2) << std::setw(8) << POS[0];  // std::left
  stream << std::fixed << std::setprecision(2) << std::setw(8) << POS[1];
  stream << std::fixed << std::setprecision(2) << std::setw(8) << POS[2];
  LINE = LINE.substr(0,30) + stream.str() + LINE.substr( 54, LINE.size() - 54 );
  return LINE;
}




bool CheckPDBLine( const std::string &LINE, const std::vector< std::string> &ATOM_TYPES, const std::vector<char> &CHAINS, const std::map<char,std::string>& ALIGNMENTS, int & count )
{
  static int prev_resid = -1000;
  static char prev_chain = ' ';

  if( LINE.size() < 5 || ( LINE.substr(0,4) != "ATOM" && LINE.substr(0,6) != "HETATM"))
    { return false;}

  //  std::cout << __FUNCTION__ << std::endl;
  //  std::cout << LINE << std::endl;
  int resid;
  auto chain = Chain(LINE);

  if( ( ATOM_TYPES.size() > 0 && std::find( ATOM_TYPES.begin() , ATOM_TYPES.end(), AtomName( LINE) ) == ATOM_TYPES.end() ) ||
      ( CHAINS.size() > 0 && std::find( CHAINS.begin(), CHAINS.end(), chain )  == CHAINS.end() ) )
    { return false;}

  //  std::cout << "accepted" << std::endl;

  if( ALIGNMENTS.size() > 0){
    //std::cout << "why?" << std::endl;
    if( chain != prev_chain)
      {
	count = -1;
	prev_chain = chain;
      }

    resid = Resid( LINE);

    if( resid != prev_resid){
      ++count;
      prev_resid = resid;
    }

    auto itr = ALIGNMENTS.find(chain);

    //    std::cout << "res: " << itr->second[count] << std::endl;

    if( itr == ALIGNMENTS.end() || itr->second[count] == '-'){ return false;}

    if( aa3to1[ Resname(LINE) ] != itr->second[count] )
      {
	std::cout << LINE << std::endl;

	std::cout << "ERROR: " << itr->second[count] << " vs " << aa3to1[ Resname(LINE) ] << std::endl;
      }
  }
  //std::cout << "accepted" << std::endl;
  return true;
}



std::pair< float, std::vector< std::vector<float> > >
ReadPDBPositionsAndEnergy( const std::string &FILE, const std::string &ENERGY_IDENTIFIER, const std::vector< std::string> &ATOM_TYPES, const std::vector<char> &CHAINS, const std::map<char,std::string>& ALIGNMENTS )
{
    //std::cout << __FUNCTION__ << std::endl;

    std::string
		file;
    float
		energy = -1e12;
    if( ENERGY_IDENTIFIER == "")
    {
    	//std::cout << "read energies from text file" << std::endl;
        std::ifstream
    		in( FILE);
        if( !in)
        { std::cerr << __FUNCTION__ << " ERROR: File not open: " << FILE << std::endl; exit(1);}
        in >> file >> energy;
        in.close(); in.clear();
    }
    else
    {
    	//std::cout << "read energies from pdb file" << std::endl;
    	file = FILE;
    }


    std::ifstream
		in( file);
    if( !in)
    { std::cerr << __FUNCTION__ << " ERROR: File not open: " << file << std::endl; exit(1);}

    std::string
		line;
    std::vector< std::vector< float > >
    	pos;
    int
		count = -1;


    while( std::getline(in, line ) )
    {
    	if( CheckPDBLine( line , ATOM_TYPES, CHAINS, ALIGNMENTS,count ))
		{
		  std::vector< float> p(3);
		  p[0] = X(line);
		  p[1] = Y(line);
		  p[2] = Z(line);
		  pos.push_back(p);
		}
        else if( ENERGY_IDENTIFIER != "" && line.find( ENERGY_IDENTIFIER) != std::string::npos)
		{
		  //	  std::cout << line << std::endl;
		  std::vector<std::string> v;
		  Split( line, v);
		  energy = std::atof( v.back().c_str() );
		  //	  std::cout << energy << std::endl;
		}
    }
    in.close();
    in.clear();
    return std::make_pair( energy, pos );
}



void WritePDB( std::ofstream &OUT, const std::string &NAME, const float &MIN, const float &ENERGY, const float &MAX)
{
	//std::cout << __FUNCTION__ << std::endl;
    std::ifstream
		in( NAME);
    if( !in)
    { std::cerr << __FUNCTION__ <<  " ERROR: File not open: " << NAME << std::endl; exit(1);}

    std::string
		line;
    int
		//count = -1,
		loc = 0;

    while( std::getline(in, line ) ){
    	if( line.size() < 3){ continue;}
    	else if( line.substr(0,4) == "ATOM" || line.substr(0,6) == "HETATM")
    	{
    		//	WritePos( line, POS[count] );
    		if( loc == 0){WriteBFactor( line, MIN);}
    		else if( loc == 1){WriteBFactor( line, MAX);}   // DOES NOT MAKE MUCH SENSE !!!!
    		else{WriteBFactor( line, ENERGY);}
    		OUT << line << std::endl;
			++loc;
			//	++count;
    	}
    	else if(  line.substr(0,3) == "END")
    	{ return;}
    }
    in.close();
    in.clear();
}

double SinoidTransition( const double &CURRENT, const double &XMIN, const double &XMAX, const double &YMIN, const double &YMAX)
{
	if( CURRENT < XMIN)
	{ return YMIN;}
	if( CURRENT > XMAX)
	{ return YMAX;}
	double
		x = M_PI * ( (CURRENT-XMIN)/(XMAX-XMIN) - 0.5 ),
		delta_y = YMAX - YMIN,
		y_shift = YMIN + 0.5 * delta_y,
		y = 0.5 * delta_y * sin(x) + y_shift;

	return y;
}

void EnergyTransition( std::ostream &OUT, const double & RMSD, const double &SUM, const double &OFFSET, const double &ENERGY_1, const double &ENERGY_2, int NR_STEPS , int &COUNT, const double &TRANS_RATIO, const double &DIST_TO_FREE)
{
	float
		delta = RMSD / float(NR_STEPS + 1),
		sum = SUM,
		energy,
		local,
		x;

	for( int i = 1; i <= NR_STEPS; ++i)
	{
		++COUNT;
		sum += delta;
		x = i * delta;
		if( x <= 0.5 * RMSD * ( 1 - TRANS_RATIO) )
		{
			local = x;
			energy = SinoidTransition( x, 0.0, DIST_TO_FREE, ENERGY_1, OFFSET );
		}
		else if( x > 0.5 * RMSD * (1 + TRANS_RATIO) )
		{
			local = RMSD-x;
			energy = SinoidTransition( RMSD - x , 0.0 , DIST_TO_FREE , ENERGY_2 , OFFSET );
		}
		else // transition
		{
			local = ( x / RMSD - 0.5 * (1 - TRANS_RATIO) ) / TRANS_RATIO;
			energy = (1-local) * SinoidTransition( x, 0.0, DIST_TO_FREE, ENERGY_1, OFFSET )
					+ local * SinoidTransition( RMSD - x , 0.0 , DIST_TO_FREE , ENERGY_2 , OFFSET );
		}
		OUT << COUNT << "\t" << local << "\t" << sum << "\t" << energy	<< "\ttrans_" + std::to_string(i) << std::endl; // "\t" << ENERGY_1 << "\t" << ENERGY_2 << "\t" << OFFSET << std::endl;
	}
}





#endif /* READER_H_ */
