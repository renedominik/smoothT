/*
 * string_functions.h
 *
 *  Created on: Oct 7, 2019
 *      Author: Rene Staritzbichler
 *  License: CC BY 4.0
 */

#ifndef STRING_FUNCTIONS_H_
#define STRING_FUNCTIONS_H_


std::string
Upper( std::string  STR)
{
  std::transform( STR.begin(), STR.end(), STR.begin(), ::toupper);
  return STR;
}

std::map<std::string, char>
aa3to1 =
{
{ Upper("Gly"), 'G'},
{ Upper("Ala"), 'A'},
{ Upper("Leu"), 'L'},
{ Upper("Met"), 'M'},
{ Upper("Phe"), 'F'},
{ Upper("Trp"), 'W'},
{ Upper("Lys"), 'K'},
{ Upper("Gln"), 'Q'},
{ Upper("Glu"), 'E'},
{ Upper("Ser"), 'S'},
{ Upper("Pro"), 'P'},
{ Upper("Val"), 'V'},
{ Upper("Ile"), 'I'},
{ Upper("Cys"), 'C'},
{ Upper("Tyr"), 'Y'},
{ Upper("His"), 'H'},
{ Upper("Arg"), 'R'},
{ Upper("Asn"), 'N'},
{ Upper("Asp"), 'D'},
{ Upper("Thr"), 'T'}};


std::string &
Strip( std::string &STR)
{
  STR.erase( std::remove_if( STR.begin(), STR.end(), isspace), STR.end());
  return STR;
}


template <class Container>
void Split(const std::string& str, Container& cont)
{
    std::istringstream iss(str);
    std::copy(std::istream_iterator<std::string>(iss),
         std::istream_iterator<std::string>(),
         std::back_inserter(cont));

}



#endif /* STRING_FUNCTIONS_H_ */
