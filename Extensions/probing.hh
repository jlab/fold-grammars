#ifndef PROBING_HH
#define PROBING_HH

#include <iostream>
#include <fstream>
#include <string>

inline double getSHAPEscore(const TUSubsequence &leftBase) {
	static bool isLoaded = false;
	static std::vector<double> probingData;

	if (!isLoaded) {
		std::string line;
		std::ifstream infile (getProbingDataFilename());
		if (infile.is_open()) {
		    while (getline (infile,line)) {
		    	probingData.push_back(atof(line.c_str()));
		    }
		    infile.close();
		}
		if (probingData.size() < (leftBase.seq->n)) {
			std::cerr << "Warning: chemical probing data file '" << getProbingDataFilename() << "' misses " << (leftBase.seq->n - probingData.size()) << " data-row(s) " << std::endl << "         compared to the number of nucleotides in your input sequence." << std::endl << "         Missing values will be set to 0.0!" << std::endl;
		}
		if (probingData.size() > (leftBase.seq->n)) {
			std::cerr << "Warning: chemical probing data file '" << getProbingDataFilename() << "' contains " << (probingData.size()-leftBase.seq->n) << " more row(s) " << std::endl << "         than there are nucleotides in your input sequence." << std::endl << "         Exceeding data lines will be ignored!" << std::endl;
		}

		isLoaded = true;
	}
	
	double score = 0.0;
	for (unsigned int i = leftBase.i; i < leftBase.j && i < probingData.size(); i++) {
		score += probingData.at(i);
	}
	
	return score;
}

inline List_Ref<std::pair<int, double> >& paretoFilter(List_Ref<std::pair<int, double> >& in)
{
	std::list <std::pair<int, double> > newListIn;

	// Computation of the Pareto front.
	for(List_Ref<std::pair<int, double> >::iterator i = in->begin(); i != in->end(); ++i) {

		bool iDominate = false;
		std::list <std::list<std::pair<int, double> >::iterator > toDelete;

		if(newListIn.empty()) {
			iDominate = true;
		} else {
			for(std::list<std::pair<int, double> >::iterator j = newListIn.begin(); j != newListIn.end(); ++j) {
				// i dominate j.
				if( ( (*i).first == (*j).first && (*i).second == (*j).second ) ||
				    ( (*i).first < (*j).first && (*i).second > (*j).second ) ||
				    ( (*i).first < (*j).first && (*i).second == (*j).second ) ||
				    ( (*i).first == (*j).first && (*i).second > (*j).second )) {
					iDominate = true;

					// Delete j from the List.
					toDelete.push_back(j);
				} else {
					// i and j co-dominate.
					if( ( (*i).first < (*j).first && (*i).second < (*j).second ) ||
					    ( (*i).first > (*j).first && (*i).second > (*j).second ) ) {
						iDominate = true;
					} else // j dominate i.
					{
						iDominate = false;
						break;
					}
				}
			}
		}

		// Add the solution i to the list.
		if(iDominate) {	newListIn.push_back(std::make_pair((*i).first,(*i).second)); }

		// Prune the newListIn of all deleted objects.
		for(std::list <std::list<std::pair<int, double> >::iterator >::iterator it = toDelete.begin(); it != toDelete.end(); ++it) { newListIn.erase((*it)); }
	}

	// Copy the solution Vector in List_Ref list
	List_Ref<std::pair<int, double> > *newList = new List_Ref<std::pair<int, double> >();

	for(std::list<std::pair<int, double> >::iterator i = newListIn.begin(); i != newListIn.end(); ++i) { (*newList)->push_back(*i); }

	return *newList;
}


inline List_Ref<std::pair<std::pair<int, double>, String > >& paretoFilterSmooth(List_Ref<std::pair<std::pair<int, double>, String > >& in)
{
	std::list <std::pair<std::pair<int, double>, String > > newListIn;
	// Computation of the Pareto front.

	for(List_Ref<std::pair<std::pair<int, double>, String > >::iterator i = in->begin(); i != in->end(); ++i) {
		if(newListIn.empty()) {
			newListIn.push_back(*i);
		} else {
			for(std::list<std::pair<std::pair<int, double>, String > >::iterator j = newListIn.begin(); j != newListIn.end(); j++)
			{
				if( (*i).first.first > (*j).first.first && (*i).first.second <= (*j).first.second )
				{
					// i is dominated by j.
					break;
				}
				else
				{
					if( (*i).first.first == (*j).first.first && (*i).first.second <= (*j).first.second )
					{
						// i is dominated by j.
						break;
					}
					else
					{
						if( (*i).first.first > (*j).first.first && (*i).first.second > (*j).first.second )
						{
							// i and j co-dominate. We insert i after j so we do nothing. Test for insertion in case of end of the list.

							std::list<std::pair<std::pair<int, double>, String > >::iterator it_test = j;
							it_test++;
							if( it_test == newListIn.end() )
							{
								newListIn.push_back(std::make_pair(std::make_pair((*i).first.first,(*i).first.second), (*i).second));
								break;
							}
						}
						else
						{
							if( (*i).first.first < (*j).first.first && (*i).first.second < (*j).first.second )
							{
								// i and j co-dominate. We insert i on the list before j.
								newListIn.insert(j, std::make_pair(std::make_pair((*i).first.first,(*i).first.second), (*i).second));
							}
							else
							{
								if(  ( (*i).first.first < (*j).first.first && (*i).first.second >= (*j).first.second ) || ( (*i).first.first == (*j).first.first && (*i).first.second > (*j).first.second ) )
								{
									// i dominate j. We place i before j and we prune the tree.
									newListIn.insert(j, std::make_pair(std::make_pair((*i).first.first,(*i).first.second), (*i).second));

									// We prune the list.
									std::list<std::pair<std::pair<int, double>, String > >::iterator it_del = j;

									bool loopOn = true;
									while( loopOn )
									{
										if( it_del != newListIn.end() )
										{
											if( (*it_del).first.first > (*i).first.first && (*it_del).first.second <= (*i).first.second )
											{
												it_del = newListIn.erase(it_del);
											}
											else
											{
												loopOn = false;
											}
										}
										else
										{
											loopOn = false;
										}
									}
									if(!loopOn)
									{
										break;
									}
								}
								else
								{
									// End of the List case.
									std::cerr<<"Error !!!!!!!!!!"<<std::endl;
									exit(0);
								}
							}
						}
					}
				}
			}
		}
	}

	// Copy the solution Vector in List_Ref list
	List_Ref<std::pair<std::pair<int, double>, String > > *newList = new List_Ref<std::pair<std::pair<int, double>, String > >();

	for(std::list<std::pair<std::pair<int, double>, String > >::iterator i = newListIn.begin(); i != newListIn.end(); ++i) { (*newList)->push_back(*i); }

	return *newList;
}

#endif


