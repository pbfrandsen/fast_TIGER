/* fast_TIGER (version 1.0), a program for computing TIGER rates
 * (TIGER: Tree Independent Generation of Evolutionary Rates).
 *
 * Copyright (C) January 2015 by Paul Frandsen and Christoph Mayer
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, see <http://www.gnu.org/licenses/>.
 *
 * For any enquiries send an email to
 * Paul Frandsen: paulbfrandsen@gmail.com
 * or
 * Christoph Mayer: c.mayer.zfmk@uni-bonn.de
 *
 * When publishing work that is based on the results of Fast_Tiger please cite:
 * 
 * Frandsen, P.B., Calcott, B., Mayer, C., Lanfear, R., 2015, Automatic selection of
 * partitioning schemes for phylogenetic analyses using iterative k-means clustering
 * of site rates, BMC Evolutionary Biology 15:13.
 *
 */
//
//  site_pattern.h
//  
//
//  Created by Paul Frandsen on 7/8/14.
//
//

#ifndef _site_pattern_h
#define _site_pattern_h

#include <vector>
#include "faststring2.h"
#include "CSplit2.h"
#include <map>
#include "fast-dynamic-bitset/fast-dynamic-bitset.h"

class SitePattern
{
    faststring pattern;
    std::map<char, std::vector<int> > partition_map;
    std::vector<CSplit> returnable_vector_of_bitsets;
    bool invariant;
    bool empty;
    
public:
    // Takes a faststring that conatins the site pattern
    bool is_invariant()
    {
        return invariant;
    }
    bool is_empty()
    {
        return empty;
    }
    
    SitePattern(faststring new_pattern):invariant(false), empty(false)
    {
		
        pattern = new_pattern;
    }
	
    SitePattern():invariant(false), empty(false)
    {
        pattern = "";
    }
	
    void reset_pattern(faststring new_pattern)
    {
        returnable_vector_of_bitsets.clear();
        pattern = new_pattern;
        invariant = false;
        empty = false;
    }
	
    // print to the stream that the user chooses
    void print(std::ostream &os)
    {
        os << pattern;
    }
    
    // create a map containing the nucleotide states as keys with a vector of ints
    std::map<char, std::vector<int> > gen_pattern_partition()
    {
        // Initiate a map to hold the partition splits
        // character state as keys and taxon numbers as values
        std::vector<int> a;
        std::vector<int> c;
        std::vector<int> g;
        std::vector<int> t;
        partition_map['A'] = a;
        partition_map['C'] = c;
        partition_map['G'] = g;
        partition_map['T'] = t;
        
        // figure out the length of the site pattern
        unsigned long length = pattern.length();

        // loop through the pattern and assign partitions to map
        for (unsigned i = 0; i < length; ++i)
        {
            if (pattern[i] == 'A')
                partition_map['A'].push_back(i);
            else if (pattern[i] == 'C')
                partition_map['C'].push_back(i);
            else if (pattern[i] == 'G')
                partition_map['G'].push_back(i);
            else if (pattern[i] == 'T')
                partition_map['T'].push_back(i);
        }
        return partition_map;
    }
	
    std::vector<CSplit> gen_pattern_bitsets()
    {
        // figure out the length of the site pattern
        unsigned long length = pattern.length();
        unsigned i;
        
        // initiate four different CSplits that we'll return in the vector
        CSplit split_A(length);
        CSplit split_C(length);
        CSplit split_G(length);
        CSplit split_T(length);
        
        for (i = 0; i < length; ++i)
        {
            if (pattern[i] == 'A' || pattern[i] == 'a')
            {
                split_A.set(i);
//                std::cout << "set as A" << std::endl;
            }
            else if (pattern[i] == 'C' || pattern[i] == 'c')
            {
                split_C.set(i);
//                std::cout << "set as C" << std::endl;
            }
            else if (pattern[i] == 'G' || pattern[i] == 'g')
            {
                split_G.set(i);
//                std::cout << "set as G" << std::endl;
            }
            else if (pattern[i] == 'T' || pattern[i] == 't')
            {
                split_T.set(i);
//                std::cout << "set as T" << std::endl;
            }
            else
            {
                split_A.set(i);
                split_C.set(i);
                split_G.set(i);
                split_T.set(i);
            }
        }
        
        // Check if site is full of ambig characters, if it is, flag it as empty
        if (split_A.none() & split_C.none() & split_G.none() & split_T.none())
            empty = true;
        
        // Create the flipped splits to check if site is invariant or empty.
        CSplit split_A_flip = split_A;
        split_A_flip.flip();
        CSplit split_C_flip = split_C;
        split_C_flip.flip();
        CSplit split_G_flip = split_G;
        split_G_flip.flip();
        CSplit split_T_flip = split_T;
        split_T_flip.flip();

        // Check to see if site is invariant, if it is, it's pa will always be 1.
        if (split_A_flip.none() || split_C_flip.none() || split_G_flip.none() || split_T_flip.none())
            invariant = true;
        
        
        returnable_vector_of_bitsets.push_back(split_A);
        returnable_vector_of_bitsets.push_back(split_C);
        returnable_vector_of_bitsets.push_back(split_G);
        returnable_vector_of_bitsets.push_back(split_T);
        return returnable_vector_of_bitsets;
    }
    
    std::vector<CSplit> gen_morph_pattern_bitsets()
    {
        // figure out the length of the site pattern
        unsigned long length = pattern.length();
        int i;
        
        // initiate ten different morphological CSplits that we'll return in the vector
        CSplit split_0(length);
        CSplit split_1(length);
        CSplit split_2(length);
        CSplit split_3(length);
        CSplit split_4(length);
        //        CSplit split_5(length);
        //        CSplit split_6(length);
        //        CSplit split_7(length);
        //        CSplit split_8(length);
        //        CSplit split_9(length);
        
        
        for (i = 0; i < length; ++i)
        {
            if (pattern[i] == '0')
            {
                split_0.set(i);
                //                std::cout << "set as A" << std::endl;
            }
            else if (pattern[i] == '1')
            {
                split_1.set(i);
                //                std::cout << "set as C" << std::endl;
            }
            else if (pattern[i] == '2')
            {
                split_2.set(i);
                //                std::cout << "set as G" << std::endl;
            }
            else if (pattern[i] == '3')
            {
                split_3.set(i);
                //                std::cout << "set as T" << std::endl;
            }
            else if (pattern[i] == '4')
            {
                split_4.set(i);
                //                std::cout << "set as T" << std::endl;
            }
            //            else if (pattern[i] == '5')
            //            {
            //                split_5.set(i);
            //                //                std::cout << "set as T" << std::endl;
            //            }
            //            else if (pattern[i] == '6')
            //            {
            //                split_6.set(i);
            //                //                std::cout << "set as T" << std::endl;
            //            }
            //            else if (pattern[i] == '7')
            //            {
            //                split_7.set(i);
            //                //                std::cout << "set as T" << std::endl;
            //            }
            //            else if (pattern[i] == '8')
            //            {
            //                split_8.set(i);
            //                //                std::cout << "set as T" << std::endl;
            //            }
            //            else if (pattern[i] == '9')
            //            {
            //                split_9.set(i);
            //                //                std::cout << "set as T" << std::endl;
            //            }
            else if (pattern[i] == '-')
            {
                split_0.set(i);
                split_1.set(i);
                split_2.set(i);
                split_3.set(i);
                split_4.set(i);
            }
            else if (pattern[i] == 'N')
            {
                split_0.set(i);
                split_1.set(i);
                split_2.set(i);
                split_3.set(i);
                split_4.set(i);
            }
            else if (pattern[i] == 'n')
            {
                split_0.set(i);
                split_1.set(i);
                split_2.set(i);
                split_3.set(i);
                split_4.set(i);
            }
            else if (pattern[i] == '?')
            {
                split_0.set(i);
                split_1.set(i);
                split_2.set(i);
                split_3.set(i);
                split_4.set(i);
            }
        }
        
        // Check if site is full of ambig characters, if it is, flag it as empty
        if (split_0.none() & split_1.none() & split_2.none() & split_3.none() & split_4.none())
            empty = true;
        
        // Create the flipped splits to check if site is invariant or empty.
        CSplit split_0_flip = split_0;
        split_0_flip.flip();
        CSplit split_1_flip = split_1;
        split_1_flip.flip();
        CSplit split_2_flip = split_2;
        split_2_flip.flip();
        CSplit split_3_flip = split_3;
        split_3_flip.flip();
        CSplit split_4_flip = split_4;
        split_4_flip.flip();
        
        // Check to see if site is invariant, if it is, it's pa will always be 1.
        if (split_0_flip.none() || split_1_flip.none() || split_2_flip.none() || split_3_flip.none() || split_4_flip.none())
            invariant = true;
        
        
        returnable_vector_of_bitsets.push_back(split_0);
        returnable_vector_of_bitsets.push_back(split_1);
        returnable_vector_of_bitsets.push_back(split_2);
        returnable_vector_of_bitsets.push_back(split_3);
        returnable_vector_of_bitsets.push_back(split_4);
        return returnable_vector_of_bitsets;
    }
    
    std::map<char, std::vector<int> > give_pattern_partition()
    {
        return partition_map;
    }
	
    std::vector<CSplit>& give_pattern_bitsets()
    {
        return returnable_vector_of_bitsets;
    }
    
    // Now we need to have a function that can calculate and return the a(x, P(i))
    // between this site pattern and another
    double calculate_pa(SitePattern & j)
    {
        double axpi = 0;
        double pa;
        int num_parts = 4;
        std::map<char, std::vector<int> > partition_to_compare = j.give_pattern_partition();
        // do a whole bunch of set comparisons, to see if any of the sets in the site pattern object are subsets of the
        // sets in the site pattern reference that is fed into the function. I use the "include" function to do this.
        // If the set is empty, then I move onto the next set.

        
        
        
        
        if (partition_to_compare['A'].empty())
            num_parts -= 1;
        else
        {
//            std::cout << "Trying A..." << std::endl;
            if ((std::includes(partition_map['A'].begin(), partition_map['A'].end(), partition_to_compare['A'].begin(), partition_to_compare['A'].end()))
                && (!partition_map['A'].empty() ))
            {
                axpi += 1;
//                std::cout << "A to A." << std::endl;
            }

            else if ((std::includes(partition_map['C'].begin(), partition_map['C'].end(), partition_to_compare['A'].begin(), partition_to_compare['A'].end()))
                     && (!partition_map['C'].empty() ))
            {
                axpi += 1;
//                std::cout << "A to C." << std::endl;
            }

            else if ((std::includes(partition_map['G'].begin(), partition_map['G'].end(), partition_to_compare['A'].begin(), partition_to_compare['A'].end()))
                     && (!partition_map['G'].empty() ))
            {
                axpi += 1;
//                std::cout << "A to G." << std::endl;
            }

            else if ((std::includes(partition_map['T'].begin(), partition_map['T'].end(), partition_to_compare['A'].begin(), partition_to_compare['A'].end()))
                     && (!partition_map['T'].empty() ))
            {
                axpi += 1;
//                std::cout << "A to T." << std::endl;
            }

        }
        if (partition_to_compare['C'].empty())
            num_parts -= 1;
        else
        {
//            std::cout << "Trying C..." << std::endl;
            if ((std::includes(partition_map['A'].begin(), partition_map['A'].end(), partition_to_compare['C'].begin(), partition_to_compare['C'].end()))
                && (!partition_map['A'].empty() ))
            {
                axpi += 1;
//                std::cout << "C to A." << std::endl;
            }

            else if ((std::includes(partition_map['C'].begin(), partition_map['C'].end(), partition_to_compare['C'].begin(), partition_to_compare['C'].end()))
                     && (!partition_map['C'].empty() ))
            {
                axpi += 1;
//                std::cout << "C to C." << std::endl;
            }

            else if ((std::includes(partition_map['G'].begin(), partition_map['G'].end(), partition_to_compare['C'].begin(), partition_to_compare['C'].end()))
                     && (!partition_map['G'].empty() ))
            {
                axpi += 1;
//                std::cout << "C to G." << std::endl;
            }

            else if ((std::includes(partition_map['T'].begin(), partition_map['T'].end(), partition_to_compare['C'].begin(), partition_to_compare['C'].end()))
                     && (!partition_map['T'].empty() ))
            {
                axpi += 1;
//                std::cout << "C to T." << std::endl;
            }

        }
        if (partition_to_compare['G'].empty())
            num_parts -= 1;
        else
        {
//            std::cout << "Trying G..." << std::endl;
            if ((std::includes(partition_map['A'].begin(), partition_map['A'].end(), partition_to_compare['G'].begin(), partition_to_compare['G'].end()))
                 && (!partition_map['A'].empty() ))
            {
                axpi += 1;
//                std::cout << "G to A." << std::endl;
            }

            else if ((std::includes(partition_map['C'].begin(), partition_map['C'].end(), partition_to_compare['G'].begin(), partition_to_compare['G'].end()))
                 && (!partition_map['C'].empty() ))
            {
                axpi += 1;
//                std::cout << "G to C." << std::endl;
            }

            else if ((std::includes(partition_map['G'].begin(), partition_map['G'].end(), partition_to_compare['G'].begin(), partition_to_compare['G'].end()))
                 && (!partition_map['G'].empty() ))
            {
                axpi += 1;
//                std::cout << "G to G." << std::endl;
            }

            else if ((std::includes(partition_map['T'].begin(), partition_map['T'].end(), partition_to_compare['G'].begin(), partition_to_compare['G'].end()))
                && (!partition_map['T'].empty() ))
            {
                axpi += 1;
//                std::cout << "G to T." << std::endl;
            }
        }
        if (partition_to_compare['T'].empty())
            num_parts -= 1;
        else
        {
//            std::cout << "Trying T..." << std::endl;
            if ((std::includes(partition_map['A'].begin(), partition_map['A'].end(), partition_to_compare['T'].begin(), partition_to_compare['T'].end()))
                && (!partition_map['A'].empty() ))
            {
                axpi += 1;
//                std::cout << "T to A." << std::endl;
            }

            else if ((std::includes(partition_map['C'].begin(), partition_map['C'].end(), partition_to_compare['T'].begin(), partition_to_compare['T'].end()))
                     && (!partition_map['C'].empty() ))
            {
                axpi += 1;
//                std::cout << "T to C." << std::endl;
            }

            else if ((std::includes(partition_map['G'].begin(), partition_map['G'].end(), partition_to_compare['T'].begin(), partition_to_compare['T'].end()))
                     && (!partition_map['G'].empty() ))
            {
                axpi += 1;
//                std::cout << "T to G." << std::endl;
            }

            else if ((std::includes(partition_map['T'].begin(), partition_map['T'].end(), partition_to_compare['T'].begin(), partition_to_compare['T'].end()))
                     && (!partition_map['T'].empty() ))
            {
                axpi += 1;
//                std::cout << "T to T." << std::endl;
            }

        }
//        std::cout << axpi << " " << num_parts << std::endl;
        pa = axpi/num_parts;
        return pa;
    }
	
    double calculate_pa_2(SitePattern & j)
    {
        double axpi = 0;
        double pa;
        int num_parts = 4;
        
//        partition_to_compare[0].print(std::cout);
//        std::cout << std::endl;
//        std::cout << !partition_to_compare[0].any() << std::endl;
//        std::cout << partition_to_compare[0].none() << std::endl;
//        std::cout << partition_to_compare[0].any() << std::endl;
//        std::cout << !partition_to_compare[0].none() << std::endl;
//        std::cout << std::endl;
        if (j.is_empty())
            return 1;
        else if (invariant)
            return 1;
        else
        {
            // now do a bunch of comparisons. Yeehaw.
            if (j.give_pattern_bitsets()[0].none())
            {
                num_parts -= 1;
    //            std::cout << "Taking away one." << std::endl;
            }
            else
            {
                if ((returnable_vector_of_bitsets[0].any()) && (returnable_vector_of_bitsets[0].is_parameter_a_subset_of_this(j.give_pattern_bitsets()[0])))
                {
                    axpi += 1;
    //                std::cout << "A to A." << std::endl;
                }
                else if ((returnable_vector_of_bitsets[1].any()) && (returnable_vector_of_bitsets[1].is_parameter_a_subset_of_this(j.give_pattern_bitsets()[0])))
                {
                    axpi += 1;
    //                std::cout << "A to C." << std::endl;
                }
                else if ((returnable_vector_of_bitsets[2].any()) && (returnable_vector_of_bitsets[2].is_parameter_a_subset_of_this(j.give_pattern_bitsets()[0])))
                {
                    axpi += 1;
    //                std::cout << "A to G." << std::endl;
                }
                else if ((returnable_vector_of_bitsets[3].any()) && (returnable_vector_of_bitsets[3].is_parameter_a_subset_of_this(j.give_pattern_bitsets()[0])))
                {
                    axpi += 1;
    //                std::cout << "A to T." << std::endl;
                }
            }
            if (j.give_pattern_bitsets()[1].none())
            {
                num_parts -= 1;
    //            std::cout << "Taking away one." << std::endl;
            }
            else
            {
                if ((returnable_vector_of_bitsets[0].any()) && (returnable_vector_of_bitsets[0].is_parameter_a_subset_of_this(j.give_pattern_bitsets()[1])))
                {
                    axpi += 1;
    //                std::cout << "C to A." << std::endl;
                }
                else if ((returnable_vector_of_bitsets[1].any()) && (returnable_vector_of_bitsets[1].is_parameter_a_subset_of_this(j.give_pattern_bitsets()[1])))
                {
                    axpi += 1;
    //                std::cout << "C to C." << std::endl;
                }
                else if ((returnable_vector_of_bitsets[2].any()) && (returnable_vector_of_bitsets[2].is_parameter_a_subset_of_this(j.give_pattern_bitsets()[1])))
                {
                    axpi += 1;
    //                std::cout << "C to G." << std::endl;
                }
                else if ((returnable_vector_of_bitsets[3].any()) && (returnable_vector_of_bitsets[3].is_parameter_a_subset_of_this(j.give_pattern_bitsets()[1])))
                {
                    axpi += 1;
    //                std::cout << "C to T." << std::endl;
                }
            }
            if (j.give_pattern_bitsets()[2].none())
            {
                num_parts -= 1;
    //            std::cout << "Taking away one." << std::endl;
            }
            else
            {
                if ((returnable_vector_of_bitsets[0].any()) && (returnable_vector_of_bitsets[0].is_parameter_a_subset_of_this(j.give_pattern_bitsets()[2])))
                {
                    axpi += 1;
    //                std::cout << "G to A." << std::endl;
                }
                else if ((returnable_vector_of_bitsets[1].any()) && (returnable_vector_of_bitsets[1].is_parameter_a_subset_of_this(j.give_pattern_bitsets()[2])))
                {
                    axpi += 1;
    //                std::cout << "G to C." << std::endl;
                }
                else if ((returnable_vector_of_bitsets[2].any()) && (returnable_vector_of_bitsets[2].is_parameter_a_subset_of_this(j.give_pattern_bitsets()[2])))
                {
                    axpi += 1;
    //                std::cout << "G to G." << std::endl;
                }
                else if ((returnable_vector_of_bitsets[3].any()) && (returnable_vector_of_bitsets[3].is_parameter_a_subset_of_this(j.give_pattern_bitsets()[2])))
                {
                    axpi += 1;
    //                std::cout << "G to T." << std::endl;
                }
            }
            if (j.give_pattern_bitsets()[3].none())
            {
                num_parts -= 1;
    //            std::cout << "Taking away one." << std::endl;
            }
            else
            {
                if ((returnable_vector_of_bitsets[0].any()) && (returnable_vector_of_bitsets[0].is_parameter_a_subset_of_this(j.give_pattern_bitsets()[3])))
                {
                    axpi += 1;
    //                std::cout << "T to A." << std::endl;
                }
                else if ((returnable_vector_of_bitsets[1].any()) && (returnable_vector_of_bitsets[1].is_parameter_a_subset_of_this(j.give_pattern_bitsets()[3])))
                {
                    axpi += 1;
    //                std::cout << "T to C." << std::endl;
                }
                else if ((returnable_vector_of_bitsets[2].any()) && (returnable_vector_of_bitsets[2].is_parameter_a_subset_of_this(j.give_pattern_bitsets()[3])))
                {
                    axpi += 1;
    //                std::cout << "T to G." << std::endl;
                }
                else if ((returnable_vector_of_bitsets[3].any()) && (returnable_vector_of_bitsets[3].is_parameter_a_subset_of_this(j.give_pattern_bitsets()[3])))
                {
                    axpi += 1;
    //                std::cout << "T to T." << std::endl;
                }
            }
            pa = axpi/num_parts;
            return pa;
        }
    }
};

#endif
