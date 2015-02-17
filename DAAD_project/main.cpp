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


#include <iostream>
#include <fstream>
#include "../CSequences2.h"
#include "../faststring2.h"
#include "../site_pattern.h"
#include "../CFile/CFile2_1.h"
#include <map>
#include <vector>
#include <iterator>
#include <cstdlib>
#include <numeric>
//#include "../kmeans.h"
//#include "../sp_container.h"

using namespace std;

CSequences2 read_fasta_file(int type, const char *fn )
{
    CFile file;
    
    cerr << "Reading fasta file:" << fn << endl;
    
    file.ffopen(fn);
    
    if (file.fail())
    {
        cerr << "Could not open specified file: " << fn << endl;
        exit(-1);
    }
    
    CSequences2 *seqs;
    if (type == 1)
        seqs = new CSequences2(CSequence_Mol::dna);
    else
        seqs = new CSequences2(CSequence_Mol::protein);
    
    CSequence_Mol::processing_flag pflag = CSequence_Mol::processing_flag(0);
    seqs->read_from_Fasta_File(file, pflag, 0, -1, true);
    file.ffclose();
    
    cerr << "Found " << seqs->GetTaxaNum() << " taxa." << endl;
    cerr << "Found " << seqs->GetPosNum() << " positions." << endl;
    
    return *seqs;
}

CSequences2 read_phylip_file(int type, const char *fn )
{
    CFile file;
    
    cerr << "Reading phylip file:" << fn <<endl;
    
    file.ffopen(fn);
    
    if (file.fail())
    {
        cerr << "Could not open specified file: " << fn << endl;
        exit(-1);
    }
    
    CSequences2 *seqs;
    if (type == 1)
        seqs = new CSequences2(CSequence_Mol::dna);
    else
        seqs = new CSequences2(CSequence_Mol::protein);
    
    seqs->read_from_Phylip_File(file);
    file.ffclose();
    
    cerr << "Found " << seqs->GetTaxaNum() << " taxa." << endl;
    cerr << "Found " << seqs->GetPosNum() << " positions." << endl;
    
    return *seqs;
}

/* This main function is for when using bitsets */
int main(int argc, const char ** argv)
{
    // open the phylip file
    if (argc != 3)
    {
        cerr << "Usage: " << argv[0] << " data-type phylip-filename" << endl;
        cerr << "Data-type must be either dna or protein" << endl;
        exit(0);
    }
    
    
    faststring filename = argv[2];
    faststring data_type = argv[1];
    
    
    if (data_type == "dna")
    {
        CSequences2 sequences = read_phylip_file(1, filename.c_str());
        unsigned alignment_length = sequences.GetPosNum();
        
        vector<faststring> seq_patterns;
        vector<SitePattern> all_patterns(alignment_length);
        
        seq_patterns.reserve(alignment_length);
        
        // TODO: Use reference for this rather than making a copy
        seq_patterns = sequences.get_pattern_vec();
        
        SitePattern temp(seq_patterns[0]);
        // TODO: Use reference here as well
        vector <CSplit> site_partition = temp.gen_pattern_bitsets();
        
        unsigned i;
        unsigned j;
        int k = 0;
        
        for (i = 0; i < alignment_length; ++i)
        {
            temp.reset_pattern(seq_patterns[i]);
            site_partition = temp.gen_pattern_bitsets();
            // Define a container object first, then use add function to add site patterns. Call the add function and the add function will call the pushback
            all_patterns[i] = (temp);
        }
        

        
        double pa_sum;
        vector<double> site_rates(alignment_length);
        
        // If the alignment is larger than 500,000 sites, this small loop
        // just makes a subset of the alignment patterns for comparison,
        // I have hard coded it to 5...we can also do more or less later.
        // by changing the "subset_size" variable.
        if (alignment_length > 500000)
        {
            unsigned subset_length = 10000;
            unsigned subset_size = alignment_length/subset_length;
            vector<SitePattern> subset_patterns(subset_length);
            // Make the subset of patterns
            for (i = 0; i < subset_length - 1; ++i)
            {
                subset_patterns[i] = all_patterns[k];
                k += subset_size;
            }
            
            for (i = 0; i < alignment_length; ++i)
            {
                if (i%1000 == 0)
                    cerr << i << " of " << alignment_length << " generations complete." << endl;
                pa_sum = 0;
                for (j = 0; j < subset_length - 1; ++j)
                {
                    pa_sum += all_patterns[i].calculate_pa_2(subset_patterns[j]);
                }
                site_rates[i] = (pa_sum/(subset_length - 1));
            }
        }
        // If the alignment is smaller than or equal to 500,000 sites than we do every comparison
        else
        {
            for (i = 0; i < alignment_length; ++i)
            {
                if (i%1000 == 0)
                    cerr << i << " of " << alignment_length << " generations complete." << endl;
                pa_sum = 0;
                for (j = 0; j < alignment_length; ++j)
                {
                    if (i == j)
                        continue;
                    else
                        pa_sum += all_patterns[i].calculate_pa_2(all_patterns[j]);
                }
                site_rates[i] = (pa_sum/(alignment_length - 1));
            }
        }
	cerr << i << " of " << alignment_length << " generations complete." << endl;
        faststring rates_name = filename;
//        int l;
//        for (l = 0; l < 3; ++l)
//        {
//            rates_name.pop_back();
//        }
        rates_name.push_back('_');
        rates_name.push_back('r');
        rates_name.push_back('8');
        rates_name.push_back('s');
        rates_name.push_back('.');
        rates_name.push_back('t');
        rates_name.push_back('x');
        rates_name.push_back('t');
        ofstream rates_file;
        rates_file.open(rates_name.c_str());
        
        copy(site_rates.begin(), site_rates.end(), ostream_iterator<double>(rates_file, "\n"));
//        vector <int> cluster_assigns;
//        cluster_assigns = kmeans(site_rates, 4);
//        copy(cluster_assigns.begin(), cluster_assigns.end(), ostream_iterator<int>(cout, ","));
//        cout << endl;
    }
    else if (data_type == "protein")
    {
        CSequences2 sequences = read_phylip_file(2, filename.c_str());
        cerr << "Sorry this program currently only accepts DNA alignments." << endl;
    }
    else if (data_type == "morphology")
    {
        CFile file;
        
        cerr << "Reading phylip file:" << filename <<endl;
        
        file.ffopen(filename.c_str());
        
        
        if (file.fail())
        {
            cerr << "Could not open specified file: " << filename << endl;
            exit(-1);
        }
        
        unsigned        all_lines_N;
        unsigned        curr_taxon;
        unsigned        posNum;
        unsigned        taxaNum;
        
        char read_mode;
        
        // The allowed symbols from the phylip manual:
        
        
        char ch;
        faststring line;
        std::vector<faststring> fsvec;
        faststring *p;
        
        ch = file.peek();
        
        // If the first char is a white space or a digit, we expect this
        // is the line with the number of taxa and seq. positions.
        if ( isspace(ch) || isdigit(ch) )
        {
            file.getline(line);
            split(fsvec, line);
            
            if (fsvec.size()==2 && fsvec[0].isAnUnsigned() && fsvec[1].isAnUnsigned() )
            {
                taxaNum = fsvec[0].ToUnsigned();
                posNum  = fsvec[1].ToUnsigned();
                
                cout << "Number of taxa     " << taxaNum << endl;
                cout << "Number of residues " << posNum  << endl;
            }
        }
        
        vector<faststring> all_seqs(taxaNum);
        
        vector<faststring*> all_lines;
        
        p = new faststring();
        file.getline(*p);
        p->removeSpacesFront();
        p->removeSpacesBack();
        
        while (!file.fail() )
        {
            all_lines.push_back(p);
            p = new faststring();
            file.getline(*p);
            p->removeSpacesFront();
            p->removeSpacesBack();
        }
        unsigned i;
        all_lines_N = all_lines.size();
        
        //     std::cerr << all_lines[all_lines_N-1]->length() << std::endl;
        
        while (all_lines[all_lines_N-1]->length() == 0)
        {
            --all_lines_N;
            delete all_lines[all_lines_N];
            std::vector<faststring*>::iterator it = all_lines.end();
            --it;
            all_lines.erase(it);
        }
        
        unsigned number_blank_lines=0;
        
        for (i=0; i<all_lines_N; ++i)
        {
            if (all_lines[i]->length() == 0)
            {
                ++number_blank_lines;
            }
        }
        
        if (all_lines_N == 0 || number_blank_lines == all_lines_N) // Nothing to do???
            return 0;
        
        if (number_blank_lines == 0) // sequential format
        {
            read_mode = 's';
            // Example:
            //a          ATGATTCAAC CTCAGACCCT TTTAAATGTA GCAGATAACA GTGGAGCTCG
            //           AAAATTGATG
            //b          ATGATTCAAC CTCAGACCCA TTTAAATGTA GCGGATAACA GCGGGGCTCG
            //           AGAATTGATG
            
            // This can only be interpreted, if the number
            // of taxa and positions has been specified in
            // the first line.
            // What makes it more complicated
            // is the fact, that the phylip format tolerates
            // numbers and other stuff at the beginning of pure sequence lines.
            
            faststring  seq_taxon_name;
            faststring  seq_in_one_line;
            //       unsigned    current_length = 0;
            unsigned    curr_line_num=0;
            
            faststring  *currLine;
            
            
            if (taxaNum == 0 || posNum == 0)
            {
                // delete all faststrings in all_lines
                for (unsigned ind=0; ind < all_lines.size(); ++ind)
                    delete all_lines[ind];
                return -3;
            }
            
            
            // For all taxa
            //
            for (curr_taxon = 0; curr_taxon < taxaNum; ++curr_taxon)
            {
                if (curr_line_num >= all_lines_N)
                {
                    std::cerr << "Parse error while reading the input file:\nUnexpected end of file in line: " << curr_line_num <<  std::endl;
                    std::cerr << "Trying to read " << taxaNum << " taxa but found only " << curr_taxon << "!" <<  std::endl;
                    
                    exit(-45);
                }
                
                currLine       = all_lines[curr_line_num];
                const char* pos = currLine->c_str();
                //	 currLineLen     = currLine->length();
                seq_taxon_name.clear();
                seq_in_one_line.clear();
                
                // Copy the sequence name
                for (;*pos != '\0' && *pos != ' '; ++pos)
                {
                    seq_taxon_name.push_back(*pos);
                }
                
                // Skip all spaces:
                while (*pos == ' ' && *pos != '\0')
                    ++pos;
                
                // Read data after sequence name:
                if (*pos)
                {
                    // Copy the sequence data:
                    {
                        // Find first character that shall be copied.
                        // pos should point to the first char of the sequence.
                        char c;
                        
                        // Skip all non-allowed symbols - since we skiped spaces, there should be none.
                        c = *pos;
                    }
                    // Now pos points to the first char that should be copied:
                    while (*pos)
                    {
                        if (*pos != ' ')
                        {
                            seq_in_one_line.push_back(*pos);
                        }
                        ++pos;
                    }
                } // END if (*pos) AND read data after sequence name.
                
                while (seq_in_one_line.length() < posNum)
                {
                    ++curr_line_num;
                    if (curr_line_num >= all_lines_N)
                    {
                        std::cerr << "Parse error while reading the input file. The file has been interpreted as being in sequential phylip format." << std::endl;
                        std::cerr << "Several problems can trigger this error: (i) Wrong number of residues in sequences compared to the number specified " << std::endl;
                        std::cerr << "in the file header. (ii) Wrong number of taxa (sequences) compared to the number specified in the file header." << std::endl;
                        std::cerr << "(iii) Missing blank line between blocks of data if this file should be in interleaved format." << std::endl;
                        exit(-44);
                    }
                    currLine       = all_lines[curr_line_num];
                    //	   const char* pos = currLine->c_str();
                    {
                        // Find first character that shall be copied.
                        char c;
                        pos = currLine->c_str();
                        c = *pos;
                    }
                    // Now pos points to the first char that should be copied:
                    // Here we are in the sequential section, reading all residues until
                    // we have filled our sequence.
                    while (*pos)
                    {
                        if (*pos != ' ')
                        {
                            seq_in_one_line.push_back(*pos);
                        }
                        ++pos;
                    }
                }
                // If we are here, the sequence should have been copied to
                // seq_in_one_line
                //                cout << "Here is the taxon: " << seq_taxon_name << endl;
                //                cout << "And here is the sequence: " << seq_in_one_line << endl;
                
                all_seqs[curr_taxon] = seq_in_one_line;
                
                // The next line should start with a new sequence:
                ++curr_line_num;
                
            } // End for loop - for all taxa
            
        }    // End sequential format
        
        std::vector<faststring> pattern_vec(posNum);
        
        unsigned k, l;
        
        faststring tmp;
        
        for (k=0; k<posNum; ++k)
        {
            tmp.clear();
            
            for (l=0; l < taxaNum; ++l)
            {
                tmp.push_back(all_seqs[l][k]);
            }
            pattern_vec[k] = tmp;
        }
        
        SitePattern temp(pattern_vec[0]); //seq_patterns
        // TODO: Use reference here as well
        vector <CSplit> site_partition = temp.gen_pattern_bitsets();
        vector<SitePattern> all_patterns(posNum);
        
        int m; //i
        int n; //j
        int o = 0; //k
        
        for (m = 0; m < posNum; ++m)
        {
            temp.reset_pattern(pattern_vec[m]);
            site_partition = temp.gen_morph_pattern_bitsets();
            // Define a container object first, then use add function to add site patterns. Call the add function and the add function will call the pushback
            all_patterns[m] = (temp);
        }
        
        
        
        double pa_sum;
        vector<double> site_rates(posNum);
        
        // If the alignment is larger than 500,000 sites, this small loop
        // just makes a subset of the alignment patterns for comparison,
        // I have hard coded it to 5...we can also do more or less later.
        // by changing the "subset_size" variable.
        if (posNum > 50000)
        {
            int subset_size = 25;
            int subset_length = posNum/subset_size;
            vector<SitePattern> subset_patterns(subset_length);
            // Make the subset of patterns
            for (m = 0; m < subset_length - 1; ++m)
            {
                subset_patterns[m] = all_patterns[o];
                o += 5;
            }
            
            for (m = 0; m < posNum; ++m)
            {
                if (i%1000 == 0)
                    cerr << m << " of " << posNum << " generations complete." << endl;
                pa_sum = 0;
                for (n = 0; n < subset_length - 1; ++n)
                {
                    pa_sum += all_patterns[m].calculate_pa_2(subset_patterns[n]);
                }
                site_rates[m] = (pa_sum/(subset_length - 1));
            }
        }
        // If the alignment is smaller than or equal to 500,000 sites than we do every comparison
        else
        {
            for (m = 0; m < posNum; ++m)
            {
                if (m%1000 == 0)
                    cerr << m << " of " << posNum << " generations complete." << endl;
                pa_sum = 0;
                for (n = 0; n < posNum; ++n)
                {
                    if (m == n)
                        continue;
                    else
                        pa_sum += all_patterns[m].calculate_pa_2(all_patterns[n]);
                }
                site_rates[m] = (pa_sum/(posNum - 1));
            }
        }
        faststring rates_name = filename;
        //        int l;
        //        for (l = 0; l < 3; ++l)
        //        {
        //            rates_name.pop_back();
        //        }
        rates_name.push_back('_');
        rates_name.push_back('r');
        rates_name.push_back('8');
        rates_name.push_back('s');
        rates_name.push_back('.');
        rates_name.push_back('t');
        rates_name.push_back('x');
        rates_name.push_back('t');
        ofstream rates_file;
        rates_file.open(rates_name.c_str());
        
        copy(site_rates.begin(), site_rates.end(), ostream_iterator<double>(rates_file, "\n"));
    }
    else
    {
        cerr << "Unknown data type: " << data_type << ". Please enter either dna or protein." << endl;
        return 0;
    }
    
    

    return 0;
}
