#ifndef DNA_ROUTINES_H
#define DNA_ROUTINES_H

/* Todo:  count DNA bases, guess data type, is AA, etc.  */

#include <cstring>
#include <cctype>

static const char DNARNA_symbols[]   = "aAcCgGtTuU";
/////////static const faststring aa_symbols     = "aArRnNdDbBcCeEqQzZgGhHiIlLkKmMfFpPsStTwWyYvVxX";
static const char aa_symbols[]       = "aArRnNdDcCeEqQzZgGhHiIlLkKmMfFpPsStTwWyYvVxX";

static const char aa_symbols_ext[]       = "aArRnNdDcCeEqQzZgGhHiIlLkKmMfFpPsStTwWyYvVxXUuOo*";

static const char non_aa_in_ABC [] = "BJOUZbjouz";

static const char DNARNA_iupac_symbols [] = "aAcCgGtTuURYSWryswKMBDkmbdHVNhvn";
static const char phylip_DNARNA_allowed_symbols [] = "?-nNaAcCgGtTuURYSWryswKMBDkmbdHVhvoOxX";

static const char phylip_neither_aa_dna_symbol [] = "jJzZ";
static const char allowed_non_ABC [] = "-?*";

// O and U are special amino acids that are usually not in the models.
// Chars that are not an aa symbol: B, J, O, U, Z
// X is the amibuity code


inline bool is_DNA_base_upper(char c)
{ 
  if (c == 'A' || c == 'C' || c == 'G' || c == 'T' )
    return true;
  return false;
}

inline bool is_DNA_base_lower(char c)
{ 
  if (c == 'a' || c == 'c' || c == 'g' || c == 't' )
    return true;
  return false;
}

inline bool is_DNA_base(char c)
{
  if (c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'a' || c == 'c' || c == 'g' || c == 't' )
    return true;
  return false;
}

inline bool is_DNA_or_RNA_base(char c)
{
  if (c == 'A' || c == 'C' || c == 'G' || c == 'T' || 
      c == 'a' || c == 'c' || c == 'g' || c == 't' || c == 'U' || c == 'u')
    return true;
  return false;
}

inline bool is_allowed_non_ABC(char c)
{
  if (c == '-' || c == '?' || c == '*')
    return true;
  else
    return false;
}

inline bool is_phylip_aa_dna_symbol(char c) // Does not include -? or similar
{
  if (c >='a' && c <= 'z')
  {  
      return true;
  }
  else if (c >='A' && c <= 'Z')
  {
      return true;
  }
  return false;
}


inline bool is_DNA_iupac_ambig(char c)
{
  if (c == 'R'  ||  c == 'Y'  || c == 'S' || c == 'W' ||
      c == 'r'  ||  c == 'y'  || c == 's' || c == 'w' ||
      c == 'K'  ||  c == 'M'  || c == 'B' || c == 'D' ||
      c == 'k'  ||  c == 'm'  || c == 'b' || c == 'd' ||
      c == 'H'  ||  c == 'V'  || c == 'N' ||
      c == 'h'  ||  c == 'v'  || c == 'n' || c == '?' )
    return true;
  else
    return false;
}

inline int is_aa_ambig(char c)
{
  if (c == 'b' || c == 'B' || c == '?' || c == 'x' || c == 'X' || 
      c == 'j' || c == 'J' || c == 'z' || c == 'Z')
    return true;
  else
    return false;
}



inline int is_aa_or_aa_ambig(char c) // includes '?' and  B, Z, J, X
{
  if (c >='a' && c <= 'z')
  {
    if (c != 'o' && c != 'u') // ou
      return true;
  }
  else if (c >='A' && c <= 'Z')
  {
    if (c != 'O' && c != 'U') // OU
      return true;
  }

  if (c == '?' || c == '*')
    return true;

  return false;
  //  return (strchr(aa_symbols, c) != NULL);
}

inline int is_aa_or_aa_ambig_extended(char c) // O and U are allowed  // includes '?' and  B, Z, J, X
{
  if (c >='a' && c <= 'z')
  {
      return true;
  }
  else if (c >='A' && c <= 'Z')
  {
      return true;
  }

  if (c == '?' || c == '*')
    return true;

  return false;
  //  return (strchr(aa_symbols, c) != NULL);
}


inline char complementBase_of_upper(char c)
{
  switch (c) {
  case 'A':  c = 'T'; break;
  case 'T':  c = 'A'; break;
  case 'G':  c = 'C'; break;
  case 'C':  c = 'G'; break;
    
  case 'R':  c = 'Y'; break;
  case 'Y':  c = 'R'; break;

  case 'M':  c = 'K'; break;
  case 'K':  c = 'M'; break;

  case 'B':  c = 'V'; break;
  case 'V':  c = 'B'; break;

  case 'D':  c = 'H'; break;
  case 'H':  c = 'D'; break;

  case 'N':  c = 'N'; break;
  case '?':  c = '?'; break;
  case '-':  c = '-'; break;

  default:
    c = '*';
    break;
  }
  return c;
}

inline char complementBase_of_lower(char c)
{
  switch (c) {
  case 'a':  c = 't'; break;
  case 't':  c = 'a'; break;
  case 'g':  c = 'c'; break;
  case 'c':  c = 'g'; break;
    
  case 'r':  c = 'y'; break;
  case 'y':  c = 'r'; break;

  case 'm':  c = 'k'; break;
  case 'k':  c = 'm'; break;

  case 'b':  c = 'v'; break;
  case 'v':  c = 'b'; break;

  case 'd':  c = 'h'; break;
  case 'h':  c = 'd'; break;

  case 'n':  c = 'n'; break;
  case '?':  c = '?'; break;
  case '-':  c = '-'; break;

  default:
    c = '*';
    break;
  }
  return c;
}

inline char complementBase(char c)
{
  if ( isupper(c) )
    return complementBase_of_upper(c);
  else
    return complementBase_of_lower(c);
}


inline void CG_AT_content_in_region(const char *beg, const char *end, unsigned &AT, unsigned &CG)
{
  AT = 0;
  CG = 0;

  char c;

  while (beg != end)
  {
    c = *beg;
    if (c == 'C' || c == 'c' || c == 'G' || c == 'g')
      ++CG;
    else if (c == 'A' || c == 'a' || c == 'T' || c == 't')
      ++AT;

    ++beg;
  }
}


#endif

