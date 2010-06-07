//////////////////////////////////////////////////////////////////////
// DuplexEngine.ipp
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Wrapper macros for certain model features.
//////////////////////////////////////////////////////////////////////

#include <list>
#include <utility>
#include <cassert>
#include <cstdlib>

//////////////////////////////////////////////////////////////////////
// DuplexEngine::IsComplementary()
//
// Determine if a pair of positions is considered "complementary."
//////////////////////////////////////////////////////////////////////

template<class RealT>
bool DuplexEngine<RealT>::IsComplementary(int i, int j) const
{
    Assert(1 <= i && i <= L1, "Index out-of-bounds.");
    Assert(1 <= j && j <= L2, "Index out-of-bounds.");

    return is_complementary[s1[i]][s2[j]];
}

//////////////////////////////////////////////////////////////////////
// DuplexEngine::DuplexEngine()
//
// Constructor
//////////////////////////////////////////////////////////////////////

template<class RealT>
DuplexEngine<RealT>::DuplexEngine(bool allow_noncomplementary) :
    allow_noncomplementary(allow_noncomplementary),
    cache_initialized(false),
    parameter_manager(NULL),
    L1(0), L2(0),
    SIZE(0)
{
    // precompute mapping from characters to index representation
    std::memset(char_mapping, BYTE(alphabet.size()), 256);
    for (size_t i = 0; i < alphabet.size(); i++)
    {
        char_mapping[BYTE(tolower(alphabet[i]))] = 
            char_mapping[BYTE(toupper(alphabet[i]))] = i;
    }
    
    // precompute complementary pairings
    for (int i = 0; i <= M; i++)
        for (int j = 0; j <= M; j++)
            is_complementary[i][j] = 0;
    
    is_complementary[char_mapping[BYTE('A')]][char_mapping[BYTE('U')]] = 
        is_complementary[char_mapping[BYTE('U')]][char_mapping[BYTE('A')]] = 
        is_complementary[char_mapping[BYTE('G')]][char_mapping[BYTE('U')]] = 
        is_complementary[char_mapping[BYTE('U')]][char_mapping[BYTE('G')]] = 
        is_complementary[char_mapping[BYTE('C')]][char_mapping[BYTE('G')]] = 
        is_complementary[char_mapping[BYTE('G')]][char_mapping[BYTE('C')]] = 1;
}

//////////////////////////////////////////////////////////////////////
// DuplexEngine::~DuplexEngine()
//
// Destructor.
//////////////////////////////////////////////////////////////////////

template<class RealT>
DuplexEngine<RealT>::~DuplexEngine()
{
}

//////////////////////////////////////////////////////////////////////
// DuplexEngine::RegisterParameters()
//
// Establish a mapping between parameters in the inference
// engine and parameters in the parameter manager.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void DuplexEngine<RealT>::RegisterParameters(ParameterManager<RealT> &parameter_manager)
{
    char buffer[1000];
    char buffer2[1000];

    cache_initialized = false;
    this->parameter_manager = &parameter_manager;
    parameter_manager.ClearParameters();
    
#if SINGLE_HYPERPARAMETER
    parameter_manager.AddParameterGroup("all_params");
#endif
    
#if PARAMS_BASE_PAIR
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("base_pair");
#endif
    for (int i = 0; i <= M; i++)
    {
        for (int j = 0; j <= M; j++)
        {
            if (i == M || j == M)
            {
                score_base_pair[i][j] = std::pair<RealT,RealT>(0, 0);
            }
            else
            {
                sprintf(buffer, "base_pair_%c%c", alphabet[i], alphabet[j]);
                sprintf(buffer2, "base_pair_%c%c", alphabet[j], alphabet[i]);
                if (strcmp(buffer, buffer2) < 0)
                    parameter_manager.AddParameterMapping(buffer, &score_base_pair[i][j]);
                else
                    parameter_manager.AddParameterMapping(buffer2, &score_base_pair[i][j]);
            }
        }
    }
#endif
    
#if PARAMS_BASE_PAIR_DIST
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("base_pair_dist_at_least");
#endif
    for (int i = 0; i < D_MAX_BP_DIST_THRESHOLDS; i++)
    {
        sprintf(buffer, "base_pair_dist_at_least_%d", BP_DIST_THRESHOLDS[i]);
        parameter_manager.AddParameterMapping(buffer, &score_base_pair_dist_at_least[i]);
    }
#endif
    
#if PARAMS_TERMINAL_MISMATCH
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("terminal_mismatch");
#endif
    for (int i1 = 0; i1 <= M; i1++)
    {
        for (int j1 = 0; j1 <= M; j1++)
        {
            for (int i2 = 0; i2 <= M; i2++)
            {
                for (int j2 = 0; j2 <= M; j2++)
                {
                    if (i1 == M || j1 == M || i2 == M || j2 == M)
                    {
                        score_terminal_mismatch[i1][j1][i2][j2] = std::pair<RealT,RealT>(0, 0);
                    }
                    else
                    {
                        sprintf(buffer, "terminal_mismatch_%c%c%c%c", alphabet[i1], alphabet[j1], alphabet[i2], alphabet[j2]);
                        parameter_manager.AddParameterMapping(buffer, &score_terminal_mismatch[i1][j1][i2][j2]);
                    }
                }
            }
        }
    }
#endif
    
#if PARAMS_HAIRPIN_LENGTH
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("hairpin_length_at_least");
#endif
    for (int i = 0; i <= D_MAX_HAIRPIN_LENGTH; i++)
    {
        sprintf(buffer, "hairpin_length_at_least_%d", i);
        parameter_manager.AddParameterMapping(buffer, &score_hairpin_length_at_least[i]);
    }
#endif
    
#if PARAMS_HAIRPIN_3_NUCLEOTIDES
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("hairpin_3_nucleotides");  
#endif
    for (int i1 = 0; i1 <= M; i1++)
    {
        for (int i2 = 0; i2 <= M; i2++)
        {
            for (int i3 = 0; i3 <= M; i3++)
            {
                if (i1 == M || i2 == M || i3 == M)
                {
                    score_hairpin_3_nucleotides[i1][i2][i3] = std::pair<RealT,RealT>(0, 0);
                }
                else
                {
                    sprintf(buffer, "hairpin_3_nucleotides_%c%c%c", alphabet[i1], alphabet[i2], alphabet[i3]);
                    parameter_manager.AddParameterMapping(buffer, &score_hairpin_3_nucleotides[i1][i2][i3]);
                }
            }
        }
    }
#endif
    
#if PARAMS_HAIRPIN_4_NUCLEOTIDES
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("hairpin_4_nucleotides");  
#endif
    for (int i1 = 0; i1 <= M; i1++)
    {
        for (int i2 = 0; i2 <= M; i2++)
        {
            for (int i3 = 0; i3 <= M; i3++)
            {
                for (int i4 = 0; i4 <= M; i4++)
                {
                    if (i1 == M || i2 == M || i3 == M || i4 == M)
                    {
                        score_hairpin_4_nucleotides[i1][i2][i3][i4] = std::pair<RealT,RealT>(0, 0);
                    }
                    else
                    {
                        sprintf(buffer, "hairpin_4_nucleotides_%c%c%c%c", alphabet[i1], alphabet[i2], alphabet[i3], alphabet[i4]);
                        parameter_manager.AddParameterMapping(buffer, &score_hairpin_4_nucleotides[i1][i2][i3][i4]);
                    }
                }
            }
        }
    }
#endif
    
#if PARAMS_HELIX_LENGTH
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("helix_length_at_least");
#endif
    for (int i = 0; i <= D_MAX_HELIX_LENGTH; i++)
    {
        if (i < 3)
        {
            score_helix_length_at_least[i] = std::pair<RealT,RealT>(0, 0);
        }
        else
        {
            sprintf(buffer, "helix_length_at_least_%d", i);
            parameter_manager.AddParameterMapping(buffer, &score_helix_length_at_least[i]);
        }
    }
#endif
    
#if PARAMS_ISOLATED_BASE_PAIR
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("isolated_base_pair");
#endif
    parameter_manager.AddParameterMapping("isolated_base_pair", &score_isolated_base_pair);
#endif
    
#if PARAMS_INTERNAL_EXPLICIT
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("internal_explicit");
#endif
    for (int i = 0; i <= D_MAX_INTERNAL_EXPLICIT_LENGTH; i++)
    {
        for (int j = 0; j <= D_MAX_INTERNAL_EXPLICIT_LENGTH; j++)
        {
            if (i == 0 || j == 0)
            {
                score_internal_explicit[i][j] = std::pair<RealT,RealT>(0, 0);
            }
            else
            {
                sprintf(buffer, "internal_explicit_%d_%d", std::min(i, j), std::max(i, j));
                parameter_manager.AddParameterMapping(buffer, &score_internal_explicit[i][j]);
            }
        }
    }
#endif

#if PARAMS_BULGE_LENGTH
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("bulge_length_at_least");
#endif
    for (int i = 0; i <= D_MAX_BULGE_LENGTH; i++)
    {
        if (i == 0)
        {
            score_bulge_length_at_least[i] = std::pair<RealT,RealT>(0, 0);
        }
        else
        {
            sprintf(buffer, "bulge_length_at_least_%d", i);
            parameter_manager.AddParameterMapping(buffer, &score_bulge_length_at_least[i]);
        }
    }
#endif
    
#if PARAMS_INTERNAL_LENGTH
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("internal_length_at_least");
#endif
    for (int i = 0; i <= D_MAX_INTERNAL_LENGTH; i++)
    {
        if (i < 2)
        {
            score_internal_length_at_least[i] = std::pair<RealT,RealT>(0, 0);
        }
        else
        {
            sprintf(buffer, "internal_length_at_least_%d", i);
            parameter_manager.AddParameterMapping(buffer, &score_internal_length_at_least[i]);
        }
    }
#endif
    
#if PARAMS_INTERNAL_SYMMETRY
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("internal_symmetric_length_at_least");
#endif
    for (int i = 0; i <= D_MAX_INTERNAL_SYMMETRIC_LENGTH; i++)
    {
        if (i == 0)
        {
            score_internal_symmetric_length_at_least[i] = std::pair<RealT,RealT>(0, 0);
        }
        else
        {
            sprintf(buffer, "internal_symmetric_length_at_least_%d", i);
            parameter_manager.AddParameterMapping(buffer, &score_internal_symmetric_length_at_least[i]);
        }
    }
#endif

#if PARAMS_INTERNAL_ASYMMETRY
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("internal_asymmetry_at_least");
#endif
    for (int i = 0; i <= D_MAX_INTERNAL_ASYMMETRY; i++)
    {
        if (i == 0)
        {
            score_internal_asymmetry_at_least[i] = std::pair<RealT,RealT>(0, 0);
        }
        else
        {
            sprintf(buffer, "internal_asymmetry_at_least_%d", i);
            parameter_manager.AddParameterMapping(buffer, &score_internal_asymmetry_at_least[i]);
        }
    }
#endif

#if PARAMS_BULGE_0x1_NUCLEOTIDES
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("bulge_0x1_nucleotides");
#endif
    for (int i1 = 0; i1 <= M; i1++)
    {
        if (i1 == M)
        {
            score_bulge_0x1_nucleotides[i1] = std::pair<RealT,RealT>(0, 0);
            score_bulge_1x0_nucleotides[i1] = std::pair<RealT,RealT>(0, 0);
        }
        else
        {
            sprintf(buffer, "bulge_0x1_nucleotides_%c", alphabet[i1]);
            parameter_manager.AddParameterMapping(buffer, &score_bulge_0x1_nucleotides[i1]);
            parameter_manager.AddParameterMapping(buffer, &score_bulge_1x0_nucleotides[i1]);
        }
    }
#endif

#if PARAMS_BULGE_0x2_NUCLEOTIDES
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("bulge_0x2_nucleotides");
#endif
    for (int i1 = 0; i1 <= M; i1++)
    {
        for (int i2 = 0; i2 <= M; i2++)
        {
            if (i1 == M || i2 == M)
            {
                score_bulge_0x2_nucleotides[i1][i2] = std::pair<RealT,RealT>(0, 0);
                score_bulge_2x0_nucleotides[i1][i2] = std::pair<RealT,RealT>(0, 0);
            }
            else
            {
                sprintf(buffer, "bulge_0x2_nucleotides_%c%c", alphabet[i1], alphabet[i2]);
                parameter_manager.AddParameterMapping(buffer, &score_bulge_0x2_nucleotides[i1][i2]);
                parameter_manager.AddParameterMapping(buffer, &score_bulge_2x0_nucleotides[i1][i2]);
            }
        }
    }
#endif
    
#if PARAMS_BULGE_0x3_NUCLEOTIDES
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("bulge_0x3_nucleotides");
#endif
    for (int i1 = 0; i1 <= M; i1++)
    {
        for (int i2 = 0; i2 <= M; i2++)
        {
            for (int i3 = 0; i3 <= M; i3++)
            {
                if (i1 == M || i2 == M)
                {
                    score_bulge_0x3_nucleotides[i1][i2][i3] = std::pair<RealT,RealT>(0, 0);
                    score_bulge_3x0_nucleotides[i1][i2][i3] = std::pair<RealT,RealT>(0, 0);
                }
                else
                {
                    sprintf(buffer, "bulge_0x3_nucleotides_%c%c%c", alphabet[i1], alphabet[i2], alphabet[i3]);
                    parameter_manager.AddParameterMapping(buffer, &score_bulge_0x3_nucleotides[i1][i2][i3]);
                    parameter_manager.AddParameterMapping(buffer, &score_bulge_3x0_nucleotides[i1][i2][i3]);
                }
            }
        }
    }
#endif

#if PARAMS_INTERNAL_1x1_NUCLEOTIDES
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("internal_1x1_nucleotides");
#endif
    for (int i1 = 0; i1 <= M; i1++)
    {
        for (int i2 = 0; i2 <= M; i2++)
        {
            if (i1 == M || i2 == M)
            {
                score_internal_1x1_nucleotides[i1][i2] = std::pair<RealT,RealT>(0, 0);
            }
            else 
            {          
                sprintf(buffer, "internal_1x1_nucleotides_%c%c", alphabet[i1], alphabet[i2]);
                sprintf(buffer2, "internal_1x1_nucleotides_%c%c", alphabet[i2], alphabet[i1]);
                if (strcmp(buffer, buffer2) < 0)
                    parameter_manager.AddParameterMapping(buffer, &score_internal_1x1_nucleotides[i1][i2]);
                else
                    parameter_manager.AddParameterMapping(buffer2, &score_internal_1x1_nucleotides[i1][i2]);
            }
        }
    }
#endif
  
#if PARAMS_INTERNAL_1x2_NUCLEOTIDES
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("internal_1x2_nucleotides");
#endif
    for (int i1 = 0; i1 <= M; i1++)
    {
        for (int i2 = 0; i2 <= M; i2++)
        {
            for (int i3 = 0; i3 <= M; i3++)
            {
                if (i1 == M || i2 == M || i3 == M)
                {
                    score_internal_1x2_nucleotides[i1][i2][i3] = std::pair<RealT,RealT>(0, 0);
                }
                else
                {
                    sprintf(buffer, "internal_1x2_nucleotides_%c%c%c", alphabet[i1], alphabet[i2], alphabet[i3]);
                    parameter_manager.AddParameterMapping(buffer, &score_internal_1x2_nucleotides[i1][i2][i3]);
                    sprintf(buffer, "internal_2x1_nucleotides_%c%c%c", alphabet[i1], alphabet[i2], alphabet[i3]);
                    parameter_manager.AddParameterMapping(buffer, &score_internal_2x1_nucleotides[i1][i2][i3]);
                }
            }
        }
    }
#endif
  
#if PARAMS_INTERNAL_2x2_NUCLEOTIDES
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("internal_2x2_nucleotides");
#endif
    for (int i1 = 0; i1 <= M; i1++)
    {
        for (int i2 = 0; i2 <= M; i2++)
        {
            for (int i3 = 0; i3 <= M; i3++)
            {
                for (int i4 = 0; i4 <= M; i4++)
                {
                    if (i1 == M || i2 == M || i3 == M || i4 == M)
                    {
                        score_internal_2x2_nucleotides[i1][i2][i3][i4] = std::pair<RealT,RealT>(0, 0);
                    }
                    else
                    {
                        sprintf(buffer, "internal_2x2_nucleotides_%c%c%c%c", alphabet[i1], alphabet[i2], alphabet[i3], alphabet[i4]);
                        sprintf(buffer2, "internal_2x2_nucleotides_%c%c%c%c", alphabet[i3], alphabet[i4], alphabet[i1], alphabet[i2]);
                        if (strcmp(buffer, buffer2) < 0)
                            parameter_manager.AddParameterMapping(buffer, &score_internal_2x2_nucleotides[i1][i2][i3][i4]);
                        else
                            parameter_manager.AddParameterMapping(buffer2, &score_internal_2x2_nucleotides[i1][i2][i3][i4]);
                    }
                }
            }
        }
    }
#endif

#if PARAMS_HELIX_STACKING
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("helix_stacking");
#endif
    for (int i1 = 0; i1 <= M; i1++)
    {
        for (int j1 = 0; j1 <= M; j1++)
        {
            for (int i2 = 0; i2 <= M; i2++)
            {
                for (int j2 = 0; j2 <= M; j2++)
                {
                    if (i1 == M || j1 == M || i2 == M || j2 == M)
                    {
                        score_helix_stacking[i1][j1][i2][j2] = std::pair<RealT,RealT>(0, 0);
                    }
                    else
                    {
                        sprintf(buffer, "helix_stacking_%c%c%c%c", alphabet[i1], alphabet[j1], alphabet[i2], alphabet[j2]);
                        sprintf(buffer2, "helix_stacking_%c%c%c%c", alphabet[j2], alphabet[i2], alphabet[j1], alphabet[i1]);
                        if (strcmp(buffer, buffer2) < 0)
                            parameter_manager.AddParameterMapping(buffer, &score_helix_stacking[i1][j1][i2][j2]);
                        else
                            parameter_manager.AddParameterMapping(buffer2, &score_helix_stacking[i1][j1][i2][j2]);
                    }
                }
            }
        }
    }
#endif

#if PARAMS_HELIX_CLOSING
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("helix_closing");
#endif
    for (int i = 0; i <= M; i++)
    {
        for (int j = 0; j <= M; j++)
        {
            if (i == M || j == M)
            {
                score_helix_closing[i][j] = std::pair<RealT,RealT>(0, 0);
            }
            else
            {
                sprintf(buffer, "helix_closing_%c%c", alphabet[i], alphabet[j]);
                parameter_manager.AddParameterMapping(buffer, &score_helix_closing[i][j]);
            }
        }
    }
#endif

#if PARAMS_MULTI_LENGTH
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("multi_length");
#endif
    parameter_manager.AddParameterMapping("multi_base", &score_multi_base);
    parameter_manager.AddParameterMapping("multi_unpaired", &score_multi_unpaired);
    parameter_manager.AddParameterMapping("multi_paired", &score_multi_paired);
#endif

#if PARAMS_DANGLE
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("dangle");
#endif
    for (int i1 = 0; i1 <= M; i1++)
    {
        for (int j1 = 0; j1 <= M; j1++)
        {
            for (int i2 = 0; i2 <= M; i2++)
            {
                if (i1 == M || j1 == M || i2 == M)
                {
                    score_dangle_left[i1][j1][i2] = std::pair<RealT,RealT>(0, 0);
                }
                else
                {
                    sprintf(buffer, "dangle_left_%c%c%c", alphabet[i1], alphabet[j1], alphabet[i2]);
                    parameter_manager.AddParameterMapping(buffer, &score_dangle_left[i1][j1][i2]);
                }
            }
        }
    }
  
    for (int i1 = 0; i1 <= M; i1++)
    {
        for (int j1 = 0; j1 <= M; j1++)
        {
            for (int j2 = 0; j2 <= M; j2++)
            {
                if (i1 == M || j1 == M || j2 == M)
                {
                    score_dangle_right[i1][j1][j2] = std::pair<RealT,RealT>(0, 0);
                }
                else
                {
                    sprintf(buffer, "dangle_right_%c%c%c", alphabet[i1], alphabet[j1], alphabet[j2]);
                    parameter_manager.AddParameterMapping(buffer, &score_dangle_right[i1][j1][j2]);
                }
            }
        }
    }
#endif

#if PARAMS_EXTERNAL_LENGTH
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("external_length");
#endif
    parameter_manager.AddParameterMapping("external_unpaired", &score_external_unpaired);
    parameter_manager.AddParameterMapping("external_paired", &score_external_paired);
#endif

}

//////////////////////////////////////////////////////////////////////
// DuplexEngine::LoadSequence()
//
// Load an RNA sequence.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void DuplexEngine<RealT>::LoadSequence(const SStruct& str1, const SStruct& str2)
{
    cache_initialized = false;
    
    // compute dimensions
    L1 = str1.GetLength();
    L2 = str2.GetLength();
    SIZE = (L1+1)*(L2+1);
    
    // allocate memory
    offset.resize(L1+1);
    s1.resize(L1+1);
    s2.resize(L2+1);
    allow_unpaired.resize(SIZE);
    allow_paired.resize(SIZE);

    // convert sequences to index representation
    const std::string &seq1 = str1.GetSequences()[0];
    s1[0] = BYTE(alphabet.size());
    for (int i = 1; i <= L1; i++)
        s1[i] = char_mapping[BYTE(seq1[i])];
    const std::string &seq2 = str2.GetSequences()[0];
    s2[0] = BYTE(alphabet.size());
    for (int i = 1; i <= L2; i++)
        s2[i] = char_mapping[BYTE(seq2[i])];

    // compute indexing scheme for two dimension array;
    for (int i = 0; i <= L1; i++) offset[i] = (L2+1)*i;

    // allow all ranges to be unpaired, and all pairs of letters
    // to be paired;
    for (int i = 0; i < SIZE; i++)
        allow_unpaired[i] = allow_paired[i] = 1;

    // prevent the non-letter before each sequence from pairing with anything;
    for (int i = 0; i <= L1; i++) allow_paired[offset[i]+0] = 0;
    for (int j = 0; j <= L2; j++) allow_paired[offset[0]+j] = 0;

    // enforce complementarity of base-pairings
    if (!allow_noncomplementary)
    {
        // for each pair of non-complementary letters in the sequence, disallow the pairing
        for (int i = 1; i <= L1; i++)
        {
            for (int j = 1; j <= L2; j++)
                if (!IsComplementary(i,j)) allow_paired[offset[i]+j] = 0;
        }
    }
}

//////////////////////////////////////////////////////////////////////
// DuplexEngine::InitializeCache()
//
// Initialize scoring cache prior to inference.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void DuplexEngine<RealT>::InitializeCache()
{
    if (cache_initialized) return;
    cache_initialized = true;

    // initialize length and distance scoring
#if PARAMS_HELIX_LENGTH
    cache_score_helix_length[0].first = score_helix_length_at_least[0].first;
    for (int i = 1; i <= D_MAX_HELIX_LENGTH; i++)
        cache_score_helix_length[i].first = cache_score_helix_length[i-1].first + score_helix_length_at_least[i].first;
#endif

#if PARAMS_BULGE_LENGTH
    RealT temp_cache_score_bulge_length[D_MAX_BULGE_LENGTH+1];
    temp_cache_score_bulge_length[0] = score_bulge_length_at_least[0].first;
    for (int i = 1; i <= D_MAX_BULGE_LENGTH; i++)
        temp_cache_score_bulge_length[i] = temp_cache_score_bulge_length[i-1] + score_bulge_length_at_least[i].first;
#endif
    
#if PARAMS_INTERNAL_LENGTH
    RealT temp_cache_score_internal_length[D_MAX_INTERNAL_LENGTH+1];
    temp_cache_score_internal_length[0] = score_internal_length_at_least[0].first;
    for (int i = 1; i <= D_MAX_INTERNAL_LENGTH; i++)
        temp_cache_score_internal_length[i] = temp_cache_score_internal_length[i-1] + score_internal_length_at_least[i].first;
#endif
    
#if PARAMS_INTERNAL_SYMMETRY
    RealT temp_cache_score_internal_symmetric_length[D_MAX_INTERNAL_SYMMETRIC_LENGTH+1];
    temp_cache_score_internal_symmetric_length[0] = score_internal_symmetric_length_at_least[0].first;
    for (int i = 1; i <= D_MAX_INTERNAL_SYMMETRIC_LENGTH; i++)
        temp_cache_score_internal_symmetric_length[i] = temp_cache_score_internal_symmetric_length[i-1] + score_internal_symmetric_length_at_least[i].first;
#endif
    
#if PARAMS_INTERNAL_ASYMMETRY
    RealT temp_cache_score_internal_asymmetry[D_MAX_INTERNAL_ASYMMETRY+1];
    temp_cache_score_internal_asymmetry[0] = score_internal_asymmetry_at_least[0].first;
    for (int i = 1; i <= D_MAX_INTERNAL_ASYMMETRY; i++)
        temp_cache_score_internal_asymmetry[i] = temp_cache_score_internal_asymmetry[i-1] + score_internal_asymmetry_at_least[i].first;
#endif
    
    // precompute score for single-branch loops of length l1 and l2
    for (int l1 = 0; l1 <= C_MAX_SINGLE_LENGTH; l1++)
    {
        for (int l2 = 0; l1+l2 <= C_MAX_SINGLE_LENGTH; l2++)
        {
            cache_score_single[l1][l2].first = RealT(0);

            // skip over stacking pairs
            if (l1 == 0 && l2 == 0) continue;

            // consider bulge loops
            if (l1 == 0 || l2 == 0)
            {
#if PARAMS_BULGE_LENGTH
                cache_score_single[l1][l2].first += temp_cache_score_bulge_length[std::min(D_MAX_BULGE_LENGTH, l1+l2)];
#endif
            }

            // consider internal loops
            else
            {
#if PARAMS_INTERNAL_EXPLICIT
                if (l1 <= D_MAX_INTERNAL_EXPLICIT_LENGTH && l2 <= D_MAX_INTERNAL_EXPLICIT_LENGTH)
                    cache_score_single[l1][l2].first += score_internal_explicit[l1][l2].first;
#endif
#if PARAMS_INTERNAL_LENGTH
                cache_score_single[l1][l2].first += temp_cache_score_internal_length[std::min(D_MAX_INTERNAL_LENGTH, l1+l2)];
#endif
#if PARAMS_INTERNAL_SYMMETRY
                if (l1 == l2)
                    cache_score_single[l1][l2].first += temp_cache_score_internal_symmetric_length[std::min(D_MAX_INTERNAL_SYMMETRIC_LENGTH, l1)];
#endif
#if PARAMS_INTERNAL_ASYMMETRY
                cache_score_single[l1][l2].first += temp_cache_score_internal_asymmetry[std::min(D_MAX_INTERNAL_ASYMMETRY, Abs(l1-l2))];
#endif
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////
// DuplexEngine::LoadValues()
//
// Load parameter values.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void DuplexEngine<RealT>::LoadValues(const std::vector<RealT> &values)
{
    if (values.size() != parameter_manager->GetNumLogicalParameters()) Error("Parameter size mismatch.");
    
    cache_initialized = false;
    for (size_t i = 0; i < values.size(); i++)
    {
        std::vector<std::pair<RealT,RealT> *> physical_parameters = parameter_manager->GetPhysicalParameters(i);
        for (size_t j = 0; j < physical_parameters.size(); j++)
            physical_parameters[j]->first = values[i];
    }
}

//////////////////////////////////////////////////////////////////////
// DuplexEngine::GetCounts()
//
// Return counts.
//////////////////////////////////////////////////////////////////////

template<class RealT>
std::vector<RealT> DuplexEngine<RealT>::GetCounts()
{
    std::vector<RealT> counts(parameter_manager->GetNumLogicalParameters());
    
    // clear counts for physical parameters
    for (size_t i = 0; i < parameter_manager->GetNumLogicalParameters(); i++)
    {
        std::vector<std::pair<RealT,RealT> *> physical_parameters = parameter_manager->GetPhysicalParameters(i);
        for (size_t j = 0; j < physical_parameters.size(); j++)
            counts[i] += physical_parameters[j]->second;
    }

    return counts;
}

//////////////////////////////////////////////////////////////////////
// DuplexEngine::ClearCounts()
//
// Set all counts to zero.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void DuplexEngine<RealT>::ClearCounts()
{
    // clear counts for physical parameters
    for (size_t i = 0; i < parameter_manager->GetNumLogicalParameters(); i++)
    {
        std::vector<std::pair<RealT,RealT> *> physical_parameters = parameter_manager->GetPhysicalParameters(i);
        for (size_t j = 0; j < physical_parameters.size(); j++)
            physical_parameters[j]->second = RealT(0);
    }

    // clear counts for cache
#if PARAMS_BASE_PAIR_DIST
    for (int i = 0; i <= BP_DIST_LAST_THRESHOLD; i++)
        cache_score_base_pair_dist[i].second = RealT(0);
#endif

#if PARAMS_HAIRPIN_LENGTH
    for (int i = 0; i <= D_MAX_HAIRPIN_LENGTH; i++)
        cache_score_hairpin_length[i].second = RealT(0);
#endif

#if PARAMS_HELIX_LENGTH
    for (int i = 0; i <= D_MAX_HELIX_LENGTH; i++)
        cache_score_helix_length[i].second = RealT(0);
#endif

    for (int l1 = 0; l1 <= C_MAX_SINGLE_LENGTH; l1++)
        for (int l2 = 0; l2 <= C_MAX_SINGLE_LENGTH; l2++)
            cache_score_single[l1][l2].second = RealT(0);
}

//////////////////////////////////////////////////////////////////////
// DuplexEngine::FinalizeCounts()
//
// Apply any needed transformations to counts.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void DuplexEngine<RealT>::FinalizeCounts()
{
    // perform transformations
#if PARAMS_BASE_PAIR_DIST
    for (int i = 0; i < D_MAX_BP_DIST_THRESHOLDS; i++)
        for (int j = BP_DIST_THRESHOLDS[i]; j <= BP_DIST_LAST_THRESHOLD; j++)
            score_base_pair_dist_at_least[i].second += cache_score_base_pair_dist[j].second;
#endif
    
#if PARAMS_HAIRPIN_LENGTH
    for (int i = 0; i <= D_MAX_HAIRPIN_LENGTH; i++)
        for (int j = i; j <= D_MAX_HAIRPIN_LENGTH; j++)
            score_hairpin_length_at_least[i].second += cache_score_hairpin_length[j].second;
#endif
    
#if PARAMS_HELIX_LENGTH
    for (int i = 0; i <= D_MAX_HELIX_LENGTH; i++)
        for (int j = i; j <= D_MAX_HELIX_LENGTH; j++)
            score_helix_length_at_least[i].second += cache_score_helix_length[j].second;
#endif

    // allocate temporary storage
#if PARAMS_BULGE_LENGTH
    RealT temp_cache_counts_bulge_length[D_MAX_BULGE_LENGTH+1];
    std::fill(temp_cache_counts_bulge_length, temp_cache_counts_bulge_length + D_MAX_BULGE_LENGTH+1, RealT(0));
#endif
    
#if PARAMS_INTERNAL_LENGTH
    RealT temp_cache_counts_internal_length[D_MAX_INTERNAL_LENGTH+1];
    std::fill(temp_cache_counts_internal_length, temp_cache_counts_internal_length + D_MAX_INTERNAL_LENGTH+1, RealT(0));
#endif
    
#if PARAMS_INTERNAL_SYMMETRY
    RealT temp_cache_counts_internal_symmetric_length[D_MAX_INTERNAL_SYMMETRIC_LENGTH+1];
    std::fill(temp_cache_counts_internal_symmetric_length, temp_cache_counts_internal_symmetric_length + D_MAX_INTERNAL_SYMMETRIC_LENGTH+1, RealT(0));
#endif
    
#if PARAMS_INTERNAL_ASYMMETRY
    RealT temp_cache_counts_internal_asymmetry[D_MAX_INTERNAL_ASYMMETRY+1];
    std::fill(temp_cache_counts_internal_asymmetry, temp_cache_counts_internal_asymmetry + D_MAX_INTERNAL_ASYMMETRY+1, RealT(0));
#endif

    // compute contributions
    for (int l1 = 0; l1 <= C_MAX_SINGLE_LENGTH; l1++)
    {
        for (int l2 = 0; l1+l2 <= C_MAX_SINGLE_LENGTH; l2++)
        {
            // skip over stacking pairs
            if (l1 == 0 && l2 == 0) continue;

            // consider bulge loops
            if (l1 == 0 || l2 == 0)
            {
#if PARAMS_BULGE_LENGTH
                temp_cache_counts_bulge_length[std::min(D_MAX_BULGE_LENGTH, l1+l2)] += cache_score_single[l1][l2].second;
#endif
            }

            // consider internal loops
            else
            {
#if PARAMS_INTERNAL_EXPLICIT
                if (l1 <= D_MAX_INTERNAL_EXPLICIT_LENGTH && l2 <= D_MAX_INTERNAL_EXPLICIT_LENGTH)
                    score_internal_explicit[l1][l2].second += cache_score_single[l1][l2].second;
#endif
#if PARAMS_INTERNAL_LENGTH
                temp_cache_counts_internal_length[std::min(D_MAX_INTERNAL_LENGTH, l1+l2)] += cache_score_single[l1][l2].second;
#endif
#if PARAMS_INTERNAL_SYMMETRY
                if (l1 == l2)
                    temp_cache_counts_internal_symmetric_length[std::min(D_MAX_INTERNAL_SYMMETRIC_LENGTH, l1)] += cache_score_single[l1][l2].second;
#endif
#if PARAMS_INTERNAL_ASYMMETRY
                temp_cache_counts_internal_asymmetry[std::min(D_MAX_INTERNAL_ASYMMETRY, Abs(l1-l2))] += cache_score_single[l1][l2].second;
#endif
            }
        }
    }

#if PARAMS_BULGE_LENGTH
    for (int i = 0; i <= D_MAX_BULGE_LENGTH; i++)
        for (int j = i; j <= D_MAX_BULGE_LENGTH; j++)
            score_bulge_length_at_least[i].second += temp_cache_counts_bulge_length[j];
#endif
    
#if PARAMS_INTERNAL_LENGTH
    for (int i = 0; i <= D_MAX_INTERNAL_LENGTH; i++)
        for (int j = i; j <= D_MAX_INTERNAL_LENGTH; j++)
            score_internal_length_at_least[i].second += temp_cache_counts_internal_length[j];
#endif
    
#if PARAMS_INTERNAL_SYMMETRY
    for (int i = 0; i <= D_MAX_INTERNAL_SYMMETRIC_LENGTH; i++)
        for (int j = i; j <= D_MAX_INTERNAL_SYMMETRIC_LENGTH; j++)
            score_internal_symmetric_length_at_least[i].second += temp_cache_counts_internal_symmetric_length[j];
#endif
    
#if PARAMS_INTERNAL_ASYMMETRY
    for (int i = 0; i <= D_MAX_INTERNAL_ASYMMETRY; i++)
        for (int j = i; j <= D_MAX_INTERNAL_ASYMMETRY; j++)
            score_internal_asymmetry_at_least[i].second += temp_cache_counts_internal_asymmetry[j];
#endif
}

// score for leaving s[i] unpaired

#define ScoreUnpairedPosition(i) (RealT(0))
#define CountUnpairedPosition(i,v)

// score for leaving s[i+1...j] unpaired

#define ScoreUnpaired(i,j) (RealT(0))
#define CountUnpaired(i,j,v)

//////////////////////////////////////////////////////////////////////
// DuplexEngine::LoopScore()
//
// Returns the score for nucleotides in a single-branch loop 
// spanning i to j and p to q.
//
// In an RNA structure, this would look like
//
//                       ...      ...
//                        |        |
//                      x[p] -- y[q]
// position p -------->  o         o  <----- position q
//                    x[p-1]     y[q+1]
//                      |           |
//                     ...         ...
//                      |           |
//                   x[i+1]       y[j-1]
// position i -------->  o         o  <----- position j
//                      x[i] -- x[j]
//                        |        |
//                     x[i-1] -- x[j+1]
//
//////////////////////////////////////////////////////////////////////

template<class RealT>
inline RealT DuplexEngine<RealT>::LoopScore(int i, int j, int p, int q) const
{
    Assert(0 < i && i < p && p <= L1, "Single-branch loop boundaries invalid.");
    Assert(0 < q && q < j && j <= L2, "Single-branch loop boundaries invalid.");

    const int l1 = p-i-1;
    const int l2 = j-q-1;
    
#if (!defined(NDEBUG) || PARAMS_BULGE_0x1_NUCLEOTIDES || PARAMS_BULGE_0x2_NUCLEOTIDES || PARAMS_BULGE_0x3_NUCLEOTIDES || PARAMS_INTERNAL_1x1_NUCLEOTIDES || PARAMS_INTERNAL_1x2_NUCLEOTIDES || PARAMS_INTERNAL_2x2_NUCLEOTIDES)
    Assert(l1 + l2 > 0 && l1 >= 0 && l2 >= 0 && l1 + l2 <= C_MAX_SINGLE_LENGTH, "Invalid single-branch loop size.");
#endif
    
    return 
        ScoreUnpaired(i,p)
        + ScoreUnpaired(q,j)
#if PARAMS_BULGE_0x1_NUCLEOTIDES
        + (l1 == 0 && l2 == 1 ? score_bulge_0x1_nucleotides[s2[j-1]].first : RealT(0))
        + (l1 == 1 && l2 == 0 ? score_bulge_1x0_nucleotides[s1[i+1]].first : RealT(0))
#endif
#if PARAMS_BULGE_0x2_NUCLEOTIDES
        + (l1 == 0 && l2 == 2 ? score_bulge_0x2_nucleotides[s2[j-2]][s2[j-1]].first : RealT(0))
        + (l1 == 2 && l2 == 0 ? score_bulge_2x0_nucleotides[s1[i+1]][s1[i+2]].first : RealT(0))
#endif
#if PARAMS_BULGE_0x3_NUCLEOTIDES
        + (l1 == 0 && l2 == 3 ? score_bulge_0x3_nucleotides[s2[j-3]][s2[j-2]][s2[j-1]].first : RealT(0))
        + (l1 == 3 && l2 == 0 ? score_bulge_3x0_nucleotides[s1[i+1]][s1[i+2]][s1[i+3]].first : RealT(0))
#endif
#if PARAMS_INTERNAL_1x1_NUCLEOTIDES
        + (l1 == 1 && l2 == 1 ? score_internal_1x1_nucleotides[s1[i+1]][s2[j-1]].first : RealT(0))
#endif
#if PARAMS_INTERNAL_1x2_NUCLEOTIDES
        + (l1 == 1 && l2 == 2 ? score_internal_1x2_nucleotides[s1[i+1]][s2[j-2]][s2[j-1]].first : RealT(0))
        + (l1 == 2 && l2 == 1 ? score_internal_2x1_nucleotides[s1[i+1]][s1[i+2]][s2[j-1]].first : RealT(0))
#endif
#if PARAMS_INTERNAL_2x2_NUCLEOTIDES
        + (l1 == 2 && l2 == 2 ? score_internal_2x2_nucleotides[s1[i+1]][s1[i+2]][s2[j-2]][s2[j-1]].first : RealT(0))
#endif
        ;
}

template<class RealT>
RealT DuplexEngine<RealT>::ComputeInside()
{
    InitializeCache();

    LogPartFunc=NEG_INF;
    // initialize the inside table
    inside.clear(); inside.resize(SIZE, RealT(NEG_INF));
    
    for (int i=1; i<=L1; ++i)
    {
        for (int j=L2; j>0; --j)
        {
            if (!allow_paired[offset[i]+j]) continue;

            {
                RealT sum = score_external_unpaired.first * (i-1 + L2-j);
                if (i>1) sum += score_dangle_right[s2[j]][s1[i]][s1[i-1]].first;
                if (j<L2) sum += score_dangle_left[s2[j]][s1[i]][s2[j+1]].first;
                sum += score_base_pair[s2[j]][s1[i]].first;
                sum += score_helix_closing[s2[j]][s1[i]].first;
                inside[offset[i]+j] = Fast_LogAdd(inside[offset[i]+j], sum);
            }

            int p_min=std::max(i-C_MAX_SINGLE_LENGTH, 1);
            for (int p=i-1; p>=p_min; --p)
            {
                int l1=i-p-1;
                int q_max=std::min(C_MAX_SINGLE_LENGTH+p-i+j, L2);
                for (int q=j+1; q<=q_max; ++q) 
                {
                    int l2=q-j-1;
                    if (!allow_paired[offset[p]+q]) continue;

                    RealT sum = inside[offset[p]+q];
                    if (l1==0 && l2==0)
                    {
                        sum += 
                            score_base_pair[s1[i]][s2[j]].first +
                            score_helix_stacking[s1[p]][s2[q]][s1[i]][s2[j]].first;
                    }
                    else
                    {
                        sum += score_terminal_mismatch[s1[p]][s2[q]][s1[p+1]][s2[q-1]].first
                            + score_terminal_mismatch[s2[j]][s1[i]][s2[j+1]][s1[i-1]].first
                            + score_base_pair[s1[i]][s2[j]].first
                            + LoopScore(p, q, i, j);
                    }
                    inside[offset[i]+j] = Fast_LogAdd(inside[offset[i]+j], sum);
                }
            }

            {
                RealT sum = inside[offset[i]+j];
                sum += score_external_unpaired.first * (L1-i + j-1);
                if (i<L1) sum += score_dangle_left[s1[i]][s2[j]][s1[i+1]].first;
                if (j>1) sum += score_dangle_right[s1[i]][s2[j]][s2[j-1]].first;
                sum += score_helix_closing[s1[i]][s2[j]].first;
                LogPartFunc = Fast_LogAdd(LogPartFunc, sum);
            }
        }
    }
    return LogPartFunc;
}

template<class RealT>
RealT DuplexEngine<RealT>::ComputeOutside()
{
    InitializeCache();

    RealT LogPartFunc=NEG_INF;
    // initialize the outside table
    outside.clear(); outside.resize(SIZE, RealT(NEG_INF));
    
    for (int i=L1; i>0; --i)
    {
        for (int j=1; j<=L2; ++j)
        {
            if (!allow_paired[offset[i]+j]) continue;

            {
                RealT sum = 0;
                sum += score_external_unpaired.first * (L1-i + j-1);
                if (i<L1) sum += score_dangle_left[s1[i]][s2[j]][s1[i+1]].first;
                if (j>1) sum += score_dangle_right[s1[i]][s2[j]][s2[j-1]].first;
                sum += score_helix_closing[s1[i]][s2[j]].first;
                outside[offset[i]+j] = Fast_LogAdd(outside[offset[i]+j], sum);
            }

            int p_min=std::max(i-C_MAX_SINGLE_LENGTH, 1);
            for (int p=i-1; p>=p_min; --p)
            {
                int l1=i-p-1;
                int q_max=std::min(C_MAX_SINGLE_LENGTH+p-i+j, L2);
                for (int q=j+1; q<=q_max; ++q) 
                {
                    int l2=q-j-1;
                    if (!allow_paired[offset[p]+q]) continue;

                    RealT sum = outside[offset[i]+j];
                    if (l1==0 && l2==0)
                    {
                        sum += 
                            score_base_pair[s1[i]][s2[j]].first +
                            score_helix_stacking[s1[p]][s2[q]][s1[i]][s2[j]].first;
                    }
                    else
                    {
                        sum += score_terminal_mismatch[s1[p]][s2[q]][s1[p+1]][s2[q-1]].first
                            + score_terminal_mismatch[s2[j]][s1[i]][s2[j+1]][s1[i-1]].first
                            + score_base_pair[s1[i]][s2[j]].first
                            + LoopScore(p, q, i, j);
                    }
                    outside[offset[p]+q] = Fast_LogAdd(outside[offset[p]+q], sum);
                }
            }

            {
                RealT sum = outside[offset[i]+j];
                sum += score_external_unpaired.first * (i-1 + L2-j);
                if (i>1) sum += score_dangle_right[s2[j]][s1[i]][s1[i-1]].first;
                if (j<L2) sum += score_dangle_left[s2[j]][s1[i]][s2[j+1]].first;
                sum += score_base_pair[s2[j]][s1[i]].first;
                sum += score_helix_closing[s2[j]][s1[i]].first;
                LogPartFunc = Fast_LogAdd(LogPartFunc, sum);
            }
        }
    }
    return LogPartFunc;
}

template<class RealT>
void DuplexEngine<RealT>::ComputePosterior()
{
    InitializeCache();

    // initialize the posterior table
    posterior.clear(); posterior.resize(SIZE, 0.0);
    posterior2.clear(); posterior2.resize(SIZE, 0.0);
    
    for (int i=1; i<=L1; ++i)
    {
        for (int j=L2; j>0; --j)
        {
            if (!allow_paired[offset[i]+j]) continue;
            posterior[offset[i]+j] = Fast_Exp(inside[offset[i]+j]+outside[offset[i]+j]-LogPartFunc);
            if (i>1 && j<L2 && allow_paired[offset[i-1]+j+1])
            {
                RealT e = score_base_pair[s1[i]][s2[j]].first +
                    score_helix_stacking[s1[i-1]][s2[j+1]][s1[i]][s2[j]].first;
                posterior2[offset[i]+j] =
                    Fast_Exp(inside[offset[i-1]+j+1]+outside[offset[i]+j]+e-LogPartFunc);
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::GetPosterior()
//
// Return posterior probability matrix, thresholded.
//////////////////////////////////////////////////////////////////////

template<class RealT>
RealT *DuplexEngine<RealT>::GetPosterior(const RealT posterior_cutoff) const
{
    RealT *ret = new RealT[SIZE];
    for (int i = 0; i < SIZE; i++)
        ret[i] = (posterior[i] >= posterior_cutoff ? posterior[i] : RealT(0));
    return ret;
}

template<class RealT>
RealT *DuplexEngine<RealT>::GetPosterior(const RealT posterior_cutoff,
					    std::vector<RealT>& p) const
{
    p.resize(SIZE);
    for (int i = 0; i < SIZE; i++)
        p[i] = (posterior[i] >= posterior_cutoff ? posterior[i] : RealT(0));
    return &p[0];
}

template<class RealT>
RealT *DuplexEngine<RealT>::GetPosterior2(const RealT posterior_cutoff) const
{
    RealT *ret = new RealT[SIZE];
    for (int i = 0; i < SIZE; i++)
        ret[i] = (posterior2[i] >= posterior_cutoff ? posterior2[i] : RealT(0));
    return ret;
}

template<class RealT>
RealT *DuplexEngine<RealT>::GetPosterior2(const RealT posterior_cutoff,
					    std::vector<RealT>& p) const
{
    p.resize(SIZE);
    for (int i = 0; i < SIZE; i++)
        p[i] = (posterior2[i] >= posterior_cutoff ? posterior2[i] : RealT(0));
    return &p[0];
}
