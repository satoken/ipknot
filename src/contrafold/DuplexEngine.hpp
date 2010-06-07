//////////////////////////////////////////////////////////////////////
// DuplexEngine.hpp
//////////////////////////////////////////////////////////////////////

#ifndef DUPLEXENGINE_HPP
#define DUPLEXENGINE_HPP

#include <queue>
#include <vector>
#include <string>
#include "Config.hpp"
#include "SStruct.hpp"
#include "ParameterManager.hpp"
#include "Utilities.hpp"
#include "LogSpace.hpp"

//////////////////////////////////////////////////////////////////////
// class DuplexEngine
//////////////////////////////////////////////////////////////////////

template<class RealT>
class DuplexEngine
{
    const bool allow_noncomplementary;
    unsigned char char_mapping[256];
    int is_complementary[M+1][M+1];
    bool cache_initialized;
    ParameterManager<RealT> *parameter_manager;
    
    // dimensions
    int L1, L2, SIZE;

    // sequence data
    std::vector<int> s1, s2;
    std::vector<int> offset;
    std::vector<int> allow_unpaired, allow_paired;

    std::vector<RealT> inside, outside;
    std::vector<RealT> posterior, posterior2;
    RealT LogPartFunc;

    // parameters

#if PARAMS_BASE_PAIR
    std::pair<RealT,RealT> score_base_pair[M+1][M+1];
#endif
#if PARAMS_BASE_PAIR_DIST
    std::pair<RealT,RealT> score_base_pair_dist_at_least[D_MAX_BP_DIST_THRESHOLDS];
    std::pair<RealT,RealT> cache_score_base_pair_dist[BP_DIST_LAST_THRESHOLD+1];
#endif
#if PARAMS_TERMINAL_MISMATCH
    std::pair<RealT,RealT> score_terminal_mismatch[M+1][M+1][M+1][M+1];
#endif
#if PARAMS_HAIRPIN_LENGTH
    std::pair<RealT,RealT> score_hairpin_length_at_least[D_MAX_HAIRPIN_LENGTH+1];
    std::pair<RealT,RealT> cache_score_hairpin_length[D_MAX_HAIRPIN_LENGTH+1];
#endif
#if PARAMS_HAIRPIN_3_NUCLEOTIDES
    std::pair<RealT,RealT> score_hairpin_3_nucleotides[M+1][M+1][M+1];
#endif
#if PARAMS_HAIRPIN_4_NUCLEOTIDES
    std::pair<RealT,RealT> score_hairpin_4_nucleotides[M+1][M+1][M+1][M+1];
#endif
#if PARAMS_HELIX_LENGTH
    std::pair<RealT,RealT> score_helix_length_at_least[D_MAX_HELIX_LENGTH+1];
    std::pair<RealT,RealT> cache_score_helix_length[D_MAX_HELIX_LENGTH+1];
#endif
#if PARAMS_ISOLATED_BASE_PAIR
    std::pair<RealT,RealT> score_isolated_base_pair;
#endif
#if PARAMS_INTERNAL_EXPLICIT
    std::pair<RealT,RealT> score_internal_explicit[D_MAX_INTERNAL_EXPLICIT_LENGTH+1][D_MAX_INTERNAL_EXPLICIT_LENGTH+1];
#endif
#if PARAMS_BULGE_LENGTH
    std::pair<RealT,RealT> score_bulge_length_at_least[D_MAX_BULGE_LENGTH+1];
#endif
#if PARAMS_INTERNAL_LENGTH
    std::pair<RealT,RealT> score_internal_length_at_least[D_MAX_INTERNAL_LENGTH+1];
#endif
#if PARAMS_INTERNAL_SYMMETRY
    std::pair<RealT,RealT> score_internal_symmetric_length_at_least[D_MAX_INTERNAL_SYMMETRIC_LENGTH+1];
#endif
#if PARAMS_INTERNAL_ASYMMETRY
    std::pair<RealT,RealT> score_internal_asymmetry_at_least[D_MAX_INTERNAL_ASYMMETRY+1];
#endif
#if PARAMS_BULGE_0x1_NUCLEOTIDES
    std::pair<RealT,RealT> score_bulge_0x1_nucleotides[M+1];
    std::pair<RealT,RealT> score_bulge_1x0_nucleotides[M+1];
#endif
#if PARAMS_BULGE_0x2_NUCLEOTIDES
    std::pair<RealT,RealT> score_bulge_0x2_nucleotides[M+1][M+1];
    std::pair<RealT,RealT> score_bulge_2x0_nucleotides[M+1][M+1];
#endif
#if PARAMS_BULGE_0x3_NUCLEOTIDES
    std::pair<RealT,RealT> score_bulge_0x3_nucleotides[M+1][M+1][M+1];
    std::pair<RealT,RealT> score_bulge_3x0_nucleotides[M+1][M+1][M+1];
#endif
#if PARAMS_INTERNAL_1x1_NUCLEOTIDES
    std::pair<RealT,RealT> score_internal_1x1_nucleotides[M+1][M+1];
#endif
#if PARAMS_INTERNAL_1x2_NUCLEOTIDES
    std::pair<RealT,RealT> score_internal_1x2_nucleotides[M+1][M+1][M+1];
    std::pair<RealT,RealT> score_internal_2x1_nucleotides[M+1][M+1][M+1];
#endif
#if PARAMS_INTERNAL_2x2_NUCLEOTIDES
    std::pair<RealT,RealT> score_internal_2x2_nucleotides[M+1][M+1][M+1][M+1];
#endif
#if PARAMS_HELIX_STACKING
    std::pair<RealT,RealT> score_helix_stacking[M+1][M+1][M+1][M+1];
#endif
#if PARAMS_HELIX_CLOSING
    std::pair<RealT,RealT> score_helix_closing[M+1][M+1];
#endif
#if PARAMS_MULTI_LENGTH
    std::pair<RealT,RealT> score_multi_base;
    std::pair<RealT,RealT> score_multi_unpaired;
    std::pair<RealT,RealT> score_multi_paired;
#endif
#if PARAMS_DANGLE
    std::pair<RealT,RealT> score_dangle_left[M+1][M+1][M+1];
    std::pair<RealT,RealT> score_dangle_right[M+1][M+1][M+1];
#endif
#if PARAMS_EXTERNAL_LENGTH
    std::pair<RealT,RealT> score_external_unpaired;
    std::pair<RealT,RealT> score_external_paired;
#endif
    
    // cache
    std::pair<RealT,RealT> cache_score_single[C_MAX_SINGLE_LENGTH+1][C_MAX_SINGLE_LENGTH+1];
#if ( PARAMS_HELIX_LENGTH || PARAMS_ISOLATED_BASE_PAIR ) && FAST_HELIX_LENGTHS 
    std::vector<std::pair<RealT,RealT> > cache_score_helix_sums;
#endif

    void FillScores(typename std::vector<std::pair<RealT, RealT> >::iterator begin, typename std::vector<std::pair<RealT, RealT> >::iterator end, RealT value);
    void FillCounts(typename std::vector<std::pair<RealT, RealT> >::iterator begin, typename std::vector<std::pair<RealT, RealT> >::iterator end, RealT value);
    bool IsComplementary(int i, int j) const;
    
    RealT LoopScore(int i, int j, int p, int q) const;

    std::vector<RealT> GetCounts();
    void ClearCounts();
    void InitializeCache();
    void FinalizeCounts();
    
public:
    // constructor and destructor
    DuplexEngine(bool allow_noncomplementary);
    ~DuplexEngine();

    // register params with the parameter manager
    void RegisterParameters(ParameterManager<RealT> &parameter_manager);
                            
    // load sequence
    void LoadSequence(const SStruct &str1, const SStruct& str2);
    
    // load parameter values                        
    void LoadValues(const std::vector<RealT> &values);

    // MEA inference
    RealT ComputeInside();
    RealT ComputeOutside();
    RealT ComputeLogPartitionCoefficient() const { return LogPartFunc; }
    void ComputePosterior();

    RealT *GetPosterior(const RealT posterior_cutoff) const;
    RealT *GetPosterior(const RealT posterior_cutoff, std::vector<RealT>& p) const;
    RealT *GetPosterior2(const RealT posterior_cutoff) const;
    RealT *GetPosterior2(const RealT posterior_cutoff, std::vector<RealT>& p) const;

    int GetOffset(int i) { return offset[i]; }
};

#include "DuplexEngine.ipp"

#endif
