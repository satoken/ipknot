/*
 *bpp.cpp*
 The main code for base pair probability calculation.

 author: He Zhang
 created by: 04/2019
*/

#include <stdio.h> 
#include <sys/time.h>
#include "LinearPartition.h"
#include "Utils/utility.h"
#include "Utils/utility_v.h"

#define SPECIAL_HP

using namespace std;
using namespace LinearPartition;

void BeamCKYParser::output_to_file(string file_name, const char * type) {

    if(!file_name.empty()) {
        printf("Outputing base pairing probability matrix to %s...\n", file_name.c_str()); 
        FILE *fptr = fopen(file_name.c_str(), type); 
        if (fptr == NULL) { 
            printf("Could not open file!\n"); 
            return; 
        }
        int turn = no_sharp_turn?3:0;
#if 0
        int turn = no_sharp_turn?3:0;
        for (int i = 1; i <= seq_length; i++) {
            for (int j = i + turn + 1; j <= seq_length; j++) {
                pair<int, int> key = make_pair(i,j);
                auto got = Pij.find(key);
                if (got != Pij.end()){
                    fprintf(fptr, "%d %d %.4e\n", i, j, got->second);
                }
            }
        }
#else
        for (int i=1; i<=seq_length; i++)
            for (const auto [j, p]: Pij[i]) 
                if (i+turn<j)
                    fprintf(fptr, "%d %d %.4e\n", i, j, p);
#endif
        fprintf(fptr, "\n");
        fclose(fptr); 
        printf("Done!\n"); 
    }

    return;
}

void BeamCKYParser::get_posterior(vector<float>& bp, vector<int>& offset) const 
{
  const int L = seq_length;
  bp.resize((L+1)*(L+2)/2);
  std::fill(bp.begin(), bp.end(), 0.0);
  offset.resize(L+1);
  for (int i=0; i<=L; ++i)
    offset[i] = i*((L+1)+(L+1)-i-1)/2;

  int turn = no_sharp_turn?3:0;
  for (int i=1; i<=L; ++i)
#if 0
    for (int j=i+turn+1; j<=L; ++j) {
        pair<int, int> key = make_pair(i,j);
        auto got = Pij.find(key);
        if (got != Pij.end()){
            //fprintf(fptr, "%d %d %.4e\n", i, j, got->second);
            bp[offset[i]+j] = got->second;
        }
    }
#else
    for (const auto [j, p]: Pij[i])
        if (i+turn<j)
            bp[offset[i]+j] = p;
#endif
}

void BeamCKYParser::get_posterior(vector<vector<pair<unsigned int, float>>>& bp) const 
{
    bp = Pij;
}

void BeamCKYParser::cal_PairProb(State& viterbi) {

    Pij.resize(seq_length+1);

    for(int j=0; j<seq_length; j++){
        for(auto &item : bestP[j]){
            int i = item.first;
            State state = item.second;
            
            float temp_prob_inside = state.alpha + state.beta - viterbi.alpha;
            if (temp_prob_inside > float(-9.91152)) {
                float prob = Fast_Exp(temp_prob_inside);
                if(prob > 1.0) prob = 1.0;
                if(prob < bpp_cutoff) continue;
                Pij[i+1].emplace_back(j+1, prob);
                Pij[j+1].emplace_back(i+1, prob);
            }
        }
    }
#if 0
    // -o mode: output to a single file with user specified name;
    // bpp matrices for different sequences are separated with empty lines
    if (!bpp_file.empty()){
        output_to_file(bpp_file, "a");
    } 

    // -prefix mode: output to multiple files with user specified prefix;
    else if (!bpp_file_index.empty()) {
        output_to_file(bpp_file_index, "w");
    }
#endif
    return;
}

template <bool LPV, class T>
void BeamCKYParser::outside(vector<int> next_pair[], const std::vector<int>* cons /*=NULL*/){
    typedef T value_type;
      
    struct timeval parse_starttime, parse_endtime;

    gettimeofday(&parse_starttime, NULL);

    bestC[seq_length-1].beta = 0.0;

    // from right to left
    value_type newscore;
    for(int j = seq_length-1; j > 0; --j) {
        int nucj = nucs[j];
        int nucj1 = (j+1) < seq_length ? nucs[j+1] : -1;

        unordered_map<int, State>& beamstepH = bestH[j];
        unordered_map<int, State>& beamstepMulti = bestMulti[j];
        unordered_map<int, State>& beamstepP = bestP[j];
        unordered_map<int, State>& beamstepM2 = bestM2[j];
        unordered_map<int, State>& beamstepM = bestM[j];
        State& beamstepC = bestC[j];

        // beam of C
        {
            // C = C + U
            if (j < seq_length-1) {
                if (! (use_constraints && !allow_unpaired_position[j+1])) {
                    if (LPV) {
                        Fast_LogPlusEquals(beamstepC.beta, (bestC[j+1].beta));
                    } else {
                        newscore = score_external_unpaired(j+1, j+1);
                        Fast_LogPlusEquals(beamstepC.beta, bestC[j+1].beta + newscore);
                    }
                }
            }
        }
    
        // beam of M
        {
            for(auto& item : beamstepM) {
                int i = item.first;
                State& state = item.second;
                if (j < seq_length-1) {
                    if (use_constraints && !allow_unpaired_position[j+1])
                        continue;
                    if (LPV) {
                        Fast_LogPlusEquals(state.beta, bestM[j+1][i].beta);
                    } else {
                        newscore = score_multi_unpaired(j + 1, j + 1);
                        Fast_LogPlusEquals(state.beta, bestM[j+1][i].beta + newscore);
                    }
                }
            }
        }

        // beam of M2
        {
            for(auto& item : beamstepM2) {
                int i = item.first;
                State& state = item.second;

                // 2. M = M2
                Fast_LogPlusEquals(state.beta, beamstepM[i].beta);

                // 1. multi-loop
                {
                    for (int p = i-1; p >= std::max(i - SINGLE_MAX_LEN, 0); --p) {
                        int nucp = nucs[p];
                        int q = next_pair[nucp][j];

                        if (q != -1 && ((i - p - 1) <= SINGLE_MAX_LEN)) {

                            if (use_constraints){
                                if (p < i - 1 && !allow_unpaired_position[p+1])
                                    break;
                                if (!allow_unpaired_position[p]){
                                    q = (*cons)[p];
                                    if (q < p) break;
                                }
                                if (q > j+1 && q > allow_unpaired_range[j])
                                    continue;
                                int nucq = nucs[q];
                                if (!allow_paired(p, q, cons, nucp, nucq))
                                    continue;
                            }

                            if (LPV) {
                                Fast_LogPlusEquals(state.beta, bestMulti[q][p].beta);
                            } else {
                                newscore = score_multi_unpaired(p+1, i-1) +
                                        score_multi_unpaired(j+1, q-1);
                                Fast_LogPlusEquals(state.beta, bestMulti[q][p].beta + newscore);
                            }
                        }
                    }
                }
            }
        }

        // beam of P
        {  
            for(auto& item : beamstepP) {
                int i = item.first;
                State& state = item.second;
                int nuci = nucs[i];
                int nuci_1 = (i-1>-1) ? nucs[i-1] : -1;

                // 2. M = P
                if(i > 0 && j < seq_length-1){
                    if (LPV) {
                        int score_M1 = - v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
                        Fast_LogPlusEquals(state.beta, (beamstepM[i].beta + score_M1/kT));
                    } else {
                        newscore = score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
                        Fast_LogPlusEquals(state.beta, beamstepM[i].beta + newscore);
                    }
                }

                // 3. M2 = M + P
                int k = i - 1;
                if ( k > 0 && !bestM[k].empty()) {
                    float m1_alpha;
                    if (LPV) {
                        int M1_score = - v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
                        m1_alpha = M1_score/kT;
                    } else {
                        newscore = score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
                        m1_alpha = newscore;
                    }
                    float m1_plus_P_alpha = state.alpha + m1_alpha;
                    for (auto &m : bestM[k]) {
                        int newi = m.first;
                        State& m_state = m.second;
                        Fast_LogPlusEquals(state.beta, (beamstepM2[newi].beta + m_state.alpha + m1_alpha));
                        Fast_LogPlusEquals(m_state.beta, (beamstepM2[newi].beta + m1_plus_P_alpha));
                    }
                }

                // 4. C = C + P
                {
                    int k = i - 1;
                    if (k >= 0) {
                        int nuck = nuci_1;
                        int nuck1 = nuci;
                        float external_paired_alpha_plus_beamstepC_beta;
                        if (LPV) {
                            int score_external_paired = - v_score_external_paired(k+1, j, nuck, nuck1,
                                                                    nucj, nucj1, seq_length);
                            external_paired_alpha_plus_beamstepC_beta = beamstepC.beta + score_external_paired/kT;
                        } else {
                            newscore = score_external_paired(k+1, j, nuck, nuck1, nucj, nucj1, seq_length);
                            external_paired_alpha_plus_beamstepC_beta = beamstepC.beta + newscore;
                        }
                        Fast_LogPlusEquals(bestC[k].beta, state.alpha + external_paired_alpha_plus_beamstepC_beta);
                        Fast_LogPlusEquals(state.beta, bestC[k].alpha + external_paired_alpha_plus_beamstepC_beta);
                    } else {
                        // value_type newscore;
                        if (LPV) {
                            int score_external_paired = - v_score_external_paired(0, j, -1, nucs[0],
                                                                    nucj, nucj1, seq_length);
                            Fast_LogPlusEquals(state.beta, (beamstepC.beta + score_external_paired/kT));
                        } else {
                            newscore = score_external_paired(0, j, -1, nucs[0],
                                                                nucj, nucj1, seq_length);
                            Fast_LogPlusEquals(state.beta, beamstepC.beta + newscore);
                        }
                    }
                }

                // 1. generate new helix / single_branch
                // new state is of shape p..i..j..q
                if (i >0 && j<seq_length-1) {
                    value_type precomputed; 
                    if (!LPV)
                        precomputed = score_junction_B(j, i, nucj, nucj1, nuci_1, nuci);
                    for (int p = i - 1; p >= std::max(i - SINGLE_MAX_LEN, 0); --p) {
                        int nucp = nucs[p];
                        int nucp1 = nucs[p + 1]; 
                        int q = next_pair[nucp][j];

                        if (use_constraints){
                            if (p < i-1 && !allow_unpaired_position[p+1]) // p+1 can be unpaired
                                break;
                            if (!allow_unpaired_position[p]){ // p must be paired
                                q = (*cons)[p];
                                if (q < p) break; // p is )
                            }
                        }

                        while (q != -1 && ((i - p) + (q - j) - 2 <= SINGLE_MAX_LEN)) {
                            int nucq = nucs[q];

                            if (use_constraints){
                                if (q>j+1 && q > allow_unpaired_range[j]) // loop
                                    break;
                                if (!allow_paired(p, q, cons, nucp, nucq)) // p q is ) (
                                    break;
                            }

                            int nucq_1 = nucs[q - 1];

                            if (p == i - 1 && q == j + 1) {
                                // helix
                                if (LPV) {
                                    int score_single = -v_score_single(p,q,i,j, nucp, nucp1, nucq_1, nucq,
                                                                nuci_1, nuci, nucj, nucj1);
                                    Fast_LogPlusEquals(state.beta, (bestP[q][p].beta + score_single/kT));
                                } else {
                                    newscore = score_helix(nucp, nucp1, nucq_1, nucq);
                                    Fast_LogPlusEquals(state.beta, bestP[q][p].beta + newscore);
                                }
                            } else {
                                // single branch
                                if (LPV) {
                                    int score_single = - v_score_single(p,q,i,j, nucp, nucp1, nucq_1, nucq,
                                                    nuci_1, nuci, nucj, nucj1);
                                    Fast_LogPlusEquals(state.beta, (bestP[q][p].beta + score_single/kT));
                                } else {
                                    newscore = score_junction_B(p, q, nucp, nucp1, nucq_1, nucq) +
                                            precomputed + 
                                            score_single_without_junctionB(p, q, i, j, nuci_1, nuci, nucj, nucj1);
                                    Fast_LogPlusEquals(state.beta, bestP[q][p].beta + newscore);
                                }
                            }
                            q = next_pair[nucp][q];
                        }
                    }
                }
            }
        }

        // beam of Multi
        {
            for(auto& item : beamstepMulti) {
                int i = item.first;
                State& state = item.second;

                int nuci = nucs[i];
                int nuci1 = nucs[i+1];
                int jnext = next_pair[nuci][j];

                // 2. generate P (i, j)
                {
                    if (LPV) {
                        int score_multi = - v_score_multi(i, j, nuci, nuci1, nucs[j-1], nucj, seq_length);
                        Fast_LogPlusEquals(state.beta, (beamstepP[i].beta + score_multi/kT));
                    } else {
                        newscore = score_multi(i, j, nuci, nuci1, nucs[j-1], nucj, seq_length);
                        Fast_LogPlusEquals(state.beta, beamstepP[i].beta + newscore);
                    }
                }

                if (jnext != -1 && use_constraints){
                    int nucjnext = nucs[jnext];
                    if (jnext > allow_unpaired_range[j] || !allow_paired(i, jnext, cons, nuci, nucjnext))
                        continue;
                }

                // 1. extend (i, j) to (i, jnext)
                {
                    if (jnext != -1) {
                        if (LPV) {
                            Fast_LogPlusEquals(state.beta, (bestMulti[jnext][i].beta));
                        } else {
                            newscore = score_multi_unpaired(j, jnext - 1);
                            Fast_LogPlusEquals(state.beta, bestMulti[jnext][i].beta + newscore);
                        }
                    }
                }
            }
        }
    }  // end of for-loo j

    gettimeofday(&parse_endtime, NULL);
    double parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec-parse_starttime.tv_usec)/1000000.0;

// #ifdef testtime
    if(is_verbose) printf("Base Pairing Probabilities Calculation Time: %.2f seconds.\n", parse_elapsed_time);
// #endif


    fflush(stdout);

    return;
}

