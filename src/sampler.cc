#include <mkl.h>
#include <omp.h>
#include "sampler.h"

#define BRNG VSL_BRNG_MT19937 
#define GAMMA_METHOD VSL_RNG_METHOD_GAMMA_GNORM
#define UNIFORM_METHOD VSL_RNG_METHOD_UNIFORM_STD
#define GAUSSIAN_METHOD VSL_RNG_METHOD_GAUSSIAN_ICDF 

#define MIN_PROB_VALUE -70000000

Sampler::Sampler() {
    // unsigned int SEED = time(0);
    unsigned int SEED = 0; 
    vslNewStream(&stream, BRNG,  SEED);
}

float Sampler::PLUT2BProb(const Ss& plu_top, \
        const vector<Ss>& plu_btm, \
        trie<S, S_FL>& rule_counts, S_FL& parent_counts) {
    float prob = 0;
    assert(plu_top.size() == plu_btm.size());
    for (size_t i = 0; i < plu_top.size(); ++i) {
        float ruleprob = 0;
        const S_FL& p_w = rule_counts.find(plu_btm[i]) -> data;
        if (p_w.find(plu_top[i]) == p_w.end()) {
            ruleprob = MIN_PROB_VALUE;
        }
        else {
            ruleprob = log(p_w.at(plu_top[i]) / parent_counts.at(plu_top[i]));
        }
        prob += ruleprob;
    }
    return prob;
}

float Sampler::IncreasePLUT2B(const Ss& plu_top, \
        const vector<Ss>& plu_btm, \
        trie<S, S_FL>& rule_counts, S_FL& parent_counts) {
    float prob = 0;
    assert(plu_top.size() == plu_btm.size());
    for (size_t i = 0; i < plu_top.size(); ++i) {
        assert(plu_btm[i].size()> 0);
        float& current_parent_count = parent_counts[plu_top[i]];
        float& current_rule_count = \
                 rule_counts[plu_btm[i]][plu_top[i]];
        float ruleprob = current_rule_count / current_parent_count; 
        prob += log(ruleprob);
        ++current_parent_count;
        ++current_rule_count;
    }
    return prob;
}

void Sampler::AddClusterAssignment(Datum* datum, \
                   map<symbol, Cluster*>& cluster_counter) {
    vector<Segment*> segments = datum -> segments();
    vector<Segment*>::iterator s_iter = segments.begin();
    for (; s_iter != segments.end(); ++s_iter) {
        if (((*s_iter) -> id_symbol()).string_reference() != "\'\'") {
            cluster_counter[(*s_iter) -> id_symbol()] -> Plus(*s_iter);
        }
    }
}

void Sampler::RemoveClusterAssignment(Datum* datum, \
                        map<symbol, Cluster*>& cluster_counter) {
    vector<Segment*> segments = datum -> segments();
    vector<Segment*>::iterator s_iter = segments.begin();
    for (; s_iter != segments.end(); ++s_iter) {
        if (((*s_iter) -> id_symbol()).string_reference() != "\'\'") {
            cluster_counter[(*s_iter) -> id_symbol()] -> Minus(*s_iter); 
        }
    }
    // datum -> ClearSegs();
}

Ss Sampler::SampleRHS(symbol plu_top, const float total_prob, \
        const Ss& plu_btm_candidates, \
        trie<S, vector<vector<float> > >& seg_prob_given_cluster, \
        trie<S, S_FL>& rules, \
        vector<Bound*>& bounds, \
        const int b1, const int b2) {

    if (RHS_DEBUG) {
        cerr << "In SampleRHS" << endl;
        cerr << "plu_top = " << plu_top << endl;
        cerr << "b1 = " << b1 << ", b2 = " << b2 << endl;
    }

    float random_unit_sample = rvg.GetUniformSample(); 
    float prob_threshold = total_prob + log(random_unit_sample);
    float prob_so_far = MIN_PROB_VALUE;
    
    // the size of seg_prob should be b x b
    size_t b = bounds.size();
    const int max_num_units = _config -> max_num_units();
    const int max_duration = _config -> max_duration();

    vector<int> accumulated_frame_num;
    int total_frame_num = 0;
    for (size_t i = 0 ; i < b; ++i) {
        total_frame_num += bounds[i] -> get_frame_num();
        accumulated_frame_num.push_back(total_frame_num);
    }
    int start_frame = b1 == 0 ? 0 : accumulated_frame_num[b1 - 1];
    int duration = accumulated_frame_num[b2] - start_frame;
    if (duration <= max_num_units * max_duration || b2 - b1 <= 1) {
        if (duration <= max_duration || b2 == b1) {
            vector<symbol>::const_iterator s_iter = plu_btm_candidates.begin();
            for (; s_iter != plu_btm_candidates.end(); ++s_iter) {
                Ss rhs(1, *s_iter);
                float probs; 
                trie<S, S_FL>::iterator trie_iter = rules.find(rhs);
                if (trie_iter == rules.end()) {
                    probs = MIN_PROB_VALUE;
                }
                else {
                    S_FL& parent_weight = trie_iter -> data;
                    probs = parent_weight.find(plu_top) == parent_weight.end() ? \
                            MIN_PROB_VALUE : parent_weight[plu_top];
                }
                if (RHS_DEBUG) {
                    cerr << "prob(" << plu_top << " --> " << rhs << ") = " << probs << endl;
                }
                if (probs > MIN_PROB_VALUE) {
                    probs += seg_prob_given_cluster[rhs][b1][b2]; 
                    if (RHS_DEBUG) {
                        cerr << "prob(" << plu_top << " --> b[" << b1 << "][" << b2 << "]) = " << probs << endl;
                    }
                    prob_so_far = ToolKit::SumLogs(probs, prob_so_far); 
                    if (RHS_DEBUG) {
                        cerr << "prob_so_far =  " << prob_so_far 
                            << ", prob threshold = " << prob_threshold 
                            << ", total prob = " << total_prob << endl;
                    }
                    if (prob_so_far >= prob_threshold) {
                        if (RHS_DEBUG) {
                            cerr << "returning " << rhs << endl;
                        }
                        return rhs;
                    }
                }
            }
            if (b2 - b1 >= 1) {
                vector<symbol>::const_iterator s0_iter = plu_btm_candidates.begin();
                for (; s0_iter != plu_btm_candidates.end(); ++s0_iter) {
                    vector<symbol>::const_iterator s1_iter = plu_btm_candidates.begin();
                    for (; s1_iter != plu_btm_candidates.end(); ++s1_iter) {
                        Ss rhs(1, *s0_iter);
                        rhs.push_back(*s1_iter);
                        trie<S, S_FL>::iterator trie_iter = rules.find(rhs);
                        float probs;
                        if (trie_iter == rules.end()) {
                            probs = MIN_PROB_VALUE;
                        }
                        else {
                            S_FL& parent_weight = trie_iter -> data;
                            probs = parent_weight.find(plu_top) == parent_weight.end() ? \
                                    MIN_PROB_VALUE : parent_weight[plu_top];
                        }
                        if (RHS_DEBUG) {
                            cerr << "prob(" << plu_top << " --> " << rhs << ") = " << probs << endl;
                        }
                        if (probs > MIN_PROB_VALUE) {
                            trie<S, vector<vector<float> > >::iterator c_iter = seg_prob_given_cluster.find(rhs);
                            if (c_iter != seg_prob_given_cluster.end()) {
                                probs += seg_prob_given_cluster[rhs][b1][b2];
                                if (RHS_DEBUG) {
                                    cerr << "prob(" << plu_top << " --> b[" << b1 << "][" << b2 << "]) = " << probs << endl;
                                }
                                prob_so_far = ToolKit::SumLogs(probs, prob_so_far);
                                if (RHS_DEBUG) {
                                    cerr << "prob_so_far =  " << prob_so_far 
                                        << ", prob threshold = " << prob_threshold 
                                        << ", total prob = " << total_prob << endl;
                                }
                                if (prob_so_far >= prob_threshold) {
                                    if (RHS_DEBUG) {
                                        cerr << "returning " << rhs << endl;
                                    }
                                    return rhs;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    cerr << "plu_top = " << plu_top << endl;
    cerr << "b1 = " << b1 << ", b2 = " << b2 << endl;
    cerr << "frame number = " << duration 
             << ", which is not smaller than " 
             << max_num_units * max_duration << " and b2 - b1 > 1" << endl;
    cerr << "Error in sample RHS" << endl;
    cerr << "total_prob = " << total_prob << ", probthreshold = " << prob_threshold << ", prob_so_far = " << prob_so_far << endl;
    cerr << "uni sample = " << random_unit_sample << endl;
    S matched_plu(plu_top.string_reference().substr(2));
    Ss matched_rhs(1, matched_plu);
    return matched_rhs;
}

void Sampler::SampleForward(const Ss& plu, Datum* datum, \
          const Ss& plu_btm_candidates, \
          map<symbol, vector<vector<float> > >& seg_prob, \
          trie<S, vector<vector<float> > >& prob_table,
          trie<S, S_FL>& rules, \
          vector<vector<ProbList<int> > >& B, \
          vector<vector<ProbList<int> > >& Bstar, \
          vector<Ss>& new_segmentation, \
          vector<Segment*>& proposed_segs) {

    S emp_sym("\'\'");
    Ss emp_rhs(1, emp_sym);

    vector<Bound*> bounds = datum -> bounds();
    int l = plu.size();
    int b = bounds.size();
    int i = 0, j = 0;

    while (i < l && j < b) {
        // Sample from B to decide which letter to go to
        if (DEBUG) {
            cout << "plu: " << plu[i] << endl;
            vector<float> next_i_dist = B[i][j].probs();
            for (int gg = 0; gg < (int) next_i_dist.size(); ++gg) {
                cout << "prob[" << gg + i << "]: " << next_i_dist[gg] << " ";
            }
            cout << endl;
            vector<float> log_probs = B[i][j].probs();
            ToolKit::MaxRemovedLogDist(log_probs);
            ToolKit::NormalizeDist(log_probs);
            for (int gg = 0; gg < (int) log_probs.size(); ++gg) {
                cout << "norm_prob[" << gg + i << "]: " << log_probs[gg] << " ";
            }
            cout << endl;
        }
        int next_i = B[i][j].index(\
                SampleIndexFromLogDistribution(B[i][j].probs()));
        /*
        if (DEBUG) {
            vector<float> dist_next_i = B[i][j].probs();
            cout << "B[" << i << "][" << j << "]: " << endl;
            for (unsigned int z = 0; z < dist_next_i.size(); ++z) {
                cout << z << ": " << dist_next_i[z] << endl;
            }
        }
        */
        // Sample from Bstar to decide which bound to go to
        /*
        if (DEBUG) {
            cout << "Sampling from Bstar to decide which bound to go to" << endl;
            vector<float> possible_next_j = Bstar[next_i][j].probs();
            for (int gg = 0; gg < (int) possible_next_j.size(); ++gg) {
                cout << "prob[" << gg + j << "]: " << possible_next_j[gg] << " ";
            }
            cout << endl;
        }
        */
        int next_j = Bstar[next_i][j].index(\
                SampleIndexFromLogDistribution(Bstar[next_i][j].probs())); 
        /*
        if (DEBUG) {
            cout << "next_j: " << next_j << endl;
            cout << "Done sampling next bound" << endl;
            vector<float> next_dist = Bstar[next_i][j].probs();
            for (unsigned int z = 0; z < next_dist.size(); ++z) {
                cout << "z: " <<  z << " " << next_dist[z] << ',';
            }
            cout << endl;
        }
        */
        S plu_top = plu[next_i - 1];
        Ss rhs = SampleRHS(plu_top, seg_prob[plu_top][j][next_j - 1], \
                plu_btm_candidates, prob_table, rules, bounds, j, next_j - 1); 
        if (DEBUG) { 
            cout << "next_i: " << next_i << endl;
            /*
            cout << "next_j: " << next_j << endl;
            cout << "rhs: ";
            for (size_t i = 0; i < rhs.size(); ++i) {
                cout << rhs[i] << " ";
            }
            cout << endl;
            */
        }
        // Done sampling (except segment state_seq) add information
        if (next_i > i + 1) {
            for (int k = i + 1; k < next_i; ++k) {
                Segment* emp_seg = new Segment(emp_sym);
                proposed_segs.push_back(emp_seg);
                new_segmentation[k - 1] = emp_rhs;
            }
        }
        // Create proper segments
        if (rhs.size() == 1) {
            Segment* segment = new Segment(rhs[0]);
            proposed_segs.push_back(segment);
            for (int p = j; p < next_j; ++p) {
                segment -> push_back(bounds[p]);
            }
        }
        else {
            Ss rhs_0(1, rhs[0]), rhs_1(1, rhs[1]);
            vector<float> composition_prob;
            for (int p = j; p < next_j - 1; ++p) {
                composition_prob.push_back(prob_table[rhs_0][j][p] + \
                                    prob_table[rhs_1][p + 1][next_j - 1]);
            }
            int boundary = SampleIndexFromLogDistribution(composition_prob) + j;
            Segment* segment1 = new Segment(rhs[0]);
            proposed_segs.push_back(segment1);
            for (int p = j; p <= boundary; ++p) {
                segment1 -> push_back(bounds[p]);
            }
            Segment* segment2 = new Segment(rhs[1]);
            proposed_segs.push_back(segment2);
            for (int p = boundary + 1; p < next_j; ++p) {
                segment2 -> push_back(bounds[p]);
            }
        }
        new_segmentation[next_i - 1] = rhs; 
        i = next_i;
        j = next_j;
    }
    if (i >= l && j < b) {
        cout << "l_ptr: " << i << " " << ", b_ptr: " << j << " out of (" << l << ", " << b << ")\n";
        cout << "Align not complete in Sampling forward, need to check " << datum -> tag() << endl;
        datum -> set_corrupted(true);
        datum -> set_corrupted_type(5);
        exit(-1);
    }
    else if (i < l && j >= b) {
        for (int k = i + 1; k <= l; ++k) {
            Segment* emp_seg = new Segment(emp_sym);
            proposed_segs.push_back(emp_seg);
            new_segmentation[k - 1] = emp_rhs;
        }
    }
}

void Sampler::MessageBackward(const Ss& plu, Datum* datum, \
        map<symbol, vector<vector<float> > >& segment_prob, \
        trie<S, S_FL>& rules, \
        vector<vector<ProbList<int> > >& B, 
        vector<vector<ProbList<int> > >& Bstar) {

    int l = plu.size(); 
    int b = (datum -> bounds()).size();
    vector<Bound*> bounds = datum -> bounds();

    // Initialization
    B[l][b].push_back(0, -1);
    int i = l - 1;
    int num_epsilons_in_a_row = 0;
    symbol emp("\'\'");
    Ss rhs(1, emp);

    while (num_epsilons_in_a_row < _config -> max_num_epsilons() && i > 0) {
        float probs;
        trie<S, S_FL>::iterator trie_iter = rules.find(rhs);
        if (trie_iter == rules.end()) {
                break;
        }
        else {
            S_FL& parent_weight = trie_iter -> data;
            probs = parent_weight.find(plu[i]) == parent_weight.end() ? \
                    MIN_PROB_VALUE : parent_weight[plu[i]];
            if (probs == MIN_PROB_VALUE) {
                break;
            }
            else {
                B[i][b].push_back(probs + B[i + 1][b].value(), -1);
                ++num_epsilons_in_a_row;
                --i;
            }
        }
    }

    for (; i > 0; --i) {
        B[i][b].push_back(MIN_PROB_VALUE, -1);
    }
    int max_num_units = _config -> max_num_units();
    int max_duration = _config -> max_duration();

    int total_frame_num = 0;
    vector<int> accumulated_frame_num(b, 0);
    for (int i = 0; i < b; ++i) {
        total_frame_num += bounds[i] -> get_frame_num();
        accumulated_frame_num[i] = total_frame_num;
    }
    // Compute B and Bstar
    for (int i = b - 1; i >= 0; --i) {
        for (int j = l; j >= 0; --j) {
            // Compute Bstar
            if (j > 0) {
                symbol plu_top = plu[j - 1];
                int start_frame = i == 0 ? 0 : accumulated_frame_num[i - 1];
                for (int k = i + 1; k <= b && \
                        ((accumulated_frame_num[k - 1] - start_frame) \
                         <= max_num_units * max_duration || k <= i + 2 ); ++k) { 
                    Bstar[j][i].push_back(B[j][k].value() + segment_prob[plu_top][i][k - 1], k);
                }
            }
            // Compute B
            if (j == l) {
                B[j][i].push_back(MIN_PROB_VALUE, -1);
            }
            else if ((j > 0 && i > 0) || (j == 0 && i == 0)){
                int k = j + 1;
                int num_epsilons_in_a_row = 0;
                float epsilon_prob = 0;
                while (k <= l && num_epsilons_in_a_row <= _config -> max_num_epsilons()) {
                    B[j][i].push_back(Bstar[k][i].value() + epsilon_prob, k);
                    ++num_epsilons_in_a_row;
                    trie<S, S_FL>::iterator trie_iter = rules.find(rhs);
                    if (trie_iter == rules.end()) {
                        break;
                    }
                    else {
                        S_FL& parent_weight = trie_iter -> data;
                        if (parent_weight.find(plu[k - 1]) == parent_weight.end()) {
                            break;
                        }
                        float probs = parent_weight[plu[k - 1]];
                        if (probs == MIN_PROB_VALUE) {
                            break;
                        } 
                        else {
                            epsilon_prob += probs;
                            ++k;
                        }
                    }
                }
            }
        }
    }
}

vector<symbol> Sampler::GetPluCandidates(\
        map<symbol, vector<vector<float> > >& table) {
    vector<symbol> plu_btm_candidates;
    map<symbol, vector<vector<float> > >::const_iterator \
        iter = table.begin(); 
    for (; iter != table.end(); ++iter) {
        plu_btm_candidates.push_back(iter -> first);
    }
    return plu_btm_candidates;
}

vector<symbol> Sampler::GetPluCandidates(\
        map<symbol, Cluster*>& clusters) {
    vector<symbol> plu_btm_candidates;
    map<symbol, Cluster*>::const_iterator c_iter = clusters.begin(); 
    for (; c_iter != clusters.end(); ++c_iter) {
        plu_btm_candidates.push_back(c_iter -> first);
    }
    return plu_btm_candidates;
}

void Sampler::ComputeSegProbGivenMultiCs(Datum* datum, \
        const Ss& plu_btm_candidates, \
        trie<S, S_FL>& rules, \
        trie<S, vector<vector<float> > >& prob_table) {

    vector<Bound*> bounds = datum -> bounds(); 
    size_t b = bounds.size();
    const int max_num_units = _config -> max_num_units();
    const int max_duration = _config -> max_duration();

    vector<int> accumulated_frame_num;
    int total_frame_num = 0;
    for (size_t i = 0 ; i < b; ++i) {
        total_frame_num += bounds[i] -> get_frame_num();
        accumulated_frame_num.push_back(total_frame_num);
    }
    // initialize space
    vector<vector<float> > min_prob (b, vector<float> (b, MIN_PROB_VALUE));

    vector<symbol>::const_iterator s1_iter = plu_btm_candidates.begin();
    for (; s1_iter != plu_btm_candidates.end(); ++s1_iter) {
        vector<symbol>::const_iterator s2_iter = plu_btm_candidates.begin();
        for (; s2_iter != plu_btm_candidates.end(); ++s2_iter) {
            Ss rhs(1, *s1_iter);
            rhs.push_back(*s2_iter);
            trie<S, S_FL>::iterator trie_iter = rules.find(rhs);
            if (trie_iter != rules.end() && (trie_iter -> data).size() != 0) {
                prob_table.insert(rhs, min_prob);
            }
        }
    }
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (size_t s0_ptr = 0; s0_ptr < plu_btm_candidates.size(); ++s0_ptr) {
        /*
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        */
        for (size_t s1_ptr = 0; s1_ptr < plu_btm_candidates.size(); ++s1_ptr) {
            Ss rhs(1, plu_btm_candidates[s0_ptr]);
            rhs.push_back(plu_btm_candidates[s1_ptr]);
            trie<S, vector<vector<float> > >::iterator trie_iter = prob_table.find(rhs);
            if (trie_iter != prob_table.end()) {
                vector<vector<float> >& table = (trie_iter -> data);
                Ss rhs_0(1, plu_btm_candidates[s0_ptr]), rhs_1(1, plu_btm_candidates[s1_ptr]);
                trie<S, vector<vector<float> > >::iterator trie_iter_0, trie_iter_1;
                trie_iter_0 = prob_table.find(rhs_0);
                trie_iter_1 = prob_table.find(rhs_1);
                assert(trie_iter_0 != prob_table.end());
                assert(trie_iter_1 != prob_table.end());
                const vector<vector<float> >& seg_prob_given_c_0 = trie_iter_0 -> data;
                const vector<vector<float> >& seg_prob_given_c_1 = trie_iter_1 -> data;
                /*
                #ifdef _OPENMP
                #pragma omp parallel for
                #endif
                */
                for (size_t b1 = 0; b1 < b; ++b1) {
                    int start_frame = b1 == 0 ? 0 : accumulated_frame_num[b1 - 1];
                    for (size_t b2 = b1 + 1; b2 < b; ++b2) {
                        int duration = accumulated_frame_num[b2] - start_frame;
                        if (duration <= max_num_units * max_duration || b2 - b1 == 1) {
                            for (size_t p = b1; p < b2; ++p) {
                                float probs = (seg_prob_given_c_0[b1][p] + \
                                    seg_prob_given_c_1[p + 1][b2]);
                                table[b1][b2] = ToolKit::SumLogs(probs, table[b1][b2]);
                            }
                        }
                    }
                }
            }
        }
    }
}

void Sampler::FastComputeSegProbGivenPluTop(Datum* datum, \
        const Ss& plu_btm_candidates, \
        const Ss& plu, \
        trie<S, S_FL>& rules, \
        trie<S, vector<vector<float> > >& prob_table, \
        map<symbol, vector<vector<float> > >& seg_prob) { 

    vector<Bound*> bounds = datum -> bounds();
    size_t b = bounds.size();
    const int max_num_units = _config -> max_num_units();
    const int max_duration = _config -> max_duration();

    vector<int> accumulated_frame_num;
    int total_frame_num = 0;
    for (size_t i = 0 ; i < b; ++i) {
        total_frame_num += bounds[i] -> get_frame_num();
        accumulated_frame_num.push_back(total_frame_num);
    }

    Ss plu_set;
    for (size_t i = 0; i < plu.size(); ++i) {
        if (seg_prob.find(plu[i]) == seg_prob.end()) {
            plu_set.push_back(plu[i]);
            seg_prob[plu[i]].resize(b, vector<float> (b, MIN_PROB_VALUE));
        }
    }

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (size_t i = 0; i < plu_set.size(); ++i) { 
        vector<vector<float> >& plu_top_seg_prob = seg_prob[plu_set[i]];
        for (size_t s0_ptr = 0; s0_ptr < plu_btm_candidates.size(); ++s0_ptr) {
            Ss rhs(1, plu_btm_candidates[s0_ptr]);
            trie<S, vector<vector<float> > >::iterator c_trie_iter = prob_table.find(rhs);
            assert(c_trie_iter != prob_table.end());
            vector<vector<float> >& single_c_prob_table = c_trie_iter -> data; 
            // unary rules
            float rule_prob;
            trie<S, S_FL>::iterator trie_iter = rules.find(rhs); 
            if (trie_iter == rules.end()) {
                rule_prob = MIN_PROB_VALUE;
            }
            else {
                S_FL& parent_weight = (trie_iter -> data);
                rule_prob = parent_weight.find(plu_set[i]) == parent_weight.end() ? \
                            MIN_PROB_VALUE : parent_weight[plu_set[i]];
            }
            if (rule_prob > MIN_PROB_VALUE) {
                /*
                #ifdef _OPENMP
                #pragma omp parallel for
                #endif
                */
                for (size_t b1 = 0; b1 < b; ++b1) {
                    int start_frame = b1 == 0 ? 0 : accumulated_frame_num[b1 - 1];
                    for (size_t b2 = b1; b2 < b; ++b2) {
                        int duration = accumulated_frame_num[b2] - start_frame;
                        if (duration <= max_duration || b2 - b1 == 0) {
                            float map_prob = rule_prob + single_c_prob_table[b1][b2]; 
                            plu_top_seg_prob[b1][b2] = \
                                ToolKit::SumLogs(plu_top_seg_prob[b1][b2], map_prob);
                        }
                    }
                }
            }
        }
        // binary rules
        for (size_t s0_ptr = 0; s0_ptr < plu_btm_candidates.size(); ++s0_ptr) {
            for (size_t s1_ptr = 0; s1_ptr < plu_btm_candidates.size(); ++s1_ptr) {
                vector<symbol> rhs(1, plu_btm_candidates[s0_ptr]);
                rhs.push_back(plu_btm_candidates[s1_ptr]);
                trie<S, vector<vector<float> > >::iterator c_trie_iter = prob_table.find(rhs);
                if (c_trie_iter != prob_table.end()) {
                    vector<vector<float> >& multi_c_prob_table = c_trie_iter -> data;
                    float rule_prob; 
                    trie<S, S_FL>::iterator trie_iter = rules.find(rhs);
                    if (rules.find(rhs) == rules.end()) {
                        rule_prob = MIN_PROB_VALUE;
                        cerr << "This should happen. If rules.find(rhs) == rules.end(), \
                            then prob_table.find(rhs) should be prob_table.end() too." << endl;
                        exit(-1);
                    }
                    else {
                        S_FL& parent_weight = trie_iter -> data;
                        rule_prob = parent_weight.find(plu_set[i]) == parent_weight.end() ? \
                                   MIN_PROB_VALUE : parent_weight[plu_set[i]]; 
                    }
                    if (rule_prob > MIN_PROB_VALUE) {
                        /*
                        #ifdef _OPENMP
                        #pragma omp parallel for
                        #endif
                        */
                        for (size_t b1 = 0; b1 < b; ++b1) {
                            int start_frame = b1 == 0 ? 0 : accumulated_frame_num[b1 - 1];
                            for (size_t b2 = b1 + 1; b2 < b; ++b2) {
                                int duration = accumulated_frame_num[b2] - start_frame;
                                if (duration <= max_num_units * max_duration || b2 - b1 == 1) {
                                    float map_prob = rule_prob + multi_c_prob_table[b1][b2];
                                    plu_top_seg_prob[b1][b2] = ToolKit::SumLogs(plu_top_seg_prob[b1][b2], map_prob);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    if (PLU_DEBUG) {
        map<S, vector<vector<float> > >::iterator m_iter = seg_prob.begin();
        for (; m_iter != seg_prob.end(); ++m_iter) { 
            float total_prob = MIN_PROB_VALUE;
            for (size_t s0_ptr = 0; s0_ptr < plu_btm_candidates.size(); ++s0_ptr) {
                Ss rhs(1, plu_btm_candidates[s0_ptr]);
                trie<S, S_FL>::iterator rule_iter = rules.find(rhs);
                float rule_prob;
                if (rule_iter == rules.end()) {
                    cerr << "Prob(" << m_iter -> first << " --> " << rhs << ") = " << MIN_PROB_VALUE << endl;
                    rule_prob = MIN_PROB_VALUE;
                }
                else {
                    S_FL& parent_weight = rule_iter -> data;
                    rule_prob = parent_weight.find(m_iter -> first) == parent_weight.end() ? \
                                      MIN_PROB_VALUE : parent_weight[m_iter -> first];
                    cerr << "Prob(" << m_iter -> first << " --> " << rhs << ") = " << rule_prob << endl;
                }
                trie<S, vector<vector<float> > >::iterator c_iter = prob_table.find(rhs);
                cerr << "Prob(" << rhs << " --> b[0][2]) = " << (c_iter -> data)[0][2] << endl;
                cerr << "total_prob_old = " << total_prob << " ";
                total_prob = ToolKit::SumLogs(total_prob, rule_prob + (c_iter -> data)[0][2]);
                cerr << ", total_prob_new = " << total_prob << endl;

                for (size_t s1_ptr = 0; s1_ptr < plu_btm_candidates.size(); ++ s1_ptr) {
                    rhs.assign(1, plu_btm_candidates[s0_ptr]);
                    rhs.push_back(plu_btm_candidates[s1_ptr]);
                    trie<S, S_FL>::iterator rule_iter = rules.find(rhs);
                    float rule_prob;
                    if (rule_iter == rules.end()) {
                        cerr << "Prob(" << m_iter -> first << " --> " << rhs << ") = " << MIN_PROB_VALUE << endl;
                        rule_prob = MIN_PROB_VALUE;
                    }
                    else {
                        S_FL& parent_weight = rule_iter -> data;
                        rule_prob = parent_weight.find(m_iter -> first) == parent_weight.end() ? \
                                          MIN_PROB_VALUE : parent_weight[m_iter -> first];
                        cerr << "Prob(" << m_iter -> first << " --> " << rhs << ") = " << rule_prob << endl;
                    }
                    trie<S, vector<vector<float> > >::iterator c_iter = prob_table.find(rhs);
                    cerr << "Prob(" << rhs << " --> b[0][2]) = " << (c_iter -> data)[0][2] << endl;

                    cerr << "total_prob_old = " << total_prob << " ";
                    total_prob = ToolKit::SumLogs(total_prob, rule_prob + (c_iter -> data)[0][2]);
                    cerr << ", total_prob_new = " << total_prob << endl;
                }
            }
            cerr << "Plu: " << m_iter -> first << " ";
            cerr << "Check whether the two are the same: " << (m_iter -> second)[0][2] << " " << total_prob << endl;
        }
    }
}


void Sampler::ComputeSegProbGivenCluster(Datum* datum, \
    const map<symbol, Cluster*>& clusters, \
    trie<S, vector<vector<float> > >& segment_prob_given_cluster) {

    map<symbol, Cluster*>::const_iterator iter = clusters.begin();
    for (; iter != clusters.end(); ++iter) {
        Ss id(1, (iter -> first));
        bool inserted ATTRIBUTE_UNUSED = \
        segment_prob_given_cluster.insert(id, \
                (iter -> second) -> ConstructSegProbTable(datum -> bounds())).second;
        assert(inserted);
        if (SEG_DEBUG) {
            cout << "For Cluster: " << iter -> first << endl;
            int b = (datum -> bounds()).size();
            const vector<vector<float> >& prob_table = \
                    segment_prob_given_cluster.find(id) -> data;
            for (int j = 0; j < b; ++j) {
                for (int k = 0; k < b; ++k) {
                   cout << "b[" << j << "][" << k << "]: " << \
                    prob_table[j][k] << " "; 
                }
                cout << endl;
            }
            cout << endl;
        }

    }
}

vector<Ss> Sampler::Resegment(Datum* datum, 
        trie<S, S_FL>& rules, \
        trie<S, S_FL>& rule_counts, \
        S_FL& parent_counts, \
        const vector<symbol>& plu, \
        map<symbol, Cluster*>& clusters, \
        map<symbol, Cluster*>& counters, \
        const float pi0, const float r0) {

    int l = plu.size(); 
    int b = (datum -> bounds()).size();
    vector<Ss> new_segmentation;
    new_segmentation.resize(plu.size());

    if (SEG_DEBUG) {
        cout << "removing assginment" << endl;
    }
    RemoveClusterAssignment(datum, counters);

    if (SEG_DEBUG) {
        cout << "segprob given clusters" << endl;
    }
    // Use trie to store prob table
    trie<symbol, vector<vector<float> > > seg_prob_given_Cs;

    // Compute P(b|c) 
    ComputeSegProbGivenCluster(datum, clusters, seg_prob_given_Cs);

    // build p(b|c1c2)
    vector<symbol> plu_btm_candidates = GetPluCandidates(clusters);
    if (SEG_DEBUG) {
        cout << "start building prob given multi cs: " << time(0) << endl;
    }
    ComputeSegProbGivenMultiCs(datum, plu_btm_candidates, rules, seg_prob_given_Cs);

    if (SEG_DEBUG) {
        cout << "stop building prob given multi cs: " << time(0) << endl;
    }

    if (SEG_DEBUG) {
        for (size_t s0_ptr = 0; s0_ptr < plu_btm_candidates.size(); ++s0_ptr) {
            for (size_t s1_ptr = 1; s1_ptr < plu_btm_candidates.size(); ++ s1_ptr) {
                Ss rhs(1, plu_btm_candidates[s0_ptr]);
                rhs.push_back(plu_btm_candidates[s1_ptr]);
                if (seg_prob_given_Cs.find(rhs) == seg_prob_given_Cs.end()) {
                    cerr << "Can't find SegClusterProb for " << rhs  << endl;
                    exit(-1);
                }
                vector<vector<float> >& prob_given_Cs = (seg_prob_given_Cs.find(rhs)) -> data;
                cout << "prob_" << rhs << "[" << 0 << "][" << 2 << "]: " << prob_given_Cs[0][2] << endl;
            }
        } 
    }

    // Compute p(b|l)
    if (SEG_DEBUG) {
        cout << "seg prob given plu_top " << endl;
    }
    // cout << "start fast compute: " << time(0) << endl;
    map<symbol, vector<vector<float> > > segment_prob;
    FastComputeSegProbGivenPluTop(datum, plu_btm_candidates, plu, rules, seg_prob_given_Cs, segment_prob);
    // cout << "stop fast compute: " << time(0) << endl;

    if (SEG_DEBUG) {
        cout << "message starts!" << endl;
    }
    // Compute B and Bstar
    vector<vector<ProbList<int> > > Bstar, B; 
    B.resize(l + 1); 
    Bstar.resize(l + 1); 
    for (int i = 0 ; i <= l; ++i) {
        B[i].resize(b + 1);
        Bstar[i].resize(b + 1);
    }
    if (MSG_DEBUG) {
        cout << "Message Backward" << endl;
    }
    MessageBackward(plu, datum, segment_prob, rules, B, Bstar);

    if (B_DEBUG) {
        for (int i = 0; i <= l; ++i) {
            for (int j = 0; j <= b; ++j) {
                cout << "B[" << i << "][" << j << "]: " << B[i][j].value() << " ";
            }
            cout << endl;
        }
    }
    if (B_DEBUG) {
        for (int i = 0; i <= l; ++i) {
            for (int j = 0; j <= b; ++j) {
                cout << "Bstar[" << i << "][" << j << "]: " << Bstar[i][j].value() << " ";
            }
            cout << endl;
        }
    }

    // Sample forward
    if (MSG_DEBUG) {
        cout << "Sample Forward" << endl;
    }
    vector<Segment*> proposed_segs;
    SampleForward(plu, datum, plu_btm_candidates, \
            segment_prob, seg_prob_given_Cs,
            rules, B, Bstar, new_segmentation, proposed_segs);

    // get pi1 and r1

    float r1 = PLUT2BProb(plu, new_segmentation, \
                           rule_counts, parent_counts);

    float pi1 = IncreasePLUT2B(plu, new_segmentation, \
                              rule_counts, parent_counts);

    float accept = exp((pi1 + r0) - (pi0 + r1));

    cerr << "accept ratio in Sampler = " << accept << endl;
    if (rvg.GetUniformSample() <= accept) {
        datum -> ClearSegs();
        datum -> set_segments(proposed_segs); 
    }
    else {
        vector<Segment*>::iterator s_iter = proposed_segs.begin();
        for (; s_iter != proposed_segs.end(); ++s_iter) {
            delete (*s_iter);
        }
    }

    // Sample State sequence and Mixture ID for segments
    vector<Segment*> segments = datum -> segments();
    if (MSG_DEBUG) {
        cout << "Sampling state and mixture seq" << endl;
        cout << "segment size: " << segments.size() << endl;
    }
    for (int i = 0; i < (int) segments.size(); ++i) {
        if (SEG_DEBUG) {
            cout << "Sampling: " << i << endl;
        }
        if ((segments[i] -> id_symbol()).string_reference() != "\'\'") {
            SampleStateMixtureSeq(segments[i], \
                    clusters[segments[i] -> id_symbol()]);
        }
    }
    // Add newly sampled alignment results
    if (SEG_DEBUG) {
        cout << "Adding cluster assignement" << endl;
    }
    AddClusterAssignment(datum, counters);
    return new_segmentation;
} 

int Sampler::SampleIndexFromLogDistribution(vector<float> log_probs) {
    ToolKit::MaxRemovedLogDist(log_probs);
    return SampleIndexFromDistribution(log_probs);
}

int Sampler::SampleIndexFromDistribution(vector<float> probs) {
   float sum = ToolKit::NormalizeDist(probs);
   // sample a random number from a uniform dist [0,1]
   float random_unit_sample = rvg.GetUniformSample(); 
   while (random_unit_sample == 0 || random_unit_sample == 1) {
       random_unit_sample = rvg.GetUniformSample(); 
   }
   // cout << "random_unit_sample = " << random_unit_sample << endl;
   // figure out the index 
   size_t index = 0; 
   sum = probs[index];
   while (random_unit_sample > sum) {
       if (++index < probs.size()) {
           sum += probs[index];
       }
       else {
           break;
       }
   }
   if (index >= probs.size()) {
       index = probs.size() - 1;
   }
   return index;
}

void Sampler::SampleStateMixtureSeq(Segment* segment, Cluster* cluster) {
    // Message Backward
    if (SEG_DEBUG) {
        cout << "segment id: " << segment -> id_symbol() << endl;
        cout << "data size: " << segment -> frame_num() << endl;
        cout << "message backward for a segment" << endl;
        cout << "cluster pointer: " << cluster << endl;
    }
    vector<vector<ProbList<int> > > B = cluster -> MessageBackwardForASegment(segment);
    // Sample Forward
    vector<int> state_seq;
    vector<int> mix_seq;
    int s;
    int m;
    for (int i = 0, j = 0; j < segment -> frame_num(); ++j) {
        s = B[i][j].index(SampleIndexFromLogDistribution(B[i][j].probs()));
        state_seq.push_back(s - 1);
        m = SampleIndexFromLogDistribution((cluster -> emission(s - 1)).\
            ComponentLikelihood(segment -> frame(j)));
        mix_seq.push_back(m);
        i = s;
    }
    segment -> set_state_seq(state_seq);
    segment -> set_mix_seq(mix_seq);
    // Delete B
}

void Sampler::SampleClusterParams(Cluster* cluster, Cluster* counter) {
    int state_num_ = cluster -> state_num(); 
    if (counter == NULL) {
        cerr << "counter should not be null for this version" << endl;
        exit(-1);
    }
    else {
        vector<vector<float> > trans_probs;
        vector<vector<float> > trans_counts = counter -> transition_probs();
        for (int i = 0; i < state_num_; ++i) {
            for (int j = i; j < state_num_ + 1; ++j) {
                trans_counts[i][j] += _config -> transition_alpha();
            }
            vector<float> trans_prob = SampleDirFromGamma(state_num_ + 1 - i, &trans_counts[i][i]);
            for (int j = 0; j < i; ++j) {
                trans_prob.insert(trans_prob.begin(), MIN_PROB_VALUE);
            }
            trans_probs.push_back(trans_prob);
        }
        cluster -> set_trans(trans_probs);
        for (int i = 0; i < state_num_; ++i) {
            int mix_num_ = cluster -> emission(i).mix_num(); 
            vector<float> weight_count = (counter -> emission(i)).weight();
            for (int j = 0; j < mix_num_; ++j) {
                // Sample Gaussian
                vector<float> mean_count = (counter -> emission(i)).mixture(j).mean();
                vector<float> pre_count = (counter -> emission(i)).mixture(j).pre();
                vector<float> pre = SampleGaussianPre(mean_count, pre_count, weight_count[j]);
                vector<float> mean = SampleGaussianMean(pre, mean_count, weight_count[j]);
                (cluster -> emission(i)).mixture(j).set_mean(mean);
                (cluster -> emission(i)).mixture(j).set_pre(pre);
                (cluster -> emission(i)).mixture(j).set_det();
                weight_count[j] += _config -> mix_alpha();
            }
            vector<float> weight = SampleDirFromGamma(mix_num_, &weight_count[0], -5);
            (cluster -> emission(i)).set_weight(weight);
        }
    }
}

vector<float> Sampler::SampleDirFromGamma(int n, float* alpha, float min) {
    vector<float> samples(n, 0);
    bool need_to_set_sparse = true;
    for (int i = 0; i < n; ++i) {
        if (vsRngGamma(GAMMA_METHOD, stream, 1, &samples[i], alpha[i], 0, 1) != VSL_STATUS_OK) {
           cout << "Error when calling SampleDirFromGamma" << endl;
           cout << "The parameters are:" << endl;
           cout << "an: " << alpha[i] << " bn: 1" << endl; 
           exit(1);
        }
        samples[i] = samples[i] == 0 ? min : log(samples[i]);
        if (samples[i] < min) {
            samples[i] = min;
        }
        if (samples[i] > min) {
            need_to_set_sparse = false;
        }
    }
    if (need_to_set_sparse) {
        samples[0] = 0;
    }
    float sum = ToolKit::SumLogs(samples); 
    for (int i = 0; i < n; ++i) {
        samples[i] -= sum; 
    }
    return samples;
}

vector<float> Sampler::SampleGaussianPre(vector<float> mean_count, \
        vector<float> pre_count, float n) {
    vector<float> gaussian_b0;
    vector<float> gaussian_u0;
    vector<float> new_pre(_config -> dim(), 1);
    gaussian_b0 = _config -> gaussian_b0();
    gaussian_u0 = _config -> gaussian_u0();    
    for (int i = 0; i < _config -> dim(); ++i) {
        float bn = UpdateGammaRate(\
                gaussian_b0[i], pre_count[i], mean_count[i], \
                n, gaussian_u0[i]);
        float an = _config -> gaussian_a0() + n / 2;
        if (vsRngGamma(GAMMA_METHOD, stream, 1, &new_pre[i], an, 0, 1 / bn) != VSL_STATUS_OK) {
            cout << "Error when calling SampleGaussianPre" << endl;
            cout << "The parameters are: " << endl;
            cout << "an: " << an << " bn: " << bn << endl;
            exit(1);
        } 
    }
    return new_pre;
}

float Sampler::UpdateGammaRate(float b0, float x2, float x, float n, float u0) {
    float k0 = _config -> gaussian_k0();
    if (n == 0) {
        return b0;
    }
    else {
        return b0 + 0.5 * (x2 - x * x / n) + (k0 * n * (x / n - u0) * (x / n - u0))/(2 * (k0 + n));
    }
}

vector<float> Sampler::SampleGaussianMean(vector<float> pre, \
        vector<float> count, float n) {
    vector<float> new_mean(_config -> dim(), 0);
    vector<float> gaussian_u0;
    float k0 = _config -> gaussian_k0();
    gaussian_u0 = _config -> gaussian_u0();
    for (int i = 0; i < _config -> dim(); ++i) {
        float un = (k0 * gaussian_u0[i] + count[i]) / (k0 + n);
        float kn = k0 + n;
        float std = sqrt(1 / (kn * pre[i]));
        vsRngGaussian(GAUSSIAN_METHOD, stream, 1, &new_mean[i], un, std); 
    }
    return new_mean;
}

void Sampler::InsertEmptyStrings(Datum* data) {
    string empty_str("\'\'");
    symbol empty_sym(empty_str);
    int num_epsilons_in_a_row = 0;
    vector<Segment*> segs = data -> segments();
    vector<Segment*>::iterator s_iter = segs.begin(); 

    int counter = 0;
    for (; s_iter != segs.end(); ++s_iter) {
        if (num_epsilons_in_a_row >= _config -> max_num_epsilons()) {
            num_epsilons_in_a_row = 0;
        }
        float uni = rvg.GetUniformSample();
        if (uni <= _config -> empty_string_prior() && counter != 0) {
            Segment* empty_seg = new Segment(empty_sym);
            s_iter = segs.insert(s_iter, empty_seg);
            ++num_epsilons_in_a_row;
        }
        else {
            num_epsilons_in_a_row = 0;
        }
        ++counter;
    }
    data -> set_segments(segs);
}

Sampler::~Sampler() {
    vslDeleteStream(&stream);
}
