#include "mylib/common.hpp"
#include "mylib/graph.hpp"

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cassert>
#include <chrono>
#include <unordered_map>
#include <unordered_set>

void print_usage(char *fil){
  fprintf(stderr, "Usage: %s [graph_file] [probability_file] [source_vertices_file] [order_file]\n", fil);
}

using State = std::pair<uint64_t, std::vector<int8_t>>;

namespace std{
  template <>
  struct hash<State>{
    public:
    uint64_t operator()(const State& s) const{
      uint64_t h = FNV_OFFSET_BASIS_64;
      h = FNV_PRIME_64* h ^ s.first;
      for(const auto& v : s.second) h = FNV_PRIME_64* h ^ v;
      return h;
    }
  };
}

class DPBlock{
public:
  int level;
  double p = 0.0;
  std::vector<double> q;
  int8_t cnum;
  size_t lo;
  size_t hi;
  std::vector<int8_t> vlo;
  std::vector<int8_t> vhi;
  
  DPBlock() {};
  DPBlock(int _level, int8_t _siz): level(_level), cnum(_siz), q(_siz), vlo(_siz), vhi(_siz) {};
};

int main(int argc, char **argv){
  if(argc < 5){
    fprintf(stderr, "ERROR: too few arguments.\n");
    print_usage(argv[0]);
    exit(EXIT_FAILURE);
  }
  
  Graph G;
  int n, m;
  std::unordered_set<int> srcs;
  std::vector<double> pi;
  
  {
    Graph H;
    if(!H.readfromFile(argv[1])){
      fprintf(stderr, "ERROR: reading graph file %s failed.\n", argv[1]);
      print_usage(argv[0]);
      exit(EXIT_FAILURE);
    }
    
    n = H.numV();
    m = H.numE();
    
    std::vector<double> prob(m);
    pi.resize(m);
    {
      FILE *fp;
      if((fp = fopen(argv[2], "r")) == NULL){
        fprintf(stderr, "ERROR: reading probability file %s failed.\n", argv[2]);
        print_usage(argv[0]);
        exit(EXIT_FAILURE);
      }
      
      for(size_t i=0; i<m; ++i){
        fscanf(fp, "%lf", &prob[i]);
      }
      fclose(fp);
    }
    
    {
      FILE *fp;
      if((fp = fopen(argv[3], "r")) == NULL){
        fprintf(stderr, "ERROR: reading source vertices file %s failed.\n", argv[3]);
        print_usage(argv[0]);
        exit(EXIT_FAILURE);
      }
      
      int src;
      while(fscanf(fp, "%d", &src) != EOF){
        srcs.emplace(src);
      }
      fclose(fp);
    }
    
    if(!G.readfromFile(argv[4])){
      fprintf(stderr, "ERROR: reading order file %s failed.\n", argv[4]);
      print_usage(argv[0]);
      exit(EXIT_FAILURE);
    }
    
    for(size_t i=0; i<m; ++i){
      pi[i] = prob[H.etovar(G.e[i].first, G.e[i].second)];
    }
  }
  
  auto cstart = std::chrono::system_clock::now();
  
  // compute frontiers and src positions
  std::vector<std::vector<int>> fro(m+1);
  std::vector<std::vector<int>> fro_med(m);
  std::vector<std::vector<int>> srclist(m);
  size_t src_final = 0;
  {
    std::vector<int> deg(n+1);
    for(size_t i=0; i<m; ++i){
      ++deg[G.e[i].first];
      ++deg[G.e[i].second];
      auto it = srcs.find(G.e[i].first);
      if(it != srcs.end()){
        srclist[i].emplace_back(*it);
        srcs.erase(it);
        src_final = i;
      }
      it = srcs.find(G.e[i].second);
      if(it != srcs.end()){
        srclist[i].emplace_back(*it);
        srcs.erase(it);
        src_final = i;
      }
    }
    
    std::vector<int> deg2(n+1);
    for(size_t i=0; i<m; ++i){
      fro_med[i].reserve(n);
      fro[i+1].reserve(n);
      ++deg2[G.e[i].first];
      ++deg2[G.e[i].second];
      for(size_t j=1; j<=n; ++j){
        if(deg[j] && deg2[j]) fro_med[i].emplace_back(j);
      }
      --deg[G.e[i].first];
      --deg[G.e[i].second];
      for(size_t j=1; j<=n; ++j){
        if(deg[j] && deg2[j]) fro[i+1].emplace_back(j);
      }
      fro_med[i].shrink_to_fit();
      fro[i+1].shrink_to_fit();
    }
  }
  
  // compute correspondence of adjacent frontiers
  std::vector<std::vector<int>> med_to_prev(m); // correspondence of positions from fro_med[i] to fro[i]
  std::vector<std::vector<int>> prev_to_med(m); // inverse of above arrays
  std::vector<std::vector<int>> next_to_med(m); // correspondence of positions from fro[i+1] to fro_med[i]
  std::vector<std::vector<int>> med_to_next(m); // inverse of above arrays
  std::vector<std::pair<int, int>> e_pos(m);    // edge e[i] positions within fro_med[i]
  for(int i=0; i<m; ++i){
    // compute correspondence between fro_med[i] and fro[i]
    {
      size_t kk = fro[i].size();
      size_t ll = fro_med[i].size();
      
      med_to_prev[i].resize(ll);
      prev_to_med[i].resize(kk);
      std::fill(med_to_prev[i].begin(), med_to_prev[i].end(), -1);
      
      size_t l = 0;
      for(size_t k=0; k<kk; ++k){
        while(fro[i][k] != fro_med[i][l]){
          ++l;
        }
        med_to_prev[i][l] = k;
        prev_to_med[i][k] = l;
      }
    }
    // compute correspondence between fro[i+1] and fro_med[i]
    {
      size_t kk = fro_med[i].size();
      size_t ll = fro[i+1].size();
      
      next_to_med[i].resize(ll);
      med_to_next[i].resize(kk);
      std::fill(med_to_next[i].begin(), med_to_next[i].end(), -1);
      
      size_t k = 0;
      for(size_t l=0; l<ll; ++l){
        while(fro_med[i][k] != fro[i+1][l]){
          ++k;
        }
        next_to_med[i][l] = k;
        med_to_next[i][k] = l;
      }
    }
    // compute edge positions within fro_med[i]
    {
      size_t kk = fro_med[i].size();
      for(size_t k=0; k<kk; ++k){
        if(G.e[i].first == fro_med[i][k]) e_pos[i].first = k;
        if(G.e[i].second== fro_med[i][k]) e_pos[i].second= k;
      }
    }
  }
  
  std::vector<DPBlock> dp;
  std::vector<std::unordered_map<State, size_t>> maps(m+1);
  size_t snum = 2;
  
  // id=0: dummy node (sink node)
  dp.emplace_back(m, 2);
  dp[0].q[0] = 0.0;
  dp[0].q[1] = 1.0;
  
  // id=1: root node
  {
    State root;
    root.first = 0ULL;
    maps[0].emplace(root, 1);
    dp.emplace_back(0, 0);
    dp[1].p = 1.0;
  }
  
  for(size_t i=0; i<m; ++i){
    size_t kk = fro[i].size();
    size_t tt = fro_med[i].size();
    size_t ll = fro[i+1].size();
    for(const auto& ent : maps[i]){
      size_t now_id = ent.second;
      const State& now_state = ent.first;
      
      // generate intermediate state
      State med_state;
      med_state.second.resize(tt);
      int8_t cc = dp[now_id].cnum;
      
      for(size_t t=0; t<tt; ++t){
        med_state.second[t] = med_to_prev[i][t] >= 0 ? now_state.second[med_to_prev[i][t]] : cc++;
      }
      med_state.first = now_state.first;
      // introduce asterisk when src enters the (intermediate-)frontier
      for(const auto& src_now : srclist[i]){
        if(src_now == G.e[i].first) med_state.first |= 1ULL << med_state.second[e_pos[i].first];
        else                        med_state.first |= 1ULL << med_state.second[e_pos[i].second];
      }
      
      // lo_state processing
      {
        State lo_state;
        lo_state.second.resize(ll);
        std::vector<int8_t> renum(cc, -1);
        // generate lo_state
        int8_t cc_new = 0;
        for(size_t l=0; l<ll; ++l){
          int8_t c_tmp = med_state.second[next_to_med[i][l]];
          if(renum[c_tmp] < 0) renum[c_tmp] = cc_new++;
          lo_state.second[l] = renum[c_tmp];
        }
        // asterisk renumbering
        bool prune = false;
        lo_state.first = 0ULL;
        {
          uint64_t astset = med_state.first;
          while(astset){
            uint64_t ast_c2ton = astset & (-astset);
            astset ^= ast_c2ton;
            int8_t ast_c = renum[log2ton(ast_c2ton)];
            if(ast_c < 0){  // asterisk component leaves without being conncted to other components
              prune = true; // if so, prune this state
              break;
            }
            lo_state.first |= 1ULL << ast_c;
          }
        }
        
        size_t lo_id = 0;
        if(prune){ // prune lo_state if prune = true
          dp[now_id].lo = lo_id;
          // if so, computing vlo is performed by the comparison between now_state and med_state
          for(size_t k=0; k<kk; ++k){
            int med_tmp = prev_to_med[i][k];
            dp[now_id].vlo[now_state.second[k]] = (i >= src_final && med_state.first == 1ULL << med_state.second[med_tmp]) ? 1 : 0;
          }
        }else{ // not pruned
          auto it = maps[i+1].find(lo_state);
          
          if(it != maps[i+1].end()){ // maps[i+1] has already had entry
            lo_id = it->second;
          }else{                     // there is no entry
            maps[i+1].emplace(lo_state, snum);
            lo_id = snum++;
            dp.emplace_back(i+1, cc_new);
          }
          dp[now_id].lo = lo_id;
          // if not pruned, computing vlo is performed by the comparison between now_state and med_state, then renumber
          for(size_t k=0; k<kk; ++k){
            dp[now_id].vlo[now_state.second[k]] = med_state.second[prev_to_med[i][k]];
          }
          for(auto&& val : dp[now_id].vlo){
            val = renum[val];
          }
        }
      }
      
      // generate hi-side intermediate state
      int8_t cat_to   = med_state.second[e_pos[i].first];
      int8_t cat_from = med_state.second[e_pos[i].second];
      for(auto&& val : med_state.second){
        val = (val != cat_from) ? val : cat_to;
      }
      if(med_state.first & (1ULL << cat_from)){
        med_state.first ^= 1ULL << cat_from;
        med_state.first |= 1ULL << cat_to;
      }
      
      // hi_state processing
      {
        State hi_state;
        hi_state.second.resize(ll);
        std::vector<int8_t> renum(cc, -1);
        
        // generate hi_state
        int8_t cc_new = 0;
        for(size_t l=0; l<ll; ++l){
          int8_t c_tmp = med_state.second[next_to_med[i][l]];
          if(renum[c_tmp] < 0) renum[c_tmp] = cc_new++;
          hi_state.second[l] = renum[c_tmp];
        }
        // asterisk renumbering
        bool prune = false;
        hi_state.first = 0ULL;
        {
          uint64_t astset = med_state.first;
          while(astset){
            uint64_t ast_c2ton = astset & (-astset);
            astset ^= ast_c2ton;
            int8_t ast_c = renum[log2ton(ast_c2ton)];
            if(ast_c < 0){  // asterisk component leaves without being conncted to other components
              prune = true; // if so, prune this state
              break;
            }
            hi_state.first |= 1ULL << ast_c;
          }
        }
        
        size_t hi_id = 0;
        if(prune){ // prune hi_state if prune = true
          dp[now_id].hi = hi_id;
          // if so, computing vhi is performed by the comparison between now_state and med_state
          for(size_t k=0; k<kk; ++k){
            int med_tmp = prev_to_med[i][k];
            dp[now_id].vhi[now_state.second[k]] = (i >= src_final && med_state.first == 1ULL << med_state.second[med_tmp]) ? 1 : 0;
          }
        }else{ // not pruned
          auto it = maps[i+1].find(hi_state);
          if(it != maps[i+1].end()){ // maps[i+1] has already had entry
            hi_id = it->second;
          }else{                     // there is no entry
            maps[i+1].emplace(hi_state, snum);
            hi_id = snum++;
            dp.emplace_back(i+1, cc_new);
          }
          dp[now_id].hi = hi_id;
          // if not pruned, computing vhi is performed by the comparison between now_state and hi_state
          for(size_t k=0; k<kk; ++k){
            dp[now_id].vhi[now_state.second[k]] = med_state.second[prev_to_med[i][k]];
          }
          for(auto&& val : dp[now_id].vhi){
            val = renum[val];
          }
        }
      }
    }
  }
  
  // dp calculation
  for(size_t ind=1; ind<snum; ++ind){
    DPBlock& dp_now = dp[ind];
    dp[dp_now.lo].p += (1.0 - pi[dp_now.level]) * dp_now.p;
    dp[dp_now.hi].p += pi[dp_now.level]         * dp_now.p;
  }
  for(size_t ind=snum-1; ind>=2; --ind){
    DPBlock& dp_now = dp[ind];
    for(size_t c=0; c<dp_now.cnum; ++c){
      dp_now.q[c]  = dp_now.vlo[c] >= 0 ? (1.0 - pi[dp_now.level]) * dp[dp_now.lo].q[dp_now.vlo[c]] : 0.0;
      dp_now.q[c] += dp_now.vhi[c] >= 0 ? pi[dp_now.level]         * dp[dp_now.hi].q[dp_now.vhi[c]] : 0.0;
    }
  }
  
  // levelwise calculation
  for(size_t i=1; i<m; ++i){
    printf("LEVEL %zu:\n", i);
    size_t kk = fro[i].size();
    
    std::vector<double> res(kk);
    
    for(const auto& ent : maps[i]){
      size_t now_id = ent.second;
      const State& now_state = ent.first;
      
      for(size_t k=0; k<kk; ++k){
        res[k] += dp[now_id].p * dp[now_id].q[now_state.second[k]];
      }
    }
    
    for(size_t k=0; k<kk; ++k){
      printf("%d : %.15lf\n", fro[i][k], res[k]);
    }
  }
  
  auto cend = std::chrono::system_clock::now();
  double ctime = std::chrono::duration_cast<std::chrono::milliseconds>(cend-cstart).count();
  
  fprintf(stderr, "calc time: %.6lf ms\n", ctime);
  fprintf(stderr, "#(states): %zu\n", snum);
  
  return 0;
}