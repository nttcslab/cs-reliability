#include "mylib/common.hpp"
#include "mylib/graph.hpp"

#include "tdzdd/DdSpec.hpp"
#include "tdzdd/DdStructure.hpp"
#include "tdzdd/DdEval.hpp"
#include "tdzdd/spec/FrontierBasedSearch.hpp"
#include "tdzdd/util/Graph.hpp"

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cassert>
#include <unordered_set>
#include <chrono>

class ProbEval : public tdzdd::DdEval<ProbEval, double> {
private:
  std::vector<double> prob_list_;
  
public:
  ProbEval(const std::vector<double>& prob_list) : prob_list_(prob_list) {}
  
  void evalTerminal(double& p, bool one) const { p = one ? 1.0 : 0.0; }
  
  void evalNode(double& p, int level,
                tdzdd::DdValues<double, 2> const& values) const {
    double pc = prob_list_[prob_list_.size() - level];
    p = values.get(0) * (1 - pc) + values.get(1) * pc;
  }
};

void print_usage(char *fil){
  fprintf(stderr, "Usage: %s [graph_file] [probability_file] [terminal_vertices_file] [order_file]\n", fil);
}

int main(int argc, char **argv){
  if(argc < 5){
    fprintf(stderr, "ERROR: too few arguments.\n");
    print_usage(argv[0]);
    exit(EXIT_FAILURE);
  }
  
  Graph G;
  int n, m;
  std::unordered_set<int> srcs;
  int dest;
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
      std::vector<int> srcvec;
      while(fscanf(fp, "%d", &src) != EOF){
        srcvec.emplace_back(src);
      }
      fclose(fp);
      for(size_t i=0; i<srcvec.size()-1; ++i){
        srcs.emplace(srcvec[i]);
      }
      dest = srcvec[srcvec.size()-1];
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
  
  tdzdd::Graph tG;
  
  for(const auto& edg : G.e){
    tG.addEdge(std::to_string(edg.first), std::to_string(edg.second));
  }
  tG.update();
  
  {
    auto pstart = std::chrono::system_clock::now();
    
    for(const auto& src : srcs){
      tG.setColor(std::to_string(src), 1);
    }
    tG.setColor(std::to_string(dest), 1);
    tG.update();
    
    tdzdd::FrontierBasedSearch fbs(tG, -1, false, false);
    
    tdzdd::DdStructure<2> dd(fbs);
    dd.useMultiProcessors(false);
    double rel = dd.evaluate(ProbEval(pi));
    
    printf("%d: %.15lf\n", dest, rel);
    
    tG.clearColors();
    
    auto pend = std::chrono::system_clock::now();
    double ptime = std::chrono::duration_cast<std::chrono::milliseconds>(pend-pstart).count();
    fprintf(stderr, "%d: %.6lf ms\n", dest, ptime);
  }
  
  auto cend = std::chrono::system_clock::now();
  double ctime = std::chrono::duration_cast<std::chrono::milliseconds>(cend-cstart).count();
  
  fprintf(stderr, "calc time: %.6lf ms\n", ctime);
  
  return 0;
}