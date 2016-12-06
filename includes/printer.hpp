#ifndef CMR_PRINTER_H
#define CMR_PRINTER_H

#include <iostream>

#include "util.hpp"

class CMR_Printer {
 public:
  CMR_Printer(std::vector<int> & _tour_nodes, std::vector<int> & _tour_edges,
	       std::vector<double> & _lp_edges, std::vector<CMR::Edge> & _edges) :
  tour_nodes(_tour_nodes), tour_edges(_tour_edges),
    m_lp_edges(_lp_edges), edges(_edges) {}

  void print_edge(const int i){
    std::cout << "Edge number " << i
	      << "[" << edges[i].end[0] << ", " << edges[i].end[1] << "] ";
  };
  
  void best_tour_nodes(){
    for(int i = 0; i < tour_nodes.size(); i++)
      std::cout << "tour_nodes[" << i << "] = "
		<< tour_nodes[i] << std::endl;
    };
  
  void best_tour_edges(){
    for(int i = 0; i < tour_edges.size(); i++)
      if(tour_edges[i] == 1){
	print_edge(i);
	std::cout << std::endl;
      }
  };
  
  void lp_edges(){
    for(int i = 0; i < m_lp_edges.size(); i++)
      if(m_lp_edges[i] > CMR::Epsilon::Zero){
	print_edge(i);
	std::cout << " LP val: " << m_lp_edges[i] << std::endl;
      }
  };

  void pivot(const CMR::LP::PivType piv_stat){
    std::cout << "Pivot status: ";
    switch(piv_stat){
    case(CMR::LP::PivType::Frac):
      std::cout << "Fractional\n";
      break;
    case(CMR::LP::PivType::Subtour):
      std::cout << "Integral subtour\n";
      break;
    case(CMR::LP::PivType::Tour):
      std::cout << "New tour\n";
      break;
    case(CMR::LP::PivType::FathomedTour):
      std::cout << "Tour fathomed optimal\n";
      break;
    }
  };

 private:
  std::vector<int> &tour_nodes;
  std::vector<int> &tour_edges;
  std::vector<double> &m_lp_edges;
  std::vector<CMR::Edge> &edges;
};


#endif
