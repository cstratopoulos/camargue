#ifndef PSEP_BBUTILS_H
#define PSEP_BBUTILS_H

class PSEP_BBNode {
 public:
  enum class NType {
    ROOT, LEFT, RIGHT
  };

  enum class NStat {
    UNVISITED, VISITED, FATHOMED
  };

 PSEP_BBNode(NType _type, int _edge) :
  node_type(_type), node_stat(NStat::UNVISITED), branch_edge(_edge) {}

  NType type() const {
    return node_type;
  }

  NStat status() const {
    return node_stat;
  }

  int edge() const {
    if(node_type != NType::ROOT)
      return branch_edge;
    else
      return -1;
  }

 private:
  NType node_type;
  NStat node_stat;
  int branch_edge;
};

#endif
