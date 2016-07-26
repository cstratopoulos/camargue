#ifndef PSEP_BBUTILS_H
#define PSEP_BBUTILS_H

class PSEP_BBNode {
 public:
 PSEP_BBNode(PSEP_BBNode::NType _type, int _edge) :
  node_type(_type), branch_edge(_edge), node_stat(PSEP_BBNode::UNVISITED) {}
  
  enum class NType {
    ROOT, LEFT, RIGHT
  };

  NType &type() const {
    return node_type;
  }

  enum class NStat {
    UNVISITED, VISITED, FATHOMED
  };

  NStat &status() const {
    return node_stat;
  }

  int edge() const {
    if(node_type != NType::ROOT)
      return branch_edge;
    else
      return -1;
  }

 private:
  PSEP_BBNode::Ntype node_type;
  PSEP_BBNode::NStat node_stat;
  int branch_edge;
};

#endif
