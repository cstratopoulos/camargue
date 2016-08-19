#include <fstream>
		   
#include <iostream>

#include "Graph.h"

using namespace std;
using namespace PSEP;


Edge::Edge(int e0, int e1, int _len):
  len(_len),
  removable(false){
  end[0] = e0;
  end[1] = e1;
}

void Graph::print_edges() {
    for (int i = 0; i < edges.size(); ++i)
    {
        cout << edges[i].end[0] << ", " << edges[i].end[1]
	     << "-(" << edges[i].len << ")" << endl;
    }
}

int GraphUtils::connected(SupportGraph *G, int *icount, std::vector<int> &island,
		       int starting_node){
  *icount = 0;
  for(int i = 0; i < G->node_count; i++)
    G->nodelist[i].mark = 0;

  dfs(starting_node, G, icount, island);

  if(*icount == G->node_count)
    return 1;
  else
    return 0;
}

void GraphUtils::dfs(int n, SupportGraph *G, int *icount, std::vector<int> &island)
{
  int neighbor;
  SNode *pn;

  island[*icount] = n;
  (*icount)++;

  pn = &G->nodelist[n];
  pn->mark = 1;

  for(int i = 0; i < pn->s_degree; i++){
    neighbor = pn->adj_objs[i].other_end;
    if(G->nodelist[neighbor].mark == 0)
      dfs(neighbor, G, icount, island);
  }
}

void GraphUtils::get_delta (std::vector<int> &nodelist, std::vector<Edge> &elist,
			 int *deltacount_p, std::vector<int> &delta,
			 std::vector<int> &marks){
  for(int i = 0; i < nodelist.size(); i++)
    marks[nodelist[i]] = 1;

  int k = 0;
  for(int i = 0; i < elist.size(); i++)
    if(marks[elist[i].end[0]] + marks[elist[i].end[1]] == 1)
      delta[k++] = i;

  *deltacount_p = k;
  
  for(int i = 0; i < nodelist.size(); i++)
    marks[nodelist[i]] = 0;
}

void GraphUtils::get_delta (int nsize, int *nlist, int ecount, int *elist,
			 int *deltacount, int *delta, int *marks){
  int i, k = 0;

  for (i = 0; i < nsize; i++) marks[nlist[i]] = 1;

  for (i = 0; i < ecount; i++) {
    if (marks[elist[2*i]] + marks[elist[2*i+1]] == 1) {
      delta[k++] = i;
    }
  }

  *deltacount = k;

  for (i = 0; i < nsize; i++) marks[nlist[i]] = 0;
}

int GraphUtils::build_s_graph (int node_count, int edge_count,
			    vector<Edge> &edges,
			    vector<int> &support_indices,
			    vector<double> &m_lp_edges, SupportGraph *G_s)
{
  int i, ind, a, b;
  SNode *n;
  s_adjobj *p;

  if(G_s->nodelist) free(G_s->nodelist);
  if(G_s->adjspace) free(G_s->adjspace);
  G_s->nodelist = (SNode *) malloc(node_count * sizeof(SNode));
  G_s->adjspace = (s_adjobj *) malloc(2 * edge_count * sizeof(SNode));
  if(!G_s->nodelist || !G_s->adjspace){
    fprintf(stderr, "Out of memory for support nodelist or adjspace\n");
    return 1;
  }

  for(i = 0; i < node_count; i++)
    G_s->nodelist[i].s_degree = 0;

  for(i = 0; i < edge_count; i++){
    ind = support_indices[i];
    a = edges[ind].end[0]; b = edges[ind].end[1];
    G_s->nodelist[a].s_degree++;
    G_s->nodelist[b].s_degree++;
  }

  p = G_s->adjspace;
  for(i = 0; i < node_count; i++){
    G_s->nodelist[i].adj_objs = p;
    p += G_s->nodelist[i].s_degree;
    G_s->nodelist[i].s_degree = 0;
  }

  for(i = 0; i < edge_count; i++){
    ind = support_indices[i];
    a = edges[ind].end[0]; b = edges[ind].end[1];
    n = &G_s->nodelist[a];
    n->adj_objs[n->s_degree].other_end = b;
    n->adj_objs[n->s_degree].edge_index = ind;
    n->adj_objs[n->s_degree].lp_weight = m_lp_edges[ind];
    n->s_degree++;
    n = &G_s->nodelist[b];
    n->adj_objs[n->s_degree].other_end = a;
    n->adj_objs[n->s_degree].edge_index = ind;
    n->adj_objs[n->s_degree].lp_weight = m_lp_edges[ind];
    n->s_degree++;
  }

  G_s->node_count = node_count; G_s->edge_count = edge_count;

  return 0;
}
