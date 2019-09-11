#ifndef TREE_H_
#define TREE_H_

#include <list>
#include <map>
using std::map;
using std::list;

template <class T>
class tree {
  int newest;
 public:
  const static int root;
  std::map<int,T> node;
  std::map<int,list<int> > children;
  int add_leaf(int);
 tree(): newest(1) {};
};

template <class T>
const int tree<T>::root = 0;

template <class T>
inline int tree<T>::add_leaf(int parent_node){
  int current = newest++;
  children[parent_node].push_back(current);
  return current;
}

#endif /* TREE_H_*/
