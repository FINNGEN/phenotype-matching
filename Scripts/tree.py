from typing import AbstractSet, List, Dict, Any

class Tree(object):
    def __init__(self, name: Any, data: Any) -> 'Tree':
        self.child=[]
        self.name=name
        self.data=data
    
    def add_child(self,ch):
        self.child.append(ch)

    def get_children(self) -> List['Tree']:
        return self.child

def get_tree_nodes(tree) -> Dict[Any,Any]:
    tempstorage = []
    for c in tree.get_children():
        tempstorage.append(get_tree_nodes(c) )
    if not tempstorage:#no children
        return {tree.name:tree.data}
    out={}
    out[tree.name]=tree.data
    for s in tempstorage:
        out.update(s)
    return out