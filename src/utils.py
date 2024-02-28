from itertools import chain, combinations
import pygraphviz as pgv 


def powerset(iterable):
        s = list(iterable)
        return list(chain.from_iterable(combinations(s, r) for r in range(len(s)+1)))

def draw(tree, fname):
    ptree = pgv.AGraph(strict=False, directed=True)
    ptree.add_edges_from(list(tree.edges))
    ptree.layout("dot")
    ptree.draw(fname)


def merge_lists(list1, list2):
        # Base case: if either list is empty, return the other list as the only merging
        if not list1:
            return [list2]
        if not list2:
            return [list1]

        # Recursive case:
        # For each possible merging, recursively merge the remaining parts of the lists
        merged_lists = []
        for merged_part1 in merge_lists(list1[1:], list2):
            merged_lists.append([list1[0]] + merged_part1)
        for merged_part2 in merge_lists(list1, list2[1:]):
            merged_lists.append([list2[0]] + merged_part2)
        return merged_lists


