from itertools import chain, combinations
import pygraphviz as pgv 
import pickle 
import pandas as pd 
import timeit 
import functools


def timeit_decorator(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        start_time = timeit.default_timer()
        result = func(*args, **kwargs)
        end_time = timeit.default_timer()
        execution_time = end_time - start_time
        print(f"Execution time of {func.__name__}: {execution_time} seconds")
        return result
    return wrapper

def powerset(iterable):
        s = list(iterable)
        return list(chain.from_iterable(combinations(s, r) for r in range(len(s)+1)))

def draw(tree, fname):
    ptree = pgv.AGraph(strict=False, directed=True)
    ptree.add_edges_from(list(tree.edges))
    ptree.layout("dot")
    ptree.draw(fname)

def concat_and_sort(mylist):
        concat_list = list(chain(*mylist))
        concat_list = sorted(concat_list, key=lambda x: x.cost)
        return concat_list
     

def get_top_n(all_trees, top_n):
        '''
        @params all_trees: list of lists of Solution Objects 
        returns the top_n minimum costs trees
        '''
    
        # concat_list = list(chain(*all_trees))
        # concat_list = sorted(concat_list, key=lambda x: x.cost)
        concat_list = concat_and_sort(all_trees)
        
        if len(concat_list) >= top_n:
            return concat_list[:top_n]
        else:
            return concat_list


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


def pickle_object(obj, file_path):
    """
    Pickle an object and save it to a file.

    Args:
    - obj: The object to pickle.
    - file_path: The file path where the pickled object will be saved.
    """
    with open(file_path, 'wb') as f:
        pickle.dump(obj, f)

def load_pickled_object(file_path):
    """
    Load a pickled object from a file.

    Args:
    - file_path: The file path from which to load the pickled object.

    Returns:
    - The unpickled object.
    """
    with open(file_path, 'rb') as f:
        obj = pickle.load(f)
    return obj

def dict_to_df(mydict, colnames):
    df = pd.DataFrame.from_dict(mydict, orient='index', columns=[colnames[1]])
    # Reset index to create a column from dictionary keys
    df.reset_index(inplace=True)
    df.columns = colnames
    return df 