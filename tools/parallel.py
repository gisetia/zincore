from typing import Callable, Iterable, List, Dict, Optional, Union
import pandas as pd
import numpy as np
import multiprocessing as mp
from functools import partial


def groupby_apply_parallel(grouped_df, func: Callable,
                           cores: int = mp.cpu_count(),
                           maxtasksperchild: Optional[int] = None, **kwargs):
    '''Paralellize pandas df.groupby.apply'''

    partial_ = partial(func, **kwargs)

    with mp.Pool(cores, maxtasksperchild=maxtasksperchild) as p:
        returned_list = p.map(partial_, [group for name, group in grouped_df])

        if type(returned_list[0]) == pd.Series:

            idx1 = returned_list[0].index

            if np.all([idx1.equals(x.index) for x in returned_list]):
                # If all series share an index, concatenate along axis 1
                # Shared indices will then become dataframe columns
                returned_concat = pd.concat(returned_list, axis=1).T
                returned_concat.index = [name for name, group in grouped_df]
                returned_concat.index.name = grouped_df.keys
            else:
                # If series have different indices, concatenate them as rows
                print(x.shape for x in returned_list)
                returned_concat = pd.concat(returned_list)
                returned_concat.index = [name for name, group in grouped_df
                                         for x in range(0, len(returned_list))]
        elif type(returned_list[0]) == pd.DataFrame:
            names = [name for name, group in grouped_df]

            if len(names[0]) == 1:
                for idx, df in enumerate(returned_list):
                    df.index = [names[idx] for x in range(len(df))]
                returned_concat = pd.concat(returned_list)
                returned_concat.index.name = grouped_df.keys
            else:
                for idx, df in enumerate(returned_list):
                    df[grouped_df.keys] = [(names[idx][0], names[idx][1])
                                           for x in range(len(df))]
                returned_concat = pd.concat(returned_list)
                returned_concat = returned_concat.set_index(grouped_df.keys)
        else:
            returned_concat = pd.Series(returned_list,
                                        index=[name for name, group
                                               in grouped_df])
            returned_concat.index.name = grouped_df.keys

    return returned_concat


def func_parallel(func: Callable, arg_list: Iterable,
                  cores: Optional[int] = mp.cpu_count(),
                  maxtasksperchild: Optional[int] = None, return_dict=True,
                #   chunk_size=1,
                  **kwargs) -> Union[list, dict]:
    '''Parallelize function 'func', given a set of arguments 'arg_list'.
    Returns dict where keys are the arguments and values are the function
    applied to the given arguments.
    '''
    partial_ = partial(func, **kwargs)

    if cores == 1:
        returned_list = list(map(partial_, arg_list))
    else:
        with mp.Pool(cores, maxtasksperchild=maxtasksperchild) as p:
            returned_list = p.map(partial_, arg_list)

    if return_dict:
        returned_dict = dict(zip(arg_list, returned_list))
        return returned_dict
    else:
        return returned_list


def func_parallel_lazy(func: Callable, arg_list: Iterable,
                       cores: Optional[int] = mp.cpu_count(),
                       maxtasksperchild: Optional[int] = None,
                       #   lazy=False,
                       chunk_size=1, **kwargs):

    partial_ = partial(func, **kwargs)
    with mp.Pool(cores, maxtasksperchild=maxtasksperchild) as p:
        returned_list = p.imap(partial_, arg_list, chunk_size)

        # returned_list = list(returned_list)
        # return returned_list

        for i in returned_list:
            yield i


def chunks(generator, chunk_size):
    '''Yield successive n-sized chunks from generator.'''
    chunk = []
    for item in generator:
        chunk.append(item)
        if len(chunk) == chunk_size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk
