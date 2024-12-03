import os
import pyBigWig
import pysam
import pandas as pd
import numpy as np
import subprocess as sp

from typing import Optional
from scipy.signal import argrelextrema

def call_peaks(wt_bam_dir: str, save_dir: str, input_file: bool,
               ctrl_bam_dir: Optional[str] = None) -> None:

    if input_file:
        ctrl_bam_dir = 'sample_data/ChIPseq/ctrl_BAM'

    macs_env_path = '/Users/gis/Documents/venvs/macs3_env/bin'
    macs_env = os.environ.copy()
    macs_env['PATH'] += os.pathsep + macs_env_path

    wt_bam_files = sorted([x for x in os.listdir(wt_bam_dir) if
                           x.endswith('.bam')])
    if input_file:
        ctrl_bam_file = [x for x in os.listdir(ctrl_bam_dir) if
                         x.endswith('.bam')][0]

    if not os.path.exists(save_dir):
        os.mkdir(save_dir)

    for bam in wt_bam_files[0:1]:
        bam_name = bam.split('.')[0]

        print('Processing:', bam_name)

        macs_cmd = (f'macs3 callpeak '
                    f'-t {wt_bam_dir}/{bam} '
                    f'-f BAMPE '
                    f'--outdir {save_dir} '
                    f'--name {bam_name} '
                    )

        if input_file:
            macs_cmd += f'-c {ctrl_bam_dir}/{ctrl_bam_file}'

        sp.call(macs_cmd, shell=True, env=macs_env)

    print('Done: peak calling...')


def find_max():
    return 1

def expand_summit():
    return 1