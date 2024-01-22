import alignmentsTemplate
import alignSAbDab
import alignOAS
from _02_rmsd_pairwise import run_dtw
from datetime import datetime

def combi_align(chain, cdr, length, num):
    alignmentsTemplate.template(chain, cdr, length)
    alignSAbDab.sabdab_align(chain, cdr, length, num)
    alignOAS.oas_align(chain, cdr, length, num)

def run_pipeline(chain, cdr, length, num):
    # Align is controlled for those structures that are already plotted.
    print("XXX - Running Align =", datetime.now())
    combi_align(chain, cdr, length, num)

    # DTW simply reads the files made from the align step.
    print("XXX - Running DTW =", datetime.now())
    run_dtw(chain, cdr, length)
    print("XXX - Finish DTW =", datetime.now())

def run_experiments(chain_cdr_length_combinations, n_samps):
    for chain, cdr, length in chain_cdr_length_combinations:
        run_pipeline(chain, cdr, length, n_samps)

if __name__ == '__main__':
    n_samps = 4200#0
    chain_cdr_length_combinations = [
        # ('L', 1, 6),
        # ('L', 1, 7),
        # ('L', 1, 8),
        # ('L', 1, 9),
        # ('L', 1, 11),
        # ('L', 1, 12),

        # ('L', 2, 3),

        # ('L', 3, 8),
        # ('L', 3, 9),
        ('L', 3, 10)#,
        # ('L', 3, 11),

        # ('H', 1, 8),
        # ('H', 1, 9),
        # ('H', 1, 10),

        # ('H', 2, 7),
        # ('H', 2, 8),
        # ('H', 2, 10)
    ]
    run_experiments(chain_cdr_length_combinations, n_samps)
