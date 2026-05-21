import pandas as pd
import os
import sys
from multiprocessing import Pool, cpu_count
import time

OUTPUT_DIR = "/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/output/GCN4_pipeline/speed/NARDINI_retry"

def process_seq(args):
    idx, seq = args
    try:
        # Each sequence gets its own subdirectory — no collisions
        seq_dir = os.path.join(OUTPUT_DIR, f"seq_{idx:04d}")
        os.makedirs(seq_dir, exist_ok=True)
        os.chdir(seq_dir)
        
        from localcider.sequenceParameters import SequenceParameters
        SeqObj = SequenceParameters(seq)
        SeqObj.save_zscoresAndPlots(num_scrambles=100000, random_seed=None)
        return (idx, seq, "ok")
    except Exception as e:
        return (idx, seq, f"ERROR: {e}")

def main():
    exponential_fit = pd.read_csv(
        '../../../output/GCN4_pipeline/speed/exponential_fit.csv', index_col=0
    )
    
    seqs = list(enumerate(exponential_fit["ADseq"]))
    n_workers = cpu_count()
    print(f"Using {n_workers} workers for {len(seqs)} sequences", flush=True)
    
    start = time.time()
    with Pool(processes=n_workers) as pool:
        for result in pool.imap_unordered(process_seq, seqs, chunksize=1):
            idx, seq, status = result
            if status != "ok":
                print(f"[{idx}] {seq[:20]}... -> {status}", flush=True)
            elif idx % 10 == 0:
                elapsed = time.time() - start
                print(f"[{idx}/{len(seqs)}] done ({elapsed:.0f}s elapsed)", flush=True)
    
    print(f"All done in {time.time() - start:.1f}s", flush=True)

if __name__ == "__main__":
    main()