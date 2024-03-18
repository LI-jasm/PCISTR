# PCISTR

> **This is the code for our paper** `Protein Complex Identification Method Based on Spatiotemporal Constrained RNA-Protein Heterogeneous Network`.

## Dependencies
Recent versions of the following packages for Python 3 are required:

* Anaconda3 (python 3.9)
* numpy==1.21.2
* torch==1.9.1
* scipy==1.7.1
* scikit-learn==0.24.2

## Reproduction

```
## Documentation
src
  │  args.py
  │  Decoupling_matrix_aggregation.py
  │  emb_train.py
  │  Model.py
  │  spatiotemporal.py
  │  Utils.py
Cluster_core_attachment.py
Compare_performance.py
MHPIN_emb.py
Update_linking_weight.py
```

## Usage
* Step 1: spatiotemporal.py: Introduce the spatiotemporal expression data of proteins from `series_matrix.txt` and `yeast_compartment_benchmark.txt`, respectively, and extract the basic interaction patterns.
* Step 2: MHPIN_emb.py: Calculate the wide-domain interaction pattern matrix and deep-domain interaction pattern matrix based on the generated basic interaction patterns. Use the pattern matrix and Protein-RNA interaction matrix as the heterogeneous network structure information, and use the dual view aggregator to jointly learn the protein embedding under the two patterns.
* Step 3: Update_linking_weight.py: Calculate the cosine similarity based on the protein embedding representation and update the weight of each edge in the PIN.
* Step 4: Cluster_core_attachment.py: Generate a set of structurally closely related protein seed cores in the PIN with the cluster mining algorithm, and filter them based on the embedding similarity. Generate complexes based on the core-attachment structure.
* Step 5: Compare_performance.py: Compare the performance of PCITSR and other classical methods based on different evaluation metrics.

## Contact

If you have any questions, please feel free to contact us!

**Zeqian Li** @github.com/LI-jasm
**Email:** [lizeqian@dlmu.edu.cn](mailto:lizeqian@dlmu.edu.cn)
**Site:** [GitHub](https://github.com/LI-jasm)
