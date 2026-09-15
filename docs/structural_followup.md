# Structural Modeling & Virtual Screening Protocol

Once the communication analysis nominates high-priority ligand–receptor pairs, the following steps can be performed on GPU-enabled infrastructure (Google Colab Pro, institutional HPC, or cloud).

## 1. Complex structure prediction

1. Retrieve protein sequences from UniProt (prefer mature / extracellular domains for transmembrane receptors).
2. Submit to **AlphaFold3** (alphafoldserver.com) or **ColabFold** multimer mode.
3. Download the top-ranked model (PDB/CIF).
4. Inspect the interface in PyMOL or ChimeraX; identify key contact residues.

## 2. Binding-site definition

- Use the predicted complex to define a docking box that covers the ligand–receptor interface or the orthosteric pocket of the receptor.
- Prepare the receptor (remove ligand atoms if co-folded, add hydrogens, assign charges) with AutoDockTools or OpenBabel → PDBQT.

## 3. Compound library

Recommended starting libraries (all publicly available):

- FDA-approved drugs (DrugBank, ZINC FDA subset, SWEETLEAD)
- Bioactive / probe libraries (ChEMBL, BindingDB)

Convert to PDBQT format.

## 4. Virtual screening

```bash
vina --receptor receptor.pdbqt \
     --ligand  library/ \
     --out     results/ \
     --exhaustiveness 16 \
     --dir     docking_out/
```

Rank poses by predicted binding affinity. Prioritize compounds that:

- Occupy the ligand-binding pocket, or
- Form steric clashes with critical interface residues of the complex.

## 5. Experimental follow-up

Top candidates should be tested in:

- Co-culture assays (leukemic cells + MSCs / endothelial cells)
- Migration / adhesion / survival read-outs under niche-protective conditions
- In vivo niche-homing or residual-disease models

## Priority pairs suggested by the demo pipeline

Typical high-scoring axes that emerge from the communication module:

| Ligand  | Receptor | Biological rationale                  |
|---------|----------|---------------------------------------|
| CXCL12  | CXCR4    | Classic HSC / LSC retention axis      |
| VEGFA   | KDR/FLT1 | Angiogenic niche remodeling           |
| IL6     | IL6R     | Inflammatory / survival signaling     |
| CSF1    | CSF1R    | Myeloid niche support                 |
| TGFB1   | TGFBR1/2 | Quiescence / immune modulation        |
| SPP1    | CD44     | Adhesion / chemoresistance            |
