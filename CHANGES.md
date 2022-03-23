
extn_coeff_fasta

* 2022-03-23 bug fixed: fixed std ouput results display. "molid" changed to "mol_name" in the code and some code polishing. 

2.2
---
* 2021-12-10 bug fixed: for e(214 nm) the contribution from Trp was almost doubled. 

2.0-2.1
-------
* changes allowing multiple sequeces procesing and JSON formatting.



pI_fasta_v1.4/

1.4
-----
* changes allowing multiple sequnce processing and JSON output formatting



rdkit_pI/

Changes
=======

2022-03-23 bug fix: the output was incomplete if the extra keys (print_fragmet_pkas, etc.) were switched on.


3.2
-----

* Removed "dimorphite_site_substructures_smarts.py". Now it is loaded from the Dimorphite_DL folder. 
* New compound smiles file in the TEST folder.




smi2scrambledfasta_v1.0/




dimorphite_dl_pka/

1.2.4_for_rdkit_pI
-----
* Code modified for the use with rdkit_pI.py code that allows isoelectric point calcvulation
  of nonatural peptides. Modified functions from Dimorphite_DL are used to find the
  substructure matches in teh molecule. Moreover, the mapping between substrutures 
  and the pKa value from "site_substructures.smarts" file is used for isoelectric point calculation
* TATA core added to "site_substructures.smarts"
* Python readable file "dimorphite_dl_site_substructures_smarts.py" created
* Folder "training_data" removed 

1.2.4
-----

* Dimorphite-DL now better protonates compounds with polyphosphate chains
  (e.g., ATP). See `site_substructures.smarts` for the rationale behind the
  added pKa values.
* Added test cases for ATP and NAD.
* `site_substructures.smarts` now allows comments (lines that start with `#`).
* Fixed a bug that affected how Dimorphite-DL deals with new protonation
    states that yield invalid SMILES strings.
  * Previously, it simply returned the original input SMILES in these rare
    cases (better than nothing). Now, it instead returns the last valid SMILES
    produced, not necessarily the original SMILES.
  * Consider `O=C(O)N1C=CC=C1` at pH 3.5 as an example.
    * Dimorphite-DL first deprotonates the carboxyl group, producing
      `O=C([O-])n1cccc1` (a valid SMILES).
    * It then attempts to protonate the aromatic nitrogen, producing
      `O=C([O-])[n+]1cccc1`, an invalid SMILES.
    * Previously, it would output the original SMILES, `O=C(O)N1C=CC=C1`. Now
      it outputs the last valid SMILES, `O=C([O-])n1cccc1`.
* Improved suport for the `--silent` option.
* Reformatted code per the [*Black* Python code
  formatter](https://github.com/psf/black).

1.2.3
-----

* Updated protonation of nitrogen, oxygen, and sulfur atoms to be compatible
  with the latest version of RDKit, which broke backwards compatibility.
* Added "silent" option to suppress all output.
* Added code to suppress unnecessary RDKit warnings.
* Updated copyright to 2020.

1.2.2
-----

* Added a new parameter to limit the number of variants per compound
  (`--max_variants`). The default is 128.

1.2.1
-----

* Corrected a bug that rarely misprotonated/deprotonated compounds with
  multiple ionization sites (e.g., producing a carbanion).

1.2
---

* Corrected a bug that led Dimorphite-DL to sometimes produce output molecules
  that are non-physical.
* Corrected a bug that gave incorrect protonation states for rare molecules
  (aromatic rings with nitrogens that are protonated when electrically
  neutral, e.g. pyridin-4(1H)-one).
* `run_with_mol_list()` now preserves non-string properties.
* `run_with_mol_list()` throws a warning if it cannot process a molecule,
  rather than terminating the program with an error.

1.1
---

* Dimorphite-DL now distinguishes between indoles/pyrroles and
  Aromatic_nitrogen_protonated.
* It is now possible to call Dimorphite-DL from another Python script, in
  addition to the command line. See the `README.md` file for instructions.

1.0
---

The original version described in:

Ropp PJ, Kaminsky JC, Yablonski S, Durrant JD (2019) Dimorphite-DL: An
open-source program for enumerating the ionization states of drug-like small
molecules. J Cheminform 11:14. doi:10.1186/s13321-019-0336-9.
