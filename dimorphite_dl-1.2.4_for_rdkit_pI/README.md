dimorphite_dl-1.2.4_for_rdkit_pI
===================

What is it?
-----------

It is custom spin off version of Dimorphite-DL 1.2.4 software that allows finding 
ionizable centers in the 2D structure of the molecule. It also provides with the set
of pKa values for all defined ionizable centres. 

Citation
--------

If you use Dimorphite-DL in your research, please cite:

Ropp PJ, Kaminsky JC, Yablonski S, Durrant JD (2019) Dimorphite-DL: An
open-source program for enumerating the ionization states of drug-like small
molecules. J Cheminform 11:14. doi:10.1186/s13321-019-0336-9.

Check this for the oridinal code:
https://github.com/UnixJunkie/dimorphite_dl

Licensing
---------

Dimorphite-DL is released under the Apache 2.0 license. See LICENCE.txt for
details.

Usage
-----

It can be used from a Python code. Example: 

def calc_pkas_dimorphite_dl(smi_list):

    from dimorphite_dl_site_substructures_smarts import data_txt
    from rdkit import Chem
    from dimorphite_dl import Protonate

    skip_site_names = ['TATA','*Amide']

    base_pkas=[]
    acid_pkas=[]
    diacid_pkas=[]

    for smiles in smi_list:
        protonation_sites = Protonate({'smiles':smiles}).get_protonation_sites()

        for sites in protonation_sites:
            #(3, 'BOTH', '*Amide') 
            site_name = sites[2]

            if site_name in skip_site_names: continue

            if site_name in D_dimorphite_dl_type_pka.keys():
                site_data = D_dimorphite_dl_type_pka[site_name]
                if site_data['type']=='base':
                        base_pkas.append((site_data['pka'],smiles))
                if site_data['type']=='acid':
                        acid_pkas.append((site_data['pka'],smiles))
                if site_data['type']=='diacid':
                        diacid_pkas.append(((site_data['pka1'],site_data['pka2']),smiles))
                
            else:
                print('Error:  not known dimorphite site type ' + site_type)
                sys.exit(1)



Caveats
-------

Dimorphite-DL deprotonates indoles and pyrroles around pH 14.5. But these
substructures can also be protonated around pH -3.5. Dimorphite does not
perform the protonation.

Authors and Contacts
--------------------

See the `CONTRIBUTORS.md` file for a full list of contributors. Please contact
Andrey Frolov (andrey.frolov@astrazeneca.com) with any questions.

