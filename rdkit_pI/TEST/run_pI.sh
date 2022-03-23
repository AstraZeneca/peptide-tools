# example usage. currently set to run from teh "rdkit_pI" folder.
# Make sure to have RDkit avaialble. 
# Make sure to set teh proper path to dimorphite_dl_pka (see below) - if intent to use it.
# Or make sure to have ACD Perceptabat avaialble - if intent to use it. 

module load rdkit
#module load acdperceptabatch
export PYTHONPATH=$PYTHONPATH:${pwd}/../dimorphite_dl_pka

i=$1
o=$2

usage='USAGE: ./run_pI.sh inputfile outpiutfile'
if [ -z "$i" ]; then echo "no input.  $usage" ; fi
if [ -z "$o" ]; then echo "no output.  $usage" ; fi

#python ../rdkit_pI.py -i ${i} --plot_titration_curve --print_fragment_pkas --use_acdlabs &> $o
#python ../rdkit_pI.py -i ${i} --plot_titration_curve --print_fragment_pkas --use_dimorphite --json  &> $o
#python ../rdkit_pI.py -i ${i} --plot_titration_curve --print_fragment_pkas --use_dimorphite  &> $o
python ../rdkit_pI.py -i ${i} --plot_titration_curve --print_fragment_pkas --use_dimorphite 


