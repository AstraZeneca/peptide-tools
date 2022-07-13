

for v in $(cat test.smarts | sed "s: ::g") ; do

    vv=$( echo $v | sed 's:\"::g' )
    smi1=$( echo $vv | awk -F ',' '{print $1}') 
    smi2=$( echo $vv | awk -F ',' '{print $2}') 
    smi3=$( echo $vv | awk -F ',' '{print $3}') 
    name=$( echo $vv | awk -F ',' '{print $4}') 

    echo $name

    p=$PWD
    mkdir $name
    cd $name
        echo $smi1 > smi1.smi 
        echo $smi2 > smi2.smi 
        echo $smi3 > smi3.smi 
        ../run_pI.sh smi1.smi smi1.out
        ../run_pI.sh smi2.smi smi2.out
        ../run_pI.sh smi3.smi smi3.out
    cd $p

done

