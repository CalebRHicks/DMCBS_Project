#batch submit

#echo -n "Node number: "
#read node

node=12
echo "running on node $node"

date
`which mpirun` -n 24 -host cop$node ./dmcbcsmpi <dmcbcs.dat> "out.001"
date

nmax=2
i=1
while [ "$i" -lt "$nmax" ]
do
    ((i +=1))
    padded=$(printf "%03d" $i)
    date
    `which mpirun` -n 24 -host cop$node ./dmcbcsmpi <dmcbcs.dat> "out.$padded"
    date
done
