for i in `seq 1 4`;
do
  ./bbcache -b -t 1 -p 1 -n 100000000 -r 0.005 ../data/CrossValidation/compas_1_train.out ../data/CrossValidation/compas_1_train.label ../data/CrossValidation/compas_1_train.minor &
done
